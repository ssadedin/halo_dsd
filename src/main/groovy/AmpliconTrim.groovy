import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.util.SequenceUtil;

import graxxia.IntegerStats
import java.util.zip.*

import org.biojava3.alignment.template.Aligner;

/**
 * A trimmer that uses the known amplicon boundaries and sizes to improve trimming 
 * accuracy and avoid trimming non-adapter sequence from the middle of amplicons.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class AmpliconTrim {
    
    /**
     * Command line options
     */
    def opts = [:] // default to empty map for unit tests
    
    /**
     * Default sequence used for matching adapters
     */
    static final DEFAULT_ADAPTER = "AGATCGGAAGAG"
    
    /**
     * All possible permutations of adapter sequences that we will check for
     */
    List<String> adapterSeqs     
      
    /**
     * Index olding references to FASTA for all the amplicons
     */
    FASTAIndex index 
      
    /**
     * Count of how many reads we identify containing adapters
     */
    int adapterCount = 0

    /**
     * Count of the total number of reads processed
     */
    int total = 0

    /**
     * Count of how many trimmed reads had amplicons identified
     */
    int adapterCountWithAmplicons = 0
    
    /**
     * How many base pairs between read end and amplicon end
     * (estimated if not provided)
     */
    IntRange ampliconReadOffset = 0..5
    
    /**
     * Seed size to use for indexing amplicons
     */
    int seedSize = 30
    
    boolean verbose = false
    
    static void main(String [] args) {
      
      def err = System.err
      
      CliBuilder cli = new CliBuilder(usage:"java -jar AmpliconTrim.jar [options]\n\n", writer: new PrintWriter(System.err))
      cli.with {
          '1' args: 1, "FASTQ for read 1 (forward)"
          '2' args: 1, "FASTQ for read 2 (reverse)"
          't1' args: 1, "Output file for trimmed reads (forward)"
          't2' args: 1, "Output file for trimmed reads (reverse)"
          's' args: 1, "Alignment score threshold for considering reads to be overlapping"
          'a' args: 1, "adapter sequence to probe for at end of reads (eg: AGATCGGAAGAGCG)"
          'f' args: 1, "FASTA for amplicons"
          'r' args: 1, "Reference fasta (if not using -f)"
          'b' args: 1, "bed file of amplicons (if not using -f)"
          'o' args: 1, "maximum offset of read end from amplicon ends (5)"
          'ao' args: 1, "Offset of read ends from amplicon end (bp). Optional, if not provided, offset will be estimated from data."
          'v' "Turn on verbose mode"
          'vv' "Turn on super verbose (David Goode) mode"
      }
      
      def opts = cli.parse(args)
      ['1','2','t1','t2'].each {
          if(!opts[it]) {
              err.println()
              err.println "Please provide option -$it"
              err.println()
              cli.usage()
              System.exit(1)
          }
      }
      
      if(!opts.f && !(opts.b && opts.r)) {
          err.println()
          err.println("Please specify either -f or both of -r and -b")
          err.println()
          System.exit(1)
      }
      
      String readFile1 = opts['1']
      String readFile2 = opts['2']
      
      def trimmer = new AmpliconTrim(opts:opts)
      if(opts.v || opts.vv)
          trimmer.verbose = true
          
      trimmer.trim(readFile1, readFile2)
    }
    
    void indexAdapter(String adapter) {
         String adapterComp = FASTA.reverseComplement(adapter)
          
          // Create a list of all hte possible terminators for reads that would be adapters
          adapterSeqs = (1..adapter.length()).collect { adapter.substring(0,it) }.reverse()
          
          if(opts.vv) {
              println "All indexed adapter sequences are: $adapterSeqs"
          }
    }
    
    void indexReference(String referencePath, String ampliconBedFilePath) {
      index = new FASTAIndex(new FASTA(referencePath), 0..(opts.o ? opts.o.toInteger():5), 0, seedSize, new BED(ampliconBedFilePath)) 
      if(opts.vv) {
          println "Indexed Amplicon Sequences are: $index.sequences"
          println "Indexed Amplicon Names are: $index.sequenceNames"
      }
    }
    
    IntRange determineOffsetRange(String readFile1, String readFile2) {
        
        int identifiedCount = 0
        IntegerStats stats = new IntegerStats(50)
        FASTQ.eachPair(readFile1, readFile2) { FASTQRead r1, FASTQRead r2 ->
            def ampliconR1 = index.querySequenceName(r1.bases)
            if(ampliconR1) {
               int seedOffset = findAmpliconOffset(r1, index.sequenceNames[ampliconR1])
               stats.addValue(seedOffset)
                ++identifiedCount
            }
            
            def ampliconR2 = index.querySequenceName(r2.bases)
            if(ampliconR2) {
               int seedOffset = findAmpliconOffset(r2, index.sequenceNames[ampliconR2])
                stats.addValue(seedOffset)
                ++identifiedCount
            }
            
            if(identifiedCount>1000) 
                throw new Abort()
        }
        
        if(identifiedCount < 50) 
            throw new RuntimeException("Less than 50 reads could be mapped to amplicons. As a result, it was not possible to estimate read / amplicon offsets. Please provide them explicitly with the -ao option.")
        
        println "Median offset from amplicon start is ${stats.median}"
        if(opts.vv) 
            println "Offset Distribution: "
            
        def offsets = []
        stats.values.eachWithIndex { freq, offset ->
            if(opts.vv)
                println "$offset: $freq"
            offsets << [freq: freq, offset: offset]
        }
        
        // Pick off the maximum and second maximum
        offsets.sort { -it.freq }
        
        println "Estimated offsets of read start from amplicon start are ${offsets[0].offset} and ${offsets[1].offset}"
        return offsets[0].offset..offsets[1].offset
    }
    
    /**
     * Examine a read belonging to an amplicon in both forward and reverse (reverse complement) 
     * orientation, and return its offset within the amplicon.
     * 
     * @param r             	the read to examine
     * @param ampliconSequence  the amplicon sequence
     * @return
     */
    int findAmpliconOffset(FASTQRead r, String ampliconSequence) {
		
		// The seed that was used for lookup - it has to be inside the amplicon somewhere
        String seed = r.bases.substring(0,index.seedSize)
		
		// Look relative to start of the samplicon
        int seedOffset = ampliconSequence.indexOf(seed)
        if(seedOffset<0 || seedOffset > index.seedSize) {
			
			// Not there? look in reverse complement
            seedOffset = ampliconSequence.size() - ampliconSequence.lastIndexOf(FASTA.reverseComplement(seed)) - index.seedSize
            if(seedOffset>index.seedSize) {
                System.err.println "ERROR: seed not found in amplicon queried from index!"
                System.err.println "ERROR: seed     = $seed"
                System.err.println "ERROR: reverse  = ${FASTA.reverseComplement(seed)}"
                System.err.println "ERROR: amplicon = ${ampliconSequence}"
                throw new IllegalStateException("Bad amplicon index state")
            }
        }
        return seedOffset
    }
    
    void trim(String readFile1, String readFile2) {

        String adapter = opts.a?:DEFAULT_ADAPTER

        this.indexAdapter(adapter)

        // Index of the amplicons - we will look up the amplicons for reference, although they aren't
        // currently used in determining the overlap
        System.err.println "Reading amplicons ..."
        if(opts.f) {
            // First read in all the amplicon FASTAs
            String amplicon_fasta = opts.f
            index = new FASTAIndex(new FASTA(amplicon_fasta),0..(opts.o ? opts.o.toInteger():5))
        }
        else {
            indexReference(opts.r, opts.b)
        }
        
        System.err.println "Inspecting amplicon / read offsets"
        
        if(opts.ao) {
            this.ampliconReadOffset = 0..(opts.ao.toInteger())
        }
        else {
            this.ampliconReadOffset = determineOffsetRange(readFile1, readFile2)
        }
        
//        System.exit(0)

        System.err.println "Inspecting reads ..."

        Writer trimOut1 = new PrintWriter(new GZIPOutputStream(new File(opts.t1).newOutputStream()))
        Writer trimOut2 = new PrintWriter(new GZIPOutputStream(new File(opts.t2).newOutputStream()))

        FASTQ.eachPair(readFile1, readFile2) { FASTQRead r1, FASTQRead r2 ->

            ++total

            try {
                // Check for reverse orientation (read 1 is on alternate strand)
                (r1,r2) = checkReadPair(r1,r2)
    
                // Check for forward orientation (read 1 is on forward strand)
                (r2,r1) = checkReadPair(r2,r1)
            }
            catch(Exception e) {
                println "WARNING: failed to process read $r1.name"
                e.printStackTrace()
            }

            r1.write(trimOut1)
            r2.write(trimOut2)
        }

        trimOut1.close()
        trimOut2.close()

        println " Summary ".center(80,"=")
        println "$adapterCount / $total reads contained adapter sequence at end"
        println "$adapterCountWithAmplicons / $adapterCount of reads containing adapters matched expected amplicon"
        println "=" * 80
      }  
    
    /**
     * Check if the second read ends with adapter sequence and the first
     * read aligns to an amplicon boundary. If so, trim the reads to remove
     * the adapter sequence and return a new read pair reflecting the 
     * trimmed reads.
     * 
     * @param r1    First read in pair
     * @param r2    Second read in pair
     * @return
     */
    List<FASTQRead> checkReadPair(FASTQRead r1, FASTQRead r2) {
        
      // Does the read end with adapter sequence?
      List trimmedReads = null
      for(adapterSeq in adapterSeqs) { 
          int index = r2.bases.lastIndexOf(adapterSeq)
//          println "Index of $adapterSeq = $index in $r2.bases"
          if(index == (r2.bases.size() - adapterSeq.size())) { // adapter sequence at end
              trimmedReads = checkAdapterForReadPair(adapterSeq, r1, r2)
          }
          else
          if(adapterSeq.size() > 8 && index > 20) { // if long enough, look for it in the  middle
              trimmedReads = checkAdapterForReadPair(adapterSeq, r1, r2)
          }
          if(trimmedReads)
              break
      }
      
      if(trimmedReads)
          return trimmedReads
      else
          return [r1,r2]
    }
    
    /**
     * Check the given read pairs for match with the specified adapter sequence and
     * a matching amplicon
     * 
     * @param endingAdapter
     * @param r1    Read 1 in the pair
     * @param r2    Read 2 in the pair
     * @return  if adapter sequence was identified, a list containing a new R1 and R2. If no
     *          adapter sequence was identified, null
     */
    List<FASTQRead> checkAdapterForReadPair(String endingAdapter, FASTQRead r1, FASTQRead r2) {
        
        if(r1.bases.size() != r2.bases.size()) {
            if(verbose)
                println "Skipping trim of $r1.readName because r1 and r2 are different lengths (${r1.bases.size()} vs ${r2.bases.size()}). Assume already trimmed?"
            return null
        }
            
       
        int adapterPos = r2.bases.lastIndexOf(endingAdapter)
        int adapterLength = r2.bases.size() - adapterPos
            
        // Ok, so check if the two reads match after a 5bp gap
        String reverseCompR2 = FASTA.reverseComplement(r2.bases)
        int adjoinMax = Math.min(reverseCompR2.size()-1, adapterLength + ampliconReadOffset.to+10)
        String adjoiningBases1 = reverseCompR2.substring(adapterLength + ampliconReadOffset.to, adjoinMax)
        String adjoiningBases2 = reverseCompR2.substring(adapterLength + ampliconReadOffset.from, 
                                                         adapterLength + ampliconReadOffset.from + 10)
        
        if(opts.ignoren) {
            
        }
        
        boolean match1 = r1.bases.startsWith(adjoiningBases1)
        boolean match2 = r1.bases.startsWith(adjoiningBases2)
        
        if(verbose) {
            println "Adapter match: $endingAdapter (index=$adapterPos, length=$adapterLength) with adjoining bases: $adjoiningBases1, $adjoiningBases2 (match = $match1,$match2)"
        }
        
        if(match1 || match2) {
                
            // Look up r1 in the amplicon index
            String ampliconStart = r1.bases.substring(0,seedSize)
            
            if(opts.vv)
                println "Looking up amplicon for sequence : $ampliconStart"
                
            String amplicon  = index.sequences[ampliconStart]
                
            if(amplicon) {
                
                if(verbose)
                    println "Amplicon matching read identified at $ampliconStart: $amplicon"
                ++adapterCountWithAmplicons
            }
                
            if(r1.bases.size() == 0) {
                System.err.println "WARNING: empty bases for read $r1.name?!"
            }
            
            if(adapterLength<0) {
                System.err.println "WARNING: adapterLength = $adapterLength for read $r1.name?!"
            }
            
            // Check with alignment score
            // Here only align to read body
            def alignment = Align.global(r1.bases.substring(0,r1.bases.size()-adapterLength), reverseCompR2.substring(adapterLength))
            
            if(opts.vv) {
                println "Align:\n$r1.bases\n$reverseCompR2"
                println "Result (score=$alignment.score):\n" + alignment.profile.toString()
                println "Distance = " + alignment.getDistance()
            }
            
            if(alignment.score > 400 || ((r1.bases.size() - adapterLength < 100) && (alignment.distance < 0.1))) {
              println "=" * 80
              def ampliconSize = (index.sequenceNames[amplicon] ? index.sequenceNames[amplicon].size() : "?")
              println "Adapter sequence match at end of $r2.name: $endingAdapter (${adapterLength} bases) pair alignment score = $alignment.score amplicon ${amplicon?:'N/A'} size = $ampliconSize"
              
              // We realign the whole read ONLY so as to display the final result accurately in the output
              def align2 = Align.global(r1.bases, reverseCompR2)
              println(("-" * (adapterLength)) + index.sequenceNames[amplicon]?:"")
              print align2.profile.toString()
              println("^" * adapterLength)
              ++adapterCount
              
              // Note: trim extra 5 bases when there is a 5 base offset from amplicon start
              FASTQRead r1TrimEnd = r1.trimEnd(match1 ? adapterLength+ampliconReadOffset.to : adapterLength)
              FASTQRead r2TrimEnd = r2.trimEnd(adapterLength)
              return [r1TrimEnd,r2TrimEnd]
            }
            else { 
                if(opts.v)
                    println "Read $r1.readName failed to trim due to insufficient alignment of adapter / read hybrid"
                return null
            }
        }
        else {
           if(verbose) {
               println "Other read starting bases do not match:\nExpected R1: $adjoiningBases2\nActual R1  : ${r1.bases.substring(0,Math.max(endingAdapter.size()*2,6))}"
               println " " * 13 + "^" * endingAdapter.size()
           }
           return null
        }
        assert false, "Should not get here"
    }
}
