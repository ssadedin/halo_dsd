import net.sf.samtools.SAMRecord;

/**
 * Combines information about HaloPlex amplicons and reads together to offer integrated
 * metrics and functions.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class HaloPlexBAM {
    
    Regions amplicons
    
    SAM reads
    
    SAM readMates
    
    /**
     * Optionally, reads that do not align to an expected position can be looked up in
     * a hash based index. This allows some read counts to be more accurate.
     */
    FASTAIndex ampliconIndex = null
    
    boolean verbose = false
    
    def log = System.err
    
    public HaloPlexBAM() {
    }

    public HaloPlexBAM(Map options=[:], Regions amplicons, SAM reads, FASTA reference) {
        this(reads,amplicons)
        this.ampliconIndex = 
            new FASTAIndex(reference, 0..(options.maxOffset?:5), 0, options.seedSize?:20, amplicons)
    }
    
    public HaloPlexBAM(Regions amplicons, SAM reads) {
        this.reads = reads
        this.amplicons = amplicons
    }
    
    /**
     * Locate the amplicon that the given read belongs to. For improvement in 
     * efficiency, the mate pair of read r can be passed as well. Otherwise it
     * will be queried, however this option is quite slow.
     */
    List<Amplicon> findAmplicons(SAMRecord r, IntRange offsets=0..5, SAMRecord r2=null) {
        
        if(readMates == null)
            readMates = new SAM(reads.samFile)
        
        // First look at amplicons starting at exactly the position
        def amplicons = []
        if(!r.readNegativeStrandFlag) {
          for(int pos in [r.alignmentStart-offsets.from-1,r.alignmentStart-offsets.to-1]) {
              int oldSize = amplicons.size()
              amplicons.addAll(this.amplicons.startingAt(r.referenceName,pos))
              
              if(verbose)
                log.println "Read $r.readName at $r.referenceName:$r.alignmentStart matched ${amplicons.size()} amplicons at $pos"
           }
          
          // If ambiguous, filter by ones best matching the other end
          if(amplicons.size()>1) {
            if(r2 == null) {
                try {
                  r2 = readMates.samFileReader.queryMate(r)
                }
                catch(net.sf.samtools.SAMFormatException ex) {
                    // Get this with some rare cases because bwa mem maps the mate to multiple "primary" positions,
                    // which Picard does not like
                }
            }
            if(r2)
              amplicons = amplicons.grep { r2.alignmentEnd in [it.to+offsets.from+1,it.to-(offsets.to-1)] }
          }
            
          if(!amplicons && verbose)
             log.println "Unmatched: $r.readName at $r.referenceName:$r.alignmentStart (mapQ = $r.mappingQuality)"
        }
        else {
          for(int pos in [r.alignmentEnd-offsets.from-1,r.alignmentEnd+offsets.to-1]) {
            amplicons.addAll(this.amplicons.endingAt(r.referenceName,pos))
            
            if(verbose)
                log.println "Read $r.readName at $r.referenceName:$r.alignmentEnd matched ${amplicons.size()} amplicons at $pos"
          }
          
          if(amplicons.size()>1) {
              amplicons = amplicons.grep { r.mateAlignmentStart in [it.from+offsets.from+1,it.from+offsets.to+1] }
          }
          
          if(!amplicons && verbose)
             log.println "Unmatched: $r.readName at $r.referenceName:$r.alignmentStart (mapQ = $r.mappingQuality)"
        }
        
        return amplicons
    }
    
    /**
     * Annotates the amplicons for this BAM with counts of the reads that map to each
     * amplicon.
     * 
     * @return
     */
    Map countReads(Map options = [:]) {
        ProgressCounter counter = options.progress 
        
        boolean verboseCount = verbose || options.verbose
        if(counter == true) 
            counter = new ProgressCounter()
            
        FASTAIndex index = null
        if(options.index)
            index = options.index
            
        int count = 0
        int countAmbiguous = 0
        int countUnmatched = 0        
        int countRescuedByIndex = 0        
        int countMatched = 0
        reads.eachPair { SAMRecord r1, SAMRecord r2 ->
            
           if(options.includeChrs && !options.includeChrs.contains(r1.referenceName))
               return
               
           if(options.excludeChrs && options.excludeChrs.contains(r1.referenceName))
               return
               
            if(verboseCount)
               log.println "Read $r1.readName starting at $r1.alignmentStart, ending at $r1.alignmentEnd"
               
            if(counter)
               counter.count()
               
            ++count
            
            def amplicons = findAmplicons(r1, 0..5, r2)
            if(amplicons.size()>1) {
                log.println "${amplicons.size()} amplicons remaining after matching mate end - multiple legit amplicons match read $r1.readName"
                ++countAmbiguous
            }
            
            if(amplicons) {
                amplicons.each { GRange a ->
                    countReadToAmplicon(a,r1)
                }
                ++countMatched
                return
            }
            
            // If we have got here we did not find any amplicon at the correct position
            // If the user provided a FASTA index, we can attempt lookup of the amplicon in there
            if(index) {
//                log.println "Checking index for read ${r1.readString.substring(0,30)+'...'}"
//                log.println "Checking index for read1 ${r1.readString}"
                String ampliconName = index.querySequenceName(r1.readString)
                if(!ampliconName) {
//                    log.println "Checking index for read2 ${r2.readString}"
                    ampliconName = index.querySequenceName(r2.readString)
                }
                if(ampliconName) {
                    // Find the amplicon - they are named in the standard format, chr:start-end
                    Region amplicon = new Region(ampliconName)
                    List candidates = this.amplicons.startingAt(amplicon.chr, amplicon.from)
                    if(candidates.size() == 1) {
                        ++countRescuedByIndex
                        countReadToAmplicon(candidates[0], r1)
                        if(verboseCount)
                            log.println "Read $r1.readName rescued by hash lookup lookup to amplicon $ampliconName"
                        return
                    }
                    else {
                        if(verboseCount) {
                            log.println "WARNING: Read $r1.readName multimaps to amplicons ${candidates*.toString()}"
                            ++countAmbiguous
                        }
                    }
                }
            }
            
            ++countUnmatched
            // Not found - add to unmatched
            if(verboseCount)
                System.err.println "Unmatched: $r1.readName at $r1.referenceName:${r1.readNegativeStrandFlag?r1.alignmentEnd:r1.alignmentStart} (mapQ = $r1.mappingQuality)"
        }
        
        return [matched: countMatched, ambiguous: countAmbiguous, rescued: countRescuedByIndex, unmatched: countUnmatched, total: count]
    }
    
    /**
     * Increment 'start' or 'end' attributes to the given GRange, depending on read orientation
     */
    private countReadToAmplicon(GRange a, SAMRecord r) {
        if(!a.extra || (a.extra instanceof String)) { a.extra = [start:0, end:0] };
        if(r.readNegativeStrandFlag)
            ++a.extra.start
        else
            ++a.extra.end 
    }
}
