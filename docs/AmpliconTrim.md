# Amplicon Based Trimming

## General Trimming Approaches

The goal of trimming methods is to remove non-biological sequences from HTS
reads to avoid those sequences appearing as artefactual variant calls and
having other adverse effects on downstream processing of the data. Since most
contamination usually arises from adapters used in the library preparation
steps leading up to sequencing, the most common approach implemented by tools
is to perform an alignment between known contamination (adapter) sequences and
the bases observed in a read. When a sufficiently close match is observed over
a long enough segment, the matching part of the sequence is interpreted as
contaminant and the corresponding bases are removed. 

## Trimming Approaches Specialised for Read-Through

The aforementioned methods are designed to detect contamination which may
theoretically appear anywhere within a read. However, one particular case
stands above others in importance: the case of "read through" into adapter
sequence at the 3' end of reads. In paired end sequencing, a characteristic
signature allows for an optimised method to detect read-through adapter
contamination. In Illumina libraries, a long universal adapter sequence is
ligated to each end of every fragment. When read-through occurs, a symmetry
appears in the arrangement of the adapter sequence and the body of the
fragment: each read contains identical sequence (after reverse complementation)
representing the entirety of the true fragment. This is then followed directly
by the known adapter sequence.  Recognising the identical sequence in the two
reads allows for a stringent test for the presence of read-through adapter
contamination Trimmomatic, for example, has implemented a mode named
"ILLUMINACLIP" to specifically take advantage of this read configuration.
Trimmomatic implements the test by appending the known adapter contamination
sequence to each end of every read.  It then performs an alignment of the
resulting sequence of R1 and R2 for each read pair. 

## Our Approach

Amplicon based data such as HaloPlex presents even more prior information than
just the known read configuration mentioned above: we know in advance the
predicted start and end of each fragment on the sequencer. Thus even more power
for detecting contamination can be obtained if knowledge of the predicted
mapping locations and the known adapter sequence is combined together. This
idea forms the inspiration for the method we have implemented, called
"AmpliconTrim".  The AmpliconTrim algorithm begins with a file containing the
known amplicons in the capture and the known adapter sequence, A, as well as
the sequencing reads themselves. The algorithm then proceeds as follows:

 1. Initially, every amplicon is indexed by the first and last 10 bases
 2. The adapter sequence A is split into a list containing every prefix of A
 3. Similarly, every prefix of the reverse complement of A is also added to the list
 4. The entire prefix list is then sorted in order of descending length
 5. Each read pair consisting of R1 and R2 is then processed in turn. For each pair:
    i.   the adapter sequence prefixes are taken in order of decreasing size
    ii.  each adapter sequence prefix is searched within the R1 to find the last index at which 
         the prefix occurs (if it is found)
    iii. the first match that is identified terminates the search, resulting in the longest matching
         prefix, and the index at which it occurred within the read
    iv.  if no match is identified, the read pair is left untrimmed
    v.   the putative adapter sequence is removed from both R1 and R2, including all subsequent bases
    vi.  after removal of adapter sequence, the body of the reads is expected to be identical 
         after taking the reverse complement of R2. The first 10 bases (by default) are tested for a
         match between the start of R1 and the end of R2
    vii. if no match exists, steps 5.ii - 5.vi are repeated, searching for progressively shorter
         adapter sequence matches
    viii. the first 10 bases of the body of the reads is searched in the amplicon index. If no amplicon
          exists starting with the read body, trimming is aborted for the read pair
    ix.  Finally, a full Needleman-Wunsch alignment is performed to align the read body with the
         identified amplicon body. If a sufficient alignment score is observed, the trimmed body of the 
         reads is written to the output instead of the original reads.

In steps 5(i) - 5(vii) this algorithm is approximately equivalent to that
implemented by Trimmomatic. The key difference is in the final two steps where
reads are only trimmed when there is a concordance between the trimmed amplicon
and a known sequence expected from the HaloPlex design.  By constraining
trimming to these instances, the trimming algorithm can be made highly
stringent without losing sensitivity. As a result, even contamination
consisting of a single base can be trimmed without resulting in accidental
trimming of non-adapter sequence. 

## Implementation

The implementation can be found in ../src/main/groovy/AmpliconTrim.groovy

