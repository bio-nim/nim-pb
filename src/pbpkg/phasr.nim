# vim: sw=4 ts=4 sts=4 tw=0 et:
import hts
import os
import ./kmers

proc foo*() =
    echo "foo"
#[
proc readaln*(bfn: string; fasta: string) =
 var b:Bam

 #open reference fasta with htsnim
 reference = open_fasta(fasta)

 open(b, bfn, index=false)
 echo "[INFO] reading bam"
 for record in b:
  var rseq: string
  discard sequence(record, rseq)
  var kmers: pot_t = dna_to_kmers(rseq, 21)
  
  # for each overlapping read
  # complement to remove kmers that are in the reference where the read maps. 
  rseq = get_ref(reference, record.ref, record.start-20, record.end+20)
  var refkmers = dna_to_kmers(rseq, 21)
  complement(kmers, refkmers)
  
  # find all read that overlap current read

  # count shared kmers between all overlapping reads

  # build weighted graph were nodes are reads and the edges are overlaps, weights are the number of shared kmers between reads

  edge between two nodes:
  - read id1, read id2
  - weight
  node:
  - read id
  - read start
  - read end

  # run unknown algorithm

  - merger function between nodes
  - dynamic programming

  ouput:   
  	   tuples:
   	   - phase block id [0,1,2 ... ]
	   - phase block start
	   - phase block end
   	   - phase [0,1]
   	   - vector of read ids
]#

proc main*(ref_fn: string, aln_fn: string) =
    echo "[INFO] input data:", ref_fn, " ", aln_fn
    #readaln(aln)

when isMainModule:
    main()
