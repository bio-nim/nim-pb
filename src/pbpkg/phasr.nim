# vim: sw=4 ts=4 sts=4 tw=0 et:
from strutils import format
import hts
import os
import ./kmers

proc foo*() =
    echo "foo"

proc readaln*(bfn: string; fasta: string) =
    const klen = 21
    var b: hts.Bam

    hts.open(b, bfn, index=true)
    echo "[INFO] reading bam"
    for record in b:
        echo format("$# $# ($#) [$# .. $#] $#", record.tid, record.chrom, record.qname, record.start, record.stop,
            ($record.cigar).substr(0, 32))
        var rseq: string
        discard hts.sequence(record, rseq)
        var kmers: pot_t = dna_to_kmers(rseq, klen)
 #[

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

proc main*(aln_fn: string, ref_fn: string) =
    echo "[INFO] input reference (fasta):", ref_fn, ", reads:", aln_fn
    if strutils.find(ref_fn, "fa") == -1:
        echo format("[WARN] Bad fasta filename? '$#'", ref_fn)
    var refx: hts.Fai
    assert hts.open(refx, ref_fn)
    assert refx.len() == 1
    let reference_dna = refx.get(refx[0])
    readaln(aln_fn, reference_dna)
    echo "bye"

when isMainModule:
    main()
