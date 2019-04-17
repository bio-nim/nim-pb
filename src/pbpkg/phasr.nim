# vim: sw=4 ts=4 sts=4 tw=0 et:
from sequtils import nil
from strutils import format
import deques
import hts
import os
from ./util import raiseEx
import ./kmers

proc foo*() =
    echo "foo"
proc showRec(record: Record) =
    echo format("$# $# ($#) [$# .. $#] $#", record.tid, record.chrom, record.qname, record.start, record.stop,
        ($record.cigar).substr(0, 32))
type
    ProcessedRecord = object
        # These are both refs already.
        rec: Record
        kmers: pot_t

proc processRecord(record: Record, klen: int): ProcessedRecord =
    var rseq: string
    discard hts.sequence(record, rseq)
    let kmers: pot_t = dna_to_kmers(rseq, klen)
    return ProcessedRecord(rec: record, kmers: kmers)

type
    Pileup = object
        # I'm not sure I've named this well. It's a list of records that overlap a given one,
        # all processed for kmers already.
        precs: seq[ProcessedRecord]
        current: int  # This is the index into precs of the record we are currently comparing with.
        # All other recs overlap this one.

iterator overlaps(b: hts.Bam, klen: int): Pileup =
    var current: Record
    var current_stack_index: int
    var stack = deques.initDeque[ProcessedRecord](64)
    for r in b:
        let record = hts.copy(r) # need a copy because iterator stores record on the Bam struct
        showRec(record)  # DEBUG

        stack.addLast(processRecord(record, klen))
        if stack.len() == 1:
            current_stack_index = 0
            current = record
        if record.start >= current.stop:
            yield Pileup(precs: sequtils.toSeq(stack), current: current_stack_index)  # YIELD

            # Switch current to next record.
            current_stack_index += 1
            current = stack[current_stack_index].rec
            let num_before = current_stack_index - 1
            for i in 0 .. num_before - 1:
                let ri = stack.peekFirst().rec
                #echo " Pop?", ri.stop, "<=?", current.start
                if ri.stop <= current.start:
                    discard stack.popFirst()
                    current_stack_index -= 1
                    assert current == stack[current_stack_index].rec

proc readaln*(bfn: string; fasta: string) =
    const klen = 21
    var b: hts.Bam

    hts.open(b, bfn, index=true)
    # We do not really need the index, but it proves that the Bam is sorted.
    # Actually, we expect it to be sorted with @SO=coordinate, so maybe we
    # should verify that. TODO(CD).
    # Note that the sort is wrong for circular genomes; the secondary key is
    # the query, not the coordinate. So we will need to do something tricky
    # for circular references. TODO(CD).
    echo "[INFO] reading bam"
    for ovlps in overlaps(b, klen=21):
        echo "len=", len(ovlps.precs)
        echo ovlps.current, " range=", ovlps.precs[0].rec.start, "..", ovlps.precs[^1].rec.stop
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
    if not hts.open(refx, ref_fn):
        raiseEx(format("Could not open '$#'", ref_fn))
    assert refx.len() == 1
    let reference_dna = refx.get(refx[0])
    readaln(aln_fn, reference_dna)
    echo "bye"

when isMainModule:
    main()
