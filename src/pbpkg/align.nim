# vim: sw=4 ts=4 sts=4 tw=0 et:
import hts
from algorithm import nil
from sequtils import nil
from strutils import format
from tables import contains, `[]`, `[]=`
from ./util import log

proc logRec(record: Record) =
    # I think len == stop-start+1, but I need to verify. ~cd
    var s: string
    hts.sequence(record, s)
    log(format("$# $# ($#) [$# .. $#] $# seqlen=$#", record.tid, record.chrom,
        record.qname, record.start, record.stop,
        ($record.cigar).substr(0, 32), s.len()))

type
    Params = object
        min_len: int
        min_frac: float

# Someday we might want to look at the CIGAR alignments.
# But for now, we call it a decent alignment if it is long enough.
proc update_counts(bam_fn: string, params: Params,
        everything: var tables.CountTable[string],
        exclusions: var tables.CountTable[string]) =
    var
        totals: tables.Table[string, int]
        qlengths: tables.Table[string, int]
        q: string # temp
        b: hts.Bam
    totals = tables.initTable[string, int]()
    qlengths = tables.initTable[string, int]()

    hts.open(b, bam_fn)
    defer: hts.close(b)
    for record in b:
        let key: string = record.qname
        tables.inc(everything, key)
        let reflen = record.stop - record.start + 1
        totals[key] = reflen + tables.getOrDefault(totals, key, 0)
        hts.sequence(record, q)
        qlengths[key] = len(q)
        if reflen >= params.min_len:
            tables.inc(exclusions, key)
            #log("exclude:" & key)
        #logRec(record)
    for key, qlen in tables.pairs(qlengths):
        let reflentotal = totals[key]
        let frac = float(reflentotal) / float(qlen)
        if frac >= params.min_frac:
            tables.inc(exclusions, key)
            #log("exclude:" & key)
        #log(format("$# $#/$#=$#", key, reflentotal, qlen, frac))

proc align_filter*(bams_fofn: string, min_len = 300000, min_frac = 0.70) =
    ## Print subreads which have decent alignments in any of the bam inputs.
    var params = Params(min_len: min_len, min_frac: min_frac)
    log(params)
    var everything = tables.initCountTable[string]()
    var exclusions = tables.initCountTable[string]()
    for fn in lines(bams_fofn):
        log("Processing ", fn)
        update_counts(fn, params, everything, exclusions)
    log("everything.len=", tables.len(everything))
    log("exclusions.len=", tables.len(exclusions))
    var sorted_exclusions: seq[string] = sequtils.toSeq(tables.keys(exclusions))
    algorithm.sort(sorted_exclusions)
    for key in sorted_exclusions:
        #log("Count of '", key, "': ", everything[key])
        echo key
