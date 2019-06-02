# vim: sw=4 ts=4 sts=4 tw=0 et:
import hts
from algorithm import nil
from sequtils import nil
from strutils import format
from tables import contains, `[]`

proc log(words: varargs[string, `$`]) =
    for word in words:
        write(stderr, word)
    write(stderr, '\l')

proc logRec(record: Record) =
    # I think len == stop-start+1, but I need to verify. ~cd
    log(format("$# $# ($#) [$# .. $#] $#", record.tid, record.chrom, record.qname, record.start, record.stop,
        ($record.cigar).substr(0, 32)))

type
    Params = object
        min_len: int

# Someday we might want to look at the CIGAR alignments.
# But for now, we call it a decent alignment if it is long enough.
proc update_counts(bam_fn: string, params: Params,
        everything: var tables.CountTable[string],
        exclusions: var tables.CountTable[string]) =
    var b: hts.Bam
    hts.open(b, bam_fn)    
    defer: hts.close(b)
    for record in b:
        let key: string = record.qname
        tables.inc(everything, key)
        let len = record.stop - record.start + 1
        if len >= params.min_len:
            tables.inc(exclusions, key)
        else:
            logRec(record)

proc align_filter*(bams_fofn: string, min_len=1000) =
    ## Print subreads which have decent alignments in any of the bam inputs.
    var params = Params(min_len:min_len)
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
