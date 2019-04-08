# vim: sw=4 ts=4 sts=4 tw=0 et:
#from strformat import fmt
import deques
from os import nil
from strutils import format
from pbpkg/kmers import nil

proc main*(args: seq[string]): int =
    var sq = "ATCGGCTACTATT"

    echo format("Starting seq: $#", sq)

    var ans_lookups = [
        "AGCCGATGATAA",
        "TAGCCGATGATA",
        "ATCGGCTACTAT",
        "TCGGCTACTATT",
    ]
    var kms: kmers.pot_t = kmers.dna_to_kmers(sq, 12)
    var qms: kmers.pot_t = kmers.dna_to_kmers(sq, 12)

    echo "kms"
    kmers.print_pot(kms)
    echo ""
    echo "qms"
    kmers.print_pot(qms)

    discard kmers.make_searchable(kms)  # also sorts

    var final_res: int = 0
    var i: uint32 = 0

    while i < kms.n:
        let tmp = kmers.bin_to_dna(kms.seeds[i].kmer, kms.word_size,
                                   kms.seeds[i].strand)
        let res = cmp(tmp, ans_lookups[i])

        echo format("[$#:DNA->BIT->DNA:TEST] kmer: $# pos:$# expecting:$# observed:$# [$#]",
            system.currentSourcePath(), kms.seeds[i].kmer, kms.seeds[i].pos, ans_lookups[i], tmp,
                if res != 0: "FAIL" else: "PASS")

        final_res = final_res or res

        inc(i)

    var hits = kmers.search(kms, qms)

    try:
        while true:
            let pair = deques.popLast(hits)
            echo format("qb:$# tb:$# qs:$# ts:$# $# $#",
                pair.a.pos, pair.b.pos,
                pair.a.strand, pair.b.strand,
                kmers.bin_to_dna(pair.a.kmer, kms.word_size, pair.a.strand),
                kmers.bin_to_dna(pair.b.kmer, kms.word_size, pair.b.strand))
    except:
        discard

    return final_res

when isMainModule:
    quit main(os.commandLineParams())
