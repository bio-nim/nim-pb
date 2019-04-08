# vim: sw=4 ts=4 sts=4 tw=0 et:
import unittest
import pbpkg/kmers
import deques
import sequtils
import sets

test "bin_to_dna":
    check bin_to_dna(0, 1, false) == "A"
    check bin_to_dna(1, 1, false) == "C"
    check bin_to_dna(2, 1, false) == "G"
    check bin_to_dna(3, 1, false) == "T"

    check bin_to_dna(0b00011011, 4, false) == "ACGT"
    check bin_to_dna(0b00011011, 4, true)  == "TGCA"

test "dna_to_kmers":
    check dna_to_kmers("AAAA", 2).seeds.len() == 6

test "sorted_kmers":
    let
        sq = "ATCGGCTACTATT"
        expected = [
            "AGCCGATGATAA",
            "TAGCCGATGATA",
            "ATCGGCTACTAT",
            "TCGGCTACTATT",
        ]
        k = 12
        kms: pot_t = dna_to_kmers(sq, k)
    discard make_searchable(kms)  # sort
    let got = sequtils.mapIt(kms.seeds, bin_to_dna(it.kmer, k.uint8, it.strand))
    check got == expected

test "search":
    let
        sq = "ATCGGCTACTATT"
        k = 12
        kms: pot_t = dna_to_kmers(sq, k)
        qms: pot_t = dna_to_kmers(sq, k)
    check make_searchable(kms) == 0
    let hits = search(kms, qms)
    check hits.len() == 4
    check sets.toHashSet(seqUtils.toSeq(hits)).len() == 4  # 4 unique items
    #check sets.len(sets.toHashSet(seqUtils.toSeq(deques.items(hits)))) == 4  # same as above
