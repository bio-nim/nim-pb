# vim: sw=4 ts=4 sts=4 tw=0 et:
from pbpkg/kmers import hash # avoiding "*" imports
import unittest
import deques
import sequtils
import sets

test "bin_to_dna":
    check kmers.bin_to_dna(0, 1, false) == "A"
    check kmers.bin_to_dna(1, 1, false) == "C"
    check kmers.bin_to_dna(2, 1, false) == "G"
    check kmers.bin_to_dna(3, 1, false) == "T"

    check kmers.bin_to_dna(0b00011011, 4, false) == "ACGT"
    check kmers.bin_to_dna(0b00011011, 4, true)  == "TGCA"

test "dna_to_kmers":
    check kmers.dna_to_kmers("AAAA", 2).seeds.len() == 6

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
        kms = kmers.dna_to_kmers(sq, k)
    discard kmers.make_searchable(kms)  # sort
    let got = sequtils.mapIt(kms.seeds, kmers.bin_to_dna(it.kmer, k.uint8, it.strand))
    check got == expected

test "search":
    let
        sq = "ATCGGCTACTATT"
        k = 12
        kms = kmers.dna_to_kmers(sq, k)
        qms = kmers.dna_to_kmers(sq, k)
    check kmers.make_searchable(kms) == 0
    let hits = kmers.search(kms, qms)
    check hits.len() == 4
    check sets.toHashSet(seqUtils.toSeq(hits)).len() == 4  # 4 unique items
    #check sets.len(sets.toHashSet(seqUtils.toSeq(deques.items(hits)))) == 4  # same as above

test "complement":
 let
  sq = "ATCGGCTACTATT"
  k = 12
  kms = kmers.dna_to_kmers(sq, k)
  qms = kmers.dna_to_kmers(sq, k)
 discard kmers.make_searchable(qms)
 discard kmers.complement(kms, qms)
 check kmers.nkmers(kms) == 0
 check kmers.nkmers(qms) == 4
