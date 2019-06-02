import unittest
import "../src/pbpkg/gc"

test "okay_count":
    var dna: string = "ATCG"
    var c = gc.countBases(dna)
    echo(c.okayBases)
    check c.okayBases == 4
    dna = "atcg"
    c = gc.countBases(dna)
    check c.okayBases == 4

test "not_okay_count":
    var dna: string = "XXXX"
    var c = gc.countBases(dna)
    check c.okayBases == 0
    dna = "ancg"
    c = gc.countBases(dna)
    check c.okayBases == 3

test "A_pos":
    var dna: string = "aAGC"
    var c = gc.countBases(dna)
    check gc.A(c) == 2

test "A_neg":
    var dna: string = "TTGC"
    var c = gc.countBases(dna)
    check gc.A(c) == 0

test "CG_0":
    var dna: string = "TTGC"
    var c = gc.countBases(dna)
    check gc.GC(c) == 2

test "GC_1":
    var dna: string = "TTgC"
    var c = gc.countBases(dna)
    check gc.GC(c) == 2

test "revComp_0":
    var dna: string = "TTgC"
    var rev = gc.revComp(dna)
    check rev == "GcAA"

test "revComp_1":
    var dna: string = "TTgCNNaa"
    var rev = gc.revComp(dna)
    check rev == "ttNNGcAA"
