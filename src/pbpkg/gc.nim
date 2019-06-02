import tables

var complement = {'A': 'T', 'a': 't', 'T': 'A', 'a': 't', 'G': 'C', 'g': 'c',
        'C': 'G', 'c': 'g'}.toTable

var seq_nt4_table: array[256, int] = [
        0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

type
    BaseCounts* = object
        okayBases*: int
        counts*: tables.CountTable[char]

proc revComp*(dna: string): string =
    var j: int = dna.len-1

    for i in 0..j:
        var c: char = dna[j - i]
        if c in complement:
            result &= complement[c]
        else:
            result &= c

proc countBases*(dna: string): BaseCounts =

    var b: BaseCounts
    b.okayBases = 0
    b.counts = initCountTable[char]()

    for c in dna:
        b.counts.inc(c)
        if seq_nt4_table[int(c)] < 4:
            b.okayBases += 1

    return b

proc A*(c: BaseCounts): int =
    if 'A' in c.counts:
        result += c.counts['A']
    if 'a' in c.counts:
        result += c.counts['a']

proc T*(c: BaseCounts): int =
    if 'T' in c.counts:
        result += c.counts['T']
    if 't' in c.counts:
        result += c.counts['t']

proc G*(c: BaseCounts): int =
    if 'G' in c.counts:
        result += c.counts['G']
    if 'g' in c.counts:
        result += c.counts['g']

proc C*(c: BaseCounts): int =
    if 'C' in c.counts:
        result += c.counts['C']
    if 'c' in c.counts:
        result += c.counts['c']

proc AT*(c: BaseCounts): int =
    return A(c) + T(c)

proc GC*(c: BaseCounts): int =
    return G(c) + C(c)

proc perGC*(c: BaseCounts): float =
    return GC(c) / c.okayBases

proc gcSkew*(c: BaseCounts): float =
    return (G(c) - C(c)) / GC(c)
