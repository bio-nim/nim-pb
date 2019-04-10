# vim: sw=4 ts=4 sts=4 tw=0 et:
import deques
import tables
from algorithm import sort
from hashes import nil
from strutils import format

type
    Dna* = string             # someday, this might be an array
    Bin* = uint64             # compact bitvector of DNA
    ##  In bitvector, A is 0, C is 1, G is two, and T is 3.

    ##  kmer - a uint64 supporting a maximum of 32 DNA bases.
    ##  pos  - position along the sequence
    ##  strand - if true, reverse kmers
    seed_t* = object
        kmer*: Bin
        pos*: uint32
        strand*: bool

    ##  a & b are two seed_t's designed for matching in the hash lookup
    seed_pair_t* = object
        a*: seed_t
        b*: seed_t

    Hash* = int

    ##  seeds - a pointer to the kmers
    ##  n  - the number of kmers in the database (h)
    pot_t* = ref object
        word_size*: uint8     # <=32
        n*: uint32            # same as len(seeds)
        seeds*: seq[seed_t]
        ht*: ref tables.Table[Bin, int]
        searchable*: bool

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


##  @return uninitialized
#
proc newDna(size: int): Dna =
    return newString(size)

# hashes for sets and tables

proc hash*(s: kmers.seed_t): hashes.Hash =
   #hashes.hash(s.pos) + hashes.hash(s.kmer shl 8) + hashes.hash(s.strand)
   hashes.hash([s.pos.uint64, s.kmer.uint64, s.strand.uint64])

proc hash*(p: kmers.seed_pair_t): hashes.Hash =
   hashes.hash([hash(p.a), hash(p.b)])

# convenience for C coders

template `<<`(a, b: uint64): uint64 =
    a shl b

template `>>`(a, b: uint64): uint64 =
    a shr b


##  Converts a char * into a set of seed_t objects.
##  @param  sq  - sequence
##  @param  k - kmer size (<=32)
##  @return pot
#
proc dna_to_kmers*(sq: Dna; k: int): pot_t =
    if sq.len == 0 or k > 32:
        return nil
    new(result)

    var
        shift1: uint64 = 2'u64 * (k - 1).uint64
        mask: uint64 = (1'u64 << (2 * k).uint64) - 1
    #echo format("shift1=$# mask=$#", shift1, mask)

    var forward_kmer: seed_t
    var reverse_kmer: seed_t

    forward_kmer.kmer = 0
    forward_kmer.pos = 0
    reverse_kmer.kmer = 0
    reverse_kmer.pos = 0
    forward_kmer.strand = false
    reverse_kmer.strand = true

    var kmer_stack = deques.initDeque[seed_t](128)

    ##  lk is the length of the kmers being built on the fly. The variable n is the total number of
    var
        i: int
        lk: int
        n: int
    i = 0
    lk = 0
    n = 0

    while i < sq.len():
        let ch = cast[uint8](sq[i])
        let c = seq_nt4_table[ch].uint64
        if c < 4:
            forward_kmer.kmer = (forward_kmer.kmer << 2 or c) and mask
            reverse_kmer.kmer = (reverse_kmer.kmer >> 2) or (3'u64 xor c) << shift1
            #echo format("[$#]=$# $#==$#($# $#) f:$# r:$#",
            #    i, sq[i], ch, c, (3'u8 xor c), (3'u8 xor c).uint64 shl shift1, forward_kmer.kmer, reverse_kmer.kmer)
            inc(lk)
        else:
            ##  advance the window beyond the unknown character
            lk = 0
            inc(i, k)
            inc(forward_kmer.pos, k)
            forward_kmer.kmer = 0
            inc(reverse_kmer.pos, k)
            reverse_kmer.kmer = 0

        if lk >= k:
            inc(n, 2)
            deques.addLast(kmer_stack, forward_kmer)
            deques.addLast(kmer_stack, reverse_kmer)
            inc(forward_kmer.pos, 1)
            inc(reverse_kmer.pos, 1)
        inc(i)

    var kmers: pot_t
    new(kmers)
    kmers.seeds = newSeq[seed_t](n)
    kmers.n = n.uint32
    kmers.word_size = k.uint8
    kmers.searchable = false
    i = 0

    while i < n:
        kmers.seeds[i] = kmer_stack.popLast()
        #kmers.seeds[i] = kmer_stack[i] # would also be fine
        #echo format("[$#]->$#", i, kmers.seeds[i].kmer)
        inc(i)

    return kmers

##  A function to convert the binary DNA back into character
##  @param kmer   up to 32 2-bit bases
##  @param k      kmer length
##  @param strand If true, start at kth bit and go backwards.
##
##  Zero is A, one is C, G is two, and T is 3
#
proc bin_to_dna*(kmer: Bin; k: uint8; strand: bool): Dna =
    var lookup: array[4, char] = ['A', 'C', 'G', 'T']
    var mask: uint64 = 3
    var i: uint8 = 0
    var tmp: uint64 = 0
    var offset: uint64 = 0

    var dna = newDna(k.int)
    i = 0
    while i < k:
        tmp = kmer
        offset = if not strand: (k - i - 1) * 2 else: (i * 2)
        tmp = tmp >> offset
        dna[i] = lookup[mask and tmp]
        inc(i)

    return dna


##  Prints the pot structure to STDOUT
##  @param pot a ref to the pot
#
proc print_pot*(pot: pot_t) =
    var i: uint32 = 0

    i = 0
    while i < pot.n:
        var dna = bin_to_dna(pot.seeds[i].kmer, pot.word_size,
                         pot.seeds[i].strand)
        echo format("pos:$# strand:$# seq:$# bin:$#",
            pot.seeds[i].pos, pot.seeds[i].strand, dna, pot.seeds[i].kmer)
        inc(i, 1)

proc cmp_seeds(a, b: seed_t): int =
    let c = a.kmer
    let d = b.kmer

    if c < d:
        return -1

    if c == d:
        if a.pos < b.pos:
            return -1
        else:
            return 0

    return 1

##  A function that sorts and loads the kmers into a hash table and sets the sort
##  flag.
##  @param  pot_t * - a pointer to a pot
##  @return 0 if a new table was created
#
proc make_searchable*(kms: pot_t): int {.discardable.} =
    if kms.searchable:
        return 1
    kms.seeds.sort(cmp_seeds)
    kms.ht = newTable[Bin, int]()

    var i: int = 0
    while i < kms.seeds.len():
        let key = kms.seeds[i].kmer
        if kms.ht.hasKeyOrPut(key, i):
            echo format("WARNING: Duplicate seed $# @$#, not re-adding @$#", key, i, kms.ht[key], i)
        inc(i)

    kms.searchable = true
    return 0

##  A function that simply checks for the presence or absense of a kmer in a
##  pot regaurdless of the position
##  @param  pot_t * - a pointer to a pot
##  @return false if kmer doesn't exist or pot is not searchable
#
proc haskmer*(target: pot_t; query: Bin): bool =
 if not target.searchable:
  return false
 if target.ht.hasKey(query):
  return true
 return false

proc search*(target: pot_t; query: pot_t): deques.Deque[seed_pair_t] =
    discard make_searchable(target)
    echo format("Searching through $# kmers", query.n)
    var hit_stack = deques.initDeque[seed_pair_t](128)
    var hit: seed_pair_t
    var hit_index: int

    var i: int = 0
    #echo format("target.ht=$#", target.ht)
    #echo format("query.ht=$#", query.ht)
    while i < query.seeds.len():
        let key = query.seeds[i].kmer
        if key in target.ht:
            hit_index = target.ht[key]
            #echo format("For $# ($#), ql=$# tl=$#, hit_index=$#", i, key, query.seeds.len(), target.seeds.len(), hit_index)
            while (hit_index < target.seeds.len() and key == target.seeds[hit_index].kmer):
                #echo format("--For $# ($#), ql=$# tl=$#, hit_index=$#", i, key, query.seeds.len(), target.seeds.len(), hit_index)
                hit.a = query.seeds[i]
                hit.b = target.seeds[hit_index]
                deques.addLast(hit_stack, hit)
                inc(hit_index, 1)
        inc(i)

    return hit_stack
