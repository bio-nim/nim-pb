import hts
import gc
import strutils
import algorithm
import random
import sets
import std/wordwrap

type

    dat* = tuple[pos: int, skew: float, accum: float]

    skewDat* = object
        data*: seq[dat]
        mini*: int
        minv*: float

proc calcSkew*(dna: string, win: int, step: int): skewDat =
    var sk: skewDat
    var accum: float = 0
    var i, j: int = 0;
    var min: float = 100000000

    while(i <= dna.len - win - 1):
        var c = gc.countBases(dna[i..i+win])
        var d: dat
        var skew = gc.gcSkew(c)

        accum += skew
        d.pos = i
        d.skew = skew
        d.accum = accum
        sk.data.add(d)

        if accum < min:
            sk.mini = j
            sk.minv = accum
            min = accum

        inc(j)
        i += step

    return sk

proc printSkew*(seqname: string, sk: skewDat) =
    var output = open(seqname & ".gc_skew.txt", fmWrite)
    for i in sk.data:
        output.write(seqname, " ", i.pos, " ", i.skew.formatFloat(ffDecimal,
                2), " ", i.accum.formatFloat(ffDecimal, 2), "\n")
    output.close()

proc randomize*(input: string, output: string) =
    ##randomly rotates left rotates the sequence and writes to the output file.
    var fai: Fai
    if not fai.open(input):
        quit "couldn't open fasta"

    var outfh = open(output, fmWrite)
    defer:
        close(outfh)
    for i in 0..<fai.len:
        let chrom_name = fai[i]
        let chrom_len = fai.chrom_len(chrom_name)
        var full_sequence = fai.get(chrom_name)
        var ran = rand(chrom_len)
        discard algorithm.rotateLeft(full_sequence, ran)
        echo "[INFO] randomizing: ", chrom_name
        outfh.write(">", chrom_name, " randomly_shifted_by_bp:-",
                 ran, "/", chrom_len, "\n")
        outfh.write(wrapWords(full_sequence), "\n")


proc loadWhiteList(fin: string): HashSet[string] =
    var whitelist = initHashSet[string]()
    if fin == "":
        return whitelist
    for l in fin.lines:
        echo l
        whitelist.incl(l)
    return whitelist

proc reorient(fin: string, fon: string, wl: string, w: int, s: int,
        print: bool) =
    var fai: Fai
    if not fai.open(fin):
        quit "couldn't open fasta"

    var whiteList = loadWhiteList(wl)

    var output = open(fon, fmWrite)
    defer:
        output.close()

    for i in 0..<fai.len:
        var chrom_name = fai[i]
        var chrom_len = fai.chrom_len(chrom_name)
        var full_sequence = fai.get(chrom_name)
        var sdf = calcSkew(full_sequence, w, s)
        if not (whiteList.contains(chrom_name)) and (whiteList.len > 0):
            sdf.data[sdf.mini].pos = 0
        if print:
            printSkew(chrom_name, sdf)
        echo "#windows:", sdf.data.len, " pivot index:", sdf.mini
        echo "pivot pos: ", sdf.data[sdf.mini].pos, " / ", chrom_len,
         " seq: ", chrom_name
        discard algorithm.rotateLeft(full_sequence, sdf.data[sdf.mini].pos)
        output.write(">", chrom_name, " shifted_by_bp:-",
                   sdf.data[sdf.mini].pos, "/", chrom_len, "\n")
        output.write(wrapWords(full_sequence), "\n")


proc main*(input: string, output: string, wl = "", window = 500, step = 200,
        print = false) =
    ##reorients circular sequences based on gc-skew distribution and writes to output.
    if input == "" or output == "":
        quit "\n\n [FATAL] missing input or output required options \n\n"
    reorient(input, output, wl, window, step, print)

when isMainModule:
    main()
