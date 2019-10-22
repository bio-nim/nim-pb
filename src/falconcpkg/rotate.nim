from ./util import nil
from strformat import nil
import hts
import gc
import strutils
import algorithm
import random
import sets
import std/wordwrap
import logging
import os


var logger* = newConsoleLogger(fmtStr = "[$time] - $levelname: ")
addHandler(logger)

type

    dat* = tuple[pos: int, skew: float, accum: float]

    skewDat* = object
        data*: seq[dat]
        mini*: int
        minv*: float

proc checkEmptyFile(fin: string): bool =
    var finfo = getFileInfo(fin)
    if finfo.size == 0:
        return true
    return false

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

proc randomize*(input: string, output: string, seed: int64 = 0) =
    ##randomly rotates left rotates the sequence and writes to the output file.
    if seed != 0:
        random.randomize(seed)
    else:
        random.randomize()
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
        var ran = random.rand(chrom_len)
        discard algorithm.rotateLeft(full_sequence, ran)
        echo "[INFO] randomizing: ", chrom_name
        outfh.write(">", chrom_name, " randomly_shifted_by_bp:-",
                 ran, "/", chrom_len, "\n")
        outfh.write(wrapWords(full_sequence), "\n")

type
    WhiteList = ref object
        specific: HashSet[string]

proc loadWhiteList*(fin: string): WhiteList =
    ## If filename is blank, return nil, to indicate that everything is whitelisted.
    var whitelist: WhiteList
    if fin == "":
        return whitelist # nil
    new(whitelist)
    whitelist.specific = initHashSet[string]()
    echo "whitelist:"
    for l in fin.lines:
        echo " whitelisted: ", l
        whitelist.specific.incl(l)
    return whitelist

proc whitelisted*(whitelist: WhiteList, chrom_name: string): bool =
    if isNil(whitelist):
        return true
    return whitelist.specific.contains(chrom_name)

type
    SimpleFastaWriter = ref object
        fout: File
        width: int
    SimpleFastqWriter = ref object
        fout: File
        width: int

proc newSimpleFastaWriter(fn: string, width: int=80): SimpleFastaWriter =
    new(result)
    result.fout = open(fn, fmWrite)
    result.width = width

proc write(obj: SimpleFastaWriter, full_sequence: string, chrom_name: string, shift: int) =
    let chrom_len = len(full_sequence)
    obj.fout.write(">", chrom_name, " shifted_by_bp:-",
                shift, "/", chrom_len, "\n")
    obj.fout.write(wrapWords(full_sequence, obj.width), "\n")

proc close(obj: SimpleFastaWriter) =
    close(obj.fout)

proc newSimpleFastqWriter(fn: string, width: int=80): SimpleFastqWriter =
    new(result)
    result.fout = open(fn, fmWrite)
    result.width = width

proc write(obj: SimpleFastqWriter, full_sequence, full_qvs: string, chrom_name: string, shift: int) =
    let chrom_len = len(full_sequence)
    obj.fout.write("@", chrom_name, " shifted_by_bp:-",
                shift, "/", chrom_len, "\n")
    obj.fout.write(wrapWords(full_sequence, obj.width), "\n")
    obj.fout.write("+\n")
    obj.fout.write(wrapWords(full_qvs, obj.width), "\n")

proc close(obj: SimpleFastqWriter) =
    close(obj.fout)

iterator FaiReader(fn: string, full_sequence: var string): string {.closure.} =
    # Yield chrom_name; modify full_sequence
    var fai: Fai
    if not fai.open(fn):
        let msg = strformat.fmt("Problem loading fasta file '{fn}'")
        util.raiseEx(msg)
    for i in 0..<fai.len:
        var chrom_name = fai[i]
        var chrom_len = fai.chrom_len(chrom_name)
        full_sequence = fai.get(chrom_name) # modify input var
        yield chrom_name

iterator FastqReader(fn: string, full_sequence, full_qvs: var string): string {.closure.} =
    # Yield chrom_name; modify full_sequence, full_qvs
    for i in 0..1:
        full_sequence = "GATTACA"
        full_qvs = "123"
        let chrom_name = "mychrom"
        yield chrom_name

proc reorientFASTA(fin: string, fon: string, wl: string, w: int, s: int,
        print: bool) =
    var writer = newSimpleFastaWriter(fon)
    defer: writer.close()

    if checkEmptyFile(fin):
        logger.log(lvlNotice, "Empty input, output will be empty.")
        return

    let whiteList = loadWhiteList(wl)

    var full_sequence: string # alg.rotate needs var, so Reader cannot just return it.

    for chrom_name in FaiReader(fin, full_sequence):
        var sdf = calcSkew(full_sequence, w, s)
        if not whitelisted(whiteList, chrom_name):
            sdf.data[sdf.mini].pos = 0
        if print:
            printSkew(chrom_name, sdf)
        #echo "#windows:", sdf.data.len, " pivot index:", sdf.mini
        #echo "pivot pos: ", sdf.data[sdf.mini].pos, " / ", len(full_sequence),
        # " seq: ", chrom_name

        discard algorithm.rotateLeft(full_sequence, sdf.data[sdf.mini].pos)

        writer.write(full_sequence, chrom_name, sdf.data[sdf.mini].pos)

proc reorientFASTQ(fin: string, fon: string, wl: string, w: int, s: int,
        print: bool) =
    var writer = newSimpleFastqWriter(fon)
    defer: writer.close()

    if checkEmptyFile(fin):
        logger.log(lvlNotice, "Empty input, output will be empty.")
        return

    let whiteList = loadWhiteList(wl)

    var full_sequence, full_qvs: string # alg.rotate needs var, so Reader cannot just return it.

    for chrom_name in FastqReader(fin, full_sequence, full_qvs):
        var sdf = calcSkew(full_sequence, w, s)
        if not whitelisted(whiteList, chrom_name):
            sdf.data[sdf.mini].pos = 0
        if print:
            printSkew(chrom_name, sdf)
        #echo "#windows:", sdf.data.len, " pivot index:", sdf.mini
        #echo "pivot pos: ", sdf.data[sdf.mini].pos, " / ", len(full_sequence),
        # " seq: ", chrom_name

        discard algorithm.rotateLeft(full_sequence, sdf.data[sdf.mini].pos)

        writer.write(full_sequence, full_qvs, chrom_name, sdf.data[sdf.mini].pos)


proc main*(input_fn: string, output_fn: string, wl_fn = "", window = 500,
        step = 200, print = false) =
    ##reorients circular sequences based on gc-skew distribution and writes to output.
    if input_fn == "" or output_fn == "":
        logger.log(lvlFatal, "Missing input or output required options.")
        quit 1
    logger.log(lvlInfo, "Reorienting.")
    reorientFASTA(input_fn, output_fn, wl_fn, window, step, print)

when isMainModule:
    main()
