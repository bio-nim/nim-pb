# vim: sw=4 ts=4 sts=4 tw=0 et:

## RaptorDB v0.2 specification.
## Each line is a separate entity. Order of lines does not matter. Lines start with a keyword character,
## depending on that encode:
## - Version: "V\t<float-version_num>"
## - Files: "F\t<int64-file_id>\t<str-file_path>\t<str-file_format>"
## - Sequences: "S\t<int64-seq_id>\t<str-header>\t<int64-seq_len>\t<int64-file_id>\t<int64-start_offset_in_file>\t<int64-data_length>"
## - Blocks: "B\t<int64-block_id>\t<int64-seq_id_start>\t<int64-seq_id_end>\t<int64-num_bases_in_block>"

from sets import nil
from strformat import fmt
from strutils import nil
from ./util import nil

# If we wanted to parse into memory, we could use these.
type
    FileRecord = tuple
        file_id: int64
        file_path: string
        file_format: string
    SequenceRecord = tuple
        seq_id: int64
        header: string
        seq_len: int64
        file_id: int64
        start_offset_in_file: int64
        data_length: int64
    BlockRecord = tuple
        block_id: int64
        seq_id_start: int64
        seq_id_end: int64
    Db = object
        version: string
        files: seq[FileRecord]
        seqs: seq[SequenceRecord]
        blocks: seq[BlockRecord]

#proc toSequenceRecord(fields: seq[string]): SequenceRecord =
#    # should use strutils.parseInt()
#    result = (fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6])

#proc writeSequence(sout: File, fields: seq[string]) =
#    sout.write(
#        'S', '\t',
#        fields[1], '\t',
#        fields[2], '\t',
#        fields[3], '\t',
#        fields[4], '\t',
#        fields[5], '\t',
#        fields[6], '\n',

type
    SeqLineWriter = object
        num_seq_lines: int64
        sout: File
        num_bases_written: int64
        block_size: int64
        block_start_bases: int64
        block_start_id: int64
        num_blocks: int64

proc initSeqLineWriter(sout: File, block_size_MB: int): SeqLineWriter =
    let block_size = block_size_MB * 1024 * 1024
    return SeqLineWriter(num_seq_lines: 0, sout: sout, num_bases_written: 0, block_size: block_size,
        block_start_bases: 0, block_start_id: 0, num_blocks: 0)
proc write_block(w: var SeqLineWriter) =
    let curr_bases_in_block = w.num_bases_written - w.block_start_bases;
    let line = fmt("B\t{w.num_blocks}\t{w.block_start_id}\t{w.num_seq_lines}\t{curr_bases_in_block}\n")
    w.sout.write(line)
    inc w.num_blocks
    w.block_start_id = w.num_seq_lines
    w.block_start_bases = w.num_bases_written
proc write(w: var SeqLineWriter, split_line: seq[string]) =
    # This function updates the sequence IDs depending on the number
    # of already output sequence lines.
    var temp = split_line
    temp[1] = $w.num_seq_lines
    w.sout.write(strutils.join(temp, "\t"))
    w.sout.write('\n')
    inc w.num_seq_lines
    let nbases: int64 = strutils.parseBiggestInt(split_line[3])
    w.num_bases_written += nbases

    # Write the block if we reached the limit.
    let curr_bases_in_block = w.num_bases_written - w.block_start_bases;
    if curr_bases_in_block >= w.block_size:
        w.write_block()
proc close(w: var SeqLineWriter) =
    w.block_size = 1 # to force the last block if anything has been written
    w.write_block()

# Instead, we will stream-filter.
#
proc stream(sin, sout: File, blacklist: sets.HashSet[string]) =
    let block_size_MB = 1
    var writer = initSeqLineWriter(sout, block_size_MB)
    defer: writer.close()

    for line in lines(sin):
        if len(line) == 0:
            continue
        let keyletter = line[0]
        case keyletter
        of 'V', 'F':
            # Version and File list information.
            sout.write(line)
            sout.write('\n')
            continue
        of 'B':
            # We need to re-block because filtering will
            # reduce the size.
            continue
        of 'S':
            var fields = strutils.splitWhitespace(line)
            assert len(fields) == 7
            if not sets.contains(blacklist, fields[2]):
                writer.write(fields)
                #sout.write(line)
                #sout.write('\n')
            continue
        else:
            let msg = fmt("Unknown line '{line}'")
            #util.raiseEx(msg)
            util.log(msg)
            continue


proc filter*(blacklist_fn: string="") =
    ## Read/write raptor-db to/from stdin/stdout.
    ## Exclude zmws in blacklist.
    util.log("filter sans ", blacklist_fn)
    var blacklist = sets.initHashSet[string]()
    if len(blacklist_fn) > 0:
        var sin: File
        if not open(sin, blacklist_fn, fmRead):
            util.raiseEx(fmt("Failed to open blacklist file '{blacklist_fn}'"))
        defer: close(sin)
        for zmw in lines(sin):
            util.log(fmt("Skipping ({zmw})"))
            if sets.contains(blacklist, zmw):
                util.raiseEx(fmt("Found a repeat in blacklist '{blacklist_fn}': '{zmw}'\n Something is wrong!"))
            sets.incl(blacklist, zmw)
    stream(stdin, stdout, blacklist)
