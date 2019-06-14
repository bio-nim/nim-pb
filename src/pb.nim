# vim: sw=4 ts=4 sts=4 tw=0 et:
#import pbpkg/zev
from pbpkg/align import nil
from pbpkg/raptor_db import nil
from pbpkg/rotate import nil
from pbpkg/phasr import nil

proc dataset(extras: seq[string]) =
    echo "pb dataset"
proc kmers(int_dummy: int = 42, string_dummy: string = "hello") =
    echo "pb kmers"
proc utils(extras: seq[string], float_req: float) =
    echo "pb utils ..."
#proc zev() =
#    echo "starting"
#    zev.foo()
#    echo "finished"

when isMainModule:
    import cligen
    dispatchMulti(
        # [zev],
        [dataset, short = {}, help = {}],
        [kmers, short = {"int_dummy": 'd'}, help = {}],
        [utils, short = {}, help = {"float_req": "special help message"}],
        [align.align_filter, cmdName="align-filter"],
        [raptor_db.filter, cmdName="raptor-db-filter"],
        [rotate.main, cmdName="circ-orient",
            help = {
       "input": "fasta file of circular sequences",
       "output": "fasta file output",
       "wl": "white list of sequences to rotate, one per line, no spaces, no trailing spaces",
       "window": "window size to caculate gc-skew",
       "step": "window step",
       "print": "print skew data to files, one per sequence"
            },
        ],
        [rotate.randomize, cmdName="circ-randomize",
            help = {
       "input": "fasta file of circular sequences",
       "output": "fasta file output",
       "seed": "set seed, if non-zero",
            },
        ],
        [phasr.main, cmdName="phasr",
          help = {
            "aln_fn": "BAM alignment, sorted on 'coordinate'",
            "ref_fn": "FASTA reference",
            "out_fn": "Output file name",
            "iterations" : "Number of phasing iterations per read",
            "kmersize"   : "Kmer size",
            "delta"      : "Frequency cut",
          },
        ],
    )
