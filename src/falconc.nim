# vim: sw=4 ts=4 sts=4 tw=0 et:
#import falconcpkg/zev
from falconcpkg/align import nil
from falconcpkg/raptor_db import nil
from falconcpkg/rotate import nil
from falconcpkg/phasr import nil
from falconcpkg/overlapFilter import nil
from falconcpkg/pbcromwell import nil
from falconcpkg/pbreports import nil
from falconcpkg/stats import nil
from falconcpkg/falcon/rr_hctg_track import nil
import cligen

const
    cNimbleData = staticRead("../falconc.nimble")
    cGitSha1 {.strdefine.} = staticExec("git describe --always --tags HEAD")
    cToolVersion = cNimbleData.fromNimble("version")

proc get_version(): auto =
    cToolVersion & "+git." & cGitSha1

proc version() =
    echo cToolVersion & "+git." & cGitSha1
proc dataset(extras: seq[string]) =
    echo "falconc dataset"
proc kmers(int_dummy: int = 42, string_dummy: string = "hello") =
    echo "falconc kmers"
proc utils(extras: seq[string], float_req: float) =
    echo "falconc utils ..."
#proc zev() =
#    echo "starting"
#    zev.foo()
#    echo "finished"

when isMainModule:
    # Show version at start-up, for now.
    stderr.writeLine("  version=", get_version(), ", nim-version=",
            system.NimVersion)

    dispatchMulti(
        [version],
        [dataset, short = {}, help = {}],
        [kmers, short = {"int_dummy": 'd'}, help = {}],
        [utils, short = {}, help = {"float_req": "special help message"}],
        [align.bam_filter_ipa, cmdName = "bam-filter-ipa"],
        [align.bam_filter_ipa, cmdName = "align-filter",
                doc = "alias for bam-filter-ipa"],
        [align.bam_tags_enrich, cmdName = "bam-tags-enrich",
         help = {
          "input-fn": "Bam or Sam filename (based on its extension), or '-'",
          "output-fn": "Bam or Sam filename (based on its extension)",
            }
        ],
        [align.bam_filter_clipped, cmdName = "bam-filter-clipped",
         help = {
          "input-fn": "Bam or Sam filename (based on its extension), or '-'",
          "output-fn": "Bam or Sam filename (based on its extension)",
          "max-clipping": "Maximum clipping on left or right of query, in basepairs",
          "end-margin": "Maximum margin on contig ends, in basepairs",
          "Flags-exclude": "Integer (0x ok) of flags to exclude, independent of other filters",
          "tags-enrich": "Also enrich the tags. (See bam-tags-enrich.)",
          "verbose": "Show each skipped alignment, and a count.",
            }
        ],
        [raptor_db.filter, cmdName = "raptor-db-filter"],
        [raptor_db.calc_length_cutoff, cmdName = "raptor-db-calc-length-cutoff",
            help = {
        "fail-low-cov": "Exit non-zero if the input coverage was insufficient to satisfy the requirement.",
        "alarms-fn": "Optional JSON file to write SMRT Link alarms. (typically used w/ -f)",
            },
        ],
        [rotate.main, cmdName = "circ-orient",
            help = {
        "input-fn": "fasta (or fastq) file of circular sequences",
        "output-fn": "fasta (or fastq) file output",
        "wl-fn": "white list of sequences to rotate, one per line, no spaces, no trailing spaces",
        "window": "window size to caculate gc-skew",
        "step": "window step",
        "print": "whether to print skew data to files ('SEQ.gc_skew.txt', one per SEQ)"
            },
        ],
        [rotate.randomize, cmdName = "circ-randomize",
            help = {
        "input": "fasta file of circular sequences",
        "output": "fasta file output",
        "seed": "set seed, if non-zero",
            },
        ],
        [phasr.main, cmdName = "phasr",
          help = {
            "aln-fn": "BAM alignment, sorted on 'coordinate'",
            "ref-fn": "FASTA reference",
            "out-fn": "Output file name",
            "iterations": "Number of phasing iterations per read",
            "kmersize": "Kmer size",
            "delta": "Frequency cut",
            },
        ],
        [overlapFilter.runMergeBlacklists, cmdName = "m4filt-merge-blacklist",
        ],
        [overlapFilter.runStage1, cmdName = "m4filt-stage1",
        ],
        [overlapFilter.runStage2, cmdName = "m4filt-stage2",
        ],
        [overlapFilter.runDumpBlacklist, cmdName = "m4filt-dump-blacklist",
        ],
        [overlapFilter.m4filtRunner, cmdName = "m4filt",
         help = {
          "idt-stage1": "Stage one percent identity filter, formatted as percentage, overlaps < %idt are skipped",
          "idt-stage2": "Stage two percent identify filter",
          "min-len": "Minimum read length, reads shorter than minLen will be discarded",
          "min-cov": "Minimum number of overlaps on either side of a read",
          "max-cov": "Maximum number of overlaps on either side of a read",
          "max-diff": "Reads are skipped is abs(5p-3p) overlap counts > maxDiff",
          "bestn": "Keep N best overlaps at 5prime AND 3prime of a read",
          "min-depth": "Depths lower than minDepth are considered gaps",
          "gap-filt": "Run depth filter, takes a little more time",
          "n-proc": "Number of processes to run locally",
          "filter-log-fn": "Write read filter stats to this file",
          "out-fn": "Final m4 overlap file",
            }
        ],
        [overlapFilter.falconRunner, cmdName = "m4filt-falconRunner",
         help = {
          "las-json-fn": "List of las files from falcon, e.g ../1-preads_ovl/las-merge-combine/las_fofn.json",
          "idt-stage1": "Stage one percent identity filter, formatted as percentage, overlaps < %idt are skipped",
          "idt-stage2": "Stage two percent identify filter",
          "min-len": "Minimum read length, reads shorter than minLen will be discarded",
          "min-cov": "Minimum number of overlaps on either side of a read",
          "max-cov": "Maximum number of overlaps on either side of a read",
          "max-diff": "Reads are skipped is abs(5p-3p) overlap counts > maxDiff",
          "bestn": "Keep N best overlaps at 5prime AND 3prime of a read",
          "min-depth": "Depths lower than minDepth are considered gaps",
          "gap-filt": "Run depth filter, takes a little more time",
          "n-proc": "Number of processes to run locally",
          "filter-log-fn": "Write read filter stats to this file",
          "debug-log-fn": "Write stderr to this file",
          "out-fn": "Final m4 overlap file",
            }
        ],
        [overlapFilter.ipaRunner, cmdName = "m4filt-ipaRunner",
         help = {
          "ovls-fofn-fn": "List of m4 files from ipa/raptor",
          "idt-stage1": "Stage one percent identity filter, formatted as percentage, overlaps < %idt are skipped",
          "idt-stage2": "Stage two percent identify filter",
          "min-len": "Minimum read length, reads shorter than minLen will be discarded",
          "min-cov": "Minimum number of overlaps on either side of a read",
          "max-cov": "Maximum number of overlaps on either side of a read",
          "max-diff": "Reads are skipped is abs(5p-3p) overlap counts > maxDiff",
          "bestn": "Keep N best overlaps at 5prime AND 3prime of a read",
          "min-depth": "Depths lower than minDepth are considered gaps",
          "gap-filt": "Run depth filter, takes a little more time",
          "n-proc": "Number of processes to run locally",
          "filter-log-fn": "Write read filter stats to this file",
          "outputFn": "Final m4 overlap file",
            }
        ],
        [pbcromwell.remove_las, cmdName = "pbcromwell-rm-las",
         short = {"dry_run": 'n'},
        ],
        [pbreports.circ, cmdName = "pbreports-circ",
         help = {
          "fasta-fn": "FASTA filename, preferably indexed (with .fai)",
          "circ-fn": "Text list of circular contigs (newline-delimited)",
            }
        ],
        [stats.assembly, cmdName = "stats-assembly",
         help = {
          "fasta-fn": "FASTA filename, preferably indexed (with .fai)",
            }
        ],
        [rr_hctg_track.run_stage1, cmdName = "rr-hctg-track1",
        ],
        [rr_hctg_track.run_stage2, cmdName = "rr-hctg-track2",
        ],
    )
