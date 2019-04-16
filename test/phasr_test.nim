# vim: sw=4 ts=4 sts=4 tw=0 et:
from strutils import format
from pbpkg/util import raiseEx
import hts
import os
import pbpkg/phasr

const
    #ref_fn = "/pbi/dept/secondary/siv/testdata/hgap/synth5k/ref.fasta"
    #aln_fn = "/pbi/dept/secondary/siv/testdata/hgap/synth5k/aln.bam"
    #ref_fn = "/pbi/dept/secondary/siv/testdata/hgap/greg200k-sv2/ref1.fasta"
    #aln_fn = "/pbi/dept/secondary/siv/testdata/hgap/greg200k-sv2/align.bam"
    ref_fn = "/pbi/dept/secondary/siv/testdata/hgap/greg200k-sv2-ccs/ref1.fasta"
    aln_fn = "/pbi/dept/secondary/siv/testdata/hgap/greg200k-sv2-ccs/align.bam"

proc main() =
    echo "[INFO] input data:", ref_fn, " ", aln_fn
    #var refx: hts.Fai
    #if not hts.open(refx, ref_fn):
    #    raiseEx(format("Could not open '$#'", ref_fn))
    #assert refx.len() ==  1  # for now
    #let reference_dna = refx.get(refx[0])
    ##assert reference_dna.len() == 5000  # that was for synth5k
    #assert reference_dna.len() > 150000
    #assert reference_dna.len() < 250000
    ##var bamx: hts.Bam
    ##hts.open(bamx, aln_fn, index=true)

    ##phasr.readaln(aln_fn, reference_dna.Dna)
    phasr.main(aln_fn, ref_fn)

main()
