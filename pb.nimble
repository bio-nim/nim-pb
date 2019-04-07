# Package

version       = "0.1.0"
author        = "Zev Kronenberg"
author        = "Christopher Dunn"
description   = "Utilities for bioinformatics"
license       = "BSD-3-Clause"
srcDir        = "src"
installDirs   = @["pbpkg"]
bin           = @["pb"]


# Dependencies

requires "nim >= 0.19.6"

task integ, "Runs integration tests":
  let cmd = "nim c -r test/kmer_test1"
  echo cmd
  exec cmd
