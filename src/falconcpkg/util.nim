# vim: sts=4:ts=4:sw=4:et:tw=0
from cpuinfo import nil
from os import nil
from threadpool import nil

type PbError* = object of Exception

proc raiseEx*(msg: string) {.discardable.} =
    raise newException(PbError, msg)

proc isEmptyFile*(fin: string): bool =
    var finfo = os.getFileInfo(fin)
    if finfo.size == 0:
        return true
    return false

template withcd*(newdir: string, statements: untyped) =
    let olddir = os.getCurrentDir()
    os.setCurrentDir(newdir)
    defer: os.setCurrentDir(olddir)
    statements

proc log*(words: varargs[string, `$`]) =
    for word in words:
        write(stderr, word)
    write(stderr, '\l')

proc adjustThreadPool*(n: int) =
    ## n==0 => use ncpus
    var size = n
    if n == 0:
        size = cpuinfo.countProcessors()
    if size > threadpool.MaxThreadPoolSize:
        size = threadpool.MaxThreadPoolSize
    log("ThreadPoolsize=#, MaxThreadPoolSize=#, NumCpus=#",
        size, threadpool.MaxThreadPoolSize, cpuinfo.countProcessors())
    threadpool.setMaxPoolSize(size)