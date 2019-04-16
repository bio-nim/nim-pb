# vim: sw=4 ts=4 sts=4 tw=0 et:
import pbpkg/phasr
proc phasr() =
    phasr.foo()
when isMainModule:
    import cligen
    dispatch(phasr)
