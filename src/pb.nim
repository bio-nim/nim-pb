# vim: sw=4 ts=4 sts=4 tw=0 et:
import pbpkg/zev
proc zev() =
    zev.foo()
when isMainModule:
    import cligen
    dispatch(zev)
