# vim: sw=4 ts=4 sts=4 tw=0 et:

proc dataset(extras: seq[string]) =
    echo "pb dataset"
proc kmers(int_dummy: int = 42, string_dummy: string = "hello") =
    echo "pb kmers"
proc utils(extras: seq[string], float_req: float) =
    echo "pb utils ..."

when isMainModule:
    import cligen
    dispatchMulti(
        [dataset, short={}, help={}],
        [kmers, short={"int_dummy": 'd'}, help={}],
        [utils, short={}, help={"float_req": "special help message"}],
    )
