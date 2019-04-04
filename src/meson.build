runtime_sources = ['kmers.c', 'stack.c', 'kmers.h', 'stack.h']

libruntime_lt_c=1
libruntime_lt_r=0
libruntime_lt_a=0

libruntime_so_version = '@0@.@1@.@2@'.format((libruntime_lt_c - libruntime_lt_a),
                  libruntime_lt_a,
                  libruntime_lt_r)


kmer_lib_static = static_library('zkmers',
                            runtime_sources,
                            install : false)
