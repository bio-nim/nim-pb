PREFIX?=${CURDIR}
NIMBLE_DIR?=${CURDIR}/nimbleDir
export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

quick:
	nim c -r tests/t_kmers.nim
integ:
	${MAKE} -C test
help:
	nimble -h
	nimble tasks
test:
	nimble test --debug # uses "tests/" directory by default
install:
	${NIMBLE_INSTALL}
pretty:
	find . -name '*.nim' | xargs -L1 nimpretty --indent=4

.PHONY: test
