PREFIX?=${CURDIR}
NIMBLE_DIR?=${CURDIR}/nimbleDir
export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

all:
	${MAKE} sub
	${MAKE} install
quick:
	nim c -r tests/t_kmers.nim
integ:
	${MAKE} -C integ-tests
	nimble integ --debug # slow, for now
help:
	nimble -h
	nimble tasks
test:
	nimble test --debug # uses "tests/" directory by default
install:
	${NIMBLE_INSTALL}
sub:
	#cd repos/cligen; nimble develop  # so we never need to reinstall after edits
	cd vendor/nim-kmers; ${NIMBLE_INSTALL}
	cd vendor/nim-networkx; ${NIMBLE_INSTALL}
pretty:
	find . -name '*.nim' | xargs -L1 nimpretty --indent=4

.PHONY: test
