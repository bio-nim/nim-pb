PREFIX?=${CURDIR}
NIMBLE_DIR?=${CURDIR}/nimbleDir
export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

test-kmers:
	nim c -r tests/t_kmers
help:
	nimble -h
	nimble tasks
test:
	nimble test --debug # uses "tests/" directory by default
integ:
	nimble integ --debug
install:
	${NIMBLE_INSTALL}

.PHONY: test
