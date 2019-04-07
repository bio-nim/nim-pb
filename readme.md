# PB Utilities
So far, these are just Nim re-writes of
http://bitbucket.pacificbiosciences.com:7990/users/zkronenberg/repos/utils/

## Set up

    source module.sh
    export NIMBLE_DIR=$(pwd)/.git/nimble  # recommended

## Test
To run unit-tests in the "tests/" directory:

    nimble test

## Test integration

    nimble integ  # redundant with our unit-tests, but you get the idea

## Install

    nimble install

## Debug and develop

    make
    make integ
