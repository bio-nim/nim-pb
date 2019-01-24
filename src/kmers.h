#ifndef KMERS_H
#define KMERS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include "stack.h"

/**
 * kmer - a uint64 supporting a maximum of 32 DNA bases.
 * pos  - position along the sequence
 */
typedef struct seed {
        uint64_t kmer;
        uint32_t pos;
        uint8_t strand;
} seed_t;

/**
 * seeds - a pointer to the kmers
 * n  - the number of kmers in the database
 */
typedef struct seed_holder {
        seed_t * seeds;
        uint32_t n;
} pot_t;


/**
 * Converts a char * into a set of seed_t objects.
 * @param  seq  - char * sequence
 * @param  len  - length of the sequence
 * @param  k - kmer size
 * @return      [description]
 */
pot_t * dna_to_kmer (char * seq, uint32_t len, uint8_t k);

/**
 * Prints the pot structure to STDOUT
 * @param pot a pointer to the pot
 */
void print_pot(pot_t * pot);

#ifdef __cplusplus
}
#endif

#endif /* KMERS_H */
