#ifndef KMERS_H
#define KMERS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include "stack.h"
#include "khash.h"
#include <stdlib.h>     /* qsort */

KHASH_MAP_INIT_INT(z, uint32_t);

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
 * a & b are two seed_t's designed for matching in the hash lookup
 */
typedef struct seed_pair {
	seed_t * a;
	seed_t * b;
} seed_pair_t;

/**
 * seeds - a pointer to the kmers
 * n  - the number of kmers in the database
 */
typedef struct seed_holder {
	size_t word_size;
	uint32_t n;
	seed_t * seeds;
	khash_t(z) *h;
	uint8_t searchable;
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

/**
 * A function to convert the binary DNA back into character
 * @param target - char * to load
 * @param kmer   [description]
 * @param k      [description]
 * @param strand [description]
 */
char * bin_to_dna(uint64_t kmer, uint8_t k, uint8_t strand);

/**
 * A function that sorts and loads the kmers into a hash table and sets the sort
 * flag.
 * @param  pot_t * - a pointer to a pot
 * @return int - zero upon success
 */
int make_searchable(pot_t *);

Stack * search(pot_t *, pot_t *);

#ifdef __cplusplus
}
#endif

#endif /* KMERS_H */
