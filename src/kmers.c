#include "kmers.h"

pot_t * dna_to_kmer(char * seq, uint32_t len, uint8_t k){

/**
 * Zero is A, one is C, G is two, and T is 3
 *
 */

        unsigned char seq_nt4_table[256] = {
                0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
        };


        if(len == 0 || seq == NULL || k > 32 ) return NULL;

        uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2*k)-1;

        seed_t * forward_kmer   = malloc(sizeof(seed_t));
        seed_t * reverse_kmer  = malloc(sizeof(seed_t));

        forward_kmer->kmer = 0; forward_kmer->pos = 0;
        reverse_kmer->kmer = 0; reverse_kmer->pos = 0;
        forward_kmer->strand = 0; reverse_kmer->strand = 1;


        Stack * kmer_stack = createStack(sizeof(seed_t), 100);

        // l is the length of the kmers being built on the fly. The variable n is the total number of
        uint32_t i, l, n; i = 0; l = 0; n = 0;

        for(; i < len; i++) {
                int c = seq_nt4_table[(uint8_t)seq[i]];
                //    fprintf(stderr, "%c ", seq[i]);
                if(c < 4) {
                        forward_kmer->kmer = (forward_kmer->kmer << 2  | c ) & mask;
                        reverse_kmer->kmer = (reverse_kmer->kmer >> 2) | (3ULL^c) << shift1;
                        l++;
                }
                else{
                        // advance the window beyond the unknown character
                        l = 0;  i += k;
                        forward_kmer->pos += k; forward_kmer->kmer = 0;
                        reverse_kmer->pos += k; reverse_kmer->kmer = 0;
                }
                if(l >= k) {
                        n+=2;

                        stackPush(kmer_stack, (void *)forward_kmer);
                        stackPush(kmer_stack, (void *)reverse_kmer);
                        forward_kmer->pos +=1; reverse_kmer->pos+=1;
                }
        }

        pot_t  * kmers   = malloc(sizeof(pot_t));
        kmers->seeds     = malloc(sizeof(seed_t)*n);
        kmers->n         = n;
        kmers->word_size = k;
        i = 0;

        while(!stackPop(kmer_stack, (void *) &kmers->seeds[i])) {
                i++;
        }

        stackDestroy(kmer_stack);

        return kmers;
}

char *  bin_to_dna( uint64_t kmer, uint8_t k, uint8_t strand){
        char lookup[4] = {'T', 'G', 'C', 'A'};
        uint8_t mask = 3;
        uint8_t i    = 0;
        uint64_t tmp = 0;
        uint64_t offset = 0;

        char * dna = malloc(sizeof(char)*k);
        for(i = 0; i < k; i++) {
                tmp = kmer;
                offset =  strand == 0 ? (k-i-1)*2 : (i * 2);
                tmp = tmp >> offset;
                dna[i] = lookup[mask & tmp];
        }

        return dna;

}

void print_pot(pot_t * pot){
        uint32_t i = 0;

        for(i = 0; i < pot->n; i += 1) {
                char * dna = bin_to_dna(pot->seeds[i].kmer, pot->word_size, pot->seeds[i].strand);
                fprintf(stderr, "%s\n",  dna );
                free(dna);
        }
}
