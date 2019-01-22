#include "kmers.h"

seed_t * dna_to_kmer(char * seq, uint32 len, uint8_t k){

        if(len == 0 || seq == NULL) return NULL;

        uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2*k)-1;

        seed_t * foward_kmer   = malloc(sizeof(seed_t));
        seed_t * reverse_kmer  = malloc(sizeof(seed_t));

        foward_kmer.kmer = 0; forward_kmer.pos = 0;
        reverse_kmer.kmer = 0; reverse_kmer.pos = 0;
        forward_kmer.strand = 0; reverse_kmer.strand = 1;


        kmer_stack * createStack(sizeof(seed_t), 100);

        // l is the length of the kmers being built on the fly. The variable n is the total number of
        uint32 i, l, n; i = 0; l = 0; n = 0;

        for(; i < len; i++) {
                int c = seq_nt4_table[(uint8_t)seq[i]];
                if(c < 4) {
                        forward_kmer.kmer = (forward_kmer.kmer << 2  | c ) & mask;
                        kmers.km[1] = (kmers.km[1] >> 2) | (3ULL^c) << shift1;
                        foward_kmer.pos++; reverse_kmer.pos++; l++;
                }
                else{
                        // advance the window beyond the unknown character
                        l = 0;  i += k;
                        foward_kmer.pos += k; forward_kmer.kmer = 0;
                        reverse_kmer.pos += k; reverse_kmer.kmer = 0;
                }
                if(l >= k) {
                        n++;
                        stackPush(kmer_stack, (void *)forward_kmer);
                        stackPush(kmer_stack, (void *)reverse_kmer);
                }
        }

        seed_t * kmers = malloc(sizeof(seed_t)*n);

        return kmers;
}
