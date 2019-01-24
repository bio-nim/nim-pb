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
                        fprintf(stderr, "WTF\n" );
                        // advance the window beyond the unknown character
                        l = 0;  i += k;
                        forward_kmer->pos += k; forward_kmer->kmer = 0;
                        reverse_kmer->pos += k; reverse_kmer->kmer = 0;
                }
                if(l >= k) {
                        n+=2;
                        fprintf(stderr, "k:%i l:%i p1:%i p2:%i\n", k, l, forward_kmer->pos, reverse_kmer->pos );
                        stackPush(kmer_stack, (void *)forward_kmer);
                        stackPush(kmer_stack, (void *)reverse_kmer);
                        forward_kmer->pos +=1; reverse_kmer->pos+=1;
                }
        }
        fprintf(stderr, "\n" );

        pot_t  * kmers = malloc(sizeof(pot_t));
        kmers->seeds = malloc(sizeof(seed_t)*n);
        kmers->n     = n;
        i = 0;

        while(!stackPop(kmer_stack, (void *) &kmers->seeds[i])) {
                i++;
        }

        stackDestroy(kmer_stack);

        return kmers;
}

void print_kmer(uint64_t kmer, uint8_t k, uint8_t strand){
        char * lookup[4] = {'T', 'G', 'C', 'A'};
        uint8_t mask = 3;
        uint8_t i    = 0;
        uint64_t tmp = 0;
        uint64_t offset = 0;

        for(i = 0; i < k; i++) {
                tmp = kmer;
                offset =  strand == 0 ? (k-i-1)*2 : (i * 2);
                tmp = tmp >> offset;
                fprintf(stderr, "%c", lookup[mask & tmp]);

        }

}

void print_pot(pot_t * pot){
        uint32_t i = 0;

        fprintf(stderr, "number of kmer seeds: %i\n", pot->n + 1);


        for(i = 0; i < pot->n - 1; i += 2) {
                fprintf(stderr, "%i kmer:%llu position:%i strands:%i%i\n", i, pot->seeds[i].kmer, pot->seeds[i].pos, pot->seeds[i].strand, pot->seeds[i+1].strand );
                fprintf(stderr, "%i kmer:%llu position:%i strands:%i%i\n", i, pot->seeds[i+1].kmer, pot->seeds[i].pos, pot->seeds[i].strand, pot->seeds[i+1].strand );
                print_kmer(pot->seeds[i].kmer, 16, pot->seeds[i].strand);
                fprintf(stderr, "\n" );
                print_kmer(pot->seeds[i+1].kmer, 16, pot->seeds[i+1].strand);
                fprintf(stderr, "\n" );
        }


}
