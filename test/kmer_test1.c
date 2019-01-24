#include "../src/kmers.h"


int main(int argc, char ** argv){

        char seq[27] = "ATGCCCATCGGGACATATATCCGATC";

        pot_t * kms =  dna_to_kmer(seq, 27, 16);

        print_pot(kms);

}
