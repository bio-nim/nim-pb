#include "../src/kmers.h"
#include <string.h>

int main(int argc, char ** argv){

	char seq[13] = "ATCGGCTACTATT";

	fprintf(stderr, "Starting seq: %s\n", seq );

	char ans_lookups[6][13] = {
		"TCGGCTACTATT",
		"ATCGGCTACTAT",
		"TAGCCGATGATA",
		"AGCCGATGATAA"
	};

	pot_t * kms =  dna_to_kmer(seq, 13, 12);
	pot_t * qms =  dna_to_kmer(seq, 13, 12);

	fprintf(stderr, "%s\n", "kms" );
	print_pot(kms);
	fprintf(stderr, "\n\n" );
	fprintf(stderr, "%s\n", "qms" );
	print_pot(qms);

	make_searchable(kms);

	uint32_t i;

	int final_res = 0;

	for(i = 0; i < kms->n; i++) {
		char * tmp = bin_to_dna(kms->seeds[i].kmer, kms->word_size, kms->seeds[i].strand);
		int res = strcmp(tmp, ans_lookups[i]);

		fprintf(stderr, "[%s:DNA->BIT->DNA:TEST] kmer: %llu pos:%i expecting:%s observed:%s [%s]\n", __FILE__, kms->seeds[i].kmer, kms->seeds[i].pos, ans_lookups[i], tmp, res != 0 ? "FAIL" : "PASS" );

		if(res != 0 ) final_res = 1;

		free(tmp);
	}


	Stack * hits = search(kms, qms);

	seed_pair_t pair;

	while(!stackPop(hits, (void*)&pair)) {
		fprintf(stderr, "qb:%i tb:%i qs:%i ts:%i %s %s\n", pair.a->pos, pair.b->pos, pair.a->strand, pair.b->strand, bin_to_dna(pair.a->kmer, kms->word_size, pair.a->strand), bin_to_dna(pair.b->kmer, kms->word_size, pair.b->strand) );
	}

	return final_res;

}
