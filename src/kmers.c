#include "kmers.h"

int _compare_kmers(const void * a, const void * b){

	uint64_t c = (*(seed_t *)a).kmer;
	uint64_t d = (*(seed_t *)b).kmer;

	if (c < d ) return -1;
	// This makes sure they are never the same, meaning stable sort!
	if ( c == d ) {
		if( (*(seed_t *)a).pos < (*(seed_t *)b).pos) {
			return -1;
		}
		else{
			return 0;
		}
	}
	return 1;
}


const unsigned char seq_nt4_table[] = {
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


pot_t * dna_to_kmer(char * seq, uint32_t len, uint8_t k){

/**
 * Zero is A, one is C, G is two, and T is 3
 */



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

	pot_t  * kmers    = malloc(sizeof(pot_t));
	kmers->seeds      = malloc(sizeof(seed_t)*n);
	kmers->n          = n;
	kmers->word_size  = k;
	kmers->searchable = 0;
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
		fprintf(stderr, "pos:%i strand:%i seq:%s\n",  pot->seeds[i].pos, pot->seeds[i].strand, dna );
		free(dna);
	}
}

int make_searchable(pot_t * kms){
	if(kms->searchable) return 1;
	//sorts the kmers in the pot, not really required and could be removed.
	qsort(kms->seeds, kms->n, sizeof(seed_t), _compare_kmers);
	kms->h = kh_init(z);

	int ret;
	uint32_t i = 0;
	khiter_t k;

	for(i = 0; i < kms->n; i++) {
		k = kh_get(z, kms->h, kms->seeds[i].kmer);
		if(k == kh_end(kms->h)) {
			k = kh_put(z, kms->h, kms->seeds[i].kmer, &ret );
			kh_value(kms->h, k) = i;
		}
	}

	kms->searchable = 1;
	return 0;
}

Stack * search(pot_t * target, pot_t * query){
	make_searchable(target);

	fprintf(stderr, "Searching through %i kmers\n", query->n );

	Stack * hit_stack = createStack(sizeof(seed_pair_t), 100);

	seed_pair_t hit;

	int ret;
	uint32_t i = 0;
	khiter_t k;
	uint32_t hit_index;

	for(i = 0; i < query->n; i++) {
		k = kh_get(z, target->h, query->seeds[i].kmer);
		if(k != kh_end(target->h)) {
			hit_index = kh_value(target->h, k);
			while(query->seeds[i].kmer == target->seeds[hit_index].kmer && hit_index < query->n) {

				hit.a = &query->seeds[i];
				hit.b = &target->seeds[hit_index];;
				stackPush(hit_stack, (void*)&hit);

				hit_index+=1;
			}
		}
	}
	return hit_stack;
}
