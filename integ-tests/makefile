default: test1 test-phasr
test1:
	nim c kmer_test1
	./kmer_test1 >| kmer_test1.result.txt
	diff kmer_test1.expect.txt kmer_test1.result.txt
test-phasr:
	nim c phasr_test
	./phasr_test
