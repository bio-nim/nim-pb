default: test1
test1:
	nim c -d:debug kmer_test1
	./kmer_test1 >| kmer_test1.result.txt
	diff kmer_test1.expect.txt kmer_test1.result.txt
