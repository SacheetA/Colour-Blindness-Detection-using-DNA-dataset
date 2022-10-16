# Colour-Blindness-Detection-using-DNA-dataset
Chromosome X of the reference DNA  is compared with the target DNA reads (small segments of original DNA of the target of length ~ 100) to detect colour blindness

Dataset Required: (Not provided here due to copyright reasons)

1. Chromosome X of reference DNA: A string of characters for eg. 'ATATGCCTCGCATCGATAAACTTTCGCGCATACGACATCGA' of length approx ~150 million characters

2. BWT transform of the dataset above (chromosome X) needs to be performed and the correspoding last column and first column of the BWT transform are required for the code to run.

3. The first column of BWT transform has characters in alphabetical order, the location of these characters in the actual reference Chromosome X is required.

4. Reads: The target DNA is divided into segments of length ~100 approx. when DNA sequence need to be generated in a laboratory. Reads are those DNA segments in digital format for eg. 'CTCGCATCGATATACGACATGCCTAACT'. There are ~3 million approx. such reads that need to be compared with the reference chromosome X.

NOTE:

The BWT transform method is used as it reduces time complexity of matching read with reference DNA chromosome X from order(len of reference chromosome ~150m) to order(len of read ~100) for each read.
