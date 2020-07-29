# btMSI-CAST

b/tMSI-CAST (blood/tissue MSI Caller Adjusted with Sequence duplicaTes), an NGS-based composite MSI detection method that has a two-mode caller compatible with both tissue and plasma samples and does not require matched normal samples.

# Installation
b/tMSI-CAST is written in Perl. Directly download the entire directory and run msicast.pl to perform MSI dectection on a single test sample. No extra compilation is required.

# Usage
Version 0.10
Usage: perl msicast.pl [options]

## Options
-`database` <string> the directory that contains baseline databases.
-`bam` <string> the path of input BAM file.
-flank <int> spanning reads are defined as fully covering the coordinates of the microsatellite in reference plus $flank bp in both 5’ and 3’ directions. default=2
-num <int> the minimum allele supporting reads post-deduplication, default=2
-minsize <int> minimum family size consdiered during deduplication, default=1
-quality <int> minimum mapping quality, default=1
-mode <string> available modes are 'DelRatio' and 'KLD', default=DelRatio
-outdir <string> output directory
-prefix <string> prefix for output files
 
# Example
An example command to run MSI detection on a single pre-deduplication BAM is:
```
perl $0 -database <database> -bam result/sample20171207-B_I25.sorted.mkdup.realign.bam -flank 2 -num 2 -minsize 1 -quality 1 -mode DelRatio -outdir MSI -prefix sample20171207-B_I25
```
The program takes a pre-deduplication BAM as input, perform duplication removal, call loci-level MSI status and call sample-level MSI status.

Output files include:
1. `prefix`.result.txt
Each entry contains spanning reads information (readID, microsatellite sequence, read sequence, CIGAR and loci name) from the input BAM file.
 
2. `prefix`.DUP.txt
Each row contains the loci name, reference name, aligned starting position, insert size, lengths of each members pertaining to a family.

3. `prefix`.stat.txt
Each row represents info. of all fragments covering a locus, including DeletionRatios/KLD score, dup-ratios, and allele length offsets compared to the germline allele and corresponding frequencies.

4. `prefix`.result.txt
One line sample-level MSI calling results.


# Contact
We can be reached by: Lili Zhao, zhao.lili@genecast.com.cn; Hongyu Xie, xie.hongyu@genecast.com.cn
