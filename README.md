# btMSI-CAST

b/tMSI-CAST (blood/tissue MSI Caller Adjusted with Sequence duplicaTes), an NGS-based composite MSI detection method that has a two-mode caller compatible with both tissue and plasma samples and does not require matched normal samples.


# Usage
Version 0.10
Usage:  msicast [options]

## options
-bin <string>
-bam <string> path of input BAM file
-flank <int> Spanning reads are defined as fully covering the coordinates of the microsatellite in reference plus $flank bp in both 5’ and 3’ directions. default=2
-num <int> the minimum allele supporting reads post-deduplication, default=2
-minsize <int> minimum family size consdiered during deduplication, default=1
-quality <int> minimum mapping quality, default=1
-mode <string> available modes are 'DelRatio' and 'KLD', default=DelRatio
-outdir <string>
-prefix <string> prefix for output files
 
# Example
take pre-deduplication BAMs as input and perform duplication removal by a proprietary majority-voting strategy.

spanning reads information (readID, microsatellite sequence, read sequence, CIGAR and loci name) from the input BAM file is incudled in "result".txt
example.result.txt

# Contact
We can be reached: Lili Zhao, zhao.lili@genecast.com.cn; Hongyu Xie, xie.hongyu@genecast.com.cn
