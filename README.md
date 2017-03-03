# cgDMR-miner
###a generalized method for the identification of differentially methylated CpG regions (DMRs) on multiple shallow sequencing datasets
## Usage:
./src/cgDMR-miner.pl samples_info cpg_sites_bed_list min_depth p_value_cutoff minimum_effect_size

## Input files:
1. samples_info: a tab separated file with three columns and one row per sample, (1) sample name, (2) sample id, (3) path to list of methylation frequency files. For each sample a list of methylation frequency file must be generated, which is a tab separated table comprised of three columns and one row per chromosome: (1) sample name matching samples_info file, (2) path to the specific methylation frequency file, (3) chromosome name. Note that both sample names and sample ids should not contain any spaces.
2. cpg_sites_bed_list: a tab separated file with two columns, (1) chromosome name, (2) path to cpg positions bed file. For each chromosome, a bed file containing the CpG positions to be considered must be provided. If chromosome bed file is missing, then that chromosome will be ignored. Note that chromosome names should not contain any spaces.
3. min_depth: an integer that is the minimum total depth required in each sample for DMR summarization. Default is 10.
4. p_value_cutoff: a floating value [0-1] that is the maximum p-value for DMRs. Default is 0.001
5. minimum_effect_size: a floating value greater or equal to 0 that is the effect size cutoff for DMR calling. This is equals to the square root of the test statistic divided by N (total depth of coverage for all samples). 
 
