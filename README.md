# cgDMR-miner
###a generalized method for the identification of differentially methylated CpG regions (DMRs) on multiple shallow sequencing datasets
## Usage:
perl ./src/cgDMR-miner.pl samplesInfo cpgBedList minDepth pValueCutoff minEffectSize segmentationMode

### Description of input options
1. samplesInfo |  a tab separated file with three columns and one row per sample, (1) sample name, (2) sample id, (3) path to list of methylation frequency files. For each sample a list of methylation frequency file must be generated, which is a tab separated table comprised of three columns and one row per chromosome: (1) sample name matching samples_info file, (2) path to the specific methylation frequency file, (3) chromosome name. Note that both sample names and sample ids should not contain any spaces.
2. cpgBedList | a tab separated file with two columns, (1) chromosome name, (2) path to cpg positions bed file. For each chromosome, a bed file containing the CpG positions to be considered must be provided. If chromosome bed file is missing, then that chromosome will be ignored. Note that chromosome names should not contain any spaces.
3. minDepth: an integer that is the minimum total depth required in each sample for DMR summarization. Default is 10.
4. pValueCutoff | a floating value [0-1] that is the maximum p-value for DMRs. Default is 0.001
5. minEffectSize | a floating value greater or equal to 0 that is the effect size cutoff for DMR calling. This is equals to the square root of the test statistic divided by N (total depth of coverage for all samples). 
6. segmentationMode | "HMM" or "CBS". HMM option will utilize 5-states Hidden Markov Model while CBS option will utilize circular binary segmentation. I recommend HMM for WGBS data and CBS for other genome-wide datasets such as RRBS or capture sequencing data.

## Information about the method

![SFigure1](https://dinhdiep.github.com/cgDMR-miner/img/SFigure1.jpg.jpg) 
![SFigure2](https://dinhdiep.github.com/cgDMR-miner/img/SFigure2.jpg.jpg) 
![SFigure3](https://dinhdiep.github.com/cgDMR-miner/img/SFigure3.jpg.jpg) 
![SFigure4](https://dinhdiep.github.com/cgDMR-miner/img/SFigure4.jpg.jpg) 
![SFigure5](https://dinhdiep.github.com/cgDMR-miner/img/SFigure5.jpg.jpg) 
 
