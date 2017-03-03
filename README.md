# cgDMR-miner
###a generalized method for the identification of differentially methylated CpG regions (DMRs) on multiple shallow sequencing datasets
## Usage:
`perl ./src/cgDMR-miner.pl samplesInfo cpgBedList minDepth pValueCutoff minEffectSize segmentationMode`

### Description of input options
Input | Explanation of values
------------|-------------------------------------
1. samplesInfo |  a tab separated file with three columns and one row per sample, (1) sample name, (2) sample id, (3) path to list of methylation frequency files. For each sample a list of methylation frequency file must be generated, which is a tab separated table comprised of three columns and one row per chromosome: (1) sample name matching samples_info file, (2) path to the specific methylation frequency file, (3) chromosome name. Note that both sample names and sample ids should not contain any spaces.
2. cpgBedList | a tab separated file with two columns, (1) chromosome name, (2) path to cpg positions bed file. For each chromosome, a bed file containing the CpG positions to be considered must be provided. If chromosome bed file is missing, then that chromosome will be ignored. Note that chromosome names should not contain any spaces.
3. minDepth: an integer that is the minimum total depth required in each sample for DMR summarization. Default is 10.
4. pValueCutoff | a floating value [0-1] that is the maximum p-value for DMRs. Default is 0.001
5. minEffectSize | a floating value greater or equal to 0 that is the effect size cutoff for DMR calling. This is equals to the square root of the test statistic divided by N (total depth of coverage for all samples). 
6. segmentationMode | "HMM" or "CBS". HMM option will utilize 5-states Hidden Markov Model while CBS option will utilize circular binary segmentation. I recommend HMM for WGBS data and CBS for other genome-wide datasets such as RRBS or capture sequencing data.

## Information about the method
### Quantifying methylation variabilities
Methylation variabilities are quantified for artificial examples of ten samples (s1 to s10) across five CpG sites (#1-5). In (a), the acceptable metrics are bolded as they exhibit decrease in methylation variabilities for #1 to #3. Furthermore, Jensen-Shannon distance metric is best able to distinguish 0.1 difference from 1.0 difference, since it has 25 folds difference between #1 and #3. In (b), the acceptable metrics should give the same methylation variabilities for examples #4 and #5. 

![SFigure1](https://github.com/dinhdiep/cgDMR-miner/blob/master/img/SFigure1.jpg.jpg) 

### Clustering five similar methods to cgDMR-miner from simulated 20X depth of coverage data
Five methods at simulated 20-fold average depth of coverage are compared using the DMRs called by each. Distances are 1 – Pearson’s correlation and is calculated on a binary matrix with “1” indicating that a CpG was identified as DMR or “0” indicating that it was not. Only CpGs which have been called DMR in a method are considered. cgDMR-miner is most similar to MethylPy and Hon et al. 2013. Ziller et al. 2013 is the most different, probably due to a higher number of false positives. DSS-general is in-between Ziller et al. 2013 and the other methods.

![SFigure2](https://github.com/dinhdiep/cgDMR-miner/blob/master/img/SFigure2.jpg.jpg) 

### Receiver operating characteristics curves for DMRs with 0.1, 0.2, and 0.4 methylation differences
Receiver operating characteristics curves for cgDMR-miner, our implementation of the Hon et al. 2013 method, and MethylPy for simulated datasets with 5X depth of coverage and DMR methylation differences of 0.10, 0.20, and 0.40.

![SFigure3](https://github.com/dinhdiep/cgDMR-miner/blob/master/img/SFigure3.jpg.jpg) 

### Differentially methylated regions overlap with transcription factor binding sites 
DMRs identified from 36 human tissues are overlapped with a list of all predicted transcription factor binding sites. The frequencies of TF binding are plotted against the distance from the center of a random sampling of 100k segments without replacement. In (a), the DMRs versus non-DMRs segments from cgDMR-miner are compared. In (b), the DMRs from cgDMR-miner are compared against the DMRs from Schultz et al. 2015 using either Ziller et al. 2013 approach or using MethylPy. For the total number of TFBS overlapping DMRs, cgDMR-miner identifies 332,473 (out of 737,084 DMRs called), MethylPy identifies 313,656 (out of 626,418 DMRs called), and Ziller et al. 2013 approach identifies 337,958 (out of 1,198,131 DMRs called). 

![SFigure4](https://github.com/dinhdiep/cgDMR-miner/blob/master/img/SFigure4.jpg.jpg) 

### DHMR example missed by MethylPy
MethylPy missed a DHMR region from the TAB-seq dataset consisting of two normal kidneys and two renal carcinoma samples. This region is found in the gene body of ZNF382, which encodes the transcription factor KZNF. KZNF have been associated with playing a critical role as tumor suppressor in multiple carcinomas. 
![SFigure5](https://github.com/dinhdiep/cgDMR-miner/blob/master/img/SFigure5.jpg.jpg) 
 
