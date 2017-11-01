## The following script constructs the BS dataset and perform smoothing
# This script does not order the positions! Make sure to order them before running this!!!

args = commandArgs(trailingOnly=T)
NAME = args[1]
TESTNS = as.numeric(args[2])
TESTH = as.numeric(args[3])
sitesName = args[4]

msg.trap <- capture.output( suppressMessages( library( bsseq ) ) )
msg.trap <- capture.output( suppressMessages( library( pryr ) ) )

## Reading information
data <- read.table(NAME, sep="\t", header=F)

if(is.unsorted(data$V2, strictly=T)){
  stop("Unsorted count data")
}

BS_all <- BSseq(pos = data$V2, 
                chr = data$V1, 
	        M = as.matrix(data$V4,ncol=1), 
                Cov = as.matrix(data$V5,ncol=1), 
                sampleNames = "V2")

#====Parameters below will be modified per run
pData(BS_all)$Rep <- "V2"
#validObject(BS_all)
#pData(BS_all)

BS_all <- BSmooth(BS_all, 
                  mc.cores = 1, 
                  ns = TESTNS, 
                  h = TESTH, 
                  maxGap=100000, 
                  parallelBy="sample", 
                  verbose=FALSE)
#print(mem_used())
CpGcovered <- granges(BS_all)
myM2 <- round(getMeth(BS_all), 3)

write.table(data.frame(chr=as.vector(CpGcovered@seqnames), start = as.vector(CpGcovered@ranges@start), methylSmoothed=myM2), file=sitesName, col.name=T,row.name=F,quote=F,sep="\t")
#print(gc())
