library(methods)

ttest <- function(a,b){
	results = try(t.test(a,b, alternative= c("two.sided"), paired=TRUE), silent=TRUE)
	if(is(results, "try-error")) return (NA) else return (results$p.value)
}

data = read.table("MethylMatrix.dmrs.freq",T, row.names=1)

using_samples = read.table("label_groups", F, sep="\t", stringsAsFactor=FALSE)

# for meandiff function
controls <- using_samples[which(grepl("ctrl", using_samples$V2)),]
cases <- using_samples[which(grepl("case", using_samples$V2)),]
controls_pos = unlist(sapply(controls$V1, grep, colnames(data)))
cases_pos = unlist(sapply(cases$V1, grep, colnames(data)))
controls_pos = controls_pos[order(controls_pos)]
cases_pos = cases_pos[order(cases_pos)]
print(controls_pos)
print(cases_pos)
data$mypval = apply(data, 1, function(x) ttest(x[controls_pos],x[cases_pos]));
data$myadjpval = p.adjust(data$mypval, method="BH")
write.table(data, file="MethylMatrix.TTest", quote=F, sep="\t");
