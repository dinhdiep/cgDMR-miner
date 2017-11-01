cgDMRminer_dir='/media/2TB_store1/cgDMR-miner_v1.0'
all_methyl_freq='mf_list'

grep -v NOPASS chr*.results.txt | grep -v index | cut -f 1 | sed 's/,/\t/g' > dmrs.seg.txt
perl $cgDMRminer_dir/src/allMethyl2Matrix.pl $all_methyl_freq 20 4 dmrs.seg.txt freq dmr dmrs
Rscript $cgDMRminer_dir/src/getSignifTTest.R


