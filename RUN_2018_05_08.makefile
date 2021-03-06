
BLK := 100

CHR := $(shell seq 1 22)

all:

queue: $(foreach chr, $(CHR), jobs/20180508/$(chr).twas.txt.gz)

long: $(foreach chr, $(CHR), jobs/20180508/$(chr).long-twas.txt.gz)

# % = 19
jobs/20180508/%.twas.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)

	awk -vB=$(BLK) -vN=$(shell zcat FQTL-v6/chr$*/50/combined.txt.gz | awk '$$NF == .9' | wc -l) 'BEGIN { nb = int(N/B); for(b = 0; b <= nb; ++b) print "./make.twas-v6.R" FS "FQTL-v6/chr$*/50/combined.txt.gz" FS "gwas_stat/hg19/ptsd_civ_ea/$*.txt.gz" FS "geno_v6/chr$*" FS ((B * b + 1) ":" (B* (1 + b))) FS ("result/20180508/civ/$*/" b ".txt.gz") }' | gzip > $@

	awk -vB=$(BLK) -vN=$(shell zcat FQTL-v6/chr$*/50/combined.txt.gz | awk '$$NF == .9' | wc -l) 'BEGIN { nb = int(N/B); for(b = 0; b <= nb; ++b) print "./make.twas-v6.R" FS "FQTL-v6/chr$*/50/combined.txt.gz" FS "gwas_stat/hg19/ptsd_mil_ea/$*.txt.gz" FS "geno_v6/chr$*" FS ((B * b + 1) ":" (B* (1 + b))) FS ("result/20180508/mil/$*/" b ".txt.gz") }' | gzip >> $@

	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N ptsd_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/20180508/%.long-twas.txt.gz: jobs/20180508/%.twas.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=12:00:00 -b y -j y -N ptsd_$*_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@
