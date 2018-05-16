
BLK := 50

CHR := $(shell seq 1 22)

all:

queue: $(foreach chr, $(CHR), jobs/20180516/$(chr).twas.txt.gz)

long: $(foreach chr, $(CHR), jobs/20180516/$(chr).long-twas.txt.gz)

# % = 19
jobs/20180516/%.twas.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)

	awk -vB=$(BLK) -vN=$(shell zcat FQTL-v8/fqtl_$*.txt.gz | wc -l) 'BEGIN { nb = int(N/B); for(b = 0; b <= nb; ++b) print "./make.twas-v8.R" FS "FQTL-v8/fqtl_$*.txt.gz" FS "gwas_stat/ptsd_civ_ea/$*.txt.gz" FS "geno/chr$*" FS ((B * b + 1) ":" (B* (1 + b))) FS ("result/20180516/civ/$*/" b ".txt.gz") }' | gzip > $@

	awk -vB=$(BLK) -vN=$(shell zcat FQTL-v8/fqtl_$*.txt.gz | wc -l) 'BEGIN { nb = int(N/B); for(b = 0; b <= nb; ++b) print "./make.twas-v8.R" FS "FQTL-v8/fqtl_$*.txt.gz" FS "gwas_stat/ptsd_mil_ea/$*.txt.gz" FS "geno/chr$*" FS ((B * b + 1) ":" (B* (1 + b))) FS ("result/20180516/mil/$*/" b ".txt.gz") }' | gzip >> $@

	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N ptsd_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/20180516/%.long-twas.txt.gz: jobs/20180516/%.twas.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=12:00:00 -b y -j y -N ptsd_$*_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@
