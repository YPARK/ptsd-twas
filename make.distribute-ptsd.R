#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
source('Util.R')

library(dplyr)
library(readr)
library(methods)

out.dir <- '/broad/hptmp/ypp/twas/gwas_stat/'
system('mkdir -p ' %&&% out.dir)

bim.cols <- c('chr', 'rs', 'snp.loc')
bim.types <- 'ic_i__'

read.bim <- function(...) read_tsv(..., col_names = bim.cols, col_types = bim.types)

read.ptsd <- function(gwas.file, ref.snps.tab) {

    ptsd.cols <- c('rs', 'gwas.a1', 'gwas.a2', 'gwas.beta', 'gwas.se', 'gwas.p')
    ptsd.types <- 'cccddd_'

    gwas.tab <- read_tsv(gwas.file, col_names = ptsd.cols, col_types = ptsd.types, skip = 1) %>%
        right_join(ref.snps.tab, by = 'rs') %>%
            na.omit() %>%
                mutate(gwas.a1 = toupper(gwas.a1), gwas.a2 = toupper(gwas.a2))

}

lift.over <- function(chr.input, gwas.tab, temp.hdr) {

    dir.create(dirname(temp.hdr), recursive = TRUE)

    temp.in.file <- temp.hdr %&&% chr.input %&&% '.bed.gz'
    temp.out.file <- temp.hdr %&&% chr.input %&&% '.lift'
    temp.null.file <- temp.hdr %&&% chr.input %&&% '.null'

    gwas.tab %>% filter(chr == chr.input) %>%
        mutate(snp.1 = as.integer(snp.loc -1)) %>%
            mutate(snp.loc = as.integer(snp.loc)) %>%
                mutate(chr = 'chr' %&&% chr) %>% mutate(prev = snp.loc) %>%
                    select(chr, snp.1, snp.loc, rs, prev) %>%
                        write_tsv(path = temp.in.file, col_names = FALSE)

    sys.lift.cmd <- './bin/liftOver' %&&% ' ' %&&% temp.in.file %&&% ' ' %&&% 'hg19ToHg38.over.chain.gz' %&&% ' ' %&&% temp.out.file %&&% ' ' %&&% temp.null.file

    sys.gzip.cmd <- 'gzip -f ' %&&% temp.out.file

    system(sys.lift.cmd %&&% '; ' %&&% sys.gzip.cmd)

    log.msg('System: ' %&&% sys.lift.cmd)
    log.msg('System: ' %&&% sys.gzip.cmd)

    lift.col <- c('chr', 'remove', 'snp.new', 'rs', 'snp.loc')
    lift.tab <- read_tsv(temp.out.file %&&% '.gz', col_names = lift.col) %>%
        select(-remove) %>%
            mutate(chr = gsub(chr, pattern = 'chr', replacement = '')) %>%
                mutate(chr = as.integer(chr))

    out.tab <- gwas.tab %>%
        filter(chr == chr.input) %>%
            right_join(lift.tab, by = c('chr', 'rs', 'snp.loc')) %>%
                na.omit() %>%
                    select(-snp.loc) %>% mutate(snp.loc = snp.new) %>%
                        mutate(gwas.z = signif(gwas.beta/gwas.se, 4)) %>%
                            select(chr,snp.loc,rs,gwas.a1,gwas.a2,gwas.beta,gwas.se,gwas.p,gwas.z)
    
    sys.rm <- 'rm -f ' %&&% out.dir %&&% gwas.name %&&% '_temp_*'
    log.msg('System: ' %&&% sys.rm)
    system(sys.rm)

    return(as.data.frame(out.tab))
}

## Full pipeline for hg38
proc.gwas <- function(gwas.file, gwas.name, ref.snps.tab) {

    out.files <- out.dir %&&% gwas.name %&&% '/' %&&% 1:22 %&&% '.txt.gz'

    dir.create(dirname(out.files[1]), recursive = TRUE)
    if(all(file.exists(out.files))) return(NULL)

    gwas.tab <- read.ptsd(gwas.file, ref.snps.tab = ref.snps.tab)

    temp <- out.dir %&&% gwas.name %&&% '_temp_'

    lift.list <- lapply(1:22, lift.over, gwas.tab = gwas.tab, temp.hdr = temp)

    dir.create(out.dir %&&% gwas.name, recursive = TRUE)

    for(.chr in 1:22) {
        if(!file.exists(out.files[.chr])) {
            write_tsv(lift.list[[.chr]], path = out.files[.chr])
        }
    }
}

proc.gwas.hg19 <- function(gwas.file, gwas.name, ref.snps.tab) {

    out.files <- out.dir %&&% gwas.name %&&% '/' %&&% 1:22 %&&% '.txt.gz'

    dir.create(dirname(out.files[1]), recursive = TRUE)

    if(all(file.exists(out.files))) return(NULL)

    gwas.tab <- read.ptsd(gwas.file, ref.snps.tab = ref.snps.tab) %>%
        select(chr, snp.loc, rs, gwas.a1, gwas.a2, gwas.beta, gwas.se, gwas.p) %>%
            mutate(gwas.z = gwas.beta / gwas.se)

    for(.chr in 1:22) {
        if(!file.exists(out.files[.chr])) {
            write_tsv(gwas.tab %>% filter(chr == .chr), out.files[.chr])
        }
    }
}

ref.snps.tab <- '1KG_EUR/chr' %&&% 1:22 %&&% '.bim' %>%
    lapply(read.bim) %>%
        bind_rows()

gwas.file <- 'PTSD/military_wave_1.5/ancestry_specific_military/EA/EA_military_7_wave1.51.txt.gz'
gwas.name <- 'ptsd_mil_ea'
proc.gwas(gwas.file, gwas.name, ref.snps.tab)

warnings()

gwas.file <- 'PTSD/civilian_wave_1.5/ancestry_specific_civilian/EA/EA_civilian_5_wave1.51.txt.gz'
gwas.name <- 'ptsd_civ_ea'
proc.gwas(gwas.file, gwas.name, ref.snps.tab)

warnings()

## just use hg19
gwas.file <- 'PTSD/military_wave_1.5/ancestry_specific_military/EA/EA_military_7_wave1.51.txt.gz'
gwas.name <- 'hg19/ptsd_mil_ea'
proc.gwas.hg19(gwas.file, gwas.name, ref.snps.tab)

warnings()

gwas.file <- 'PTSD/civilian_wave_1.5/ancestry_specific_civilian/EA/EA_civilian_5_wave1.51.txt.gz'
gwas.name <- 'hg19/ptsd_civ_ea'
proc.gwas.hg19(gwas.file, gwas.name, ref.snps.tab)

warnings()
