#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 5) {
    q()
}

fqtl.stat.file <- argv[1]             # e.g., 'FQTL-v6/chr1/50/combined.txt.gz'
gwas.stat.file <- argv[2]             # e.g., 'gwas_stat/hg19/ptsd_mil_ea/1.txt.gz'
geno.hdr <- argv[3]                   # e.g., 'geno_v6/chr1'
gg.vec <- eval(parse(text = argv[4])) # e.g., 1:5
out.file <- argv[5]                   # e.g., ''

dir.create(dirname(out.file), recursive = TRUE)

################################################################

library(dplyr)
library(readr)
library(tidyr)
library(zqtl)
source('Util.R')
source('NWAS.R')

if(file.exists(out.file)) {
    log.msg('file exists: %s', out.file)
    q()
}

################################################################

cis.dist <- 1e6
n.cutoff <- 10
n.perm <- 5e7
n.blk <- 1024
n.round <- ceiling(n.perm/n.blk)

temp.dir <- system('mkdir -p /broad/hptmp/ypp/ptsd-twas/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/ptsd-twas/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

################################################################
fqtl.cols <- c('ensg', 'chr', 'tss', 'tis.idx', 'tis.name', 'tis.beta', 'tis.se', 'tis.lo',
               'snp.name', 'snp.beta', 'snp.se', 'snp.lo', 'factor', 'lo.cutoff')

fqtl.stat <- read_tsv(fqtl.stat.file, col_names = fqtl.cols) %>%
    filter(lo.cutoff == 0.9)

snp.stat <- fqtl.stat %>%
    select(ensg, factor, snp.name, snp.beta, snp.se)

gwas.stat <- read_tsv(gwas.stat.file)

take.twas <- function(gg) {

    tss <- fqtl.stat[gg, 'tss'] %>% .unlist() %>% as.integer()
    info <- fqtl.stat %r% gg %>% select(ensg, factor)

    .split <- function(...) .unlist(...) %>% strsplit(split = '[|]') %>% (function(x) x[[1]])

    qtl.info <-
        data.frame(snp = snp.stat[gg, 'snp.name'] %>% .split(),
                   qtl.beta = snp.stat[gg, 'snp.beta'] %>% .split() %>% as.numeric(),
                   qtl.se = snp.stat[gg, 'snp.se'] %>% .split() %>% as.numeric()) %>%
                       separate(snp, c('chr', 'snp.loc', 'qtl.a1', 'qtl.a2', 'ver'), sep = '_') %>%
                           select(-ver) %>%
                               mutate(snp.loc = as.integer(snp.loc),
                                      chr = as.integer(chr))

    plink <- subset.plink(geno.hdr,
                          chr = qtl.info[1,'chr'],
                          plink.lb = min(qtl.info$snp.loc),
                          plink.ub = max(qtl.info$snp.loc),
                          temp.dir)

    if(is.null(plink)) return(NULL)

    matched <- match.qtl.gwas(plink, gwas.stat, qtl.info, tss, cis.dist)
    if(nrow(matched) < 1) return(NULL)

    qtl.z.tab <- matched %>% filter(!is.na(qtl.z))
    if(nrow(qtl.z.tab) < 1) return(NULL)

    svd.out <- zqtl::take.ld.svd(plink$BED, options = list(eig.tol = 1e-2, do.stdize = TRUE))

    V.t <- svd.out$V.t %c% matched$x.pos
    D <- svd.out$D

    V.t.obs <- svd.out$V.t %c% qtl.z.tab$x.pos
    gwas.z <- qtl.z.tab$gwas.z
    qtl.z <- qtl.z.tab$qtl.z
    gwas.z.tot <- matched$gwas.z

    obs.stat <- func.NWAS(qtl.z, gwas.z, V.t.obs, D)

    blk.mat <- func.blk.ind(nrow(qtl.z.tab), n.blk)

    ## adaptive permutation
    c.tot <- 0
    p.tot <- 0
    z.obs.abs <- abs(obs.stat[, 'z'])
    for(b in seq(1, n.round)){
        set.seed(b)

        stat.z <- func.NWAS.qtl.perm(qtl.z, gwas.z.tot, V.t, D, blk.mat) %>%
            select(z) %>% .unlist()

        c.tot <- c.tot + sum((abs(stat.z) + 1e-8) >= z.obs.abs)
        p.tot <- p.tot + length(stat.z)

        if(c.tot > n.cutoff){
            break;
        }
        cat('\n', c.tot, '/', p.tot, '... ')
    }

    log.msg('Finished: %d, %s', gg, info[1])

    ret <- bind_cols(info, obs.stat) %>% mutate(p.val = (1 + c.tot)/(1 + p.tot))
    return(ret)
}

gg.vec <- pmin(gg.vec, nrow(fqtl.stat)) %>% unique()

out.tab <- lapply(gg.vec, take.twas) %>% bind_rows() %>%
    mutate(z = signif(z, 4),
           theta = signif(theta, 4),
           theta.se = signif(theta.se, 4),
           p.val = signif(p.val, 4))

if(dir.exists(temp.dir)) { system('rm -r ' %&&% temp.dir) }
write_tsv(out.tab, path = out.file)
log.msg('Finished')
