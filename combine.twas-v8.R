#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)
source('Util.R')

library(dplyr)
library(tidyr)
library(readr)

read.twas <- function(chr, data, result.dir) {
    result.files <- list.files(path = result.dir %&&% '/' %&&% data %&&% '/' %&&% chr %&&% '/',
                               pattern = '.txt.gz', full.names = TRUE)
    twas.tab <- result.files %>% lapply(FUN = function(f) { suppressMessages(read_tsv(f)) }) %>%
        bind_rows()
    ret <- twas.tab %>%
        simplify.ensg() %>%
            mutate(data)
}

read.tis.tab <- function(chr) {
    fqtl.stat.file <- 'FQTL-v8/fqtl_' %&&% chr %&&% '.txt.gz'
    tis.stat <- suppressMessages(read_tsv(fqtl.stat.file)) %>%
        simplify.ensg() %>%
            select(ensg, factor, tis, tis.theta, tis.sd, tis.lodds)
}

twas.mil.tab <- lapply(1:22, read.twas, data = 'mil', result.dir = 'result/20180516') %>%
    bind_rows() %>%
        select(ensg, factor, z, theta, theta.se, p.val) %>%
            rename(z.mil = z, twas.mil = theta, se.mil = theta.se, pval.mil = p.val)

twas.civ.tab <- lapply(1:22, read.twas, data = 'civ', result.dir = 'result/20180516') %>%
    bind_rows() %>%
        select(ensg, factor, z, theta, theta.se, p.val) %>%
            rename(z.civ = z, twas.civ = theta, se.civ = theta.se, pval.civ = p.val)

tis.stat <- lapply(1:22, read.tis.tab) %>% bind_rows()

twas.tab <- twas.mil.tab %>%
    left_join(twas.civ.tab) %>%
        left_join(tis.stat) %>%
            mutate(z.delta = (z.mil - z.civ) / sqrt(se.mil^2 + se.civ^2)) %>%
                mutate(p.delta = 2 * pnorm(abs(z.delta), lower.tail = FALSE))

################################################################
## protein coding
twas.coding <- read.coding.genes() %>%
    left_join(twas.tab) %>%
        left_join(tis.stat) %>%
            na.omit()

tis.names <- read_tsv('FQTL-v8-tissues.txt',
                      col_names = c('tis.name','tis.size'),
                      col_types = '_ci_') %>%
                          mutate(tis = 1:n())

## lineararized data
.split.bar <- function(s) strsplit(s, split = '[|]')

twas.coding.linear <-
    twas.coding %>%
        unnest(tis = .split.bar(tis),
               tis.theta = .split.bar(tis.theta),
               tis.sd = .split.bar(tis.sd),
               tis.lodds = .split.bar(tis.lodds)) %>%
                   mutate(tis = as.integer(tis),
                          tis.theta = as.numeric(tis.theta),
                          tis.sd = as.numeric(tis.sd),
                          tis.lodds = as.numeric(tis.lodds)) %>%
                   left_join(tis.names) %>%
                       filter(tis.lodds > 0)

.mean <- function(...) mean(...) %>% signif(digits = 2)
.min <- function(...) min(...) %>% signif(digits = 2)

twas.coding.out <- twas.coding.linear %>%    
    mutate(tis.z = signif(tis.theta / tis.sd, 2)) %>%
        mutate(tis.pip = signif(1/(1+exp(-tis.lodds)), 2)) %>%
            group_by(chr, lb, ub, strand, ensg, factor, hgnc) %>%
                summarize(z.mil = .mean(z.mil),
                          twas.mil = .mean(twas.mil),
                          se.mil = .mean(se.mil),
                          pval.mil = .min(pval.mil),
                          z.civ = .mean(z.civ),
                          twas.civ = .mean(twas.civ),
                          se.civ = .mean(se.civ),
                          pval.civ = .min(pval.civ),
                          tis.name = paste(tis.name, collapse = '|'),
                          tis.z = paste(tis.z, collapse = '|'),
                          tis.pip = paste(tis.pip, collapse = '|'))

################################################################
write_tsv(twas.coding.out, path = 'result/twas.gtex_v8.txt.gz')
