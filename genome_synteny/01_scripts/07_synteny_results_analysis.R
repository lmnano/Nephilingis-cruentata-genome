#synteny analysis for all species

library(dplyr)
library(stringr)
library(ggplot2)
library(readxl)
library(tidyr)


chrData = "./00_data/all_sp.xlsx"
tableFile = "../00_data/06.1_cactus_synteny_tables_full/"

tspeciesDir = "../00_data/06.1_cactus_synteny_full/"
tspeciesDirs = list.dirs(path = tspeciesDir, full.names = FALSE, recursive = FALSE)

lapply(tspeciesDirs, function(tspecies){
  #get files
  
  pathSp = paste0("./00_data/06.1_cactus_synteny_full/", tspecies, "/")
  
  pslList = list.files(pathSp, full.names = TRUE)
  pslNames = list.files(pathSp)
  
  #psl format col names
  psl_colnames <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", 
                    "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", 
                    "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", 
                    "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts")
  
  #psl import
  psls = lapply(pslList, read.table, sep = "\t", col.names = psl_colnames)
  names(psls) = pslNames
  
  #import chr data
  
  chrs = read_excel(chrData, sheet = "chr_only", col_names = FALSE)
  colnames(chrs) = c("chr", "len", "species")
  chrs_type = read_excel(chrData, sheet = "chr_type_no_synteny", col_names = TRUE)
  colnames(chrs) = c("chr", "len", "species")
  chrs_num = read_excel(chrData, sheet = "chr_num_no_synteny", col_names = FALSE)
  colnames(chrs_num) = c("chr", "len", "species", "chr_num")
  
  percChrNoLenFilter = lapply(1:length(psls), function(x, psllst, chrdf, chrt){

    pslname = names(psllst[x])
    
    sps = gsub(pattern = "06.1_synteny_", replacement = "", x = pslname)
    sps = gsub(pattern = ".psl$", replacement = "", x = sps)
    target = gsub(pattern = "(.*)_\\w+$", replacement = "\\1", x = sps)
    target = gsub(pattern = "(.*)_\\w+$", replacement = "\\1", x = target)
    query = gsub(pattern = "^\\w+_(.*_.*$)", replacement = "\\1", x = sps)
    
    chrdfFilter = chrdf %>%
      mutate(chr = gsub(pattern = " .*$", replacement = "", x = chr)) %>%
      filter(species == target | species == query)
    
    chrtFilter = chrt %>%
      mutate(chr = gsub(pattern = " .*$", replacement = "", x = chr)) %>%
      filter(species == target | species == query)
    
    chrtTarget = chrtFilter %>%
      filter(species == target) %>%
      select(-c("species", "len")) %>%
      rename("chrt" = "chr_type")
    
    chrtQuery = chrtFilter %>%
      filter(species == query) %>%
      select(-c("species", "len")) %>%
      rename("chrq" = "chr_type")
    
    df = psllst[[x]] %>%
      filter(tName %in% chrdfFilter$chr) %>%
      filter(qName %in% chrdfFilter$chr) %>%
      mutate(tLength = tEnd - tStart) %>%
      mutate(tPerc = tLength / tSize) %>%
      group_by(tName, qName) %>%
      summarise(percSum = sum(tPerc)) %>%
      arrange(percSum, .by_group = TRUE)
    
    df2tname = psllst[[x]] %>%
      filter(tName %in% chrdfFilter$chr) %>%
      filter(qName %in% chrdfFilter$chr) %>%
      count(tName, name = "tcount")
    
    df2 = psllst[[x]] %>%
      filter(tName %in% chrdfFilter$chr) %>%
      filter(qName %in% chrdfFilter$chr) %>%
      count(tName, qName, name = "tqcount") %>%
      left_join(df2tname) %>%
      mutate(percCount = tqcount / tcount) %>%
      select(-c("tcount", "tqcount"))
    
    statsdf1 = df %>%
      left_join(df2) %>%
      left_join(chrtTarget, by = c("tName" = "chr")) %>%
      left_join(chrtQuery, by = c("qName" = "chr")) %>%
      unite(chrt, chrq, sep = "", col="chrsyn") %>%
      group_by(chrsyn)
    
    df = df %>%
      left_join(df2) %>%
      filter(percCount >= 0.1)
    
    return(df)
    
    
  }, psllst = psls, chrdf = chrs, chrt = chrs_type)
  names(percChrNoLenFilter) = pslNames
  
  chrMatch = lapply(1:length(percChrNoLenFilter), function(x, dfs, chrdf){
    
    pslname = names(dfs[x])
    
    sps = gsub(pattern = "06.1_synteny_", replacement = "", x = pslname)
    sps = gsub(pattern = ".psl$", replacement = "", x = sps)
    target = gsub(pattern = "(.*)_\\w+$", replacement = "\\1", x = sps)
    target = gsub(pattern = "(.*)_\\w+$", replacement = "\\1", x = target)
    query = gsub(pattern = "^\\w+_(.*_.*$)", replacement = "\\1", x = sps)
    
    chrData = chrdf %>%
      mutate(chrCode = gsub(pattern = " .*$", replacement = "", x = chr)) %>%
      select(c("species", "chr_num", "chrCode"))
    
    
    df = dfs[[x]] %>%
      group_by(tName) %>%
      mutate(maxMatch = max(percCount)) %>%
      filter(percCount == maxMatch) %>%
      left_join(chrData, by = c("tName" = "chrCode")) %>%
      rename(tSpecies = species, tNum = chr_num) %>%
      left_join(chrData, by = c("qName" = "chrCode")) %>%
      rename(qSpecies = species, qNum = chr_num) %>%
      ungroup() %>%
      select(c("tNum", "qNum", "maxMatch"))
    
    matchName = paste(target, query, "percCountMax", sep = "-")
    dfNames = c(target, query, matchName)
    colnames(df) = dfNames
    
    return(df)
    
  }, dfs = percChrNoLenFilter, chrdf = chrs_num)
  
  chrMatchAll = chrMatch %>%
    purrr::reduce(left_join) %>%
    t()
  outMatchPath = paste0(tableFile, tspecies, "/", tspecies, "_chr_match_table.tsv")
  write.table(chrMatchAll, file = outMatchPath, row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  
})



