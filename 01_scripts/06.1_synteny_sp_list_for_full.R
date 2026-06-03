#data preparation for synteny analysis

library(readxl)
library(dplyr)

chrData = "../00_data/all_sp.xlsx"

chrs = read_excel(chrData, sheet = "chr_only", col_names = FALSE)
colnames(chrs) = c("chr", "len", "species")

species = chrs$species
spUniq = unique(species)
spUniq = gsub(pattern = "\\..*$", replacement = "", x = spUniq)
spdf = data.frame(sp1 = spUniq, sp2 = spUniq)

spAll = expand.grid(spdf) %>%
  filter(sp1 != sp2) %>%
  group_by(sp2) %>%
  group_split()

names(spAll) = spdf$sp2

lapply(1:length(spAll), function(x, sps){
  sp = names(sps[x])
  
  df = sps[[x]]
  out = as.character(df$sp1)
  
  outFile = paste0("06_synteny_list_full_", sp, ".txt")
  
  write.table(out, file = outFile, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}, sps = spAll)
