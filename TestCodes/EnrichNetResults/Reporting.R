library(dplyr)
library(stringr)
library(xtable)
E.Net.cc <- read.table("OC_enrichnet_ranking_table.txt",
                                    header = T)


E.Net.cc <- as_data_frame(E.Net.cc)

E.Net.cc %<>% mutate(.,ID = str_sub(Annotation..pathway.process.,
                                   start = 4, end = 8),
                      Pathway_Name =
                           str_sub(Annotation..pathway.process.,start = 10)) %>%
             select(.,c("Pathway_Name","ID","XD.score","Fisher.q.value" )) %>%
             filter(., XD.score >= 1.12)

E.Net.report     <- E.Net.cc
E.Net.report[,4] <- mapply(formatC,E.Net.report[,4],
                             MoreArgs = list(format = "e", digits = 2))

print(xtable(E.Net.report), include.rownames = FALSE)
