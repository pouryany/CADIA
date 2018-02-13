install.packages("rvest")
library(rvest)


webpage <- read_html("http://www.genome.jp/dbget-bin/get_linkdb?-t+pathway+gn:T01001")

webpage %>% html_nodes("pre") %>% html_nodes("a") %>% html_text(trim = T)
