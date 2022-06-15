path <- "C:\\Users\\dango\\OneDrive\\Робочий стіл\\StatMagData\\CompStat\\csdata"
paths.to.tables <- list.files(path, full.names = TRUE)

list.of.tables <- lapply(paths.to.tables,
                         function(f.tbl) {
                           tbl <- read.csv(
                             file = f.tbl, header = F
                           )[,c(1,6)]
                           colnames(tbl) <- f1c(
                             "date", unlist(strsplit(f.tbl, "[_.]"))[2]
                           )
                           tbl
                         })

merged.table <- Reduce(function(x, y) {merge(x, y, by = "date")},
                       list.of.tables)

#write.csv(merged.table, paste(path, "merged.csv", sep="\\"))
