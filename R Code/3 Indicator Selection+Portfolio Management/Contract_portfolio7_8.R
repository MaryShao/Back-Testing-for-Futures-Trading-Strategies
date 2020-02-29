setwd("C:/Users/m8sha/Desktop/DATA")
perf <- read.csv("performance.txt", sep="", stringsAsFactors=FALSE)

perf.order <- c()
for (i in (1:nrow(perf))){
  ord = sort(as.numeric(perf[i,]), decreasing = TRUE, index.return = TRUE)$ix
  tmp = colnames(perf)[ord]
  perf.order = rbind(perf.order,tmp)
}

rownames(perf.order) = rownames(perf)
colnames(perf.order) = c(1:96)

setwd("C:/Users/m8sha/Desktop")
write.table(perf.order,'perf.order.txt')