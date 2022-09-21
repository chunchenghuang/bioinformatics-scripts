#Installation
#install.packages("ggplot2")
#install.packages("hrbrthemes")
#https://ggplot2.tidyverse.org/reference/scale_manual.html
#https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html

# Library
library(ggplot2)
library(hrbrthemes)

# Dummy data
x <- LETTERS[1:20]
y <- paste0("gene", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- sample(0:5, 400, replace=T, prob=c(0.8,0.04,0.04,0.04,0.04,0.04))
col=c("grey","lightblue","darkblue","orange","red","purple")
ggplot(data, aes(X, Y)) + 
  geom_tile(aes(fill=factor(data$Z)),colour="white") +
  scale_fill_manual(
    values=col,
    breaks=c("0","1","2","3","4","5"),
    labels=c("none","missense","duplication","deletion","substitution","insertion"))

