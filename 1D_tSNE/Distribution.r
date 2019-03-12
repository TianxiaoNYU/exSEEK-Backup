
library(Rtsne)
library(ggplot2)
library(dplyr)

library(FNN)

predata <- read.table("filter.scimpute_count.Norm_RLE.Batch_null.domains_combined.txt", header = T, row.names = 1, sep = "\t")
postdata <- read.table("filter.scimpute_count.Norm_RLE.Batch_RUV.domains_combined.txt", header = T, row.names = 1, sep = "\t")

batch_info <- read.table("~/Bioinfos/data/scirep_batch.txt", sep = ",", header = T, row.names = 1, stringsAsFactors = T)
batch_info$RNA.Isolation.batch <- factor(batch_info$RNA.Isolation.batch)
batch_temp <- batch_info
batch_temp$names <- rownames(batch_info)
batch_temp <- arrange(batch_temp, names)
rownames(batch_temp) <- batch_temp$names
batch_info <- as.data.frame(batch_temp[,-ncol(batch_temp)])

sample_class <- as.factor(c(rep("Colorectal Cancer", times = 100), rep("Normal", times = 50), 
                                          rep("Normal", times = 6), rep("Prostate Cancer", times = 36)))

#batch_info[which(batch_info$RNA.Isolation.batch == 7), ]

temp <- as.data.frame(t(predata))
temp$names <- rownames(temp)
temp <- arrange(temp, names)
rownames(temp) <- temp$names
predata <- as.data.frame(temp[,-ncol(temp)])

temp <- as.data.frame(t(postdata))
temp$names <- rownames(temp)
temp <- arrange(temp, names)
rownames(temp) <- temp$names
postdata <- as.data.frame(temp[,-ncol(temp)])



p <- Rtsne(predata, dims = 1)
q <- Rtsne(postdata, dims = 1)

pre_p <- cbind(batch_info$RNA.Isolation.batch, p$Y)
post_q <- cbind(batch_info$RNA.Isolation.batch, q$Y)
sum_data <- as.data.frame(rbind(pre_p, post_q))
foo <- c(rep("Before Batch Removal", times = length(p$Y)), rep("After Batch Removal", times = length(q$Y)))
sum_data$Type <- sample_class
sum_data$State <- foo

names(sum_data) <- c("Batch", "t_SNE", "Type", "State")
sum_data$Batch <- as.factor(sum_data$Batch)
sum_data$Type <- as.factor(sum_data$Type)
sum_data$State <- as.factor(sum_data$State)

sum_data

s <- ggplot(data = sum_data, aes(x = t_SNE, col = Batch))  + geom_density() + facet_wrap(~State)
s <- s + scale_fill_brewer(palette="Set2") + 
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black")) 
pdf("test_for_batch.pdf", 10, 6)
s
dev.off()

s <- ggplot(data = sum_data, aes(x = t_SNE, col = Type))  + geom_density() + facet_wrap(~State)
s <- s + scale_fill_brewer(palette="Set2") + 
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black")) 
pdf("test_for_type.pdf", 10, 6)
s
dev.off()



p <- Rtsne(postdata, dims = 1)
tSNE_Vis_before <- as.data.frame(cbind(batch_info$RNA.Isolation.batch, p$Y))
names(tSNE_Vis_before) <- c("Batch", "Y")
s <- ggplot(data = tSNE_Vis_before, aes(x = Y, col = Batch)) + geom_density() 
s <- s + scale_fill_brewer(palette="Set2") + 
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black")) 
png("before.png")
s
dev.off()

p <- Rtsne(predata, dims = 1)
tSNE_Vis_after <- as.data.frame(cbind(batch_info$RNA.Isolation.batch, p$Y))
names(tSNE_Vis_after) <- c("Batch", "Y")
tSNE_Vis_after$Batch <- factor(tSNE_Vis_after$Batch)
s <- ggplot(data = tSNE_Vis_after, aes(x = Y, col = Batch)) + geom_density()
s <- s + scale_fill_brewer(palette="Set2") + 
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black")) 
png("after.png")
s
dev.off()

## Test for the cancer/normal -- after
p <- Rtsne(postdata, dims = 1)
tSNE_Vis_after <- as.data.frame(cbind(sample_class, p$Y))
names(tSNE_Vis_after) <- c("Type", "Y")
s <- ggplot(data = tSNE_Vis_after, aes(x = Y, fill = Type)) + geom_density()
s <- s + scale_fill_brewer(palette="Set2") + 
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black")) 
png("after_test.png")
s
dev.off()

p <- Rtsne(predata, dims = 1)
tSNE_Vis_before <- as.data.frame(cbind(sample_class, p$Y))
names(tSNE_Vis_before) <- c("Type", "Y")
s <- ggplot(data = tSNE_Vis_before, aes(x = Y, fill = Type)) + geom_density()
s <- s + scale_fill_brewer(palette="Set2") + 
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black")) 
png("before_test.png")
s
dev.off()

t1 <- tSNE_Vis_after$Y[which(tSNE_Vis_after$Type == 1)]
t0 <- tSNE_Vis_after$Y[which(tSNE_Vis_after$Type == 0)]
t2 <- tSNE_Vis_after$Y[which(tSNE_Vis_after$Type == 2)]

b1 <- tSNE_Vis_before$Y[which(tSNE_Vis_after$Type == 1)]
b0 <- tSNE_Vis_before$Y[which(tSNE_Vis_after$Type == 0)]
b2 <- tSNE_Vis_before$Y[which(tSNE_Vis_after$Type == 2)]

KL.divergence(t1, t0, k = 1)
KL.divergence(b1, b0, k = 1)




