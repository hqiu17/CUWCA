#!/usr/bin/env Rscript
# Wthin corresponding analysis

library(ade4)
library(MASS)
options(max.print = 9999999)

# Parse command line argument
cmd_args = commandArgs(T);
infile = cmd_args[1];

# Read codon data and transpose it
table<-read.table(cmd_args[1], header=T)
table<-as.data.frame(t(as.matrix(table)))

# Make a factor from synonymous codons
aa=c(rep("F", 2), rep("L", 6), rep("S", 6),
      rep("Y", 2), rep("C", 2), rep("P", 4),
      rep("H", 2), rep("Q", 2), rep("R", 6),
      rep("I", 3), rep("T", 4), rep("N", 2),
      rep("K", 2), rep("V", 4), rep("A", 4),
      rep("D", 2), rep("E", 2), rep("G", 4)
     )
aafac=as.factor(aa)

# Correspondence analysis
table.ca<-dudi.coa(table, scan=FALSE, nf=4)

# Within-group correspondence analysis
table.wt<-wca(table.ca, aafac, scan=FALSE, nf=3)

# List the coordinates of each gene on the first 4 axis
outfile = paste(infile, "gene_coordinates.txt", sep=".");
sink(outfile)   
table.wt$co
sink()

# Percentage of total variance accounted for by each axis
outfile = paste (infile, "inertia.txt", sep=".");
sink(outfile)   
inertia.dudi(table.wt, row.inertia=FALSE, col.inertia=FALSE)
sink()

# Plot the genes along the first and second axis
outfile = paste (infile, ".scat_comp1-2.pdf", sep="");
pdf(file=outfile, width=9, height=9, pointsize=16)
plot(table.wt$co$Comp1, table.wt$co$Comp2, pch=21, col="black")  
abline (h=seq(-0.5,0.5,0.05), v=seq(-0.5,0.5,0.05), col="grey")
abline (h=seq(-0.5,0.5,0.1),  v=seq(-0.5,0.5,0.1),  col="blue")
dev.off()

# Plot the genes along the second  and third axis
outfile = paste (infile, ".scat_comp2-3.pdf", sep="");
pdf(file=outfile, width=9, height=9, pointsize=16)
plot(table.wt$co$Comp2, table.wt$co$Comp3, pch=21, col="black")
abline (h=seq(-0.5,0.5,0.05),v=seq(-0.5,0.5,0.05), col="grey")
abline (h=seq(-0.5,0.5,0.1), v=seq(-0.5,0.5,0.1),  col="blue")
dev.off()

# Make histogram of the data along the first axes
outfile = paste (infile, ".his_comp1.pdf", sep="")
pdf(file=outfile, width=9, height=9, pointsize=25)
plot(hist(table.wt$co$Comp1, breaks=40), col="blue")  
dev.off()

# Make histogram of the data along the second axes
outfile = paste (infile, ".his_comp2.pdf", sep="")
pdf(file=outfile, width=9, height=9, pointsize=25)
plot(hist(table.wt$co$Comp2, breaks=40), col="blue")
dev.off()

# Contour plot of the data in 2D
outfile = paste (infile, ".kde_comp1-2.pdf", sep="")
pdf(file=outfile, width=9, height=9, pointsize=25)
f1=kde2d(table.wt$co$Comp1, table.wt$co$Comp2,n=25)
contour(f1)
abline(h=seq(-0.5,0.5,0.05), v=seq(-0.5,0.5,0.05), col="grey")
abline(h=seq(-0.5,0.5,0.1),  v=seq(-0.5,0.5,0.1),  col="blue")
dev.off()

