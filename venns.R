install.packages("VennDiagram")
library(VennDiagram)

# Reference four-set diagram with my numbers
venn.plot <- draw.quad.venn(
area1 = 276,
area2 = 103,
area3 = 72,
area4 = 64,
n12 = 103,
n13 = 70,
n14 = 63,
n23 = 53,
n24 = 50,
n34 = 41,
n123 = 53,
n124 = 50,
n134 = 41,
n234 = 37,
n1234 = 37,
category = c("Binder", "High Affinity", "SMALI", "Scansite"),
fill = c("orange", "red", "green", "blue"),
lty = "solid",
cex = 2,
cat.cex = 2,
cat.col = c("orange", "red", "green", "blue")
);

venn.plot <- draw.triple.venn(
area1 = 103,
area2 = 72,
area3 = 64,
n12 = 53,
n23 = 41,
n13 = 50,
n123 = 37,
category = c("High Affinity", "SMALI", "Scansite"),
fill = c("orange", "green", "blue"),
lty = "blank",
cex = 2,
cat.cex = 2,
cat.col = c("orange", "green", "blue")
);

# Reference four-set diagram with high affinity SH2 numbers
venn.plot <- draw.quad.venn(
area1 = 39,
area2 = 28,
area3 = 20,
area4 = 15,
n12 = 14,
n13 = 14,
n14 = 9,
n23 = 4,
n24 = 4,
n34 = 1,
n123 = 3,
n124 = 3,
n134 = 1,
n234 = 0,
n1234 = 0,
category = c("ABL1", "GRB2", "NCK1", "SHP2"),
fill = c("orange", "red", "green", "blue"),
lty = "blank",
cex = 2,
cat.cex = 2,
cat.col = c("orange", "red", "green", "blue")
);

# Reference four-set diagram showing binary overlap of SH2 affinity numbers
venn.plot <- draw.quad.venn(
area1 = 39,
area2 = 28,
area3 = 20,
area4 = 15,
n12 = 14,
n13 = 14,
n14 = 9,
n23 = 4,
n24 = 4,
n34 = 1,
n123 = 3,
n124 = 3,
n134 = 1,
n234 = 0,
n1234 = 0,
category = c("ABL1", "GRB2", "NCK1", "SHP2"),
fill = c("orange", "red", "green", "blue"),
lty = "blank",
cex = 2,
cat.cex = 2,
cat.col = c("orange", "red", "green", "blue")
);
