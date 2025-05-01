# Simon P. Sherman III PhD 
# April 21 2025

# PCA 
# Note this PCA is based on a combination of the two tutorials
# source 1 : https://www.r-bloggers.com/2021/05/principal-component-analysis-pca-in-r/#google_vignette
# source 2 : https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#google_vignette

#load packages
library(factoextra)
library(devtools)
library(usethis)
library(psych) # for Correlation and Orthogonal plots
library(ggrepel)
library(ggbiplot)
library(corrr)
library(ggcorrplot)
library(FactoMineR)

#use RStudio Session Directory

df <- read.csv("morph_knowns_all.csv") #knowns active variables
df.X <- read.csv("morph_all_pca.csv") #supplementary individuals
# note later PCAs will observe the supplementary variables by a new categorical
# as.factor column called 'Group' which is the geological groups
summary(df)

df$PPK <- as.factor(df$PPK)

str(df)
# be careful some morphometrics will convert from num to character. Inspect the data!

# if PPK or Group column is not a factor and rather a character or levels are empty
# you know it is the wrong spreadsheet used

#IMPORTANT : remove categorical columns so PCA can work
row.names(df) <- df$name
df <- df[, -1]

#color palette if more than 3 groups like tutorial in finnstats
color6 <- Polychrome::palette36.colors(n = 4) # n needs to match # of groups/PPks

pairs.panels(df[,-20], # not that 20 is rep of # of columns removing name and PPK
             gap = 0,
             bg = color6[as.numeric(as.factor(df$PPK))],
             pch=21)

# Note the very last column will be a categorical variable need to temporarily remove
#### known <- df[, 1:19] # if below doesnt work

#PCA pre normalizing for multicollinearity
pc <- prcomp(df[,-20],
             center = TRUE,
             scale. = TRUE)
attributes(pc)

pc$center #center

pc$scale #scale

print(pc)

summary(pc)

# PCA: Orthogonality of PCs normalizing multicollinearity
pairs.panels(pc$x,
             gap=0,
             bg = color6[as.numeric(as.factor(df$PPK))],
             pch=21)

# First Biplot to see what it looks like
library(devtools)
library(ggbiplot)
g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = df$PPK, #or Group
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68) # can increase to 95% confidence
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)

# NOTE: Do the same steps for df.X (unknowns)
df.X <- read.csv("morph_all_pca.csv") #supplementary individuals

summary(df.X)

df.X$PPK <- as.factor(df.X$PPK)

str(df.X)
# be careful some morphometrics will convert from num to character. Inspect the data!

#IMPORTANT : remove categorical columns so PCA can work
row.names(df.X) <- df.X$name
df.X <- df.X[, -1]

known <- df[, 1:19]

# Use the same variable names used in PCA
new_data <- df.X[, colnames(known)]

#Scale the unknowns using the known PCAs center and scale
new_scaled <- scale(new_data,
                    center = pc$center,
                    scale = pc$scale)

# Project the unknowns into PCA space
new_proj <- as.data.frame(new_scaled %*% pc$rotation)

# new PCA biplot
library(ggrepel)

# Create label names for the new points
unknown_labels <- rownames(df.X)

# Plot PCA with original data
g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = df$PPK,
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal',
        legend.position = 'top') +
  geom_point(data = new_proj, aes(x = PC1, y = PC2),
             shape = 17, color = "black", size = 1) +
  geom_text_repel(data = new_proj, aes(x = PC1, y = PC2, label = unknown_labels),
                  size = 3, color = "black")

print(g)

# To see what fits within the ellipsoids
library(car)
library(dplyr)
library(sp)

# Create a data frame for knowns in PCA space
known_scores <- as.data.frame(pc$x[, 1:2])  # PC1 and PC2
known_scores$Group <- df$PPK

# Unknowns are in new_proj (already in PC1/PC2 space)
unknown_scores <- new_proj[, 1:2]
unknown_scores$Label <- rownames(new_proj)

# Function to test point-in-ellipse for each group
point_in_ellipse <- function(unknowns, knowns_group, level = 0.68) {
  ellipse_pts <- data.frame(ellipse(center = colMeans(knowns_group[, 1:2]),
                                    shape = cov(knowns_group[, 1:2]),
                                    radius = sqrt(qchisq(level, df = 2))))
  
  inside <- apply(unknowns[, 1:2], 1, function(row) {
    sp::point.in.polygon(row[1], row[2], ellipse_pts$x, ellipse_pts$y) == 1
  })
  
  return(inside)
}

# build matric rows = unknowns 
group_names <- unique(known_scores$Group)
assignment_matrix <- sapply(group_names, function(grp) {
  grp_data <- known_scores %>% filter(Group == grp)
  point_in_ellipse(unknowns = unknown_scores, knowns_group = grp_data)
})

assigned_group <- apply(assignment_matrix, 1, function(x) {
  if (sum(x) == 1) {
    return(group_names[which(x)])
  } else {
    return(NA)  # ambiguous or outside all ellipses
  }
})

ellipse_table <- data.frame(
  Unknown_Label = rownames(unknown_scores),
  Assigned_Group = assigned_group
) %>% filter(!is.na(Assigned_Group))

# Create a proportional summary table
prop_table <- ellipse_table %>%
  group_by(Assigned_Group) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / nrow(unknown_scores))

#write .csv
write.csv(prop_table, file = "prop_table_SAA_2025.csv", row.names = TRUE)
write.csv(ellipse_table, file = "ellipse_table_predictions_SAA_2025.csv", row.names = TRUE)

#PCA Scores .csv
write.csv(known_scores, file = "Known_PCA_Scores_SAA_2025.csv", row.names = TRUE)
write.csv(unknown_scores, file = "Unknown_PCA_Scores_SAA_2025.csv", row.names = TRUE)

#For overlapping PCA prediciton ellipses to tell similarities in type by knowns compared
write.csv(assignment_matrix, file = "Assingment_Matrix_PCA_SAA_2025.csv", row.names = TRUE)

### PCA Results for variables

# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- res.pca$rotation
sdev <- res.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord[, 1:4])
