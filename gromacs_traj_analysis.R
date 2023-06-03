library(bio3d)
# locate working directory
> getwd()
#change the directory
> setwd("path to directory")

#prepare dcd file on vmd
#or
mdconvert md_0_1_noPBC.xtc -t md_0_1.gro -o traj.dcd

# Read trajectory file
> dcd <- read.dcd("traj.dcd")

# Read the pdb file
> pdb <- read.pdb(wt.pdb)

#checking rows and columns in pdb and dcd
> dim(dcd)
> ncol(trj) == length(pdb$xyz)
[1] TRUE
# Trajectory Superposition on Calpha atoms
> ca.inds <- atom.select(pdb, elety = "CA")
> xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, 
 	               fixed.inds = ca.inds$xyz, 
	               mobile.inds = ca.inds$xyz)

# Principal Component Analysis
> pc <- pca.xyz(xyz[, ca.inds$xyz])
> plot(pc, col = bwr.colors(nrow(xyz)))

# Clustering
> hc <- hclust(dist(pc$z[, 1:2]))
> grps <- cutree(hc, k = 2)
> plot(pc, col = grps)

# DCCM
> cij <- dccm(xyz[, ca.inds$xyz])
> plot(cij)
