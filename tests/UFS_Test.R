
#------------------------------------------------------------------------------#
# Importing UFS
#------------------------------------------------------------------------------#

library("UFS")

#------------------------------------------------------------------------------#
# PCA Testing
#------------------------------------------------------------------------------#

#Political party test
out = plot_reduced(senate, mag=party, gradient="~RdBu")

#Party support test
votes_by_party = t(senate) %*% diag(as.vector(t(party)))
party_popularity = apply(votes_by_party, 1, sum)
out = plot_reduced(t(senate), mag=party_popularity, gradient="~RdBu")

#Bipartisan popularity test
bipartisan_popularity = apply(senate, 2, sum)
out = plot_reduced(t(senate), mag=bipartisan_popularity, gradient="RdYlGn")

#------------------------------------------------------------------------------#
# Nearest Vectors Test
#------------------------------------------------------------------------------#

#Getting the nearest vectors to a single vector
out = nearest_vectors(t(senate), ix=1)

#Getting the nearest vectors to the first principle component
pca_senate = pca(t(senate), 2)
pc1 = pca_senate$transformation_matrix[1,] * -10
out = nearest_vectors(senate, vect=pc1)

#------------------------------------------------------------------------------#
# GVM Test
#------------------------------------------------------------------------------#

#Running on senate
ufs = gvm(senate, -1)
ufs = gvm(senate, 0)
ufs = gvm(senate, 2)
ufs = gvm(senate, 5)

#Running on bills
ufs = gvm(t(senate), -1)
ufs = gvm(t(senate), 0)
ufs = gvm(t(senate), 2)
ufs = gvm(t(senate), 5)

#------------------------------------------------------------------------------#
# PCA Test
#------------------------------------------------------------------------------#

#Running on senate
ufs = pca(senate, 1)
ufs = pca(senate, 2)
ufs = pca(senate, 5)

#Running on bills
ufs = pca(t(senate), 1)
ufs = pca(t(senate), 2)
ufs = pca(t(senate), 5)

#------------------------------------------------------------------------------#
# Var Explained Test
#------------------------------------------------------------------------------#

gvm(senate, n=5, show=FALSE)$var_explained
pca(senate, n=5, show=FALSE)$var_explained

gvm(t(senate), n=5, show=FALSE)$var_explained
pca(t(senate), n=5, show=FALSE)$var_explained

#------------------------------------------------------------------------------#
# Nearest Vectors Test
#------------------------------------------------------------------------------#

#Where is the vect vector?
#t(nearest_vectors(t(senate), vect=vect[,1]*100, n=5))
#t(nearest_vectors(t(senate), vect=vect[,1]/15.7, n=5))

#------------------------------------------------------------------------------#
# 
#------------------------------------------------------------------------------#
