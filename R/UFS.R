#Will Booker
#3/26/2018
#Unsupervised Feature Selection Package

#------------------------------------------------------------------------------#
# Helper Functions
#------------------------------------------------------------------------------#

#Converts a df/matrix/vector to a vector
convert_to_vector = function(xs, name=NULL) {
  if(!is.null(name) && name %in% names(xs)) {xs = xs$var_explained}
  return(unlist(xs))}


#Gets the running sum of a vector
cumulate_vector = function(xs) {
  running = 0; cumulative = c()
  for(x in xs) {
    running = running + x
    cumulative = c(cumulative, running)}
  return(cumulative)
}


#Returns the trace of a matrix
tr = function(mat) {return(sum(diag(mat)))}


#Creating rejection vectors from a start vector and a matrix (MxN & MxP)
calc_rejection_vectors = function(q, b) {
  return(q - t(t(q) %*% b %*% t(b) / (norm(b, type="2") ^ 2)))}


#Calculating projection ammountss
calc_projection = function(q, b) {
  return(t(t(q) %*% b / (norm(b, type="2"))))}


#Computes a matrix's displacement vectors
calc_displacement_vectors = function(df) {
  d = data.frame(matrix(nrow=dim(df)[1], ncol=dim(df)[2]))
  for(i in 1:dim(df)[2]) {
    d[,i] = df[,i] - mean(df[,i])}
  return(data.matrix(d))
}


#Normalize a set of vectors (columns of dataframe)
normalize_vectors = function(df) {
  magnitudes = apply(df, 1, norm, "2")
  df = df / magnitudes
  df[is.na(df)] = 0
  return(df)}

#------------------------------------------------------------------------------#
# Model Functions
#------------------------------------------------------------------------------#

#Getting the distribution of the unfixed values of a partially fixed random distribution
var_given_fixed = function(df, fixed_indices) {
  
  #Getting the covariance of the dataframe
  sig = cov(df)
  
  #Getting variances and covariances
  sig11 = sig[-fixed_indices, -fixed_indices]
  sig22 = sig[fixed_indices, fixed_indices]
  sig12 = sig[-fixed_indices, fixed_indices]
  sig21 = sig[fixed_indices, -fixed_indices]
  
  #Getting new variance matrix
  return(sig11 - sig12 %*% solve(sig22) %*% sig21)
}


#Runs all the analysis for PCA
model_pca = function(df, dimensions=2, use_disp_vectors=TRUE) {
  
  #Getting the displacement vectors
  df = cbind(data.matrix(df))
  if(use_disp_vectors) {disp = calc_displacement_vectors(df)}
  else{disp = df}
  
  #Getting maximum magnitudes (for plotting)
  max_magnitudes = diag(t(df) %*% df) / (dim(df)[1] - 1)
  
  #Getting the eigen vectors and values of the covariance matrix
  e = eigen(cov(df))
  
  #Getting the first two principle components
  ps = e$vectors[,1:dimensions]
  
  #Translating the matrix using the first two principle components
  ys = data.matrix(df) %*% data.matrix(ps)
  
  #Getting variance explained
  var_explained = (e$values[1:dimensions] / sum(e$values))
  
  #Getting PCA rejection vectors and magnitudes
  rejection = t(calc_rejection_vectors(t(disp), ps))
  magnitudes = diag(t(rejection) %*% rejection) / (dim(df)[1] - 1)
  
  #Returning the matrix
  return(list(df=df, rejection=rejection, magnitudes=magnitudes, max_magnitudes=max_magnitudes,
              transformed=ys, transformation_matrix=t(ps), var_explained=var_explained))
}


#Running GVM over a set of data
model_gvm = function(df, num_vectors, use_disp_vectors=TRUE) {
  
  #Getting the displacement vectors
  df = data.matrix(df)
  saved = cbind(df)
  if(use_disp_vectors) {disp = calc_displacement_vectors(df)}
  else{disp = df}
  
  #Getting the maximum magnitudes of the dataframe (for plotting)
  max_magnitudes = diag(t(df) %*% df) / (dim(df)[1] - 1)
  
  #Initializing the indices data value
  total_var = tr(cov(df))
  indices = c()
  var_remaining = c()
  
  #Feature selecting
  rejection = cbind(disp)
  if(num_vectors > 0) {
    for(i in 1:num_vectors) {
      max = which.max(apply(rejection, 2, norm, "2"))[1]
      indices = c(indices, max)
      rejection = calc_rejection_vectors(rejection, rejection[,max])
      var_remaining = c(var_remaining, tr(var_given_fixed(df, indices)) / total_var)
    }
    vectors = data.matrix(df[,indices])
  }
  else{vectors = NULL}
  
  #Getting magnitudes and selected feature vectors
  #magnitudes = apply(rejection, 2, norm, "2")
  magnitudes = diag(t(rejection) %*% rejection) / (dim(df)[1] - 1)
  
  #Getting the variance explained
  total_explained = 1 - var_remaining
  var_explained = c(total_explained[1])
  if(num_vectors > 1) { 
    for(i in 2:length(total_explained)) {
      var_explained = c(var_explained, total_explained[i] - total_explained[i-1])}}
  
  #Returning the indices of the most important vectors
  return(list(df=saved, rejection=rejection, magnitudes=magnitudes, max_magnitudes=max_magnitudes, 
              vectors=vectors, indices=indices, var_explained=var_explained))
}

#------------------------------------------------------------------------------#
# View Functions
#------------------------------------------------------------------------------#

#Plotting function from the PCA model
view_pca = function(model, use_disp_vectors=TRUE, scale=9, ...) {
  
  #Getting variables from model
  df = model$df
  mag = model$magnitudes / model$max_magnitudes
  vects = model$df %*% t(model$transformation_matrix)
  vects = scale * t(normalize_vectors(t(vects)))
  
  #Plotting the color-coded variance from PCA
  plot_reduced(t(df), mag=mag, vects=vects, normalize=FALSE, ...)
}


#Plotting the GVM data
view_gvm = function(gvm_output, use_disp_vectors=TRUE, ...) {
  
  #Extracting GVM data
  max_magnitudes = gvm_output$max_magnitudes
  magnitudes = gvm_output$magnitudes
  indices = gvm_output$indices
  df = gvm_output$df
  
  #Plotting GVM data
  out = plot_reduced(t(df), mag=magnitudes/max_magnitudes, idxs=indices, normalize=FALSE, ...)
}

#------------------------------------------------------------------------------#
# Final Functions
#------------------------------------------------------------------------------#

#' Convert Color
#' @description Casts a set of magnitudes to a color gradient.
#' @param mag The magnitudes to convert
#' @param gradient The gradient to use for coloring. "~" reverses gradient direction. Defaults to "RdYlGn".
#' @param normalize Normalizes the magnitudes to 1-9.9 before color casting. Defaults to TRUE.
#' @return Returns a vector containing the color magnitudes. Has a side-effect that changes the color palette.
#' @keywords color, casting
#' @export
#' @examples
#' convert_color(party)
convert_color = function(mag=NULL, gradient="RdYlGn", normalize=TRUE) {
  
  #If there are magnitudes to plot
  colors = "black"
  if(!is.null(mag)) {
    
    #Data manipulation
    mag = unlist(mag)
    
    #Normalizing magnitudes
    if(normalize) {
      range = max(mag) - min(mag)
      mag = mag - min(mag)
      mag = mag / range
    }
    
    #Importing color library
    library(RColorBrewer)
    
    #Getting color direction
    if(substr(gradient, 1, 1) == "~") {
      reverse = TRUE; gradient=substring(gradient,2)}
    else{reverse = FALSE}
    
    #Chosing colors
    palette_choice = brewer.pal(n=9, name=gradient)
    if(reverse){palette_choice = rev(palette_choice)}
    palette(palette_choice)
    colors = 1 + mag * 8.9
  }
  
  #Returning the color values
  return(as.vector(unlist(colors)))
}


#' Plot Reduced
#' @description Dimensionally reduces the columns of a dataset, plots the results, and colors by an input magnitude.
#' @param mag The color scale of each of the vectors. 
#' @param idxs Specific vectors to color from the dataframe.
#' @param vects Specific vectors to color from outside the dataframe.
#' @param normalize Automatically normalizes the magnitudes to a range from 0-1. Defaults to TRUE.
#' @param gradient Color gradient to use. Using a "~" reverses the direction. Defaults to "YlOrRd". 
#' @param color Color to use for plotting points. Defaults to "forestgreen".
#' @keywords plot, reduced
#' @export
#' @examples 
#' plot_reduced(senate, mag=party, gradient="~RdBu", main="Senators: PCA Plot")
#' plot_reduced(t(senate), mag=apply(t(senate) %*% diag(as.vector(t(party))), 1, sum), gradient="~RdBu", main="Bills By Party")
#' plot_reduced(t(senate), mag=apply(senate, 2, sum), gradient="RdYlGn", main="Bills By Popularity")
#' plot_reduced(t(senate), main="Senate Bills: PCA Plot")
plot_reduced = function(df, mag=NULL, idxs=NULL, vects=NULL, reduce=model_pca, 
                        normalize=TRUE, gradient="YlOrRd", color="forestgreen", ...) {
  
  #Getting the colors
  colors = convert_color(mag, gradient, normalize)
  
  #Running the reduction analysis
  pcomp = reduce(df, 2)
  
  #Plotting all vectors
  plot(pcomp$transformed, pch=19, col=colors, xlab=NA, ylab=NA, ...)
  
  #Plotting hit indices
  if(!is.null(idxs)) {
    hits = t(pcomp$transformation_matrix %*% t(df[c(idxs),,drop=FALSE]))
    points(hits, pch=19, col=color)
  }
  
  #Plotting hit vectors
  if(!is.null(vects)) {
    hits = t(pcomp$transformation_matrix %*% data.matrix(vects))
    points(hits, pch=19, col=color)
  }
}


#' Scree Plot
#' @description Plots a scree plot of one or two vectors
#' @param xs The first series to plot. 
#' @param ys An additional second series to plot. 
#' @param cumulate Plots a cumulative scree plot. Defaults to false.
#' @keywords scree, plot
#' @export
#' @examples
#' scree_plot(pca(senate, 10, FALSE), gvm(senate, 10, FALSE), TRUE)
#' scree_plot(pca(t(senate), 10, FALSE), gvm(t(senate), 10, FALSE), TRUE)
scree_plot = function(xs, ys=NULL, cumulate=FALSE) {
  
  #Normalizing inputs
  title = "Scree Plot"
  xs = convert_to_vector(xs, "var_explained")
  ys = convert_to_vector(ys, "var_explained")
  
  #Cumulating inputs
  if(cumulate) {
    title = "Cumulative Scree Plot"
    xs = c(0, cumulate_vector(xs))
    if(!is.null(ys)) {ys = c(0, cumulate_vector(ys))}}
  
  #Plotting
  plot(xs, pch=19, ylim=c(0,1), col="black", type="p", 
       main=title, xlab="Components", ylab="Variance Explained")
  lines(xs, col="black")
  lines(ys, pch=19, col="red", type="p")
  lines(ys, col="red")
}


#' Principle Component Analysis
#' @description Runs PCA and plots the variance explained by a set of column vectors.
#' @param df The dataframe to plot.
#' @param n The number of dimensions to keep.
#' @param show Displays the plot. Defaults to TRUE.
#' @return \item{transformation_matrix}{The transformation matrix.}
#' @return \item{var_explained}{The explained variance by each component.}
#' @keywords pca
#' @export
#' @examples
#' pca(senate, 2, main="Bill PCA Variance Explained")
#' pca(t(senate), 2, main="Senator PCA Variance Explained")
pca = function(df, n, show=TRUE, ...) {
  
  #Setting input vectors to 0
  if(n < 0) {n=0}
  
  #Getting outputs from the model
  model = model_pca(df, dimensions=n)
  
  #Showing the plot if required
  if(show) {view_pca(model, ...)}
  
  #Returning the required outputs
  return(list(transformed=model$transformed, transformation_matrix=model$transformation_matrix, 
              var_explained=model$var_explained))
}


#' Greedy Variance Maximization
#' @description Runs GVM and plots the variance explained by a set of column vectors.
#' @param df The dataframe to plot.
#' @param num_vectors The number of vectors to keep
#' @param show Displays the plot. Defaults to TRUE.
#' @return \item{indices}{The indices of the most important vectors.}
#' @return \item{var_explained}{The variance explained by each of the vectors.}
#' @keywords gvm
#' @export
#' @examples 
#' gvm(senate, 2, main="Bill GVM Variance Explained")
#' gvm(t(senate), 2, main="Senator GVM Variance Explained")
gvm = function(df, n, show=TRUE, ...) {
  
  #Setting input vectors to 0
  use_disp = TRUE
  if(n < 0) {n=0; use_disp=FALSE}
  
  #Getting outputs from the model
  model = model_gvm(df, num_vectors=n, use_disp_vectors=use_disp)
  
  #Showing the plot if required
  if(show) {view_gvm(model, ...)}
  
  #Returning the required outputs
  return(list(transformed=model$vectors, indices=model$indices, 
              var_explained=model$var_explained))
}


#' Nearest Vectors
#' @description Computes the nearest vectors to a starting vector and plots the display.
#' @param df The dataframe to plot.
#' @param ix Specific index vector to search from the dataframe.
#' @param vect Specific vectors to search from outside the dataframe.
#' @param n Number of neighbors to return. Defaults to 20.
#' @param show Shows the plot. Defaults to TRUE.
#' @return Returns the indices of the n closest vectors.
#' @keywords nearest, vectors
#' @export
#' @examples
#' nearest_vectors(senate, idx=1, n=5, main="Nearest Vectors to Index 1")
#' nearest_vectors(t(senate), vect=t(senate[1,]), n=5, main="Nearest Vectors to Index 1")
nearest_vectors = function(df, idx=NULL, vect=NULL, n=20, show=TRUE, ...) {
  
  #Variable conversion
  vector = vect
  index = idx
  
  #Getting percent rejection magnitude and plotting the maximum projection
  if(!is.null(index)) {vector = df[,index]}
  
  #Getting percent rejection magnitude and plotting the maximum projection
  if(!is.null(vector)) {
    percent_projection = abs(calc_projection(df, vector)) / apply(df, 2, norm, "2")
    if(show) {
      plot_reduced(t(df), mag=percent_projection, vects=vector, normalize=FALSE, gradient="Blues", color="darkorange", ...)}
    
  #Returning the nearest vectors
  sort_order = rev(order(percent_projection))
  return_data = cbind(sort_order, percent_projection[sort_order])
  if(is.null(n)) {n = dim(df)[2]}
  return(return_data[1:n,])}
}


# TODO: Allow for arbitrary vector inputs
#' Variance Explained
#' @description Calculates the variance explained by a set of vectors.
#' @param df The dataframe to plot.
#' @param idxs Specific indices to search from the dataframe.
#' @param show Shows the plot. Defaults to TRUE.
#' @keywords variance, explained
#' @export
#' @examples
#' variance_explained(senate, c(1,2,3), main="Variance Explained By [1,2,3]")
#' variance_explained(t(senate), c(1,2,3), main="Variance Explained By [1,2,3]")
variance_explained = function(df, idxs, show=TRUE, ...) {
  
  #Converting variable names
  df = t(df)
  show_plot = show
  
  #Converting function inputs
  indices = idxs
  
  #If the number of vectors is less than zero, revert symmetry breaking
  df = data.matrix(t(df))
  
  #Getting the maximum magnitudes (using the original df to lessen low variance terms)
  max_magnitudes = diag(t(df) %*% df) / (dim(df)[1] - 1)
  
  #Getting remaining magnitudes
  remaining_variance = var_given_fixed(df, indices)
  magnitudes = diag(remaining_variance)
  
  #Inserting dropped elements
  corrected_magnitudes = c()
  counter = 1
  for(i in 1:dim(df)[2]) {
    if(i %in% indices) {corrected_magnitudes = c(corrected_magnitudes, 0)}
    else {
      corrected_magnitudes = c(corrected_magnitudes, magnitudes[counter])
      counter = counter + 1}}
  magnitudes = corrected_magnitudes
  
  #Normalizing magnitudes
  normalized_magnitudes = magnitudes / max_magnitudes
  
  #Plotting normalized magnitudes
  plot_reduced(t(df), mag=normalized_magnitudes, idxs=indices, normalize=FALSE, ...)
  
  #Remaining Variance
  remaining_variance = sum(magnitudes) / sum(tr(cov(df)))
  
  #Returning the remaining variance
  return(1 - remaining_variance)
}

#------------------------------------------------------------------------------#
# 
#------------------------------------------------------------------------------#
