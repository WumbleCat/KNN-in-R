##############################################################################################################################
##############################################################################################################################
# Part I
##############################################################################################################################
##############################################################################################################################

# x and y are vectors (1 row, n columns)
dist_vect <- function(x, y) {
  distance <- sum(abs(x - y))
  return(distance)
}

# X is a matrix and y is a vector
# Basically applies dist_vect for every row of X and y to get a scalar value
# Then stores it in a vector
dist_mat <- function(X,y){
  V <- (apply(X, 1, dist_vect,y))
  return(V)
} 

# Vectorised version uses sweep function
dist_mat_fast <- function(X,y){
  D <- sweep(X, 2, y, "-")
  V <- rowSums(abs(D))
  return(V)
} 

# Demo test for dist_mat and dist_mat_fast
set.seed(414)
X <- matrix(rnorm(10, 3), 10, 3)
y <- rnorm(3)

dist_mat(X, y)
dist_mat_fast(X, y)


##############################################################################################################################
##############################################################################################################################
# Part II
##############################################################################################################################
##############################################################################################################################

# Set up variables: mainly
# n = number of rows
# X = matrix of 200x2 (original dataset)
# C = first 4 rows of X (cluster data)
set.seed(4184)
n <- 200
X <- matrix(rnorm(n*2, sample(c(-2, 2), n*2, replace = TRUE), 1), n, 2)

k <- 4
id_C <- 1:k
C <- X[id_C, ]

# find_nearest takes Matrix X (original dataset) and C (cluster center)
find_nearest <- function(X, C) {
  # cl contains information about the points in each cluster
  # It is a list of 4 vectors, where each index stands for the cluster number
  # and they contain the indices of its clusters
  cl <- vector("list", nrow(C)) 
  
  M <- X[1:nrow(X), ]
  
  # dists is a matrix (200x4) where each row represents the indices of matrix X 
  # and each column represents the distance from the center cluster to the index in X
  dists <- apply(C, 1, function(row) dist_mat(M, row))
  
  # row_mins gets the cluster center that has the closest distance
  row_mins <- apply(dists, 1, which.min)
  
  # indices contains the indices of the closest distance in X
  indices <- seq_along(row_mins)
  cl <- split(indices, row_mins)
  
  return(cl)
}


# Calculates the total distance of data points from their cluster centers
dist_tot <- function(X, C, cl) {
  dtot <- 0
  
  # For loop over each cluster
  for (j in seq_along(cl)) {
    # Indices of data points in cluster j
    cluster_indices <- cl[[j]]  
    # The center for cluster j
    cluster_center <- C[j, ]    
    
    # Sum distances of all points in this cluster to the cluster center
    for (i in cluster_indices) {
      dtot <- dtot + dist_vect(X[i, ], cluster_center)
    }
  }
  return(dtot)
}

plot_clustering <- function(X, C, cl) {
  # Define the limits of the plot
  x_limits <- c(min(X[, 1]), max(X[, 1]))
  y_limits <- c(min(X[, 2]), max(X[, 2]))
  
  # Set up the plot with defined limits and no points initially (type = "n")
  plot(X, type = "n", xlab = colnames(table)[1], ylab = colnames(table)[2], main = paste0(colnames(table)[1], "vs", colnames(table)[2]),
       asp = 1, xlim = x_limits, ylim = y_limits)
  
  # Assign unique colors to each cluster
  colors <- rainbow(length(cl))
  
  # Plot each cluster with a different color
  for (i in 1:length(cl)) {
    # Draw the cluster center first (as a larger point)
    points(C[i, 1], C[i, 2], col = colors[i], pch = 8, cex = 3)
    
    # Draw each point assigned to the cluster
    points(X[cl[[i]], 1], X[cl[[i]], 2], col = colors[i], pch = 19)
  }
}

cl <- find_nearest(X,C)
plot_clustering(X, C, cl)
png("./points.png", width = 500, height = 500)
#create the plot
plot_clustering(X, C, cl)
#close the file
dev.off()

##############################################################################################################################
##############################################################################################################################
# Part III
##############################################################################################################################
##############################################################################################################################
# How it basically works:
# FIND THE POINT CLOSEST TO THE CENTER IN EACH CLUSTER 
# WHICHEVER HAS THE LOWEST DISTANCE OUT OF EACH CLUSTER, IS THE NEW CENTER
# MAKE NEW MATRIX C WITH THE NEW CENTERS
# RERUN max_iter TIMES
clu_algo <- function(X, k, max_iter) {
  # Get random center clusters
  random_indices <- sample(nrow(X), k, replace = FALSE)
  C <- X[random_indices, ]
  
  # initialise the first C and cl
  C_min <- C
  cl_min <- find_nearest(X, C)
  cost_min <- Inf
  
  # Initialise values if costs have a convergence
  tol <- 1e-10
  converged <- FALSE
  
  # For loop to adjust the centers
  for (iter in 1:max_iter) {
    # check for each cluster
    for (i in 1:k) {
      # Get the points assigned to the current cluster
      cluster_indices <- cl_min[[i]]
      cluster_points <- X[cluster_indices, ]
      
      # Calculate the sum of distances from each point to all other points in the cluster
      distances <- apply(cluster_points, 1, function(point) sum(dist_mat_fast(cluster_points, point)))
      
      # Find the index of the point that is closest (minimum distance)
      min_idx <- which.min(distances)
      
      # Get the x,y coordinates of the center
      center <- cluster_points[min_idx, ]
      
      # Make a temporary copy to compare distances without changing the main C
      C_temp <- C_min
      
      # Set the row of the C matrix to the point that minimises the distances
      C_temp[i, ] <- center
      
      # Assign points to the new cluster centers
      cl_temp <- find_nearest(X, C_temp)
      
      # Calculate the cost
      cost <- dist_tot(X, C_temp, cl_temp)
      
      # Update if a better cluster assignment is found
      if (cost < cost_min) {
        C_min <- C_temp
        cl_min <- cl_temp
        # Check if cost has converged to find the best solution
        # Skip the first iteration
        if (k != 1){
          if (abs(cost_min-cost)<tol){
            converged <- TRUE
          }
        }
        cost_min <- cost
      }
      
      
    }
    if (converged == TRUE) {
      return(list(C_min, cl_min, cost_min))
    }
  }
  return(list(C_min, cl_min, cost_min))
}

final <- clu_algo(X,4,50)
C_min <- final[[1]]
cl_min <- final[[2]]
plot_clustering(X, C_min, cl_min)
png("./clust.png", width = 500, height = 500)
#create the plot
plot_clustering(X, C_min , cl_min)
#close the file
dev.off()

##############################################################################################################################
##############################################################################################################################
# Part IV
##############################################################################################################################
##############################################################################################################################

par(mfrow = c(2, 3))
data("iris")

# Create list of vectors that contains the matrix of each unique pair
my_list = vector("list", 0)
idx <- 1
for (i in 1:4){
  for (j in (i+1):4){
    # 
    a <- iris[i]
    b <- iris[j]
    M <- cbind(a,b)
    my_list[[idx]] <- M
    idx <- idx + 1
  }
}

for (i in 1:6){
  table <- my_list[[i]]
  final_list <- clu_algo(table,2,10)
  C_min <- final_list[[1]]
  cl_min <- final_list[[2]]
  cost_min <- final_list[[3]]
  plot_clustering(table, C_min, cl_min)
}

png("./iris.png", width = 500, height = 500)
#create the plot
for (i in 1:6){
  table <- my_list[[i]]
  final_list <- clu_algo(table,2,10)
  C_min <- final_list[[1]]
  cl_min <- final_list[[2]]
  cost_min <- final_list[[3]]
  plot_clustering(table, C_min, cl_min)
}
#close the file
dev.off()
