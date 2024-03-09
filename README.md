# Lagrangian
An R package that implements transport covariance
  functions belonging to the class of Lagrangian Spatio-Temporal
  Covariance Function

---

### License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>

---

## Sample R codes using the DiffOp package

```
library(Lagrangian)
library(dplyr)
library(fields)

x <- seq(0, 1, length.out = 10)
y <- seq(0, 1, length.out = 10)
loc2d <- expand.grid(x, y) %>% as.matrix()

t <- seq(0, 10, by = 1)
loc3d <- cbind(rep(loc2d[, 1], each = length(t)),
               rep(loc2d[, 2], each = length(t)), t)

cov_mat <- cov_uni_lagrangian(location = loc3d, sigma2 = 1, scale = 0.5, nu = 1, 
                              vx = 0.1001, vy = 0.1001)

cov_mat <- cov_uni_lagrangian(location = loc3d, sigma2 = 1, scale = 1, nu = 1, 
                              vx = 0.1001, vy = 0.1001, nonstat = 1) #deformation

cov_mat <- cov_uni_lagrangian(location = loc3d, sigma2 = 1, scale = 1, nu = 1, 
                              vx = 0.1001, vy = 0.1001, nonstat = 2) #spatially varying

plot(cov_mat[1,1:100])

chol(cov_mat)

ind <- which(loc3d[, 3] == 0)
quilt.plot(loc3d[ind, 1], loc3d[ind, 2], cov_mat[1, ind], nx = 10, ny = 10)

ind <- which(loc3d[, 3] == 3)
quilt.plot(loc3d[ind, 1], loc3d[ind, 2], cov_mat[1, ind], nx = 10, ny = 10)


## SIMULATION

library(MASS)

set.seed(1235)
Z <- mvrnorm(1, mu = rep(0, ncol(cov_mat)), Sigma = cov_mat)

quilt.plot(loc3d[ind, 1], loc3d[ind, 2], Z[ind], nx = length(x), ny = length(y), xlab = "Longitude", ylab = "", zlim = range(Z), xaxt = 'n')

```
