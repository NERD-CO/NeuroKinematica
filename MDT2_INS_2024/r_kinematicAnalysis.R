# devtools::install_github("PabRod/kinematics")


library(kinematics)

mov <- data.frame(t = c(0, 1, 2, 3, 4), 
                  x = c(0, 1, 2, 3, 4), 
                  y = c(0, 1, 4, 9, 16))

plot(mov$x, mov$y, xlab = "x", ylab = "y")


mov_analyzed <- append_dynamics(mov)

# vx, vy and aspeed: horizontal, vertical and absolute speed.
# ax, ay, and aaccel: horizontal, vertical and absolute acceleration.
# curv and curv_radius: curvature and curvature radius.
# disp_x, disp_y and adisp: horizontal, vertical and absolute displacement 
# .....(since previous time step).

# Speeds’ units are: x unit divided by t unit.
# Acceleration’s units are: x unit divided by t unit squared.
# Curvature unit is 1 divided by x unit.
# Curvature radius’ unit is the same as x.
# Displacements’ units are the same as x.

library(ggplot2)
ggplot(data = mov_analyzed, 
       mapping = aes(x = x, y = y, col = curv_radius, size = aspeed)) +
  geom_point() +
  scale_color_gradient(low="blue", high="red")



#### A second example

# Generate the times
ts <- seq(0, 10, by = 0.05)

# Calculate the positions using a function
xs <- ts * cos(ts)
ys <- ts * sin(ts)

# Store as data frame
mov <- data.frame(t = ts, 
                  x = xs, 
                  y = ys)


plot(mov$x, mov$y, xlab = "x", ylab = "y", asp = 1)


mov_analyzed <- append_dynamics(mov)

ggplot(data = mov_analyzed, 
       mapping = aes(x = x, y = y, col = curv_radius, size = aspeed)) +
  geom_point(alpha = 0.5) +
  coord_fixed() +
  scale_color_gradient(low="blue", high="red")



## 3rd example

mov <- kinematics::example_mov

plot(mov$x, mov$y, xlab = "x", ylab = "y", asp = 1)

mov_analyzed <- append_dynamics(mov)

ggplot(data = mov_analyzed, 
       mapping = aes(x = x, y = y, col = aaccel, size = aspeed)) +
  geom_point(alpha = 0.1) +
  coord_fixed() +
  scale_color_gradient(low="blue", high="red")


hist(mov_analyzed$aaccel, 
     breaks = 500, 
     xlab = 'Accelerations', 
     main = 'Acceleration histogram')

## in matlab save T, X, Y for each movement into CSV
### LOAD in 



