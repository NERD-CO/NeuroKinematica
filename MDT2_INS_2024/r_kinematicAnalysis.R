devtools::install_github("PabRod/kinematics")


library(kinematics)
library(ggplot2)

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

#library(ggplot2)
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


## CD to Computer location
setwd("E:\\Dropbox\\PowerPoint_Meta\\2024_INS_Vancouver\\Data\\finalResults")
csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty list to store data frames
list_of_dataframes <- list()

# Loop through the .csv files and read each one
for (file_name in csv_files) {
  # Read the CSV file into a data frame
  data_frame <- read.csv(file_name)
  
  # DO SOMETHING 
  mov_analyzed <- append_dynamics(data_frame)
  
  # Modify the file name to include '_EN' before the '.csv'
  new_file_name <- sub(".csv$", "_EN.csv", basename(file_name))
  
  # Save the modified data frame back to disk with the new file name
  write.csv(mov_analyzed, file.path(dirname(file_name), new_file_name), row.names = FALSE)
  
  # Add the data frame to the list using the filename as the key
  list_of_dataframes[[file_name]] <- mov_analyzed
}



ggplot(data = mov_analyzed, 
       mapping = aes(x = x, y = y, col = aaccel, size = aspeed)) +
  geom_point(alpha = 0.1) +
  coord_fixed() +
  scale_color_gradient(low="blue", high="red")


hist(mov_analyzed$aaccel, 
     breaks = 500, 
     xlab = 'Accelerations', 
     main = 'Acceleration histogram')

# Access a specific data frame by its file name
specific_df <- list_of_dataframes[["path/to/your/file.csv"]]




