## Function to interpolate y value for given x values
interpolate_y <- function(x_values, tam_table) {
  ### Initialize an empty vector for interpolated y values
  interpolated_y_values <- numeric(length(x_values))
  
  ### Pre-compute x and y vectors and slopes from tam_table
  x_col <- as.numeric(unlist(tam_table[, 1]))
  y_col <- as.numeric(unlist(tam_table[, 2]))
  slopes <- diff(y_col) / diff(x_col)
  
  for (i in seq_along(x_values)) {
    x_value <- x_values[i]
    
    #### Check if x_value is outside the range of tam_table$x
    if (is.na(x_value) || x_value < min(x_col) || x_value > max(x_col)) {
      interpolated_y_values[i] <- 0  # Set interpolated y value to zero
    } else {
      ##### Find the interval where x_value lies
      idx <- findInterval(x_value, x_col)
      
      ##### Linear interpolation
      y1 <- y_col[idx]
      slope <- slopes[idx]
      interpolated_y_values[i] <- y1 + slope * (x_value - x_col[idx])
    }
  }
  
  return(interpolated_y_values)
}
