## Function to interpolate y value for given x values
interpolate_y <- function(x_values, tam_table) {
  ### Initialize an empty vector for interpolated y values
  interpolated_y_values <- numeric(length(x_values))
  
  for (i in seq_along(x_values)) {
    x_value <- x_values[i]
    
    #### Check if x_value is outside the range of tam_table$x
    if (x_value < min(as.numeric(unlist(tam_table[,1]))) || 
        x_value > max(as.numeric(unlist(tam_table[,1]))) || 
        is.na(x_value) == TRUE ) {
      interpolated_y_values[i] <- 0  # Set interpolated y value to zero
    } else {
      ##### Find the interval where x_value lies
      idx <- findInterval(x_value, as.numeric(unlist(tam_table[,1])))
      
      ##### Linear interpolation
      y1 <- as.numeric(unlist(tam_table[,2]))[idx]
      y2 <- as.numeric(unlist(tam_table[,2]))[idx + 1]
      slope <- slopes[idx]
      interpolated_y_values[i] <- y1 + slope * (x_value - as.numeric(unlist(tam_table[,1]))[idx])
    }
  }
  
  return(interpolated_y_values)
}
