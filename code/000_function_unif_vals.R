## Function to generate user-defined number of columns that use index column
# to generate new values with a certain percentage from the uniform distribution
uniform_values <- function(df, column_name, num_columns, percent_var) {
  
  df_new <- df 
  
  for (i in 1:num_columns) {
    new_col <- apply(df, 1, function(row) {
      original_value <- row[column_name]
      
      if (is.na(original_value)) {
        return(NA)
        
      } else {
        
       return(original_value + runif(1, min = -percent_var * abs(original_value), max = percent_var * abs(original_value)))
      }
    })
    
    df_new[[paste0(column_name, "_mod_", i)]] <- new_col
    
  }
  
  return(df_new)
  
}





