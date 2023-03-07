library(httr) # get requests
library(jsonlite) # convert data types

# results categories that we can search
search_results <- httr::GET("https://www.ebi.ac.uk/ena/portal/api/results?dataPortal=ena&format=json")
search_results_content <- jsonlite::fromJSON(rawToChar(search_results$content))
View(search_results_content) # we want 'sequence'

# fields within those categories we can search
search_fields <- httr::GET("https://www.ebi.ac.uk/ena/portal/api/searchFields?result=sequence&format=json")
search_fields_content <- jsonlite::fromJSON(rawToChar(search_fields$content))
View(search_fields_content) # we want 'host'

# the host species we're interested in
species_list <- c("Picea engelmannii",
                  "Dryas octopetala")

# formatted for html encoding
species_list_formatted <- species_list
for (i in 1:length(species_list_formatted)) {
  species_list_formatted[[i]] <- sub(" ", "%20", species_list_formatted[[i]])
}

# place holder for the returned results
df <- data.frame(matrix(nrow = 0, ncol = 3))

# looping through the API query per species
for (i in 1:length(species_list_formatted)) {
  query <- httr::GET(paste0("https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&query=host=%22",species_list_formatted[[i]],'%22&format=json'))
  query_content <- jsonlite::fromJSON(rawToChar(query$content))
  query_content$host <- species_list[[i]]
  query_content <- query_content[, c(3, 1, 2)]
  df <- rbind(df, query_content)
  Sys.sleep(sample(1:5, 1))
}
