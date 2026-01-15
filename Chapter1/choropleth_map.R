library(readxl)
library(dplyr)
library(ggplot2)
library(sf)
library(tigris)

file_path <- "ohio_county_data.xlsx"

data <- read_excel(file_path)

data <- data %>%
  mutate(county_clean = tolower(gsub("\\s+county$", "", County)))

oh_counties <- counties(state = "OH", cb = TRUE, year = 2022) %>%
  st_as_sf() %>%
  mutate(county_clean = tolower(NAME))

map_data <- oh_counties %>%
  left_join(data, by = "county_clean")

add_title <- TRUE

find_col <- function(pattern) {
  hits <- grep(pattern, names(map_data), ignore.case = TRUE, value = TRUE)
  if (length(hits) == 0) {
    stop(paste0("No column matched pattern: ", pattern))
  }
  hits[1]
}

make_map <- function(var_name, title_text, file_name, fill_label) {
  p <- ggplot(map_data) +
    geom_sf(aes(fill = .data[[var_name]]), color = "white", size = 0.2) +
    scale_fill_viridis_c(option = "C", na.value = "gray90") +
    labs(fill = fill_label) +
    theme_minimal()

  if (add_title) {
    p <- p + labs(title = title_text)
  }

  ggsave(file_name, plot = p, width = 8, height = 6, units = "in")
  p
}

col_unemployment <- find_col("unemployment")
col_pm25 <- find_col("pm2\\.5|pm25")
col_some_college <- find_col("some\\s+college")

make_map(col_unemployment, "Unemployment by County (Ohio)", "unemployment_by_county.pdf", "Unemployment")
make_map(col_pm25, "Average Daily PM2.5 by County (Ohio)", "pm25_by_county.pdf", "Average Daily PM2.5")
make_map(col_some_college, "% Some College by County (Ohio)", "some_college_by_county.pdf", "% Some College")
