#Use libraries:
library(ggplot2)
library(tidyr)
library(dplyr)

#Read in table
data = read.csv("filec", sep="\t")

# Create the percentage identity brackets down to "0-5"
data <- data %>%
  mutate(perc_bracket = case_when(
    perc_identity >= 95 ~ "95-100",
    perc_identity >= 90 ~ "90-95",
    perc_identity >= 85 ~ "85-90",
    perc_identity >= 80 ~ "80-85",
    perc_identity >= 75 ~ "75-80",
    perc_identity >= 70 ~ "70-75",
    perc_identity >= 65 ~ "65-70",
    perc_identity >= 60 ~ "60-65",
    perc_identity >= 55 ~ "55-60",
    perc_identity >= 50 ~ "50-55",
    TRUE ~ "0-50"
  ))

#Change name 
colnames(data)[colnames(data) == "X.syntenic_blocks"] <- "blocks"

# Summarize the data
summary_table <- data %>%
  group_by(perc_bracket) %>%
  summarise(
    avg_syntenic_blocks = mean(`blocks`, na.rm = TRUE),
    avg_translocation_junctions = mean(Translocation_junctions, na.rm = TRUE),
    avg_inversion_junctions = mean(Inversion_junctions, na.rm = TRUE),
    avg_same_direction_duplication_junctions = mean(Same_direction_duplication_junctions, na.rm = TRUE),
    avg_loop_direction_duplication_junctions = mean(Loop_direction_duplication_junctions, na.rm = TRUE)
  ) %>%
  mutate(starting_value = case_when(
    perc_bracket == "100" ~ 100,
    perc_bracket == "95-100" ~ 95,
    perc_bracket == "90-95" ~ 90,
    perc_bracket == "85-90" ~ 85,
    perc_bracket == "80-85" ~ 80,
    perc_bracket == "75-80" ~ 75,
    perc_bracket == "70-75" ~ 70,
    perc_bracket == "65-70" ~ 65,
    perc_bracket == "60-65" ~ 60,
    perc_bracket == "55-60" ~ 55,
    perc_bracket == "50-55" ~ 50,
    TRUE ~ 0
  ))

#Convert output to data frame
summary_table_df <- as.data.frame(summary_table)

# Reshape the summary_table_df to long format and calculate proportions
summary_long <- summary_table_df %>%
    pivot_longer(cols = c(avg_translocation_junctions, avg_inversion_junctions, avg_same_direction_duplication_junctions, avg_loop_direction_duplication_junctions), 
                 names_to = "junction_type", 
                 values_to = "value") %>%
    group_by(perc_bracket) %>%
    mutate(proportion = value / sum(value) * 100) %>%  # Normalize to 100% proportions
    ungroup()

# Create pie charts for each percentage bracket
pie = ggplot(summary_long, aes(x = "", y = proportion, fill = junction_type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    facet_wrap(~ perc_bracket) +
    labs(title = "Proportions of Junctions by Percentage Bracket",
         x = NULL,
         y = NULL) +
    theme_void() +  # Removes axes and background
    theme(legend.position = "right")  # Adjust legend position

pdf("Chart_of_break_types.pdf", width=9, height=9)
pie
dev.off()