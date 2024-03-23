corr_fetal <- table(fetal$timepoint, fetal$EK_PB_annov1)
df <- as.data.frame.matrix(corr_fetal)

# Calculate percentages row-wise
df_perc <- prop.table(as.matrix(df), 1) * 100

# Convert to long format for ggplot
df_long <- pivot_longer(as.data.frame(df_perc), cols = everything(), names_to = "Cell_Type", values_to = "Percentage")
df_long$Timepoint <- rep(rownames(df_perc), times = ncol(df_perc))
# Add a column for timepoint extracted from rownames
df_long$Timepoint <- factor(df_long$Timepoint, levels = unique(df_long$Timepoint))

df_long$Timepoint <- factor(df_long$Timepoint, levels = c('Day 59','Week 9','Day 80','Day 82','Week 11','Week 12','Week 13','Week 14','Day 105','Week 15','Week 16','Day 125','Week 17','Week 18','Week 19','Week 20','Week 22','Week 24','Week 27'))
# Create dot plot with ggplot2
ggplot(df_long, aes(x = Timepoint, y = Percentage, group = Cell_Type)) +
  geom_point(aes(color = Cell_Type), size = 3) +  # Plot dots
  geom_line(aes(linetype = Cell_Type)) +    # Connect dots for each cell type with lines
  theme_minimal() +
  labs(title = "Percentage of Cell Types Over Time",
       x = "Timepoint",
       y = "Percentage",
       color = "Cell Type",
       linetype = "Cell Type")


df <- table(fetal$timepoint, fetal$EK_PB_annov1)
df <- as.data.frame(df)
df$Timepoint <- df$Var1
head(df)
df_summarized <- df %>%
  mutate(Timepoint = case_when(
    Timepoint %in% c('Day 59', 'Week 9', 'Day 80', 'Day 82', 'Week 11') ~ 'Week 8-11',
    Timepoint %in% c('Week 12', 'Week 13', 'Week 14', 'Day 105', 'Week 15') ~ 'Week 12-15',
    Timepoint %in% c('Week 16', 'Day 125', 'Week 17', 'Week 18', 'Week 19') ~ 'Week 16-19',
    Timepoint %in% c('Week 20', 'Week 22') ~ 'Week 20-22',
    Timepoint %in% c('Week 24', 'Week 27') ~ 'Week 24-27',
    TRUE ~ Timepoint
  ))
head(df_summarized)
df_summarized$Cell_Type <- df_summarized$Var2
df_summarized$Count <- df_summarized$Freq
df_grouped <- df_summarized %>%
  group_by(Timepoint, Cell_Type) %>%
  summarise(Count = sum(Count)) %>%
  ungroup()

# Calculate percentages row-wise
df_grouped <- df_grouped %>%
  group_by(Timepoint) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  ungroup()

# Convert Timepoint to a factor to ensure correct ordering
df_grouped$Timepoint <- factor(df_grouped$Timepoint, levels = c('Week 8-11', 'Week 12-15', 'Week 16-19', 'Week 20-22', 'Week 24-27'))

ggplot(df_grouped, aes(x = Timepoint, y = Percentage, group = Cell_Type)) +
  geom_point(aes(color = Cell_Type), size = 3) +  # Plot dots
  geom_line(aes(color = Cell_Type)) +    # Connect dots for each cell type with lines
  theme_minimal() +
  labs(title = "Percentage of Cell Types Over Summarized Timepoints",
       x = "Timepoint",
       y = "Percentage",
       color = "Cell Type",
       linetype = "Cell Type")

ggplot(df_grouped, aes(x = Timepoint, y = Percentage, group = Cell_Type)) +
  geom_point(aes(color = Cell_Type), size = 3) +  # Plot dots
  geom_line(aes(color = Cell_Type)) +    # Connect dots for each cell type with lines
  theme_minimal() +
  labs(title = "Percentage of Cell Types Over Summarized Timepoints",
       x = "Timepoint",
       y = "Percentage",
       color = "Cell Type",
       linetype = "Cell Type") + facet_wrap(~Cell_Type)
