############################
# Import libraries 
############################

library(vegan)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggbeeswarm)
library(gridExtra)


############################
# Set working directory
############################

setwd("D:/D-Holocene T2024/data")


############################
# File names and data loading
############################

fillname <- "Liang et al., 2024 Lacutrine brGDGT and BGC brGDGTs202407.xlsx"

GDGTdata <- readxl::read_excel(fillname, sheet = "Table S1. ss_brGDGTs")
GDGTdata <- na.omit(GDGTdata)
# Print column names and unique predicted labels
print(colnames(GDGTdata))

print(unique(GDGTdata$Predicted_Labels))

# Subset GDGT data
GDGTdata_com <- GDGTdata[, c(3:17, 27)]
colnames(GDGTdata_com)

############################
# Transform data to long format
############################

GDGTdata_long <- pivot_longer(
  GDGTdata_com,
  cols = -Predicted_Labels,
  names_to = "variable",
  values_to = "value"
)
head(GDGTdata_long)

############################
# Calculate summary statistics
############################

GDGTdata_summary <- GDGTdata_long %>%
  group_by(variable, Predicted_Labels) %>%
  summarize(
    mean_value = mean(value),
    se_value = sd(value) / sqrt(n())
  )

# Recode Predicted_Labels for readability
GDGTdata_summary$Predicted_Labels <- factor(
  GDGTdata_summary$Predicted_Labels,
  levels = c("C1","C2"),
  labels = c("Cluster 1", "Cluster 2")
)
head(GDGTdata_summary)


############################
# Generate bar plot for GDGT data summary
############################
colors <- c("#68a9cb", "#fb6e55")
brGDGT_com_plot <- ggplot(GDGTdata_summary,
                          aes(x = variable, y = mean_value, fill = Predicted_Labels)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "(A) Barplot of brGDGTs", x = "", y = "Relative abundance (%)") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.9), legend.direction = "vertical",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(fill = NULL)

print(brGDGT_com_plot)

############################
# Load and inspect environmental data
############################

envdf <- readxl::read_excel(fillname, sheet = "Table S2. ss_env")
colnames(envdf)
envdata=envdf[,c(9,10,12,18,23,24,25)]
colnames(envdata)[2:5]=c("Salinity","DO bottom", "LST","MAF")

########################################################
# Convert wide data to long format for environmental data
########################################################

envdata_long <- pivot_longer(
  envdata[,],
  cols = -Predicted_Labels,
  names_to = "variable",
  values_to = "value"
)
head(envdata_long)


selected_data_long <- envdata_long %>%
  filter(variable %in% c("DO bottom", "Salinity", "pH")) %>%
  mutate(variable = factor(variable, levels = c("Salinity", "pH", "DO bottom")))
colnames(selected_data_long)
tail(selected_data_long)
selected_data_long$Predicted_Labels<- factor(selected_data_long$Predicted_Labels,
                                             levels = c("C1","C2"),
                                             labels = c("Cluster 1", "Cluster 2"))
####################################################################################
# Generate violin and box plots for selected environmental data
####################################################################################


plot <- ggplot(selected_data_long, aes(x = Predicted_Labels, y = value, 
                                       fill = Predicted_Labels)) +
              geom_boxplot(notch = F,  color = "black", lwd = 0.2, alpha = 0.7, show.legend = FALSE) +
               geom_violin(data = selected_data_long, alpha = 0.5, position = position_dodge(width = 0.75),size = 1, color = NA, trim = FALSE) +
               ggbeeswarm::geom_quasirandom(shape = 21, size = 2, dodge.width = 0.75, color = "black", alpha = 0.5, show.legend = FALSE) +
    labs(title = "(B) Violin plot of environmental variables", x = "Cluster", y = "Value") +
  theme_bw() +
 scale_fill_manual(values = colors) +
  facet_wrap(~variable, scales = "free_y") +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(plot)


#####################################
## Perform Kruskal-Wallis test ##
#####################################


perform_kruskal_test <- function(data, group_col_name) {
  results <- lapply(data[,-which(names(data) == group_col_name)], function(var) {
    kruskal.test(var ~ data[[group_col_name]])
  })
  df_results <- data.frame(
    Environment = names(results),
    p_value = sapply(results, function(x) x$p.value),
    chi_squared = sapply(results, function(x) x$statistic),
    df = sapply(results, function(x) x$parameter),
    significance = sapply(results, function(x) ifelse(x$p.value < 0.001, "***", ifelse(x$p.value < 0.01, "**",
                                                                                       ifelse(x$p.value < 0.05, "*", ifelse(x$p.value < 0.1, ".", "")))))
  )
  df_results[order(nchar(as.character(df_results$significance)), decreasing = TRUE), ]
}

# brGDGT DATA to the Kruskal-Wallis test function
df_results_GDGTdata <- perform_kruskal_test(GDGTdata_com, "Predicted_Labels")
print(df_results_GDGTdata)

df_results_envdata <- perform_kruskal_test(envdata, "Predicted_Labels")
print(df_results_envdata)

############################
# SEM ANLYSIS #
############################
library(nlme)
library(piecewiseSEM)
library(car)
library(lmtest)
library(semTools)
library(bestNormalize)


dataset=data.frame(envdata,GDGTdata)
colnames(dataset)

model1 <- gls(MBT.5me ~ MAF+Salinity+pH+LDA.score+LST+DO.bottom, data = dataset)
model2 <- gls(IR ~MAF+Salinity+pH+LDA.score+LST+DO.bottom, data = dataset)
model3 <- gls(LDA.score ~ MAF+Salinity+pH+LST+DO.bottom, data = dataset)
model4 <- gls(Salinity ~MAF+LST, data = dataset)
model5 <-gls(pH ~  MAF+LST, data = dataset)
model6 <- gls(DO.bottom ~MAF+Salinity+LST, data = dataset)
model7 <- gls(LST~ MAF, data = dataset)


sem_model_gls <- psem(model1, model2, model3,model4,model5,model6,model7,
                      MBT.5me%~~%IR,pH%~~%Salinity)
LLchisq(sem_model_gls)

AIC(sem_model_gls)
AIC(sem_model_gls, AIC.type = "dsep")
#get the basis set
basisSet(sem_model_gls)

summary(sem_model_gls,standardize="scale")
plot(sem_model_gls,title="gls SEM")

summary_model1 <- summary(sem_model_gls)
coefficients_model1 <- summary_model1$coefficients


significant_coefficients <- coefficients_model1[coefficients_model1[, "P.Value"] <= 0.1, c("Response","Predictor", "Std.Estimate","P.Value")]



print(significant_coefficients)


significant_coefficients <- coefficients_model1[coefficients_model1[, "P.Value"] <= 0.1, c("Response","Predictor", "Std.Estimate","P.Value")]

# 在显著系数中添加星号
significant_coefficients$Significance <- ifelse(significant_coefficients$P.Value < 0.001, "***",
                                                ifelse(significant_coefficients$P.Value < 0.01, "**",
                                                       ifelse(significant_coefficients$P.Value < 0.05, "*",
                                                              ifelse(significant_coefficients$P.Value < 0.1, ".", ""))))


significant_coefficients$Std.Estimate <- sprintf("%.2f", significant_coefficients$Std.Estimate)

print(significant_coefficients)


# Clean up environment before session ends
rm(list = ls())
