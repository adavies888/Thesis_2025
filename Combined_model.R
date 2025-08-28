# Combined LASSO

# Libraries
library(ggplot2)
library(stats)
library(tidyverse)
library(readxl)
library(factoextra)
library(ggfortify)
library(ggpubr)
library(rstatix)
library(finalfit)
library(mlbench)
library(DALEX)
library(DALEXtra)
library(modelStudio)
library(recipes)
library(naniar)
library(mice)
library(VIM)
library(glmnet)
library(caret)
library(pROC)
library(data.table)
library(tableone)
library(ggfortify)
library(devtools)
library(themis)
library(pROC)
library(tidymodels)
library(tictoc)
library(yardstick)
library(vip)

# Load data
bile_acids <- read_csv("Bile_acid_reduced.csv")
mg_data <- read_csv("Metagenomics_reduced.csv")
mt_data <- read_csv("Metatranscriptomics_reduced.csv")

# remove column
bile_acids <- bile_acids[, -1]
mg_data <- mg_data[, -1]
mt_data <- mt_data[, -1]

# add the sample ID back in to perform merge
bile_acids$ID <- paste(bile_acids$Number, bile_acids$Letter)
mg_data$ID <- paste(mg_data$Number, mg_data$Letter)
mt_data$ID <- paste(mt_data$Number, mt_data$Letter)

# Merge the dfs
data1 <- merge(bile_acids, mg_data, on = ID)
data_final <- merge(data1, mt_data, on = ID)

# remove ID column
data_final <- data_final[, -3]


### Logistic regression model with LASSO features
test <- data_final %>%    
  filter(Number %in% c(10,14))

train <- data_final %>%
  setdiff(test)

test <- test %>% dplyr::select(-Number) 
train <- train %>% dplyr::select(-Number) 

lasso_factors <- c(A = 0, C = 1) 
train$Letter <- lasso_factors[train$Letter]
test$Letter <- lasso_factors[test$Letter]

x_train <- model.matrix(Letter~., train)[, -1]
x_test <- model.matrix(Letter~., test)[,-1] 

y_train <- train$Letter
y_test <- test$Letter

# Create the model
lasso_mod <- glmnet(x_train, 
                    y_train,
                    alpha = 1, 
                    family = "binomial", standardize = TRUE) # Fit lasso model on training data

plot(lasso_mod, xvar = "lambda")  # Draw plot of coefficients

# Prepare training control for LOOCV
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE)

# Convert x to data frame and combine with y
data <- as.data.frame(x_train)
data <- as.data.frame(scale(data))
data$y <- factor(y_train)  
levels(data$y) <- c("Class0", "Class1")  

# Train the model using caret
set.seed(123)
cv.out <- train(
  y ~ ., 
  data = data,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl,
  tuneGrid = expand.grid(
    alpha = 1,                        # LASSO
    lambda = 10^seq(-4, 1, length = 50)  # custom lambda grid
  ))

print(cv.out)
plot(cv.out)

bestlam <- cv.out$bestTune$lambda

lasso_coef = predict(lasso_mod, type = "coefficients", s = bestlam) # Display coefficients using lambda chosen by CV

myCoefs <- coef(lasso_mod, s=bestlam);
myCoefs[which(myCoefs != 0 ) ]              

myCoefs@Dimnames[[1]][which(myCoefs != 0 ) ] 

myResults <- data.frame(
  features = myCoefs@Dimnames[[1]][ which(myCoefs != 0 ) ], #intercept included
  coefs    = myCoefs              [ which(myCoefs != 0 ) ]  #intercept included
)

myResults$features <- gsub("`|'", "", myResults$features)

myResults <- myResults[-1, ]
myResults$sign <- ifelse(myResults$coefs >= 0, "Positive", "Negative")
myResults$features <- factor(myResults$features, levels = myResults$features[order(abs(myResults$coefs), decreasing = FALSE)])

# Plot the features
options(repr.plot.width = 100, repr.plot.height = 3)
ggplot(myResults, aes(x = features, y = abs(coefs), fill = sign)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = signif(coefs, 3)), 
            hjust = -0.1, 
            size = 3.5) +
  coord_flip(clip = "off") + 
  labs(x = NULL, y = "Weight", title = "Weight of LASSO features") +
  scale_fill_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F8766D")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Predict with the logistic regression model 
lasso_pred = predict(lasso_mod, s = bestlam, newx = x_test, type = "response")

lasso_pred

pred <- as.factor(ifelse(lasso_pred >= 0.5, 1, 0))

conf_matrix <- confusionMatrix(as.factor(y_test), pred)

# Print confusion matrix
print(conf_matrix)

# AUC
roc_obj <- roc(as.numeric(y_test), as.numeric(pred))
auc(roc_obj)


data_list <- colnames(data_final)[-c(1, 2)]

data_long <- data_final %>%
  pivot_longer(
    cols = -c(Number, Letter),
    values_to = "Value",
    names_to = "Feature"
  )


# T-test to see change between week 0 and 4
data_long %>%
  group_by(Letter, Feature) %>%
  get_summary_stats(Value, type = "mean_sd")

pvals <- c()

for (ba in data_list) {
  test_df <- data_long %>%
    filter(Feature %in% ba)
  
  res <- t.test(Value ~ Letter, data = test_df)
  
  pvals[ba] <- res$p.value
}

pvals_adj <- p.adjust(pvals, method = "fdr", n = length(pvals))

data_table <- data.frame(
  Feature = names(pvals),
  P_Value = as.numeric(pvals),
  P_Adj = as.numeric(pvals_adj)
)

# Optional
data_table <- data_table[order(data_table$P_Adj), ]

# View the table
print(data_table)


lasso_features <- as.vector(myResults$features)

data_scaled <- data_final %>% 
  mutate(across(-c(1, 2), scale))

data_long <- data_scaled %>%
  pivot_longer(
    cols = -c(Number, Letter),
    values_to = "Value",
    names_to = "Feature"
  )

data_features <- data_long  %>%
  filter(Feature %in% lasso_features)

letter_to_time <- c(A = "Baseline", C = "Week 4") 
data_features$Time <- letter_to_time[data_features$Letter]

label_positions <- data_features %>%
  group_by(Feature) %>%
  summarise(y_pos = max(Value, na.rm = TRUE) * 1.05)

annotation_plot_df <- data_table %>%
  merge(label_positions, by = "Feature")

bxp <- ggboxplot(
  data_features, x = "Time", y = "Value",
  color = "Time", palette = "jco",
  facet.by = "Feature", short.panel.labs = FALSE
)
bxp + 
  facet_wrap(~Feature) + 
  labs(
    x = "Baseline vs Week 4", 
    y = "Feature distribution") +
  geom_text(
    data = annotation_plot_df,
    aes(x = 2.2, y = 4, label = paste("p-value:", round(P_Adj, 3))),  
    inherit.aes = FALSE,
    size = 4,
    hjust = 0.5
  )

annotation_plot_df