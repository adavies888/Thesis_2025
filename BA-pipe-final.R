# Analysis of bile acids

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
library(tidymodels)
library(tictoc)
library(yardstick)
library(vip)
library(verification)
library(dplyr)
library(Boruta)
library(mlbench)
library(randomForest)

# Load data
bile_acids <- read_excel("0_Bile acid metabolomics.xlsx")

bile_acids <- bile_acids %>%
  mutate(across(-SampleID, as.numeric))

bile_acids <- bile_acids %>% 
  extract(SampleID, into = c("Number", "Letter"), regex = "([0-9]+)([A-Z])") %>%
  mutate(Number = factor(Number, levels = sort(as.numeric(unique(Number)))))

# Reduce dataset to just weeks 0 and 4
bile_acids <- bile_acids %>%
  filter(Letter %in% c("A", "C")) %>%
  mutate(Letter = factor(Letter, levels = sort(unique(Letter))))

bile_acid_list <- colnames(bile_acids)[-c(1,2)] # a list of the bile acid features

letter_to_time <- c(A = "Baseline", C = "Week 4") 
bile_acids$Time <- letter_to_time[bile_acids$Letter]

# Pivot longer
bile_acids_long <- bile_acids %>%
  pivot_longer(
    cols = -c(Number, Letter, Time),
    values_to = "Value",
    names_to = "bile_acid"
  )

# Visualise the data
bxp <- ggboxplot(
  bile_acids_long, x = "Time", y = "Value",
  color = "Time", palette = "jco",
  facet.by = "bile_acid", short.panel.labs = FALSE
)
bxp + 
  facet_wrap(~bile_acid, scales = "free_y") +
  labs(
    x = "Baseline vs Week 4",
    y = "Recorded value"
  )

# PCA
pca.res <- prcomp(bile_acids[,-c(1,2,52)])

summary(pca.res)

loadings <- pca.res$rotation

pca_var <- pca.res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)

x_label <- paste0("PC1 (", pca_var_perc[1], "%)")
y_label <- paste0("PC2 (", pca_var_perc[2], "%)")

# Create a biplot to visualise the PCA
fviz_pca_biplot(pca.res,
                repel = TRUE,
                label = "none",
                invisible = "var",
                #col.ind = bile_acids$Time,
                habillage = bile_acids$Time,
                palette = c("#00AFBB", "#B6ACE7"), 
                addEllipses = TRUE, ellipse.level = 0.9) +
  labs(x = x_label, y = y_label,  title = "(a) Bile Acid: PC 1 against 2") + 
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18))


# T-test to see change between week 0 and 4
bile_acids_long %>%
  group_by(Letter, bile_acid) %>%
  get_summary_stats(Value, type = "mean_sd")

pvals <- c()

for (ba in bile_acid_list) {
  test_df <- bile_acids_long %>%
    filter(bile_acid %in% ba)
  
  res <- t.test(Value ~ Letter, data = test_df)
  
  pvals[ba] <- res$p.value
}

pvals_adj <- p.adjust(pvals, method = "fdr", n = length(pvals)) # adjust the p-values

bile_acid_table <- data.frame(
  bile_acid = names(pvals),
  P_Value = as.numeric(pvals),
  P_Adj = as.numeric(pvals_adj)
)

bile_acid_table <- bile_acid_table[order(bile_acid_table$P_Adj), ] # sort by p-value

print(bile_acid_table) # View the table

# Plot the t-test ones < 0.05
ba_ttest <- bile_acid_table %>%
  filter(P_Adj < 0.05) %>%
  pull(bile_acid) %>%
  as.list()

bile_acids_scaled <- bile_acids %>% 
  mutate(across(-c(1, 2, 52), scale))

bile_acids_long <- bile_acids_scaled %>%
  pivot_longer(
    cols = -c(Number, Letter, Time),
    values_to = "Value",
    names_to = "bile_acid"
  )

bile_acid_features <- bile_acids_long  %>%
  filter(bile_acid %in% ba_ttest)

label_positions <- bile_acid_features %>%
  group_by(bile_acid) %>%
  summarise(y_pos = max(Value, na.rm = TRUE) * 1.05)

annotation_plot_df <- bile_acid_table %>%
  merge(label_positions, by = "bile_acid")

bxp <- ggboxplot(
  bile_acid_features, x = "Time", y = "Value",
  color = "Time", palette = "jco",
  facet.by = "bile_acid", short.panel.labs = FALSE
)
bxp + 
  facet_wrap(~bile_acid) + 
  labs(
    x = "Baseline vs Week 4",
    y = "Bile Acid distribution" ) +
  geom_text(
    data = annotation_plot_df,
    aes(x = 2.2, y = 4, label = paste("p-value:", round(P_Adj, 3))),
    inherit.aes = FALSE,
    size = 4,
    hjust = 0.5
  )


# Logisic regression model
test <- bile_acids %>%
  dplyr::select(-Time) %>%
  filter(Number %in% c(1,3)) # Test train-split

train <- bile_acids %>%
  dplyr::select(-Time) %>% 
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

# Prepare training control for LOOCV
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE)

# Convert x to data frame and combine with y
data <- as.data.frame(x_train)
data <- as.data.frame(scale(data))
data$y <- factor(y_train)
levels(data$y) <- c("Class0", "Class1")  

set.seed(92)
lasso_mod <- glmnet(x_train, 
                    y_train,
                    alpha = 1, 
                    family = "binomial", standardize = TRUE)

# Train the model using caret
set.seed(123)
cv.out <- train(
  y ~ ., 
  data = data,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl,
  tuneGrid = expand.grid(
    alpha = 1,
    lambda = 10^seq(-4, 1, length = 50)
  ))

print(cv.out)
plot(cv.out)

bestlam <- cv.out$bestTune$lambda

lasso_coef = predict(lasso_mod, 
                     type = "coefficients", 
                     s = bestlam) # Display coefficients using lambda chosen by CV

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

lasso_features <- as.vector(myResults$features) # list of lasso features
lasso_features

# Plot the coefficients
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

# Now train the lasso model and predict
lasso_mod <- glmnet(x_train, 
                    y_train,
                    alpha = 1, 
                    family = "binomial", standardize = TRUE) # Fit lasso model on training data

plot(lasso_mod, xvar = "lambda")# label = TRUE)    # Draw plot of coefficients

lasso_pred = predict(lasso_mod, s = bestlam, newx = x_test, type = "response")

lasso_pred

pred <- as.factor(ifelse(lasso_pred >= 0.5, 1, 0))

conf_matrix <- confusionMatrix(as.factor(y_test), pred)

# Print confusion matrix
print(conf_matrix)

# AUC
roc_obj <- roc(as.numeric(y_test), as.numeric(lasso_pred))
auc(roc_obj)


# Random forest Model

# Boruta feature selection
train$Letter <- as.factor(train$Letter)
levels(train$Letter) <- make.names(levels(train$Letter))

colnames(train) <- make.names(colnames(train))

set.seed(111)
boruta <- Boruta(Letter ~ ., data = train, doTrace = 2, maxRuns = 500)

print(boruta)

plot(boruta, las = 2, cex.axis = 0.7)

bor <- TentativeRoughFix(boruta)
print(bor)
boruta_stats <- attStats(boruta)

# Extract important boruta features
confirmed_features <- rownames(boruta_stats[boruta_stats$decision == 'Confirmed', ])
print(confirmed_features)
attStats(boruta)


selected_formula <- getNonRejectedFormula(boruta)

train_boruta <- train[, all.vars(selected_formula)]

# tune the model and predict on the testing data
control <- trainControl(method="LOOCV", classProbs = TRUE, summaryFunction = defaultSummary, savePredictions = "all")
set.seed(43)
tunegrid <- expand.grid(.mtry=c(3:15))
rf_gridsearch <- train(Letter~., 
                       data=train_boruta, 
                       preProcess = "scale", 
                       method="rf", 
                       metric="Accuracy", 
                       tuneGrid=tunegrid, 
                       trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)
ggplot(rf_gridsearch) + 
  theme_bw() +
  labs(title = "Random Forest ROC vs. mtry",
       x = "mtry", 
       y = "ROC")

set.seed(333) 

mtry = rf_gridsearch$bestTune$mtry
rf_model <- randomForest(Letter~., data = train_boruta, ntree=500, mtry = mtry)

colnames(test) <- make.names(colnames(test))
test_boruta <- test[, all.vars(selected_formula)]
test_boruta <- test_boruta %>% 
  dplyr::select(-c(Letter))

p <- as.numeric(predict(rf_model, test_boruta)) -1
p <- as.factor(p)
y_test <- as.factor(y_test)
confusionMatrix(p, y_test)

roc_obj <- roc(as.numeric(y_test), as.numeric(p))
auc(roc_obj)


# Find the overlap between each dataset
print("Intersect: boruta & Lasso")
intersect(confirmed_features, make.names(lasso_features))

print("Intersect: Lasso & ttest")
intersect(lasso_features, ba_ttest)

print("Intersect: boruta & ttest")
intersect(confirmed_features, make.names(ba_ttest))

as.vector(ba_ttest)



# Overlapping features:

# Manually make a list of all the overlapping features to extract
overlap_features <- c("Glycodeoxycholic Acid",
                      "Glycochenodeoxycholic Acid-3-Sulfate",
                      "6-Oxolithocholic Acid",
                      "5-beta-Cholanic Acid 12-alpha-ol-3-one",
                      "12-Ketochenodeoxycholic acid|7-Ketodeoxycholic acid",
                      "Chenodeoxycholic Acid",
                      "Deoxycholic Acid",
                      "Lithocholic Acid",
                      "Cholic acid|Ursocholic acid",
                      "Isolithocholic Acid",
                      "5-beta-Cholanic Acid-3-beta, 12-alpha-diol",
                      "Glycocholic Acid-3-Sulfate",
                      "3-alpha-Hydroxy-12 Ketolithocholic Acid",
                      "Hyodeoxycholic Acid") 


bile_acid_features <- bile_acids_long  %>%
  filter(bile_acid %in% overlap_features)


label_positions <- bile_acid_features %>%
  group_by(bile_acid) %>%
  summarise(y_pos = max(Value, na.rm = TRUE) * 1.05)

annotation_plot_df <- bile_acid_table %>%
  merge(label_positions, by = "bile_acid")

bxp <- ggboxplot(
  bile_acid_features, x = "Time", y = "Value",
  color = "Time", palette = "jco",
  facet.by = "bile_acid", short.panel.labs = FALSE
)
bxp + 
  facet_wrap(~bile_acid) +
  labs(
    x = "Baseline vs Week 4",
    y = "Bile Acid distribution"
  ) +
  geom_text(
    data = annotation_plot_df,
    aes(x = 2, y = 4, label = paste("p-value:", round(P_Adj, 3))),  
    inherit.aes = FALSE,
    size = 4,
    hjust = 0.5
  )

annotation_plot_df

# Extract the features into another df and save it
bile_acid_reduced <- bile_acids %>%
  dplyr::select(c("Number", "Letter", all_of(overlap_features)))


write.csv(bile_acid_reduced, "Bile_acid_reduced.csv")
