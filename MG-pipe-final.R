# Meta-genomics pipe

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
library(Boruta)
library(mlbench)
library(caret)
library(randomForest)


# Load data
mg_data <- read_excel("M7 alternate/Metagenomics dataset (species).xls")

mg_data <- mg_data %>%
  dplyr::select(-c("Patient", "Timepoint")) %>%
  rename(SampleID = Sample_ID) %>%
  mutate(across(-SampleID, as.numeric))

mg_data <- mg_data %>% 
  extract(SampleID, into = c("Number", "Letter"), regex = "([0-9]+)([A-Z])") %>%
  mutate(Number = factor(Number, levels = sort(as.numeric(unique(Number)))))

# Reduce dataset to just weeks 0 and 4
mg_data <- mg_data %>%
  filter(Letter %in% c("A", "C")) %>%
  mutate(Letter = factor(Letter, levels = sort(unique(Letter))))

letter_to_time <- c(A = "Baseline", C = "Week 4") 
mg_data$Time <- letter_to_time[mg_data$Letter]

# Investigate missingness
f <- mg_data %>%
  dplyr::select( -c("Number", "Letter", "Time"))

f[] <- lapply(f, function(x) { x[x == 0] <- NA; return(x) })

ll <- data.frame(is.na(f))

cols <- sapply(ll, is.logical)

ll[, cols] <- lapply(ll[, cols], as.numeric)

Miss1 <- UpSetR::upset(ll,
                       
                       nsets = 20, number.angles = 10, point.size = 3.5, line.size = 2,
                       
                       mainbar.y.label = "Missing Values", sets.x.label = "Total Number Missing Values",
                       
                       text.scale = c(2.3, 2.3, 2, 2, 2, 1.75), order.by = "freq", sets.bar.color = "red3"
                       
)

Miss1

print(vis_miss(f, warn_large_data = FALSE) + 
        theme(axis.text.x =  element_blank()))

# Remove missing values
f_2 <- mg_data
f_2[] <- lapply(f_2, function(x) { x[x == 0] <- NA; return(x) })

missing_values_per_col <- colSums(is.na(f_2)) < 21 # 30%

mg_data_2 <- mg_data %>% 
  dplyr::select(which(missing_values_per_col)) 

f <- mg_data_2 %>%
  dplyr::select( -c("Number", "Letter", "Time"))

f[] <- lapply(f, function(x) { x[x == 0] <- NA; return(x) })

ll <- data.frame(is.na(f))

cols <- sapply(ll, is.logical)

ll[, cols] <- lapply(ll[, cols], as.numeric)

Miss1 <- UpSetR::upset(ll,
                       
                       nsets = 20, number.angles = 10, point.size = 3.5, line.size = 2,
                       
                       mainbar.y.label = "Missing Values", sets.x.label = "Total Number Missing Values",
                       
                       text.scale = c(2.3, 2.3, 2, 2, 2, 1.75), order.by = "freq", sets.bar.color = "red3"
                       
)

Miss1

print(vis_miss(f, warn_large_data = FALSE) + 
        theme(axis.text.x =  element_blank()))

f <- mg_data_2 %>%
  dplyr::select( -c("Number", "Letter", "Time"))

f[] <- lapply(f, function(x) { x[x == 0] <- NA; return(x) })

missing_values_per_row <- as.numeric(round((rowSums(is.na(f))/dim(f)[2])*100,2))

missing_values_per_row

MissingIDs <- mg_data_2 %>% 
  add_column(Missing = missing_values_per_row)# %>%

ggplot(MissingIDs, aes(Missing)) + 
  geom_histogram(stat = "bin", binwidth = 1) + theme_bw() + 
  labs(x = "Missing (%)")


# Create final dataset
FinalDatasetMissing <- mg_data_2[-c(16, 26, 29), ] # 1C, 11C and 14C have the most amount missing  

mg_data <- FinalDatasetMissing

mg_data_list <- colnames(mg_data)[-c(1, 2, 50)]

mg_data_long <- mg_data %>%
  pivot_longer(
    cols = -c(Number, Letter, Time),
    values_to = "Value",
    names_to = "mg"
  )


# PCA
pca.res <- prcomp(mg_data[,-c(1,2,50)])

summary(pca.res)

loadings <- pca.res$rotation

pca_var <- pca.res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)

x_label <- paste0("PC1 (", pca_var_perc[1], "%)")
y_label <- paste0("PC2 (", pca_var_perc[2], "%)")

# plot
fviz_pca_biplot(pca.res,
                repel = TRUE,
                label = "none",
                invisible = "var",
                habillage = mg_data$Time,
                palette = c("#00AFBB", "#B6ACE7"), 
                addEllipses = TRUE, ellipse.level = 0.9) +
  labs(x = x_label, y = y_label,  title = "(c) Meta-genomics: PC 1 against 2") + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18)
  )

# T-test to see change between week 0 and 4
mg_data_long %>%
  group_by(Letter, mg) %>%
  get_summary_stats(Value, type = "mean_sd")

pvals <- c()

for (mg_item in mg_data_list) {
  test_df <- mg_data_long %>%
    filter(mg %in% mg_item)
  
  res <- t.test(Value ~ Letter, data = test_df)
  
  pvals[mg_item] <- res$p.value
}

pvals_adj <- p.adjust(pvals, method = "fdr", n = length(pvals))

mg_data_table <- data.frame(
  mg = names(pvals),
  P_Value = as.numeric(pvals),
  P_Adj = as.numeric(pvals_adj)
)

# sort by p-value
mg_data_table <- mg_data_table[order(mg_data_table$P_Adj), ]

# View the table
print(mg_data_table)

# Plot the t-test ones < 0.05
mg_ttest <- mg_data_table %>%
  filter(P_Adj < 0.05) %>%
  pull(mg) %>%
  as.list()

mg_data_scaled <- mg_data %>% 
  mutate(across(-c(1, 2, 50), scale))

mg_data_long <- mg_data_scaled %>%
  pivot_longer(
    cols = -c(Number, Letter, Time),
    values_to = "Value",
    names_to = "mg"
  )

mg_data_features <- mg_data_long  %>%
  filter(mg %in% mg_ttest)

label_positions <- mg_data_features %>%
  group_by(mg) %>%
  summarise(y_pos = max(Value, na.rm = TRUE) * 1.05)

annotation_plot_df <- mg_data_table %>%
  merge(label_positions, by = "mg")

bxp <- ggboxplot(
  mg_data_features, x = "Time", y = "Value",
  color = "Time", palette = "jco",
  facet.by = "mg", short.panel.labs = FALSE
)
bxp + 
  facet_wrap(~mg) + 
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


# Logistic regression: LASSO
test <- mg_data %>%
  dplyr::select(-Time) %>%       
  filter(Number %in% c(1,3)) # test/train split

train <- mg_data %>%
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

# Create the model
lasso_mod <- glmnet(x_train, 
                    y_train,
                    alpha = 1, 
                    family = "binomial", standardize = TRUE) # Fit lasso model on training data

plot(lasso_mod, xvar = "lambda")    # Draw plot of coefficients

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

# List the lasso features
lasso_features <- as.vector(myResults$features)
lasso_features


# Plot the lasso featues
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

# Predict
lasso_pred = predict(lasso_mod, s = bestlam, newx = x_test, type = "response")

lasso_pred

pred <- as.factor(ifelse(lasso_pred >= 0.5, 1, 0))

conf_matrix <- confusionMatrix(as.factor(y_test), pred)

# Print confusion matrix
print(conf_matrix)

# AUC
roc_obj <- roc(as.numeric(y_test), as.numeric(pred))
auc(roc_obj)


# Random forest 
# Boruta
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

confirmed_features <- rownames(boruta_stats[boruta_stats$decision == 'Confirmed', ])
print(confirmed_features)
attStats(boruta)

selected_formula <- getNonRejectedFormula(boruta)

train_boruta <- train[, all.vars(selected_formula)]

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

roc_obj <- roc(as.numeric(y_test), as.numeric(pred))
auc(roc_obj)

# Identify overlapping features
print("Intersect: boruta & Lasso")
intersect(confirmed_features, lasso_features)

print("Intersect: Lasso & ttest")
intersect(lasso_features, mg_ttest)

print("Intersect: boruta & ttest")
intersect(confirmed_features, mg_ttest)


# Overlapping features:
# Manually make a list of all the overlapping features to extract
overlap_features <- c(
  "Megasphaera",
  "Bacteroides",
  "Roseburia",
  "Blautia",
  "Faecalibacterium",
  "Clostridium",
  "Pediococcus",
  "Lachnospiraceae_unclassified",
  "Phocaeicola",
  "Ruminococcus",
  "Bifidobacterium",
  "Parabacteroides",
  "Dorea",
  "Lacrimispora",
  "Collinsella",
  "Anaerostipes",
  "Lachnospira",
  "Anaerobutyricum",
  "Mediterraneibacter",
  "Dysosmobacter",
  "Flavonifractor",
  "Eubacteriales_unclassified",
  "Agathobaculum",
  "Fusobacterium",
  "Lacticaseibacillus",
  "Lactiplantibacillus"
)

mg_data_features <- mg_data_long  %>%
  filter(mg %in% overlap_features)


label_positions <- mg_data_features %>%
  group_by(mg) %>%
  summarise(y_pos = max(Value, na.rm = TRUE) * 1.05)

annotation_plot_df <- mg_data_table %>%
  merge(label_positions, by = "mg")

bxp <- ggboxplot(
  mg_data_features, x = "Time", y = "Value",
  color = "Time", palette = "jco",
  facet.by = "mg", short.panel.labs = FALSE
)
bxp + 
  facet_wrap(~mg) + #, scales = "free_y")
  labs(
    x = "Baseline vs Week 4",
    y = "Feature distribution"
  ) +
  geom_text(
    data = annotation_plot_df,
    aes(x = 2.2, y = 4, label = paste("p-value:", round(P_Adj, 3))),  
    inherit.aes = FALSE,
    size = 4,
    hjust = 0.5
  )

annotation_plot_df

# Extract the features into another df and save it
mg_data_reduced <- mg_data %>%
  dplyr::select(c("Number", "Letter", all_of(overlap_features)))

write.csv(mg_data_reduced, "Metagenomics_reduced.csv")
