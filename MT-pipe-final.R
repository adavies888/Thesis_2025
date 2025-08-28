# Meta-transcriptomics analysis

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


### Load data
mt_data <- read_excel("M7 alternate/Metatranscriptomics (KEGG).xlsx")

mt_data <- mt_data %>%
  dplyr::select(-c("Patient", "Timepoint")) %>%
  rename(SampleID = Sample_ID) %>%
  mutate(across(-SampleID, as.numeric))

mt_data <- mt_data %>% 
  extract(SampleID, into = c("Number", "Letter"), regex = "([0-9]+)([A-Z])") %>%
  mutate(Number = factor(Number, levels = sort(as.numeric(unique(Number)))))


# Reduce dataset to just weeks 0 and 4
mt_data <- mt_data %>%
  filter(Letter %in% c("A", "C")) %>%
  mutate(Letter = factor(Letter, levels = sort(unique(Letter))))

letter_to_time <- c(A = "Baseline", C = "Week 4") 
mt_data$Time <- letter_to_time[mt_data$Letter]

# Investigate missingness
f <- mt_data %>%
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
f_2 <- mt_data
f_2[] <- lapply(f_2, function(x) { x[x == 0] <- NA; return(x) })

missing_values_per_col <- colSums(is.na(f_2)) < 17 # 40%

mt_data_2 <- mt_data %>% 
  dplyr::select(which(missing_values_per_col)) 

f <- mt_data_2 %>%
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

f <- mt_data_2 %>%
  dplyr::select( -c("Number", "Letter", "Time"))

f[] <- lapply(f, function(x) { x[x == 0] <- NA; return(x) })

missing_values_per_row <- as.numeric(round((rowSums(is.na(f))/dim(f)[2])*100,2))

missing_values_per_row

MissingIDs <- mt_data_2 %>% 
  add_column(Missing = missing_values_per_row)# %>%

ggplot(MissingIDs, aes(Missing)) +
  geom_histogram(stat = "bin", binwidth = 1) + theme_bw() + 
  labs(x = "Missing (%)")

# Create the final dataset
FinalDatasetMissing <- mt_data_2[-c(1, 6), ] # 1A and 6A are the missing samples

mt_data <- FinalDatasetMissing

mt_list <- colnames(mt_data)[-c(1, 2, 1174)]

mt_data_long <- mt_data %>%
  pivot_longer(
    cols = -c(Number, Letter, Time),
    values_to = "Value",
    names_to = "mt"
  )


# PCA
pca.res <- prcomp(mt_data[,-c(1,2,1174)])

summary(pca.res)

loadings <- pca.res$rotation

pca_var <- pca.res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)

x_label <- paste0("PC1 (", pca_var_perc[1], "%)")
y_label <- paste0("PC2 (", pca_var_perc[2], "%)")

# Plot
fviz_pca_biplot(pca.res,
                repel = TRUE,
                label = "none",
                invisible = "var",
                habillage = mt_data$Time,
                palette = c("#00AFBB", "#B6ACE7"), 
                addEllipses = TRUE, ellipse.level = 0.9) +
  labs(x = x_label, y = y_label,  title = "(b) Meta-transcriptomics: PC 1 against 2") + 
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18))


# T-test to see change between week 0 and 4
mt_data_long %>%
  group_by(Letter, mt) %>%
  get_summary_stats(Value, type = "mean_sd")

pvals <- c()

for (mt_item in mt_list) {
  test_df <- mt_data_long %>%
    filter(mt %in% mt_item)
  
  res <- t.test(Value ~ Letter, data = test_df)
  
  pvals[mt_item] <- res$p.value
}

pvals_adj <- p.adjust(pvals, method = "fdr", n = length(pvals))

mt_table <- data.frame(
  mt = names(pvals),
  P_Value = as.numeric(pvals),
  P_Adj = as.numeric(pvals_adj)
)

# sort by p-value
mt_table <- mt_table[order(mt_table$P_Adj), ]

# View the table
print(mt_table)

# Plot the t-test ones < 0.05
mt_ttest <- mt_table %>%
  filter(P_Adj < 0.05) %>%
  pull(mt) %>%
  as.list()

mt_data_scaled <- mt_data %>% 
  mutate(across(-c(1, 2, 1174), scale))

mt_data_long <- mt_data_scaled %>%
  pivot_longer(
    cols = -c(Number, Letter, Time),
    values_to = "Value",
    names_to = "mt"
  )

mt_features <- mt_data_long  %>%
  filter(mt %in% mt_ttest)

label_positions <- mt_features %>%
  group_by(mt) %>%
  summarise(y_pos = max(Value) * 1.05)

annotation_plot_df <- mt_table %>%
  merge(label_positions, by = "mt")

bxp <- ggboxplot(
  mt_features, x = "Time", y = "Value",
  color = "Time", palette = "jco",
  facet.by = "mt", short.panel.labs = FALSE
)
bxp + 
  facet_wrap(~mt) + 
  labs(
    x = "Baseline vs Week 4",
    y = "MT distribution") +
  geom_text(
    data = annotation_plot_df,
    aes(x = 2.2, y = 4, label = paste("p-value:", round(P_Adj, 3))),  
    inherit.aes = FALSE,
    size = 4,
    hjust = 0.5
  )

# LASSO
test <- mt_data %>%
  dplyr::select(-Time) %>%       
  filter(Number %in% c(2,3)) # test/train split

train <- mt_data %>%
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

plot(lasso_mod, xvar = "lambda")# label = TRUE)    # Draw plot of coefficients

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
                     s = bestlam) 

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
myResults$features <- factor(myResults$features, 
                             levels = myResults$features[order(abs(myResults$coefs), 
                                                               decreasing = FALSE)])


lasso_features <- as.vector(myResults$features) # extract the lasso features
lasso_features


# Plot the lasso features
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

# Extract important features
confirmed_features <- rownames(boruta_stats[boruta_stats$decision == 'Confirmed', ])
print(confirmed_features)

selected_formula <- getNonRejectedFormula(boruta)

train_boruta <- train[, all.vars(selected_formula)]

# train the model and predict
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

# Compare for overlap
print("Intersect: boruta & Lasso")
intersect(confirmed_features, make.names(lasso_features))

print("Intersect: Lasso & ttest")
intersect(lasso_features, mt_ttest)

print("Intersect: boruta & ttest")
intersect(confirmed_features, make.names(mt_ttest))


# Overlapping features

# Manually make a list of all the overlapping features to extract
overlap_features <- c(
  "K01588: 5-(carboxyamino)imidazole ribonucleotide mutase [EC:5.4.99.18]",
  "K03088: RNA polymerase sigma-70 factor, ECF subfamily",
  "K03152: 4-methyl-5(b-hydroxyethyl)-thiazole monophosphate biosynthesis",
  "K09014: Fe-S cluster assembly protein SufB",
  "K07473: DNA-damage-inducible protein J",
  "K10914: CRP/FNR family transcriptional regulator, cyclic AMP receptor protein",
  "K00041: tagaturonate reductase [EC:1.1.1.58]",
  "K00648: 3-oxoacyl-[acyl-carrier-protein] synthase III [EC:2.3.1.180]",
  "K01262: X-Pro aminopeptidase",
  "K01681: aconitate hydratase [EC:4.2.1.3]",
  "K01736: chorismate synthase [EC:4.2.3.5]",
  "K01808: ribose 5-phosphate isomerase B [EC:5.3.1.6]",
  "K01835: phosphoglucomutase [EC:5.4.2.2]",
  "K02874: large subunit ribosomal protein L14",
  "K03154: sulfur carrier protein",
  "K03324: phosphate:Na+ symporter, PNaS family",
  "K09794: hypothetical protein",
  "K11717: cysteine desulfurase / selenocysteine lyase [EC:2.8.1.7 4.4.1.16]",
  "K14652: 3,4-dihydroxy 2-butanone 4-phosphate synthase / GTP cyclohydrolase II [EC:4.1.99.12 3.5.4.25]",
  "K15633: 2,3-bisphosphoglycerate-independent phosphoglycerate mutase [EC:5.4.2.12]",
  "K01284: peptidyl-dipeptidase Dcp",
  "K01809: mannose-6-phosphate isomerase [EC:5.3.1.8]",
  "K00286: pyrroline-5-carboxylate reductase [EC:1.5.1.2]",
  "K01719: uroporphyrinogen-III synthase [EC:4.2.1.75]",
  "K03701: excinuclease ABC subunit A",
  "K18929: NO_NAME",
  "K01662: 1-deoxy-D-xylulose-5-phosphate synthase [EC:2.2.1.7]",
  "K01265: methionyl aminopeptidase",
  "K03303: lactate transporter, LctP family"
)

mt_features <- mt_data_long  %>%
  filter(mt %in% overlap_features)


label_positions <- mt_features %>%
  group_by(mt) %>%
  summarise(y_pos = max(Value, na.rm = TRUE) * 1.05)

annotation_plot_df <- mt_table %>%
  merge(label_positions, by = "mt")

bxp <- ggboxplot(
  mt_features, x = "Time", y = "Value",
  color = "Time", palette = "jco",
  facet.by = "mt", short.panel.labs = FALSE
)
bxp + 
  facet_wrap(~mt) + 
  labs(
    x = "Baseline vs Week 4",
    y = "Feature distribution"
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
mt_reduced <- mt_data %>%
  dplyr::select(c("Number", "Letter", all_of(overlap_features)))

write.csv(mt_reduced, "Metatranscriptomics_reduced.csv")