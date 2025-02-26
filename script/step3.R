
rm(list = ls())

set.seed(52)

library(glmnet)
load("Rdata/mydata_DMP_DMR.Rdata")
myDMR = myDMR$BumphunterDMR
rm(myDMP)

myDMR$group = ifelse(myDMR$value > 0, "Hypermethylated", "Hypomethylated")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

write.table(myDMR,file = "myDMR.txt",row.names = F,sep = "\t",quote = F)

myDMR <- annotatePeak("./myDMR.txt", tssRegion = c(-3000, 3000), TxDb = txdb,annoDb = 'org.Hs.eg.db')

myDMR = as.data.frame(myDMR@anno@elementMetadata@listData)

myDMR %>% 
  select(SYMBOL) %>% 
  distinct()

myDMR %>%
  count(SYMBOL) %>%
  filter(n > 1) 

myDMR_anno = myDMR

load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM
rm(myDMR)

dmp = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% myDMR_anno$SYMBOL)

rm(myDMP)

load("Rdata/mydata_beta_pd.Rdata")

mydata_beta = as.data.frame(mydata_beta)
mydata_beta = mydata_beta[dmp$cgname,]
mydata_beta = as.data.frame(t(mydata_beta))
mydata_beta = scale(mydata_beta)

mydata_pd$Sample_Group = factor(mydata_pd$Sample_Group, levels = c("Control", "GDM"))

# Model
la.eq <- glmnet(mydata_beta, mydata_pd$Sample_Group, 
                family='binomial', 
                intercept = F, alpha=1) 

la.eq$beta[,1]
tiff(file = "./Figure/step3_1.tiff", 
     width = 8 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
plot(la.eq,xvar = "lambda", label = F)
dev.off()

mod_cv <- cv.glmnet(mydata_beta, mydata_pd$Sample_Group, family="binomial", # 默认nfolds = 10
                    intercept = F, alpha=1)

tiff(file = "./Figure/step3_2.tiff", 
     width = 8 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
plot(mod_cv)
dev.off()

print(paste(mod_cv$lambda.min,
            log(mod_cv$lambda.min)))
print(paste(mod_cv$lambda.1se,
            log(mod_cv$lambda.1se)))

best_lambda <- mod_cv$lambda.min
best_lambda

best_model <- glmnet(mydata_beta, mydata_pd$Sample_Group, alpha = 1, lambda = best_lambda, family='binomial')
coef(best_model)
coef_best_model = coef(best_model)
coef_best_model_df <- as.data.frame(as.matrix(coef_best_model))
hubcg = coef_best_model_df %>% 
  rownames_to_column("cgnames") %>% 
  filter(s0 != 0)
hubcg = hubcg[-1,]

dmp %>% 
  filter(cgname %in% hubcg$cgnames) %>% 
  select(gene)

###
require(readr)
require(ggplot2)
require(dplyr)
require(tidyr)
require(caret)
require(corrplot)
require(Hmisc)
require(parallel)
require(doParallel)
require(ggthemes)
library(ggsci)

set.seed(52)

hub_cg = c("cg19037167", "cg22998811", "cg02249039", "cg15150348", "cg17537719", "cg26538349", "cg26792694")

load("Rdata/mydata_beta_pd.Rdata")

mydata_beta = as.data.frame(mydata_beta)
mydata_beta = mydata_beta[hub_cg,]
mydata_beta = as.data.frame(t(mydata_beta))

data = mydata_beta
data$group = mydata_pd$Sample_Group
data$group = factor(data$group, levels = c("Control", "GDM"))

### separate dataset into training and testing sets
sample_Index <- createDataPartition(data$group, p=0.7,list=FALSE)
voice_Train <- data[sample_Index,]
voice_Test <- data[-sample_Index,]

### preprocess factors for further modeling
#pp <- preProcess(voice_Train,method=c("scale","center","pca"))
pp <- preProcess(voice_Train,method=c("scale","center"))
voice_Train <- predict(pp,voice_Train)
voice_Test <- predict(pp,voice_Test)

### define formula
model_Formula <- group~cg19037167+cg22998811+cg02249039+cg15150348+cg17537719+cg26538349+cg26792694
#model_Formula <- group~PC1

###set cross-validation parameters
modelControl <- trainControl(method="repeatedcv",number=5,
                             repeats=5,allowParallel=TRUE, 
                             classProbs = TRUE)


### logistic regression
glm_Model <- train(model_Formula,
                   data=voice_Train,
                   method="glm",
                   trControl=modelControl)

#glm_Coefficients <- coef(glm_Model$finalModel)

### linear discrimant analysis
lda_Model <- train(model_Formula,
                   data=voice_Train,
                   method="lda",
                   trControl=modelControl)

### random forrest
rf_Model <- train(model_Formula,
                  data=voice_Train,
                  method="rf",
                  trControl=modelControl,
                  ntrees=500)

### Naive Bayes
nb_Model <- train(model_Formula,
                  data=voice_Train,
                  method="nb",
                  trControl=modelControl)

### SVM
svm_Model <- train(model_Formula,
                   data=voice_Train,
                   method="svmLinear",
                   trControl=modelControl,
)

### xgboost
xgboost_Model <- train(model_Formula,
                       data=voice_Train,
                       method="xgbLinear",
                       trControl=modelControl,
)

### nnet
nnet_Model <- train(model_Formula,
                    data=voice_Train,
                    method="nnet",
                    trControl=modelControl,
)

### knn
knn_Model <- train(model_Formula,
                   data=voice_Train,
                   method="knn",
                   trControl=modelControl,
)

### rpart
rpart_Model <- train(model_Formula,
                     data=voice_Train,
                     method="C5.0",
                     trControl=modelControl,
)

### AdaBoost
AdaBoost_Model <- train(model_Formula,
                        data=voice_Train,
                        method="AdaBoost.M1",
                        trControl=modelControl,
)

library(pROC)

a = voice_Test
a$glm_Prob <- predict(glm_Model, voice_Test, type = "prob")[, 2]
a$lda_Prob <- predict(lda_Model, voice_Test, type = "prob")[, 2]
a$rf_Prob <- predict(rf_Model, voice_Test, type = "prob")[, 2]
a$nb_Prob <- predict(nb_Model, voice_Test, type = "prob")[, 2]
a$svm_Prob <- predict(svm_Model, voice_Test, type = "prob")[, 2]
a$xgboost_Prob <- predict(xgboost_Model, voice_Test, type = "prob")[, 2]
a$nnet_Prob <- predict(nnet_Model, voice_Test, type = "prob")[, 2]
a$knn_Prob <- predict(knn_Model, voice_Test, type = "prob")[, 2]
a$rpart_Prob <- predict(rpart_Model, voice_Test, type = "prob")[, 2]
a$AdaBoost_Prob <- predict(AdaBoost_Model, voice_Test, type = "prob")[, 2]

colors=pal_futurama()(11)

tiff(file = "./Figure/step3_3.tiff", 
     width = 8 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
plot.roc(
  a$group,              
  a$glm_Prob,       
  col = colors[1],            
  percent = TRUE,             
  lwd = 2,                   
  print.auc = TRUE,           
  print.auc.cex = 1,          
  print.auc.pattern = "Logistic Regression: AUC=%.1f%%",      
  print.auc.y = 45
)

# 继续添加其他 ROC 曲线，参数与上述相似
plot.roc(a$group,    
         a$lda_Prob,
         col = colors[2],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Linear Regression: AUC=%.1f%%",
         print.auc.y = 40,
         add = T)

plot.roc(a$group,    
         a$rf_Prob,
         col = colors[3],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Random Forest: AUC=%.1f%%",
         print.auc.y = 35,
         add = T)

plot.roc(a$group,    
         a$nb_Prob,
         col = colors[4],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Naive Bayes: AUC=%.1f%%",
         print.auc.y = 30,
         add = T)

plot.roc(a$group,    
         a$svm_Prob,
         col = colors[5],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Support Vector Machine: AUC=%.1f%%",
         print.auc.y = 25,
         add = T)

plot.roc(a$group,    
         a$xgboost_Prob,
         col = colors[7],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Extreme Gradient Boosting: AUC=%.1f%%",
         print.auc.y = 20,
         add = T)

plot.roc(a$group,    
         a$nnet_Prob,
         col = colors[8],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Artificial Neural Network: AUC=%.1f%%",
         print.auc.y = 15,
         add = T)

plot.roc(a$group,    
         a$knn_Prob,
         col = colors[9],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "k-Nearest Neighbors: AUC=%.1f%%",
         print.auc.y = 10,
         add = T)

plot.roc(a$group,    
         a$rpart_Prob,
         col = colors[10],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Decision Tree: AUC=%.1f%%",
         print.auc.y = 5,
         add = T)

plot.roc(a$group,    
         a$AdaBoost_Prob,
         col = colors[11],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Adaptive Boosting: AUC=%.1f%%",
         print.auc.y = 0,
         add = T)
dev.off()

model_coef <- coef(glm_Model$finalModel)
print(model_coef)

logit(P)= 12.09 - 7.26 * cg19037167 + 0.16 * cg22998811 - 6.89 * cg02249039 - 6.22 * cg15150348 +
  2.20 * cg17537719 - 11.46 * cg26538349 - 0.96 * cg26792694

###
rm(list = ls())

library(xlsx)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggExtra)

rdata <- read.xlsx(file = "./shifuyou_sample_data.xlsx",
                   2,header = TRUE)
rdata <- as.data.frame(t(rdata))
colnames(rdata) <- rdata[1,]
rdata <- rdata[-1,]

rdata[c(1:13,16)] <- as.data.frame(map(rdata[c(1:13,16)],as.numeric))

hub_cg = c("cg19037167", "cg22998811", "cg02249039", "cg15150348", "cg17537719", "cg26538349", "cg26792694")

load("Rdata/mydata_beta_pd.Rdata")

mydata_beta = mydata_beta[hub_cg,]
mydata_beta = as.data.frame(t(mydata_beta))

data = mydata_beta
data$sampleID = mydata_pd$Sample_Name
data$sampleID = substr(data$sampleID,16,18)

data = data %>% 
  as_tibble() %>% 
  inner_join(rdata, by = "sampleID")

data_long <- data %>%
  select(c(hub_cg, 'Gestational_weeks')) %>%
  pivot_longer(cols = all_of(hub_cg), names_to = 'hub_cg', values_to = 'value')

# Create the plot
tiff(file = "./Figure/step3_4_Gestational weeks.tiff", 
     width = 27 * 300,
     height = 4 * 300, 
     res = 300, 
     compression = "lzw")
ggplot(data_long, aes(x = value, y = Gestational_weeks)) +
  geom_point(size = 2, color = '#EC0101', alpha = 0.5) +
  theme_bw() +  
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()) +
  geom_smooth(method = 'lm', se = T, color = '#F9B208', size = 1.5, fill = '#FEA82F') +
  stat_cor(method = "pearson", digits = 3, size = 6) +
  facet_wrap(~hub_cg, scales = 'free', ncol = 7)+
  ylab("Gestational weeks") +
  xlab("Methylation value")
dev.off()