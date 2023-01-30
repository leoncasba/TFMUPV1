####################################
# Cargar librerias
####################################

library(readxl)
library(writexl)
library(FactoMineR)
library(factoextra)
library(caret)
library(mice)
library(knitr)
library(dplyr)
library(clusterSim)
library(ggplot2)
library(e1071)
library(tidyr)
library(randomForest)
library(smotefamily)
library(kernlab)
library(DMwR)
library(class)
library(FNN)
library(nnet)
library(ROCit)
library(rpart.plot)
library(lime)

####################################
# Preprocesamiento
####################################

data <- read_xlsx("sinanomalos.xlsx")
dim(data)


cuentaNAporcent <- function(x, ndec=2){
  porcent=(sum(is.na(x))/length(x))*100
  p2 = round(porcent, digits=ndec)
}


apply(data, 2, cuentaNAporcent)

summary(data)

#Tabla de rangos para verificar las variables esten dentro de rango consistente
rangos <- t(apply(data, 2, range, na.rm = TRUE))
colnames(rangos) <- c("Min","Max")
kable(rangos, caption="Tabla de rangos de las variables num?ricas", digits=0)


###Missing data
imputacion = mice(data[,-c(1:3)], seed = 123, m = 5, print = FALSE, method = NULL) #numero de imputaciones multiples=5
mice::stripplot(imputacion) #Grafico de imputacion
imputaciones <- mice::complete(imputacion)
data2 <- imputaciones

apply(data2, 2, cuentaNAporcent) #validacion imputacion

write_xlsx(data2, "tfm_imputada.xlsx") #exportar datos imputados Excel

####################################
# PCA - Validacion informacion similar/igual ASPEN
####################################



data2 <- read_xlsx("tfm_imputada_ASPEN.xlsx")
dim(data)

#PCA1 - Funcion PCA
res.pca = PCA(data2[,c(4:14)], scale.unit=TRUE, ncp=5, graph=T)

barplot(res.pca$eig[,1],main="Autovalores", names.arg=paste("dim",1:nrow(res.pca$eig)), ylim=c(0,20))
abline(h = 1, col="red")

#PCA2 - Funcion ACP
acp<-princomp(data2[,c(4:14)], scale=TRUE, cor = T)
acp$center
summary(acp)
loadings <- acp$loadings
scores <- acp$scores

#screen plot
fviz_eig(acp, choice="eigenvalue")

#grafico contribuciones
fviz_contrib(acp, axes=1, choice="var")
fviz_contrib(acp, axes=2, choice="var")
fviz_contrib(acp, axes=3, choice="var")
fviz_contrib(acp, axes=4, choice="var")

#biplot
fviz_pca_biplot(acp, label="var", labelsize=5, addEllipses=TRUE, ellipse.level=0.95)+theme_minimal()
fviz_pca_biplot(acp, label="var", labelsize=5, addEllipses=TRUE, ellipse.level=0.95, axes=c(1,2))

####################################
# Analisis Clustering
####################################

#escalar varianza unitaria
scaled_values <- apply(data2[,c(4:11, 13)], 2, scale, center = TRUE, scale = TRUE)

#Verificar si estan correctamente escalados
paste("Scaled values of dataset have mean equal", mean(round(apply(scaled_values, 2, mean,na.rm=TRUE))), "and standard deviation equal", mean(apply(scaled_values, 2, sd, na.rm=TRUE)))

#Heatmap
mat_dist <- get_dist(x = scaled_values, method = "euclidean")
fviz_dist(dist.obj = mat_dist, lab_size = 5) +
  theme(legend.position = "none")

#Elbow method

# Initializing total within sum of squares error: wss
wss <- 0

# For 1 to 10 cluster centers, save total within sum of squares to wss variable
for (i in 1:10) {
  test <- kmeans(scaled_values2, centers = i, nstart = 20, met)
  # 
  wss[i] <- test$tot.withinss
}

# Plot total within sum of squares vs. number of clusters
plot(1:10, wss, type = "b", 
     xlab = "N?mero de clusters", 
     ylab = "Suma de cuadrados dentro de clusters (WCSS)")


#Indice Davies Bouldin
set.seed(1995)
lst1 <- lapply(2:6, function(i) kmeans(scaled_values, centers=i, nstart = 20))
names(lst1) <- paste0("cluster_", 2:6)
k_ <- as_tibble(scaled_values) %>% mutate(k2=unlist(lst1$cluster_2[1]), k3=unlist(lst1$cluster_3[1]), 
                                         k4=unlist(lst1$cluster_4[1]), k5=unlist(lst1$cluster_5[1]), k6=unlist(lst1$cluster_6[1])) 
index_DB <- lapply(8:12, function(i) index.DB(k_, k_[,i], d=NULL, centrotypes="centroids", p=2, q=2))
names(index_DB) <- paste0("index_cluster_k", 2:6)
DB_index <- rbind(index_DB$index_cluster_k2[1], index_DB$index_cluster_k3[1],index_DB$index_cluster_k4[1], index_DB$index_cluster_k5[1], index_DB$index_cluster_k6[1])
DB_index_ <- as_tibble(unlist(DB_index))
DB_index_$k <- c(2, 3, 4, 5, 6)
kable(DB_index_, digits=2, col.names=c("DB index", "Number of clusters"), caption="**Table 4.** Davies-Bouldin Index Values.")
index_DB2 <- index.DB(k_[, -c(9:12)], k_$k2, d=NULL, centrotypes="centroids", p=2, q=2)
index_DB3 <- index.DB(k_[, -c(8, 10:12)], k_$k3, d=NULL, centrotypes="centroids", p=2, q=2)
index_DB4 <- index.DB(k_[, -c(8, 9, 11, 12)], k_$k4, d=NULL, centrotypes="centroids", p=2, q=2)
index_DB5 <- index.DB(k_[, -c(8:10, 12)], k_$k5, d=NULL, centrotypes="centroids", p=2, q=2)
index_DB6 <- index.DB(k_[, -c(8:11)], k_$k6, d=NULL, centrotypes="centroids", p=2, q=2)
index_DB2$DB
index_DB3$DB
index_DB4$DB
index_DB5$DB

#Analisis silueta
fviz_nbclust(scaled_values, pam, method = "silhouette")+ theme_classic()

#2 parece el optimo

#K-means algorythm
set.seed(1995)
kmeans2 <- kmeans(scaled_values2, centers=2, nstart = 5)
kmeans3 <- kmeans(scaled_values2, centers=3, nstart = 5)

fviz_cluster(k_definitive, data = scaled_values2,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FF0000"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_classic(base_size = 18),
             main = "",
             ylab= "P2",
             xlab = "P1")
             + theme (plot.title = element_text (hjust = 0.5))+ scale_colour_viridis_d() +
  scale_fill_viridis_d()

#Fuzzy clustering

fuzzy2 <- cmeans(x=scaled_values2, centers = 2, iter.max=100, dist='euclidean', method='cmeans')
fuzzy2$centers
fuzzy2$membership

#Grafico de kmeans+Proyeccion PCA

fviz_pca_biplot(PCAdata, label="n", habillage=kmeans2$cluster,
                addEllipses=TRUE, ellipse.level=0.95, title ="K-means + Proyecci?n del PCA", axis.title = element_text(size = 2.5), 
                fill.ind = "red")

fviz_pca_biplot(PCAdata, label="n", habillage=fuzzy2$cluster,
                addEllipses=TRUE, ellipse.level=0.95, title ="Fuzzy-Clustering + Proyecci?n del PCA", 
                palette = "Set2", ggtheme = theme_minimal())

#Describiendo los cluster

# Extraer el vector de cluster del modelo
clust_customers <- fuzzy2$cluster
# Construir el data frame con el cluster asignado en cada empresa
empresas_k_clusters <- mutate(data2[, c(4:11, 13)], Cluster = clust_customers)

# Subset the data to select the observations in each cluster
cluster_1 <- empresas_k_clusters[empresas_k_clusters$Cluster == 1, ]
cluster_2 <- empresas_k_clusters[empresas_k_clusters$Cluster == 2, ]

# Calcular el tamano de cada cluster
count <- count(empresas_k_clusters, Cluster)
# Calcular el promedio de variables en cada categoria de cluster + anadir columna conteo y porcentaje
Clusters_final <- empresas_k_clusters %>% 
  group_by(Cluster) %>% 
  summarise_all(list(mean)) %>% mutate(N=count$n, Percentage=((count$n)/sum(count$n)*100))

# Realizar el t-test para cada variable y encontrar diferencias entre cluster
for(variable in names(empresas_k_clusters)) {
  if(variable != "Cluster") {
    result <- t.test(cluster_1[, variable], cluster_2[, variable])
    print(paste(variable, ": p-value =", result$p.value))
  }
}

#Clusters
write_xlsx(empresas_k_clusters, "tfmobjetivo1.xlsx") #exportar datos  Excel


###########################################
#Modelos sin cribar
##########################################

huella <- read_excel("tfmobjetivo2 - Copy.xlsx")
huella$Cluster <- factor(huella$Cluster)
huella <- huella[, -1]

#Balanceo de clases 

#Particion de los datos

trainIndex <- createDataPartition(huella$Cluster, p = 0.7, 
                                  list = FALSE, 
                                  times = 1)

#Bases de datos de test y training son creadas:
huellatrain <- huella[trainIndex,]
huellatest  <- huella[-trainIndex,]

#Oversampling con SMOTE
balanced.data <- smotefamily::SMOTE(huellatrain[,-1],huellatrain[,1],K=3,dup_size=2)
balance <- balanced.data$data
colnames(balance)[44] <- "Cluster"

table(huellatrain$Cluster) #unbalanced
table(balance$Cluster) #balanced

write_xlsx(balance, "tfm_balanceadas_.xlsx") #exportar datos  Excel para Aspen posteriormente


#Hold out repetido modelos clasificacion

tasa.aciertos.rf <- NA
tasa.aciertos.regresion<-NA
tasa.aciertos.svp<-NA
tasa.aciertos.plsda1<-NA
tasa.aciertos.tree<-NA
tasa.aciertos.knn<-NA
sensitivity.plsda1 <- NA
specificity.plsda1 <- NA
sensitivity.tree <- NA
specificity.tree <- NA
sensitivity.rf <- NA
specificity.rf <- NA
sensitivity.regresion <- NA
specificity.regresion <- NA
sensitivity.svp <- NA
specificity.svp <- NA
sensitivity.knn <- NA
specificity.knn <- NA

for(i in 1:100){
  library(smotefamily)
  cjo<-sample(1:nrow(huella),nrow(huella)*0.7, FALSE) 
  train <-huella[cjo,]
  test<-huella[-cjo,]
  train_balanced <- smotefamily::SMOTE(train[,-1],train[,1],K=3)
  train_balanced_ <- (train_balanced$data)
  colnames(train_balanced_)[44] <- "Cluster"
  train_balanced_$Cluster <- ifelse(train_balanced_$Cluster == 1, 0, 1)
  test$Cluster<- ifelse(test$Cluster == 1, 0, 1)
  test$Cluster <- factor(test$Cluster)
  train_balanced_$Cluster <- factor(train_balanced_$Cluster)
  
  #PLSDA
  plsda1 <- plsda(train_balanced_[,-44], as.factor(train_balanced_$Cluster), ncomp = 3, cv = 50)
  pred0 <- predict(plsda1, test[,-1], type = "class",ncomp = 3)
  tab_plsda <- table(pred0, test$Cluster) 
  #Tasa aciertos plsda1
  tasa.aciertos.plsda1[i]<-sum(tab_plsda[row(tab_plsda)==col(tab_plsda)])/sum(tab_plsda)
  sensitivity.plsda1[i] <- sensitivity(tab_plsda)
  specificity.plsda1[i] <- specificity(tab_plsda)
 
  
  ##Random forest
  rf<-randomForest(factor(train_balanced_$Cluster)~ .,data=as.data.frame(train_balanced_[,-44]),mtry=4,method="class",importance=TRUE, ntree=5)
  pred2 <- predict(rf, newdata=test[,-1])
  tab_rf <- table(pred2, test$Cluster)
  #Tasa aciertos rf
  tasa.aciertos.rf[i]<-sum(tab_rf[row(tab_rf)==col(tab_rf)])/sum(tab_rf)
  sensitivity.rf[i] <- sensitivity(tab_rf)
  specificity.rf[i] <- specificity(tab_rf)
  
  ##Suport Vector Machines (kernlab)
  svp <- ksvm(factor(train_balanced_$Cluster)~ ., data=train_balanced_[,-44], type = "C-svc", kernel = "rbfdot",kpar = "automatic")
  pred3 <- predict(svp, test[,-1]) 
  tab_svp <- table(pred3, test$Cluster) 
  #Tasa aciertos svp
  tasa.aciertos.svp[i]<-sum(tab_svp[row(tab_svp)==col(tab_svp)])/sum(tab_svp)
  sensitivity.svp[i] <- sensitivity(tab_svp)
  specificity.svp[i] <- specificity(tab_svp)
  
  ##Arbol decision podado
  arbol<-rpartXse(train_balanced_$Cluster ~ ., se = 0.5, data=train_balanced_[,-44]) #arbol podado segun la regla 1-SE
  pred4 <- predict(arbol, test[,-1],type = "class")
  tab_tree <- table(pred4, test$Cluster) 
  #Tasa aciertos arbol
  tasa.aciertos.tree[i]<-sum(tab_tree[row(tab_tree)==col(tab_tree)])/sum(tab_tree)
  sensitivity.tree[i] <- sensitivity(tab_tree)
  specificity.tree[i] <- specificity(tab_tree)
  
  ##Regresion Logistica
  regresion <- glm(formula = train_balanced_$Cluster ~ ., family = binomial(link = "logit"), data = train_balanced_[,-44])
  pred5<-predict(regresion, test[,-1], type = "response")
  round(pred5)
  tab_regresion <- table(round(pred5), test$Cluster)
  #Tasa aciertos regresion
  tasa.aciertos.regresion[i]<-sum(tab_regresion[row(tab_regresion)==col(tab_regresion)])/sum(tab_regresion)
  sensitivity.regresion[i] <- sensitivity(tab_regresion)
  specificity.regresion[i] <- specificity(tab_regresion)
  
  ##Vecino mas proximo
  ctrl <- trainControl(method="cv",number = 10) 
  knnFit <- train(Cluster ~ ., data = train_balanced_, method = "knn",
                  trControl = ctrl, preProcess = c("center","scale"), tuneLength = 20)
  pred6 = predict(knnFit, test, type = "raw")
  tab_knn <- table(pred6, test$Cluster)
  tasa.aciertos.knn[i]<-sum(tab_knn[row(tab_knn)==col(tab_knn)])/sum(tab_knn)
  sensitivity.knn[i] <- sensitivity(tab_knn)
  specificity.knn[i] <- specificity(tab_knn)
}

######
#Evaluacion Modelos Tasa aciertos
######

#Medias del hold-out repetido Tasa Aciertos
mae.todos<-cbind(tasa.aciertos.plsda1, tasa.aciertos.regresion, tasa.aciertos.rf, tasa.aciertos.svp, tasa.aciertos.tree, tasa.aciertos.knn)
mae.todos.mean<-apply(mae.todos,2, mean)
mae.todos.mean
mae.todos.int<-apply(mae.todos,2, quantile, probs=c(0.025,0.975))
mae.todos.int

#Comparacion medias Tasa Aciertos
vector.mae<-as.vector(mae.todos)
#al convertirlo en vector pone una columna debajo de la otra
cbind(vector.mae[1:100], mae.todos[,1])
modelo<-rep(c("PLS-DA", "RL", "RF", "SVM", "DT", "KNN"), each=100)
bloque<-factor(rep(1:100,6))
plot(factor(modelo),vector.mae, xlab="Modelo", ylab="Tasa de Acierto")

#ANOVA
summary(aov(vector.mae ~ modelo))
fm<-aov(vector.mae ~ modelo + bloque)
summary(fm)
intervals = TukeyHSD(fm, "modelo", las=2, cex.axis=0.35)
plot(intervals, las=2, cex.axis=0.6)


#Confusion Matrix-Especificidad-Sensibilidad
confusionMatrix(as.factor(pred0), as.factor(test$Cluster)) #pls1
confusionMatrix(pred2, as.factor(test$Cluster)) #rf
confusionMatrix(pred3, as.factor(test$Cluster)) #svp
confusionMatrix(pred4, as.factor(test$Cluster)) #arbol
confusionMatrix(as.factor(round(pred5)), as.factor(test$Cluster))#regresion
confusionMatrix(as.factor(pred6), as.factor(test$Cluster))#KNN

####
#CURVAS ROC/AUC
###
ROCit_pls1 <- rocit(score=as.numeric(pred0), class=test$Cluster)
ROCit_pls1$AUC
ROCit_rf <- rocit(score=as.numeric(pred2), class=test$Cluster)
ROCit_rf$AUC
ROCit_svm <- rocit(score=as.numeric(pred3), class=test$Cluster)
ROCit_svm$AUC
ROCit_dt <- rocit(score=as.numeric(pred4), class=test$Cluster)
ROCit_dt$AUC
ROCit_rl <- rocit(score=as.numeric(round(pred5)), class=test$Cluster)
ROCit_rl$AUC
ROCit_knn <- rocit(score=as.numeric(pred6), class=test$Cluster)
ROCit_knn$AUC

#PLOT
par(mfrow=c(3,2))
plot(ROCit_pls1)
text(0.75, 0.5, c("AUC=0.73"), cex=1.2)
text(0.4, 0.2, c("PLS-DA"), cex=1.5, col="red")
plot(ROCit_knn)
text(0.75, 0.5, c("AUC=0.78"), cex=1.2)
text(0.4, 0.2, c("KNN"), cex=1.5, col="red")
plot(ROCit_rf)
text(0.75, 0.5, c("AUC=0.74"), cex=1.2)
text(0.4, 0.2,c("RF"), cex=1.5, col="red")
plot(ROCit_svm)
text(0.75, 0.5, c("AUC=0.68"), cex=1.2)
text(0.4, 0.2, c("SVM"), cex=1.5, col="red")
plot(ROCit_dt)
text(0.75, 0.5, c("AUC=0.80"), cex=1.2)
text(0.4, 0.2, c("DT"), cex=1.5, col="red")
plot(ROCit_rl)
text(0.75, 0.5, c("AUC=0.80"), cex=1.2)
text(0.4, 0.2, c("R.Logistica"), cex=1.5, col="red")


#RF
varImpPlot(rf, main="",col="dark blue")

summary(regresion)
###########################################
#Modelos cribados - variables reducidas
##########################################

huella_cribado <- read_excel("tfmobjetivo2.xlsx")
huella_cribado$Cluster <- factor(huella_cribado$Cluster)


#Hold out repetido modelos clasificacion

tasa.aciertos.rf1 <- NA
tasa.aciertos.regresion1<-NA
tasa.aciertos.svp1<-NA
tasa.aciertos.plsda11<-NA
tasa.aciertos.tree1<-NA
tasa.aciertos.knn1<-NA
sensitivity.plsda11 <- NA
specificity.plsda11 <- NA
sensitivity.tree1 <- NA
specificity.tree1 <- NA
sensitivity.rf1 <- NA
specificity.rf1 <- NA
sensitivity.regresion1 <- NA
specificity.regresion1 <- NA
sensitivity.svp1 <- NA
specificity.svp1 <- NA
sensitivity.knn1 <- NA
specificity.knn1 <- NA

for(i in 1:100){
  library(smotefamily)
  cjo1<-sample(1:nrow(huella_cribado),nrow(huella_cribado)*0.75, FALSE) 
  train_cribado <-huella_cribado[cjo1,]
  test_cribado<-huella_cribado[-cjo1,]
  train_balanced__ <- smotefamily::SMOTE(train_cribado[,-1],train_cribado[,1],K=3)
  train_balanced_cribado <- (train_balanced__$data)
  colnames(train_balanced_cribado)[19] <- "Cluster"
  train_balanced_cribado$Cluster <- ifelse(train_balanced_cribado$Cluster == 1, 0, 1)
  test_cribado$Cluster<- ifelse(test_cribado$Cluster == 1, 0, 1)
  test_cribado$Cluster <- factor(test_cribado$Cluster)
  train_balanced_cribado$Cluster <- factor(train_balanced_cribado$Cluster)
  
  #PLSDA
  plsda11 <- plsda(train_balanced_cribado[,-19], as.factor(train_balanced_cribado$Cluster), ncomp = 3, cv = 50)
  pred01 <- predict(plsda11, test_cribado[,-1], type = "class",ncomp = 3)
  tab_plsda1 <- table(pred01, test_cribado$Cluster) 
  #Tasa aciertos plsda1
  tasa.aciertos.plsda11[i]<-sum(tab_plsda1[row(tab_plsda1)==col(tab_plsda1)])/sum(tab_plsda1)
  sensitivity.plsda11[i] <- sensitivity(tab_plsda1)
  specificity.plsda11[i] <- specificity(tab_plsda1)
  
  
  ##Random forest
  rf1<-randomForest(factor(train_balanced_cribado$Cluster)~ .,data=as.data.frame(train_balanced_cribado[,-19]),mtry=4,method="class",importance=TRUE, ntree=5)
  pred21 <- predict(rf1, newdata=test_cribado[,-1])
  tab_rf1 <- table(pred21, test_cribado$Cluster)
  #Tasa aciertos rf
  tasa.aciertos.rf1[i]<-sum(tab_rf1[row(tab_rf1)==col(tab_rf1)])/sum(tab_rf1)
  sensitivity.rf1[i] <- sensitivity(tab_rf1)
  specificity.rf1[i] <- specificity(tab_rf1)
  
  ##Suport Vector Machines (kernlab)
  svp1 <- ksvm(factor(train_balanced_cribado$Cluster)~ ., data=train_balanced_cribado[,-19], type = "C-svc", kernel = "rbfdot",kpar = "automatic")
  pred31 <- predict(svp1, test_cribado[,-1]) 
  tab_svp1 <- table(pred31, test_cribado$Cluster) 
  #Tasa aciertos svp
  tasa.aciertos.svp1[i]<-sum(tab_svp1[row(tab_svp1)==col(tab_svp1)])/sum(tab_svp1)
  sensitivity.svp1[i] <- sensitivity(tab_svp1)
  specificity.svp1[i] <- specificity(tab_svp1)
  
  ##Arbol decision podado
  arbol1<-rpartXse(train_balanced_cribado$Cluster ~ ., se = 0.5, data=train_balanced_cribado[,-19]) #arbol podado segun la regla 1-SE
  pred41 <- predict(arbol1, test_cribado[,-1],type = "class")
  tab_tree1 <- table(pred41, test_cribado$Cluster) 
  #Tasa aciertos arbol
  tasa.aciertos.tree1[i]<-sum(tab_tree1[row(tab_tree1)==col(tab_tree1)])/sum(tab_tree1)
  sensitivity.tree1[i] <- sensitivity(tab_tree1)
  specificity.tree1[i] <- specificity(tab_tree1)
  
  ##Regresion Logistica
  regresion1 <- glm(formula = train_balanced_cribado$Cluster ~ ., family = binomial(link = "logit"), data = train_balanced_cribado[,-19])
  #model <- glm(response ~ ., data = train, family = "binomial")
  pred51<-predict(regresion1, test_cribado[,-1], type = "response")
  round(pred51)
  tab_regresion1 <- table(round(pred51), test_cribado$Cluster)
  #Tasa aciertos regresion
  tasa.aciertos.regresion1[i]<-sum(tab_regresion1[row(tab_regresion1)==col(tab_regresion1)])/sum(tab_regresion1)
  sensitivity.regresion1[i] <- sensitivity(tab_regresion1)
  specificity.regresion1[i] <- specificity(tab_regresion1)
  
  ##Vecino mas proximo
  ctrl1 <- trainControl(method="cv",number = 10) 
  knnFit1 <- train(Cluster ~ ., data = train_balanced_cribado, method = "knn",
                  trControl = ctrl1, preProcess = c("center","scale"), tuneLength = 20)
  pred61 = predict(knnFit1, test_cribado, type = "raw")
  tab_knn1 <- table(pred61, test_cribado$Cluster)
  tasa.aciertos.knn1[i]<-sum(tab_knn1[row(tab_knn1)==col(tab_knn1)])/sum(tab_knn1)
  sensitivity.knn1[i] <- sensitivity(tab_knn1) 
  specificity.knn1[i] <- specificity(tab_knn1)
}

######
#Evaluacion Modelos Tasa aciertos
######


#Medias del hold-out repetido Tasa Aciertos
mae.todos1<-cbind(tasa.aciertos.plsda11, tasa.aciertos.regresion1, tasa.aciertos.rf1, tasa.aciertos.svp1, tasa.aciertos.tree1, tasa.aciertos.knn1)
mae.todos.mean1<-apply(mae.todos1,2, mean)
mae.todos.mean1
mae.todos.int1<-apply(mae.todos1,2, quantile, probs=c(0.025,0.975))
mae.todos.int1

#Comparacion medias Tasa Aciertos
vector.mae1<-as.vector(mae.todos1)
#al convertirlo en vector pone una columna debajo de la otra
cbind(vector.mae1[1:100], mae.todos1[,1])
modelo1<-rep(c("PLS-DA", "RL", "RF", "SVM", "DT", "KNN"), each=100)
bloque1<-factor(rep(1:100,6))
plot(factor(modelo1),vector.mae1, xlab="Modelo", ylab="Tasa de Acierto")

#ANOVA
summary(aov(vector.mae1 ~ modelo1))
fm1<-aov(vector.mae1 ~ modelo1 + bloque1)
summary(fm1)
intervals1 = TukeyHSD(fm, "modelo", las=2, cex.axis=0.35)
plot(intervals1, las=2, cex.axis=0.6)


#Confusion Matrix-Especificidad-Sensibilidad
confusionMatrix(as.factor(pred01), as.factor(test_cribado$Cluster)) #pls1
confusionMatrix(pred21, as.factor(test_cribado$Cluster)) #rf
confusionMatrix(pred31, as.factor(test_cribado$Cluster)) #svp
confusionMatrix(pred41, as.factor(test_cribado$Cluster)) #arbol
confusionMatrix(as.factor(round(pred51)), as.factor(test_cribado$Cluster))#regresion
confusionMatrix(as.factor(pred61), as.factor(test_cribado$Cluster))#KNN


####
#CURVAS ROC/AUC
###
ROCit_pls1 <- rocit(score=as.numeric(pred01), class=test_cribado$Cluster)
ROCit_pls1$AUC
ROCit_rf <- rocit(score=as.numeric(pred21), class=test_cribado$Cluster)
ROCit_rf$AUC
ROCit_svm <- rocit(score=as.numeric(pred31), class=test_cribado$Cluster)
ROCit_svm$AUC
ROCit_dt <- rocit(score=as.numeric(pred41), class=test_cribado$Cluster)
ROCit_dt$AUC
ROCit_rl <- rocit(score=as.numeric(round(pred51)), class=test_cribado$Cluster)
ROCit_rl$AUC
ROCit_knn <- rocit(score=as.numeric(pred61), class=test_cribado$Cluster)
ROCit_knn$AUC

#PLOT
par(mfrow=c(3,2))
plot(ROCit_pls1)
text(0.75, 0.5, c("AUC=0.61"), cex=1.2)
text(0.4, 0.2, c("PLS-DA"), cex=1.5, col="red")
plot(ROCit_knn)
text(0.75, 0.5, c("AUC=0.69"), cex=1.2)
text(0.4, 0.2, c("KNN"), cex=1.5, col="red")
plot(ROCit_rf)
text(0.75, 0.5, c("AUC=0.75"), cex=1.2)
text(0.4, 0.2,c("RF"), cex=1.5, col="red")
plot(ROCit_svm)
text(0.75, 0.5, c("AUC=0.67"), cex=1.2)
text(0.4, 0.2, c("SVM"), cex=1.5, col="red")
plot(ROCit_dt)
text(0.75, 0.5, c("AUC=0.80"), cex=1.2)
text(0.4, 0.2, c("DT"), cex=1.5, col="red")
plot(ROCit_rl)
text(0.75, 0.5, c("AUC=0.80"), cex=1.2)
text(0.4, 0.2, c("R.Logistica"), cex=1.5, col="red")






















