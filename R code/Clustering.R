1. Fit Lee-Carter model on each country to extract the time trend
overall_lc_kapa_v3 <- matrix(0, nrow = 260, ncol = length(read_excel_name_v3))
for(i in 1:length(read_excel_name_v3)){
  overall_lc_kapa_v3[,i] <- lc(log(data_list_v3[[i]][,3:6]))$kapa
}
lc_kapa_ts_v3 <- ts(overall_lc_kapa_v3, frequency = 52)



2. Clutsering Method 1 (based on raw time trends):
raw_lc_kapa_ts_v3 <- t(lc_kapa_ts_v3)
raw_lc_kapa_ts_NH_v3 <- raw_lc_kapa_ts_v3[-c(1,18),]


# Find the optimal no. of cluster
## Elbow method
set.seed(1019)
wcss <- numeric(10)
for (k in 1:10) {
  kmeans_model <- kmeans(raw_lc_kapa_ts_NH_v3, centers = k, nstart = 25)
  wcss[k] <- kmeans_model$tot.withinss
}

plot(1:10, wcss, type = "b", pch = 19, xlab = "Number of Clusters", 
     ylab = "Within-cluster sum of squares", main = "Elbow Method")



## Silhouette Method
sil_width <- numeric(10)
for (k in 2:10) {
  kmeans_model <- kmeans(raw_lc_kapa_ts_NH_v3, centers = k, nstart = 25)
  sil <- silhouette(kmeans_model$cluster, dist(raw_lc_kapa_ts_NH_v3))
  sil_width[k] <- mean(sil[, 3])
}

plot(2:10, sil_width[-1], type = "b", pch = 19, xlab = "Number of Clusters", 
     ylab = "Average Silhouette Width", main = "Silhouette Method")


# Apply the optimal no. of cluster
raw_lc_kapa_ts_NH_v3_clusters  <- raw_lc_kapa_ts_NH_v3  %>%
  kmeans(centers = 2, nstart = 25)

raw_lc_kapa_ts_NH_v3_clusters$cluster["Series 1"] <- 3
raw_lc_kapa_ts_NH_v3_clusters$cluster["Series 18"] <- 3
ordered_raw_lc_kapa_ts_NH_v3 <- raw_lc_kapa_ts_NH_v3_clusters$cluster[order(as.numeric(gsub("Series ", "", names(raw_lc_kapa_ts_NH_v3_clusters$cluster))))]

# cluster result
lapply(1:3, \(i) read_excel_name_v3[ordered_raw_lc_kapa_ts_NH_v3 == i])




3. Clustering Method 2 (based on time trend slopes):
lc_kapa_stl_trend_v3 <- matrix(0, nrow = length(read_excel_name_v3), ncol = 1)
for(i in 1:length(read_excel_name_v3)){
  trend <- stl(lc_kapa_ts_v3[,i], s.window="periodic")$time.series[,"trend"]
  fit.trend <- data.frame(trend=trend, t=1:length(trend))
  lc_kapa_stl_trend_v3[i,1] <- lm(trend ~ t, data = fit.trend)$coefficients[2]
}

lc_kappa_trend_slope <- lc_kapa_stl_trend_v3

# Find the optimal no. of cluster
## Elbow method
set.seed(1019)
wcss <- numeric(10)
for (k in 1:10) {
  kmeans_model <- kmeans(lc_kapa_stl_trend_v3, centers = k, nstart = 25)
  wcss[k] <- kmeans_model$tot.withinss
}

plot(1:10, wcss, type = "b", pch = 19, xlab = "Number of Clusters", 
     ylab = "Within-cluster sum of squares", main = "Elbow Method")


## Silhouette Method
sil_width <- numeric(10)
for (k in 2:10) {
  kmeans_model <- kmeans(lc_kapa_stl_trend_v3, centers = k, nstart = 25)
  sil <- silhouette(kmeans_model$cluster, dist(lc_kapa_stl_trend_v3))
  sil_width[k] <- mean(sil[, 3])
}

plot(2:10, sil_width[-1], type = "b", pch = 19, xlab = "Number of Clusters", 
     ylab = "Average Silhouette Width", main = "Silhouette Method")


# Apply the optimal no. of cluster
lc_kapa_stl_trend_v3 <- lc_kapa_stl_trend_v3 %>%
  kmeans(centers = 3, nstart=25)

# cluster result
lapply(1:3, \(i) read_excel_name_v3[lc_kapa_stl_trend_v3$cluster == i])


