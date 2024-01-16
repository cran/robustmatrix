## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- results='hide', message = FALSE, warning=FALSE--------------------------
library(robustmatrix)
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(ggrepel)
library(tidyr)
library(tibble)
library(forcats)

## -----------------------------------------------------------------------------
data(weather)
X <- weather[,,!apply(weather,3, function(x) any(is.na(x)))]
n<- dim(X)[3]
p<- dim(X)[1]
q<- dim(X)[2]

## -----------------------------------------------------------------------------
par_MMLE <- mmle(X, lambda = 0)
set.seed(1)
par_MMCD <- mmcd(X, alpha = 0.5 ,lambda = 0, nsamp = 5000, nthreads = 1)

## -----------------------------------------------------------------------------
MD <- as.numeric(mmd(X = X, 
                     mu = par_MMLE$mu, 
                     cov_row = par_MMLE$cov_row_inv, 
                     cov_col = par_MMLE$cov_col_inv, 
                     inverted = TRUE))
MD_rob <- as.numeric(mmd(X = X, 
                         mu = par_MMCD$mu, 
                         cov_row = par_MMCD$cov_row_inv, 
                         cov_col = par_MMCD$cov_col_inv, 
                         inverted = TRUE))
out_quant <- qchisq(0.95, p*q)
outliers <- which(MD_rob > out_quant)

## ---- fig.dim = c(7,4)--------------------------------------------------------
names(MD_rob) <- dimnames(X)[[3]]
subs <- MD_rob < out_quant
plt_dd <- ggplot(data = NULL, aes(x = sqrt(MD), y = sqrt(MD_rob))) +
  geom_point() +
  geom_hline(yintercept = sqrt(out_quant)) +
  geom_vline(xintercept = sqrt(out_quant)) +
  geom_label_repel(aes(x = sqrt(MD[!subs]), 
                       y = sqrt(MD_rob[!subs]), 
                       label = names(MD_rob[!subs])), 
                   size = 2.5, max.overlaps = 200, nudge_y = 0.5) +
  theme_classic() +
  xlim(c(min(sqrt(MD)),max(sqrt(MD))+1)) +
  labs(x = "MMD", y = "Robust MMD")
plt_dd

## -----------------------------------------------------------------------------
shv_cell <- array(dim = dim(X))
shv_cell[] <- matrixShapley(X, 
                            mu = par_MMCD$mu, 
                            cov_row = par_MMCD$cov_row_inv, 
                            cov_col = par_MMCD$cov_col_inv, 
                            inverted = TRUE, type = "cell")
dimnames(shv_cell) <- dimnames(X)

## -----------------------------------------------------------------------------
long_weather <- function(x, val_to, var_levels = NULL) {
  y <- x %>%
    rownames_to_column(var = "var") %>%
    pivot_longer(cols = -var, names_to = "parameter", values_to = val_to)
  if(!is.null(var_levels)){
    y <- y %>%
      mutate(var = factor(var, levels = var_levels))
  }
  y
}

var_names <- c("T [°C]" = "T",
               "SH [h]" = "SH",
               "P [mm/m²]" = "P",
               "AP [hPa]" = "AP",
               "SP [%]" = "SP")

set <- cbind("var" = 1891:2022,
             "Temperature [°C]" = 0,
             "Preciptiation [mm/m²]" = 0,
             "Solid Preciptiation [%]" = 0,
             "Sunshine hours [h]" = 0,
             "Air pressure [hPa]" = 0)
set[!(set[,1] %in% dimnames(X)[[3]]),-1] <- 1
set <- data.frame(set)
colnames(set) <- c("var",
                   "T [°C]",
                   "P [mm/m²]",
                   "SP [%]",
                   "SH [h]",
                   "AP [hPa]")
set_long <- set %>% pivot_longer(-c(var))

## ---- message = FALSE, warning=FALSE------------------------------------------
shv_year <- apply(shv_cell, c(1,3),sum)
X_center <- array(dim = dim(X), dimnames = dimnames(X))
X_center[] <- (apply(X,3,function(x) x - par_MMCD$mu))
X_center_sign <- array(dim = dim(X), dimnames = dimnames(X))
X_center_sign[] <- sign(X_center)
sign_mat <- apply(X_center,c(1,3),function(x) sign(mean(x)))
sign_mat[,-outliers] <- 0
shv_rel <- apply(shv_year, 2, function(x) x/sum(x))
shv_rel_outliers <- shv_rel
shv_rel_outliers[,-outliers] <- 0
shv_rel_outliers[shv_rel_outliers < 0] <- 0 # ONLY PLOT POSITIVE SHAPLEY VALUES


shv_years_long <- long_weather(data.frame(t(shv_rel_outliers)), val_to = "shv_rel")
sign_years_long <- long_weather(data.frame(t(sign_mat)), val_to = "sign")
plt_data_years <- inner_join(shv_years_long,sign_years_long) %>% 
  mutate(val = sign*abs(shv_rel), 
         sign = as.character(sign), 
         var = as.numeric(var),
         parameter = fct_recode(factor(parameter),!!! var_names))

## ---- message = FALSE, warning=FALSE------------------------------------------
shv_long_1895 <- long_weather(data.frame(t(shv_cell[,,outliers[1]])), val_to = "shv", var_levels = dimnames(X)[[2]])
X_long_1895 <- long_weather(data.frame(t(X[,,outliers[1]])) %>% mutate(SP = SP*100), val_to = "x", var_levels = dimnames(X)[[2]])
sign_long_1895 <- long_weather(data.frame(t(sign(X[,,outliers[1]] - par_MMCD$mu))), val_to = "sign", var_levels = dimnames(X)[[2]])
plt_data_1895 <- inner_join(shv_long_1895,sign_long_1895) %>% 
  mutate(val = sign*abs(shv), 
         sign = as.character(sign)) %>%
  inner_join(X_long_1895) %>%
  mutate(parameter = fct_recode(factor(parameter),!!! var_names))

shv_long_2022 <- long_weather(data.frame(t(shv_cell[,,n])), val_to = "shv", var_levels = dimnames(X)[[2]])
X_long_2022 <- long_weather(data.frame(t(X[,,n])) %>% mutate(SP = SP*100), val_to = "x", var_levels = dimnames(X)[[2]])
sign_long_2022 <- long_weather(data.frame(t(sign(X[,,n] - par_MMCD$mu))), val_to = "sign", var_levels = dimnames(X)[[2]])
plt_data_2022 <- inner_join(shv_long_2022,sign_long_2022) %>% 
  mutate(val = sign*abs(shv), 
         sign = as.character(sign)) %>%
  inner_join(X_long_2022) %>%
  mutate(parameter = fct_recode(factor(parameter),!!! var_names))

## ---- fig.dim = c(7,4), message = FALSE, warning=FALSE------------------------
ggplot(plt_data_years, aes(x = var, y = parameter, fill = val)) +
  geom_tile(color = "lightgray") +
  theme_minimal() +
  scale_fill_gradient2(low = "#1F78B4", mid = "white", high = "#E31A1C", 
                       midpoint = 0, guide = "none")  +
  scale_x_continuous(breaks = seq(from = 1880, to = 2020, by = 5), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1),
        axis.title = element_blank(), title = element_blank()) +
  geom_tile(data = set_long, aes(x = var, y = name, alpha = factor(value)), fill = "gray") +
  scale_alpha_manual(values = c(0,1)) +
  guides(alpha = FALSE, fill = FALSE)

## ---- fig.dim = c(7,4), message = FALSE, warning=FALSE------------------------
gridExtra::grid.arrange(ggplot(plt_data_1895, aes(x = var, y = parameter, fill = val, label = round(x))) +
                          geom_tile(color = "lightgray") +
                          theme_minimal() +
                          geom_text(size = 3) +
                          theme(axis.text.x = element_text(angle = 45, vjust = 1.4, hjust=1.4),
                                axis.title = element_blank())+
                          labs(title = "1895") + 
                          scale_fill_gradient2(low = "#1F78B4", mid = "white", high = "#E31A1C", 
                                               midpoint = 0, guide = "none"),
                        ggplot(plt_data_2022, aes(x = var, y = parameter, fill = val, label = round(x))) +
                          geom_tile(color = "lightgray") +
                          theme_minimal() +
                          geom_text(size = 3) +
                          theme(axis.text.x = element_text(angle = 45, vjust = 1.4, hjust=1.4),
                                axis.title = element_blank())+
                          labs(title = "2022") + 
                          scale_fill_gradient2(low = "#1F78B4", mid = "white", high = "#E31A1C", 
                                               midpoint = 0, guide = "none"),
                        nrow = 1)

## -----------------------------------------------------------------------------
data(darwin)
darwin_class <- dimnames(darwin)[[3]]
sub_H <- which(darwin_class == "H")
sub_P <- which(darwin_class == "P")

## -----------------------------------------------------------------------------
darwin_array_clean <- darwin[-c(1,9,18),,]
X_H <- darwin_array_clean[,,sub_H] 
X_P <- darwin_array_clean[,,sub_P] 
dimnames(X_H)[[3]] <- sub_H
dimnames(X_P)[[3]] <- sub_P

p <- dim(X_H)[1]
q <- dim(X_H)[2]
n <- dim(X_H)[3]

## -----------------------------------------------------------------------------
set.seed(1)
par_mmle <- robustmatrix::mmle(X_H)
par_mmcd <- robustmatrix::mmcd(X_H, nsamp = 500, alpha = 0.5, 
                               nthreads = 1, scale_consistency = "quant")

## ---- fig.dim = c(7,4)--------------------------------------------------------
plt_data <-   data.frame(patient_id = 1:dim(darwin_array_clean)[[3]], 
                         class = c(rep("AD", length(sub_P)),
                                   rep("H", length(sub_H))),
                         mmd = sqrt(mmd(darwin_array_clean, 
                                        par_mmcd$mu, par_mmcd$cov_row, 
                                        par_mmcd$cov_col)))

p_mmd <- ggplot(data = plt_data, 
                mapping = aes(x = patient_id, y = mmd, color = class)) + 
  geom_point(show.legend = FALSE) + 
  geom_hline(yintercept = sqrt(qchisq(0.975, p*q))) + 
  labs(x = "Subject ID", y = "MMD", color = element_blank()) + 
  theme_classic() +
  facet_grid(~ class, scales = "free")
p_mmd

## ---- message = FALSE, warning=FALSE------------------------------------------
shv_cell_P <- array(dim = dim(X_P), dimnames = dimnames(X_P))
shv_cell_P[] <- matrixShapley(X_P, 
                              mu = par_mmcd$mu, 
                              cov_row = par_mmcd$cov_row_inv, 
                              cov_col = par_mmcd$cov_col_inv, 
                              inverted = TRUE, type = "cell")

shv_cell_H <- array(dim = dim(X_H), dimnames = dimnames(X_H))
shv_cell_H[] <- matrixShapley(X_H, 
                              mu = par_mmcd$mu, 
                              cov_row = par_mmcd$cov_row_inv, 
                              cov_col = par_mmcd$cov_col_inv, 
                              inverted = TRUE, type = "cell")


shv_row_mean_P <- data.frame("val" = apply(shv_cell_P,c(1),
                                           function(x)mean(x))) %>% 
  rownames_to_column()
shv_row_mean_H <- data.frame("val" = apply(shv_cell_H,c(1),
                                           function(x)mean(x))) %>% 
  rownames_to_column()

shv_row_mean <- rbind(cbind("class" = "AD", shv_row_mean_P),
                      cbind("class" = "H", shv_row_mean_H))


shv_row_mean_P_prop <- data.frame("val" = apply(shv_cell_P,c(1),
                                                function(x)sum(x)/sum(shv_cell_P))) %>% 
  rownames_to_column()
shv_row_mean_H_prop <- data.frame("val" = apply(shv_cell_H,c(1),
                                                function(x)sum(x)/sum(shv_cell_H))) %>% 
  rownames_to_column()

shv_row_mean_prop<- rbind(cbind("class" = "AD", shv_row_mean_P_prop),
                          cbind("class" = "H", shv_row_mean_H_prop))
fct_order <- order(shv_row_mean_P_prop$val)
shv_row_mean_prop <- shv_row_mean_prop %>%
  mutate(rowname = factor(rowname, levels = shv_row_mean_P_prop$rowname[fct_order]))

## ---- fig.dim = c(7,4)--------------------------------------------------------
p_shv_features <- ggplot(data = shv_row_mean_prop, 
                         aes(x = val, y = rowname, fill = class)) + 
  geom_bar(stat = "identity", position=position_dodge(), color = "black") + 
  labs(x = "Proportional contribution to MMD²", 
       y = element_blank(), fill = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme_classic() + 
  theme(legend.position = c(0.95,0.5))
p_shv_features

## ---- eval=FALSE--------------------------------------------------------------
#  # Loading the data
#  load(url("https://wis.kuleuven.be/stat/robust/Programs/DO/do-video-data-rdata"))
#  
#  # Creating the overview plot
#  ind_overview <- c(1,487,491,495,500)
#  plts_overview <- list()
#  for(i in seq_along(ind_overview)){
#    dim_Video <- dim(Video)
#    video_grid <- expand.grid("x" = 1:dim_Video[2], "y" = 1:dim_Video[3])
#    video_long <- data.frame(cbind(video_grid,
#                                   "r" = as.vector(Video[ind_overview[i],,,1]/255),
#                                   "g" = as.vector(Video[ind_overview[i],,,2]/255),
#                                   "b" = as.vector(Video[ind_overview[i],,,3]/255))) %>% mutate(rgb.val=rgb(r,g,b))
#    plts_overview[[i]] <- ggplot(video_long, aes(y,x)) +
#      geom_raster(aes(fill=rgb.val)) +
#      scale_fill_identity() +
#      theme_void()+
#      coord_fixed() +
#      theme(plot.title = element_text(hjust = 0.5),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            panel.background = element_rect(fill = "transparent" ,colour = "black"),
#            panel.ontop = TRUE,
#            plot.margin = margin(1, 1, 1, 1, "pt")) +
#      scale_x_continuous(expand = c(0,0)) +
#      scale_y_continuous(expand = c(0,0), trans = "reverse") +
#      labs(title = paste("Frame", ind_overview[i]))
#  }
#  video_overview <- do.call(gridExtra::grid.arrange, args = list("grobs" = plts_overview, "nrow" = 1))
#  
#  # Transform to grayscale image
#  video_grayscale <- apply(Video, 1:3, mean)
#  X <- aperm(video_grayscale,c(2,3,1))
#  n <- dim(X)[3]
#  p <- dim(X)[1]
#  q <- dim(X)[2]
#  
#  # MMLE and MMCD parameter estimation
#  par_MMLE <- mmle(X, lambda = 0)
#  set.seed(1)
#  par_MMCD <- mmcd(X, alpha = 0.5,lambda = 0, nsamp = 500, nthreads = 1)
#  
#  # Compute squared Mahlanobis distances
#  MD <- mmd(X, mu = par_MMLE$mu, cov_row = par_MMLE$cov_row_inv,
#            cov_col = par_MMLE$cov_col_inv, inverted = TRUE)
#  MD_rob <- mmd(X, mu = par_MMCD$mu, cov_row = par_MMCD$cov_row_inv,
#                cov_col = par_MMCD$cov_col_inv, inverted = TRUE)
#  out_quant <- qchisq(0.99, p*q)
#  
#  # Plot MMD of all observations
#  data_plt_md <- data.frame(cbind("Frame" = 1:length(MD_rob), "MD" = sqrt(MD_rob)))
#  ggplot(data_plt_md, aes(x = Frame, y = MD)) +
#    geom_point() +
#    labs(x = "Frame", y = "MMD") +
#    scale_x_continuous(limits = range(1,n), breaks = seq(from = 50, to = n, by = 50), expand = c(0,3)) +
#    theme_classic()
#  
#  # Plot MMD of last observations
#  sub_ind <- (1:length(MD_rob))[-(1:475)]
#  frame_labels <- rep(NA,length(sub_ind))
#  frame_labels[sub_ind %in% c(487,491,495,500)] <- c("Frame 487: man left of tree", "Frame 491: man behind tree", "Frame 495: man right of tree", "Frame 500: man fully visible")
#  data_plt_md1 <- data.frame(cbind("Frame" = sub_ind, "MD" = sqrt(MD_rob[sub_ind]), "labels" = frame_labels))
#  ggplot(data_plt_md1, aes(x = as.numeric(Frame), y = as.numeric(MD), label = labels)) +
#    geom_point() +
#    geom_path() +
#    geom_label_repel(nudge_x = 75, nudge_y = -80, color = "darkblue") +
#    labs(x = "Frame", y = "MMD") +
#    scale_x_continuous(limits = range(sub_ind), breaks = seq(from = 480, to = 630, by = 10), expand = c(0,1)) +
#    theme_classic()
#  
#  # Compute Shapley values
#  shv <- array(dim = dim(X))
#  shv[,,] <- matrixShapley(X[,,], mu = par_MMCD$mu,
#                           cov_row = par_MMCD$cov_row_inv, cov_col = par_MMCD$cov_col_inv,
#                           inverted = TRUE, type = "cell")
#  
#  # Function to plot Shapley values for image data
#  plot_shapley_image <- function(image,shapley_values, lower = -100, upper = 100,
#                                 title = NULL, positive_only = FALSE, border = TRUE){
#    shapley_values[shapley_values < lower] <- lower
#    shapley_values[shapley_values > upper] <- upper
#    if(positive_only){
#      shapley_values[shapley_values < 0] <- 0
#    }
#  
#    image_grid <- expand.grid("x" = 1:nrow(shapley_values), "y" = 1:ncol(shapley_values))
#    image_data <- cbind(image_grid, "z" = as.vector(image), "shv" = as.vector(shapley_values))
#    plt <- ggplot()
#    if(any(image != 0)){
#      plt <- plt +
#        geom_raster(data = image_data, aes(x, y, fill = z)) +
#        scale_fill_gradient(low = "black", high = "white",guide=FALSE)
#    }
#    plt <- plt +
#      new_scale("fill") +
#      geom_raster(data = image_data, aes(x, y, fill = shv), na.rm = TRUE) +
#      scale_fill_gradientn(colours= c("blue", "transparent", "red"),limits=c(lower, upper),
#                           guide=FALSE, na.value = "transparent") +
#      labs(title = title) +
#      theme_void()+
#      coord_fixed() +
#      theme(plot.title = element_text(hjust = 0.5))
#    if(border){
#      plt <- plt  +
#        theme(panel.grid.major = element_blank(),
#              panel.grid.minor = element_blank(),
#              panel.background = element_rect(fill = "transparent" ,colour = "black"),
#              panel.ontop = TRUE,
#              plot.margin = margin(1, 1, 1, 1, "pt")) +
#        scale_y_continuous(expand = c(0,0)) +
#        scale_x_continuous(expand = c(0,0))
#    }
#    plt
#  }
#  
#  # Plot Shapley values for selected frames
#  lims <- 100
#  gridExtra::grid.arrange(
#    plot_shapley_image(image = t(apply(X[,,487],2,rev)), shapley_values = (t(apply(shv[,,487],2,rev))),
#                       lower = -lims, upper = lims, title = "Frame 487", positive_only = TRUE),
#    plot_shapley_image(image = t(apply(X[,,491],2,rev)), shapley_values = (t(apply(shv[,,491],2,rev))),
#                       lower = -lims, upper = lims, title = "Frame 491", positive_only = TRUE),
#    plot_shapley_image(image = t(apply(X[,,495],2,rev)), shapley_values = (t(apply(shv[,,495],2,rev))),
#                       lower = -lims, upper = lims, title = "Frame 495", positive_only = TRUE),
#    plot_shapley_image(image = t(apply(X[,,500],2,rev)), shapley_values = (t(apply(shv[,,500],2,rev))),
#                       lower = -lims, upper = lims, title = "Frame 500", positive_only = TRUE),
#    nrow = 1
#  )

