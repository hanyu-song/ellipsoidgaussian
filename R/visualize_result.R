#library(rgl)
#library(ggplot2)
#library(reshape2)
#library(dplyr)
#library(GGally)
#library(scales)
#library(grDevices)
#library(RColorBrewer)

#' ggplot2 pairs plot with correlation information
#'
#' @description
#' `pairsPlot_onedata` makes a matrix of plots for a given data set with correlation in the upper
#' triangular part, density on the diagonal, and scatter plots on the lower triangular part.
#'
#' @details
#' Arguments \code{dir} and \code{name} provide the option to save
#' the plot at location \code{dir} with name \code{name}.
#' \code{rotate_x_text} also allows users to rotate the texts on the
#' x-axis in case of space constraint.
#'
#' @param dat_new A matrix of data.
#' @param dir String, the path to where the plot should be saved. Default is NULL, meaning
#' the plot is not saved.
#' @param name String, the name of the plot to be saved. When \code{dir} is NULL, meaning
#' the plot will not be saved, the name can be left blank.
#' @param rotate_x_text Logical, whether to rotate the texts on the x-axis or not.
#' If so, they are positioned at a 60-degree angle.
#' @returns ggmatrix object from GGally, a matrix of plots
#' @export
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 theme ggsave
#' @examples
#' pairsPlot_onedata(shell, NULL)
pairsPlot_onedata <- function(dat_new, dir = NULL, name, rotate_x_text = F) {
  df <- data.frame(dat_new)
  #  df2 <- data.frame(dat_orig, type = 'o')
  #  df <- rbind(df, df2)
  num <- ncol(df)
  # select sub columns to better visualization if there are too many vars
 # if (num > 6) {
    # idx <- sample(1:num, 4)
    # df <- df %>% dplyr::select(eval(paste0('X',idx)), type)
   # num <- 6
 # }
  q <- GGally::ggpairs(df, columns = 1:num,ggplot2::aes(alpha = 0.5),
               #  upper = list(continuous = wrap(cor_func,
               #     method = 'spearman')),
               upper = list(continuous = GGally::wrap(cor_func,method = 'spearman', symbol = '\u03C1 =')),
               #     upper = list(continuous = ggally_density, combo = ggally_box_no_facet),
               #   lower = list(continuous = wrap('points',size = 0.3)),
               lower = list(continuous = myscatter),
               diag = list(continuous = mydensity, axisLabels='show'))
  if(rotate_x_text) {
    q <- q +  ggplot2::theme(axis.text.x = element_text(angle = 60))
  }
  #q <- q +
    #theme_minimal()
  #  theme(legend.position = "none",
       #   panel.grid.major = element_blank(),
        #  axis.ticks = element_blank())
  if (!is.null(dir)) {
    if(!dir.exists(dir)) {dir.create(dir)}
    ggplot2::ggsave(plot = q, filename = paste0(dir,'/', name,'.pdf'), width = 4.4, height = 4.4,device = cairo_pdf)
  }
  return(q)
}

#' ggplot2 pairs plot with two data sets in juxtaposition
#'
#' @description
#' `pairsPlot` makes a matrix of plots with two data sets in juxtaposition with
#' correlation in the upper triangular part, density on the diagonal, and
#' scatter plots on the lower triangular part.
#'
#' @details
#' Arguments \code{dir} and \code{name} provide the option to save
#' the plot at location \code{dir} with name \code{name}.
#'
#' @param data_orig A matrix of data, e.g. the original data set.
#' @param dat_new A matrix of data, e.g. the data generated from the posterior
#' predictive distribution.
#' @param dir String, the path to where the plot should be saved. Default is NULL, meaning
#' the plot is not saved.
#' @param name String, the name of the plot to be saved. When \code{dir} is NULL, meaning
#' the plot will not be saved, the name can be left blank.
#' @returns ggmatrix object from GGally, a matrix of plots
#' @export
#' @importFrom GGally ggpairs
#' @importFrom ggplot2 ggsave
#' @examples
#' dat_new <- gen_posterior_predictive(1000, samples, 2500)
#' pairsPlot(shell, dat_new, NULL)
pairsPlot <- function(dat_orig, dat_new, dir, name) {
  colnames(dat_new) <- colnames(dat_orig)
  df <- data.frame(dat_new,type = 'p')
    df2 <- data.frame(dat_orig, type = 'o')
    df <- rbind(df, df2)
  num <- ncol(df) -1
  # select sub columns to better visualization if there are too many vars
  if (num > 4) {
    # idx <- sample(1:num, 4)
    # df <- df %>% dplyr::select(eval(paste0('X',idx)), type)
    num <- 4
  }
  q <- GGally::ggpairs(df, columns = 1:num,ggplot2::aes(group = type,col = type,fill = type,alpha = 0.5),
               #  upper = list(continuous = wrap(cor_func,
               #     method = 'spearman')),
               upper = list(continuous = mycorrelations),
         #     upper = list(continuous = wrap("cor", title = '\u03C1')),
       #    upper = list(continuous = wrap(cor_func,
                         #         method = 'spearman', symbol = '\u03C1 =')),
               #     upper = list(continuous = ggally_density, combo = ggally_box_no_facet),
                #  lower = list(continuous = wrap('points',size = 0.3,alpha = 0.5)),
              lower = list(continuous = myscatter),
               diag = list(continuous = mydensity, axisLabels='show'))

  #q <- q +
   # theme_light() +
   # theme(legend.position = "none",
       #   panel.grid.major = element_blank(),
       #   panel.grid.minor = element_line(color = 'grey'))
       #   axis.ticks = element_blank())
       #  panel.background = element_blank())
  #plot.brackground = element_rect(fill = '#BFD5E3'))
          #element_rect(fill = "#BFD5E3"))
  if (!is.null(dir)) {
    if(!dir.exists(dir)) {dir.create(dir)}
    ggplot2::ggsave(plot = q, filename = paste0(dir,'/', name,'.pdf'), width = 4.4, height = 4.4,device = cairo_pdf)
  }
#  print(q)
  return(q)
}
#

#' @importFrom GGally ggally_text eval_data_col
#' @importFrom ggplot2 theme_classic
# Defines function to color according to correlation
cor_func <- function(data, mapping, method, symbol, ...){
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  corr <- cor(x, y, method=method, use='complete.obs')
  #colFn <- colorRampPalette(c("brown1", "white", "dodgerblue"),
  #                          interpolate ='spline')
 # fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]

  GGally::ggally_text(
    label = paste(symbol, as.character(round(corr, 2))),
    mapping = ggplot2::aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  ) + ggplot2::theme_classic()
  #+ #removed theme_void()
   # theme(panel.background = element_rect(fill = fill))
}


#### correlation
#' @importFrom dplyr bind_rows mutate group_by filter summarize %>% as_label
#' @importFrom ggplot2 ggplot geom_text theme_classic
mycorrelations <- function(data,mapping,...){
  data2 = data
  data2$x = as.numeric(data[,dplyr::as_label(mapping$x)])
  data2$y = as.numeric(data[,dplyr::as_label(mapping$y)])
  data2$group = data[,dplyr::as_label(mapping$colour)]

  correlation_df = data2 %>%
    dplyr::bind_rows(data2 %>% dplyr::mutate(group="Overall Corr")) %>%
    dplyr::group_by(group) %>%
    dplyr::filter(sum(!is.na(x),na.rm=T)>1) %>%
    dplyr::filter(sum(!is.na(y),na.rm=T)>1) %>%
    dplyr::summarize(estimate = round(as.numeric(cor.test(x,y,method="spearman")$estimate),2)) %>%
       #       pvalue = cor.test(x,y,method="spearman")$p.value,
         #     pvalue_star = as.character(symnum(pvalue, corr = FALSE, na = FALSE,
                                          #      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                            #    symbols = c("***", "**", "*", "'", " "))))%>%
    dplyr::group_by() %>%
    dplyr::mutate(group = factor(group, levels=c(as.character(unique(sort(data[,dplyr::as_label(mapping$colour)]))), "Overall Corr")))
  correlation_df <- correlation_df %>% dplyr::filter(!group %in% 'Overall Corr')
  ggplot2::ggplot(data=correlation_df, ggplot2::aes(x=1,y=group,color=group))+
    ggplot2::geom_text(ggplot2::aes(label=paste0(group,": ",estimate)), size = 6) + ggplot2::theme_classic()
}

#' @importFrom ggplot2 ggplot geom_point scale_x_continuous scale_y_continuous theme_minimal
#' @importFrom scales pretty_breaks
myscatter <- function(data, mapping, ...) {
  ggplot2::ggplot(data = data, mapping = mapping, ...) +
    ggplot2::geom_point(size = 0.3, alpha = 0.3,...) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +ggplot2::theme_minimal()
}
#' @importFrom ggplot2 ggplot geom_density scale_x_continuous scale_y_continuous theme_minimal
#' @importFrom scales pretty_breaks
mydensity <- function(data, mapping, ...) {
  ggplot2::ggplot(data = data, mapping = mapping, ...) +
    ggplot2::geom_density(alpha= 0.5,...) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=3)) + ggplot2::theme_minimal()

}

visualizeTwoData <- function(orig, new, col1, col2, dir, method) {
  # TO BE MODIFIED
  stopifnot(ncol(new)  == ncol(orig))
  df1 <- data.frame(new)
  df0 <- data.frame(orig)
  names(df1) <- names(df0) <- paste0('X', 1:ncol(df1))
  df1$label <- 'posterior pred.'
  df0$label <- 'original'
  df <- rbind(df0, df1)
  df$color <- ifelse(df$label %in% 'original', 'royalblue1','red')
  q <- ggplot2::ggplot(df) +
    ggplot2::geom_point(ggplot2::aes(x = get(paste0('X',col1)), y = get(paste0('X',col2)),col = label, group = label), alpha = 0.1) +
    ggplot2::theme_minimal() +
    ggplot2::xlab(paste0('X',col1)) + ggplot2::ylab(paste0('X',col2)) +
   # theme(legend.position="center right") +
    ggplot2::scale_color_manual(values=c('royalblue1', 'red')) +
    ggplot2::theme(legend.position="none") + ggplot2::coord_fixed()
  if (!is.null(dir)) {
    ggplot2::ggsave(plot = q, filename = paste0(dir, method,'_2dscatter.pdf'), width = 2.5, height = 4)
  }
  return(q)
}

findOrig <- function(vec) {
  mean(vec[length(vec)])
}

visualizeTwoData_3d <- function(dat_orig, dat_new, directory, name, axesName = NULL) {
  # axesName: a vector of names for x, y, and z axes
  data <- rbind(cbind(data.frame(dat_orig, group = "Original")), cbind(data.frame(dat_new, group = "Posterior Pred.")))


  # Add a new column with color
 # mycolors <- c('royalblue1', 'darkcyan', 'oldlace')
  data$color <- ifelse(data$group %in% 'Original', 'royalblue1','red')
  #  data$color <- ifelse(data$group %in% 'Original', 'red','green')
  # Plot
  if (!is.null(axesName)) {
   rgl::plot3d(
    x=data$X1, y=data$X2, z=data$X3,
    col = data$color,
    type = 'p',
    radius = .1,
    xlab= axesName[1], ylab= axesName[2], zlab= axesName[3], aspect = F)
  }
  else {
    rgl::plot3d(
      x=data$X1, y=data$X2, z=data$X3,
      col = data$color,
      type = 'p',
      radius = .1,
      xlab = '', ylab = '', zlab = '', aspect = F)
  }
  rgl::legend3d("top", legend =  c("Original", "Posterior Pred."), pch = 16, col =  c('royalblue1','red'), cex=1, inset=c(0.02))
  rgl::rgl.snapshot(filename = paste0(directory, name,'.png') , fmt = 'png')
  rgl::rgl.postscript(paste0(directory, name,'.pdf'),fmt="pdf")
  # To display in an R Markdown document:
  # rglwidget()

  # To save to a file:
  htmlwidgets::saveWidget(rglwidget(width = 520, height = 520),
                          file = paste0(directory,  name,'.html'),
                          libdir = "libs",
                          selfcontained = FALSE
  )
  return(q)

}
tracePlots <- function(resList,fname, dir2Save) {
  # trace plots and the associated ACF plots
  p <- dim(resList$EG$lambda)[2]
  Lambda_mu <- matrix(NA, nrow = nrow(resList$EG$tau), ncol = dim(resList$EG$lambda)[2])
  lambda2 <- matrix(NA, nrow(resList$EG$tau), ncol = p*(p+1) / 2)
  for (j in 1:nrow(Lambda_mu)) {
    Lambda_mu[j,] <- resList$EG$lambda[j,,] %*% resList$EG$mu[j,]
    mat <- tcrossprod(resList$EG$lambda[j,,])
    lambda2[j,] <- mat[lower.tri(mat, diag = T)]
  }
  df3 <- data.frame(tau = as.numeric(resList$EG$tau), iteration = 1:nrow(resList$EG$tau),
                    invsig1 = resList$EG$invsig[,1], lambda_mu = Lambda_mu[,1],
                    lambda2 = lambda2[,2])
  dirL <- paste0(dir2Save, '/lambda2/')
  dirLM <- paste0(dir2Save, '/lambda_mu/')
  dirTau <- paste0(dir2Save, '/tau/')
  dirNoise <- paste0(dir2Save, '/noise/')
  if (!dir.exists(dirL)) {
    dir.create(dirL)
  }
  if (!dir.exists(dirLM)) {
    dir.create(dirLM)
  }
  if (!dir.exists(dirTau)) {
    dir.create(dirTau)
  }
  if (!dir.exists(dirNoise)) {
    dir.create(dirNoise)
  }

  df3 <- df3 %>% dplyr::filter(iteration > floor(nrow(df3) / 2))
  q <- ggplot2::ggplot() + ggplot2::geom_line(data = df3, ggplot2::aes(x = iteration, y = lambda2)) + ggplot2::ylab(latex2exp::TeX("$(\\Lambda\\Lambda^T)_{12}$"))
  ggplot2::ggsave(q, filename = paste0(dirL, fname,'.pdf'), width = 4, height =4 )
  q <- ggplot2::ggplot() + ggplot2::geom_line(data = df3, ggplot2::aes(x = iteration, y = lambda_mu)) + ggplot2::ylab(latex2exp::TeX("$(\\Lambda\\mu)_1"))
  ggplot2::ggsave(q, filename = paste0(dirLM, fname,'.pdf'), width = 4, height =4 )
  q <- ggplot() + ggplot2::geom_line(data = df3, ggplot2::aes(x = iteration, y = tau)) + ggplot2::ylab(latex2exp::TeX("$\\tau$"))
  ggplot2::ggsave(q, filename = paste0(dirTau, fname,'.pdf'), width = 4, height =4 )
  q <- ggplot() + ggplot2::geom_line(data = df3, ggplot2::aes(x = iteration, y = invsig1))+   ggplot2::ylab(latex2exp::TeX("1/$\\sigma_1^2$"))
  ggplot2::ggsave(q, filename = paste0(dirNoise, fname,'.pdf'), width = 4, height =4 )
  pdf(file = paste0(dir2Save, '/lambda2/acf_', fname, '.pdf'), width = 4, height = 4)
  acf(df3$lambda_mu, main = latex2exp::TeX("Series: $(\\Lambda\\Lambda^T)_{12}$"))
  dev.off()
  pdf(file = paste0(dir2Save, '/lambda_mu/acf_', fname, '.pdf'), width = 4, height = 4)
  acf(df3$lambda_mu, main = latex2exp::TeX("Series: ${(\\Lambda\\mu)_1}$"))
  dev.off()
  pdf(file = paste0(dir2Save, '/noise/acf_', fname, '.pdf'), width = 4, height = 4)
  acf(df3$invsig1, main = latex2exp::TeX("Series: ${1}/{\\sigma_1^2}$"))
  dev.off()
  pdf(file = paste0(dir2Save ,'/tau/acf_', fname, '.pdf'), width = 4, height = 4)
  acf(df3$tau, main = latex2exp::TeX("Series: $\\tau$"))
  dev.off()
  return(q)
}
