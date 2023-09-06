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
#' @param dir A string, a directory where the resulting image is saved.
#' @param name A string, the name of the resulting image.
#' @param rotate_x_text Logical, whether to rotate the texts on the x-axis or not.
#' If so, they are positioned at a 60-degree angle.
#' @returns ggmatrix object from GGally, a matrix of plots
#' @export
#' @examples
#' if (require('GGally') && require('ggplot2')  &&
#' require('scales')) {
#'   pairsPlot_onedata(shell)
#'   }
pairsPlot_onedata <- function(dat_new,
                              dir = NULL,
                              name = NULL,
                              rotate_x_text = FALSE) {
  if (!requireNamespace("GGally", quietly = TRUE) ||
      !requireNamespace('ggplot2', quietly = TRUE) ||

      !requireNamespace('scales', quietly = TRUE) ) {
    stop(
      "Package \"GGally\", \"ggplot2\",\"scales\"
      must be installed to use this function.",
      call. = FALSE
    )
  }
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
               upper = list(continuous =
                              GGally::wrap(cor_func,method = 'spearman', symbol = '\u03C1 =')),
               #     upper = list(continuous = ggally_density, combo = ggally_box_no_facet),
               #   lower = list(continuous = wrap('points',size = 0.3)),
               lower = list(continuous = myscatter),
               diag = list(continuous = mydensity, axisLabels='show'))
  if(rotate_x_text) {
    q <- q +  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60))
  }
  #  axis.ticks = element_blank())
  if (!is.null(dir)) {
    if (!requireNamespace('grDevices', quietly = TRUE)) {
      stop(
        "Package \"grDevices\"
      must be installed to save the image.",
        call. = FALSE
      )
    }
    if(!dir.exists(dir)) {dir.create(dir)}
    ggplot2::ggsave(plot = q, filename = paste0(dir,'/', name,'.pdf'), width = 4.4,
           height = 4.4,device = grDevices::cairo_pdf)
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
#' @param dat_orig A matrix of data, e.g. the original data set.
#' @param dat_new A matrix of data, e.g. the data generated from the posterior
#' predictive distribution.
#' @param dir A string, a directory where the resulting image is saved.
#' @param name A string, the name of the resulting image.
#' @returns ggmatrix object from GGally, a matrix of plots
#' @export
#' @importFrom rlang .data
#' @examples
#' dat_new <- gen_posterior_predictive(1000, samples, 2500)
#' if (require('GGally') && require('ggplot2') && require('dplyr') &&
#' require('scales')) {
#'  pairsPlot(shell, dat_new)
#'   }
pairsPlot <- function(dat_orig, dat_new, dir = NULL,name = NULL) {
  if (!requireNamespace("GGally", quietly = TRUE) ||
      !requireNamespace('ggplot2', quietly = TRUE) ||
    #  !requireNamespace('dplyr', quietly = TRUE) ||
      !requireNamespace('scales', quietly = TRUE) ) {
    stop(
      "Package \"GGally\", \"ggplot2\", and \"scales\"
      must be installed to use this function.",
      call. = FALSE
    )
  }
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
  q <- GGally::ggpairs(df, columns = 1:num,ggplot2::aes(group = .data$type,
                                                        col = .data$type,
                                                        fill = .data$type,alpha = 0.5),
               #  upper = list(continuous = wrap(cor_func,
               #     method = 'spearman')),
               upper = list(continuous = suppressWarnings(mycorrelations)),
         #     upper = list(continuous = wrap("cor", title = '\u03C1')),
       #    upper = list(continuous = wrap(cor_func,
                         #         method = 'spearman', symbol = '\u03C1 =')),
               #     upper = list(continuous = ggally_density, combo = ggally_box_no_facet),
                #  lower = list(continuous = wrap('points',size = 0.3,alpha = 0.5)),
              lower = list(continuous = myscatter),
               diag = list(continuous = mydensity, axisLabels='show'))
  if (!is.null(dir)) {
    if (!requireNamespace('grDevices', quietly = TRUE)) {
      stop(
        "Package \"grDevices\"
      must be installed to save the image.",
        call. = FALSE
      )
    }
    if(!dir.exists(dir)) {dir.create(dir)}
    ggplot2::ggsave(plot = q, filename = paste0(dir,'/', name,'.pdf'), width = 4.4,
           height = 4.4,device = grDevices::cairo_pdf)
  }
  return(q)
}
#

#' Correlation for the scatter plot matrix
#'
#'@description
#' `cor_func` calculates correlation used for the scatter plot matrix
#'  constructed by GGally::ggpairs().
#'
#'@param data The data set to be plotted.
#'@param mapping The mapping.
#'@param method The method of correlation calculation.
#'@param symbol Symbol.
#'@param ... ...
#'
#'@importFrom stats cor
# Defines function to color according to correlation
cor_func <- function(data, mapping, method, symbol, ...){
  if (!requireNamespace("GGally", quietly = TRUE) ||
      !requireNamespace('ggplot2', quietly = TRUE)) {
    stop(
      "Package \"GGally\" and \"ggplot2\"
      must be installed to use this function.",
      call. = FALSE
    )
  }
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  corr <- stats::cor(x, y, method=method, use='complete.obs')
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
#' Correlation for the GGally::ggpairs
#'
#' @description
#' `mycorrelations` supplies the correlation used for the GGally::ggpairs.
#'
#' @param data A matrix of data.
#' @param mapping Mapping.
#' @param ... ...
#'
#' @importFrom dplyr bind_rows mutate group_by filter summarize %>% as_label
#' @importFrom stats cor.test
#' @importFrom rlang .data
mycorrelations <- function(data,mapping,...){
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop(
      "Package \"ggplot2\"
      must be installed to use this function.",
      call. = FALSE
    )
  }
  data2 <- data
  data2$x <- as.numeric(data[,dplyr::as_label(mapping$x)])
  data2$y <- as.numeric(data[,dplyr::as_label(mapping$y)])
  data2$group <- data[,dplyr::as_label(mapping$colour)]
  correlation_df <- suppressWarnings(data2 %>%
    dplyr::bind_rows(data2 %>% dplyr::mutate(group="Overall Corr")) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::filter(sum(!is.na(.data$x),na.rm=TRUE)>1) %>%
    dplyr::filter(sum(!is.na(.data$y),na.rm=TRUE)>1) %>%
    dplyr::summarize(estimate = round(as.numeric(stats::cor.test(.data$x,.data$y,
                                                          method="spearman")$estimate),2)) %>%
       #       pvalue = cor.test(x,y,method="spearman")$p.value,
         #     pvalue_star = as.character(symnum(pvalue, corr = FALSE, na = FALSE,
                                          #      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                            #    symbols = c("***", "**", "*", "'", " "))))%>%
    dplyr::group_by() %>%
    dplyr::mutate(group = factor(.data$group, levels=c(as.character(unique(sort(data[,dplyr::as_label(mapping$colour)]))), "Overall Corr")))
)
  correlation_df <- correlation_df %>% dplyr::filter(!.data$group %in% 'Overall Corr')
  ggplot2::ggplot(data=correlation_df, ggplot2::aes(x=1,y=.data$group,color=.data$group))+
    ggplot2::geom_text(ggplot2::aes(label=paste0(.data$group,": ",.data$estimate)), size = 6) + ggplot2::theme_classic()

}

#' Scatter plots for the GGally::ggpairs plot
#'
#' @description
#' `myscatter` provides the scatter plots used in the GGally::ggpairs plot.
#'
#'@param data A matrix of data.
#'@param mapping Mapping.
#'@param ... ...
#'
myscatter <- function(data, mapping, ...) {
  if (!requireNamespace('ggplot2', quietly = TRUE) ||
      !requireNamespace('scales', quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" and \"scales\"
      must be installed to use this function.",
      call. = FALSE
    )
  }
  ggplot2::ggplot(data = data, mapping = mapping, ...) +
    ggplot2::geom_point(size = 0.3, alpha = 0.3,...) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +ggplot2::theme_minimal()
}

#' Density plots for the GGally::ggpairs plot
#'
#' @description
#' `mydensity` provides the density plots used in the GGally::ggpairs plot.
#'
#' @param data A data frame of data.
#' @param mapping Mapping.
#' @param ... ...
#'
mydensity <- function(data, mapping, ...) {
  if (!requireNamespace('ggplot2', quietly = TRUE) ||
      !requireNamespace('scales', quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" and \"scales\"
      must be installed to use this function.",
      call. = FALSE
    )
  }
  ggplot2::ggplot(data = data, mapping = mapping, ...) +
    ggplot2::geom_density(alpha= 0.5,...) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=3)) + ggplot2::theme_minimal()

}

#' Visualize the factor loadings
#'
#' @description
#' `plot_factor_loadings` visualizes the posterior mean of the factor loadings matrix. We
#' recommend post processing the samples first using [draw_latent_factors()] and
#' [postprocess()]. The function uses (post-processed) samples of factor loadings
#' directly.
#'
#' @param lambda_samples A list of the samples, length = number of samples,
#' each element is p by k.
#' @param variable_names A vector of strings, length = p, variable names. Default
#' is NULL, meaning unspecified.
#' @param row_idx A vector row indices of the factor loadings to be plotted. Default
#' is all the rows are included.
#' @param saveDir A string, the path to where the image should be saved. Default
#' is NULL, meaning the image is not saved. When the image is saved, its name is
#' loadings.pdf.
#' @export
#'
#' @examples
#' lambda_samples <- list(a = matrix(rnorm(6),3,2), b = matrix(rnorm(6),3,2))
#' plot_factor_loadings(lambda_samples = lambda_samples)
plot_factor_loadings <- function(lambda_samples,
                                 variable_names = NULL,
                                 row_idx = NULL,
                                 saveDir = NULL) {
  if (!requireNamespace("infinitefactor", quietly = TRUE)) {
    stop(
      "Package \"infinitefactor\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!is.null(variable_names)) {
    stopifnot(nrow(lambda_samples[[1]]) == length(variable_names))
  }

  if (is.null(row_idx)) {
    row_idx <- seq_len(nrow(lambda_samples[[1]]))
  }

  mat <- infinitefactor::lmean(lambda_samples)[row_idx,]
  rownames(mat) <- variable_names[row_idx]
  p1 <- plotmat(mat = mat, saveDir = saveDir)
  return(p1)
}


#' Visualize the factor loadings
#'
#' @description
#' `plotmat` builds on the function [infinitefactor::plotmat()] and includes
#' additional features, including an option to save the figure and automatic
#' inclusion of the variable names in the plot.
#'
#' @param mat A matrix of numeric values.
#' @param color The color scheme of the plot, "green", "red" or "wes".
#' @param title A string, optional plot title.
#' @param args Optional additional ggplots arguments
#' @param saveDir The path to where the image is saved. Default is NULL, meaning
#' the figure is not saved. When it is saved, it is named as loadings.pdf.
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(6), 3, 2)
#' plotmat(mat = mat)
#'
plotmat <- function (mat,
                     color = "green",
                     title = NULL,
                     args = NULL,
                     saveDir = NULL){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }
  mat = apply(mat, 2, rev)
  longmat = melt(mat)
  Var1 = Var2 = value = NULL
  p = ggplot2::ggplot(longmat, ggplot2::aes(x = Var2, y = Var1)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value),colour = "grey20")
  if (color == "green")
    p = p + ggplot2::scale_fill_gradient2(low = "#3d52bf", high = "#33961b",
                                 mid = "white")
  if (color == "red")
    p = p + ggplot2::scale_fill_gradient2(low = "#191970", high = "#800000",
                                 mid = "white")
  if (color == "wes")
    p = p + ggplot2::scale_fill_gradient2(low = "#046C9A", high = "#D69C4E",
                                 mid = "white")
  p = p + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank(),
                         panel.grid.major = ggplot2::element_blank(),
                         panel.border = ggplot2::element_blank(),
                         panel.background = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(),
                         legend.title = ggplot2::element_text(),
                         plot.title = ggplot2::element_text(hjust = 0.5),
                         axis.text.y = ggplot2::element_text(angle = 60)) +
    ggplot2::labs(fill = " ")
  if (!is.null(title))
    p = p + ggplot2::ggtitle(title)
  if (!is.null(args))
    p = p + args
  if (!is.null(saveDir)) {
    ggplot2::ggsave(filename = paste0(saveDir,'/','loadings.pdf'),
                    plot = p,
                    width = 3,
                    height =3)
  }
  return(p)
}
