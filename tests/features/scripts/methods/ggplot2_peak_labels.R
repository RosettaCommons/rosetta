require(ggplot2)

GeomPeakLabels <- proto(Geom, {
  objname <- "peak_labels"
  draw <- function(., data, scales, coordinates, ...){

  
    peakfun <- function(x) {
      d <- density(x$length)
      peaks <- which(diff(sign(diff(d$y)))==-2)
      data.frame(x=d$x[peaks],y=d$y[peaks])
    }
    
  }
  desc_params <- list()
  default_stat <- function(.) StatIdentity
  icon <- function(.) textGrob("text", gp=gpar(cex=1.2))
  desc <- "Count Instances"
  guide_geom <- function(x) "text"
  example <- function(.){
    data <- rbind(
      data.frame(x= c(rnorm(200)+5,rnorm(170)+2),strata=factor(1),slice =factor(1)),
      data.frame(x= c(rnorm(500)+5,rnorm(32 )+2),strata=factor(2),slice =factor(1)),
      data.frame(x= c(rnorm(356)+5,rnorm(120)+2),strata=factor(2),slice =factor(2)))
    p <- ggplot(data=data, aes(x=x, colour=strata))
    p <- p + geom_density()
    p <- p + geom_indicator()
    p <- p + facet_wrap( ~ slice )
    print(p)
  }
})
geom_indicator <- GeomIndicator$build_accessor()