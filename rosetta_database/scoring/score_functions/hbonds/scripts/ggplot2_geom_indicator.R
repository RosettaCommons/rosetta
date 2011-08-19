require(ggplot2)

GeomIndicator <- proto(Geom, {
  objname <- "indicator"
  draw <- function(., data, scales, coordinates, ...){
    indicator <- prettyNum(data$indicator[1], big.mark=",")    
    if(!is.na(indicator) && !is.null(indicator)){
    	level <- data$group[1] - 1
        colour <- data$colour[1]
    	textGrob(indicator,
             unit(.97, "npc"), unit(.97, "npc") - unit(level, "line"),
             just=c("right", "top"), hjust=1, vjust=1,
             gp=gpar(col=colour, cex=.75))
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
