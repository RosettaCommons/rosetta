#creating a new scals

require(ggplot2)


#TransLogPosNeg <- Trans$new(
#  name="LogPosNeg",
#  f=function(x) sign(x)*log(abs(x)+1),
#  inverse=function(x) sign(x)*(exp(abs(x)-1)),
#  labels=function(x) sign(x)*(exp(abs(x))-1))
#
#ScaleLogPosNeg <- proto(ScaleContinuous,
#  desc = "Position scale, log transformed with negative values log transformed symmetric across.",
#  tr_default = Trans$find("LogPosNeg"),
#  objname = "log_pos_neg",
#  doc=FALSE,
#  examples=function(.) {}
#)
#
#scale_x_log_pos_neg <- ScaleLogPosNeg$build_accessor(list(variable = "\"x\""))
#scale_y_log_pos_neg <- ScaleLogPosNeg$build_accessor(list(variable = "\"y\""))
