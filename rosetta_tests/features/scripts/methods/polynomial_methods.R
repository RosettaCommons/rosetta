
require(polynom)

get_1d_polynomial_by_name <- function(polynomials,
                                      name){
  df <- polynomials[polynomials$name==name,]
  if (nrow(df) !=1){
    print( paste("Polynomial", name," could not be found"))
    stop(1)
  }
  #we don't know the order in which the rows columns selected will be returned!  
  switch(df$degree,
         p <- df[1,c("c_a")],
         p <- df[1,c("c_b","c_a")],
         p <- df[1,c("c_c","c_b","c_a")],
         p <- df[1,c("c_d","c_c","c_b","c_a",)],
         p <- df[1,c("c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_k","c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")])
  model <- list(poly=polynomial(p), name=df$name)
#  print(paste("Polynomial ",model$name,": ", paste(p,sep=","), sep=""))
  model
}

  

get_1d_polynomial <- function(polynomials,
                              don_chem_type=NULL,
                              acc_chem_type=NULL,
                              seq_sep=NULL
                              ){
  if(is.null(don_chem_type)) don_chem_type <- "hbdon_NONE"
  if(is.null(acc_chem_type)) acc_chem_type <- "hbacc_NONE"
  if(is.null(seq_sep)) seq_sep <- "seq_sep_other"
  
  df <- polynomials[
	as.character(polynomials$don_chem_type)==as.character(don_chem_type) &
        as.character(polynomials$acc_chem_type)==as.character(acc_chem_type) &
	as.character(polynomials$seq_sep)==as.character(seq_sep),]

  if (nrow(df) !=1){
    print(paste("Polynomial for (",
                don_chem_type,", ",
                acc_chem_type,", ",
                seq_sep,") is not found",
                sep=""))
    stop(0)
  }

  #we don't know the order in which the rows columns selected will be returned!  
  switch(df$degree,
         p <- df[1,c("c_a")],
         p <- df[1,c("c_b","c_a")],
         p <- df[1,c("c_c","c_b","c_a")],
         p <- df[1,c("c_d","c_c","c_b","c_a",)],
         p <- df[1,c("c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")],
         p <- df[1,c("c_k","c_j","c_i","c_h","c_g","c_f","c_e","c_d","c_c","c_b","c_a")])
  model <- list(poly=polynomial(p), xmin=df$xmin, xmax=df$xmax, name=df$name)
  #print(paste("Polynomial ",model$name,": ", paste(p,sep=","), sep=""))
  model
}

add_1d_polynomial_models <- function(polynomials,
                                     densities,
                                     variables,
                                     model_transform=I,
                                     don_chem_type=NULL,
                                     acc_chem_type=NULL,
                                     seq_sep=NULL){

  densp <- ddply(densities,
               .variables=variables,
               function(df){
                 if( any(variables == "acc_chem_type")){
                   acc_chem_type=df$acc_chem_type[1]
                 }
                 if( any(variables == "don_chem_type")){
                   don_chem_type=df$don_chem_type[1]
                 }
                 if( any(variables == "seq_sep")){
                   don_chem_type=df$seq_sep[1]
                 }

                 model <- get_1d_polynomial(polynomials,
                                            don_chem_type=don_chem_type,
                                            acc_chem_type=acc_chem_type,
                                            seq_sep=seq_sep)
                 df_fit <- df
                 df_fit$y <-model_transform(predict(model$poly, df$x))
                 df_fit$sample_source <- model$name
                 rbind(df, df_fit)
               })
  return(densp)              
}
