library("inline")
  
eval_poly_sig <- signature(variable="numeric", params="array")
eval_poly_body <- '
  double const v(Rcpp::as<double>(variable));
  Rcpp::NumericVector p(params);


  double const & degree(p[0]);
  double const & xmin(p[1]);
  double const & xmax(p[2]);
  double value(0.0);

  if((v <= xmin) || (v >= xmax)){
    return(Rcpp::wrap(value));
  }
  value = p[3];
  for(unsigned int i=4; i < 3+degree; i++){
    (value *= v) += p[i];
  }
  return(Rcpp::wrap(value));
'
eval_poly <- cfunction(eval_poly_sig, eval_poly_body, Rcpp=TRUE)

eval_fade_sig <- signature(variable="numeric", params="array")
eval_fade_body <- '
  double const v(Rcpp::as<double>(variable));
  Rcpp::NumericVector p(params);
  bool const & s(p(0));
  double const & min0(p(1));
  double const & fmin(p(2));
  double const & fmax(p(3));
  double const & max0(p(4));

  double value;
//  if(s){
//    std::cout << "smooth" << std::endl;
//  }
//  if(!s){
//    std::cout << "not smooth" << std::endl;
//  }

  if(true){
    if (v <= fmax) {
      if (v <= min0) { value = 0.0;
      } else if (v >= fmin) { value = 1.0;
      } else { value = (v - min0) / static_cast<double>(fmin-min0);
      }
    } else {
      if (v >= max0){ value = 0.0;
      } else { value = (max0 - v) / static_cast<double>(max0-fmax);
      }
    }
  } else {
    if (v <= fmax) {
      if (v <= min0) { value = 0.0;
      } else if (v >= fmin) { value = 1.0;
      } else {
	double const z((v - min0)* static_cast<double>(fmin-min0));
        value = z*z*(3-2*z);
      }
    } else {
      if (v >= max0){ value = 0.0;
      } else {
        double const z((v - fmax) * static_cast<double>(max0-fmax));
        value = z*z*(2*z-3) + 1;
      }
    }
  }  
  return(Rcpp::wrap(value));
'
eval_fade <- cfunction(eval_fade_sig, eval_fade_body, Rcpp=TRUE)

eval_energy_sig <- signature(dofs="array", polys="matrix", fades="matrix")
eval_energy_body <- '

  Rcpp::NumericVector d(dofs);
  double const & AHdist(d(0));
  double const & cosBAH(d(1));
  double const & cosAHD(d(2));
  Rcpp::NumericMatrix p(polys);
  Rcpp::NumericMatrix f(fades);

  Rcpp::NumericVector AHdist_poly(p(0,Rcpp::_));
  Rcpp::NumericVector cosBAH_short_poly(p(1,Rcpp::_));
  Rcpp::NumericVector cosBAH_long_poly(p(2,Rcpp::_));
  Rcpp::NumericVector cosAHD_short_poly(p(3,Rcpp::_));
  Rcpp::NumericVector cosAHD_long_poly(p(4,Rcpp::_));

  Rcpp::NumericVector AHdist_short_fade(f(0,Rcpp::_));
  Rcpp::NumericVector AHdist_long_fade(f(1,Rcpp::_));
  Rcpp::NumericVector cosBAH_fade(f(2,Rcpp::_));
  Rcpp::NumericVector cosAHD_fade(f(3,Rcpp::_));

  double energy =
    Rcpp::as<double>(eval_poly(Rcpp::wrap(AHdist), AHdist_poly)) *
    Rcpp::as<double>(eval_fade(Rcpp::wrap(cosBAH), cosBAH_fade)) *
    Rcpp::as<double>(eval_fade(Rcpp::wrap(cosAHD), cosAHD_fade)) +
    Rcpp::as<double>(eval_fade(Rcpp::wrap(AHdist), AHdist_short_fade)) *
    Rcpp::as<double>(eval_poly(Rcpp::wrap(cosBAH), cosBAH_short_poly)) *
    Rcpp::as<double>(eval_fade(Rcpp::wrap(cosAHD), cosAHD_fade)) +
    Rcpp::as<double>(eval_fade(Rcpp::wrap(AHdist), AHdist_long_fade)) *
    Rcpp::as<double>(eval_poly(Rcpp::wrap(cosBAH), cosBAH_long_poly)) *
    Rcpp::as<double>(eval_fade(Rcpp::wrap(cosAHD), cosAHD_fade)) +
    Rcpp::as<double>(eval_fade(Rcpp::wrap(AHdist), AHdist_short_fade)) *
    Rcpp::as<double>(eval_fade(Rcpp::wrap(cosBAH), cosBAH_fade)) *
    Rcpp::as<double>(eval_poly(Rcpp::wrap(cosAHD), cosAHD_short_poly)) +
    Rcpp::as<double>(eval_fade(Rcpp::wrap(AHdist), AHdist_long_fade)) *
    Rcpp::as<double>(eval_fade(Rcpp::wrap(cosBAH), cosBAH_fade)) *
    Rcpp::as<double>(eval_poly(Rcpp::wrap(cosAHD), cosAHD_long_poly));

  if(energy >= 0){
    return(Rcpp::wrap(1));
  }
  double const density(exp(-1*energy)); 
  return(Rcpp::wrap(density));
'
eval_energy <- cfunction(
  list(
    eval_poly=eval_poly_sig,
    eval_fade=eval_fade_sig,
    eval_energy=eval_energy_sig),
  list(
    eval_poly_body,
    eval_fade_body,
    eval_energy_body),
  includes="#include<math.h>",
  Rcpp=TRUE)

eval_hbond_energy <- eval_energy[["eval_energy"]]
