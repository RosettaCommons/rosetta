// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  A class to wrap libsvm for rosetta
/// @author TJ Brunette


// Utility Headers
#include <utility/libsvm/Svm.hh>
#include <utility/libsvm/Svm_rosetta.hh>
#include <platform/types.hh>

#include <string>
#include <iostream>

namespace utility {
namespace libsvm {

using utility::vector1;
using std::string;

Svm_node_rosetta::Svm_node_rosetta(platform::Size index, platform::Real value){
	index_ = index;
	value_ = value;
}

Svm_node_rosetta::~Svm_node_rosetta()= default;


Svm_rosetta::Svm_rosetta(string model_filename){
	const char* model_filename_c= model_filename.c_str();
	svm_model_ = svm_load_model(model_filename_c);
}

Svm_rosetta::~Svm_rosetta(){
	svm_free_model_content( svm_model_ );
	free( svm_model_ );
}

platform::Size Svm_rosetta::get_nr_class(){
	return((platform::Size)svm_get_nr_class(svm_model_));
}

vector1 < platform::Real > Svm_rosetta::predict_probability(vector1 <Svm_node_rosettaOP> & features){
	// TL 5/2013: Changed to use new and delete[] to avoid problems
	//struct svm_node *x = (struct svm_node *) malloc((features.size()+1)*sizeof(struct svm_node));
	auto *x = new svm_node[features.size()+1];
	for ( platform::Size ii=1; ii<=features.size(); ++ii ) {
		x[ii-1].index = (int)features[ii]->index();
		x[ii-1].value = (double)features[ii]->value();
	}
	x[features.size()].index = -1;
	int nr_class = svm_get_nr_class(svm_model_);
	//double *prob_estimates = (double *) malloc(nr_class*sizeof(double));
	auto *prob_estimates = new double[nr_class];
	/*double predict_label =*/ svm_predict_probability(svm_model_,x,prob_estimates);
	vector1 <platform::Real> probs_to_return;
	for ( int ii=0; ii<nr_class; ++ii ) {
		probs_to_return.push_back(prob_estimates[ii]);
	}
	delete[] prob_estimates;
	delete[] x;
	return(probs_to_return);
}

platform::Real Svm_rosetta::predict( const vector1 <Svm_node_rosettaOP> & features)
{
	// TL 5/2013: Changed to use new and delete[] to avoid problems
	//struct svm_node *x = (struct svm_node *) malloc((features.size()+1)*sizeof(struct svm_node));
	auto *x = new svm_node[features.size()+1];
	for ( platform::Size ii=1; ii<=features.size(); ++ii ) {
		x[ii-1].index = (int)features[ii]->index();
		x[ii-1].value = (double)features[ii]->value();
	}
	x[features.size()].index = -1;
	platform::Real predict_value( svm_predict(svm_model_,x) );
	delete[] x;
	return(predict_value);
}

}//libsvm
}//utility
