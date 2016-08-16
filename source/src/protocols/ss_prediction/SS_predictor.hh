// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ss_predictor.hh
/// @brief
/// @details
///
/// @author TJ Brunette


#ifndef INCLUDED_protocols_ss_prediction_SS_predictor_hh
#define INCLUDED_protocols_ss_prediction_SS_predictor_hh

//// utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <utility/libsvm/Svm_rosetta.hh>
#include <utility/libsvm/Svm_rosetta.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace ss_prediction {

using std::string;
using core::Real;
using core::Size;
using core::SSize;
using utility::vector1;
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace utility::libsvm;

class SS_predictor : public utility::pointer::ReferenceCount {
public:
	/// @brief Reads in models for SS prediction etc.
	SS_predictor(string type);
	~SS_predictor();
	vector1<vector1 <Real> > predict_ss(string fasta);

private:
	Svm_rosettaOP rd1_model;
	Svm_rosettaOP rd2_model;
	string ss_type;
	void load_models(string rd1_model_fl, string rd2_model_fl);
	string get_window_aa(string fasta,SSize position);
	vector1 <Real> predict_pos_rd1(string window_aa);
	vector1 <Real> predict_pos_rd2(vector1< vector1 <Real> > rd1_preds,SSize position,string fasta);
	vector1< vector1<Real> > predict_rd1(string fasta);
	vector1< vector1<Real> > predict_rd2(vector1<vector1<Real> > rd1_preds, string fasta);
	static const Size WINDOW_SIZE = 15;
}; // class SS_predictor

// helper functions

/// @brief helper function to get SS char at a position from a vector of reals
char get_label( utility::vector1< core::Real > const & ss_pred_pos );

/// @brief helper function to get SS char at a position from a vector of reals
core::Real get_prob( char wanted_ss, utility::vector1< core::Real > const & ss_pred_pos );

} //ss_prediction
} // protocols

#endif
