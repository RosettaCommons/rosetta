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

class SS_predictor : public utility::pointer::ReferenceCount {
public:
	/// @brief Reads in models for SS prediction etc.
	SS_predictor(std::string type);
	~SS_predictor() override;
	utility::vector1< utility::vector1 <core::Real> > predict_ss(std::string fasta);

private:
	utility::libsvm::Svm_rosettaOP rd1_model;
	utility::libsvm::Svm_rosettaOP rd2_model;
	std::string ss_type;
	void load_models(std::string rd1_model_fl, std::string rd2_model_fl);
	std::string get_window_aa(std::string fasta,core::SSize position);
	utility::vector1 <core::Real> predict_pos_rd1(std::string window_aa);
	utility::vector1 <core::Real> predict_pos_rd2(utility::vector1< utility::vector1 <core::Real> > rd1_preds, core::SSize position, std::string fasta);
	utility::vector1< utility::vector1<core::Real> > predict_rd1(std::string fasta);
	utility::vector1< utility::vector1<core::Real> > predict_rd2(utility::vector1< utility::vector1<core::Real> > rd1_preds, std::string fasta);
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
