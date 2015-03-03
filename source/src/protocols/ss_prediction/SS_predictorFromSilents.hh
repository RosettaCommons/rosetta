// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ss_predictorFromSilents.hh
/// @brief
/// @detailed
///
/// @author TJ Brunette


#ifndef INCLUDED_protocols_ss_prediction_SS_predictorFromSilents_hh
#define INCLUDED_protocols_ss_prediction_SS_predictorFromSilents_hh

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
using namespace utility::libsvm;

class SS_predictorFromSilents : public utility::pointer::ReferenceCount {
 public:
	/// @brief Reads in models for SS prediction etc.
	SS_predictorFromSilents(string type);
	~SS_predictorFromSilents();
	vector1<vector1 <Real> > predict_ss();

 private:
}; 
// helper functions

/// @brief helper function to get SS char at a position from a vector of reals
char get_label( utility::vector1< core::Real > const & ss_pred_pos );

/// @brief helper function to get SS char at a position from a vector of reals
core::Real get_prob( char wanted_ss, utility::vector1< core::Real > const & ss_pred_pos );

} //ss_prediction
} // protocols

#endif
