// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MembranePotential.hh
/// @brief  ProQembrane Potential
/// @author Bjorn Wallner


#ifndef INCLUDED_core_scoring_ProQPotential_hh
#define INCLUDED_core_scoring_ProQPotential_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/ProQPotential.fwd.hh>
//#include <core/scoring/EnvPairPotential.hh>

// Package headers
//#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

//#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
//#include <ObjexxFCL/FArray3D.hh>
//#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

class ProQPotential  {

public:
	ProQPotential();

	inline
	~ProQPotential()
	{}

	void
	setup_for_scoring( pose::Pose & pose ) const;

	void
	score(pose::Pose & pose,
		ObjexxFCL::FArray2D< Real > & feature_vector,
		ObjexxFCL::FArray1D< Real > & score,
		bool ProQ2=false
	) const;

	Size
	num_features_proqm() const;

	Size
	num_features_proq2() const;

private: // data

	ObjexxFCL::FArray2D< Real > linear_weights_; //ProQM
	ObjexxFCL::FArray1D< Real > b_; // ProQM
	ObjexxFCL::FArray2D< Real > linear_weights_proq2_;
	ObjexxFCL::FArray1D< Real > b_proq2_;

	//static Size const num_features_=260;
	static Size const num_features_proqm_=260;
	static Size const num_features_proq2_=174;
	static Size const num_models_=5;
	Size svm_model_;  //if prediction only should use one specific svm model, typically in cross validation.
	bool cross_val_;


};


} // ns scoring
} // ns core

#endif
