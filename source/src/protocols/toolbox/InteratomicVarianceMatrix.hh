// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/Job.hh
/// @brief  header file for ThreadingJob classes, part of August 2008 job distributor as planned at RosettaCon08.  This file is responsible for three ideas: "inner" jobs, "outer" jobs (with which the job distributor works) and job container (currently just typdefed in the .fwd.hh)
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_toolbox_InteratomicVarianceMatrix_hh
#define INCLUDED_protocols_toolbox_InteratomicVarianceMatrix_hh

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <core/types.hh>

#include <ObjexxFCL/FArray3.fwd.hh>


namespace protocols {
namespace toolbox {

class InteratomicVarianceMatrix {
public:

	InteratomicVarianceMatrix( core::Size n_atoms, ObjexxFCL::FArray2D_double& ivm )
	: n_atoms_( n_atoms ), ivm_( ivm ) { };

	InteratomicVarianceMatrix( core::Size n_atoms = 0 )
	: n_atoms_( n_atoms ) { };

	core::Size n_atoms() const {
		return n_atoms_;
	}

	ObjexxFCL::FArray2D_double& ivm() {
		return ivm_;
	}

	//interatomic distance matrix --- no superposition needed for this!
	void init( core::Size n_atoms, core::Size n_decoys, ObjexxFCL::FArray3_double const& coords );


	void order_parameter( core::Real epsilon, ObjexxFCL::FArray1_double& T );

	core::Real kurtosis( ObjexxFCL::FArray1_double& T );

	void optimize_kurtosis( core::Size ngrid, core::Real lb, core::Real ub );

private:
	core::Size n_atoms_;
	ObjexxFCL::FArray2D_double ivm_;

	/// @brief order parameter to given value epsilon
	ObjexxFCL::FArray1D_double T_;

	/// @brief
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real epsilon_;


};

}
}

#endif
