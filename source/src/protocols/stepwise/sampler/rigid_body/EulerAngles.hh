// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rigid_body/EulerAngles.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_EulerAngles_HH
#define INCLUDED_protocols_swa_rna_EulerAngles_HH

#include <core/types.hh>
#include <protocols/stepwise/sampler/rigid_body/EulerAngles.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rigid_body {

class EulerAngles: public utility::pointer::ReferenceCount {

public:

	//Should make sure that alpha and gamma lies in the [-Pi:Pi] range.
	//constructor
	EulerAngles();

	EulerAngles( numeric::xyzMatrix< core::Real > const & rotation_matrix );

	//destructor
	~EulerAngles();

public:

	void
	initialize_from_rotation_matrix( numeric::xyzMatrix< core::Real > const & rotation_matrix );

	void
	convert_to_rotation_matrix( numeric::xyzMatrix< core::Real > & rotation_matrix );

	void
	set_alpha( Real const setting ){ alpha_ = setting; }

	void
	set_beta( Real const setting );

	void
	set_z( Real const setting );

	void
	set_gamma( Real const setting ){ gamma_ = setting; }

	Real const & alpha() const { return alpha_; }
	Real const & beta() const { return beta_; }
	Real const & gamma() const { return gamma_; }
	Real const & z() const { return z_; }

private:

	core::Real alpha_; //phi
	core::Real beta_; //theta
	core::Real gamma_; //psi
	core::Real z_; //z=cos(beta)

};

} //rigid_body
} //sampler
} //stepwise
} //protocols

#endif
