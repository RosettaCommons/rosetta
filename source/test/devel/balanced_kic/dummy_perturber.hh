// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/balanced_kic/assertion_helpers.hh
/// @brief  Provide a dummy kinematic perturber implementation.
/// @author Kale Kundert
///
/// The dummy perturber provided by this header is meant to test code that 
/// invokes a perturber, but that doesn't care exactly what that perturber 
/// does.  Testing real perturbers directly is difficult because they are often 
/// stochastic.  Testing this perturber is trivial, because it simply changes 
/// the coordinates of the first atom to (-1, -1, -1).

#ifndef INCLUDED_devel_balanced_kic_dummy_perturer_HH
#define INCLUDED_devel_balanced_kic_dummy_perturer_HH

#include <devel/balanced_kic/algorithms.hh>
#include <devel/balanced_kic/KinematicPerturber.hh>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>

using namespace core;
using namespace devel::balanced_kic;

class DummyPerturber : public KinematicPerturber {

public:

	DummyPerturber() {
		exhausted_ = false;
	}

	std::string perturber_type() const {
		return "TestKinematicPerturber";
	}

	void perturb_chain(
			pose::Pose const &pose,
			algorithms::ParameterList &torsion_angles,
			algorithms::ParameterList &bond_angles,
			algorithms::ParameterList &bond_lengths) const {

		torsion_angles[1] = -1;
		bond_angles[1] = -1;
		bond_lengths[1] = -1;

		exhausted_ = true;
	}

	bool perturber_exhausted() const {
		return exhausted_;
	}

private:
	mutable bool exhausted_;

};

#endif
