// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/OmegaTether.hh
/// @brief  OmegaTether potential class delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_OmegaTether_hh
#define INCLUDED_core_scoring_OmegaTether_hh

// Unit Headers
#include <core/scoring/OmegaTether.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ProteinTorsion.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray4D.hh>

//Auto Headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.fwd.hh>
#include <iostream>


namespace core {
namespace scoring {



class OmegaTether : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	OmegaTether();
	~OmegaTether() {}

	Real
	eval_omega_score_residue(
		AA const res_aa,
		Real const omega
	) const;

	void
	eval_omega_score_residue(
		conformation::Residue const & res,
		Real & energy,
		Real & denergy_domega
	) const;

	void
	eval_omega_score_residue(
		AA const res_aa,
		Real const omega,
		Real & energy,
		Real & denergy_domega
	) const;


	void
	eval_omega_score_all(
		Pose & pose,
		ScoreFunction const & scorefxn
	) const;


/*
	void
	write_omega_score_all(
		Pose const & pose
	) const;

	//Real get_omega_score_residue_deriv( int res, Pose const & a_pose, ProteinTorsion torsion ) const;
	void eval_procheck_omega( Pose const & a_pose,
		Real & favorable, Real & allowed, Real & generous ) const;
*/

private:


};

}
}

#endif
