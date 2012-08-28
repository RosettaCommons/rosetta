// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran.hh
/// @brief  Ramachandran potential class delcaration
/// @author Guoli Wang

#ifndef INCLUDED_core_scoring_Ramachandran2B_hh
#define INCLUDED_core_scoring_Ramachandran2B_hh

// Unit Headers
#include <core/scoring/Ramachandran2B.fwd.hh>
//#include <core/scoring/ProteinTorsion.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {


class Ramachandran2B : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	Ramachandran2B();
	virtual ~Ramachandran2B() ; // auto-removing definition from header{}

	Real
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi
	) const;

	// Guoli Wang
	void
	eval_rama_score_residue(
		conformation::Residue const & res,
		chemical::AA const left_aa,
		chemical::AA const right_aa,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	Real
	RamaE_Lower(
		conformation::Residue const &rsd,
		chemical::AA const &neighbor,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE_Lower(
		conformation::Residue const &rsd,
		chemical::AA const &neighbor
	) const;

	Real
	RamaE_Upper(
		conformation::Residue const & rsd,
		chemical::AA const &neighbor,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE_Upper(
		conformation::Residue const & rsd,
		chemical::AA const &neighbor
	) const;

	Real
	RamaE(
		conformation::Residue const & rsd,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE(
		conformation::Residue const & rsd
	) const;
	// finished

	void
	eval_rama_score_residue(
		conformation::Residue const & res,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	IdealizeRamaEnergy(
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi,
		Real const entropy,
		ObjexxFCL::
		FArray2A< Real > const & rama_for_res
	) const;

	///////////////////////////////
	// unused??
	void
	eval_rama_score_all(
		Pose & pose,
		ScoreFunction const & scorefxn
	) const;

	void
	write_rama_score_all(
		Pose const & pose
	) const;


private:

	void read_rama();

	//static bool rama_initialized_;
	ObjexxFCL::FArray3D< Real > ram_energ_;
	ObjexxFCL::FArray1D< Real > ram_entropy_;
	ObjexxFCL::FArray4D< Real > ram_energ_left_;
	ObjexxFCL::FArray2D< Real > ram_entropy_left_;
	ObjexxFCL::FArray4D< Real > ram_energ_right_;
	ObjexxFCL::FArray2D< Real > ram_entropy_right_;

	static int const n_phi_ = 36;
	static int const n_psi_ = 36;
	static Real const binw_; // 360 / n_phi_ = 10;
	static int const n_aa_ = 20; // Ramachandran score defined for the cananical AAs only.
	static int const nullaa = 21; // Guoli Wang

	Real const rama_score_limit_;
};

}
}

#endif
