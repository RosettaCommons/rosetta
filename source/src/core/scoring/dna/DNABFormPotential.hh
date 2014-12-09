// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dna/DNABFormPotential.hh
/// @brief  DNA B-form specific torsion potential class delcaration
/// @author Jim Havranek

#ifndef INCLUDED_core_scoring_dna_DNABFormPotential_HH
#define INCLUDED_core_scoring_dna_DNABFormPotential_HH

// Unit Headers
#include <core/scoring/dna/DNABFormPotential.fwd.hh>
#include <core/scoring/ProteinTorsion.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray4D.hh>

namespace core {
namespace scoring {
namespace dna {

class TorsionFourierComponent : public utility::pointer::ReferenceCount
{

public:
	TorsionFourierComponent( Real factor, Real periodicity, Real phase ) :
			factor_( factor ), periodicity_( periodicity ), phase_( phase ) {}

	Real compute( Real const torsion_angle, Real & deriv ) const;

	Real factor() const { return factor_; }
	Real periodicity() const { return periodicity_; }
	Real phase() const { return phase_; }

private:
	Real factor_;				// The value for Vn/2 - analog to spring constant
	Real periodicity_;	// How many minima for torsion - i.e. two SP3 atoms would have 3
	Real phase_;				// Offset to minimum in potential

};


class DNABFormPotential : public utility::pointer::ReferenceCount
{

public:
	DNABFormPotential();
	~DNABFormPotential() {}

	void
	eval_dna_bform_bb_torsion_score_residue(
		conformation::Residue const & res,
		Real & score,
		Real & dscore_dchi,
		Size const torsion_id
	) const;

	void
	eval_dna_bform_chi_torsion_score_residue(
		conformation::Residue const & res,
		Real & score,
		Real & dscore_dchi
	) const;

private:

	void init_dna_bform_data();
	int dummy_;
	utility::vector1< utility::vector1< TorsionFourierComponentCOP > > bb_fourier_data;

};

}
}
}

#endif
