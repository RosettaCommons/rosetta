// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/RNA_SuitePotential.hh
/// @brief  RNA_SuitePotential potential class delcaration
/// @author Fang-Chieh Chou

#ifndef INCLUDED_core_scoring_rna_RNA_SuitePotential_HH
#define INCLUDED_core_scoring_rna_RNA_SuitePotential_HH

// Unit Headers
#include <core/scoring/rna/RNA_SuitePotential.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// ublas Headers
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace core {
namespace scoring {
namespace rna {

class RNA_SuitePotential : public utility::pointer::ReferenceCount {

public:
	RNA_SuitePotential();

	virtual ~RNA_SuitePotential();

	bool eval_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose
	) const;

	Real get_score() const { return score_; }

	utility::vector1<Real> get_deriv() const { return deriv_; }

	utility::vector1<id::TorsionID> get_torsion_ids() const {
		return torsion_ids_;
	}

private:
	Size const n_torsions_;
	utility::vector1<Real> weights_;
	utility::vector1< boost::numeric::ublas::vector<Real> > centers_;
	boost::numeric::ublas::matrix<Real> inv_cov_;
	Real offset_;
	mutable Real score_;
	mutable utility::vector1<Real> deriv_;
	mutable utility::vector1<id::TorsionID> torsion_ids_;

	void eval_score(
			boost::numeric::ublas::vector<Real> const & torsions ) const;

	void regularize_torsions(
			boost::numeric::ublas::vector<Real> & torsions ) const;
};

} //rna
} //scoring
} //core

#endif
