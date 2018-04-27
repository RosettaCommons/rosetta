// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/TNA_SuitePotential.hh
/// @brief  TNA_SuitePotential potential class delcaration
/// @author Andy Watkins

#ifndef INCLUDED_core_scoring_rna_TNA_SuitePotential_HH
#define INCLUDED_core_scoring_rna_TNA_SuitePotential_HH

// Unit Headers
#include <core/scoring/rna/TNA_SuitePotential.fwd.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// ublas Headers
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace core {
namespace scoring {
namespace rna {

class TNA_SuitePotential : public utility::pointer::ReferenceCount {

public:

	/// @details This constructor reads in data from disk and should only 
	/// be called from the ScoringManager
	TNA_SuitePotential();

	virtual ~TNA_SuitePotential();

	bool eval_score(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose
	) const;

	Real get_score() const { return score_; }

	utility::fixedsizearray1<Real,4> get_deriv() const { return deriv_; }

	utility::vector1<id::TorsionID> get_torsion_ids() const {
		return torsion_ids_;
	}

private:

	void eval_likelihood_potential(
		utility::fixedsizearray1<Real, 4> const & torsions ) const;

	void regularize_torsions(
		boost::numeric::ublas::vector<Real> & torsions ) const;

	void figure_out_offset();

private:
	Size const n_torsions_ = 0;
	utility::vector1<Real> weights_;
	utility::vector1< boost::numeric::ublas::vector<Real> > centers_;
	utility::vector1<std::string> tags_;
	boost::numeric::ublas::matrix<Real> inv_cov_{ 4, 4 };
	Real offset_ = 0;
	mutable Real score_ = 0;
	mutable utility::fixedsizearray1<Real,4> deriv_ = 0;
	mutable utility::vector1<id::TorsionID> torsion_ids_;
};

} //rna
} //scoring
} //core

#endif
