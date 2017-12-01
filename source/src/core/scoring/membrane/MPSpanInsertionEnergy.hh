// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane/SpanInsertionEnergy.hh
/// @brief  penalize spans with unnatural dG of insertion
/// @author Jonathan Weinstein


#ifndef INCLUDED_core_scoring_membrane_MPSpanInsertionEnergy_hh
#define INCLUDED_core_scoring_membrane_MPSpanInsertionEnergy_hh

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>
#include <core/conformation/membrane/Span.hh>
#include <map>

//Objexx headers

// Utility headers


namespace core {
namespace scoring {
namespace membrane {


class MPSpanInsertionEnergy : public methods::WholeStructureEnergy  {
public:
	typedef methods::WholeStructureEnergy  parent;

public:

	MPSpanInsertionEnergy();

	//clone
	methods::EnergyMethodOP
	clone() const override;


	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const override;

	void
	setup_for_derivatives(
		pose::Pose &,
		ScoreFunction const &
	)
	const override;

	void
	eval_atom_derivative(
		id::AtomID const &,
		pose::Pose const &,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const &,
		Vector &,// F1,
		Vector & // F2
	) const override;

	core::Real
	compute(
		pose::Pose const & pose
	) const;

	utility::vector1< core::conformation::membrane::Span >
	create_updated_span(
		pose::Pose const & pose
	) const;

	core::Real
	calc_span_score(
		pose::Pose const & pose, core::Size start, core::Size end
	) const;


	void
	parse_elazar_spline_table();

	numeric::interpolation::spline::CubicSpline
	split_line_to_spline(
		std::string line
	) const;

	core::Real
	spline_by_z(
		char const & res,
		core::Real const & z
	) const;

	core::Size version() const override;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:
	std::map< char, numeric::interpolation::spline::CubicSpline > res_spline_map_;
	core::Real avg_dg_span_single_;
	core::Real std_dg_span_single_;

	core::Real avg_dg_span_multi_;
	core::Real std_dg_span_multi_;

};

} //membrane
} //scoring
} //core

#endif // INCLUDED_core_scoring_membrane_MPSpanInsertionEnergy_HH
