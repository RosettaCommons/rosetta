// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPResidueLipophilicityEnergy.hh
///
/// @brief  Fullatom Smoothed Membrane Environment Energy
/// @details residue speicific enrgy by membrane depth, according to the Elazar
/// hydrophobicity scale
///    @FlesihmanLab.
///    Last Modified: 4/4/16
///
/// @author  Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
/// @author  Assaf Elazar
/// @author Sarel Fleishman

#ifndef INCLUDED_core_scoring_membrane_MPResidueLipophilicityEnergy_hh
#define INCLUDED_core_scoring_membrane_MPResidueLipophilicityEnergy_hh

// Unit Headers
#include <core/scoring/membrane/MPResidueLipophilicityEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
//#include <utility/graph/Graph.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>
#include <map>
#include <numeric/interpolation/spline/CubicSpline.hh>

namespace core {
namespace scoring {
namespace membrane {


/// @brief Fullatom Membrane Solvation Energy - Statistically Derived,
/// and smoothed derivatives
class MPResidueLipophilicityEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy  {

public:

	typedef ContextDependentOneBodyEnergy  parent;

public:
	// Constructors ////////////////////

	/// @brief Default Constructor
	MPResidueLipophilicityEnergy();

	/// @brief Create a clone of this energy method
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	// Scoring Methods ////////////////

	/// @brief Computes dScore/dNumNeighbors for all residues for rapid use in later
	/// atom derivate calculations
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;

	/// @brief Causes a neighbor graph update
	void
	setup_for_derivatives( pose::Pose &, ScoreFunction const & ) const override;

	/// @brief Evaluates the one-body energy for a residue
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		EnergyMap &
	) const override;

	utility::vector1< core::Size >
	centroid_neighbors(
		pose::Pose const & pose,
		conformation::Residue const & rsd
	) const;

	/// @brief Increments the F1 and F2 derivative vectors for an atom
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const &,
		Vector &,
		Vector &
	) const override;

	/// @brief Unused by the MPResidueLipophilicityEnergy class, returns 0
	Distance
	atomic_interaction_cutoff() const;

	core::Real
	report_ressolv(
		std::ostream & out,
		pose::Pose const & pose,
		bool print_splines
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override;

	/// @brief Tells the scoring function to maintain the TwelveANeighborGraph
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

	/// @brief convert the neighbor counts for a residue and its aa type into
	/// a score and a score derivative.
	void
	calc_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		Real & score
	) const;

	void
	parse_elazar_polynom_table() ;

	void
	parse_elazar_spline_table() ;

	numeric::interpolation::spline::CubicSpline
	split_line_to_spline(
		std::string line_
	) const;

	utility::vector1< Real >
	split_line_to_polyval(
		std::string line_
	) const;

	Real
	polynom_by_z(
		char const & res,
		Real const & z,
		bool & full_atom
	) const;

	Real
	spline_by_z(
		char const & res,
		Real const & z,
		bool & full_atom
	) const;

	utility::vector1 < core::Size >
	neighboring_atoms(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		Real const & cutoff_1,
		Real const & cutoff_2
	) const;

	Real
	burial_sigmoid(
		core::Size n_atoms,
		Real const slope,
		Real const offset
	) const;

	core::Real
	calc_residue_burial(
		pose::Pose const & pose,
		conformation::Residue const & rsd
	) const;


private:

	std::map< char,  utility::vector1< Real > > res_polynom_map_;
	std::map< char,  utility::vector1< Real > > cen_res_polynom_map_;

	std::map< char,  numeric::interpolation::spline::CubicSpline > res_spline_map_;
	std::map< char,  numeric::interpolation::spline::CubicSpline > cen_res_spline_map_;

	core::Size version() const override;

};

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPResidueLipophilicityEnergy_hh
