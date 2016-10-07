// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/symmetry/SymmetricScoreFunction.hh
/// @brief  Symmetric Score function class
/// @author Ingemar Andre


#ifndef INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_hh
#define INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_hh

// Unit headers
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace symmetry {

class SymmetricScoreFunction : public ScoreFunction
{
public:
	typedef ScoreFunction parent;
	typedef conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef conformation::symmetry::SymmetryInfoCOP SymmetryInfoCOP;

public:

	/// ctor
	SymmetricScoreFunction();

	/// @brief Constructor initializing base class data from a (possibly local) options collection
	SymmetricScoreFunction( utility::options::OptionCollection const & );

private:

	// Copy constructors and assignment operators may discard subtype information.
	// Do not use or implement
	SymmetricScoreFunction &
	operator=( SymmetricScoreFunction const & );

	SymmetricScoreFunction( SymmetricScoreFunction const & );

public:

	/// @brief NOT FOR GENERAL USE
	/// Copy the information about src into the current score function.
	/// There are deep interactions with subclasses,
	/// (the subclass information doesn't necessarily get copied)
	/// so this is primarily for advanced scorefunction manipulation.
	/// Normal usage should just use clone() and replace the OP.
	virtual void
	assign( ScoreFunction const & src);

	/// @brief NOT FOR GENERAL USE
	virtual void
	assign( SymmetricScoreFunction const & src);

	ScoreFunctionOP clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// score
	/////////////////////////////////////////////////////////////////////////////

	virtual Real
	operator ()( pose::Pose & pose ) const;

	/// @brief Initialize a MinimizationGraph and cache it in the pose's Energies object
	/// for use during minimization -- only add edges to the asymmetric unit and within it
	/// are added to the MinimizationGraph.
	virtual
	void
	setup_for_minimizing(
		pose::Pose & pose,
		kinematics::MinimizerMapBase const & min_map
	) const;


	void
	eval_twobody_neighbor_energies( pose::Pose & pose ) const;

	void
	eval_long_range_twobody_energies( pose::Pose & pose ) const;


	void
	eval_onebody_energies( pose::Pose & pose ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose ) const;

	virtual
	void
	eval_npd_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		Vector & F1,
		Vector & F2
	) const;

	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose
	) const;

	//void create_intersubunit_hbonds( pose::Pose & pose, hbonds::HBondSetOP hbond_set_subunit ) const;

	void
	intersubunit_hbond_energy( pose::Pose & pose, EnergyMap & intersubunit_energy  ) const;

	void
	symmetrical_allow_hbonds( pose::Pose & pose ) const;

	void
	set_symmetric_residue_neighbors_hbonds( pose::Pose & pose ) const;

	void
	set_symmetric_cenlist( pose::Pose & pose ) const;

	void
	correct_arrays_for_symmetry( pose::Pose & pose ) const;

	void
	correct_finalize_score( pose::Pose & pose ) const;

	static
	void
	list_options_read( utility::options::OptionKeyList & options_read );

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

inline core::scoring::ScoreFunctionOP
asymmetrize_scorefunction( ScoreFunction const & src ) {
	//Should there be some test to see if the ScoreFunction is actually symmetric first?
	return src.clone_as_base_class();
}

inline core::scoring::symmetry::SymmetricScoreFunctionOP
symmetrize_scorefunction( ScoreFunction const & src ) {
	//Should there be some test to see if the ScoreFunction is actually symmetric first?
	SymmetricScoreFunctionOP new_scorefxn( new SymmetricScoreFunction );
	new_scorefxn->assign( src );
	return new_scorefxn;
}

} // symmetry
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_symmetry_SymmetricScoreFunction )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_SymmetricScoreFunction_HH
