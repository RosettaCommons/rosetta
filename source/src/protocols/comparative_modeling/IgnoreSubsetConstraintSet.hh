// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_IgnoreSubsetConstraintSet_hh
#define INCLUDED_protocols_comparative_modeling_IgnoreSubsetConstraintSet_hh

#include <protocols/comparative_modeling/IgnoreSubsetConstraintSet.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// C++ headers
#include <set>
#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace comparative_modeling {

class IgnoreSubsetConstraintSet : public core::scoring::constraints::ConstraintSet {
public:
	typedef core::scoring::constraints::ConstraintSetOP ConstraintSetOP;
	typedef core::scoring::constraints::ConstraintSet ConstraintSet;

public:
	//IgnoreSubsetConstraintSet( ConstraintSet const & other );

	IgnoreSubsetConstraintSet( IgnoreSubsetConstraintSet const &other );

	IgnoreSubsetConstraintSet(
		std::set< int > residues_to_ignore,
		ConstraintSet const & other
	);

	virtual
	ConstraintSet const &
	operator = ( ConstraintSet const & rhs );

	virtual
	ConstraintSetOP
	clone() const;

	virtual
	void
	detached_copy( ConstraintSet const & src );

	virtual
	ConstraintSetOP
	detached_clone() const;

	virtual
	bool
	same_type_as_me( ConstraintSet const & other, bool recurse = true ) const;

	void
	residue_pair_energy(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn,
		core::scoring::EnergyMap & emap
	) const;

	/// @brief Allow the parent class implementation to add the residue constraints
	/// for this residue to the res_data_cache if this residue is not being ignored.
	virtual
	void
	setup_for_minimizing_for_residue(
		core::conformation::Residue const & rsd,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::kinematics::MinimizerMapBase const & minmap,
		core::scoring::ResSingleMinimizationData & res_data_cache
	) const;


	/// @brief Allow the parent class implenetation to add the residue-pair constraints
	/// for this residue pair to the respair_data_cache if neither residues are being ignored.
	virtual
	void
	setup_for_minimizing_for_residue_pair(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::kinematics::MinimizerMapBase const & minmap,
		core::scoring::ResSingleMinimizationData const & res1_data_cache,
		core::scoring::ResSingleMinimizationData const & res2_data_cache,
		core::scoring::ResPairMinimizationData & respair_data_cache
	) const;

	/// @brief Returns true if we're supposed to ignore this sequence position,
	/// false otherwise.
	bool ignore( int const pos ) const;

	void ignore_residue( int const pos );

	std::set< int > ignore_list() const;


	bool
	residue_pair_constraint_exists( int const pos1, int const pos2 ) const
	{
		if ( ignore(pos1) || ignore(pos2) ) return false;
		return ConstraintSet::residue_pair_constraint_exists( pos1, pos2 );
	}

protected:

	/*void
	eval_atom_derivative_for_residue_pairs(
	core::id::AtomID const & atom_id,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & weights,
	core::Vector & F1,
	core::Vector & F2
	) const;*/

private:
	std::set< int > ignore_list_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	IgnoreSubsetConstraintSet();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class IgnoreSubsetConstraintSet

} // comparative_modeling
} // protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_comparative_modeling_IgnoreSubsetConstraintSet )
#endif // SERIALIZATION


#endif
