// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Helper class for FoldConstraints Protocol, filters constraints by sequence separation
/// @details
///
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_constraints_additional_MaxSeqSepConstraintSet_hh
#define INCLUDED_protocols_constraints_additional_MaxSeqSepConstraintSet_hh


// Unit Headers
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.fwd.hh>

// Package Headers


// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/MinimizationData.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace constraints_additional {

class MaxSeqSepConstraintSet : public core::scoring::constraints::ConstraintSet
{
public:
	typedef ConstraintSet Parent;
	typedef core::scoring::constraints::ConstraintSetOP ConstraintSetOP;

public:
	MaxSeqSepConstraintSet( ConstraintSet const & other, core::kinematics::FoldTree const&f );

	MaxSeqSepConstraintSet( MaxSeqSepConstraintSet const &other );
	~MaxSeqSepConstraintSet( );
protected:
	MaxSeqSepConstraintSet( ConstraintSet const &other, core::kinematics::ShortestPathInFoldTreeCOP sp );

public:

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

	virtual ConstraintSetOP remapped_clone(
		core::pose::Pose const& src,
		core::pose::Pose const& dest,
		core::id::SequenceMappingCOP smap = NULL
	) const;

	Size max_seq_sep () const {
		return max_seq_sep_;
	}

	void set_max_seq_sep( Size setting ) {
		max_seq_sep_=setting;
		mark_revision_id_expired(); //force recompute of cached energies
	}


	void
	residue_pair_energy(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn,
		core::scoring::EnergyMap & emap
	) const;

	bool too_far( int const pos1, int const pos2 ) const;

	Size
	largest_possible_sequence_sep( core::pose::Pose const& pose ) const;


	bool
	residue_pair_constraint_exists( int const pos1, int const pos2 ) const
	{
		if ( too_far( pos1, pos2 ) ) return false;
		return Parent::residue_pair_constraint_exists( pos1, pos2 );
	}


	/// @brief If a pair of residues are not too far, then let the
	/// parent class implementation add the constraints for this
	/// respair to the respair_data_cache.  This is the only opportunity
	/// the MaxSeqSepConstraint set has to influence how constraints
	/// are evaluated during minimization.
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

	using Parent::show_violations;

	Size
	show_violations( std::ostream& out, core::pose::Pose&, Size verbose_level, core::Real threshold = 1  );

	core::kinematics::ShortestPathInFoldTree const&
	shortest_path() const {
		return *shortest_path_;
	}

	void
	virtual eval_non_residue_pair_energy(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const;


protected:

	/*
	void
	eval_atom_derivative_for_residue_pairs(
	core::id::AtomID const & atom_id,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & weights,
	core::Vector & F1,
	core::Vector & F2
	) const;*/


private:
	MaxSeqSepConstraintSet(); // private default constructor for use in the call to detached_clone

	Size max_seq_sep_;
	core::kinematics::ShortestPathInFoldTreeCOP shortest_path_;
#ifdef    SERIALIZATION
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // constraints_additional
} // protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_constraints_additional_MaxSeqSepConstraintSet )
#endif // SERIALIZATION


#endif
