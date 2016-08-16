// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/constraints/MultiConstraint.hh
/// @brief contains declarations for a type of constraint that holds multiple
/// other constrains that belong to each other and are all evaluate at once
/// @author Florian Richter (floric@u.washington.edu, march 2008)


#ifndef INCLUDED_core_scoring_constraints_MultiConstraint_hh
#define INCLUDED_core_scoring_constraints_MultiConstraint_hh

#include <core/scoring/constraints/MultiConstraint.fwd.hh>

// Unit header
#include <core/scoring/constraints/Constraint.hh>

#include <core/scoring/ScoreType.hh>

#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
//Utility Headers
#include <numeric/xyzVector.fwd.hh>

// STL Headers
#include <map>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


class MultiConstraint : public Constraint {
public:

	/// @brief default Constructor
	MultiConstraint( ScoreType const & t = dof_constraint );

	/// @brief Constructor -- performs a shallow copy of the input ConstraintCOPs.
	MultiConstraint( ConstraintCOPs const & cst_in, ScoreType const & t = dof_constraint );

	/// @brief Creates a deep copy of this %MultiConstaint, cloning all constraints that it holds.
	virtual
	ConstraintOP clone() const;

	virtual
	MultiConstraintOP empty_clone() const;

	/// @brief number of atoms involved in this MultiConstraint container
	Size natoms() const;

	/// @brief number of constraints that are held by this %MultiConstraint
	Size size() const;

	virtual std::string type() const;

	/// @brief read in constraint defiinition
	virtual
	void read_def( std::istream& data, pose::Pose const& pose,func::FuncFactory const& func_factory );

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	virtual
	bool operator == ( Constraint const & other ) const;

	virtual
	bool same_type_as_me( Constraint const & other ) const;

	/// @brief compute score
	virtual
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;

	virtual
	AtomID const & atom( Size const n ) const{
		debug_assert( n <= member_atoms_.size() );
		return member_atoms_[n];
	}

	virtual
	utility::vector1< core::Size >
	residues() const { return member_residues_; }

	//@brief translates the atom-names into numbers
	virtual void setup_for_scoring( func::XYZ_Func const &, ScoreFunction const & ) const;

	/// @brief add individual constraint into MultiConstraint
	virtual
	void
	add_individual_constraint( ConstraintCOP cst_in );

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	/// @brief compute atom deriv
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	virtual
	void show( std::ostream& out ) const;

	virtual
	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	virtual
	Size show_violations( std::ostream & out, pose::Pose const & pose, Size verbose_level, Real threshold = 1.0 ) const;

	ConstraintCOPs const &
	member_constraints() const {
		return member_constraints_;
	}

	virtual ConstraintOP remapped_clone(
		pose::Pose const& /*src*/,
		pose::Pose const& /*dest*/,
		id::SequenceMappingCOP map=NULL ) const;

	void set_effective_sequence_separation( core::Size setting ) {
		report_this_as_effective_sequence_separation_ = setting;
	}

	virtual
	core::Size choose_effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp, numeric::random::RandomGenerator& );

	virtual
	core::Size effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& ) const {
		return report_this_as_effective_sequence_separation_;
	}

protected:
	/// @brief Explicit copy constructor so that derived classes will recieve a deep copy
	/// of the member constraints held in this class.
	MultiConstraint( MultiConstraint const & src );

	/// @brief Return a vector of Constraints that are clones of the member constraints.
	ConstraintCOPs cloned_member_constraints() const;

private:

	//vector that holds the constraints
	ConstraintCOPs member_constraints_;

	//data structure that holds the atoms and atom numbers
	utility::vector1< core::Size > member_residues_;
	utility::vector1< AtomID > member_atoms_;
	std::map< AtomID, ConstraintCOPs > atomid_to_csts_;


	core::Size report_this_as_effective_sequence_separation_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //MultiConstraint

} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_MultiConstraint )
#endif // SERIALIZATION


#endif
