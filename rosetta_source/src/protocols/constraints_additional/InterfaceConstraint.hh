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

#ifndef INCLUDED_protocols_constraints_additional_InterfaceConstraint_HH
#define INCLUDED_protocols_constraints_additional_InterfaceConstraint_HH

#include <protocols/constraints_additional/InterfaceConstraint.fwd.hh>

#include <protocols/scoring/InterfaceInfo.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>

//Auto Headers
#include <numeric/xyzVector.hh>



// C++ Headers
//#include <cstdlib>
//#include <iostream>
//#include <map>
//#include <utility>


namespace protocols {
namespace constraints_additional {


///

class InterfaceConstraint : public core::scoring::constraints::Constraint {
public:

	// default c-tor
	InterfaceConstraint() : core::scoring::constraints::Constraint( core::scoring::atom_pair_constraint ) {}

	///c-tor
	InterfaceConstraint(
		AtomID const & a1,
	 	core::scoring::constraints::FuncOP func,
		core::scoring::ScoreType scoretype = core::scoring::atom_pair_constraint
	):
		core::scoring::constraints::Constraint( scoretype ),
		atom1_( a1 ),
		func_( func )
	{}

	virtual core::scoring::constraints::ConstraintOP clone() const {
		return new InterfaceConstraint( atom1_, func_, score_type() );
	}

	///
	virtual
	core::scoring::constraints::ConstraintOP clone( core::scoring::constraints::FuncOP func ) const {
		return new InterfaceConstraint( atom1_, func, score_type() );
	}


	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual
	core::scoring::constraints::ConstraintOP remapped_clone(
		core::pose::Pose const & src,
		core::pose::Pose const & dest,
		core::id::SequenceMappingCOP map = NULL
	) const;

	virtual void
	setup_for_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction const & ) const;

	///
	void
	score( core::scoring::constraints::XYZ_Func const &, core::scoring::EnergyMap const &, core::scoring::EnergyMap & emap ) const
	{
		emap[ this->score_type() ] += func_->func( dist_ );
	}

	// atom deriv
	virtual
	void
	fill_f1_f2(
		core::id::AtomID const & atom,
		core::scoring::constraints::XYZ_Func const &,
		core::Vector & F1,
		core::Vector & F2,
		core::scoring::EnergyMap const & weights
	) const;

	///
	Size
	natoms() const
	{
		return 1;
	}

	virtual
	core::scoring::constraints::ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	///
	AtomID const &
	atom( Size const n ) const
	{
		switch( n ) {
		case 1:
			return atom1_;
		default:
			utility_exit_with_message( "InterfaceConstraint::atom() bad argument" );
		}
		return atom1_;
	}


	std::string type() const {
		return "Interface";
	}

	void show( std::ostream& out ) const;
	void show_def( std::ostream& out, core::pose::Pose const& pose ) const;

	void read_def( std::istream& in, core::pose::Pose const& pose, core::scoring::constraints::FuncFactory const& func_factory );

	protocols::scoring::InterfaceInfo const & interface_from_pose( core::pose::Pose const & ) const;
	protocols::scoring::InterfaceInfo & nonconst_interface_from_pose( core::pose::Pose & ) const;


private:
	utility::vector1 < Size > exclude_chains_;
	AtomID atom1_;
	core::scoring::constraints::FuncOP func_;
	mutable AtomID atom2_;
	mutable core::Real dist_;
	mutable core::Vector f1_, f2_;
}; // class InterfaceConstraint

} // constraints
} // protocols

#endif
