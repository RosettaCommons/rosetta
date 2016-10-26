// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#ifndef INCLUDED_core_scoring_constraints_LocalCoordinateConstraint_hh
#define INCLUDED_core_scoring_constraints_LocalCoordinateConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>


// C++ Headers
#include <cstdlib>
#include <iostream>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


//#include <map>
#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

class LocalCoordinateConstraint;
typedef utility::pointer::shared_ptr< LocalCoordinateConstraint > LocalCoordinateConstraintOP;
typedef utility::pointer::shared_ptr< LocalCoordinateConstraint const > LocalCoordinateConstraintCOP;


class LocalCoordinateConstraint : public Constraint {
public:
	LocalCoordinateConstraint();

	///c-tor
	LocalCoordinateConstraint(
		id::AtomID const & a1,
		id::StubID const & fixed_stub_in,
		Vector const & xyz_target_in,
		func::FuncOP func,
		ScoreType scotype = coordinate_constraint
	);

	~LocalCoordinateConstraint();

	virtual std::string type() const;

	virtual ConstraintOP clone() const;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP map=NULL ) const;

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	virtual bool operator == ( Constraint const & rhs ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	///
	void show( std::ostream& out ) const;

	// @brief Reads the definition of a Constraint from the given std::istream,
	// using the given Pose, and the given func::FuncFactory. This method is intended
	// to be overridden by derived classes if they'd like to use the
	// ConstraintIO machinery.
	virtual void read_def( std::istream &, pose::Pose const &,func::FuncFactory const & );


	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	// @brief take coordinates, distances, angles, etc from given pose
	///
	virtual void steal_def( pose::Pose const& );

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	Real
	score(
		Vector const & xyz, //target
		Vector const & s1, //fixed_stub.a
		Vector const & s2, //fixed_stub.b
		Vector const & s3 //fixed_stub.c
	) const;

	void
	score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const;

	virtual
	core::Real
	dist( core::scoring::func::XYZ_Func const & xyz ) const;

	// atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	///
	Size
	natoms() const;

	AtomID const &
	atom( Size const n ) const;

	virtual Size show_violations(
		std::ostream& out,
		pose::Pose const& pose,
		Size verbose_level,
		Real threshold = 1
	) const;

	void set_fixed_stub( id::StubID new_stub );

	Vector xyz_target( core::pose::Pose const& local_frame_pose ) const;

	void set_xyz_target( Vector const& xyz_in, core::pose::Pose const& local_frame_pose );

private:

	// functions
	Real
	func( Real const theta ) const;

	// deriv
	Real
	dfunc( Real const theta ) const;

protected:
	/// @brief Explicit copy constructor so that derived classes will recieve a deep copy
	/// of the Func this class contains.
	LocalCoordinateConstraint( LocalCoordinateConstraint const & src );

private:
	// data
	id::AtomID atom_;
	id::StubID fixed_stub_;
	Vector xyz_target_;
	func::FuncOP func_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_LocalCoordinateConstraint )
#endif // SERIALIZATION


#endif
