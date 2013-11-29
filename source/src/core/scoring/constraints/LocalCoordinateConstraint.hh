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

// AUTO-REMOVED #include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/id/NamedStubID.hh>

// C++ Headers
#include <cstdlib>
#include <iostream>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>



//#include <map>

namespace core {
namespace scoring {
namespace constraints {

class LocalCoordinateConstraint;
typedef utility::pointer::owning_ptr< LocalCoordinateConstraint > LocalCoordinateConstraintOP;
typedef utility::pointer::owning_ptr< LocalCoordinateConstraint const > LocalCoordinateConstraintCOP;

///

class LocalCoordinateConstraint : public Constraint {
public:


	LocalCoordinateConstraint() :
		Constraint( coordinate_constraint ),
		atom_( id::BOGUS_ATOM_ID ),
		fixed_stub_( id::BOGUS_STUB_ID ),
		func_( NULL ) {}

	///c-tor
	LocalCoordinateConstraint(
	  id::AtomID const & a1,
		id::StubID const & fixed_stub_in,
		Vector const & xyz_target_in,
	 	FuncOP func,
		ScoreType scotype = coordinate_constraint
	):
		Constraint( scotype ),
		atom_(a1),
		fixed_stub_( fixed_stub_in ),
		xyz_target_( xyz_target_in ),
		func_( func )
	{
		runtime_assert( fixed_stub_.atom( 1 ) == fixed_stub_.center() || !fixed_stub_.center().valid() );
		//don't allow 4-atom stubs, because that changes other functions
	}

	~LocalCoordinateConstraint() {};

	virtual std::string type() const {
		return "LocalCoordinateConstraint";
	}

	virtual ConstraintOP clone() const {
		return new LocalCoordinateConstraint( *this );
	}

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP map=NULL ) const;

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

 	///
	void show( std::ostream& out ) const
	{
		out << "LocalCoordinateConstraint ("
				<< atom_.atomno() << "," << atom_.rsd() << "-"
				<< fixed_stub_ << ")" << std::endl;
		func_->show( out );
	}

	// @brief Reads the definition of a Constraint from the given std::istream,
	// using the given Pose, and the given FuncFactory. This method is intended
	// to be overridden by derived classes if they'd like to use the
	// ConstraintIO machinery.
	virtual void read_def( std::istream &, pose::Pose const &, FuncFactory const & );

	///
	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	// @brief take coordinates, distances, angles, etc from given pose
	///
	virtual void steal_def( pose::Pose const& );

	using Constraint::score;

	///
	Real
	score(
		Vector const & xyz, //target
		Vector const & s1, //fixed_stub.a
		Vector const & s2, //fixed_stub.b
		Vector const & s3 //fixed_stub.c
	) const;

	///
	void
	score( XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
	{
		emap[ this->score_type() ] += score( xyz( atom_ ),
			xyz( fixed_stub_.atom( 1 ) ),
			xyz( fixed_stub_.atom( 2 ) ),
			xyz( fixed_stub_.atom( 3 ) )
		);
	}

	// atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		XYZ_Func const & xyz,
		Vector & F1,
 		Vector & F2,
		EnergyMap const & weights
	) const
	{
		utility_exit_with_message( " derivative of LocalCoordinateConstraint not supported yet " );
		if ( atom != atom_ ) return;

		Vector const & xyz1( xyz( atom_ ) ), xyz2( xyz_target_ );

		Vector const f2( xyz1 - xyz2 );
		Real const dist( f2.length() ), deriv( dfunc( dist ) );
		if ( deriv != 0.0 && dist != 0.0 ) {
			Vector const f1( xyz1.cross( xyz2 ) );
			// jk: double F1 and F2 because the target is fixed
			// (matches deriv_check, and minimizes faster)
			// rhiju: No, JK, this isn't working...
			F1 += ( ( deriv / dist ) * f1 ) * weights[ this->score_type() ];
			F2 += ( ( deriv / dist ) * f2 ) * weights[ this->score_type() ];
		}
	}



	///
	Size
	natoms() const
	{
		return 4;
	}

	///
	AtomID const &
	atom( Size const n ) const
	{
		if ( n == 1 ) {
			return atom_;
		} else if ( n <= 4 ) {
			return fixed_stub_.atom( n - 1 );
		} else {
			utility_exit_with_message( "LocalCoordinateConstraint::atom() bad argument" );
		}
		return atom_;
	}

	Real
	dist( pose::Pose const & pose ) const;

	virtual Size show_violations(
		std::ostream& out,
		pose::Pose const& pose,
		Size verbose_level,
		Real threshold = 1
	) const;

	void set_fixed_stub( id::StubID new_stub ) {
		fixed_stub_ = new_stub;
	}

	Vector xyz_target( core::pose::Pose const& local_frame_pose ) const;

	void set_xyz_target( Vector const& xyz_in, core::pose::Pose const& local_frame_pose );

private:

	// functions
	Real
	func( Real const theta ) const
	{
		return func_->func( theta );
	}

	// deriv
	Real
	dfunc( Real const theta ) const
	{
		return func_->dfunc( theta );
	}

private:
	// data
	id::AtomID atom_;
	id::StubID fixed_stub_;
	Vector xyz_target_;
	FuncOP func_;
};

}
}
}

#endif
