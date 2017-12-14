// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/DOF_Node.hh
/// @brief  Kinematics
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_DOF_Node_hh
#define INCLUDED_core_optimization_DOF_Node_hh

// Unit headers
#include <core/optimization/DOF_Node.fwd.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>

#include <core/types.hh> // Vector

// // Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {


class DOF_Node : public utility::pointer::ReferenceCount
{
public:
	//typedef numeric::xyzVector< Real > Vector;
	typedef id::AtomID AtomID;
	typedef utility::vector1< AtomID > AtomIDs;
	typedef id::DOF_ID DOF_ID;
	typedef id::TorsionID TorsionID;
	typedef id::DOF_Type DOF_Type;


public:
	inline
	std::string
	to_string() const {
		std::ostringstream oss;
		oss << id << " f1=" << F1_[0]<<","<< F1_[1]<<","<< F1_[2]<<"  f2="
			<<F2_[0]<<","<< F2_[1]<<","<< F2_[2];
		return oss.str();
	}

	inline
	Vector &
	F1() { return F1_; };

	inline
	Vector &
	F2() { return F2_; }

	inline
	Vector const &
	F1() const { return F1_; };

	inline
	Vector const &
	F2() const { return F2_; }

	inline
	int
	rsd() const { return id.atom_id().rsd(); }

	inline
	int
	atomno() const { return id.atom_id().atomno(); }

	inline
	AtomID const &
	atom_id() const { return id.atom_id(); }

	inline
	DOF_Type
	type() const { return id.type(); }

	inline
	DOF_ID const &
	dof_id() const { return id; }

	inline
	int depth() const;

	inline
	AtomIDs const &
	atoms() const
	{
		return atoms_;
	}

	inline
	DOF_NodeCOP
	parent() const
	{
		return parent_;
	}

	inline
	void
	clear_atoms() {
		atoms_.clear(); // don't deallocate space -- makes DOF_Nodes reusable.
	}


	inline
	void
	add_atom( AtomID const & atom )
	{
		atoms_.push_back( atom );
	}


	/// get the rosetta torsion id for this DOF
	/**
	This may not exist, of course. But it's useful to know what it
	is when calculating derivatives of terms like rama/dunbrack/paa
	**/
	TorsionID const &
	torsion_id() const
	{
		return torsion_id_;
	}


	/// set the rosetta torsion id for this DOF
	/**
	This may not exist, of course. But it's useful to know what it
	is when calculating derivatives of terms like rama/dunbrack/paa
	**/
	void
	torsion_id(
		id::TorsionID const & id_in
	)
	{
		torsion_id_ = id_in;
	}


	/// sum derivative contributions down the tree
	inline
	void
	link_vectors()
	{
		if ( parent_ ) {
			parent_->F1() += F1_;
			parent_->F2() += F2_;
		}
	}

	// constructor
	DOF_Node(
		DOF_ID const & id_in,
		DOF_NodeOP parent_in
	):
		utility::pointer::ReferenceCount(),
		F1_(0.0),
		F2_(0.0),
		depth_(-1),
		id( id_in ),
		parent_(std::move( parent_in )),
		torsion_id_( id::TorsionID::BOGUS_TORSION_ID() ),
		dependent_( false )
	{}

	void
	set_id( DOF_ID const & setting ) {
		id = setting;
	}

	void
	set_parent( DOF_NodeOP setting )
	{
		debug_assert( setting.get() != this ); // an object in an OP should never point to itself
		parent_ = setting;
	}

	bool
	dependent() const { return dependent_; }

	void
	dependent( bool const setting ) { dependent_ = setting; }

private:
	Vector F1_;
	Vector F2_;
	mutable int depth_;
	DOF_ID id;
	AtomIDs atoms_;
	DOF_NodeOP parent_;
	TorsionID torsion_id_;
	bool dependent_; // only used/set/checked in symmetric minimization

public:

	friend
	inline
	bool
	operator< ( DOF_Node const & t1, DOF_Node const & t2 ) {
		return ( t1.depth() > t2.depth() ); // check that this gives correct sort
	}

}; // DOF_Node


inline
int
DOF_Node::depth() const
{
	if ( parent_ == nullptr ) {
		depth_ = 0;
	} else if ( depth_ < 0 ) {
		depth_ = parent_->depth() + 1;
	}
	debug_assert( depth_ >= 0 );
	return depth_;
}


} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_min_HH
