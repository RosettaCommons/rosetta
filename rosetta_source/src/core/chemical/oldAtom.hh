// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/Atom.hh
/// @brief  Energy graph class declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_chemical_Atom_hh
#define INCLUDED_core_chemical_Atom_hh

// Unit Headers
#include <core/chemical/Atom.fwd.hh>
#include <core/chemical/AtomICoor.hh>

// Project Headers
//#include <core/chemical/EnergyMap.hh>
#include <core/graph/Graph.hh>
#include <core/graph/ArrayPool.hh>
//#include <core/chemical/ScoreType.hh>
#include <numeric/xyzVector.hh>
// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>



namespace core {
namespace chemical {

/// Class Atom holds the result of a domainmap update from the
/// Conformation object held by a pose; if the internal degrees of freedom
/// for a residue (corresponding to a node in this graph) have changed
/// (and are marked with color "0" in the domainmap), then the Atom
/// object will hold that information for the ScoringFunction to retrieve
class Atom : public graph::Node
{
public:
	typedef graph::Node parent;

public:
	Atom( graph::Graph * owner, Size index );
	Atom(
			Graph * owner,
			Size index,
			std::string const & name_in,
		//	std::string const type_name,
			std::string const mm_name,
			Size const atom_type_index,
			Size const mm_atom_type_index,
			Real const charge,
			Vector const ideal_xyz,
			AtomICoor const icoor

	);
	virtual ~Atom();
	virtual void copy_from( parent const * source );

	virtual void print() const;
	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

// Const Getters
	std::string const& name() const { return name_; };
	//std::string const& type_name() const { return type_name_; };
	std::string const& mm_name() const { return mm_name_; };
	Size const& atom_type_index() const { return atom_type_index_; };
	Size const& mm_atom_type_index() const { return mm_atom_type_index_; };
	Real const& charge() const { return charge_; };
	Vector const& ideal_xyz() const { return ideal_xyz_; };
	AtomICoor const& icoor() const { return icoor_; };
// Non-const getters
	AtomICoor & icoor() { return icoor_; };
// Setters
	void name( std::string const & name ) { name_ = name; };
	//std::string const& type_name() const { return type_name_; };
	void mm_name( std::string const & name ) { mm_name_ = name; };
	void atom_type_index( Size const & atom_type_index ) { atom_type_index_ = atom_type_index; };
	void mm_atom_type_index( Size const & mm_atom_type_index ) { mm_atom_type_index_ = mm_atom_type_index; };
	void charge( Real const & charge ) { charge_ = charge; };
	void ideal_xyz( Vector const & ideal_xyz) { ideal_xyz_= ideal_xyz; };
	void icoor( AtomICoor const & icoor) { icoor_ = icoor; };


private:

	// Primary data
	std::string name_;
	//std::string const type_name_;

	// Secondary data
	std::string mm_name_;
	//	std::string const csd_atom_name_;

	Size atom_type_index_;
	/// MM atom-type index
	Size mm_atom_type_index_;
	Real charge_;
	Vector ideal_xyz_;
	AtomICoor icoor_;

};

} //namespace chemical
} //namespace core

#endif

