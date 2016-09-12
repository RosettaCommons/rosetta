// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/Edge.hh
/// @brief  Fold tree edge class
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_Edge_hh
#define INCLUDED_core_kinematics_Edge_hh


// Unit headers
#include <core/kinematics/Edge.fwd.hh>

// In a perfect world, I would ifdef for Windows PyRosetta. A perfect world this isn't.
#include <string>
#include <utility>


namespace core {
namespace kinematics {

///////////////////////////////////////////////////////////////////////////////
/// \brief single edge of the fold_tree
///
/// an edge is a path between two vertices(start and end residues). it can be
/// either a continuous segment like a normal piece of polymer ("PEPTIDE" edge,
/// index label as "-1"), a chemical connection between two residues ("CHEMICAL
/// edge), or a rigid-body transformation between two residues ("JUMP" edge,
/// index label as "1", "2",...). The edge is the basic unit of the fold tree
/// as it stores info on how to build coordinates of the end residue given that
/// of the starting residue and degrees of freedom between these two
/// vertices.
class Edge
{
public:

	/// APL -- CODE DUPLICATION -- FIX THIS IN A BETTER WAY TO RESOLVE THE CIRCULAR DEPENDENCY
	const static int PEPTIDE{-1}; // must be negative, see Edge::is_jump()
		const static int CHEMICAL{-2}; // for fold-tree edges that connect two chemically-bound residues

		public:

		/////////////////////////////////////////////////////////////////////////////
		// member access

		/// @brief start vertex, return by value
		inline
			int
			start() const
			{
			return start_;
	}

	/// @brief start vertex, return by reference
	inline
	int &
	start()
	{
		return start_;
	}


	/// @brief stop vertex, return by value
	inline
	int
	stop() const
	{
		return stop_;
	}

	/// @brief stop vertex, return by reference
	inline
	int &
	stop()
	{
		return stop_;
	}


	/// @brief start_atom, return by value
	inline
	std::string
	start_atom() const
	{
		return start_atom_;
	}

	/// @brief start atom, return by reference
	inline
	std::string&
	start_atom()
	{
		return start_atom_;
	}


	/// @brief stop_atom, return by value
	inline
	std::string
	stop_atom() const
	{
		return stop_atom_;
	}

	/// @brief stop_atom, return by reference
	inline
	std::string &
	stop_atom()
	{
		return stop_atom_;
	}

	/// @brief start-atom, alt name, return by value
	inline
	std::string
	upstream_atom() const
	{
		return start_atom_;
	}

	/// @brief start-atom, alt name, return by reference
	inline
	std::string &
	upstream_atom()
	{
		return start_atom_;
	}

	/// @brief stop-atom, alt name, return by value
	inline
	std::string
	downstream_atom() const
	{
		return stop_atom_;
	}

	/// @brief stop-atom, alt name, return by reference
	inline
	std::string &
	downstream_atom()
	{
		return stop_atom_;
	}

	/// @brief label (edge type), return by value
	inline
	int
	label() const
	{
		return label_;
	}

	/// @brief label (edge type), return by reference
	inline
	int &
	label()
	{
		return label_;
	}


	// properties

	/// @brief edge is a jump?
	inline
	bool
	is_jump() const
	{
		return ( label_ > 0 );
	}


	inline
	bool
	is_chemical_bond() const
	{
		return ( label_ == CHEMICAL );
	}

	/// @brief  Edge is peptide edge?
	inline
	bool
	is_polymer() const
	{
		return ( label_ == PEPTIDE );
	}

	/// @brief  Edge is peptide edge?
	/// deprecated
	inline
	bool
	is_peptide() const
	{
		return ( label_ == PEPTIDE );
	}

	/// @brief edge has start and stop atoms?
	inline
	bool
	has_atom_info() const
	{
		return ( start_atom_.size() && stop_atom_.size() );
	}

	inline
	bool
	keep_stub_in_residue() const
	{
		return bKeepStubInResidue_;
	}

	inline
	bool&
	keep_stub_in_residue()
	{
		return bKeepStubInResidue_;
	}

	// only one use in all the code
	// returns 1 if start<stop, -1 if stop<start, dies if jump or a chemical edge
	//
	/// @brief direction for a continuous-segement edge. 1 if start residue number < stop residue number and -1 otherwise
	inline
	int
	polymer_direction() const
	{
		debug_assert( label_ == PEPTIDE );
		return ( start_ < stop_ ? 1 : -1 );
	}


	/// @brief  Is this edge valid (false for default-constructed edges)
	inline
	bool
	valid() const
	{
		return ( start_ > 0 && stop_ > 0 && label_ != 0 );
	}

	/////////////////////////////////////////////////////////////////////////////
	// construction

	///  @brief default constructor
	Edge():
		start_(0),
		stop_(0),
		label_(0),
		start_atom_(""),
		stop_atom_(""),
		bKeepStubInResidue_( false )
	{}

	/// @brief constructor without atomno info
	Edge( int const start_in, int const stop_in, int const label_in):
		start_( start_in ),
		stop_( stop_in ),
		label_( label_in ),
		start_atom_(""),
		stop_atom_(""),
		bKeepStubInResidue_( false )
	{}


	/// @brief CHEMICAL Edge constructor (requires atomno info) -- note: a chemical
	/// edge may be built from any constructor, this one is for convenience only
	Edge( int const start_in, int const stop_in, std::string  start_atom, std::string  stop_atom ):
		start_( start_in ),
		stop_( stop_in ),
		label_( CHEMICAL ),
		start_atom_(std::move( start_atom )),
		stop_atom_(std::move( stop_atom )),
		bKeepStubInResidue_( false )
	{}

	/// @brief JUMP Edge constructor (requires atomno info) -- note: a chemical
	/// edge may be built from any constructor, this one is for convenience only
	Edge( int const start_in, int const stop_in, int label,
		std::string  start_atom, std::string  stop_atom,
		bool bKeepStubInResidue ):
		start_( start_in ),
		stop_( stop_in ),
		label_( label ),
		start_atom_(std::move( start_atom )),
		stop_atom_(std::move( stop_atom )),
		bKeepStubInResidue_( bKeepStubInResidue )
	{}

	// stream I/O ////////////////////////
	// these two should be inverses!

	/// @brief input operator
	friend std::istream & operator >>(std::istream & is, Edge & e);

	/// @brief output operator
	friend std::ostream & operator <<(std::ostream & os, const Edge & e);

	/// @brief less than operator
	friend bool operator <( Edge const & a, Edge const & b );

	/// @brief equal to operator
	friend bool operator ==( Edge const & a, Edge const & b );

	/// @brief not equal to operator
	friend bool operator !=( Edge const & a, Edge const & b );

private:

	///////
	// data

	/// start vertex (residue)
	int start_;
	/// stop vertex (residue)
	int stop_;
	/// type of the edge, continuous segement(-1) or rigid-body jump(1,2,...)
	int label_;
	/// start atom
	std::string start_atom_;
	/// stop atom
	std::string stop_atom_;

	/// STUB Info for jumps
	bool bKeepStubInResidue_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // Edge

} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_Edge_HH
