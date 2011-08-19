// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_frags_TorsionFragment_hh
#define INCLUDED_protocols_frags_TorsionFragment_hh

// Unit Headers
#include <protocols/frags/TorsionFragment.fwd.hh>

// Rosetta Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <utility/vector1_bool.hh>


/// ******************************************************************************************************
///  Code duplication alert:  TorsionFragment is going to be phased out
///   Please avoid writing any new code using these classes: look in core/fragment/ instead
///   look in /protocols/abinitio/FragmentMover  for usage examples
///
/// ******************************************************************************************************

// C++ Headers
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// #include <fstream>
// #include <sstream>

namespace protocols {
namespace frags {

using core::Real;
using core::Size;

///\brief a class for single piece of torsion fragment
///
///It stores torsion angles and secondary structure.
///Torsions are stored as vector of vector, such as TorsionFragment[frag_size][n_bb_torsion];
///SS are stored as vector, such as TorsionFragment[frag_length]
class TorsionFragment : public utility::pointer::ReferenceCount {

public:

	/// constructor, fragment_size (3mer or 9mer) and number of backbone torsions(3 for protein)
	TorsionFragment( Size const size_in, Size const nbb_in ):
		torsions_( size_in ),
		secstruct_( size_in )
	{
		for ( Size k=1; k<= size_in; ++k ) {
			torsions_[k].resize( nbb_in );
		}
	}

	/// fragment size, 3mer or 9mer?
	inline
	Size
	size() const
	{
		return torsions_.size();
	}

	/// number of backbone torsions.
	inline
	Size
	nbb() const
	{
		return ( torsions_.empty() ? 0 : torsions_[1].size() );
	}

	// insert this piece of fragment to a pose at position "begin"
	void
	insert( core::pose::Pose & pose, Size const begin ) const;

	/// set value for specific torsion in this piece of fragment.
	void
	set_torsion( Size const pos, Size const tor, Real const setting )
	{
		torsions_[pos][tor] = setting;
	}

	/// set secstruct for this position
	void
	set_secstruct( Size const pos, char const setting )
	{
		secstruct_[pos] = setting;
	}

	/// get value for specific torsion in this piece of fragment
	Real
	get_torsion( Size const pos, Size const tor ) const
	{
		return torsions_[pos][tor];
	}

	char
	get_secstruct( Size const pos ) const
	{
		return secstruct_[pos];
	}

private:
	/// torsion angles, the first dimension is fragment size, the second dimension is number of backbone torsions.
	utility::vector1< utility::vector1< Real > > torsions_;
	/// secstruct, dimensioned as fragment size
	utility::vector1< char > secstruct_;
};

typedef utility::pointer::owning_ptr< TorsionFragment > TorsionFragmentOP;

/////////////////////////////////////////////////////////////////////////////
///\brief a class for collection of fragments for a single residue position
///
///essentially a vector of TorsionFragment (owning pointers)
class SingleResidueTorsionFragmentLibrary : public utility::pointer::ReferenceCount {
/// ******************************************************************************************************
///  Code duplication alert:  TorsionFragment is going to be phased out
///   Please avoid writing any new code using these classes: look in core/fragment/ instead
///   look in /protocols/abinitio/FragmentMover  for usage examples
///
/// ******************************************************************************************************

public:
	/// insert one piece of fragment in the front
	void
	insert_fragment( TorsionFragmentOP fragment )
	{
		fragments_.insert( fragments_.begin(), fragment );
	}

	/// append one piece of fragment at the end
	void
	append_fragment( TorsionFragmentOP fragment )
	{
		fragments_.push_back( fragment );
	}

	/// number of available fragment pieces for this position
	Size
	size() const
	{
		return fragments_.size();
	}

	/// overloaded [] operator to get single piece of fragment (actual object, not owning pointers)
	TorsionFragment const &
	operator[]( Size const index ) const
	{
		return *( fragments_[ index ] );
	}


private:
	/// a collection of TorsionFragments (by owning pointers)
	utility::vector1< TorsionFragmentOP > fragments_;
};

typedef utility::pointer::owning_ptr< SingleResidueTorsionFragmentLibrary > SingleResidueTorsionFragmentLibraryOP;

/////////////////////////////////////////////////////////////////////////////
///\brief a class for classic Rosetta fragment library
///
///essentially a collection of SingleResidueTorsionFragmentLibrary (indexed by residue position)
///
class	TorsionFragmentLibrary : public utility::pointer::ReferenceCount{
/// ******************************************************************************************************
///  Code duplication alert:  TorsionFragment is going to be phased out
///   Please avoid writing any new code using these classes: look in core/fragment/ instead
///   look in /protocols/abinitio/FragmentMover  for usage examples
///
/// ******************************************************************************************************

public:
	/// default constructor
	TorsionFragmentLibrary() {}

	/// constructor with size (number of residue positions)
	TorsionFragmentLibrary( Size const size_in )
	{
		resize( size_in );
	}

	/// number of residues
	inline
	Size
	size() const {
		return fragments_.size();
	}

	/// overloaded [] operator to get framgnet pieces for a certain residue position
	SingleResidueTorsionFragmentLibrary &
	operator[]( Size const pos )
	{
		return *(fragments_[pos]);
	}

	/// overloaded [] operator to get framgnet pieces for a certain residue position (access only)
	SingleResidueTorsionFragmentLibrary const &
	operator[]( Size const pos ) const
	{
		return *(fragments_[pos]);
	}

	/// change the size of fragment library
	void
	resize( Size const new_size )
	{
		Size const old_size( fragments_.size() );
		fragments_.resize( new_size );
		for ( Size i=old_size+1; i<= new_size; ++i ) {
			fragments_[i] = new SingleResidueTorsionFragmentLibrary();
		}
	}

	// initialize fragment data from a classic Rosetta fragment library
	bool
	read_file( std::string const name, Size const frag_size, Size const nbb );

	// extract a fragment library with smaller fragment size from the one with larger lize
	bool
	derive_from_src_lib( Size const my_size, Size const src_size, TorsionFragmentLibraryCOP src_lib_op );

	/// @brief  Show some info -- right now just for debugging, ie not a full dump of the library
	void
	print( std::ostream & os ) const;

private:
	///a collection of SingleResidueTorsionFragmentLibrary owning pointers
	utility::vector1< SingleResidueTorsionFragmentLibraryOP > fragments_;

};


} // ns frags
} // ns protocols

#endif
