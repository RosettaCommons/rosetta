// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_src_devel_blab_classic_frags_TorsionFragment_HH
#define INCLUDED_src_devel_blab_classic_frags_TorsionFragment_HH

// Unit Headers
#include <protocols/frags/TorsionFragment.fwd.hh>

// Rosetta Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/moves/Mover.hh>

// ObjexxFCL Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>

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

	//using core::Real;
	//using core::Size;

///\brief a class for single piece of torsion fragment
///
///It stores torsion angles and secondary structure.
///Torsions are stored as vector of vector, such as TorsionFragment[frag_size][n_bb_torsion];
///SS are stored as vector, such as TorsionFragment[frag_length]
class TorsionFragment : public utility::pointer::ReferenceCount {

public:
	typedef core::Real Real;
	typedef core::Size Size;

public:

	/// constructor, fragment_size (3mer or 9mer) and number of backbone torsions(3 for protein)
	TorsionFragment(){}

	TorsionFragment( TorsionFragment const & src );

	TorsionFragmentOP
	clone() const;

	TorsionFragment( Size const size_in, Size const nbb_in )
	{
		set_size_and_nbb( size_in, nbb_in );
	}

	void
	set_size_and_nbb( Size const size, Size const nbb )
	{
		torsions_.resize( size );
		secstruct_.resize( size );
		sequence_.resize( size );
		for ( Size k=1; k<= size; ++k ) {
			torsions_[k].resize( nbb );
		}
	}
	virtual ~TorsionFragment();

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
	char
	get_secstruct( Size const pos ) const
	{
		return secstruct_[pos];
	}

	/// set seq for this position
	void
	set_sequence( Size const pos, char const setting )
	{
		sequence_[pos] = setting;
	}
	char
	get_sequence( Size const pos ) const
	{
		return sequence_[pos];
	}

	/// get value for specific torsion in this piece of fragment
	Real
	get_torsion( Size const pos, Size const tor ) const
	{
		return torsions_[pos][tor];
	}


private:
	/// torsion angles, the first dimension is fragment size, the second dimension is number of backbone torsions.
	utility::vector1< utility::vector1< Real > > torsions_;
	/// secstruct, dimensioned as fragment size
	utility::vector1< char > secstruct_;
	/// original fragment sequence, dimensioned as fragment size
	utility::vector1< char > sequence_;
};

std::ostream &
operator << ( std::ostream & out, TorsionFragment const & f );
std::istream &
operator >> ( std::istream & data, TorsionFragment & f );

/////////////////////////////////////////////////////////////////////////////
///\brief a class for collection of fragments for a single residue position
///
///essentially a vector of TorsionFragment (owning pointers)
class SingleResidueTorsionFragmentLibrary : public utility::pointer::ReferenceCount {
public:
	typedef core::Real Real;
	typedef core::Size Size;
/// ******************************************************************************************************
///  Code duplication alert:  TorsionFragment is going to be phased out
///   Please avoid writing any new code using these classes: look in core/fragment/ instead
///   look in /protocols/abinitio/FragmentMover  for usage examples
///
/// ******************************************************************************************************

public:
	virtual ~SingleResidueTorsionFragmentLibrary();

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

	void
	copy_fragments( SingleResidueTorsionFragmentLibrary const & src );

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

	/// erase a single fragment
	void
	erase( Size const index )
	{
		assert( index >= 1 && index <= size() );
		fragments_.erase( fragments_.begin() + (index-1) );
	}

	void
	clear()
	{
		fragments_.clear();
	}

private:
	/// a collection of TorsionFragments (by owning pointers)
	utility::vector1< TorsionFragmentOP > fragments_;
};

std::ostream &
operator << ( std::ostream & out, SingleResidueTorsionFragmentLibrary const & f );

std::istream &
operator >> ( std::istream & data, SingleResidueTorsionFragmentLibrary & lib );

/////////////////////////////////////////////////////////////////////////////
/// \brief a class for classic Rosetta fragment library
///
///essentially a collection of SingleResidueTorsionFragmentLibrary (indexed by residue position)
///
class	TorsionFragmentLibrary : public utility::pointer::ReferenceCount{
public:
	typedef core::Real Real;
	typedef core::Size Size;
/// ******************************************************************************************************
///  Code duplication alert:  TorsionFragment is going to be phased out
///   Please avoid writing any new code using these classes: look in core/fragment/ instead
///   look in /protocols/abinitio/FragmentMover  for usage examples
///
/// ******************************************************************************************************

public:
	/// default constructor
	TorsionFragmentLibrary() {}
	virtual ~TorsionFragmentLibrary();

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

	void
	delete_residue( Size const seqpos ); // invalidates all frags that overlap this position, shifts frags afterward

	void
	shift( int const current2desired_offset );

	void
	copy_fragments( TorsionFragmentLibrary const & src );

	// initialize fragment data from a classic Rosetta fragment library
	bool
	read_file( std::string const name, Size const frag_size, Size const nbb );

	// extract a fragment library with smaller fragment size from the one with larger lize
	bool
	derive_from_src_lib( Size const my_size, Size const src_size, TorsionFragmentLibraryCOP src_lib_op );

	bool
	derive_from_src_lib( Size const my_size, Size const src_size, TorsionFragmentLibrary const & src_lib );

	/// @brief  Show some info -- right now just for debugging, ie not a full dump of the library
	void
	print( std::ostream & os ) const;

	void
	clear()
	{
		fragments_.clear();
	}

private:
	///a collection of SingleResidueTorsionFragmentLibrary owning pointers
	utility::vector1< SingleResidueTorsionFragmentLibraryOP > fragments_;

};

std::ostream &
operator << ( std::ostream & out, TorsionFragmentLibrary const & f );

std::istream &
operator >> ( std::istream & data, TorsionFragmentLibrary & lib );

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FragLib : public utility::pointer::ReferenceCount {
public:
	typedef core::Real Real;
	typedef core::Size Size;

public:

	typedef std::map< Size, TorsionFragmentLibraryOP > FragMap;

public:

	FragLib( FragMap const & frag_map ):
		frag_map_( frag_map )
	{}

	FragLib() {}

	FragLib( FragLib const & src ): ReferenceCount() { copy_fragments( src ); }

	TorsionFragmentLibrary const &
	library( Size const size ) const;

	TorsionFragmentLibrary &
	library( Size const size );

	/// slow
	utility::vector1< Size >
	frag_sizes() const;

	void
	clear()
	{
		frag_map_.clear();
	}

	void
	delete_residue( Size const seqpos ); // invalidates all frags that overlap this position, shifts frags afterward

	void
	shift( int const current2desired_offset );

	void
	copy_fragments( FragLib const & src );

private:
	FragMap frag_map_;

};

std::ostream &
operator << ( std::ostream & out, FragLib const & f );

std::istream &
operator >> ( std::istream & data, FragLib & lib );

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern utility::vector1< core::Size > const empty_size_list;
	//extern std::string const empty_string;

void
add_vall_fragments(
									 utility::vector1< core::Size > const & frag_sizes,
									 core::Size const nfrags,
									 core::pose::Pose const & pose,
									 core::kinematics::MoveMap const & mm,
									 std::string const & secstruct,
									 core::Real const seq_weight,
									 core::Real const ss_weight,
									 FragLib & frag_lib,
									 utility::vector1< core::Size > const & homs_to_exclude = empty_size_list,
									 core::Real const bb_weight = 0.0,
									 std::string const & bigbins = std::string(),
									 std::string const & inputseq = std::string()
									 );
void
add_vall_fragments(
									 utility::vector1< core::Size > const & frag_sizes,
									 core::Size const nfrags,
									 core::pose::Pose const & pose,
									 core::kinematics::MoveMap const & mm,
									 utility::vector1< std::map< char, core::Real > > const & target_ss, //// THIS IS THE DIFFERENCE
									 core::Real const seq_weight,
									 core::Real const ss_weight,
									 FragLib & frag_lib,
									 utility::vector1< core::Size > const & homs_to_exclude = empty_size_list,
									 core::Real const bb_weight = 0.0,
									 std::string const & bigbins = std::string(),
									 std::string const & inputseq = std::string()
									 );

FragLibOP
setup_vall_fragments(
										 utility::vector1< core::Size > const & frag_sizes,
										 core::Size const nfrags,
										 core::pose::Pose const & pose,
										 core::kinematics::MoveMap const & mm,
										 std::string const & secstruct,
										 core::Real const seq_weight,
										 core::Real const ss_weight,
										 utility::vector1< core::Size > const & homs_to_exclude = empty_size_list
										 );

void
add_vall_cheating_fragments(
														utility::vector1< core::Size > const & frag_sizes,
														core::Size const nfrags,
														core::pose::Pose const & pose,
														core::kinematics::MoveMap const & mm,
														std::string const & secstruct,
														core::Real const seq_weight,
														core::Real const ss_weight,
														core::Real const torsion_weight,
														core::Real const min_torsion_dev,
														core::Real const max_torsion_dev,
														FragLib & frag_lib,
														utility::vector1< core::Size > const & homs_to_exclude = empty_size_list
														);

FragLibOP
setup_vall_cheating_fragments(
															utility::vector1< core::Size > const & frag_sizes,
															core::Size const nfrags,
															core::pose::Pose const & pose,
															core::kinematics::MoveMap const & mm,
															std::string const & secstruct,
															core::Real const seq_weight,
															core::Real const ss_weight,
															core::Real const torsion_weight,
															core::Real const min_torsion_dev,
															core::Real const max_torsion_dev,
															utility::vector1< core::Size > const & homs_to_exclude = empty_size_list
															);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
fill_in_gaps(
						 core::Size const nfrags,
						 core::pose::Pose const & pose,
						 std::string const & secstruct,
						 core::Real const seq_weight,
						 core::Real const ss_weight,
						 FragLib & frag_lib,
						 utility::vector1< core::Size > const & homs_to_exclude = empty_size_list,
						 bool const allow_uninitialized_secstruct = false
						 );

bool
ss_length_check(
								core::Size const min_len_helix,
								core::Size const min_len_strand,
								core::pose::Pose const & pose
								);

class TorsionFragmentMover : public protocols::moves::Mover {
public:
	typedef core::Real Real;
	typedef core::Size Size;

public:
	TorsionFragmentMover( FragLibCOP fraglib, core::kinematics::MoveMapCOP mm ):
		protocols::moves::Mover( "TorsionFragment" ),
		fraglib_( fraglib ),
		mm_( mm ),
		frag_size_( 0 ),
		check_ss_lengths_( false ),
		min_len_helix_( 1 ),
		min_len_strand_( 1 ),
		last_inserted_frag_begin_( 0 )
	{}

	std::string
	get_name() const { return "TorsionFragmentMover"; }

	void
	check_ss_lengths( bool const setting )
	{
		check_ss_lengths_ = setting;
	}

	void
	min_len_strand( Size const setting )
	{
		min_len_strand_ = setting;
	}

	void
	min_len_helix( Size const setting )
	{
		min_len_helix_ = setting;
	}

	Size
	frag_size() const { return frag_size_; }

	Size
	last_inserted_frag_begin() const { return last_inserted_frag_begin_; }

	void
	frag_size( Size const setting )
	{
		frag_size_ = setting;
	}

	bool
	pose_passes_ss_length_check( core::pose::Pose const & pose ) const;

	void
	apply( core::pose::Pose & pose );


private:
	FragLibCOP fraglib_;
	core::kinematics::MoveMapCOP mm_;
	Size frag_size_;
	bool check_ss_lengths_;
	Size min_len_helix_;
	Size min_len_strand_;
	Size last_inserted_frag_begin_;

};


void
insert_fragment(
								int const begin,
								int const end,
								core::pose::Pose & pose,
								TorsionFragmentLibrary const & lib,
								int const desired_insert_pos = 0
								);


void
insert_random_fragments_in_flexible_protein_regions(
																										utility::vector1< core::Size > const & flex_protein,
																										FragLib const & frag_lib,
																										core::pose::Pose & pose
																										);

} // ns frags
} // ns protocols
#endif
