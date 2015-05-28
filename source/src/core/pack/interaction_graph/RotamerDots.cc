// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/RotamerDots.cc
/// @brief  RotamerDots classes files - ported from rosetta++
/// @author Andrew Leaver-Fay
/// @author Ron Jacak

// Unit Headers
#include <core/pack/interaction_graph/RotamerDots.hh>

// Project headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Atom.hh>
//#include <core/scoring/sasa.hh>
#include <core/scoring/sasa/util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>


#include <core/conformation/Residue.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/ubyte.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1.hh>

// Numeric Headers
#include <numeric/constants.hh> // pi
#include <numeric/xyzVector.hh> // to get distance

// Utility Headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/options/BooleanVectorOption.hh>

// C++ Headers
#include <vector>
#include <iostream>

#include <ObjexxFCL/FArray2D.hh>

// Boost headers
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

static thread_local basic::Tracer TR_DS( "core.pack.interaction_graph.RotamerDots.DotSphere" );
static thread_local basic::Tracer TR_RD( "core.pack.interaction_graph.RotamerDots.RotamerDots" );
static thread_local basic::Tracer TR_RDC( "core.pack.interaction_graph.RotamerDots.RotamerDotsCache" );
static thread_local basic::Tracer TR_RDRD( "core.pack.interaction_graph.RotamerDots.RotamerDotsRadiusData" );


using namespace ObjexxFCL::format;

// Singleton instance and mutex static data members
namespace utility {

using core::pack::interaction_graph::RotamerDotsRadiusData;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< RotamerDotsRadiusData >::singleton_mutex_{};
template <> std::atomic< RotamerDotsRadiusData * > utility::SingletonBase< RotamerDotsRadiusData >::instance_( 0 );
#else
template <> RotamerDotsRadiusData * utility::SingletonBase< RotamerDotsRadiusData >::instance_( 0 );
#endif

}

namespace core {
namespace pack {
namespace interaction_graph {

bool unpack_ubyte( ObjexxFCL::ubyte const & value, core::Size which_bit ) {

	switch ( which_bit ) {
		case 0 :
			return value & static_cast< ObjexxFCL::ubyte >(0x01);
		case 1 :
			return value & static_cast< ObjexxFCL::ubyte >(0x02);
		case 2 :
			return value & static_cast< ObjexxFCL::ubyte >(0x04);
		case 3 :
			return value & static_cast< ObjexxFCL::ubyte >(0x08);
		case 4 :
			return value & static_cast< ObjexxFCL::ubyte >(0x10);
		case 5 :
			return value & static_cast< ObjexxFCL::ubyte >(0x20);
		case 6 :
			return value & static_cast< ObjexxFCL::ubyte >(0x40);
		case 7 :
			return value & static_cast< ObjexxFCL::ubyte >(0x80);
		default :
			utility_exit(); // this should never happen
	}

	return false;
}

void write_sphere_list_header( std::ostream & ostr, std::string const & color, bool off = false ) {

	ostr << "@spherelist color= " << color << " radius= 0.1";
	if ( off ) {
		ostr << " off";
	}
	ostr << "\n";
}

void write_dot( std::ostream & ostr, core::Vector const & center, core::Real radius, Size const dot_index, std::string const & dot_name ) {
	Vector coord = center + radius * RotamerDots::dot_coord(  dot_index );
	ostr << "{" << dot_name << "} P " << coord.x() << " " << coord.y() << " " << coord.z() << "\n";
}


void
write_sphere_list_uv1(
	std::ostream & ostr,
	std::string const & label,
	std::string const & color,
	core::Vector const & center,
	core::Real radius,
	utility::vector1< ObjexxFCL::ubyte > const & dot_masks
)
{
	bool first = true;
	Size count = 1;
	std::string name = label;
	for ( Size ii = 1; ii <= 21; ++ii ) {
		for ( Size jj = 0; jj < 8; ++jj ) {
			if ( unpack_ubyte( dot_masks[ ii ], jj ) ) {
				if ( first ) {
					first = false;
					write_sphere_list_header( ostr, color );
					write_dot( ostr, center, radius, count, name );
					name = "\"";
				} else {
					write_dot( ostr, center, radius, count, name );
				}
			}
			++count;
			if ( count == 163 ) break;
		}
		if ( count == 163 ) break;
	}
}

void write_sphere_list_farray( std::ostream & ostr, std::string const & label, std::string const & color,
	core::Vector const & center, core::Real radius, ObjexxFCL::FArray1< ObjexxFCL::ubyte > const & dot_masks ) {

	bool first = true;
	Size count = 1;
	std::string name = label;
	for ( Size ii = 1; ii <= 21; ++ii ) {
		for ( Size jj = 0; jj < 8; ++jj ) {
			if ( unpack_ubyte( dot_masks( ii ), jj ) ) {
				if ( first ) {
					first = false;
					write_sphere_list_header( ostr, color );
					write_dot( ostr, center, radius, count, name );
					name = "\"";
				} else {
					write_dot( ostr, center, radius, count, name );
				}
			}
			++count;
			if ( count == 163 ) break;
		}
		if ( count == 163 ) break;
	}
}


//----------------------------------------------------------------------------//
//----------------------------- Dot Sphere Class -----------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// default constructor, initializes all dot counts to zero
///
DotSphere::DotSphere() :
	dots_coverage_count_(),
	num_covered_(0),
	num_covered_current_( false )
{
	zero();
}

///
DotSphere::~DotSphere() {}

///
/// @brief
/// copy constructor
///
/// @details
/// memcpy is much faster than the FArray operator =
///
DotSphere::DotSphere( DotSphere const & rhs ) {
	num_covered_ = rhs.num_covered_;
	num_covered_current_= rhs.num_covered_current_;
	memcpy( dots_coverage_count_, rhs.dots_coverage_count_, NUM_COUNTS_TO_ALLOCATE );
}

///
/// @brie=
/// assignment operator
///
DotSphere const & DotSphere::operator= ( DotSphere const & rhs ) {
	num_covered_ = rhs.num_covered_;
	num_covered_current_ = rhs.num_covered_current_;
	memcpy( dots_coverage_count_, rhs.dots_coverage_count_, NUM_COUNTS_TO_ALLOCATE );
	return *this;
}

///
/// @brief
/// Comparison operator. Using this in debugging. When alternate state rotamer dots is not equal to current state
/// rotamer dots, then I print some extra debugging information. But this could be useful for other purposes, too.
/// Since RotamerDots objects contain DotSphere objects, to compare RD objects, this class needs its own comparison
/// operator, too.
///
bool DotSphere::operator!= ( DotSphere const & rhs ) {

	if ( num_covered_ != rhs.num_covered_ || num_covered_current_ != rhs.num_covered_current_ )
		return true;
	for ( Size ii=0; ii < NUM_BYTES_IN_DOT_SPHERE_OVERLAP_ARRAYS; ++ii ) {
		for ( Size jj=0; jj < 8; ++jj ) {
			Size ii8_jj = ii*8 + jj;
			if ( dots_coverage_count_[ ii8_jj ] != rhs.dots_coverage_count_[ ii8_jj ] ) {
				return true;
			}
		}
	}

	return false;
}

///
/// @details
/// sets the dot coverage counts to zero for all dots
/// memset is fast -- a lot of time is spent in this function so I'm using c-style arrays instead of the FArrays
///
void DotSphere::zero() {
	if ( num_covered_current_ && num_covered_ == 0 )
		return;
	memset( dots_coverage_count_, 0, NUM_COUNTS_TO_ALLOCATE );
	num_covered_ = 0;
	num_covered_current_ = true;
}

///
/// @brief
/// increment the count for the dots using an input ubyte array.
///
/// @details
/// Each bit in this ubyte array corresponds to a single dot on the surface of this sphere. Dot coverage counts are
/// incremented for dots whose corresponding bits in the ubyte array are '1'. dots_coverage_count_ is C-style array
/// of unsigned chars. So, index 0 returns the first char (or one byte), index 20*8+7, or 167, is the last char/byte.
///
/// @param
/// overlap_mask - a utility::vector1 of size 21 that holds ObjexxFCL ubytes. overlap_mask[ bb ] will return a single ubyte
/// which will determine whether the dots for that region of the vector should be incremented. Because it's a vector1
/// and this method uses a 0-based array, we have to remember to convert the array index to 1-based before looking at
/// what's in overlap_mask.
///
void DotSphere::increment_count( utility::vector1< ObjexxFCL::ubyte > const & overlap_mask ) {

	// for ( j = 1; val && j <= 0x80; j <<= 1 )
	// In this example, j is going to be 1, then 2, 4, 8, 16, 32, 64 and finally 128 as the loop progresses.
	// In binary, that's 0000:0001, 0000:0010, 0000:0100, 0000:1000, 0001:0000, 0010:0000, 0100:0000 and 1000:0000.
	// Since there's no way to specify binary constants in C++, it's clearer to use Hex: 0x01, 0x02, 0x04, 0x08, 0x10,
	// 0x20, 0x40, and 0x80.
	// So bitwise AND with 0x04, is bitwise AND with 0000:0100 which returns either 0000:0100 or 0000:0000. We case that
	// to a bool (which is a byte), which gives 0000:0001 or 0000:0000, which we then cast back to an unsigned char
	// because that's what the type is in dots_coverage_count_.  The reason we do this is because we want to increment
	// the count at the bit location of dots_coverage_count regardless if the value is 0010:0000 or 0000:0001. If we tried
	// to add 0000:1000 to the dots_coverage_count bit location, we'd be adding '8' and not '1' to that location.
	// This is just one way to do this; there are definitely other ways.

	for ( Size bb = 0; bb < NUM_BYTES_IN_DOT_SPHERE_OVERLAP_ARRAYS; ++bb ) {
		const Size bb8 = bb*8;
		dots_coverage_count_[ bb8 + 0 ] += static_cast <unsigned char> ( static_cast<bool> (overlap_mask[ bb+1 ] & 0x01));
		dots_coverage_count_[ bb8 + 1 ] += static_cast <unsigned char> ( static_cast<bool> (overlap_mask[ bb+1 ] & 0x02));
		dots_coverage_count_[ bb8 + 2 ] += static_cast <unsigned char> ( static_cast<bool> (overlap_mask[ bb+1 ] & 0x04));
		dots_coverage_count_[ bb8 + 3 ] += static_cast <unsigned char> ( static_cast<bool> (overlap_mask[ bb+1 ] & 0x08));
		dots_coverage_count_[ bb8 + 4 ] += static_cast <unsigned char> ( static_cast<bool> (overlap_mask[ bb+1 ] & 0x10));
		dots_coverage_count_[ bb8 + 5 ] += static_cast <unsigned char> ( static_cast<bool> (overlap_mask[ bb+1 ] & 0x20));
		dots_coverage_count_[ bb8 + 6 ] += static_cast <unsigned char> ( static_cast<bool> (overlap_mask[ bb+1 ] & 0x40));
		dots_coverage_count_[ bb8 + 7 ] += static_cast <unsigned char> ( static_cast<bool> (overlap_mask[ bb+1 ] & 0x80));
	}

	//std::cout << "dots_coverage_count_: ";
	//print( std::cout );

	num_covered_current_ = false;
}

///
/// @brief
/// returns the total number of dots on this atom whose coverage count is 0
///
/// @details
/// if the coverage count has not been modified since the last time the number of covered dots was counted, then the
/// method uses the cached result.
///
Size DotSphere::get_num_uncovered() const {
	if ( ! num_covered_current_ ) {
		count_num_covered();
	}
	return NUM_DOTS_TOTAL - num_covered_;
}

///
/// @brief
/// returns the total number of dots on this atom with a non-zero coverage count
///
/// @details
/// if the coverage count has not been modified since the last time the number of covered dots was counted, then the
/// method uses the cached result
///
Size DotSphere::get_num_covered() const {
	if ( ! num_covered_current_ )
		count_num_covered();
	return num_covered_;
}

///
/// @brief
/// iterates across all dots and stores the number with a non-zero coverage count for later use
///
/// @details
/// both num_covered_ and num_covered_current_ are declared mutable so that they may be modified in this const method
///
void DotSphere::count_num_covered() const {

	num_covered_ = 0;
	//TR_DS << "count_num_covered(): counts not current. coverage_counts_: " << std::endl;
	for ( Size ii=0; ii < NUM_DOTS_TOTAL; ++ii ) {
		//if ( ii % 8 == 0 ) TR_DS << " ";
		if ( dots_coverage_count_[ ii ] != 0 ) {
			num_covered_++;
			//TR_DS << "1";
		}
	}
	//TR_DS << std::endl;
	num_covered_current_ = true;
}

///
/// @brief
/// decrements the coverage count for this sphere by the coverage count of the rhs sphere
///
DotSphere const & DotSphere::operator -= ( DotSphere const & rhs ) {
	num_covered_current_ = false;
	for ( Size ii = 0; ii < NUM_COUNTS_TO_ALLOCATE; ++ii ) {
		dots_coverage_count_[ ii ] -= rhs.dots_coverage_count_[ ii ];
	}
	return *this;
}

///
/// @brief
/// increments the coverage count for this sphere by the coverage count of the rhs sphere
///
DotSphere const & DotSphere::operator += ( DotSphere const & rhs ) {
	num_covered_current_ = false;
	for (Size ii = 0; ii < NUM_COUNTS_TO_ALLOCATE; ++ii) {
		dots_coverage_count_[ ii ] += rhs.dots_coverage_count_[ ii ];
	}
	return *this;
}

///
/// @brief
/// Returns a boolean indicating whether the given dot is covered. Note, this function takes in a 1-based dot-index
/// and converts that to 0-based for the C-style array used by this class.
///
bool DotSphere::get_dot_covered( Size dot_index ) const {
debug_assert( dot_index > 0 && dot_index <= NUM_DOTS_TOTAL );
	return ( dots_coverage_count_[ dot_index - 1 ] != 0 );
}

///
/// @brief
/// note, this method results in loss of information; counts > 1 are truncated to 1.
/// bitwise OR with 0000:0001 results in 0000:0001.
///
void DotSphere::write_to_compact_array( utility::vector1< ObjexxFCL::ubyte > & compact ) const {

	for ( Size ii=0; ii < NUM_BYTES_IN_DOT_SPHERE_OVERLAP_ARRAYS; ++ii ) {
		Size bb8 = ii*8;
		if ( dots_coverage_count_[ bb8 + 0 ] > 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x01);
		if ( dots_coverage_count_[ bb8 + 1 ] > 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x02);
		if ( dots_coverage_count_[ bb8 + 2 ] > 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x04);
		if ( dots_coverage_count_[ bb8 + 3 ] > 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x08);
		if ( dots_coverage_count_[ bb8 + 4 ] > 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x10);
		if ( dots_coverage_count_[ bb8 + 5 ] > 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x20);
		if ( dots_coverage_count_[ bb8 + 6 ] > 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x40);
		if ( dots_coverage_count_[ bb8 + 7 ] > 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x80);
	}
}

void DotSphere::invert_to_compact_array( utility::vector1< ObjexxFCL::ubyte > & compact ) const
{
	/// convention; the dots that do not exist are not to be considered exposed.
	/// the first two bytes in the 21st position describe dots 161 and 162; the remaining
	/// dots describe 6 dots that do not exist... even if their coverage counts are 0, do not
	/// report that they are exposed!
	std::fill( compact.begin(), compact.end(), static_cast< ObjexxFCL::ubyte > (0) );
	for ( Size ii=0; ii < NUM_BYTES_IN_DOT_SPHERE_OVERLAP_ARRAYS - 1; ++ii ) {
		Size bb8 = ii*8;
		if ( dots_coverage_count_[ bb8 + 0 ] == 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x01);
		if ( dots_coverage_count_[ bb8 + 1 ] == 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x02);
		if ( dots_coverage_count_[ bb8 + 2 ] == 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x04);
		if ( dots_coverage_count_[ bb8 + 3 ] == 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x08);
		if ( dots_coverage_count_[ bb8 + 4 ] == 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x10);
		if ( dots_coverage_count_[ bb8 + 5 ] == 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x20);
		if ( dots_coverage_count_[ bb8 + 6 ] == 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x40);
		if ( dots_coverage_count_[ bb8 + 7 ] == 0 ) compact[ ii+1 ] |= static_cast< ObjexxFCL::ubyte > (0x80);
	}
	Size const last = NUM_BYTES_IN_DOT_SPHERE_OVERLAP_ARRAYS - 1;
	Size const lastx8 = 8 * last;
	if ( dots_coverage_count_[ lastx8 + 0 ] == 0 ) compact[ last+1 ] |= static_cast< ObjexxFCL::ubyte > (0x01);
	if ( dots_coverage_count_[ lastx8 + 1 ] == 0 ) compact[ last+1 ] |= static_cast< ObjexxFCL::ubyte > (0x02);

}

///
/// @brief
/// Writes coverage counts to the output stream.  if a dot is covered by 10 or more residues, prints 9 to the output stream instead.
/// Useful for debugging.
///
void DotSphere::print( std::ostream & os ) const {

	for ( Size ii=0; ii < NUM_COUNTS_TO_ALLOCATE; ++ii ) {
		if ( ii % 16 == 0 )
			os << ii+1 << ":";
		if ( dots_coverage_count_[ii] < 10 ) {
			os << (Size)dots_coverage_count_[ii];
		} else os << "9";

		if (ii % 8 == 7)
			os << " ";
	}
	os << std::endl;
	return;
}

///
/// @brief
/// invokes print on the input DotSphere object
///
std::ostream & operator<< ( std::ostream & os, DotSphere const & ds ) {
	ds.print(os);
	return os;
}


//----------------------------------------------------------------------------//
//---------------------------- Rotamer Dots Class ----------------------------//
//----------------------------------------------------------------------------//

Size const RotamerDots::num_bytes_ = 21;
Real RotamerDots::probe_radius_ = 1.4;

bool RotamerDots::sasa_arrays_initialized_ = false;
utility::vector1< core::Vector > RotamerDots::dot_coords_( 0 );

ObjexxFCL::FArray2D_int const *   RotamerDots::lg_angles_( 0 );
ObjexxFCL::FArray2D_ubyte const * RotamerDots::lg_masks_( 0 );

///
RotamerDots::RotamerDots():
	rotamer_(/* 0 */),
	num_atoms_(0),
	sasa_(0),
	sasa_is_current_(false)
{}

///
/// @brief
/// Custom constructor for a RotamerDots object
///
/// @details
/// One RotamerDots object get allocated for every state of a first class IG Node, for all first class IG Nodes of a
/// protein being designed. That's potentially a very large number of states. This class should only hold the information
/// it needs to hold to do its job.
///
RotamerDots::RotamerDots(
	conformation::ResidueCOP rotamer,
	bool exclude_hydrogen_atoms /*=false*/,
	bool use_expanded_polar_atom_radii /*=false*/
) :
	rotamer_(rotamer),
	sasa_(0),
	sasa_is_current_(false)
{
	if ( exclude_hydrogen_atoms ) {
		num_atoms_ = rotamer->nheavyatoms();
	} else {
		num_atoms_ = rotamer->natoms();
	}

	atom_counts_.resize( num_atoms_ );
	atom_sasa_.resize( num_atoms_ );

	if ( ! sasa_arrays_initialized_ ) {
		initialize_sasa_arrays();
		//that function will set sasa_arrays_initialized_ to true;
	}

	// every instance needs to have the radii_ pointer set. radii_ can't be static because two RotamerDots objects may
	// want to use different atom radii for calculating SASA.
	if ( use_expanded_polar_atom_radii ) {
		radii_ = RotamerDotsRadiusData::get_instance()->get_NACCESS_SASA_radii_with_expanded_polars();
	} else {
		radii_ = RotamerDotsRadiusData::get_instance()->get_NACCESS_SASA_radii();
	}

}

///
RotamerDots::~RotamerDots() {
	//TR_RD << "called destructor" << std::endl;
}

///
/// @brief
/// copy constructor
///
RotamerDots::RotamerDots( RotamerDots const & rhs ) :
	utility::pointer::ReferenceCount(),
	rotamer_( rhs.rotamer_ ),
	num_atoms_( rhs.num_atoms_ ),
	atom_counts_( rhs.atom_counts_ ),

	sasa_( rhs.sasa_ ),
	sasa_is_current_( rhs.sasa_is_current_ ),
	atom_sasa_( rhs.atom_sasa_ ),
	radii_( rhs.radii_ )
{}

///
/// @brief
/// Copy method for the RotamerDots class. Also used by the assignment operator.
///
void RotamerDots::copy( RotamerDots const & rhs ) {
	rotamer_ = rhs.rotamer_;
	num_atoms_ = rhs.num_atoms_;
	atom_counts_ = rhs.atom_counts_;

	sasa_ = rhs.sasa_;
	atom_sasa_ = rhs.atom_sasa_;
	sasa_is_current_ = rhs.sasa_is_current_;
	radii_ = rhs.radii_;
}

///
RotamerDots const & RotamerDots::operator=( RotamerDots const & rhs ) {
	copy( rhs );
	return *this;
}

///
/// @brief
/// Used during debugging of the HPatchIG.  Some extra information is printed if current state dots is NOT EQUAL to
/// alternate state dots at a Node/BGNode.
///
bool RotamerDots::operator!=( RotamerDots const & rhs ) {

	if ( state_unassigned() || rhs.state_unassigned() ) {
		return true; // I guess they could both be unassigned, but better to return not equal than equal
	}

	if ( rotamer_ != rhs.rotamer_ || num_atoms_ != rhs.num_atoms_ ) return true;
	if ( sasa_ != rhs.sasa_ || sasa_is_current_ != rhs.sasa_is_current_ ) return true;
	if ( atom_counts_.size() != rhs.atom_counts_.size() ) return true;
	if ( atom_sasa_.size() != rhs.atom_sasa_.size() ) return true;
	if ( radii_ != rhs.radii_ ) return true;

	//TR_RD << "operator!=() score info same. checking if dot sphere objects match." << std::endl;
	for ( Size ii=1; ii <= atom_counts_.size(); ++ii ) {
		if ( atom_counts_[ ii ] != rhs.atom_counts_[ ii ] )
			return true;
	}

	return false;
}

///
/// @brief
/// Zeros out all of the contained data except the rotamer pointer and the radii array.
///
/// @details
/// So far, this function only gets called by the BGNode::prep_for_simA() call so that multiple runs through an
/// interaction graph can be done. If the rotamer dots object on the BGNodes isn't "cleared" after a run, then the run
/// immediately following will have the incorrect counts.
///
void RotamerDots::zero() {

	// leave the rotamer and num_atoms_ variables untouched
	//rotamer_ = 0;
	//num_atoms_ = 0;

	sasa_ = 0.0;
	sasa_is_current_ = false;

	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
		atom_sasa_[ ii ] = 0.0;
		atom_counts_[ ii ].zero(); // calls zero() on each DotSphere instance
	}

}

///
/// @brief
/// Returns true if this RotamerDots object has any sphere overlap with the passed in RotamerDots object.
///
/// @details
/// This method only checks to see if two RotamerDots objects are within touching distance of each other. It is used
/// to determine whether Edges or BGEdges should be created in the IG. Calculate this using the expanded polar atom
/// radii. If we don't, there's a chance that a state substitution on a Node may cause SASA changes (when expanded polars
/// are used) on a BGNode, but if we didn't assume expanded radii in this method, there would be no edge between the two
/// nodes.
///
bool RotamerDots::overlaps( RotamerDots const & other ) const {

	Real distance_squared, atom1_radius, atom2_radius;

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		for ( Size jj = 1; jj <= other.get_num_atoms(); ++jj ) {

			distance_squared = get_atom_coords_xyz( ii ).distance_squared( other.get_atom_coords_xyz( jj ) );

			atom1_radius = get_atom_radius( ii ) + probe_radius_;
			atom2_radius = other.get_atom_radius( jj ) + probe_radius_;

			// exit if large probe radii touch
			if ( distance_squared <= (atom1_radius + atom2_radius) * (atom1_radius + atom2_radius) )
				return true;
		}
	}
	return false;
}

///
core::conformation::ResidueCOP
RotamerDots::rotamer() const {
	return rotamer_;
}

///
/// @brief
/// Is the state of this RotamerDots object unassigned?
///
bool RotamerDots::state_unassigned() const {
	if ( rotamer_ == 0 )
		return true;
	return false;
}

///
/// @brief
/// Returns the number of atoms this RotamerDots object is keeping SASA for.
///
Size RotamerDots::get_num_atoms() const {
	return num_atoms_;
}

///
/// @brief
/// Return the xyz coordinates of an atom in this RotamerDots instance.
///
numeric::xyzVector< Real > RotamerDots::get_atom_coords_xyz( Size atom_index ) const {
	if ( rotamer_ == 0 )
		return numeric::xyzVector< Real >(0,0,0);

	return rotamer_->xyz( atom_index );
}

///
/// @brief
/// Returns the SASA radius for the passed in atom type. The DB file should have been read in at construct time.
///
/// @details
/// Many of the functions in this class iterate over 1 .. num_atoms_.
/// That's not the same thing as an atom type index which is what the radii vector is indexed with. So before we can return
/// the radius, we have to convert the passed in atom_index into the right atom in the residue and then use that to get the
/// right type.
///
Real RotamerDots::get_atom_radius( Size atom_index ) const {
	if ( rotamer_ == 0 )
		return 0.0;

	conformation::Atom const & atom( rotamer_->atom( atom_index ) );
	return (*radii_)[ atom.type() ];

}

///
/// @brief
/// Same as the above, but skips the conversion from atom index to atom type index.
///
core::Real
RotamerDots::radius_for_attype( Size const attype_index ) {
	return (*radii_)[ attype_index ];
}

///
/// @brief
/// Returns the maximum atom radius. Used only by the SurfacePotential class.
///
core::Real
RotamerDots::max_atom_radius() {
	return utility::max( *radii_ );
}

///
/// @brief
/// Returns a pointer to the radii vector. Used only by the InvRotamerDots class.
///
utility::vector1< Real > const *
RotamerDots::get_radii() const {
	return radii_;
}

///
/// @brief
/// Inverts the current dot counts and saves them to the passed in vector.
///
void RotamerDots::invert_to_boolmasks( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & inv_dots ) const {

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		atom_counts_[ ii ].invert_to_compact_array( inv_dots[ ii ] );
	}
}

/// @brief invert the current dot counts for a subset of the atoms in this rotamer.
void
RotamerDots::invert_to_boolmasks(
	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & inv_dots,
	utility::vector1< Size > const & ats_to_update
) const
{
	for ( Size ii = 1, iiend = ats_to_update.size(); ii <= iiend; ++ii ) {
		atom_counts_[ ats_to_update[ ii ] ].invert_to_compact_array( inv_dots[ ats_to_update[ ii ] ] );
	}
}


core::Vector
RotamerDots::dot_coord( Size index ) {
	if ( !sasa_arrays_initialized_ )
		initialize_sasa_arrays();
	return dot_coords_[ index ];
}

///
/// @brief
/// Initializes the pointers to the angles and masks FArrays used by sasa.cc and inits the dot sphere coordinates.
///
/// @details
/// This call should only occur once (when the first RotamerDots object get constructed) and never again.
///
void RotamerDots::initialize_sasa_arrays() {

	if ( sasa_arrays_initialized_ ) return;
	sasa_arrays_initialized_ = true;

	lg_angles_ = ( & core::scoring::sasa::get_legrand_sasa_angles() );
	lg_masks_  = ( & core::scoring::sasa::get_legrand_sasa_masks()  );

	initialize_dot_coords( dot_coords_ );

	return;
}

///
/// @brief
/// computes and stores self-induced dot coverage. uses a vector1 of vector1 of ubytes to store the calculated overlap information.
///
/// @details
/// uses get_atom_atom_coverage() which in turn uses get_overlap() and get_orientation() method calls in sasa.cc to get
/// the right overlap "masks". when all atom pairs are complete, converts the masks into coverage "counts" which are stored
/// in the DotSphere class (member variable atom_counts_).
///
/// to handle the possibility of keeping expanded_polar SASA, we have to run the nested loop over atoms a second time through
/// using the expanded_polar version of get_atom_radius().
///
void RotamerDots::increment_self_overlap() {

	using namespace utility; // for utility::vector1

	vector1< vector1< ObjexxFCL::ubyte > > self_overlap( num_atoms_, vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		Vector atom1_xyz = get_atom_coords_xyz( ii );
		Real atom1_radius = get_atom_radius( ii );

		// only have to iterate over the higher indexed atoms for INTRA res overlaps
		for ( Size jj = ii+1; jj <= num_atoms_; ++jj ) {
			Vector atom2_xyz = get_atom_coords_xyz( jj );
			Real atom2_radius = get_atom_radius( jj );

			Real dist = 0.0f;
			get_atom_atom_coverage( atom1_xyz, atom1_radius, atom2_xyz, atom2_radius, self_overlap[ii], self_overlap[jj], dist );
		}
	}

	//TR_RD << "increment_self_overlap(): incrementing with counts: " << std::endl;
	//for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
	//	RotamerDotsCache rdc;
	//	rdc.print_bit_string( self_overlap[ ii ] );
	//}

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		//TR_RD << "increment_self_overlap(): calling increment_count() on atom dotsphere " << ii << std::endl;
		atom_counts_[ii].increment_count( self_overlap[ ii ] );
	}

	sasa_is_current_ = false;
}

///
/// @brief
/// Add rotamer coverage counts for dots on this object only, leaving rhs unaltered.
///
/// @details
/// In the context of the HPatchIG, this method is called by all FCNodes to increment the overlap a BG residue has on
/// the FCNode. It is called by all BG Edges, to make sure that all FCNodes that are connected to a BG residue get this
/// method called by them. 'other' in this case is a BG residue, and 'this_overlap_on_other' is the overlap that is
/// caused by all-states-possible-at-this-node on the BG node. Yes, that's keeping the same information in two places,
/// but it makes updating hpatch score calculations later faster. (ronj)
///
void RotamerDots::increment_this_and_cache(
	RotamerDots const & other,
	RotamerDotsCache & this_overlap_on_other,
	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
) {

	RotamerDotsCache this_dots_covered_by_other( num_atoms_ );
	// can't make the variable above be static because this function gets called by all kinds of FC nodes, which will
	// have varying numbers of atoms.  Then the member vector inside RDC would not be sized correctly and errors would
	// occur.

	/*TR_RD << "atom_atom_overlaps_cache.size(): " << atom_atom_overlaps_cache.size() << std::endl;
	for ( Size ii=1; ii <= atom_atom_overlaps_cache.size(); ++ii ) {
		TR_RD << "atom_atom_overlaps_cache[ " << ii << " ].size(): " << atom_atom_overlaps_cache[ ii ].size() << std::endl;
	}
	TR_RD << std::endl;*/

	get_overlap_cache( other, this_overlap_on_other, this_dots_covered_by_other, atom_atom_overlaps_cache );

	// 'increment_count_for_some' is a class method in RotamerDotsCache objects
	//TR_RD << "increment_this_and_cache(): this_dots_covered_by_other: " << std::endl;
	//this_dots_covered_by_other.print( std::cout );

	//TR_RD << "increment_this_and_cache(): others_dots_covered_by_this: " << std::endl;
	//this_overlap_on_other.print( std::cout );

	// increment_from_cached will update both regular and expanded polar SASA. the hard part is getting get_overlap_cache()
	// to do it right.
	increment_from_cached( this_dots_covered_by_other );

}

///
/// @brief
/// computes the overlap each rotamer (this & other) induce on each other and stores that information in the RotamerDotsCache objects
///
/// @param
/// other - [in] - the other RotamerDots object
/// other_dots_covered_by_this - [out] - the Cache for the dots on the surface of other that are covered by the atoms on this rotamer
/// this_dots_covered_by_other - [out] - the Cache for the dots on the surface of this that are covered by the atoms on the other rotamer
/// atom_atom_overlaps_cache   - [out] - holds a boolean indicating whether two atoms have overlapping, solvent-exposed surface area
///
void RotamerDots::get_overlap_cache(
	RotamerDots const & other,
	RotamerDotsCache & others_dots_covered_by_this,
	RotamerDotsCache & this_dots_covered_by_other,
	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
) const {

	using namespace utility;

	//apl static so that they will be allocated exactly once - we don't want to use new and delete in the inner most loop.
	static vector1< vector1< ObjexxFCL::ubyte > > vv_this_covered_by_other;
	static vector1< vector1< ObjexxFCL::ubyte > > vv_other_covered_by_this;

	// clear and resize takes alot of time; just make resize calls instead and zero out indices as necessary
	//vv_this_covered_by_other.clear();
	//vv_this_covered_by_other.resize( num_atoms_, vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
	//vv_other_covered_by_this.clear();
	//vv_other_covered_by_this.resize( other.get_num_atoms(), vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );

	if ( vv_this_covered_by_other.size() < num_atoms_ ) {
		vv_this_covered_by_other.resize( num_atoms_, vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
	}
	if ( vv_other_covered_by_this.size() < other.get_num_atoms() ) {
		vv_other_covered_by_this.resize( other.get_num_atoms(), vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
	}

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		for ( Size jj = 1; jj <= num_bytes_; ++jj ) {
			vv_this_covered_by_other[ ii ][ jj ] = 0;
		}
	}
	for ( Size ii = 1; ii <= other.num_atoms_; ++ii ) {
		for ( Size jj = 1; jj <= num_bytes_; ++jj ) {
			vv_other_covered_by_this[ ii ][ jj ] = 0;
		}
	}

	// calls the zero() class method on both of the RotamerDotsCache references passed in.  I guess we do that to invalidate the cache
	// before storing the new cache values.
	others_dots_covered_by_this.zero();
	this_dots_covered_by_other.zero();

	get_res_res_overlap( other, vv_this_covered_by_other, vv_other_covered_by_this, atom_atom_overlaps_cache );

	this_dots_covered_by_other.increment_count( vv_this_covered_by_other );
	others_dots_covered_by_this.increment_count( vv_other_covered_by_this );

}

///
/// @brief
/// Calls get_atom_atom_coverage for all atom pairs between res1 and res2. This method gets called by RotamerDots::get_overlap_cache().
///
void RotamerDots::get_res_res_overlap(
	RotamerDots const & other_res,
	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & res1_covered_by_res2, // may include more entries than atoms
	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & res2_covered_by_res1, // may include more entries than atoms
	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
) const {

	Real square_distance = 0.0f;
	bool overlap;
	for ( Size ii = 1; ii <= num_atoms_; ++ii ) { // 'this'/self residues atoms

		Vector ii_atom_xyz = get_atom_coords_xyz( ii );
		Real ii_atom_radius;
		ii_atom_radius = get_atom_radius( ii );

		for ( Size jj = 1; jj <= other_res.get_num_atoms(); ++jj ) {

			Vector jj_atom_xyz = other_res.get_atom_coords_xyz( jj );
			Real jj_atom_radius;
			jj_atom_radius = other_res.get_atom_radius( jj );

			overlap = get_atom_atom_coverage( ii_atom_xyz, ii_atom_radius, jj_atom_xyz, jj_atom_radius, res1_covered_by_res2[ ii ], res2_covered_by_res1[ jj ], square_distance );

			// set the overlaps bool if overlap was found. the outer vector holds the changing node's atoms (or other_res, as it
			// is referred to in this function) and the inner vector is for this RD object's atoms.
			if ( overlap )
				atom_atom_overlaps_cache[ jj ][ ii ] = true;

			if ( square_distance > (ii_atom_radius + jj_atom_radius) * (ii_atom_radius + jj_atom_radius) ) {
				break;
			}
		}
	}
}

///
/// @brief
/// returns false if the two spheres do not overlap at all. otherwise, saves the overlap masks to the input vectors.
///
bool RotamerDots::get_atom_atom_coverage( Vector const & at1_xyz, Real at1_base_radius,
	Vector const & at2_xyz, Real at2_base_radius,
	utility::vector1< ObjexxFCL::ubyte > & at1_sphere_covered,
	utility::vector1< ObjexxFCL::ubyte > & at2_sphere_covered, Real dist_sq ) const {

	int degree_of_overlap;
	int aphi_1_2, aphi_2_1;
	int theta_1_2, theta_2_1;
	int masknum;

	Real at1_radius = at1_base_radius + probe_radius_;
	Real at2_radius = at2_base_radius + probe_radius_;

	// exit if large probe radii do not touch
	dist_sq = at1_xyz.distance_squared( at2_xyz );
	if ( dist_sq > (at1_radius + at2_radius) * (at1_radius + at2_radius) ) {
		return false;
	}
	Real const distance = std::sqrt( dist_sq );

	//ronj this block represents the amount of surface area covered up on atom1 by atom2
	core::scoring::sasa::get_legrand_atomic_overlap( at1_radius, at2_radius, distance, degree_of_overlap );
	core::scoring::sasa::get_legrand_2way_orientation( at1_xyz, at2_xyz, aphi_1_2, theta_1_2, aphi_2_1, theta_2_1, distance );

	Size closest_dot1 = (*lg_angles_)( aphi_1_2, theta_1_2 );
	masknum = ( closest_dot1 * 100 ) + degree_of_overlap;
	for ( Size bb = 1, bbli = (*lg_masks_).index( bb, masknum ); bb <= num_bytes_; ++bb, ++bbli ) {
		at1_sphere_covered[ bb ] |= (*lg_masks_)[ bbli ];
	}

	//ronj the amount of surface area covered up on atom2 by atom1
	core::scoring::sasa::get_legrand_atomic_overlap( at2_radius, at1_radius, distance, degree_of_overlap );

	Size closest_dot2 = (*lg_angles_)( aphi_2_1, theta_2_1 );
	masknum = ( closest_dot2 * 100 ) + degree_of_overlap;
	for ( Size bb = 1, bbli = (*lg_masks_).index( bb, masknum ); bb <= num_bytes_; ++bb, ++bbli ) {
		at2_sphere_covered[ bb ] |= (*lg_masks_)[ bbli ];
	}

	return true;
}

///
/// @brief
/// Increments the dot coverage count for this rotamer from a coverage cache
///
void RotamerDots::increment_from_cached( RotamerDotsCache const & cache ) {

	if ( cache.atom_counts_.size() != atom_counts_.size() )
		TR_RD << "increment_from_cached(): cache.size(): " << cache.atom_counts_.size() << ", atom_counts_.size(): " << atom_counts_.size() << std::endl;
debug_assert( cache.atom_counts_.size() == atom_counts_.size() );
	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		atom_counts_[ ii ] += cache.atom_counts_[ ii ];
	}

	sasa_is_current_ = false;
}

///
/// @brief
/// decrements the dot coverage count for this by the coverage stored in the input RotamerDotsCache object
///
void RotamerDots::decrement_from_cached( RotamerDotsCache const & cache ) {

debug_assert( cache.atom_counts_.size() == atom_counts_.size() );
	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		atom_counts_[ ii ] -= cache.atom_counts_[ ii ];
	}

	sasa_is_current_ = false;
}

///
/// @brief
/// Add rotamer coverage counts for dots on both this and other. sets sasa_is_current_ to false on both this and rhs
///
/// @details
/// One use case involves BGNodes initializing overlap with other BGNodes. This is a brute force all BGNode v all BGNode
/// pairwise increment on the RotamerDots objects each Node holds. That's why we call increment *both*.
/// We don't care about the cache values in this case, but to use increment_both_and_cache, we have to create Cache variables
/// to use as references.
///
void RotamerDots::increment_both( RotamerDots & other ) {

	RotamerDotsCache others_dots_covered_by_this( other.get_num_atoms() );
	others_dots_covered_by_this.zero();

	RotamerDotsCache this_dots_covered_by_other( num_atoms_ );
	this_dots_covered_by_other.zero();

	// can't make the variables above be static because this function gets called by all kinds of FC nodes, which will
	// have varying numbers of atoms.  Then the member vector inside RDC would not be sized correctly and errors would
	// occur. But RDC object are pretty lightweight, so it shouldn't result in a big performance hit.
	utility::vector1< utility::vector1< bool > > dummy( other.get_num_atoms(), utility::vector1< bool >( num_atoms_, false ) );
	increment_both_and_cache( other, others_dots_covered_by_this, this_dots_covered_by_other, dummy );
}

///
/// @details
/// Add rotamer coverage counts for dots on both this and other. Cache's the overlap this and other have with each other for greater efficiency.
/// The second parameter are the dots on the surface of other_rotamer that this overlaps. The third parameter are the dots on the surface of this
/// that other_rotamer overlaps. The fourth parameter lives on the Edges of the IG and stores a boolean indicating whether two atoms have
/// overlapping, solvent-exposed surface area. Instead of recalculating this for all residue pairs every substitution, keep it in this cache.
/// The vectors are already sized. The outer vector has the "other_rotamer"'s atoms and the inner vector is for the atoms in this class's residue.
/// The structure just gets passed down to get_overlap_cache() and it eventually gets filled in get_atom_atom_coverage().
///
/// If the class is keeping the expanded polar atom SASA also, then this function makes an extra call to get_overlap_cache()
/// to get the extra overlap that happens when polar atom radii are expanded.
///
void RotamerDots::increment_both_and_cache( RotamerDots & other_rotamer, RotamerDotsCache & others_dots_covered_by_this,
	RotamerDotsCache & this_dots_covered_by_other, utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache ) {

	get_overlap_cache( other_rotamer, others_dots_covered_by_this, this_dots_covered_by_other, atom_atom_overlaps_cache );

	//TR_RD << "increment_both_and_cache(): overlap found:" << std::endl;
	//this_dots_covered_by_other.print( std::cout );
	//others_dots_covered_by_this.print( std::cout );

	increment_from_cached( this_dots_covered_by_other );
	other_rotamer.increment_from_cached( others_dots_covered_by_this );

}

///
/// @brief
/// Given the current dot coverage counts, returns the total SASA of this residue.
///
/// @details
/// This method does not do any work figuring out who is overlapping with this rotamer. It assumes that work has been
/// done. Instead, it returns the SASA of the dot counts currently held.  If the dot counts have not changed since the
/// last time get_sasa() got called then sasa_is_current_ will be true.  In that case, the method will just return the
/// the variable sasa_. If the counts have changed, it iterates over all the atoms and recalculates the total SASA.
/// That value is stored in sasa_ and sasa_is_current_ is set to true.
///
/// The reason both sasa_ and sasa_is_current_ are "mutable" is so that they can be modified inside this const function.
///
Real RotamerDots::get_sasa() const {

	if ( state_unassigned() )
		return 0.0f;

	if ( sasa_is_current_ )
		return sasa_;

	Real fraction_uncovered = 0.0;
	Real atom_radius = 0.0;
	Real const four_pi = 4.0 * Real( numeric::constants::d::pi );
	Real atom_area_exposed = 0.0;

	sasa_ = 0.0;
	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
		fraction_uncovered = static_cast< Real >( get_num_uncovered( ii ) ) / atom_counts_[ii].get_total_dots();
		atom_radius = get_atom_radius( ii ) + probe_radius_;
		atom_area_exposed = four_pi * ( atom_radius * atom_radius ) * fraction_uncovered;
		//std::cout << "get_sasa(): atom: " << ii << ", rad: " << atom_radius << ", %-uncovered: " << fraction_uncovered << ", exposed: " << atom_area_exposed << std::endl;

		atom_sasa_[ ii ] = atom_area_exposed;
		sasa_ += atom_area_exposed;
	}

	sasa_is_current_ = true;

	return sasa_;
}

///
/// @brief
/// Given the current dot coverage counts, returns the total SASA for a particular atom index.
/// Assumes that get_atom_sasa() will never be called when the object is in the unassigned state.
///
Real RotamerDots::get_atom_sasa( Size atom_index ) const {
	if ( ! sasa_is_current_ )
		get_sasa();

	return atom_sasa_[ atom_index ];
}

///
/// @brief
/// Returns the number of uncovered dots on the given atom, when using standard SASA radii.
/// Note: no expanded polars version of this method.
///
Size RotamerDots::get_num_uncovered( Size atom ) const {
	return atom_counts_[ atom ].get_num_uncovered();
}

///
/// @brief
/// Note: no expanded polars version of this method.
///
Size RotamerDots::get_num_covered_total() const {

	Size total_num_covered = 0;
	for ( Size ii=1; ii <= atom_counts_.size(); ++ii ) {
		total_num_covered += atom_counts_[ ii ].get_num_covered();
	}

	return total_num_covered;
}

///
/*void RotamerDots::write_dot_kinemage( std::ofstream & kinfile ) {

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		write_dotlist_header( kinfile, "1.4 exposed dots", "red");
		for ( Size jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj ) {
			if ( ! atom_counts_[ii].get_dot_covered( jj ) ) {
				write_dot( kinfile, ii, jj, probe_radius_ );
			}
		}
		write_dotlist_header( kinfile, "SA", "red");
		for ( Size jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj ) {
			if ( ! atom_counts_[ii].get_dot_covered( jj ) ) {
				write_dot( kinfile, ii, jj, 0 );
			}
		}
		write_dotlist_header( kinfile, "1.0 A probe Accessible", "green");
		for ( Size jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj ) {
			if ( ! atom_counts_[ii].get_dot_covered( jj ) ) {
				write_dot( kinfile, ii, jj, 0 );
			}
		}
		//write_dotlist_header( kinfile, "void surface", "blue");
		//for (int jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj) {
		//	if ( atom_counts_[ii].get_dot_covered( jj ) && ! atom_counts_small_probe_[ii].get_dot_covered( jj ) ) {
		//		write_dot( kinfile, ii, jj, 0 );
		//	}
		//}
		//write_dotlist_header( kinfile, "all dots", "white");
		//for (int jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj)
		//{
		//	write_dot( kinfile, ii, jj, 0 );
		//}
	}
}

///
void RotamerDots::write_dotlist_header( std::ofstream & kinfile, std::string master_name, std::string color ) {
	kinfile << "@dotlist master= {" << master_name << "} color= " << color << "\n";
}

///
void RotamerDots::write_dot( std::ofstream & kinfile, Size atom, Size dot, Real radius ) {
	static numeric::xyzVector< Real > coord;
	coord = get_atom_coords_xyz( atom );
	coord += (radius + get_atom_radius( atom )) * get_dot_coord( dot );

	write_dot( kinfile, coord, "dot" );
}

///
void RotamerDots::write_dot( std::ofstream & kinfile, numeric::xyzVector< Real > const & coord, std::string atname ) {
	static std::string last_atname = "";
	if ( last_atname == atname ) {
		kinfile << "{\"} P " << coord.x() << " " << coord.y() << " " << coord.z() << "\n";
	} else {
		kinfile << "{" << atname << "} P " << coord.x() << " " << coord.y() << " " << coord.z() << "\n";
		last_atname = atname;
	}
}

///
numeric::xyzVector<Real> const & RotamerDots::get_dot_coord( Size dot_id ) {
	if ( ! dot_sphere_coordinates_initialized_ ) {
		initialize_dot_sphere_coordinates_from_file();
	}
	return dot_sphere_coordinates_[ dot_id ];
}

///
void RotamerDots::initialize_dot_sphere_coordinates_from_file() {

	dot_sphere_coordinates_.resize( DotSphere::NUM_DOTS_TOTAL );

	std::ifstream dotfile("sphere.txt");
	for ( Size ii = 1; ii <= DotSphere::NUM_DOTS_TOTAL; ++ii ) {
		Real x, y, z;
		dotfile >> x >> y >> z;
		dot_sphere_coordinates_[ ii ] = -1 * numeric::xyzVector< Real >( x,y,z );
	}
	dotfile.close();

	dot_sphere_coordinates_initialized_ = true;
}
*/

///
void RotamerDots::print( std::ostream & os ) const {

	if ( state_unassigned() ) {
		os << "dots: unassigned" << std::endl;
		return;
	}

	os << "dots: " << rotamer_->name3() << rotamer_->seqpos() << std::endl;
	for ( Size ii=1; ii <= num_atoms_; ++ii ) {
		os << "atom " << rotamer_->atom_name( ii ) << ", rad: " << get_atom_radius( ii )
			<< ", covered: " << ObjexxFCL::format::I(3,atom_counts_[ii].get_num_covered()) << ", counts: ";
		atom_counts_[ ii ].print( os );
	}
	os << "num_covered: " << get_num_covered_total() << ", sasa_is_current_: " << sasa_is_current_ << ", sasa: " << get_sasa() << std::endl;

}


std::string RotamerDots::name3() const { return rotamer_->name3(); } // for operator<<
core::Size  RotamerDots::seqpos() const { return rotamer_->seqpos(); } // for operator<<

///
/// @brief
/// invokes print on the input RotamerDots object
///
std::ostream & operator<< ( std::ostream & os, RotamerDots const & rd ) {
	if ( rd.state_unassigned() ) {
		os << "unassigned";
		return os;
	}
	os << "rotamer: " << rd.name3() << rd.seqpos() << ", sasa: " << rd.get_sasa();
	return os;
}

/// @details The contents of "sphere.txt", now hard coded.
/// The coordinates in the original file are all on the opposite size of the
/// unit sphere from where they genuinely appear to be.
void
RotamerDots::initialize_dot_coords( utility::vector1< core::Vector > & dot_coords )
{
	dot_coords.resize( 162 );
	dot_coords[1] = -1* Vector(0,0,1);
	dot_coords[2] = -1* Vector(0.276393,0.850651,0.447214);
	dot_coords[3] = -1* Vector(0.894427,0,0.447214);
	dot_coords[4] = -1* Vector(0.16246,0.5,0.850651);
	dot_coords[5] = -1* Vector(0.525731,0,0.850651);
	dot_coords[6] = -1* Vector(0.0776089,0.238856,0.967949);
	dot_coords[7] = -1* Vector(0.251148,0,0.967949);
	dot_coords[8] = -1* Vector(0.361803,0.262866,0.894427);
	dot_coords[9] = -1* Vector(0.688191,0.5,0.525731);
	dot_coords[10] = -1* Vector(0.232827,0.716567,0.657513);
	dot_coords[11] = -1* Vector(0.447214,0.525731,0.723607);
	dot_coords[12] = -1* Vector(0.483974,0.716567,0.502295);
	dot_coords[13] = -1* Vector(0.638197,0.262866,0.723607);
	dot_coords[14] = -1* Vector(0.753443,0,0.657513);
	dot_coords[15] = -1* Vector(0.831052,0.238856,0.502295);
	dot_coords[16] = -1* Vector(-0.723607,0.525731,0.447214);
	dot_coords[17] = -1* Vector(-0.425325,0.309017,0.850651);
	dot_coords[18] = -1* Vector(-0.203183,0.147621,0.967949);
	dot_coords[19] = -1* Vector(-0.138197,0.425325,0.894427);
	dot_coords[20] = -1* Vector(-0.262866,0.809017,0.525731);
	dot_coords[21] = -1* Vector(-0.609548,0.442863,0.657513);
	dot_coords[22] = -1* Vector(-0.361803,0.587785,0.723607);
	dot_coords[23] = -1* Vector(-0.531939,0.681718,0.502295);
	dot_coords[24] = -1* Vector(-0.0527864,0.688191,0.723607);
	dot_coords[25] = -1* Vector(0.029644,0.864188,0.502295);
	dot_coords[26] = -1* Vector(-0.723607,-0.525731,0.447214);
	dot_coords[27] = -1* Vector(-0.425325,-0.309017,0.850651);
	dot_coords[28] = -1* Vector(-0.203183,-0.147621,0.967949);
	dot_coords[29] = -1* Vector(-0.447214,0,0.894427);
	dot_coords[30] = -1* Vector(-0.850651,0,0.525731);
	dot_coords[31] = -1* Vector(-0.609548,-0.442863,0.657513);
	dot_coords[32] = -1* Vector(-0.67082,-0.16246,0.723607);
	dot_coords[33] = -1* Vector(-0.812731,-0.295242,0.502295);
	dot_coords[34] = -1* Vector(-0.67082,0.16246,0.723607);
	dot_coords[35] = -1* Vector(-0.812731,0.295242,0.502295);
	dot_coords[36] = -1* Vector(0.276393,-0.850651,0.447214);
	dot_coords[37] = -1* Vector(0.16246,-0.5,0.850651);
	dot_coords[38] = -1* Vector(0.0776089,-0.238856,0.967949);
	dot_coords[39] = -1* Vector(-0.138197,-0.425325,0.894427);
	dot_coords[40] = -1* Vector(-0.262866,-0.809017,0.525731);
	dot_coords[41] = -1* Vector(0.232827,-0.716567,0.657513);
	dot_coords[42] = -1* Vector(-0.0527864,-0.688191,0.723607);
	dot_coords[43] = -1* Vector(0.029644,-0.864188,0.502295);
	dot_coords[44] = -1* Vector(-0.361803,-0.587785,0.723607);
	dot_coords[45] = -1* Vector(-0.531939,-0.681718,0.502295);
	dot_coords[46] = -1* Vector(0.361803,-0.262866,0.894427);
	dot_coords[47] = -1* Vector(0.688191,-0.5,0.525731);
	dot_coords[48] = -1* Vector(0.638197,-0.262866,0.723607);
	dot_coords[49] = -1* Vector(0.831052,-0.238856,0.502295);
	dot_coords[50] = -1* Vector(0.447214,-0.525731,0.723607);
	dot_coords[51] = -1* Vector(0.483974,-0.716567,0.502295);
	dot_coords[52] = -1* Vector(0.723607,0.525731,-0.447214);
	dot_coords[53] = -1* Vector(0.951057,0.309017,0);
	dot_coords[54] = -1* Vector(0.956626,0.147621,0.251148);
	dot_coords[55] = -1* Vector(0.861803,0.425325,0.276393);
	dot_coords[56] = -1* Vector(0.587785,0.809017,0);
	dot_coords[57] = -1* Vector(0.67082,0.688191,0.276393);
	dot_coords[58] = -1* Vector(0.436009,0.864188,0.251148);
	dot_coords[59] = -1* Vector(0.809017,0.587785,0);
	dot_coords[60] = -1* Vector(0.860696,0.442863,-0.251148);
	dot_coords[61] = -1* Vector(0.687157,0.681718,-0.251148);
	dot_coords[62] = -1* Vector(-0.276393,0.850651,-0.447214);
	dot_coords[63] = -1* Vector(0,1,0);
	dot_coords[64] = -1* Vector(0.155218,0.955423,0.251148);
	dot_coords[65] = -1* Vector(-0.138197,0.951057,0.276393);
	dot_coords[66] = -1* Vector(-0.587785,0.809017,0);
	dot_coords[67] = -1* Vector(-0.447214,0.850651,0.276393);
	dot_coords[68] = -1* Vector(-0.687157,0.681718,0.251148);
	dot_coords[69] = -1* Vector(-0.309017,0.951057,0);
	dot_coords[70] = -1* Vector(-0.155218,0.955423,-0.251148);
	dot_coords[71] = -1* Vector(-0.436009,0.864188,-0.251148);
	dot_coords[72] = -1* Vector(-0.894427,0,-0.447214);
	dot_coords[73] = -1* Vector(-0.951057,0.309017,0);
	dot_coords[74] = -1* Vector(-0.860696,0.442863,0.251148);
	dot_coords[75] = -1* Vector(-0.947214,0.16246,0.276393);
	dot_coords[76] = -1* Vector(-0.951057,-0.309017,0);
	dot_coords[77] = -1* Vector(-0.947214,-0.16246,0.276393);
	dot_coords[78] = -1* Vector(-0.860696,-0.442863,0.251148);
	dot_coords[79] = -1* Vector(-1,0,0);
	dot_coords[80] = -1* Vector(-0.956626,0.147621,-0.251148);
	dot_coords[81] = -1* Vector(-0.956626,-0.147621,-0.251148);
	dot_coords[82] = -1* Vector(-0.276393,-0.850651,-0.447214);
	dot_coords[83] = -1* Vector(-0.587785,-0.809017,0);
	dot_coords[84] = -1* Vector(-0.687157,-0.681718,0.251148);
	dot_coords[85] = -1* Vector(-0.447214,-0.850651,0.276393);
	dot_coords[86] = -1* Vector(0,-1,0);
	dot_coords[87] = -1* Vector(-0.138197,-0.951057,0.276393);
	dot_coords[88] = -1* Vector(0.155218,-0.955423,0.251148);
	dot_coords[89] = -1* Vector(-0.309017,-0.951057,0);
	dot_coords[90] = -1* Vector(-0.436009,-0.864188,-0.251148);
	dot_coords[91] = -1* Vector(-0.155218,-0.955423,-0.251148);
	dot_coords[92] = -1* Vector(0.723607,-0.525731,-0.447214);
	dot_coords[93] = -1* Vector(0.587785,-0.809017,0);
	dot_coords[94] = -1* Vector(0.436009,-0.864188,0.251148);
	dot_coords[95] = -1* Vector(0.67082,-0.688191,0.276393);
	dot_coords[96] = -1* Vector(0.951057,-0.309017,0);
	dot_coords[97] = -1* Vector(0.861803,-0.425325,0.276393);
	dot_coords[98] = -1* Vector(0.956626,-0.147621,0.251148);
	dot_coords[99] = -1* Vector(0.809017,-0.587785,0);
	dot_coords[100] = -1* Vector(0.687157,-0.681718,-0.251148);
	dot_coords[101] = -1* Vector(0.860696,-0.442863,-0.251148);
	dot_coords[102] = -1* Vector(0.262866,0.809017,-0.525731);
	dot_coords[103] = -1* Vector(0.531939,0.681718,-0.502295);
	dot_coords[104] = -1* Vector(0.447214,0.850651,-0.276393);
	dot_coords[105] = -1* Vector(0.309017,0.951057,0);
	dot_coords[106] = -1* Vector(0.138197,0.951057,-0.276393);
	dot_coords[107] = -1* Vector(-0.029644,0.864188,-0.502295);
	dot_coords[108] = -1* Vector(-0.688191,0.5,-0.525731);
	dot_coords[109] = -1* Vector(-0.483974,0.716567,-0.502295);
	dot_coords[110] = -1* Vector(-0.67082,0.688191,-0.276393);
	dot_coords[111] = -1* Vector(-0.809017,0.587785,0);
	dot_coords[112] = -1* Vector(-0.861803,0.425325,-0.276393);
	dot_coords[113] = -1* Vector(-0.831052,0.238856,-0.502295);
	dot_coords[114] = -1* Vector(-0.688191,-0.5,-0.525731);
	dot_coords[115] = -1* Vector(-0.831052,-0.238856,-0.502295);
	dot_coords[116] = -1* Vector(-0.861803,-0.425325,-0.276393);
	dot_coords[117] = -1* Vector(-0.809017,-0.587785,0);
	dot_coords[118] = -1* Vector(-0.67082,-0.688191,-0.276393);
	dot_coords[119] = -1* Vector(-0.483974,-0.716567,-0.502295);
	dot_coords[120] = -1* Vector(0.262866,-0.809017,-0.525731);
	dot_coords[121] = -1* Vector(-0.029644,-0.864188,-0.502295);
	dot_coords[122] = -1* Vector(0.138197,-0.951057,-0.276393);
	dot_coords[123] = -1* Vector(0.309017,-0.951057,0);
	dot_coords[124] = -1* Vector(0.447214,-0.850651,-0.276393);
	dot_coords[125] = -1* Vector(0.531939,-0.681718,-0.502295);
	dot_coords[126] = -1* Vector(0.850651,0,-0.525731);
	dot_coords[127] = -1* Vector(0.812731,-0.295242,-0.502295);
	dot_coords[128] = -1* Vector(0.947214,-0.16246,-0.276393);
	dot_coords[129] = -1* Vector(1,0,0);
	dot_coords[130] = -1* Vector(0.947214,0.16246,-0.276393);
	dot_coords[131] = -1* Vector(0.812731,0.295242,-0.502295);
	dot_coords[132] = -1* Vector(0,0,-1);
	dot_coords[133] = -1* Vector(0.425325,0.309017,-0.850651);
	dot_coords[134] = -1* Vector(0.609548,0.442863,-0.657513);
	dot_coords[135] = -1* Vector(0.361803,0.587785,-0.723607);
	dot_coords[136] = -1* Vector(-0.16246,0.5,-0.850651);
	dot_coords[137] = -1* Vector(0.0527864,0.688191,-0.723607);
	dot_coords[138] = -1* Vector(-0.232827,0.716567,-0.657513);
	dot_coords[139] = -1* Vector(0.138197,0.425325,-0.894427);
	dot_coords[140] = -1* Vector(0.203183,0.147621,-0.967949);
	dot_coords[141] = -1* Vector(-0.0776089,0.238856,-0.967949);
	dot_coords[142] = -1* Vector(-0.447214,0.525731,-0.723607);
	dot_coords[143] = -1* Vector(-0.525731,0,-0.850651);
	dot_coords[144] = -1* Vector(-0.638197,0.262866,-0.723607);
	dot_coords[145] = -1* Vector(-0.753443,0,-0.657513);
	dot_coords[146] = -1* Vector(-0.361803,0.262866,-0.894427);
	dot_coords[147] = -1* Vector(-0.251148,0,-0.967949);
	dot_coords[148] = -1* Vector(-0.638197,-0.262866,-0.723607);
	dot_coords[149] = -1* Vector(-0.16246,-0.5,-0.850651);
	dot_coords[150] = -1* Vector(-0.447214,-0.525731,-0.723607);
	dot_coords[151] = -1* Vector(-0.232827,-0.716567,-0.657513);
	dot_coords[152] = -1* Vector(-0.361803,-0.262866,-0.894427);
	dot_coords[153] = -1* Vector(-0.0776089,-0.238856,-0.967949);
	dot_coords[154] = -1* Vector(0.0527864,-0.688191,-0.723607);
	dot_coords[155] = -1* Vector(0.425325,-0.309017,-0.850651);
	dot_coords[156] = -1* Vector(0.361803,-0.587785,-0.723607);
	dot_coords[157] = -1* Vector(0.609548,-0.442863,-0.657513);
	dot_coords[158] = -1* Vector(0.138197,-0.425325,-0.894427);
	dot_coords[159] = -1* Vector(0.203183,-0.147621,-0.967949);
	dot_coords[160] = -1* Vector(0.67082,-0.16246,-0.723607);
	dot_coords[161] = -1* Vector(0.67082,0.16246,-0.723607);
	dot_coords[162] = -1* Vector(0.447214,0,-0.894427);
}


//----------------------------------------------------------------------------//
//------------------------- RotamerDotsRadiusData  ---------------------------//
//----------------------------------------------------------------------------//

RotamerDotsRadiusData *
RotamerDotsRadiusData::create_singleton_instance()
{
	return new RotamerDotsRadiusData;
}

/// @brief private constructor to guarantee the singleton
RotamerDotsRadiusData::RotamerDotsRadiusData()
{
	// Previously this was a static data member of class RotamerDots
	polar_expansion_radius_ = 1.0;

	// 1. Initialize SASA radii
	using namespace core::chemical;
	TR_RDRD << "RotamerDotsRadiusData::RotamerDotsRadiusData(): reading in sasa radii database file" << std::endl;

	//j setup the radii array, indexed by the atom type int. atom index for looking up an extra data type stored in the AtomTypes
	//ronj reads the values out of the database file sasa_radii.txt in the extras folder of atom_type_sets and stores the values
	//ronj for each atom type into the radii array. each index of the radii array corresponds to some atom type.
	AtomTypeSet const & atom_type_set = * ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	ROSETTA_SASA_radii_.resize( atom_type_set.n_atomtypes(), 0.0 );

	core::Size SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "SASA_RADIUS_LEGACY" );

	TR_RDRD << "ROSETTA_SASA_radii_: [ ";
	for ( core::Size ii=1; ii <= atom_type_set.n_atomtypes(); ++ii ) {
		ROSETTA_SASA_radii_[ ii ] = atom_type_set[ ii ].extra_parameter( SASA_RADIUS_INDEX );
		TR_RDRD << ROSETTA_SASA_radii_[ ii ] << ", ";
	}
	TR_RDRD << "]" << std::endl;

	// 2. Initialize the Naccess radii
	TR_RDRD << "RotamerDotsRadiusData::RotamerDotsRadiusData(): reading in sasa radii database file" << std::endl;

	//AtomTypeSet const & atom_type_set = * ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	NACCESS_SASA_radii_.resize( atom_type_set.n_atomtypes(), 0.0 );

	core::Size NACCESS_SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "NACCESS_SASA_RADIUS" );

	TR_RDRD << "NACCESS_SASA_radii_: [ ";
	for ( core::Size ii=1; ii <= atom_type_set.n_atomtypes(); ++ii ) {
		NACCESS_SASA_radii_[ ii ] = atom_type_set[ ii ].extra_parameter( NACCESS_SASA_RADIUS_INDEX );
		TR_RDRD << NACCESS_SASA_radii_[ ii ] << ", ";
	}
	TR_RDRD << "]" << std::endl;

	// 3. Initialize the expanded-polar Naccess radii
	TR_RDRD << "get_NACCESS_SASA_radii_with_expanded_polars(): reading in sasa radii database file" << std::endl;

	NACCESS_SASA_radii_with_expanded_polars_.resize( NACCESS_SASA_radii_.size() );

	// used to figure out what atom_type a given index is
	//AtomTypeSet const & atom_type_set = * ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );

	TR_RDRD << "NACCESS_SASA_radii_with_expanded_polars_: [ ";
	for ( Size ii=1; ii <= NACCESS_SASA_radii_.size(); ++ii ) {
		NACCESS_SASA_radii_with_expanded_polars_[ ii ] = NACCESS_SASA_radii_[ ii ];

		core::chemical::AtomType const & at( atom_type_set[ ii ] );
		if ( at.element() == "N" || at.element() == "O" ) {
			NACCESS_SASA_radii_with_expanded_polars_[ ii ] += polar_expansion_radius_;
		}
		TR_RDRD << NACCESS_SASA_radii_with_expanded_polars_[ ii ] << ", ";
	}
	TR_RDRD << "]" << std::endl;

}

utility::vector1< Real > const *
RotamerDotsRadiusData::get_ROSETTA_SASA_radii() const
{
	return &ROSETTA_SASA_radii_;
}

utility::vector1< Real > const *
RotamerDotsRadiusData::get_NACCESS_SASA_radii() const
{
	return &NACCESS_SASA_radii_;
}

utility::vector1< Real > const *
RotamerDotsRadiusData::get_NACCESS_SASA_radii_with_expanded_polars() const
{
	return &NACCESS_SASA_radii_with_expanded_polars_;

}


//----------------------------------------------------------------------------//
//------------------------- Rotamer Dots Cache Class -------------------------//
//----------------------------------------------------------------------------//

///
RotamerDotsCache::RotamerDotsCache() {}

///
RotamerDotsCache::RotamerDotsCache( Size num_atoms ) {
	atom_counts_.resize( num_atoms );
}

///
/// @brief
/// copy constructor
///
RotamerDotsCache::RotamerDotsCache( RotamerDotsCache const & rhs ) :
	atom_counts_( rhs.atom_counts_ )
{}

///
RotamerDotsCache::~RotamerDotsCache() {}

///
/// @brief
/// assignment operator
///
RotamerDotsCache const & RotamerDotsCache::operator=( RotamerDotsCache const & rhs ) {
	atom_counts_ = rhs.atom_counts_;
	return *this;
}

///
void RotamerDotsCache::resize( Size num_atoms ) {
	//TR_RDC << "resize() called with num_atoms: " << num_atoms << std::endl;
	atom_counts_.clear();
	atom_counts_.resize( num_atoms );
}

///
/// @brief
/// sets the dot counts to zero for all atoms
///
/// @details
/// if the cache already knows the atom's counts are uniformly 0, it skips it
///
void
RotamerDotsCache::zero() {
	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
		//if ( non_zero_overlap_[ ii ] ) {
			atom_counts_[ ii ].zero(); // calls zero() on each DotSphere instance
		//}
	}
	//atom_counts_.clear(); // leaving this in causes the vector to be sized down to 0, causing problems
}

///
/// @brief
/// increments the dot coverage counts for all atoms in the cache
///
/// @param
/// covered - [in] - compact ubyte array of dots; '1' for any covered dot for the vdw + 1.4 A sphere
///
void RotamerDotsCache::increment_count( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & covered ) {

	//TR_RDC << "increment_count(): atom_counts_.size(): " << atom_counts_.size() << ", covered.size(): " << covered.size() << std::endl;
	// do not assert this -- let there be more entries, possibly, in covered array;
	//debug_assert( atom_counts_.size() == covered.size() );
debug_assert( atom_counts_.size() <= covered.size() );
	for ( Size ii = 1, ii_end = atom_counts_.size(); ii <= ii_end; ++ii ) {
		//utility::vector1< ObjexxFCL::ubyte > const & atom_mask = covered[ ii ];
		atom_counts_[ ii ].increment_count( covered[ ii ] );
	}

}

///
/// @details
/// Called by BGEdges with a reference to a vector1 of vector1 of ubytes representing where to put the compact (count)
/// based representation of all the atom overlap masks.
///
/// Note: This method only writes the standard SASA dot counts to the passed in array. No version of this function
/// exists for expanded polar atom SASA dot counts.
///
void RotamerDotsCache::write_to_compact_array( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & compact ) const {
debug_assert( compact.size() != 0 );
	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
		atom_counts_[ ii ].write_to_compact_array( compact[ ii ] );
	}
}

///
void RotamerDotsCache::print_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const {
	for ( Size bb = 1; bb <= values.size(); ++bb ) {
		if ( (bb-1)*8 % 16 == 0 ) std::cout << (bb-1) * 8 << ":";
		for ( int index=7; index >= 0; index-- ) {
			std::cout << ( ( (int)values[ bb ] >> index ) & 1 );
		}
		std::cout << " ";
	}
	std::cout << std::endl;
}

///
void RotamerDotsCache::print( std::ostream & os ) const {
	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
		os << "atom " << I(2,ii) << ": ";
		atom_counts_[ ii ].print( os );
	}
}

//----------------------------------------------------------------------------//
//------------------------------ Inverted Dots Class -------------------------//
//----------------------------------------------------------------------------//

Real const InvRotamerDots::max_dist_from_dot_to_intersection = 0.8;

InvRotamerDots::InvRotamerDots() :
	rotamer_( /* 0 */ )
{}


InvRotamerDots::InvRotamerDots( InvRotamerDots const & src ) :
	utility::pointer::ReferenceCount(),
	rotamer_( src.rotamer_ ),
	inv_dots_( src.inv_dots_ ),
	radii_( src.radii_ )
{}

InvRotamerDots::~InvRotamerDots() {}

InvRotamerDots const &
InvRotamerDots::operator= ( InvRotamerDots const & rhs ) {
	if ( this != & rhs ) {
		rotamer_  = rhs.rotamer_;
		inv_dots_ = rhs.inv_dots_;
		radii_ = rhs.radii_;
	}
	return *this;
}

void
InvRotamerDots::setup_from_rotamer_dots( RotamerDots const & rdots ) {
	rotamer_ = rdots.rotamer();
	if ( ! rotamer_ ) return;
	if ( inv_dots_.size() < rdots.get_num_atoms() ) {
		inv_dots_.resize( rdots.get_num_atoms(), utility::vector1< ObjexxFCL::ubyte >( RotamerDots::num_bytes_, ObjexxFCL::ubyte(0) )  );
	}
	rdots.invert_to_boolmasks( inv_dots_ );
	radii_ = rdots.get_radii();
}

void
InvRotamerDots::setup_from_rotamer_dots(
	RotamerDots const & rdots,
	utility::vector1< Size > const & ats_to_update
)
{
	rotamer_ = rdots.rotamer();
	if ( ! rotamer_ ) return;
	if ( inv_dots_.size() < rdots.get_num_atoms() ) {
		inv_dots_.resize( rdots.get_num_atoms(), utility::vector1< ObjexxFCL::ubyte >( RotamerDots::num_bytes_, ObjexxFCL::ubyte(0) )  );
	}
	rdots.invert_to_boolmasks( inv_dots_, ats_to_update );
	radii_ = rdots.get_radii();
}


core::conformation::ResidueCOP
InvRotamerDots::rotamer() const
{
	return rotamer_;
}

/// @brief Is the intersection between two atoms on this inv-rotamer-dots object exposed?
bool
InvRotamerDots::atom_overlap_is_exposed( Size at1, Size at2 ) const {
debug_assert( rotamer_ );
	return overlap_exposed( rotamer_->atom( at1 ), inv_dots_[ at1 ], rotamer_->atom( at2 ), inv_dots_[ at2 ] );
}

/// @brief Is the intersection between one atom on this inv-rotamer-dots object,
/// and one atom on another inv-rotamer-dots object exposed?
bool
InvRotamerDots::atom_overlap_is_exposed(
	Size at_this,
	InvRotamerDots const & other,
	Size at_other
) const
{
	return overlap_exposed( rotamer_->atom( at_this ), inv_dots_[ at_this ], other.rotamer_->atom( at_other ), other.inv_dots_[ at_other ] );
}


bool InvRotamerDots::dot_exposed( Size atomid, Size dot_index ) const {
debug_assert( dot_index > 0 && dot_index <= 162 );
	dot_index -= 1;
	Size const which_byte = dot_index / 8;
	Size const which_bit  = dot_index - which_byte * 8 ;
	return unpack_ubyte( inv_dots_[ atomid ][ which_byte + 1 ], which_bit );
}

void InvRotamerDots::write_exposed_dots_to_kinemage(
	std::ostream & ostr,
	bool group
) const
{
	if ( ! rotamer_ ) return;

	if ( group ) {
		ostr << "@group { invdots " << rotamer_->seqpos() << "} dominant\n";
	}

	for ( Size ii = 1; ii <= rotamer_->nheavyatoms(); ++ii ) {
		Real const ii_rad = (*radii_)[ rotamer_->atom(ii).type() ] + RotamerDots::probe_radius_;
		write_sphere_list_uv1( ostr, rotamer_->atom_name( ii ), "gray", rotamer_->xyz( ii ), ii_rad, inv_dots_[ ii ] );
	}
}

void
InvRotamerDots::write_circle_intersection_mask_to_kinemage(
	std::ostream & ostr,
	Size const atom_this,
	InvRotamerDots const & invdots_other,
	Size const atom_other,
	bool group
) const
{
	core::conformation::Atom const & at1( rotamer_->atom( atom_this ) );
	core::conformation::Atom const & at2( invdots_other.rotamer_->atom( atom_other ) );

	Real const rad1 = (*radii_)[ at1.type() ] + RotamerDots::probe_radius_;
	Real const rad2 = (*radii_)[ at2.type() ] + RotamerDots::probe_radius_;

	Real const dist_sq = at1.xyz().distance_squared( at2.xyz() );
debug_assert( dist_sq < (rad1 + rad2) * (rad1 + rad2) );
	Real const distance = std::sqrt( dist_sq );

	Real const step_size1 = rad1 * 0.02;
	Real const step_size2 = rad2 * 0.02;

	Size const nsteps1 = (Size)(ceil( max_dist_from_dot_to_intersection / step_size1 ));
	Size const nsteps2 = (Size)(ceil( max_dist_from_dot_to_intersection / step_size2 ));

	int degree_of_overlap1, degree_of_overlap1_stepped, degree_of_overlap2, degree_of_overlap2_stepped;
	int aphi_1_2, aphi_2_1;
	int theta_1_2, theta_2_1;

	//ronj this block represents the amount of surface area covered up on atom1 by atom2
	core::scoring::sasa::get_legrand_atomic_overlap( rad1, rad2, distance, degree_of_overlap1 );
	core::scoring::sasa::get_legrand_2way_orientation( at1.xyz(), at2.xyz(), aphi_1_2, theta_1_2, aphi_2_1, theta_2_1, distance );

	utility::vector1< ObjexxFCL::ubyte > ring1( 21, ObjexxFCL::ubyte(0) );
	utility::vector1< ObjexxFCL::ubyte > ring2( 21, ObjexxFCL::ubyte(0) );

	utility::vector1< ObjexxFCL::ubyte > hit_ring1( 21, ObjexxFCL::ubyte(0) );
	utility::vector1< ObjexxFCL::ubyte > hit_ring2( 21, ObjexxFCL::ubyte(0) );

	Size closest_dot1 = (*RotamerDots::lg_angles_)( aphi_1_2, theta_1_2 );
	if ( degree_of_overlap1 + nsteps1 > 100 ) {
		degree_of_overlap1_stepped = 100;
	} else {
		degree_of_overlap1_stepped = degree_of_overlap1 + nsteps1;
	}

	int const masknum1a = ( closest_dot1 * 100 ) + degree_of_overlap1;
	int const masknum1b = ( closest_dot1 * 100 ) + degree_of_overlap1_stepped;

	// so we take two "offsets" into the "masks" table: 1) the one that normally gets used to figure out which dots are covered
	// by the neighboring atom, and 2) one that's one "step" in, which represents the "ring" of dots that's just past the ones
	// that are covered by the neighboring atom. if we then take the inverse (negate) the normally used mask, we'll get 0's 
	// whereever there are dots covered by the other atom (instead of 1's) and 1's everywhere else. if we logical AND that
	// result with the dots that one "step" in, we'll get 1's at just the ring of dots that's next to the ones that are covered.

	// if we then logical AND the ring of dots with all of the exposed dots on this atom, we can determine if there are any 
	// exposed dots adjacent to the intersection circle.
	
	for ( Size bb = 1, bblia = (*RotamerDots::lg_masks_).index( bb, masknum1a ), bblib = (*RotamerDots::lg_masks_).index( bb, masknum1b );
			bb <= RotamerDots::num_bytes_; ++bb, ++bblia, ++bblib ) {
		ring1[ bb ] = (*RotamerDots::lg_masks_)[ bblib ] & ~ (*RotamerDots::lg_masks_)[ bblia ];
		hit_ring1[ bb ] = inv_dots_[ atom_this ][ bb ] & ( (*RotamerDots::lg_masks_)[ bblib ] & ~ (*RotamerDots::lg_masks_)[ bblia ] );
	}
	std::cout << "exposed dots: ";
	utility::vector1< ObjexxFCL::ubyte > exposed_copy1 = inv_dots_[ atom_this ];
	print_dot_bit_string( exposed_copy1 );
	std::cout << "ring1: ";
	print_dot_bit_string( ring1 );
	std::cout << "hit_ring1: ";
	print_dot_bit_string( hit_ring1 );

	//ronj the amount of surface area covered up on atom2 by atom1
	core::scoring::sasa::get_legrand_atomic_overlap( rad2, rad1, distance, degree_of_overlap2 );

	Size closest_dot2 = (*RotamerDots::lg_angles_)( aphi_2_1, theta_2_1 );
	if ( degree_of_overlap2 + nsteps2 > 100 ) {
		degree_of_overlap2_stepped = 100;
	} else {
		degree_of_overlap2_stepped = degree_of_overlap2 + nsteps2;
	}

	int const masknum2a = ( closest_dot2 * 100 ) + degree_of_overlap2;
	int const masknum2b = ( closest_dot2 * 100 ) + degree_of_overlap2_stepped;

	for ( Size bb = 1, bblia = (*RotamerDots::lg_masks_).index( bb, masknum2a ), bblib = (*RotamerDots::lg_masks_).index( bb, masknum2b );
			bb <= RotamerDots::num_bytes_; ++bb, ++bblia, ++bblib ) {
		ring2[ bb ] = (*RotamerDots::lg_masks_)[ bblib ] & ~ (*RotamerDots::lg_masks_)[ bblia ];
		hit_ring2[ bb ] = invdots_other.inv_dots_[ atom_other ][ bb ] & ( (*RotamerDots::lg_masks_)[ bblib ] & ~ (*RotamerDots::lg_masks_)[ bblia ] );
	}
	std::cout << "exposed dots: ";
	utility::vector1< ObjexxFCL::ubyte > exposed_copy2 = invdots_other.inv_dots_[ atom_other ];
	print_dot_bit_string( exposed_copy2 );
	std::cout << "ring2: ";
	print_dot_bit_string( ring2 );
	std::cout << "hit_ring2: ";
	print_dot_bit_string( hit_ring2 );

	if ( group ) {
		ostr << "@group { invdots " << rotamer_->seqpos() << " " << rotamer_->atom_name( atom_this ) << "} dominant\n";
	}
	write_sphere_list_uv1( ostr, rotamer_->atom_name( atom_this ), "gray", at1.xyz(), rad1, ring1 );

	if ( group ) {
		ostr << "@group { invdots " << invdots_other.rotamer_->seqpos() << " " << invdots_other.rotamer_->atom_name( atom_other ) << "} dominant\n";
	}
	write_sphere_list_uv1( ostr,  invdots_other.rotamer_->atom_name( atom_other ), "gray", at2.xyz(), rad2, ring2 );

	if ( group ) {
		ostr << "@group { olap " << rotamer_->seqpos() << " " << rotamer_->atom_name( atom_this ) << "} dominant\n";
	}
	write_sphere_list_uv1( ostr, rotamer_->atom_name( atom_this ), "blue", at1.xyz(), rad1, hit_ring1 );

	if ( group ) {
		ostr << "@group { olap " << invdots_other.rotamer_->seqpos() << " " << invdots_other.rotamer_->atom_name( atom_other ) << "} dominant\n";
	}
	write_sphere_list_uv1( ostr,  invdots_other.rotamer_->atom_name( atom_other ), "blue", at2.xyz(), rad2, hit_ring2 );

}


bool
InvRotamerDots::overlap_exposed(
	core::conformation::Atom const & at1,
	utility::vector1< ObjexxFCL::ubyte > const & at1exposed_dots,
	core::conformation::Atom const & at2,
	utility::vector1< ObjexxFCL::ubyte > const & at2exposed_dots
) const
{

	Real const rad1 = (*radii_)[ at1.type() ] + RotamerDots::probe_radius_;
	Real const rad2 = (*radii_)[ at2.type() ] + RotamerDots::probe_radius_;

	Real const dist_sq = at1.xyz().distance_squared( at2.xyz() );
debug_assert( dist_sq <= (rad1 + rad2) * (rad1 + rad2) );
	Real const distance = std::sqrt( dist_sq );

	Real const step_size1 = rad1 * 0.02;
	Real const step_size2 = rad2 * 0.02;

	Size const nsteps1 = (Size)(ceil( max_dist_from_dot_to_intersection / step_size1 ));
	Size const nsteps2 = (Size)(ceil( max_dist_from_dot_to_intersection / step_size2 ));

	int degree_of_overlap1, degree_of_overlap2;
	int aphi_1_2, aphi_2_1;
	int theta_1_2, theta_2_1;

	//ronj this block represents the amount of surface area covered up on atom1 by atom2
	core::scoring::sasa::get_legrand_atomic_overlap( rad1, rad2, distance, degree_of_overlap1 );
	core::scoring::sasa::get_legrand_2way_orientation( at1.xyz(), at2.xyz(), aphi_1_2, theta_1_2, aphi_2_1, theta_2_1, distance );

	bool at1_intersection_exposed = false;
	Size closest_dot1 = (*RotamerDots::lg_angles_)( aphi_1_2, theta_1_2 );
	if ( degree_of_overlap1 + nsteps1 > 100 ) {
		for ( Size bb = 1; bb <= RotamerDots::num_bytes_; ++bb ){
			if ( at1exposed_dots[ bb ] ) {
				at1_intersection_exposed = true;
				break;
			}
		}
	} else {
		int masknum = ( closest_dot1 * 100 ) + degree_of_overlap1 + nsteps1;
		for ( Size bb = 1, bbli = (*RotamerDots::lg_masks_).index( bb, masknum ); bb <= RotamerDots::num_bytes_; ++bb, ++bbli ){
			if ( at1exposed_dots[ bb ] & (*RotamerDots::lg_masks_)[ bbli ] ) {
				at1_intersection_exposed = true;
				break;
			}
		}
	}

	if ( ! at1_intersection_exposed ) return false;

	//ronj the amount of surface area covered up on atom2 by atom1
	core::scoring::sasa::get_legrand_atomic_overlap( rad2, rad1, distance, degree_of_overlap2 );

	bool at2_intersection_exposed = false;

	Size closest_dot2 = (*RotamerDots::lg_angles_)( aphi_2_1, theta_2_1 );
	if ( degree_of_overlap2 + nsteps2 > 100 ) {
		for ( Size bb = 1; bb <= RotamerDots::num_bytes_; ++bb ){
			if ( at2exposed_dots[ bb ] ) {
				at2_intersection_exposed = true;
				break;
			}
		}
	} else {
		int masknum = ( closest_dot2 * 100 ) + degree_of_overlap2 + nsteps2;
		for ( Size bb = 1, bbli = (*RotamerDots::lg_masks_).index( bb, masknum ); bb <= RotamerDots::num_bytes_; ++bb, ++bbli ){
			if ( at2exposed_dots[ bb ] & (*RotamerDots::lg_masks_)[ bbli ] ) {
				at2_intersection_exposed = true;
				break;
			}
		}
	}

	return at2_intersection_exposed;
}

///
/// @brief
/// Helper method I am using to confirm that the dots are being overlapped and bits are being set correctly.
///
void InvRotamerDots::print_dot_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const {
	for ( Size bb = 1; bb <= RotamerDots::num_bytes_; ++bb ) {
		if ( (bb-1)*8 % 16 == 0 ) std::cout << (bb-1) * 8 << ":";
		for ( int index=7; index >= 0; index-- ) {
			int bit = ( ( (int)values[ bb ] >> index ) & 1 );
			std::cout << bit;
		}
		std::cout << " ";
	}
	std::cout << std::endl;
}


} // interaction_graph
} // pack
} // core


//  LocalWords:  endl
