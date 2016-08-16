// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/VDW_Grid.cc
/// @brief
/// @details
/// @author Caleb Geniesse, geniesse@stanford.edu


#include <core/pose/rna/VDW_Grid.hh>
#include <utility/vector1.hh>


#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.rna.VDW_Grid" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace rna {


// Constructor
VDW_Grid::VDW_Grid():
	is_occupied_( false ),
	bin_max_( 0 ),
	atom_bin_size_( 0 ),
	bin_offset_( 0 ),
	ref_xyz_( 0 )
{
	bins_.clear();
}

// Copy Constructor
VDW_Grid::VDW_Grid( VDW_Grid const & src ):
	bins_( src.bins_ ),
	is_occupied_( src.is_occupied_ ),
	bin_max_( src.bin_max_ ),
	atom_bin_size_( src.atom_bin_size_ ),
	bin_offset_( src.bin_offset_ ),
	ref_xyz_( src.ref_xyz_ )
{
}


//Destructor
VDW_Grid::~VDW_Grid()
{}


void
VDW_Grid::setup( int const & bin_max ) const
{
	bin_max_ = bin_max;
	utility::vector1< bool > one_dim_bin( bin_max*2, false );
	utility::vector1< utility::vector1< bool > > two_dim_bin( bin_max*2, one_dim_bin );
	bins_.assign( bin_max*2, two_dim_bin );
}


void
VDW_Grid::reset() const
{
	for ( core::Size x = 1; x <= bins_.size(); x++ ) {
		for ( core::Size y = 1; y <= bins_[x].size(); y++ ) {
			for ( core::Size z = 1; z <= bins_[x][y].size(); z++ ) {
				set_bin( x, y, z, false );
			}
		}
	}
	is_occupied_ = false;
}


void
VDW_Grid::reset( utility::vector1< Atom_Bin > const & occupied_xyz_bins_ ) const
{
	for ( core::Size n = 1; n <= occupied_xyz_bins_.size(); n++ ) {
		set_xyz_bin( occupied_xyz_bins_[ n ], false );
	}
	is_occupied_ = false;
}


core::Size
VDW_Grid::size() const
{
	//return bins_.size();
	return bins_.size() * bins_.size() * bins_.size();
}


bool
VDW_Grid::get_bin( int const & x, int const & y, int const & z ) const
{
	return bins_[x][y][z];
}


bool
VDW_Grid::get_xyz_bin( Atom_Bin const & xyz_bin ) const
{
	return get_bin( xyz_bin.x, xyz_bin.y, xyz_bin.z );
}


void
VDW_Grid::set_bin( int const & x, int const & y, int const & z, bool const & value ) const
{
	bins_[x][y][z] = value;
	if ( value && !is_occupied_ ) is_occupied_ = value;
}


void
VDW_Grid::set_xyz_bin( Atom_Bin const & xyz_bin, bool const & value ) const
{
	set_bin( xyz_bin.x, xyz_bin.y, xyz_bin.z, value );
}

bool
VDW_Grid::is_occupied() const
{
	return is_occupied_;
}

void
VDW_Grid::set_bin_max( int const & value ) const {
	bin_max_ = value;
}

int
VDW_Grid::get_bin_max() const {
	return bin_max_;
}

void
VDW_Grid::set_atom_bin_size( core::Real const & value ) const {
	atom_bin_size_ = value;
}

core::Real
VDW_Grid::get_atom_bin_size() const {
	return atom_bin_size_;
}

void
VDW_Grid::set_bin_offset( int const & value ) const {
	bin_offset_ = value;
}

int
VDW_Grid::get_bin_offset() const {
	return bin_offset_;
}

void
VDW_Grid::set_ref_xyz( numeric::xyzVector< core::Real > const & value ) const {
	ref_xyz_ = value;
}

numeric::xyzVector< core::Real >
VDW_Grid::get_ref_xyz() const {
	return ref_xyz_;
}


} //rna
} //pose
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::rna::VDW_Grid::save( Archive & arc ) const {
	arc( CEREAL_NVP( bins_ ) ); // utility::vector1<utility::vector1<utility::vector1<_Bool> > >
	arc( CEREAL_NVP( is_occupied_ ) ); // _Bool
	arc( CEREAL_NVP( bin_max_ ) ); // int
	arc( CEREAL_NVP( atom_bin_size_ ) ); // core::Real
	arc( CEREAL_NVP( bin_offset_ ) ); // int
	arc( CEREAL_NVP( ref_xyz_ ) ); // numeric::xyzVector<core::Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::rna::VDW_Grid::load( Archive & arc ) {
	arc( bins_ ); // utility::vector1<utility::vector1<utility::vector1<_Bool> > >
	arc( is_occupied_ ); // _Bool
	arc( bin_max_ ); // int
	arc( atom_bin_size_ ); // core::Real
	arc( bin_offset_ ); // int
	arc( ref_xyz_ ); // numeric::xyzVector<core::Real>
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::rna::VDW_Grid );
CEREAL_REGISTER_TYPE( core::pose::rna::VDW_Grid )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_rna_VDW_Grid )
#endif // SERIALIZATION
