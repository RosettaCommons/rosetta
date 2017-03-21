// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/ActiveSiteGrid.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/downstream/ActiveSiteGrid.hh>

// Package headers
#include <protocols/match/BumpGrid.hh>

// Project headers
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

using namespace ObjexxFCL;

ActiveSiteGrid::~ActiveSiteGrid() {}

ActiveSiteGrid::ActiveSiteGrid() :
	bin_width_( 0.0 ),
	bb_( Vector( 0.0 ), Vector( 0.0 ) ),
	grid_( Bool3DGridOP( new Bool3DGrid ) ),
	reset_grid_bb_( false )
{
}

ActiveSiteGrid::ActiveSiteGrid( ActiveSiteGrid const & other ) :
	utility::pointer::ReferenceCount(),
	bin_width_( other.bin_width_ ),
	bb_( other.bb_ ),
	grid_( Bool3DGridOP( new Bool3DGrid( *other.grid_ ) )),
	reset_grid_bb_( other.reset_grid_bb_ )
{
}

ActiveSiteGrid &
ActiveSiteGrid::operator = ( ActiveSiteGrid const & rhs )
{
	if ( this != & rhs ) {
		bin_width_ = rhs.bin_width_;
		bb_ = rhs.bb_;
		grid_ = Bool3DGridOP( new Bool3DGrid( *rhs.grid_ ) );
		reset_grid_bb_ = rhs.reset_grid_bb_;
	}
	return *this;
}

void
ActiveSiteGrid::initialize_from_gridlig_file( std::string const & fname )
{
	std::string filename = fname;
	utility::io::izstream istr( filename.c_str() );

	if ( ! istr ) {
		utility_exit_with_message( "Gridlig file " + fname + " not found" );
	}

	std::string name, liggrid;
	istr >> name >> liggrid;
	runtime_assert( name == "NAME:" );
	runtime_assert( liggrid == "gridlig" );
	std::string base;
	Real xbase( 0.0 ), ybase( 0.0 ), zbase( 0.0 );
	istr >> base;
	runtime_assert( base == "BASE:" );
	istr >> xbase; runtime_assert( ! istr.bad() );
	istr >> ybase; runtime_assert( ! istr.bad() );
	istr >> zbase; runtime_assert( ! istr.bad() );

	std::string size;
	istr >> size;
	runtime_assert( size == "SIZE:" );
	Size xsize( 0 ), ysize( 0 ), zsize( 0 );
	istr >> xsize; runtime_assert( ! istr.bad() );
	istr >> ysize; runtime_assert( ! istr.bad() );
	istr >> zsize; runtime_assert( ! istr.bad() );

	runtime_assert( xsize != 0 );
	runtime_assert( ysize != 0 );
	runtime_assert( zsize != 0 );

	std::string length;
	istr >> length;
	runtime_assert( length == "LENGTH:");
	Real xwidth( 0.0 ), ywidth( 0.0 ), zwidth( 0.0 );

	istr >> xwidth; runtime_assert( ! istr.bad() );
	istr >> ywidth; runtime_assert( ! istr.bad() );
	istr >> zwidth; runtime_assert( ! istr.bad() );

	runtime_assert( xwidth != 0 );
	runtime_assert( ywidth != 0 );
	runtime_assert( zwidth != 0 );

	/// the widths must be the same for all three dimensions.
	runtime_assert( xwidth == ywidth );
	runtime_assert( xwidth == zwidth );

	Vector lower_corner( xbase, ybase, zbase );
	Vector upper_corner( lower_corner );
	upper_corner.x() += xwidth * xsize;
	upper_corner.y() += ywidth * ysize;
	upper_corner.z() += zwidth * zsize;

	bin_width_ = xwidth;
	bb_ = BoundingBox( lower_corner, upper_corner );
	grid_ = Bool3DGridOP( new Bool3DGrid );
	grid_->set_bin_width( bin_width_ );
	grid_->set_bounding_box( bb_ );
	reset_grid_bb_ = false;


	numeric::geometry::hashing::Bin3D bin; bin[ 1 ] = 0; bin[ 2 ] = 0; bin[ 3 ] = 0;
	std::string line;
	std::istringstream line_stream;
	Size occupied;

	Size linenumber = 4; // already read first three lines.

	getline(istr, line); /// clear the EOL from the last >>

	while ( istr ) {
		getline(istr, line);
		if ( ! istr ) break;
		++linenumber;

		if ( is_blank(line) ) {
			++bin[ 1 ];
			bin[ 2 ] = 0;
			bin[ 3 ] = 0;
			continue;
		}

		line_stream.clear();
		line_stream.str(line);
		line_stream.seekg( std::ios::beg );

		if ( bin[ 1 ] >= xsize ) {
			utility_exit_with_message( "Expected active site definition to end at line " + utility::to_string( linenumber ) + ".  X dim is out of range" );
		}

		if ( bin[ 2 ] >= ysize ) {
			utility_exit_with_message( "Expected blank line to separate y/z slices at line " + utility::to_string( linenumber ) + ". Y dim is out of range" );
		}

		for ( bin[ 3 ] = 0; bin[ 3 ] < zsize; ++bin[ 3 ] ) {
			line_stream >> occupied;
			if ( line_stream.bad() ) {
				std::cerr << "ERROR: liggrid file " << filename << " near line " << linenumber << std::endl;
				std::cerr << "ERROR: Expected to read " << zsize << " entries, but found only " << bin[ 3 ] << " entries" << std::endl;
				utility_exit_with_message( "Could not read active-site gridlig definition file" );
			}
			if ( occupied ) {
				grid_->set_value_for_bin( bin, true );
			}
		}
		++bin[ 2 ];
	}
}


/// @brief Set the bounding box for this grid
void
ActiveSiteGrid::set_bounding_box( BoundingBox const & bb )
{
	bb_ = bb;
	grid_->set_bounding_box( bb_ );
}

void
ActiveSiteGrid::set_bin_width( Real width )
{
	bin_width_ = width;
	grid_->set_bin_width( bin_width_ );
}

Bool3DGrid const &
ActiveSiteGrid::grid() const
{
	debug_assert( ! reset_grid_bb_ );
	return *grid_;
}


/// @brief Is a point in this grid active?  False for a point outside the bounding box.
bool
ActiveSiteGrid::occupied( Vector const & p ) const
{
	return grid_->occupied( p );
}

/// @brief Reset all the voxels to false
void
ActiveSiteGrid::clear()
{
	grid_->clear();
}

void
ActiveSiteGrid::enlargen_to_capture_volume_within_radius_of_residue(
	core::conformation::Residue const & res,
	Real radius
)
{
	if ( res.nheavyatoms() < 1 ) return;

	Vector lower( res.xyz( 1 )), upper( lower );
	for ( Size ii = 2; ii <= res.nheavyatoms(); ++ii ) {
		lower.min( res.xyz( ii ) );
		upper.max( res.xyz( ii ) );
	}
	swell_bounding_box( lower, upper, radius );
}

void
ActiveSiteGrid::enlargen_to_capture_volume_within_radius_of_sidechain(
	core::conformation::Residue const & res,
	Real radius
)
{
	if ( res.first_sidechain_atom() > res.nheavyatoms() ) return;

	Vector lower( res.xyz( res.first_sidechain_atom() )), upper( lower );
	for ( Size ii = res.first_sidechain_atom() + 1; ii <= res.nheavyatoms(); ++ii ) {
		lower.min( res.xyz( ii ) );
		upper.max( res.xyz( ii ) );
	}
	swell_bounding_box( lower, upper, radius );
}

void
ActiveSiteGrid::enlargen_to_capture_volume_within_radius_of_backbone(
	core::conformation::Residue const & res,
	Real radius
)
{
	if ( res.nheavyatoms() < 1 || res.first_sidechain_atom() == 1 ) return;

	Vector lower( res.xyz( 1 )), upper( lower );
	for ( Size ii = 2; ii < res.first_sidechain_atom(); ++ii ) {
		lower.min( res.xyz( ii ) );
		upper.max( res.xyz( ii ) );
	}
	swell_bounding_box( lower, upper, radius );
}

void
ActiveSiteGrid::or_within_radius_of_residue(
	core::conformation::Residue const & res,
	Real radius
)
{
	prep_grid();
	for ( Size ii = 1; ii <= res.nheavyatoms(); ++ii ) {
		grid_->or_by_sphere_liberal( res.xyz( ii ), radius );
	}
}

/// @brief Set all the voxels within a certain radius of the sidechain atoms to true.
void
ActiveSiteGrid::or_within_radius_of_sidechain(
	core::conformation::Residue const & res,
	Real radius
)
{
	for ( Size ii = res.first_sidechain_atom(); ii <= res.nheavyatoms(); ++ii ) {
		grid_->or_by_sphere_liberal( res.xyz( ii ), radius );
	}
	prep_grid();
}

/// @brief Set all the voxels within a certain radius of the backbone atoms to true.
void
ActiveSiteGrid::or_within_radius_of_backbone(
	core::conformation::Residue const & res,
	Real radius
)
{
	prep_grid();
	for ( Size ii = 1; ii < res.first_sidechain_atom(); ++ii ) {
		grid_->or_by_sphere_liberal( res.xyz( ii ), radius );
	}
}

void ActiveSiteGrid::initialize()
{
	prep_grid();
}


void
ActiveSiteGrid::swell_bounding_box(
	Vector lower,
	Vector upper,
	Real radius
)
{
	lower -= radius;
	upper += radius;

	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( lower( ii ) < bb_.lower()( ii ) ) {
			lower.min( bb_.lower() );
			bb_.set_lower( lower );
			reset_grid_bb_ = true;
		}
	}

	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( upper( ii ) > bb_.upper()( ii ) ) {
			upper.max( bb_.upper() );
			bb_.set_upper( upper );
			reset_grid_bb_ = true;
		}
	}
}

void
ActiveSiteGrid::prep_grid()
{
	if ( ! reset_grid_bb_ ) return;

	grid_ = Bool3DGridOP( new Bool3DGrid );
	grid_->set_bin_width( bin_width_ );

	Vector snapped_lower( bb_.lower() );
	for ( Size ii = 1; ii <= 3; ++ii ) {
		/// align the lower boundary to an integer
		snapped_lower( ii ) = static_cast< Real > ( static_cast< int > ( snapped_lower( ii ) ));
	}
	bb_.set_lower( snapped_lower );

	grid_->set_bounding_box( bb_ );
	reset_grid_bb_ = false;

}

}
}
}
