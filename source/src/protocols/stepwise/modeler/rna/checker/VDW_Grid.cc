// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/checker/VDW_Grid.cc
/// @brief
/// @detailed
/// @author Caleb Geniesse, geniesse@stanford.edu


#include <protocols/stepwise/modeler/rna/checker/VDW_Grid.hh>
#include <utility/vector1.hh>


#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.modeler.rna.checker.VDW_Grid" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {


// Constructor
VDW_Grid::VDW_Grid():
	is_occupied_( false )
{
	bins_.clear();
}


// Copy Constructor
VDW_Grid::VDW_Grid( VDW_Grid const & src ):
	bins_( src.bins_ ),
	is_occupied_( src.is_occupied_ )
{
}


//Destructor
VDW_Grid::~VDW_Grid()
{}


void
VDW_Grid::setup( int const & bin_max ) const
{
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

} //checker
} //rna
} //modeler
} //stepwise
} //protocols
