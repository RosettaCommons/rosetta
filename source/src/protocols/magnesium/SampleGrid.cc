// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/SampleGrid.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/SampleGrid.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.magnesium.SampleGrid" );

using namespace core;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// SampleGrid allows predefinition of positions where metal ions will be sampled, based on neighborhoods
//  of hydrogen-bond acceptor atoms in pose.
//
// There's probably a nice object somewhere in Rosetta (e.g., map) for doing this.
// Anyway, FArray, with helper classes, will work for now.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace magnesium {

//Constructor
SampleGrid::SampleGrid( pose::Pose const & pose )
{
	figure_out_box_bounds( pose );
}

//Destructor
SampleGrid::~SampleGrid()
{}

///////////////////////////////////////////////////////////////////////////////
utility::vector1< Vector >
SampleGrid::get_mg_positions( pose::Pose const & pose ) {

	using namespace core::scoring::rna;

	scan_res_ = figure_out_scan_res( input_scan_res_, pose );

	create_grid();
	Size xgridsize = min_distance_grid_.size1();
	Size ygridsize = min_distance_grid_.size2();
	Size zgridsize = min_distance_grid_.size3();
	ObjexxFCL::FArray3D< Size > closest_res_grid( xgridsize, ygridsize, zgridsize );
	closest_res_grid = 0;

	// THIS IS RIDICULOUS.
	scoring::rna::RNA_ScoringInfo const & rna_scoring_info( scoring::rna::rna_scoring_info_from_pose( pose ) );
	utility::vector1< utility::vector1< Size > > const
		atom_numbers_for_mg_calculation( rna_scoring_info.atom_numbers_for_mg_calculation() );

	for ( Size q = 1; q <= scan_res_.size(); q++ ) {

		Size const n = scan_res_[ q ];
		if ( n > pose.total_residue() ) continue;

		utility::vector1< Size > const & atom_numbers1   ( atom_numbers_for_mg_calculation[ n ]  );
		for ( Size m = 1; m <= atom_numbers1.size(); ++m ) {

			core::conformation::Residue const & rsd1 = pose.residue( n );
			Size const i = atom_numbers1[ m ];
			if ( rsd1.is_virtual(i) ) continue;
			if ( !rsd1.heavyatom_is_an_acceptor(i) ) continue;

			Vector const i_xyz( rsd1.xyz(i) );

			Real const subgrid_radius = 3.0;
			Real const min_radius = 1.5;

			// defines the neighborhood around each acceptor in which to scan for Mg(2+) positions.
			Size xbinmin, xbinmax, ybinmin, ybinmax, zbinmin, zbinmax;
			define_bins( i_xyz.x(), subgrid_radius, xmin_, xgridsize, xyz_step_, xbinmin, xbinmax );
			define_bins( i_xyz.y(), subgrid_radius, ymin_, ygridsize, xyz_step_, ybinmin, ybinmax );
			define_bins( i_xyz.z(), subgrid_radius, zmin_, zgridsize, xyz_step_, zbinmin, zbinmax );

			for ( Size xbin = xbinmin; xbin <= xbinmax; xbin++ ) {
				for ( Size ybin = ybinmin; ybin <= ybinmax; ybin++ ) {
					for ( Size zbin = zbinmin; zbin <= zbinmax; zbin++ ) {

						Real const x = get_position( xbin, xmin_, xyz_step_ );
						Real const y = get_position( ybin, ymin_, xyz_step_ );
						Real const z = get_position( zbin, zmin_, xyz_step_ );
						Vector const mg_position = Vector(x,y,z);
						Real const r = ( mg_position - i_xyz ).length() ;
						if ( r > subgrid_radius ) continue;
						if ( r < min_radius ) continue;

						if ( closest_res_grid( xbin, ybin, zbin ) == 0 ||
								r < min_distance_grid_( xbin, ybin, zbin ) ) {
							closest_res_grid ( xbin, ybin, zbin ) = n;
							min_distance_grid_( xbin, ybin, zbin ) = r;
						}

					} // z
				} // y
			} // x
		} //  m = atom_numbers
	} // q = scan_res


	// Now go back through entire grid, and pick out bins that were within acceptable distances
	// to acceptors on desired residues.
	utility::vector1< Vector > mg_positions;
	for ( Size xbin = 1; xbin <= xgridsize; xbin++ ) {
		for ( Size ybin = 1; ybin <= ygridsize; ybin++ ) {
			for ( Size zbin = 1; zbin <= zgridsize; zbin++ ) {
				Size const closest_res = closest_res_grid( xbin, ybin, zbin );
				if ( closest_res == 0 ) continue;
				if ( tether_to_closest_res_ && !input_scan_res_.has_value( closest_res ) ) continue;
				Real const x = get_position( xbin, xmin_, xyz_step_ );
				Real const y = get_position( ybin, ymin_, xyz_step_ );
				Real const z = get_position( zbin, zmin_, xyz_step_ );
				mg_positions.push_back( Vector( x, y, z) );
			}
		}
	}

	return mg_positions;
}

///////////////////////////////////////////////////////////
utility::vector1< Size >
SampleGrid::figure_out_scan_res( utility::vector1< Size > const & input_scan_res,
	pose::Pose const & pose ) {

	utility::vector1< Size > scan_res = input_scan_res;

	if ( input_scan_res.size() == 0 ) {
		for ( Size n = 1; n <= pose.total_residue(); n++ ) scan_res.push_back( n );
	}

	if ( tether_to_closest_res_ ) {
		// expand neighborhood to scan, but later
		// only sample grid positions that have closest residue in input_res.
		static Distance NBR_DIST_CUTOFF( 12.0 );
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			if ( scan_res.has_value( n ) ) continue;
			for ( Size q = 1; q <= input_scan_res.size(); q++ ) {
				if ( ( pose.residue( input_scan_res[q] ).nbr_atom_xyz() -
						pose.residue( n ).nbr_atom_xyz() ).length() < NBR_DIST_CUTOFF ) {
					scan_res.push_back( n );
					break;
				}
			}
		}
	}

	return scan_res;
}

///////////////////////////////////////////////////////////////////////////////
void
SampleGrid::figure_out_box_bounds( pose::Pose const & pose )
{

	// determine bounds of scan.
	bool init( false );
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {

		if ( pose.residue( i ).is_virtual_residue() ) continue;
		if ( pose.residue( i ).name3() == " MG" ) continue;

		for ( Size j = 1; j <= pose.residue( i ).natoms(); j++ ) {
			Vector pos = pose.residue( i ).xyz( j );

			if ( !init ) {
				xmin_ = pos.x();
				xmax_ = pos.x();
				ymin_ = pos.y();
				ymax_ = pos.y();
				zmin_ = pos.z();
				zmax_ = pos.z();
				init = true;
			}

			if ( pos.x() < xmin_ ) xmin_ = pos.x();
			if ( pos.x() > xmax_ ) xmax_ = pos.x();
			if ( pos.y() < ymin_ ) ymin_ = pos.y();
			if ( pos.y() > ymax_ ) ymax_ = pos.y();
			if ( pos.z() < zmin_ ) zmin_ = pos.z();
			if ( pos.z() > zmax_ ) zmax_ = pos.z();

		}
	}

	// a little padding to be safe, and round to integer.
	xmin_ = int( xmin_ ) - 1.0;
	ymin_ = int( ymin_ ) - 1.0;
	zmin_ = int( zmin_ ) - 1.0;

	xmax_ = int( xmax_ ) + 1.0;
	ymax_ = int( ymax_ ) + 1.0;
	zmax_ = int( zmax_ ) + 1.0;

}

///////////////////////////////////////////////////////////////////////////////
void
SampleGrid::create_grid()
{
	Size xgridsize( ( xmax_ - xmin_ ) / xyz_step_ + 0.5 );
	Size ygridsize( ( ymax_ - ymin_ ) / xyz_step_ + 0.5 );
	Size zgridsize( ( zmax_ - zmin_ ) / xyz_step_ + 0.5 );
	min_distance_grid_.dimension( xgridsize, ygridsize, zgridsize );
	min_distance_grid_ = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
void
SampleGrid::define_bins( Real const x,
	Real const subgrid_radius,
	Real const xmin,
	Size const xgridsize,
	Real const xyz_increment,
	Size & xbinmin,
	Size & xbinmax ) const
{

	xbinmin = std::max( 1,              int( ( x - subgrid_radius - xmin )/xyz_increment ) + 1 ); // The +1 is because we are indexing by 1.
	xbinmax = std::min( int(xgridsize), int( ( x + subgrid_radius - xmin )/xyz_increment ) + 1 );
}

///////////////////////////////////////////////////////////////////////////////
Real
SampleGrid::get_position( Size const xbin, Real const xmin, Real const xyz_increment ) const
{
	return (  xmin + ( xbin - 0.5 ) * xyz_increment ); // The -0.5 is to get to the center of the bin.
}


} //magnesium
} //protocols
