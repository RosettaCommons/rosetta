// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/make_rot_lib/MakeRotLibJobInputter.hh
/// @brief  Implementation file for MakeRotLibOptionsData class
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// unit headers
#include <protocols/make_rot_lib/MakeRotLibOptionsData.hh>

// core headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

// basic headers
#include <basic/Tracer.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/exit.hh>

// c++ headers
#include <string>
#include <fstream>
#include <sstream>

namespace protocols {
namespace make_rot_lib {

static basic::Tracer TR( "protocols.make_rot_lib.MakeRotLibOptionsData" );

/// @brief Reads in options file. Parses it twice. Stores data.
MakeRotLibOptionsData::MakeRotLibOptionsData( std::string filename ) :
	n_chi_( 0 ), n_centroids_( 0 ), semirotameric_( false ), KbT_( 1.4 )
{

	// init stuff to hold lines
	utility::vector1< std::string > lines;
	std::string line;

	// check to see if options file exists
	if ( !utility::file::file_exists( filename.c_str() ) ) {
		utility_exit_with_message( "Cannot find options file"+filename );
	}

	// open file
	std::ifstream options( filename.c_str() );

	// read in each line, ignore comments
	while ( getline( options, line ) ) {
		//std::istringstream l( line );
		if ( line.size() < 1 || line[0] == '#' ) continue;
		lines.push_back( line );
	}

	// close options file
	options.close();

	// first pass - get AA name, number of chi angles, number of centroids
	for ( Size i( 1 ); i<= lines.size(); ++i ) {
		std::string const & line( lines[i] );
		std::istringstream l( line );
		std::string tag;

		l >> tag;

		// get number of chi
		if ( tag == "NUM_CHI" ) {
			l >> n_chi_;
		} else if ( tag == "NUM_BB" ) {
			// get number of BB
			l >> n_bb_;
		} else if ( tag == "CENTROID" ) {
			// get number of clusters
			++n_centroids_;
		} else if ( tag == "AA_NAME" ) {
			// get amino acid name
			l >> name_;
		} else if ( tag == "SEMIROTAMERIC" ) {
			semirotameric_ = true;
		} else if ( tag == "ROTAMERIC" ) {
			semirotameric_ = false;
		} else if ( tag == "TEMPERATURE" ) {
			l >> KbT_;
		}

	}
	//TR << "read first pass" << std::endl;
	// resize arrays based on inputs
	chi_ranges_.resize( n_chi_ );
	bb_ids_.resize( n_bb_ );
	bb_ranges_.resize( n_bb_ );

	utility::vector1< utility::vector1< core::Real > > rotwells_for_chi;
	Size const nrotchi = semirotameric_ ? n_chi_ - 1 : n_chi_;
	rotwells_for_chi.resize( nrotchi );
	bool rotwells_specified = false;

	//TR << "resized arrays " << std::endl;
	// second pass - get omg, phi, psi, eps torsion ranges, chi torsion ranges, and centroid data
	core::Size bb_i = 1;
	for ( core::Size i( 1 ); i <= lines.size(); ++i ) {
		std::string const & line( lines[ i ] );
		std::istringstream l( line );
		std::string tag;

		l >> tag;

		// get phi range
		if ( tag == "PHI_RANGE" ) {
			//l >> phi_range_.low >> phi_range_.high >> phi_range_.step;
			l >> bb_ranges_[bb_i].low >> bb_ranges_[bb_i].high >> bb_ranges_[bb_i].step;
			bb_ids_[bb_i] = 2;
			++bb_i;
		} else if ( tag == "PSI_RANGE" ) {
			// get psi range
			l >> bb_ranges_[bb_i].low >> bb_ranges_[bb_i].high >> bb_ranges_[bb_i].step;
			bb_ids_[bb_i] = 3;
			++bb_i;
		} else if ( tag == "BB_RANGE" ) {
			l >> bb_ranges_[bb_i].low >> bb_ranges_[bb_i].high >> bb_ranges_[bb_i].step >> bb_ids_[bb_i];
			++bb_i;
		} else if ( tag == "OMG_RANGE" ) {
			// get omg range
			l >> omg_range_.low >> omg_range_.high >> omg_range_.step;
		} else if ( tag == "EPS_RANGE" ) {
			// get eps range
			l >> eps_range_.low >> eps_range_.high >> eps_range_.step;
		} else if ( tag == "CHI_RANGE" ) {
			// get chi range
			core::Size chi_num( 0 );
			l >> chi_num;
			l >> chi_ranges_[ chi_num ].low >> chi_ranges_[ chi_num ].high >> chi_ranges_[ chi_num ].step;
		} else if ( tag == "CENTROID" ) {
			// get centroids
			CentroidRotNumVec temp_crnv;
			if ( semirotameric_ ) {
				temp_crnv.resize( n_chi_-1 );
				for ( core::Size i( 1 ); i <= n_chi_-1; ++i ) l >> temp_crnv[ i ].angle >> temp_crnv[ i ].rot_num;
			} else {
				temp_crnv.resize( n_chi_ );
				for ( core::Size i( 1 ); i <= n_chi_; ++i ) l >> temp_crnv[ i ].angle >> temp_crnv[ i ].rot_num;
			}
			centroid_data_.push_back( temp_crnv );
		} else if ( tag == "ROTWELLS" ) {
			rotwells_specified = true;
			Size chi_num, n_rotwells;
			l >> chi_num >> n_rotwells;
			TR << chi_num << " " << n_rotwells << std::endl;
			for ( Size ii = 1; ii <= n_rotwells; ++ii ) {
				Size rotwell_angle;
				l >> rotwell_angle;
				TR << rotwell_angle << " " << std::endl;
				rotwells_for_chi[ chi_num ].push_back( rotwell_angle );
			}
		}
	}

	if ( centroid_data_.size() > 0 && rotwells_for_chi[ 1 ].size() > 0 ) {
		utility_exit_with_message( "Warning: you specified both centroids and rotamer well combinations. We can't do both." );
	}

	//using the centroids designation (especially PRO/B3P) we want to skip this
	if ( rotwells_specified ) {
		// Okay, we need to generate all combinations of rotamer wells.
		// indices is a vector over all chi that says which rotwell we're working with
		// for each chi...
		utility::vector1< Size > indices;
		indices.resize( nrotchi + 1, 1 );
		indices[ nrotchi + 1 ] = 0;
		rotwells_for_chi.push_back( utility::vector1< core::Real >( 1, 0 ) );

		Size p = 1;
		while ( indices[ nrotchi + 1 ] == 0 ) {
			CentroidRotNumVec temp_crnv;
			temp_crnv.resize( nrotchi );
			TR << "Pushing back ";
			for ( Size i( 1 ); i <= nrotchi; ++i ) {
				temp_crnv[ i ].angle   = rotwells_for_chi[ i ][ indices[ i ] ];
				temp_crnv[ i ].rot_num = indices[ i ];
				TR << " " << temp_crnv[ i ].angle << " ";
			}
			centroid_data_.push_back( temp_crnv );
			TR << std::endl;
			indices[ 1 ] += 1;
			while ( indices[ p ] > rotwells_for_chi[ p ].size() ) {
				if ( p <= nrotchi ) indices[ p ] = 1;
				else                indices[ p ] = 0;
				p = p + 1;
				indices[ p ] += 1;
				if ( indices[ p ] <= rotwells_for_chi[ p ].size() ) p = 1;
			}
		}
	}
	TR << "Resulting centroids:" << std::endl;
	for ( auto const & crnv : centroid_data_ ) {
		for ( auto const & dat : crnv ) TR << dat.angle << " " << dat.rot_num << " ";
		TR << std::endl;
	}
	n_centroids_ = centroid_data_.size();


	//TR << "Second pass done "<<std::endl;
	// lastly get the polymer type
	using namespace core::chemical;
	ResidueTypeSetCOP RTS( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	ResidueType const & RT( RTS->name_map( name_ ) );

	if ( RT.is_protein() ) {
		polymer_type_ = PEPTIDE;
	} else if ( RT.is_peptoid() ) {
		polymer_type_ = PEPTOID;
	} else {
		utility_exit_with_message("The MakeRotLib protocol only works for peptides and peptoids.");
	}

}

} // namespace make_rot_lib
} // namespace protocols
