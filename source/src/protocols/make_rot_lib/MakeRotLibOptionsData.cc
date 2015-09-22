// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

static THREAD_LOCAL basic::Tracer TR( "protocols.make_rot_lib.MakeRotLibOptionsData" );

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
		}
	}
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
