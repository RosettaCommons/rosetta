// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNP_LowResPotential.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Kalli Kappel


// Unit headers
#include <core/scoring/rna/RNP_LowResPotential.hh>

// Package headers
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/kinematics/MinimizerMapBase.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <utility/io/izstream.hh>
#include <basic/database/open.hh>
#include <core/chemical/rna/util.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is a statistically derived low-resolution potential for RNA/protein interactions
// For RNA/protein modeling, this is meant to supplement the RNA low-res and protein low-res score
// functions
//
///////////////////////////////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "core.scoring.rna.RNP_LowResPotential" );
using namespace core::chemical;
using namespace core::chemical::rna;
using namespace basic::options;
using namespace basic::options::OptionKeys::score;

namespace core {
namespace scoring {
namespace rna {

/// c-tor
RNP_LowResPotential::RNP_LowResPotential():
	max_aa_( 20 ),
	max_base_( 4 ),
	num_xbins_( 10 ),
	num_ybins_( 10 ),
	num_dbins_( 10 ),
	num_backbone_dbins_( 7 )
{
	// I don't want all of these arrays initialized by default
	// only initialize them if there are nonzero wts for these score terms
	// (otherwise bound to have memory issues)
	// rnp_basepair_xy_.dimension( num_xbins_, num_ybins_, max_base_ /*4*/, max_aa_ /*20*/ );
	// rnp_stack_xy_.dimension( num_xbins_, num_ybins_, max_base_ /*4*/, max_aa_ /*20*/ );
	//rnp_pair_.dimension( num_dbins_, max_base_, max_aa_ );
	// rnp_pair_base_interface_protein_buried_.dimension( num_dbins_, max_base_, max_aa_ );
	// rnp_pair_base_interface_protein_notburied_.dimension( num_dbins_, max_base_, max_aa_ );

	// initialize the first time that evaluate is called
	initialize_rnp_base_pair();
	// initialize_rnp_pair();
	// initialize_rnp_stack_xy();
	initialize_rnp_aa_rna_backbone();
}

// Define copy constructor and clone?

/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResPotential::initialize_rnp_base_pair() {
	rnp_basepair_xy_.dimension( num_xbins_, num_ybins_, max_base_ /*4*/, max_aa_ /*20*/ );

	// use different file depending on whether we're scoring from actual centroid positions or "CEN"
	// positions
	std::string filename;
	if ( basic::options::option[ basic::options::OptionKeys::score::FA_low_res_rnp_scoring ]() ) { //default false
		filename = "scoring/rna/rnp_base_pair_actual_centroid.txt";
	} else {
		filename = "scoring/rna/rnp_base_pair.txt";
	}

	utility::io::izstream data_stream( basic::database::full_name( filename ) );

	if ( !data_stream ) {
		std::cerr << "Can't file specified RNP basepair potential file: " << filename << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
	}
	TR << "Reading RNP basepair x - y potential file: " << filename << std::endl;

	// read in the data
	// just copying what's been done in other potentials, look up base by number, AA by name
	chemical::AA aa;
	char base;
	Size base_num, xbin, ybin;
	Real potential;
	while ( data_stream >> base ) {
		base_num = core::chemical::rna::convert_acgu_to_1234( base );
		//data_stream >> aa >> xbin >> ybin >> potential >> skip;
		data_stream >> aa >> xbin >> ybin >> potential;
		rnp_basepair_xy_( xbin, ybin, base_num, aa ) = potential;
	}
	// std::string aa_str( "VAL" );
	// chemical::AA aa_val( chemical::aa_from_name( aa_str ) );
	// std::cout << "Some example values!" << std::endl;
	// std::cout << rnp_basepair_xy_(9,2,2,aa_val ) << std::endl;
	// std::cout << rnp_basepair_xy_(9,9,2,aa_val ) << std::endl;

	TR << "Finished reading RNP basepair x - y potential file: " << filename << std::endl;
	data_stream.close();

}


/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResPotential::initialize_rnp_aa_rna_backbone() {
	rnp_aa_rna_backbone_.dimension( num_backbone_dbins_ /*7*/, max_aa_ /*20*/ );

	std::string filename;
	if ( basic::options::option[ basic::options::OptionKeys::score::FA_low_res_rnp_scoring ]() ) { //default false
		filename = "scoring/rna/rnp_backbone_potential_actual_centroid.txt";
	} else {
		filename = "scoring/rna/rnp_backbone_potential.txt";
	}

	utility::io::izstream data_stream( basic::database::full_name( filename ) );

	if ( !data_stream ) {
		std::cerr << "Can't file specified RNP stack potential file: " << filename << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
	}
	TR << "Reading RNP AA RNA backbone potential file: " << filename << std::endl;

	// read in the data
	// just copying what's been done in other potentials, look up base by number, AA by name
	chemical::AA aa;
	Size dbin;
	// bins go from 3 - 10 (i.e. 3-4, 4-5, 5-6, 6-7, 7-8, 8-9, 9-10)
	Real potential;
	while ( data_stream >> aa ) {
		data_stream >> dbin >> potential;
		dbin -= 2; // b/c bins start at 3
		rnp_aa_rna_backbone_( dbin, aa ) = potential;
	}

	// for ( Size j =1; j <= max_aa_; ++j ) {
	//  for ( Size i = 1; i <= num_backbone_dbins_; ++i ) {
	//   std::cout << j << " " << i << " " << rnp_aa_rna_backbone_( i, j) << std::endl;
	//  }
	// }

	TR << "Finished reading RNP AA RNA backbone potential file: " << filename << std::endl;
	data_stream.close();

}


/////////////////////////////////////////////////////////////////////////////
//void
//RNP_LowResPotential::initialize_rnp_stack_xy() {
// rnp_stack_xy_.dimension( num_xbins_, num_ybins_, max_base_ /*4*/, max_aa_ /*20*/ );
//
// std::string const filename( "scoring/rna/rnp_stack.txt" );
// utility::io::izstream data_stream( basic::database::full_name( filename ) );
//
// if ( !data_stream ) {
//  std::cerr << "Can't file specified RNP stack potential file: " << filename << std::endl;
//  utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
// }
// TR << "Reading RNP stack potential file: " << filename << std::endl;
//
// // read in the data
// // just copying what's been done in other potentials, look up base by number, AA by name
// chemical::AA aa;
// char base;
// Size base_num, xbin, ybin;
// Real potential;
// while ( data_stream >> base ) {
//  base_num = core::chemical::rna::convert_acgu_to_1234( base );
//  //data_stream >> aa >> xbin >> ybin >> potential >> skip;
//  data_stream >> aa >> xbin >> ybin >> potential;
//  rnp_stack_xy_( xbin, ybin, base_num, aa ) = potential;
// }
//// std::string aa_str( "VAL" );
//// chemical::AA aa_val( chemical::aa_from_name( aa_str ) );
//// std::cout << "Some example values!" << std::endl;
//// std::cout << rnp_basepair_xy_(9,2,2,aa_val ) << std::endl;
//// std::cout << rnp_basepair_xy_(9,9,2,aa_val ) << std::endl;
//
// TR << "Finished reading RNP stack potential file: " << filename << std::endl;
// data_stream.close();
//
//
//}
/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResPotential::initialize_rnp_pair() {

	rnp_pair_base_interface_protein_buried_.dimension( num_dbins_, max_base_, max_aa_ );
	rnp_pair_base_interface_protein_notburied_.dimension( num_dbins_, max_base_, max_aa_ );

	std::string filename;
	if ( basic::options::option[ basic::options::OptionKeys::score::FA_low_res_rnp_scoring ]() ) { //default false
		filename = "scoring/rna/rnp_pair_protein_buried_actual_centroid.txt";
	} else {
		filename = "scoring/rna/rnp_pair_protein_buried.txt";
	}

	utility::io::izstream data_stream( basic::database::full_name( filename ) );

	if ( !data_stream ) {
		std::cerr << "Can't file RNP pair potential file: " << filename << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
	}
	TR << "Reading RNP pair potential file: " << filename << std::endl;

	// read in the data
	chemical::AA aa;
	char base;
	Size d_bin, base_num;
	Real potential;
	while ( data_stream >> base ) {
		base_num = core::chemical::rna::convert_acgu_to_1234( base );
		data_stream >> aa >> d_bin >> potential;
		rnp_pair_base_interface_protein_buried_( d_bin, base_num, aa ) = potential;
	}

	TR << "Finished reading RNP pair potential file: " << filename << std::endl;


	std::string filename2;
	if ( basic::options::option[ basic::options::OptionKeys::score::FA_low_res_rnp_scoring ]() ) { //default false
		filename2 = "scoring/rna/rnp_pair_protein_notburied_actual_centroid.txt";
	} else {
		filename2 = "scoring/rna/rnp_pair_protein_notburied.txt";
	}

	utility::io::izstream data_stream2( basic::database::full_name( filename2 ) );

	if ( !data_stream2 ) {
		std::cerr << "Can't file RNP pair potential file: " << filename2 << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
	}
	TR << "Reading RNP pair potential file: " << filename2 << std::endl;

	// read in the data
	while ( data_stream2 >> base ) {
		base_num = core::chemical::rna::convert_acgu_to_1234( base );
		data_stream >> aa >> d_bin >> potential;
		rnp_pair_base_interface_protein_notburied_( d_bin, base_num, aa ) = potential;
	}

	TR << "Finished reading RNP pair potential file: " << filename2 << std::endl;


}
/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResPotential::evaluate_rnp_base_pair_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const & x,
	Real const & y,
	Real & rnp_bp_score
) const {

	if ( rsd1.has_variant_type( REPLONLY ) || rsd2.has_variant_type( REPLONLY ) ) {
		rnp_bp_score = 0.0;
		return;
	}
	
	// Only evaluate these score terms between RNA and protein residues
	if ( !(( rsd1.is_RNA() && rsd2.is_protein() ) || ( rsd1.is_protein() && rsd2.is_RNA() )) ) {
		rnp_bp_score = 0.0;
		return;
	}

	chemical::AA aa;
	Size base_num;

	if ( rsd1.is_protein() ) {
		aa = rsd1.aa();
		base_num = core::chemical::rna::convert_acgu_to_1234( rsd2.name1() );
	} else {
		aa = rsd2.aa();
		base_num = core::chemical::rna::convert_acgu_to_1234( rsd1.name1() );
	}
	// probably need to do some interpolation, this is not smooth...
	Size xbin = Size( ( x + 10 ) / 2 ) + 1;
	Size ybin = Size( ( y + 10 ) / 2 ) + 1;

	if ( xbin < 1 ) xbin = 1;
	if ( ybin < 1 ) ybin = 1;
	if ( xbin > 10 ) xbin = 10;
	if ( ybin > 10 ) ybin = 10;

	rnp_bp_score = rnp_basepair_xy_( xbin, ybin, base_num, aa );

}
/////////////////////////////////////////////////////////////////////////////
//
//void
//RNP_LowResPotential::evaluate_rnp_stack_xy_score(
// conformation::Residue const & rsd1,
// conformation::Residue const & rsd2,
// Real const & x,
// Real const & y,
// Real & rnp_stack_score
//) const {
//
// // Only evaluate these score terms between RNA and protein residues
// if (!(( rsd1.is_RNA() && rsd2.is_protein() ) || ( rsd1.is_protein() && rsd2.is_RNA() ))) {
//  rnp_stack_score = 0.0;
//  return;
// }
//
// chemical::AA aa;
// Size base_num;
//
// if ( rsd1.is_protein() ) {
//  aa = rsd1.aa();
//  base_num = core::chemical::rna::convert_acgu_to_1234( rsd2.name1() );
// } else {
//  aa = rsd2.aa();
//  base_num = core::chemical::rna::convert_acgu_to_1234( rsd1.name1() );
// }
// // probably need to do some interpolation, this is not smooth...
// Size xbin = Size( ( x + 10 ) / 2 ) + 1;
// Size ybin = Size( ( y + 10 ) / 2 ) + 1;
//
// if ( xbin < 1 ) xbin = 1;
// if ( ybin < 1 ) ybin = 1;
// if ( xbin > 10 ) xbin = 10;
// if ( ybin > 10 ) ybin = 10;
//
// rnp_stack_score = rnp_stack_xy_( xbin, ybin, base_num, aa );
//
//}
/////////////////////////////////////////////////////////////////////////////
Real
RNP_LowResPotential::evaluate_rnp_aa_rna_backbone_score(
	conformation::Residue const & protein_rsd,
	Real const & dist_to_backbone
) const {

	chemical::AA aa;
	aa = protein_rsd.aa();

	// get the dbin
	Size d_bin = Size( dist_to_backbone ) + 1 - 2; // subtract 2 b/c bins span from 3-10
	if ( d_bin < 1 ) d_bin = 1;
	if ( d_bin > num_backbone_dbins_ ) d_bin = num_backbone_dbins_;

	return rnp_aa_rna_backbone_( d_bin, aa );

}


/////////////////////////////////////////////////////////////////////////////

void
RNP_LowResPotential::evaluate_rnp_pair_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool const & rsd1_is_interface,
	bool const & rsd1_is_buried,
	bool const & rsd2_is_interface,
	bool const & rsd2_is_buried,
	Real const & d,
	Real & rnp_pair_score
) const {

	// Right now just do this calculation if the base is at the interface
	// and the protein is at the interface
	// and the protein is either buried or not buried
	if ( !rsd1_is_interface || !rsd2_is_interface ) return; // both need to be at the interface

	chemical::AA aa;
	Size base_num;
	bool is_protein_buried;

	if ( rsd1.is_protein() ) {
		aa = rsd1.aa();
		base_num = core::chemical::rna::convert_acgu_to_1234( rsd2.name1() );
		is_protein_buried = rsd1_is_buried;
	} else {
		aa = rsd2.aa();
		base_num = core::chemical::rna::convert_acgu_to_1234( rsd1.name1() );
		is_protein_buried = rsd2_is_buried;
	}

	// get the dbin
	Size d_bin = Size( d ) + 1;
	if ( d_bin < 1 ) d_bin = 1;
	if ( d_bin > num_dbins_ ) d_bin = num_dbins_;

	if ( is_protein_buried ) {
		//if ( rsd1_is_buried ) {
		rnp_pair_score = rnp_pair_base_interface_protein_buried_( d_bin, base_num, aa );
	} else {
		rnp_pair_score = rnp_pair_base_interface_protein_notburied_( d_bin, base_num, aa );
	}

}

/////////////////////////////////////////////////////////////////////////////

} //rna
} //scoring
} //core
