// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNP_LowResStackData.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Kalli Kappel


// Unit headers
#include <core/scoring/rna/RNP_LowResStackData.hh>

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

static THREAD_LOCAL basic::Tracer TR( "core.scoring.rna.RNP_LowResStackData" );
using namespace core::chemical;
using namespace core::chemical::rna;
using namespace basic::options;
using namespace basic::options::OptionKeys::score;

namespace core {
namespace scoring {
namespace rna {

/// c-tor
RNP_LowResStackData::RNP_LowResStackData():
	max_aa_( 20 ),
	max_base_( 4 ),
	num_xbins_( 10 ),
	num_ybins_( 10 )
{

	// initialize the first time that evaluate is called
	initialize_rnp_stack_xy();
}

// Define copy constructor and clone?

/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResStackData::initialize_rnp_stack_xy() {
	rnp_stack_xy_.dimension( num_xbins_, num_ybins_, max_base_ /*4*/, max_aa_ /*20*/ );

	std::string const filename( "scoring/rna/rnp_stack.txt" );
	utility::io::izstream data_stream( basic::database::full_name( filename ) );

	if ( !data_stream ) {
		std::cerr << "Can't file specified RNP stack potential file: " << filename << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
	}
	TR << "Reading RNP stack potential file: " << filename << std::endl;

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
		rnp_stack_xy_( xbin, ybin, base_num, aa ) = potential;
	}
	// std::string aa_str( "VAL" );
	// chemical::AA aa_val( chemical::aa_from_name( aa_str ) );
	// std::cout << "Some example values!" << std::endl;
	// std::cout << rnp_basepair_xy_(9,2,2,aa_val ) << std::endl;
	// std::cout << rnp_basepair_xy_(9,9,2,aa_val ) << std::endl;

	TR << "Finished reading RNP stack potential file: " << filename << std::endl;
	data_stream.close();


}
/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResStackData::evaluate_rnp_stack_xy_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const & x,
	Real const & y,
	Real & rnp_stack_score
) const {

	if ( rsd1.has_variant_type( REPLONLY ) || rsd2.has_variant_type( REPLONLY ) ) {
		rnp_stack_score = 0.0;
		return;
	}

	// Only evaluate these score terms between RNA and protein residues
	if ( !(( rsd1.is_RNA() && rsd2.is_protein() ) || ( rsd1.is_protein() && rsd2.is_RNA() )) ) {
		rnp_stack_score = 0.0;
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

	rnp_stack_score = rnp_stack_xy_( xbin, ybin, base_num, aa );

}
/////////////////////////////////////////////////////////////////////////////

} //rna
} //scoring
} //core
