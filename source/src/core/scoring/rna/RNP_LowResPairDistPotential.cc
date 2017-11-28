// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNP_LowResPairDistPotential.cc
/// @brief  Statistically derived RNP low-res potential, simple distance-based
/// @author Kalli Kappel


// Unit headers
#include <core/scoring/rna/RNP_LowResPairDistPotential.hh>

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

static basic::Tracer TR( "core.scoring.rna.RNP_LowResPairDistPotential" );
using namespace core::chemical;
using namespace core::chemical::rna;
using namespace basic::options;
using namespace basic::options::OptionKeys::score;

namespace core {
namespace scoring {
namespace rna {

/// c-tor
RNP_LowResPairDistPotential::RNP_LowResPairDistPotential():
	max_aa_( 25 ),
	max_base_( 8 ),
	num_dbins_( 10 )
{

	rnp_pair_dist_potential_.dimension( num_dbins_, max_base_, max_aa_ );
	rnp_pair_dist_potential_ = 0.0; // zero everything out

	initialize_rnp_pair_dist();

}

// Define copy constructor and clone?

/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResPairDistPotential::initialize_rnp_pair_dist() {
	//rnp_pair_dist_potential_.dimension( num_dbins_, max_base_, max_aa_ );
	//rnp_pair_dist_potential_ = 0.0; // zero everything out

	use_actual_centroid_ = basic::options::option[ basic::options::OptionKeys::score::rna::FA_low_res_rnp_scoring ]();

	std::string filename;
	if ( use_actual_centroid_ ) { //default false
		filename = "scoring/rna/rnp_pair_dist_potential_renorm_actual_centroid.txt";
	} else {
		filename = "scoring/rna/rnp_pair_dist_potential_renorm.txt";
	}

	utility::io::izstream data_stream( basic::database::full_name( filename ) );

	if ( !data_stream ) {
		std::cerr << "Can't file specified RNP pair dist file: " << filename << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
	}
	TR << "Reading RNP pair distance potential file: " << filename << std::endl;

	// read in the data
	// just copying what's been done in other potentials, look up base by number, AA by name
	Size aa, base;
	Size dbin;
	// bins go from 3 - 10 (i.e. 3-4, 4-5, 5-6, 6-7, 7-8, 8-9, 9-10)
	Real potential;
	while ( data_stream >> dbin ) {
		data_stream >> base >> aa >> potential;
		dbin += 1; // numbering of bins starts at 0 in the data file
		//  dbin -= 2; // b/c bins start at 3
		rnp_pair_dist_potential_( dbin, base, aa ) = potential;
	}

	// for ( Size j =1; j <= max_aa_; ++j ) {
	//  for ( Size i = 1; i <= num_backbone_dbins_; ++i ) {
	//   std::cout << j << " " << i << " " << rnp_aa_rna_backbone_( i, j) << std::endl;
	//  }
	// }

	TR << "Finished reading RNP pair distance potential file: " << filename << std::endl;
	data_stream.close();

	// Initialize vectors of RNA and protein atoms
	// Ordered list of RNA atoms
	RNA_atoms_.push_back( " P  "); // 1
	RNA_atoms_.push_back( " C5'"); // 2
	RNA_atoms_.push_back( " C1'"); // 3
	RNA_atoms_.push_back( " C3'"); // 4
	RNA_atoms_.push_back( "RCEN"); // get the base centroid

	protein_atoms_.push_back(" CA "); // 1
	protein_atoms_.push_back(" C  "); // 2
	protein_atoms_.push_back(" O  "); // 3
	protein_atoms_.push_back(" N  "); // 4
	protein_atoms_.push_back(" CB "); // 5
	protein_atoms_.push_back("CEN");


}


/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResPairDistPotential::evaluate_rnp_pair_dist_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real & rnp_pair_dist_score
) const {

	if ( rsd1.has_variant_type( REPLONLY ) || rsd2.has_variant_type( REPLONLY ) ) {
		rnp_pair_dist_score = 0.0;
		return;
	}

	// Only evaluate these score terms between RNA and protein residues
	if ( !(( rsd1.is_RNA() && rsd2.is_protein() ) || ( rsd1.is_protein() && rsd2.is_RNA() )) ) {
		rnp_pair_dist_score = 0.0;
		return;
	}


	if ( rsd1.is_protein() ) {
		rnp_pair_dist_score = calc_rnp_pair_dist_score( rsd1 /*protein residue*/, rsd2 /*RNA residue*/ );
	} else {
		rnp_pair_dist_score = calc_rnp_pair_dist_score( rsd2 /*protein residue*/, rsd1 /*RNA residue*/ );
	}

}

/////////////////////////////////////////////////////////////////////////////
Real
RNP_LowResPairDistPotential::calc_rnp_pair_dist_score(
	conformation::Residue const & protein_rsd,
	conformation::Residue const & RNA_rsd
) const {

	Real total_rnp_pair_dist_score = 0.0;

	bool is_centroid = protein_rsd.type().mode() == core::chemical::CENTROID_t;
	if ( !is_centroid && !use_actual_centroid_ ) {
		//tr.Warning << "rnp vdw energy not computed b/c protein is not centroid" << std::endl;
		return total_rnp_pair_dist_score;
	}

	// Loop through the RNA and protein atoms
	// each atom is going to have a specified index
	for ( core::Size R_ind = 1; R_ind <= RNA_atoms_.size(); ++R_ind ) {
		// figure out the RNA atom index for the potential
		Size R_ind_potential;
		if ( RNA_atoms_[R_ind] != "RCEN" ) R_ind_potential = R_ind;
		else R_ind_potential = core::chemical::rna::convert_acgu_to_1234( RNA_rsd.name1() ) + (R_ind - 1);

		for ( core::Size P_ind = 1; P_ind <= protein_atoms_.size(); ++P_ind ) {
			// figure out the protein atom index for the potential
			if ( (protein_atoms_[P_ind] == " CB ") && protein_rsd.name1() == 'G' ) continue;
			Size P_ind_potential;
			if ( protein_atoms_[P_ind] != "CEN" ) P_ind_potential = P_ind;
			else  P_ind_potential = convert_aa_to_index( protein_rsd.name1() ) + (P_ind -1);

			Real dist;
			Vector protein_xyz;
			if ( use_actual_centroid_ ) {
				protein_xyz = protein_rsd.actcoord();
			} else {
				protein_xyz = protein_rsd.xyz( protein_atoms_[P_ind] );
			}
			if ( RNA_atoms_[R_ind] != "RCEN" ) {
				// get the distance
				Vector diff_xyz = RNA_rsd.xyz( RNA_atoms_[R_ind] ) - protein_xyz;
				dist = diff_xyz.length();
			} else { // Use base centroid
				Vector RNA_base_centroid = chemical::rna::get_rna_base_centroid( RNA_rsd );
				Vector diff_xyz = RNA_base_centroid - protein_xyz;
				dist = diff_xyz.length();
			}
			// which bin does this correspond to?
			Size dbin = Size ( dist/2. ) + 1;
			if ( dbin > num_dbins_ ) dbin = num_dbins_;
			// get the potential
			// probably need to do some interpolation, this is not smooth...
			//   std::cout << "Distance: " << dist;
			//   std::cout << " dbin: " << dbin;
			//   std::cout << " RNA atom " << RNA_atoms_[R_ind];
			//   std::cout << " Protein rsd/atom " << protein_rsd.name1() << "/" << protein_atoms_[P_ind];
			//   std::cout << " RNA index " << R_ind_potential;
			//   std::cout << " Protein index " << P_ind_potential;
			//   std::cout << " SCORE: " << rnp_pair_dist_potential_( dbin, R_ind_potential, P_ind_potential ) << std::endl;;
			total_rnp_pair_dist_score += rnp_pair_dist_potential_( dbin, R_ind_potential, P_ind_potential );
		}
	}

	return total_rnp_pair_dist_score;

}
/////////////////////////////////////////////////////////////////////////////

Size
RNP_LowResPairDistPotential::convert_aa_to_index( char const c ) const
{
	// just copying what was done for RNA, but I'm sure
	// there's a much nicer way of doing this
	// also similar to another index conversion in util, should unify!
	if ( c == 'A' ) return 1;
	if ( c == 'R' ) return 2;
	if ( c == 'N' ) return 3;
	if ( c == 'D' ) return 4;
	if ( c == 'C' ) return 5;
	if ( c == 'E' ) return 6;
	if ( c == 'Q' ) return 7;
	if ( c == 'G' ) return 8;
	if ( c == 'H' ) return 9;
	if ( c == 'I' ) return 10;
	if ( c == 'L' ) return 11;
	if ( c == 'K' ) return 12;
	if ( c == 'M' ) return 13;
	if ( c == 'F' ) return 14;
	if ( c == 'P' ) return 15;
	if ( c == 'S' ) return 16;
	if ( c == 'T' ) return 17;
	if ( c == 'W' ) return 18;
	if ( c == 'Y' ) return 19;
	if ( c == 'V' ) return 20;
	TR << "What is this? " << c << std::endl;
	utility_exit_with_message( "Asked for protein aa index for unknown residue_name" );
	return 0;


}

} //rna
} //scoring
} //core
