// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/PairEPotential.cc
/// @brief  pairE knowledge-based potential class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Original function authors of attributed functions noted below
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/PairEPotential.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>

// ObjexxFCL headers

// Utility headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>

// Numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/interpolation/interpolation.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>

// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>




namespace core {
namespace scoring {

PairEPotential::PairEPotential() :
	pair_score_min_sep_( 1 ), // from pdbstatistics_pack
	pair_score_cb_thresh_( 16 ),
	pair_score_bin_range_( 1.5 ),
	pair_score_bin_base_( 3.0 ),
	max_bin_( 3 ) // APL from 3 to 2; defines range from 0 A to 6 A instead of all the way out to 7.5 A.
{
	//initialize all the data
	// Constants
	int const max_aa( 20 ); // Only set up for CAAs: NCAAs will share values via "cat_key()" lookup

	// Open the residue pair statistics file
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/PairEPotential/pdb_pair_stats_fine" );

	// Dimension/allocate the pair_corr array
	pair_corr_.dimension( max_aa, max_aa, 2, 2, 5, TableProbability( 0.0 ) );

	// Read the file and assign the array elements
	int aa1, aa2, e1, e2, r12_bin;
	TableProbability pair_probability;
	while ( stream ) {
		using namespace ObjexxFCL::format;
		stream
			>> bite( 2, aa1 ) >> skip( 1 )
			>> bite( 2, aa2 ) >> skip( 1 )
			>> bite( 2, e1 ) >> skip( 1 )
			>> bite( 2, e2 ) >> skip( 1 )
			>> bite( 2, r12_bin ) >> skip( 1 )
			>> bite( 12, pair_probability ) >> skip;

		//iwd  A few entries have 0 observations, which leads to a likelihood ratio of 0,
		//iwd  which causes INF and NAN values when you take the log of it.
		//iwd  It's also unrealistic, as the value is just due to statistics of small numbers.
		//iwd  This (conservative) replacement is the minimum (non-zero) probability in the file.
		TableProbability const min_prob = 0.05;
		if( pair_probability < min_prob ) pair_probability = min_prob;

		if ( stream ) {
			pair_corr_( aa1, aa2, e1, e2, r12_bin ) = pair_probability;
		}
	}
	stream.close();


	//bk One interesting property of the pair statistics are that they favor unlike charges being
	//bk near each other, but only very close range interactions between like charges are
	//bk disfavored.  This is presumably because like charged residues frequently come togethor to
	//bk bind charged ligands.  When performing protein design in the absence of extra ligands, it is
	//bk therefore probably more correct to penalize like charges being near each other.
	//bk The following section makes it unfavorable to have like charged residues near each other
	//bk if the -use_electrostic_repulsion flag is turned on.  The penalty is set to be roughly equal
	//bk but opposite to the favorable energy given to unlike charges.
	{ // Electrostatic repulsion
		//using namespace core::conformation::amino::AminoAcidKeys;
		using namespace core::chemical;//conformation;
		using namespace basic::options;
		using namespace OptionKeys::packing;

		if ( option[ use_electrostatic_repulsion ] ) {
			for ( int e1 = 1; e1 <= 2; ++e1 ) {
				for ( int e2 = 1; e2 <= 2; ++e2 ) {
					pair_corr_( aa_asp, aa_asp, e1, e2, 1 ) = 0.3;   // 3-4.5 angstroms
					pair_corr_( aa_asp, aa_asp, e1, e2, 2 ) = 0.5;   // 4.5-6 angstroms
					pair_corr_( aa_asp, aa_asp, e1, e2, 3 ) = 0.75;  // 6-7.5 angstroms

					pair_corr_( aa_asp, aa_glu, e1, e2, 1 ) = 0.3;
					pair_corr_( aa_asp, aa_glu, e1, e2, 2 ) = 0.5;
					pair_corr_( aa_asp, aa_glu, e1, e2, 3 ) = 0.75;

					pair_corr_( aa_glu, aa_asp, e1, e2, 1 ) = 0.3;
					pair_corr_( aa_glu, aa_asp, e1, e2, 2 ) = 0.5;
					pair_corr_( aa_glu, aa_asp, e1, e2, 3 ) = 0.75;

					pair_corr_( aa_glu, aa_glu, e1, e2, 1 ) = 0.3;
					pair_corr_( aa_glu, aa_glu, e1, e2, 2 ) = 0.5;
					pair_corr_( aa_glu, aa_glu, e1, e2, 3 ) = 0.75;

					pair_corr_( aa_lys, aa_lys, e1, e2, 1 ) = 0.3;
					pair_corr_( aa_lys, aa_lys, e1, e2, 2 ) = 0.5;
					pair_corr_( aa_lys, aa_lys, e1, e2, 3 ) = 0.75;

					pair_corr_( aa_arg, aa_arg, e1, e2, 1 ) = 0.3;
					pair_corr_( aa_arg, aa_arg, e1, e2, 2 ) = 0.5;
					pair_corr_( aa_arg, aa_arg, e1, e2, 3 ) = 0.75;

					pair_corr_( aa_arg, aa_lys, e1, e2, 1 ) = 0.3;
					pair_corr_( aa_arg, aa_lys, e1, e2, 2 ) = 0.5;
					pair_corr_( aa_arg, aa_lys, e1, e2, 3 ) = 0.75;

					pair_corr_( aa_lys, aa_arg, e1, e2, 1 ) = 0.3;
					pair_corr_( aa_lys, aa_arg, e1, e2, 2 ) = 0.5;
					pair_corr_( aa_lys, aa_arg, e1, e2, 3 ) = 0.75;
				}
			}
		}
	}
}

bool
PairEPotential::pair_term_energy_exists( conformation::Residue const & rsd ) const
{
	return ( (rsd.is_polar() || rsd.is_aromatic() ) && rsd.is_protein() );
}

Energy
PairEPotential::pair_term_energy(
	conformation::Residue const & res1,
	int res1_num_10A_neighbors,
	conformation::Residue const & res2,
	int res2_num_10A_neighbors
) const
{
	Probability temp1, temp2, temp3;
	return pair_term_energy( res1, res1_num_10A_neighbors, res2, res2_num_10A_neighbors, temp1, temp2, temp3 );
}

Energy
PairEPotential::pair_term_energy(
	conformation::Residue const & res1,
	int res1_num_10A_neighbors,
	conformation::Residue const & res2,
	int res2_num_10A_neighbors,
	Probability & pair_lhood_ratio,
	Probability & pair_lhood_ratio_high,
	Probability & pair_lhood_ratio_low
) const
{
	//using namespace pdbstatistics_pack; // Various constants and tables
	//using namespace pdbstatistics_pack::pdbstatistics; // Energy lookup tables
	using namespace basic::options;
	using namespace core::chemical; //conformation;
	using namespace basic::options::OptionKeys;
	//using namespace core::conformation::amino::AminoAcidKeys;
	using numeric::abs_difference;
	using numeric::interpolation::interpolated;

	assert( res1.seqpos() != res2.seqpos() ); // Only call for distinct residues
	assert( res1.is_polar() || res1.is_aromatic() );
	assert( res2.is_polar() || res2.is_aromatic() ); // Only for polar amino acids: Caller does exclusion (prevents call overhead)
	//assert( pair_corr_.I1() == pair_corr_.I2() ); // First 2 index ranges should match //this is a silly assert

//	if ( !is_protein(aa1)  ||  !is_protein(aa2) ) ) return; //dr okay for dupes but not other nnaa where we won't have this info     //! Change to an NCAA exclusion????

	//jk option to suppress computing pair term for histidine (numbers are skewed due to metal-binding sites)
	if ( ( option[ corrections::score::no_his_his_pairE ] ) &&
		( ( res1.aa() == aa_his ) &&
		( res2.aa() == aa_his ) ) )
	{
		return Energy( 0.0 );
	}

	if ( ( option[ corrections::score::no_his_DE_pairE ] ) &&
		((( res1.aa() == aa_his ) &&
		( res2.aa() == aa_asp || res2.aa() == aa_glu ) ) ||
		( ( res1.aa() == aa_asp || res1.aa() == aa_glu ) &&
		( res2.aa() == aa_his ))) )
	{
		return Energy( 0.0 );
	}



	if ( pair_score_min_sep_ > 1 ) { // Short-circuit for speed
		if ( res1.polymeric_sequence_distance( res2 ) < pair_score_min_sep_ ) return Energy( 0.0 );
	}

	// Amino acid indexes for lookup
	int const pair_corr_n_AA( pair_corr_.u1() );
	AA const aa1n( res1.aa() );
	if ( aa1n > pair_corr_n_AA ) return Energy( 0.0 ); // Unsupported amino acid
	AA const aa2n( res2.aa() );
	if ( aa2n > pair_corr_n_AA ) return Energy( 0.0 ); // Unsupported amino acid

	// Number of neighbor bins
	int const e1( ( res1_num_10A_neighbors <= pair_score_cb_thresh_ ) ? 1 : 2 );
	int const e2( ( res2_num_10A_neighbors <= pair_score_cb_thresh_ ) ? 1 : 2 );

	// Action center distance
	Distance const r12( res1.actcoord().distance( res2.actcoord() ) );
	// If r12 == 0, one of two things has happened:
	// (1) You're comparing a residue to itself -- don't.
	// (2) The residues actcoords aren't being updated and are (0,0,0)
	// -- see ResidueType.requires_actcoord() and the appropriate .params files.
	assert( r12 > 0 );

	// 1. Get the lower of the two bin averages bracketing r12 and set r12_bin to this value.
	//    Consider 'bin average' to exist halfway through bin's range.
	//
	// 2. Find the difference between the low_bin average and the actual value of r12 in bin units
	//
	//int const max_bin( 3 ); //apl should this be hard coded here?

	Distance const r12_bin_real( std::max(
	 ( r12 / pair_score_bin_range_ ) + 1 - pair_score_bin_base_,
	 Distance( 0.5 ) ) ); // First bin is [ 3 A, 4.5 A ] but we use it to represent r12 down to 0 A
	int const r12_bin( std::max( numeric::nint( r12_bin_real ), 1 ) );
	if ( r12_bin > max_bin_ ) return Energy( 0.0 );
	Distance const r12_alpha( r12_bin_real - ( r12_bin - Distance( 0.5 ) ) );
	assert( ( r12_alpha >= Distance( 0.0 ) ) && ( r12_alpha <= Distance( 1.0 ) ) );

	// Set the low and high bin averages for interpolation

	pair_lhood_ratio_low =
	 ( pair_corr_(aa1n,aa2n,e1,e2,r12_bin) + pair_corr_(aa2n,aa1n,e2,e1,r12_bin) ) * Probability( 0.5 );

	pair_lhood_ratio_high =  r12_bin == max_bin_ ? Probability( 1.0 ) :
	 ( r12_bin_real == Distance( 0.5 ) ? pair_lhood_ratio_low :
	 ( pair_corr_(aa1n,aa2n,e1,e2,r12_bin+1) + pair_corr_(aa2n,aa1n,e2,e1,r12_bin+1) ) * Probability( 0.5 ) );

	//std::cout << "pairE potential: residues " << res1.seqpos() << " & " << res2.seqpos() << " with pairE = ";
	//std::cout << -std::log( interpolated( r12_alpha, pair_lhood_ratio_low, pair_lhood_ratio_high ) ) << std::endl;


	// Return the energy
	pair_lhood_ratio = interpolated( r12_alpha, pair_lhood_ratio_low, pair_lhood_ratio_high );
	return -std::log( pair_lhood_ratio );
}


Energy
PairEPotential::pair_term_energy_and_deriv(
	conformation::Residue const & res1,
	int res1_num_10A_neighbors,
	conformation::Residue const & res2,
	int res2_num_10A_neighbors,
	EnergyDerivative & dpairE_dr
) const
{
	Probability pair_lhood_ratio( 1.0 ), pair_lhood_ratio_low( 1.0 ), pair_lhood_ratio_high( 1.0 );

	Energy pairE = pair_term_energy( res1, res1_num_10A_neighbors, res2, res2_num_10A_neighbors,
		pair_lhood_ratio, pair_lhood_ratio_high, pair_lhood_ratio_low);

	dpairE_dr = -(1/pair_lhood_ratio) *
	 ( pair_lhood_ratio_high - pair_lhood_ratio_low ) / pair_score_bin_range_;

	 return pairE;

}

}
}
