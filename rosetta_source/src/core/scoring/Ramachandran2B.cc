// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran2B.cc
/// @brief  Neighbor Dependent Ramachandran potential class implementation
/// @author Guoli Wang

// Unit Headers
#include <core/scoring/Ramachandran2B.hh>

// Package Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ProteinTorsion.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/interpolation/periodic_range/half/interpolation.hh>

// Utility Headers
#include <utility/io/izstream.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2A.hh>



using namespace ObjexxFCL;

namespace core {
namespace scoring {

basic::Tracer T("core.scoring.Ramachandran2B");

Real const Ramachandran2B::binw_( 10.0 );

Ramachandran2B::Ramachandran2B() :
	ram_energ_( n_phi_, n_psi_, n_aa_, 0.0 ),
	ram_entropy_( n_aa_, 0.0 ),
	ram_energ_left_( n_phi_, n_psi_, n_aa_, n_aa_, 0.0 ),
	ram_entropy_left_( n_aa_, n_aa_, 0.0 ),
	ram_energ_right_( n_phi_, n_psi_, n_aa_, n_aa_, 0.0 ),
	ram_entropy_right_( n_aa_, n_aa_, 0.0 ),
	rama_score_limit_( 20 ) // customizable in the future, possibly by command line flags.
{
	read_rama();
}


///////////////////////////////////////////////////////////////////////////////

/// @brief evaluate rama score for each (protein) residue and store that score
/// in the pose.energies() object
void
Ramachandran2B::eval_rama_score_all(
	pose::Pose & pose,
	ScoreFunction const & scorefxn
) const
{
	if ( scorefxn.has_zero_weight( rama ) ) return; // unnecessary, righ?

	//double rama_sum = 0.0;

	// in pose mode, we use fold_tree.cutpoint info to exclude terminus
	// residues from rama calculation. A cutpoint could be either an artificial
	// cutpoint such as loop cutpoint or a real physical chain break such as
	// multiple-chain complex. For the artificial cutpoint, we may need to
	// calculate rama scores for cutpoint residues, but for the real chain break
	// cutpoint, we don't want to do that. So here we first loop over all the
	// residue in the protein and exclude those ones which are the cutpoints.
	// Then we loop over the cutpoint residues and add rama score for residues
	// at artificial cutpoints, i.e., cut_weight != 0.0, which means that
	// jmp_chainbreak_score is also calculated for this cutpoint. Note that the
	// default value for cut_weight here is dependent on whether
	// jmp_chainbreak_weight is set. This is to ensure that rama score for
	// termini residues are not calculated when jmp_chainbreak_weight is 0.0,
	// e.g normal pose docking.

	int const total_residue = pose.total_residue();

	// retrieve cutpoint info // apl do we actually need this data?
	// if so, Pose must provide it 'cause we're offing all global data
	//
	//kinematics::FoldTree const & fold_tree(
	//		pose.fold_tree() );
	//int const n_cut( fold_tree.num_cutpoint() );

	//FArray1D< Real > cut_weight( n_cut,
	//	scorefxns::jmp_chainbreak_weight == 0.0 ? 0.0 : 1.0 ); // apl need to handle

	//if( cut_weight.size1() == scorefxns::cut_weight.size1() )
	//	cut_weight = scorefxns::cut_weight;

	// exclude chain breaks

	Energies & pose_energies( pose.energies() );

	for ( int ii = 1; ii <= total_residue; ++ii )
	{
		if ( pose.residue(ii).is_protein()  && ! pose.residue(ii).is_terminus()  )
		{
			Real rama_score,dphi,dpsi;
			eval_rama_score_residue(pose.residue(ii), pose.residue(ii-1).aa(), pose.residue(ii+1).aa(), rama_score, dphi, dpsi);
			T << "Rama:eval_all: residue " << ii << " " << pose.residue(ii).name() <<
				" " << ii-1 << " " << pose.residue(ii-1).name() << " " << ii+1 << " " <<
				pose.residue(ii+1).name() << " = " << rama_score << std::endl;
			pose_energies.onebody_energies( ii )[rama] = rama_score;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran2B::write_rama_score_all( Pose const & /*pose*/ ) const
{}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran2B::eval_rama_score_residue(
	conformation::Residue const & rsd,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	//assert( pose.residue(res).is_protein() );
	assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 || rsd.is_terminus() ) { // begin or end of chain
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	eval_rama_score_residue( rsd.aa(), phi, psi, rama, drama_dphi, drama_dpsi );
}

///////////////////////////////////////////////////////////////////////////////
// modified by GL
void
Ramachandran2B::eval_rama_score_residue(
	conformation::Residue const &center,
	chemical::AA const left_aa,
	chemical::AA const right_aa,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	assert( center.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( center.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( center.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	if( ! basic::options::option[ basic::options::OptionKeys::score::ramaneighbors ] ) {
		eval_rama_score_residue( center.aa(), phi, psi, rama, drama_dphi, drama_dpsi );
	} else {
		//if(left.seqpos() == center.seqpos()) {
		//	Ramachandran::RamaE_Upper(center, right_aa, drama_dphi, drama_dpsi);
		//} else if(right.seqpos() == center.seqpos()) {
		//	Ramachandran::RamaE_Lower(center, left_aa, drama_dphi, drama_dpsi);
		//} else {
		Real rama_L(0.0), drama_dphi_L(0.0), drama_dpsi_L(0.0);
		Real rama_R(0.0), drama_dphi_R(0.0), drama_dpsi_R(0.0);
		Real rama_0(0.0), drama_dphi_0(0.0), drama_dpsi_0(0.0);
		rama_L = RamaE_Lower(center, left_aa, drama_dphi_L, drama_dpsi_L);
		rama_R = RamaE_Upper(center, right_aa, drama_dphi_R, drama_dpsi_R);
		rama_0 = RamaE(center, drama_dphi_0, drama_dpsi_0);

		rama = rama_L + rama_R - rama_0;
		drama_dphi = drama_dphi_L + drama_dphi_R - drama_dphi_0;
		drama_dpsi = drama_dpsi_L + drama_dpsi_R - drama_dpsi_0;
		//}
	}
}

void
Ramachandran2B::IdealizeRamaEnergy(
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi,
	Real const entropy,
	FArray2A< Real > const & rama_for_res
) const
{
	using namespace numeric::interpolation::periodic_range::half;
	Real interp_E = bilinearly_interpolated( phi, psi, binw_, n_phi_, rama_for_res, drama_dphi, drama_dpsi );
	// rama = IdealizeRamaEnergy(ram_entropy_(center.aa(), leftIndex, rightIndex), interp_E, drama_dphi, drama_dpsi);
	rama = entropy + interp_E;
	//	std::cout << "Rama::eval_res: " <<  interp_E << " rama " << rama << std::endl;

	if ( ! basic::options::option[basic::options::OptionKeys::corrections::score::rama_not_squared] ) {
		if ( rama > 1.0 ) {
			Real rama_squared = rama * rama;
			if ( rama_squared > rama_score_limit_ ) {
				drama_dphi = 0.0;
				drama_dpsi = 0.0;
				rama = rama_score_limit_;
			} else {
				drama_dphi *= 2.0 * rama;
				drama_dpsi *= 2.0 * rama;
				rama = rama_squared;
			}
		}
	}

	// std::cout << " rama: " << rama << " dphi " << drama_dphi << " dpsi " << drama_dpsi << std::endl;
}

// end modification

///////////////////////////////////////////////////////////////////////////////
// modified by GL according to Andrew's suggestion
Real
Ramachandran2B::RamaE_Lower(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor
) const
{
	Real drama_dphi, drama_dpsi;
	return RamaE_Lower( rsd, neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Lower(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		return 0.0;
	}

	// if neighbor independent protocol is selected, return 0.0
	if( ! option[ score::ramaneighbors ] ) {
		return 0.0;
	}

	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	FArray2A< Real > const rama_for_res( ram_energ_left_(1, 1, rsd.aa(), neighbor), zero_index, zero_index );
	Real entropy = ram_entropy_left_(rsd.aa(), neighbor);

	Real rama;
	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
	return rama;
}
Real
Ramachandran2B::RamaE_Upper(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor
) const
{
	Real drama_dphi, drama_dpsi;
	return RamaE_Upper( rsd, neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Upper(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		return 0.0;
	}

	// if neighbor independent protocol is selected, return 0.0
	if( ! basic::options::option[ basic::options::OptionKeys::score::ramaneighbors ] )
	{
		return 0.0;
	}

	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	FArray2A< Real > const rama_for_res( ram_energ_right_(1, 1, rsd.aa(), neighbor), zero_index, zero_index );
	Real entropy = ram_entropy_right_(rsd.aa(), neighbor);

	Real rama;
	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
	return rama;
}

Real
Ramachandran2B::RamaE(
	conformation::Residue const &rsd
) const
{
	Real drama_dphi(0.0), drama_dpsi(0.0);
	return RamaE( rsd, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE(
	conformation::Residue const &rsd,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{

	using namespace numeric;

	assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		return 0.0;
	}

	Real ramaE(0.0);
	eval_rama_score_residue( rsd.aa(), phi, psi, ramaE, drama_dphi, drama_dpsi );
	return ramaE;
}

///////////////////////////////////////////////////////////////////////////////
///
Real
Ramachandran2B::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi
) const
{

	Real rama, drama_dphi, drama_dpsi;
	eval_rama_score_residue( res_aa, phi, psi, rama, drama_dphi, drama_dpsi );
	return rama;
}

///////////////////////////////////////////////////////////////////////////////
///
void
Ramachandran2B::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	FArray2A< Real > const rama_for_res( ram_energ_(1, 1, res_aa ), zero_index, zero_index );
	Real entropy = ram_entropy_(res_aa);

	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
}


// Guoli Wang
void
Ramachandran2B::read_rama()
{
	using namespace basic::options;

	using namespace basic::options::OptionKeys;

	int aa_num( 0 ), aa_num_left( 0 ), aa_num_right( 0 );
	int phi_bin( 0 ), psi_bin( 0 ), ss_type( 0 );
	int tCounts( 0 );
	Real tProb( 0.0 ), tEnergy( 0.0 );

	Size line_count( 0 );


	//utility::io::izstream  iunit;
#ifndef WIN32
#ifndef __CYGWIN__
	clock_t starttime = clock();
#endif
#endif
	std::string energyFileName = basic::options::option[ in::file::rama2b_map ]().name() ; // "wrapHDPprobs36.both";
	T << "Read in ramachandran map: " <<  energyFileName << std::endl;
	utility::io::izstream iRamaEnergy;
	basic::database::open( iRamaEnergy, energyFileName );
	while( ! iRamaEnergy.eof() ) {
		++line_count;
		iRamaEnergy >> aa_num >> aa_num_left >> aa_num_right >> ss_type >> phi_bin >> psi_bin >> tCounts >> tProb >> tEnergy;
		// std::cout << " aa_num " << aa_num << " aa_num_left " << aa_num_left << " aa_num_right " << aa_num_right << " ss_type " << ss_type <<
		// 			" phi_bin " << phi_bin << " psi_bin " << psi_bin << " tProb " << tProb << " tEnergy " << tEnergy << std::endl;
		if(aa_num > n_aa_) continue;

		int phiIndex = phi_bin / 10 + 1;
		int psiIndex = psi_bin / 10 + 1;
		Real entropy = -1.0 * tProb * tEnergy;

		if( aa_num_left == nullaa && aa_num_right == nullaa ) {
			ram_energ_( phiIndex, psiIndex, aa_num ) = tEnergy;
			ram_entropy_( aa_num ) += entropy;
		} else if( aa_num_left != nullaa ) {
			ram_energ_left_( phiIndex, psiIndex, aa_num, aa_num_left ) = tEnergy;
			ram_entropy_left_( aa_num, aa_num_left ) += entropy;
		} else if( aa_num_right != nullaa ) {
			ram_energ_right_( phiIndex, psiIndex, aa_num, aa_num_right ) = tEnergy;
			ram_entropy_right_( aa_num, aa_num_right ) += entropy;
		}
	}

	iRamaEnergy.close();
#ifndef WIN32
#ifndef __CYGWIN__
	clock_t stoptime = clock();
	T << "Reading Rama from database took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << " seconds" << std::endl;
#endif
#endif
}


} // scoring
} // core
