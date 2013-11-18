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
/// @author Amelie Stein (amelie.stein@ucsf.edu) Oct 2012 -- rama2b lookup table for loop modeling

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
#include <numeric/random/random.hh>

// Utility Headers
#include <utility/io/izstream.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

// AS -- to get access to get_torsion_bin()
#include <core/conformation/util.hh>

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

// @brief Auto-generated virtual destructor
Ramachandran2B::~Ramachandran2B() {}

basic::Tracer T("core.scoring.Ramachandran2B");

Real const Ramachandran2B::binw_( 10.0 );
//AS
Real const Ramachandran2B::rama_sampling_thold_(0.00075 ); // only sample torsions with Rama prob above this value -- note that values are directly copied from the Ramachandran.cc implementation, might need tweaking
Real const Ramachandran2B::rama_sampling_factor_( 10.0 ); // factor for increased precision of Rama sampling table
ObjexxFCL::FArray4D< Real > Ramachandran2B::left_ram_probabil_( Ramachandran2B::n_phi_, Ramachandran2B::n_psi_, Ramachandran2B::n_aa_, Ramachandran2B::n_aa_ );
ObjexxFCL::FArray4D< Real > Ramachandran2B::right_ram_probabil_( Ramachandran2B::n_phi_, Ramachandran2B::n_psi_, Ramachandran2B::n_aa_, Ramachandran2B::n_aa_ );

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
			left_ram_probabil_( phiIndex, psiIndex, aa_num_left, aa_num ) = tProb; 
		} else if( aa_num_right != nullaa ) {
			ram_energ_right_( phiIndex, psiIndex, aa_num, aa_num_right ) = tEnergy;
			ram_entropy_right_( aa_num, aa_num_right ) += entropy;
			right_ram_probabil_( phiIndex, psiIndex, aa_num, aa_num_right ) = tProb; 

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

	///////////////////////////////////////////////////////////////////////////////
    /// Initialize the table holding the sample-able torsion space for each residue and its left neighbor
    /// with each torsion given indices proportionate to its probability
    /// @author Amelie Stein amelie.stein@ucsf.edu
    /// @brief based on the corresponding function for the Ramachandran class, but with adapted dimensions to accommodate the two neighbors
    void
    Ramachandran2B::init_rama_sampling_table_left( const char torsion_bin )  // to be adapted !! 
    {
        //rama_sampling_table_.resize(n_aa_);
        utility::vector1< utility::vector1< utility::vector1< utility::vector1< Real > > > >  current_rama_sampling_table;
        current_rama_sampling_table.resize(n_aa_);
        //int ss_type=3; // WARNING -- data for neighbor-dependent Rama is only available for 3 in the current file (Rama08.dat) -- and thus currently this information isn't even encoded in the rama_probabil_ table
        FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
        for (int left_aa=1; left_aa<=n_aa_; left_aa++) { // loop over all residue types
            current_rama_sampling_table[left_aa].resize(n_aa_);
            for (int aa=1; aa<=n_aa_; aa++) { // loop over all residue types
                current_rama_sampling_table[left_aa][aa].resize(n_aa_);
                FArray2A< Real > const rama_for_res( left_ram_probabil_(1, 1, left_aa, aa), zero_index, zero_index ); // is this correct? Why are we not accessing each one for the specific phi/psi combination?
                Size max_allowed = n_phi_ * n_psi_;
                Size actual_allowed = 0;
                Real min_val = 1.0; // minimum probability (above rama_sampling_thold_) observed for this residue
                Real max_val = 0.0; // maximum probability (above rama_sampling_thold_) observed for this residue
                utility::vector1< utility::vector1< Real> > res_torsions( max_allowed ); // double vector of allowed torsions for this residue
                utility::vector1< Real > res_probs( max_allowed ); // rama probs of allowed torsions for this residue (coupled to res_torsions by index)
				// current_rama_sampling_table[left_aa][aa][right_aa].resize(max_allowed); // I think this is resized later anyway
                for (int i=0; i<n_phi_; i++) {
                    for (int j=0; j<n_psi_; j++) {
                        Real res_prob = rama_for_res(i,j);
						
						//std::cerr << res_prob << std::endl;
						
                        if ( res_prob > rama_sampling_thold_ ) {
                            actual_allowed++; 
                            if (res_prob < min_val) min_val = res_prob;
                            else if (res_prob > max_val) max_val = res_prob;
							//std::cout << res_prob << std::endl;
							//res_probs[actual_allowed] = res_prob;
							
                            utility::vector1< Real > torsion(2);
                            Real cur_phi, cur_psi;
                            if (i <= n_phi_ / 2) {
                                cur_phi = i;
                            }
                            else {
                                cur_phi = 0 - (n_phi_ - i);
                            }
                            if (j <= n_psi_ / 2) {
                                cur_psi = j;
                            }
                            else {
                                cur_psi = 0 - (n_psi_ - j);
                            }
                            
							
                            char cur_tb = ' ';
                            if (torsion_bin != 'X') 
                                cur_tb = core::conformation::get_torsion_bin(cur_phi * 10, cur_psi * 10); //  AS -- how can we get the factor properly / without hard-coding? - also: this takes very long... 
                            if (torsion_bin == 'X' || cur_tb == torsion_bin) { 
                                actual_allowed++;
                                res_probs[actual_allowed] = res_prob;
                                torsion[1] = static_cast <Real> (/*i*/ cur_phi * binw_);  // phi
                                torsion[2] = static_cast <Real> (/*j*/ cur_psi * binw_);  // psi
                                res_torsions[actual_allowed] = torsion;
                            }
                            //rama_sampling_table_[aa][++actual_allowed] = torsion;
                            //} else {
                            //std::cerr << "warning -- discarding phi/psi " << cur_phi << "/" << cur_psi << " because they're not in torsion bin " << torsion_bin << std::endl;
                        }
                    }
                }
				
                /* not sure how to adapt this... 
				 if( ((int)aa < (int)1) || ((int)aa > (int)current_rama_sampling_table.size()) ){
				 std::cerr << "AA exceeded size of rama_sampling_table_ AA=" + ObjexxFCL::string_of( aa ) + " size()=" + ObjexxFCL::string_of( current_rama_sampling_table.size() );
				 continue; // Avoid death.
				 }
                 */
				
				
				
				// now populate the rama_sampling_table_ so the torsions are given index space proporionate to their probs
                current_rama_sampling_table[left_aa][aa].resize(Size(actual_allowed * (max_val / min_val) * rama_sampling_factor_));
                Size index=0; // to increment the index into the aa's rama_sampling_table_ vector
				//std::cout << "for aa " << aa << ":" << std::endl;
                for (Size tor = 1; tor <= actual_allowed; tor++) {
                    Size n_indices = Size(( res_probs[tor] / min_val ) * rama_sampling_factor_);
					//std::cout << "n_indices for torsion " << tor << ": " << n_indices << std::endl;
                    for (Size ind=1; ind<=n_indices; ind++) {
                        index++;
                        if( (int(index) < 1) || (int(index) > (int)current_rama_sampling_table[left_aa][aa].size()) ){
                            std::cerr <<  "index exceeded size of rama_sampling_table_[aa] index=" +
                            ObjexxFCL::string_of( (index) ) +
                            " rama_sampling_table_[aa].size()=" +
                            ObjexxFCL::string_of( current_rama_sampling_table[left_aa][aa].size() ) +
                            " AA=" + ObjexxFCL::string_of( aa );
							
                            continue; // avoid certain death - we dont yet understand why its failing here occasionally
                        }
                        if( ((int)tor < 1) || ((int)tor > (int)res_torsions.size()) ){
                            std::cerr << "tor exceeded size of rama_sampling_table_[aa] index=" +
                            ObjexxFCL::string_of( tor ) +
                            " res_torsions.size()=" + ObjexxFCL::string_of( res_torsions.size() ) +
                            " AA=" + ObjexxFCL::string_of( aa );
                            continue; // avoid certain death - we dont yet understand why its failing here occasionally
                        }
						
                        current_rama_sampling_table[left_aa][aa][index] = res_torsions[tor];
                    }
                }
                current_rama_sampling_table[left_aa][aa].resize(index);
				//std::cerr << "populating (left) table for torsion bin " << torsion_bin << " and " << AA(left_aa) << " " << AA(aa) << " -- " << index << " entries" << std::endl; // AS debug
				
            } // loop over aa
        } // loop over left_aa
        
        if (torsion_bin == 'X') {
            left_rama_sampling_table_ = current_rama_sampling_table;
        }
        core::Size tb_index = get_torsion_bin_index(torsion_bin);
        //std::cerr << "storing table for torsion bin " << torsion_bin << " --> " << tb_index << std::endl;
        if (left_rama_sampling_table_by_torsion_bin_.size() < tb_index)	
            left_rama_sampling_table_by_torsion_bin_.resize(tb_index);
        //std::cerr << "table resized" << std::endl;
        left_rama_sampling_table_by_torsion_bin_[tb_index] = current_rama_sampling_table;
        
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    /// Sample phi/psi torsions with probabilities proportionate to their
    /// Ramachandran probabilities
    /// Note -- this function had previously required that the option
    /// loops::nonpivot_torsion_sampling be active.  This function now
    /// performs a just-in-time check to initialize these tables the first
    /// time they are requested -- To properly multi-thread this code, the
    /// function should nab a mutex so that no two threads try to execute
    /// the code at once.
    void
    Ramachandran2B::random_phipsi_from_rama_left(
												 AA const left_aa,
												 AA const pos_aa,
												 Real & phi,
												 Real & psi
												 ) const    {
        
        if ( left_rama_sampling_table_.size() == 0 ) {
            /// Danger -- not threadsafe.
            const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_left('X');
        }
        
        Size n_torsions = left_rama_sampling_table_[left_aa][pos_aa].size();
        Size index = numeric::random::random_range(1, n_torsions);
        
        // following lines set phi and set to values drawn proportionately from Rama space
        // AS Nov 2013 - note that phi/psi bins are lower left corners, so we only want to add "noise" from a uniform distribution
        phi = left_rama_sampling_table_[left_aa][pos_aa][index][1] + numeric::random::uniform() * binw_;
        psi = left_rama_sampling_table_[left_aa][pos_aa][index][2] + numeric::random::uniform() * binw_;
    }
    
	
    
    
    ///////////////////////////////////////////////////////////////////////////////
    /// Initialize the table holding the sample-able torsion space for each residue and its right neighbor
    /// with each torsion given indices proportionate to its probability
    /// @author Amelie Stein amelie.stein@ucsf.edu
    /// @brief based on the corresponding function for the Ramachandran class, but with adapted dimensions to accommodate the two neighbors
    void
    Ramachandran2B::init_rama_sampling_table_right( const char torsion_bin )  
    {
        //rama_sampling_table_.resize(n_aa_);
        utility::vector1< utility::vector1< utility::vector1< utility::vector1< Real > > > >  current_rama_sampling_table;
        current_rama_sampling_table.resize(n_aa_);
        //int ss_type=3; // WARNING -- data for neighbor-dependent Rama is only available for 3 in the current file (Rama08.dat) -- and thus currently this information isn't even encoded in the rama_probabil_ table
        FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
        for (int aa=1; aa<=n_aa_; aa++) { // loop over all residue types
            current_rama_sampling_table[aa].resize(n_aa_);
            for (int right_aa=1; right_aa<=n_aa_; right_aa++) { // loop over all residue types
                FArray2A< Real > const rama_for_res( right_ram_probabil_(1, 1, aa, right_aa), zero_index, zero_index ); // is this correct? Why are we not accessing each one for the specific phi/psi combination?
                Size max_allowed = n_phi_ * n_psi_;
                Size actual_allowed = 0;
                Real min_val = 1.0; // minimum probability (above rama_sampling_thold_) observed for this residue
                Real max_val = 0.0; // maximum probability (above rama_sampling_thold_) observed for this residue
                utility::vector1< utility::vector1< Real> > res_torsions( max_allowed ); // double vector of allowed torsions for this residue
                utility::vector1< Real > res_probs( max_allowed ); // rama probs of allowed torsions for this residue (coupled to res_torsions by index)
				// current_rama_sampling_table[left_aa][aa][right_aa].resize(max_allowed); // I think this is resized later anyway
                for (int i=0; i<n_phi_; i++) {
                    for (int j=0; j<n_psi_; j++) {
                        Real res_prob = rama_for_res(i,j);
						
						//std::cerr << res_prob << std::endl;
						
                        if ( res_prob > rama_sampling_thold_ ) {
                            actual_allowed++; 
                            if (res_prob < min_val) min_val = res_prob;
                            else if (res_prob > max_val) max_val = res_prob;
							//std::cout << res_prob << std::endl;
							//res_probs[actual_allowed] = res_prob;
							
                            utility::vector1< Real > torsion(2);
                            Real cur_phi, cur_psi;
                            if (i <= n_phi_ / 2) {
                                cur_phi = i;
                            }
                            else {
                                cur_phi = 0 - (n_phi_ - i);
                            }
                            if (j <= n_psi_ / 2) {
                                cur_psi = j;
                            }
                            else {
                                cur_psi = 0 - (n_psi_ - j);
                            }
							
                            char cur_tb = ' ';
                            if (torsion_bin != 'X') 
                                cur_tb = core::conformation::get_torsion_bin(cur_phi * 10, cur_psi * 10); //  AS -- how can we get the factor properly / without hard-coding? - also: this takes very long... 
                            if (torsion_bin == 'X' || cur_tb == torsion_bin) { 
                                actual_allowed++;
                                res_probs[actual_allowed] = res_prob;
                                torsion[1] = static_cast <Real> (/*i*/ cur_phi * binw_);  // phi
                                torsion[2] = static_cast <Real> (/*j*/ cur_psi * binw_);  // psi
                                res_torsions[actual_allowed] = torsion;
                            }
							//rama_sampling_table_[aa][++actual_allowed] = torsion;
							//} else {
							//std::cerr << "warning -- discarding phi/psi " << cur_phi << "/" << cur_psi << " because they're not in torsion bin " << torsion_bin << std::endl;
                        }
                    }
                }
				
				/* not sure how to adapt this... 
				 if( ((int)aa < (int)1) || ((int)aa > (int)current_rama_sampling_table.size()) ){
				 std::cerr << "AA exceeded size of rama_sampling_table_ AA=" + ObjexxFCL::string_of( aa ) + " size()=" + ObjexxFCL::string_of( current_rama_sampling_table.size() );
				 continue; // Avoid death.
				 }
				 */
				
				
				
				// now populate the rama_sampling_table_ so the torsions are given index space proporionate to their probs
                current_rama_sampling_table[aa][right_aa].resize(Size(actual_allowed * (max_val / min_val) * rama_sampling_factor_));
                Size index=0; // to increment the index into the aa's rama_sampling_table_ vector
				//std::cout << "for aa " << aa << ":" << std::endl;
                for (Size tor = 1; tor <= actual_allowed; tor++) {
                    Size n_indices = Size(( res_probs[tor] / min_val ) * rama_sampling_factor_);
					//std::cout << "n_indices for torsion " << tor << ": " << n_indices << std::endl;
                    for (Size ind=1; ind<=n_indices; ind++) {
                        index++;
                        if( (int(index) < 1) || (int(index) > (int)current_rama_sampling_table[aa][right_aa].size()) ){
                            std::cerr <<  "index exceeded size of rama_sampling_table_[aa] index=" +
                            ObjexxFCL::string_of( (index) ) +
                            " rama_sampling_table_[aa].size()=" +
                            ObjexxFCL::string_of( current_rama_sampling_table[aa][right_aa].size() ) +
                            " AA=" + ObjexxFCL::string_of( aa );
							
                            continue; // avoid certain death - we dont yet understand why its failing here occasionally
                        }
                        if( ((int)tor < 1) || ((int)tor > (int)res_torsions.size()) ){
                            std::cerr << "tor exceeded size of rama_sampling_table_[aa] index=" +
                            ObjexxFCL::string_of( tor ) +
                            " res_torsions.size()=" + ObjexxFCL::string_of( res_torsions.size() ) +
                            " AA=" + ObjexxFCL::string_of( aa );
                            continue; // avoid certain death - we dont yet understand why its failing here occasionally
                        }
						
                        current_rama_sampling_table[aa][right_aa][index] = res_torsions[tor];
                    }
                }
                current_rama_sampling_table[aa][right_aa].resize(index);
				
				//std::cerr << "populating (right) table for torsion bin " << torsion_bin << " and " << AA(aa) << " " << AA(right_aa) << " -- " << index << " entries" << std::endl;
				
				// problem: some of these (highly specialized) bins are empty, but it can still happen that we request the respective bin
				// -- workaround 1: use the data across all bins (X) instead
				// -- workaround 2: use the minimum fraction from both the left and the right side, to make sure we don't run into this
				
				
            } // loop over right_aa
        } // loop over aa
        
        if (torsion_bin == 'X') {
            right_rama_sampling_table_ = current_rama_sampling_table;
        }
        core::Size tb_index = get_torsion_bin_index(torsion_bin);
        //std::cerr << "storing table for torsion bin " << torsion_bin << " --> " << tb_index << std::endl;
        if (right_rama_sampling_table_by_torsion_bin_.size() < tb_index)	
            right_rama_sampling_table_by_torsion_bin_.resize(tb_index);
        //std::cerr << "table resized" << std::endl;
        right_rama_sampling_table_by_torsion_bin_[tb_index] = current_rama_sampling_table;
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    /// Sample phi/psi torsions with probabilities proportionate to their
    /// Ramachandran probabilities
    /// Note -- this function had previously required that the option
    /// loops::nonpivot_torsion_sampling be active.  This function now
    /// performs a just-in-time check to initialize these tables the first
    /// time they are requested -- To properly multi-thread this code, the
    /// function should nab a mutex so that no two threads try to execute
    /// the code at once.
    void
    Ramachandran2B::random_phipsi_from_rama_right(
                                                  AA const pos_aa,
                                                  AA const right_aa,
                                                  Real & phi,
                                                  Real & psi
                                                  ) const    {
        
        if ( right_rama_sampling_table_.size() == 0 ) {
            /// Danger -- not threadsafe.
            const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_right('X');
        }
        
        Size n_torsions = right_rama_sampling_table_[pos_aa][right_aa].size();
        Size index = numeric::random::random_range(1, n_torsions);
        
        // following lines set phi and set to values drawn proportionately from Rama space
        // AS Nov 2013 - note that phi/psi bins are lower left corners, so we only want to add "noise" from a uniform distribution
        phi = right_rama_sampling_table_[pos_aa][right_aa][index][1] + numeric::random::uniform() * binw_;
        psi = right_rama_sampling_table_[pos_aa][right_aa][index][2] + numeric::random::uniform() * binw_;
    }
    
    
    core::Size Ramachandran2B::get_torsion_bin_index(char torsion_bin) const 
    { 
        return toupper(torsion_bin) - toupper('A') + 1;
    }
	
    // just to avoid code duplication
    void
    Ramachandran2B::init_rama_sampling_tables_by_torsion_bin() 
    {
		const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_left( 'X' ); // to allow wildcards in the torsion string, and to have "backup" data to access in case the selected bin is empty for a given combination (e.g., V-Y (G) doesn't have any entries, but if we filled it based on the right side we may request a G bin anyway)
        const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_left( 'A' );
        const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_left( 'B' );
        const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_left( 'E' );
        const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_left( 'G' );
		
		const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_right( 'X' ); // to allow wildcards in the torsion string
		const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_right( 'A' );
        const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_right( 'B' );
        const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_right( 'E' );
        const_cast< Ramachandran2B * > (this)->init_rama_sampling_table_right( 'G' );
		
		//std::cerr << " rama2b initialization by torsion bin done" << std::endl;
        
    }
	
	
	
	///////////////////////////////////////////////////////////////////////////////
	/// Sample phi/psi torsions with probabilities proportionate to their
	/// Ramachandran probabilities -- this version performs lookup restricted to specified torsion bins
	/// based on random_phipsi_from_rama and has the same issue for parallel running	
	
	/// @author Amelie Stein (amelie.stein@ucsf.edu)
	/// @date Fri May 11 15:52:01 PDT 2012
	/// @details returns a random phi/psi combination within the given torsion bin -- WARNING: this will only work for the torsion bins that are currently implemented
	
	void
	Ramachandran2B::random_phipsi_from_rama_by_torsion_bin_left(
																AA const left_aa,
																AA const pos_aa,
																Real & phi,
																Real & psi, 
																char const torsion_bin 
																) const
	{
		if (left_rama_sampling_table_by_torsion_bin_.size() == 0) {
			const_cast< Ramachandran2B * > (this)->init_rama_sampling_tables_by_torsion_bin(); // covers both left and righ
			// not threadsafe either
		}
		
		core::Size tb_index = get_torsion_bin_index(torsion_bin);
		
		Size n_torsions = left_rama_sampling_table_by_torsion_bin_[tb_index][left_aa][pos_aa].size();
		
		
		//if (n_torsions == 0) { // debugging
		//	std::cerr << " error -- no entries found for (left) " << torsion_bin << " -- " << AA(left_aa) << " " << AA(pos_aa) << std::endl;
		//}
		
		Size index = numeric::random::random_range(1, n_torsions);
		
		// following lines set phi and set to values drawn proportionately from Rama space
        // AS Nov 2013 - note that phi/psi bins are lower left corners, so we only want to add "noise" from a uniform distribution
		phi = left_rama_sampling_table_by_torsion_bin_[tb_index][left_aa][pos_aa][index][1] + numeric::random::uniform() * binw_;
		psi = left_rama_sampling_table_by_torsion_bin_[tb_index][left_aa][pos_aa][index][2] + numeric::random::uniform() * binw_;
		
	} // random_phipsi_from_rama_by_torsion_bin_left
	
	void
	Ramachandran2B::random_phipsi_from_rama_by_torsion_bin_right(
																 AA const pos_aa,
																 AA const right_aa,
																 Real & phi,
																 Real & psi, 
																 char const torsion_bin 
																 ) const
	{
		if (right_rama_sampling_table_by_torsion_bin_.size() == 0) {
			const_cast< Ramachandran2B * > (this)->init_rama_sampling_tables_by_torsion_bin(); // covers both left and righ
			// not threadsafe either
		}
		
		core::Size tb_index = get_torsion_bin_index(torsion_bin);
		
		Size n_torsions = right_rama_sampling_table_by_torsion_bin_[tb_index][pos_aa][right_aa].size();
		//if (n_torsions == 0) { // debugging
		//	std::cerr << " error -- no entries found for (right) " << torsion_bin << " -- "  << AA(pos_aa) << " " << AA(right_aa) << std::endl;
		//}
		
		Size index = numeric::random::random_range(1, n_torsions);
		
		// following lines set phi and set to values drawn proportionately from Rama space
        // AS Nov 2013 - note that phi/psi bins are lower left corners, so we only want to add "noise" from a uniform distribution
		phi = right_rama_sampling_table_by_torsion_bin_[tb_index][pos_aa][right_aa][index][1] + numeric::random::uniform() * binw_;
		psi = right_rama_sampling_table_by_torsion_bin_[tb_index][pos_aa][right_aa][index][2] + numeric::random::uniform() * binw_;
		
	} // random_phipsi_from_rama_by_torsion_bin_right
	
	
	
    void
    Ramachandran2B::get_entries_per_torsion_bin_left( 
													 AA const left_aa,
													 AA const pos_aa,
													 std::map< char, core::Size > & tb_frequencies ) const
    {
		// check if the tables are initialized, and if not, do so
        if (left_rama_sampling_table_by_torsion_bin_.size() == 0) 
            const_cast< Ramachandran2B * > (this)->init_rama_sampling_tables_by_torsion_bin();	
		tb_frequencies['A'] = left_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('A')][left_aa][pos_aa].size();
		tb_frequencies['B'] = left_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('B')][left_aa][pos_aa].size();
		tb_frequencies['E'] = left_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('E')][left_aa][pos_aa].size();
		tb_frequencies['G'] = left_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('G')][left_aa][pos_aa].size();
		tb_frequencies['X'] = left_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('X')][left_aa][pos_aa].size();
    }
	
	
    void
    Ramachandran2B::get_entries_per_torsion_bin_right(
													  AA const pos_aa,
													  AA const right_aa,
													  std::map< char, core::Size > & tb_frequencies ) const
    {
		// check if the tables are initialized, and if not, do so
        if (right_rama_sampling_table_by_torsion_bin_.size() == 0) 
            const_cast< Ramachandran2B * > (this)->init_rama_sampling_tables_by_torsion_bin();	
		tb_frequencies['A'] = right_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('A')][pos_aa][right_aa].size();
		tb_frequencies['B'] = right_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('B')][pos_aa][right_aa].size();
		tb_frequencies['E'] = right_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('E')][pos_aa][right_aa].size();
		tb_frequencies['G'] = right_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('G')][pos_aa][right_aa].size();
		tb_frequencies['X'] = right_rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('X')][pos_aa][right_aa].size();
    }
	
	


} // scoring
} // core
