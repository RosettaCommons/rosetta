// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran.cc
/// @brief  Ramachandran potential class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Modified by Vikram K. Mulligan (vmullig@uw.edu) for D-amino acids, noncanonical alpha-amino acids, etc.

// Unit Headers
#include <core/scoring/Ramachandran.hh>

// Package Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/interpolation/periodic_range/half/interpolation.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/string.functions.hh>

// AS -- to get access to get_torsion_bin()
#include <core/conformation/util.hh>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <utility/vector1.hh>
#include <sstream>



using namespace ObjexxFCL;

namespace core {
namespace scoring {

// @brief Auto-generated virtual destructor
Ramachandran::~Ramachandran() {}

typedef Ramachandran R;

bool R::rama_initialized_( false );
Real const R::binw_( 10.0 );
Real const R::rama_sampling_thold_( 0.00075 ); // only sample torsions with Rama prob above this value
Real const R::rama_sampling_factor_( 10.0 ); // factor for increased precision of Rama sampling table
ObjexxFCL::FArray4D< Real > R::ram_probabil_( R::n_phi_, R::n_psi_, 3, R::n_aa_ );
ObjexxFCL::FArray4D_int R::ram_counts_( R::n_phi_, R::n_psi_, 3, R::n_aa_ );
ObjexxFCL::FArray4D< Real > R::ram_energ_( R::n_phi_, R::n_psi_, 3, R::n_aa_);
ObjexxFCL::FArray2D< Real > R::ram_entropy_( 3, R::n_aa_ );

Ramachandran::Ramachandran()
{
	using namespace basic::options;
	read_rama(
		option[ OptionKeys::corrections::score::rama_map ]().name(),
		option[ OptionKeys::corrections::score::use_bicubic_interpolation ]);
}


Ramachandran::Ramachandran(
	std::string const & rama_map_filename,
	bool use_bicubic_interpolation
) {
	read_rama(rama_map_filename, use_bicubic_interpolation);
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Returns true if passed a core::chemical::AA corresponding to a
/// D-amino acid, and false otherwise.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
Ramachandran::is_canonical_d_aminoacid(
	AA const res_aa
) const {
	return core::chemical::is_canonical_D_aa(res_aa);
}

///////////////////////////////////////////////////////////////////////////////

/// @brief When passed a d-amino acid, returns the l-equivalent.  Returns
/// aa_unk otherwise.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::chemical::AA
Ramachandran::get_l_equivalent(
	AA const d_aa
) const {
	return core::chemical::get_L_equivalent(d_aa);
}


///////////////////////////////////////////////////////////////////////////////

/// @brief evaluate rama score for each (protein) residue and store that score
/// in the pose.energies() object
void
Ramachandran::eval_rama_score_all(
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
		if ( pose.residue(ii).is_protein()  && ! pose.residue(ii).is_terminus() && ! pose.residue(ii).is_virtual_residue() )
		{
			Real rama_score,dphi,dpsi;
			if(is_normally_connected(pose.residue(ii))) {
				eval_rama_score_residue(pose.residue(ii),rama_score,dphi,dpsi);
				//printf("Residue %i is normal.\n", ii); fflush(stdout); //DELETE ME -- FOR TESTING ONLY
				//std::cout << "Rama: residue " << ii << " = " << rama_score << std::endl;
				pose_energies.onebody_energies( ii )[rama] = rama_score;
			} else {
				//printf("Residue %i: THIS SHOULD HAPPEN ONLY IF THIS RESIDUE HAS WEIRD CONNECTIONS.", ii); fflush(stdout); //DELETE ME -- FOR TESTING ONLY
				eval_rama_score_residue_nonstandard_connection(pose, pose.residue(ii),rama_score,dphi,dpsi);
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran::write_rama_score_all( Pose const & /*pose*/ ) const
{}

	///////////////////////////////////////////////////////////////////////////////
	/// Initialize the table holding the sample-able torsion space for each residue
	/// with each torsion given indices proportionate to its probability
	void
	Ramachandran::init_rama_sampling_table(
										   char const torsion_bin) const // torsion_bin defaults to 'X' unless specified
	{
		//rama_sampling_table_.resize(n_aa_);
		utility::vector1< utility::vector1< utility::vector1< Real > > > current_rama_sampling_table;
		current_rama_sampling_table.resize(n_aa_);
		int ss_type=3;
		FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
		for (int aa=1; aa<=n_aa_; aa++) { // loop over all residue types
			FArray2A< Real > const rama_for_res( ram_probabil_(1, 1, ss_type, aa), zero_index, zero_index );
			Size max_allowed = n_phi_ * n_psi_;
			Size actual_allowed = 0;
			Real min_val = 1.0; // minimum probability (above rama_sampling_thold_) observed for this residue
			Real max_val = 0.0; // maximum probability (above rama_sampling_thold_) observed for this residue
			utility::vector1< utility::vector1< Real> > res_torsions( max_allowed ); // double vector of allowed torsions for this residue
			utility::vector1< Real > res_probs( max_allowed ); // rama probs of allowed torsions for this residue (coupled to res_torsions by index)
			//rama_sampling_table_[aa].resize(max_allowed);
			for (int i=0; i<n_phi_; i++) {
				for (int j=0; j<n_psi_; j++) {
					Real res_prob = rama_for_res(i,j);
					if ( res_prob > rama_sampling_thold_ ) {
						//actual_allowed++; // AS: moved down into if statement, otherwise some entries are empty -> segfault
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
					}
				}
			}

			if( ((int)aa < (int)1) || ((int)aa > (int)current_rama_sampling_table.size()) ){
				std::cerr << "AA exceeded size of rama_sampling_table_ AA=" + ObjexxFCL::string_of( aa ) + " size()=" + ObjexxFCL::string_of( current_rama_sampling_table.size() );
				continue; // Avoid death.
			}

			// now populate the rama_sampling_table_ so the torsions are given index
			// space proporionate to their probs
			current_rama_sampling_table[aa].resize(Size(actual_allowed * (max_val / min_val) * rama_sampling_factor_));
			Size index=0; // to increment the index into the aa's rama_sampling_table_ vector
			for (Size tor = 1; tor <= actual_allowed; tor++) {
				Size n_indices = Size(( res_probs[tor] / min_val ) * rama_sampling_factor_);
				for (Size ind=1; ind<=n_indices; ind++) {
					index++;
					if( (int(index) < 1) || (int(index) > (int)current_rama_sampling_table[aa].size()) ){
						std::cerr <<  "index exceeded size of rama_sampling_table_[aa] index=" +
						ObjexxFCL::string_of( (index) ) +
						" rama_sampling_table_[aa].size()=" +
						ObjexxFCL::string_of( current_rama_sampling_table[aa].size() ) +
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

					current_rama_sampling_table[aa][index] = res_torsions[tor];
				}
			}
			current_rama_sampling_table[aa].resize(index);
		}

		if (torsion_bin == 'X') {
			rama_sampling_table_ = current_rama_sampling_table;
		}
		core::Size tb_index = get_torsion_bin_index(torsion_bin);
		if (rama_sampling_table_by_torsion_bin_.size() < tb_index)
			rama_sampling_table_by_torsion_bin_.resize(tb_index);
		rama_sampling_table_by_torsion_bin_[tb_index] = current_rama_sampling_table;
	}

	void
	Ramachandran::init_uniform_sampling_table() const
	{
		uniform_sampling_table_.resize(n_aa_);

		for (int aa = 1; aa <= n_aa_; aa++) {
			FArray2A< Real >::IR const zero_index(0, n_phi_ - 1);
			FArray2A< Real > const rama_for_res(
					ram_probabil_(1, 1, 3, aa), zero_index, zero_index );

			for (int i = 0; i < n_phi_; i++) {
				for (int j = 0; j < n_psi_; j++) {

					// We want the sampling table to include each allowed phi/psi bin
					// exactly once, to ensure uniform sampling.  Forbidden phi/psi bins
					// should not appear in the table:

					if (rama_for_res(i, j) < rama_sampling_thold_) continue;

					utility::vector1<Real> torsion(2);
					torsion[1] = i * binw_;
					torsion[2] = j * binw_;

					uniform_sampling_table_[aa].push_back(torsion);
				}
			}
		}
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
	Ramachandran::random_phipsi_from_rama(
										  AA const res_aa,
										  Real & phi,
										  Real & psi
										  ) const
	{
		using numeric::random::uniform;

		if ( rama_sampling_table_.size() == 0 ) {
			/// Danger -- not threadsafe.
			init_rama_sampling_table('X');
		}

		// Use the equivalent L-amino acid if this is a D-amino acid.
		AA res_aa2 = res_aa;
		if (is_canonical_d_aminoacid(res_aa)) res_aa2=get_l_equivalent(res_aa);

		Size n_torsions = rama_sampling_table_[res_aa2].size();
		Size index = numeric::random::random_range(1, n_torsions);

		// Select a random phi/psi bin from the rama_sampling_table_.  This table
		// has been setup such that the number of times each psi/psi bin appears is
		// proportional to the rama weight of that bin.  Once a bin is chosen, a
		// specific phi/psi pair is chosen using a uniform distribution and no
		// interpolation.

		// AS Nov 2013 - note that phi/psi bins are lower left corners, so we only
		// want to add "noise" from a uniform distribution

		phi = rama_sampling_table_[res_aa2][index][1] + numeric::random::uniform() * binw_;
		psi = rama_sampling_table_[res_aa2][index][2] + numeric::random::uniform() * binw_;

		// Invert phi and psi if this is a D-amino acid.
		if (is_canonical_d_aminoacid(res_aa)) {
			phi=360.0-phi;
			psi=360.0-psi;
		}
	}

	void
	Ramachandran::uniform_phipsi_from_allowed_rama(
											AA const res_aa,
											Real & phi,
											Real & psi
											) const
	{
		using numeric::random::uniform;
		using numeric::random::random_range;

		if (uniform_sampling_table_.empty()) {
			init_uniform_sampling_table();
		}

		Size n_torsions = uniform_sampling_table_[res_aa].size();
		Size index = random_range(1, n_torsions);

		phi = uniform_sampling_table_[res_aa][index][1] + binw_ * uniform();
		psi = uniform_sampling_table_[res_aa][index][2] + binw_ * uniform();
	}

	bool
	Ramachandran::phipsi_in_allowed_rama(
											AA const aa,
											Real phi,
											Real psi
											) const
	{
		// Figure out which torsion bin the given phi/psi pair is in.  Angles must be
		// given in degrees, but need not be given in any specific range.

		phi = numeric::nonnegative_principal_angle_degrees(phi);
		psi = numeric::nonnegative_principal_angle_degrees(psi);

		Size phi_bin = floor(phi / binw_);
		Size psi_bin = floor(psi / binw_);

		// Return true if this bin has a non-zero rama probability.

		FArray2A<Real>::IR const zero_index(0, n_phi_ - 1);
		FArray2A<Real> const rama_for_res(
				ram_probabil_(1, 1, 3, aa), zero_index, zero_index);

		return rama_for_res(phi_bin, psi_bin) > rama_sampling_thold_;
	}

	bool
	Ramachandran::phipsi_in_forbidden_rama(
											AA const aa,
											Real phi,
											Real psi
											) const
	{
		return ! phipsi_in_allowed_rama(aa, phi, psi);
	}

	///////////////////////////////////////////////////////////////////////////////
	/// Sample phi/psi torsions with probabilities proportionate to their
	/// Ramachandran probabilities -- this version performs lookup restricted to specified torsion bins
	/// based on random_phipsi_from_rama and has the same issue for parallel running

	/// @author Amelie Stein (amelie.stein@ucsf.edu)
	/// @date Fri May 11 15:52:01 PDT 2012
	/// @details returns a random phi/psi combination within the given torsion bin -- WARNING: this will only work for the torsion bins that are currently implemented

	void
	Ramachandran::random_phipsi_from_rama_by_torsion_bin(
														 AA const res_aa,
														 Real & phi,
														 Real & psi,
														 char const torsion_bin
														 ) const
	{


		//utility::vector1< utility::vector1< utility::vector1< Real > > > current_rama_sampling_table; // depends on the torsion bin

		//if (rama_sampling_table_by_torsion_bin_.find(torsion_bin) == rama_sampling_table_by_torsion_bin_.end()) {
		if (rama_sampling_table_by_torsion_bin_.size() == 0) {
			init_rama_sampling_tables_by_torsion_bin();
			// not threadsafe either
			// initialize the table for this torsion bin and residue
			// TODO: store the original sampling table before doing this, and then write it back after storing the result of this function, so that we don't lose or overwrite this information
			//std::cerr << " generating rama sampling table for torsion bin " << torsion_bin << std::endl; // make sure that all of this is only done once

		}

		core::Size tb_index = get_torsion_bin_index(torsion_bin);
		/*
		 std::map< char, utility::vector1< utility::vector1< utility::vector1< Real > > > >::const_iterator m_iter = rama_sampling_table_by_torsion_bin_.find(torsion_bin);
		 if (m_iter != rama_sampling_table_by_torsion_bin_.end()) {
		 //std::cerr << "trying to fetch the table for torsion bin " << torsion_bin << std::endl;
		 current_rama_sampling_table = m_iter->second;
		 }

		 std::cerr << " rama sampling table size / # torsion bins" << rama_sampling_table_by_torsion_bin_.size() << std::endl;
		 std::cerr << " rama sampling table size " << rama_sampling_table_by_torsion_bin_[tb_index].size() << std::endl;
		 std::cerr << " -- res: " << res_aa << " " << rama_sampling_table_by_torsion_bin_[tb_index][res_aa].size() << std::endl;
		 */

		Size n_torsions = rama_sampling_table_by_torsion_bin_[tb_index][res_aa].size();
		Size index = numeric::random::random_range(1, n_torsions);

		/*
		 // check if for some reason the bins are not populated
		 std::cerr << current_rama_sampling_table[res_aa][index].size() << " index " << index << std::endl;
		 for (int i = 0; i < current_rama_sampling_table[res_aa][index].size(); i++)
		 std::cerr << current_rama_sampling_table[res_aa][index][i] << " - " << i << std::endl;

		 std::cerr << current_rama_sampling_table[res_aa][index][1] << " index " << index << " / 1" << std::endl;
		 std::cerr << current_rama_sampling_table[res_aa][index][2] << " index " << index << " / 2" << std::endl;
		 */
		// following lines set phi and set to values drawn proportionately from Rama space
        // AS Nov 2013 - note that phi/psi bins are lower left corners, so we only want to add "noise" from a uniform distribution
		phi = rama_sampling_table_by_torsion_bin_[tb_index][res_aa][index][1] + numeric::random::uniform() * binw_;
		psi = rama_sampling_table_by_torsion_bin_[tb_index][res_aa][index][2] + numeric::random::uniform() * binw_;
		// DJM: debug
		//std::cout << "res_aa: " << res_aa << std::endl;
		//std::cout << "phi: " << phi << std::endl;
		//std::cout << "psi: " << psi << std::endl;


	} // random_phipsi_from_rama_by_torsion_bin



	core::Size Ramachandran::get_torsion_bin_index(char torsion_bin) const
	{
		return toupper(torsion_bin) - toupper('A') + 1;
	}


	// just to avoid code duplication
	void
	Ramachandran::init_rama_sampling_tables_by_torsion_bin() const
	{
		init_rama_sampling_table( 'A' );
		init_rama_sampling_table( 'B' );
		init_rama_sampling_table( 'E' );
		init_rama_sampling_table( 'G' );
		init_rama_sampling_table( 'X' ); // to allow wildcards in the torsion string
	}

	void
	Ramachandran::get_entries_per_torsion_bin( AA const res_aa, std::map< char, core::Size > & tb_frequencies ) const
	{
		// check if the tables are initialized, and if not, do so
		if (rama_sampling_table_by_torsion_bin_.size() == 0) {
			init_rama_sampling_tables_by_torsion_bin();
		}
		tb_frequencies['A'] = rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('A')][res_aa].size();
		tb_frequencies['B'] = rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('B')][res_aa].size();
		tb_frequencies['E'] = rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('E')][res_aa].size();
		tb_frequencies['G'] = rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('G')][res_aa].size();
		tb_frequencies['X'] = rama_sampling_table_by_torsion_bin_[get_torsion_bin_index('X')][res_aa].size();
	}

///////////////////////////////////////////////////////////////////////////////
void
Ramachandran::eval_rama_score_residue_nonstandard_connection(
	core::pose::Pose const & mypose,
	conformation::Residue const & res,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {

	if(!basic::options::option[basic::options::OptionKeys::score::rama_score_nonstandard_connections] || res.connection_incomplete(1) || res.connection_incomplete(2)) { //If the rama_score_nonstandard_connections flag is not set, OR this is an open terminus, don't score this residue.
		rama=0.0;
		drama_dphi=0.0;
		drama_dpsi=0.0;
		return;
	}

	const core::conformation::Residue &lowerres=mypose.residue(res.residue_connection_partner(1));
	const core::conformation::Residue &upperres=mypose.residue(res.residue_connection_partner(2));

	//NOTE: The following assumes that we're scoring something with an alpha-amino acid backbone!  At the time of this writing (7 Feb 2013), rama checks that the residue being scored is either a standard L- or D-amino acid.
	Real phi = numeric::nonnegative_principal_angle_degrees(
				numeric::dihedral_degrees(
					lowerres.xyz( lowerres.residue_connect_atom_index(res.residue_connection_conn_id(1)) ),
					res.xyz("N"), //Position of N
					res.xyz("CA"), //Position of CA
					res.xyz("C") //Position of C
				)
			);
	Real psi=numeric::nonnegative_principal_angle_degrees(
				numeric::dihedral_degrees(
					res.xyz("N"), //Position of N
					res.xyz("CA"), //Position of CA
					res.xyz("C"), //Position of C
					upperres.xyz( upperres.residue_connect_atom_index(res.residue_connection_conn_id(2)) )
				)
			);

	//printf("rsd %lu phi=%.3f psi=%.3f\n", res.seqpos(), phi, psi); fflush(stdout); //DELETE ME

	core::chemical::AA ref_aa = res.aa();
	if(res.backbone_aa() != core::chemical::aa_unk) ref_aa = res.backbone_aa(); //If this is a noncanonical that specifies a canonical to use as a Rama template, use the template.

	eval_rama_score_residue( ref_aa, phi, psi, rama, drama_dphi, drama_dpsi );

	return;
}

///////////////////////////////////////////////////////////////////////////////
Real
Ramachandran::eval_rama_score_residue(
	conformation::Residue const & rsd
) const
{
	Real rama, drama_dphi, drama_dpsi;
	eval_rama_score_residue( rsd, rama, drama_dphi, drama_dpsi );
	return rama;
}

///////////////////////////////////////////////////////////////////////////////
void
Ramachandran::eval_rama_score_residue(
	conformation::Residue const & rsd,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

	//assert( pose.residue(res).is_protein() );
	assert( rsd.is_protein() );

	//Determine whether connections are incomplete:
	//bool incomplete_connections = (rsd.has_incomplete_connection(rsd.lower_connect_atom()) || rsd.has_incomplete_connection(rsd.upper_connect_atom()));

	if ( 0.0 == nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)) || 0.0 == nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)) || rsd.is_terminus() || rsd.is_virtual_residue() /*|| incomplete_connections*/) { // begin or end of chain -- don't calculate rama score
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	Real phi=0.0;
	Real psi=0.0;


	if(is_normally_connected(rsd)) { //If this residue is conventionally connected
		phi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1));
		psi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2));
		core::chemical::AA ref_aa = rsd.aa();
		if(rsd.backbone_aa() != core::chemical::aa_unk) ref_aa = rsd.backbone_aa(); //If this is a noncanonical that specifies a canonical to use as a Rama template, use the template.

		eval_rama_score_residue( ref_aa, phi, psi, rama, drama_dphi, drama_dpsi );
	} else { //If this residue is unconventionally connected (should be handled elsewhere)
		//printf("Residue %lu: THIS SHOULD NEVER OCCUR!\n", rsd.seqpos()); fflush(stdout); //DELETE ME -- FOR TESTING ONLY
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	return;
}


///////////////////////////////////////////////////////////////////////////////
///
Real
Ramachandran::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi
) const
{

	Real rama, drama_dphi, drama_dpsi;
	eval_rama_score_residue( res_aa, phi, psi, rama, drama_dphi, drama_dpsi );
	return rama;
}

void
Ramachandran::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace basic::options;
	eval_rama_score_residue(
		option[ OptionKeys::corrections::score::use_bicubic_interpolation ],
		option[ OptionKeys::corrections::score::rama_not_squared ],
		res_aa, phi, psi, rama, drama_dphi, drama_dpsi);
}

///////////////////////////////////////////////////////////////////////////////
///
void
Ramachandran::eval_rama_score_residue(
	bool use_bicubic_interpolation,
	bool rama_not_squared,
	AA const res_aa,
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const
{
	using namespace numeric;

//db
//db secondary structure dependent tables favor helix slightly.
//db only use if have predicted all alpha protein
//db
// rhiju and db: no longer use alpha-specific rama, after
//  tests on 1yrf and other all alpha proteins. 2-8-07

// apl -- removing ss dependence on rama in first implementation of mini
// after reading rhiju and david's comment above.  We will need a structural annotation
// obect (structure.cc, maybe a class SecStruct) at some point.  The question
// remains whether a pose should hold that object and be responsible for its upkeep,
// or whether such an object could be created as needed to sit alongside a pose.
//
//	std::string protein_sstype = get_protein_sstype();
//	int ss_type;
//	if ( use_alpha_rama_flag() && get_protein_sstype() == "a" ) {
//		ss_type = ( ( ss == 'H' ) ? 1 : ( ( ss == 'E' ) ? 2 : 3 ) );
//	} else {
	int ss_type = 3;
//	}

//     do I (?cems/cj/cdb?) want to interpolate probabilities or log probs???
//     currently am interpolating  probs then logging.

	//int const res_aa( rsd.aa() );
	// 	int const res_aa( pose.residue( res ).aa() );

	core::chemical::AA res_aa2 = res_aa;
	core::Real phi2 = phi;
	core::Real psi2 = psi;
	core::Real d_multiplier=1.0; //A multiplier for derivatives: 1.0 for L-amino acids, -1.0 for D-amino acids.

	if(is_canonical_d_aminoacid(res_aa)) { //If this is a D-amino acid, invert phi and psi and use the corresponding L-amino acid for the calculation
		res_aa2 = get_l_equivalent(res_aa);
		phi2 = -phi;
		psi2 = -psi;
		d_multiplier = -1.0;
	}

	if ( use_bicubic_interpolation ) {

		rama = rama_energy_splines_[ res_aa2 ].F(phi2,psi2);
		drama_dphi = d_multiplier*rama_energy_splines_[ res_aa2 ].dFdx(phi2,psi2);
		drama_dpsi = d_multiplier*rama_energy_splines_[ res_aa2 ].dFdy(phi2,psi2);
		//printf("drama_dphi=%.8f\tdrama_dpsi=%.8f\n", drama_dphi, drama_dpsi); //DELETE ME!
		//if(is_canonical_d_aminoacid(res_aa)) { printf("rama = %.4f\n", rama); fflush(stdout); } //DELETE ME!
		return; // temp -- just stop right here
	} else {

		FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
		FArray2A< Real > const rama_for_res( ram_probabil_(1, 1, ss_type, res_aa2), zero_index, zero_index );
		Real interp_p,dp_dphi,dp_dpsi;

		using namespace numeric::interpolation::periodic_range::half;
		interp_p = bilinearly_interpolated( phi2, psi2, binw_, n_phi_, rama_for_res, dp_dphi, dp_dpsi );

		if ( interp_p > 0.0 ) {
			rama = ram_entropy_(ss_type, res_aa2 ) - std::log( static_cast< double >( interp_p ) );
			double const interp_p_inv_neg = -1.0 / interp_p;
			drama_dphi = interp_p_inv_neg * d_multiplier * dp_dphi;
			drama_dpsi = interp_p_inv_neg * d_multiplier * dp_dpsi;
		} else {
			//if ( runlevel > silent ) { //apl fix this
			//	std::cout << "rama prob = 0. in eval_rama_score_residue!" << std::endl;
			//	std::cout << "phi" << SS( phi ) << " psi" << SS( psi ) <<
			//	 " ss " << SS( ss ) << std::endl;
			//}
			drama_dphi = 0.0;
			drama_dpsi = 0.0;
			rama = 20.0;
		}

		if ( ! rama_not_squared ) {
			if ( rama > 1.0 ) {
				Real const rama_squared = rama * rama;
				if ( rama_squared > 20.0 ) {
					////  limit the score, but give the true derivative
					////   as guidance out of the flat section of map
					drama_dphi = 0.0;
					drama_dpsi = 0.0;
					rama = 20.0;
				} else {
					drama_dphi = 2 * rama * drama_dphi;
					drama_dpsi = 2 * rama * drama_dpsi;
					rama = rama_squared;
				}
			}
		}
	}
}



///////////////////////////////////////////////////////////////////////////////
void Ramachandran::eval_procheck_rama(
	Pose const & /*pose*/,
	Real & /*favorable*/,
	Real & /*allowed*/,
	Real & /*generous*/
) const
{}

bool
Ramachandran::is_normally_connected (
	conformation::Residue const & res
) const {
	if(!basic::options::option[basic::options::OptionKeys::score::rama_score_nonstandard_connections]) return true;
	if (res.connection_incomplete(1) || res.connection_incomplete(2)) return true;
	return ( (res.residue_connection_partner(1) == res.seqpos()-1) && (res.residue_connection_partner(2) == res.seqpos()+1) );
}

void
Ramachandran::read_rama_map_file (
	utility::io::izstream * iunit
) {
	int aa_num,phi_bin,psi_bin,ss_type;
	Real check,min_prob,max_prob;
	double entropy;
	char line[60];
	int scan_count;
	float pval, eval; // vars for sscanf float I/O

	for ( int i = 1; i <= n_aa_ ; ++i ) {
		for ( int ii = 1; ii <= 3; ++ii ) {
			entropy = 0.0;
			check = 0.0;
			min_prob = 1e36;
			max_prob = -min_prob;
			for ( int j = 1; j <= 36; ++j ) {
				for ( int k = 1; k <= 36; ++k ) {
					iunit->getline( line, 60 );
					if ( iunit->eof() ) {
						return;
					} else if ( iunit->fail() ) { // Clear and continue: NO ERROR DETECTION
						iunit->clear();
					}
					std::sscanf( line, "%5d", &aa_num );
					std::sscanf( line+6, "%5d", &ss_type );
					std::sscanf( line+12, "%5d", &phi_bin );
					std::sscanf( line+18, "%5d", &psi_bin );
					std::sscanf( line+24, "%5d", &ram_counts_(j,k,ii,i) );
					std::sscanf( line+30, "%12f", &pval );
					ram_probabil_(j,k,ii,i) = pval;
					scan_count = std::sscanf( line+43, "%12f", &eval );
					ram_energ_(j,k,ii,i) = eval;

					if ( scan_count == EOF ) continue; // Read problem: NO ERROR DETECTION

// This is the Slick & Slow (S&S) stream-based method that is too slow for large
// files like this one, at least under the GCC 3.3.1 stream implementation.
// It should be retried on future releases and target compilers because there is
// no reason it cannot be competitive with good optimization and inlining.
// If this is used the <cstdio> can be removed.
//
//					iunit >> bite( 5, aa_num ) >> skip( 1 ) >>
//					 bite( 5, ss_type ) >> skip( 1 ) >>
//					 bite( 5, phi_bin ) >> skip( 1 ) >>
//					 bite( 5, psi_bin ) >> skip( 1 ) >>
//					 bite( 5, ram_counts(j,k,ii,i) ) >> skip( 1 ) >>
//					 bite( 12, ram_probabil(j,k,ii,i) ) >> skip( 1 ) >>
//					 bite( 12, ram_energ(j,k,ii,i) ) >> skip;
//					if ( iunit.eof() ) {
//						goto L100;
//					} else if ( iunit.fail() ) { // Clear and continue: NO ERROR DETECTION
//						iunit.clear();
//						iunit >> skip;
//					}

					check += ram_probabil_(j,k,ii,i);
					entropy += ram_probabil_(j,k,ii,i) *
					 std::log( static_cast< double >( ram_probabil_(j,k,ii,i) ) );
					min_prob = std::min(ram_probabil_(j,k,ii,i),min_prob);
					max_prob = std::max(ram_probabil_(j,k,ii,i),max_prob);
				}
			}
			ram_entropy_(ii,i) = entropy;
		}
//cj		std::cout << SS( check ) << SS( std::log(min_prob) ) <<
//cj		 SS( std::log(max_prob) ) << SS( entropy ) << std::endl;
	}

}

void
Ramachandran::read_rama(
	std::string const & rama_map_filename,
	bool use_bicubic_interpolation
) {

	utility::io::izstream  iunit;

  // search in the local directory first
  iunit.open( rama_map_filename );

  if ( !iunit.good() ) {
    iunit.close();
    if(!basic::database::open( iunit, rama_map_filename )){
			std::stringstream err_msg;
			err_msg << "Unable to open Ramachandran map '" << rama_map_filename << "'.";
			utility_exit_with_message(err_msg.str());
		}
  }

//cj      std::cout << "index" << "aa" << "ramachandran entropy" << std::endl;
//KMa add_phospho_ser 2006-01
	read_rama_map_file (&iunit);

	iunit.close();
	iunit.clear();

	if ( use_bicubic_interpolation ) {
		using namespace numeric;
		using namespace numeric::interpolation::spline;
		rama_energy_splines_.resize( chemical::num_canonical_aas );
		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			BicubicSpline ramaEspline;
			MathMatrix< Real > energy_vals( 36, 36 );
			for ( Size jj = 0; jj < 36; ++jj ) {
				for ( Size kk = 0; kk < 36; ++kk ) {
					energy_vals( jj, kk ) = -std::log( ram_probabil_(jj+1,kk+1,3,ii )) + ram_entropy_(3,ii) ;
				}
			}
			BorderFlag periodic_boundary[2] = { e_Periodic, e_Periodic };
			Real start_vals[2] = {5.0, 5.0}; // grid is shifted by five degrees.
			Real deltas[2] = {10.0, 10.0}; // grid is 10 degrees wide
			bool lincont[2] = {false,false}; //meaningless argument for a bicubic spline with periodic boundary conditions
			std::pair< Real, Real > unused[2];
			unused[0] = std::make_pair( 0.0, 0.0 );
			unused[1] = std::make_pair( 0.0, 0.0 );
			ramaEspline.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
			rama_energy_splines_[ ii ] = ramaEspline;
		}
	}
//cj      std::cout << "========================================" << std::endl;
}



}
}
