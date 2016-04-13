// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/bin_transitions/BinTransitionData.cc
/// @brief  Class member functions for BinTransitionData class.
/// @details This class stores data associated with transitions from one mainchain torsion bin to another (e.g. ABEGO bins, OO-ABBA bins, etc.)
/// for ONE specific type of transition (ith residue has certain properties, i+1st residue has certain other properties).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <core/scoring/bin_transitions/BinTransitionData.hh>

// File I/O
#include <basic/database/open.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

//Random
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

// Other Headers
#include <basic/Tracer.hh>

namespace core {
namespace scoring {
namespace bin_transitions {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.bin_transitions.BinTransitionData" );

/// @brief Default constructor for BinTransitionData
///
BinTransitionData::BinTransitionData(): //TODO -- initialize variables here:
	n_mainchain_torsions_i_(0),
	n_mainchain_torsions_iplus1_(0),
	n_bins_i_(0),
	binranges_i_(),
	binnames_i_(),
	n_bins_iplus1_(0),
	binranges_iplus1_(),
	binnames_iplus1_(),
	matrix_initialized_(false),
	matrix_finalized_(false),
	probability_matrix_(),
	properties_i_(),
	prohibited_properties_i_(),
	properties_iplus1_(),
	prohibited_properties_iplus1_(),
	res_identities_i_(),
	prohibited_res_identities_i_(),
	res_identities_iplus1_(),
	prohibited_res_identities_iplus1_(),
	binsums_i_(),
	binsums_iplus1_(),
	total_counts_(0.0),
	binsums_i_cdf_(),
	binsums_iplus1_cdf_(),
	subbin_type_i_(BTSB_NONE),
	subbin_type_iplus1_(BTSB_NONE),
	subbin_cdf_i_(),
	subbin_ranges_i_(),
	subbin_cdf_iplus1_(),
	subbin_ranges_iplus1_(),
	bin_iplus1_cdf_given_i_(),
	bin_i_cdf_given_iplus1_()
{}

/// @brief Copy constructor for BinTransitionData
///
BinTransitionData::BinTransitionData( BinTransitionData const &src ): //TODO -- copy variables here:
	utility::pointer::ReferenceCount(),
	n_mainchain_torsions_i_(src.n_mainchain_torsions_i_),
	n_mainchain_torsions_iplus1_(src.n_mainchain_torsions_iplus1_),
	n_bins_i_(src.n_bins_i_),
	binranges_i_( src.binranges_i_ ),
	binnames_i_(src.binnames_i_),
	n_bins_iplus1_(src.n_bins_iplus1_),
	binranges_iplus1_( src.binranges_iplus1_ ),
	binnames_iplus1_(src.binnames_iplus1_),
	matrix_initialized_(src.matrix_initialized_),
	matrix_finalized_(src.matrix_finalized_),
	probability_matrix_(src.probability_matrix_),
	properties_i_( src.properties_i_ ),
	prohibited_properties_i_( src.prohibited_properties_i_ ),
	properties_iplus1_( src.properties_iplus1_ ),
	prohibited_properties_iplus1_( src.prohibited_properties_iplus1_ ),
	res_identities_i_( src.res_identities_i_ ),
	prohibited_res_identities_i_( src.prohibited_res_identities_i_ ),
	res_identities_iplus1_( src.res_identities_iplus1_ ),
	prohibited_res_identities_iplus1_( src.prohibited_res_identities_iplus1_ ),
	binsums_i_( src.binsums_i_ ),
	binsums_iplus1_( src.binsums_iplus1_ ),
	total_counts_( src.total_counts_ ),
	binsums_i_cdf_( src.binsums_i_cdf_ ),
	binsums_iplus1_cdf_( src.binsums_iplus1_cdf_ ),
	subbin_type_i_( src.subbin_type_i_ ),
	subbin_type_iplus1_( src.subbin_type_iplus1_ ),
	subbin_cdf_i_( src.subbin_cdf_i_ ),
	subbin_ranges_i_( src.subbin_ranges_i_ ),
	subbin_cdf_iplus1_( src.subbin_cdf_iplus1_ ),
	subbin_ranges_iplus1_( src.subbin_ranges_iplus1_ ),
	bin_iplus1_cdf_given_i_( src.bin_iplus1_cdf_given_i_ ),
	bin_i_cdf_given_iplus1_( src.bin_i_cdf_given_iplus1_ )
{}

/// @brief Default destructor for BinTransitionData
///
BinTransitionData::~BinTransitionData() {}

/// @brief Clone operation for BinTransitionData.
/// @details Returns an owning pointer to a copy of this object.
BinTransitionDataOP BinTransitionData::clone() const
{ return BinTransitionDataOP( new BinTransitionData( *this ) ); }

/// @brief Set up the sub-bins and their cumulative probability distributions (if appropriate).
///
void BinTransitionData::set_up_subbins() {
	if ( !matrix_initialized() || matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::set_up_subbins(): The matrix must be initialized but not finalized before calling this function.\n" );

	/////////// FOR RESIDUE I ///////////
	if ( subbin_type_i() == BTSB_L_AA || subbin_type_i() == BTSB_D_AA || subbin_type_i() == BTSB_L_PRO || subbin_type_i() == BTSB_D_PRO || subbin_type_i() == BTSB_GLY ) {
		runtime_assert_string_msg( n_mainchain_torsions_i()==3, "In BinTransitionData::set_up_subbins(): The " + get_subbin_type_name(subbin_type_i()) + " sub-bin type requires that there be 3 mainchain torsions." );

		core::Real phipsi_multiplier(1.0);
		core::chemical::AA aatype = core::chemical::aa_ala;
		if ( subbin_type_i() == BTSB_L_AA ) {
			aatype=core::chemical::aa_ala;
		} else if ( subbin_type_i() == BTSB_D_AA ) {
			aatype=core::chemical::aa_ala;
			phipsi_multiplier=-1.0;
		} else if ( subbin_type_i() == BTSB_L_PRO ) {
			aatype=core::chemical::aa_pro;
		} else if ( subbin_type_i() == BTSB_D_PRO ) {
			aatype=core::chemical::aa_pro;
			phipsi_multiplier=-1.0;
		} else if ( subbin_type_i() == BTSB_GLY ) {
			aatype=core::chemical::aa_gly;
		}

		// A Ramachandran object used for torsion bins.
		core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		core::Size const n_phi_bins( rama.n_phi_bins() ); // The number of phi bins in the rama object
		core::Size const n_psi_bins( rama.n_phi_bins() ); // The number of psi bins in the rama object

		subbin_cdf_i_.clear();
		subbin_ranges_i_.clear();
		for ( core::Size bin_i=1, binimax=n_bins_i(); bin_i<=binimax; ++bin_i ) { //Loop through bins for the ith residue.
			utility::vector1 /*sub-bins*/ < core::Real > curbin_cdf; //The cumulative distribution function for the sub-bins of the current bin
			utility::vector1 /*sub-bins*/ < utility::vector1 /*mainchain torsions*/ < std::pair /*start and end vals*/ < core::Real, core::Real > > > curbin_subbin_defs; //The definitions of the sub-bins of the current bin.
			core::Real const phimin (binranges_i_[bin_i][1].first);
			core::Real const phimax (binranges_i_[bin_i][1].second);
			core::Real const psimin (binranges_i_[bin_i][2].first);
			core::Real const psimax (binranges_i_[bin_i][2].second);
			core::Real const omegamin (binranges_i_[bin_i][3].first);
			core::Real const omegamax (binranges_i_[bin_i][3].second);

			bool cis ( in_bin( omegamin, omegamax, 0.0 ) );

			core::Real phi=-180.0;
			core::Real psi=-180.0;
			core::Real const phi_increment = 360.0/static_cast<core::Real>(n_phi_bins);
			core::Real const psi_increment = 360.0/static_cast<core::Real>(n_psi_bins);

			for ( core::Size i=1; i<=n_phi_bins; ++i ) { //Loop though all phi bins
				psi=-180.0;
				for ( core::Size j=1; j<=n_psi_bins; ++j ) { //Loop though all psi bins
					//if(TR.visible()) TR << binnames_i_[bin_i] << "\tphi=" << phi << "\tpsi=" << psi << "\tphimin=" << phimin << "\tpsimin=" << psimin; //DELETE ME
					if ( !(
							( in_bin(phimin, phimax, phi) || in_bin(phimin, phimax, phi+phi_increment) )
							&&
							( in_bin(psimin, psimax, psi) || in_bin(psimin, psimax, psi+psi_increment)) )
							) {
						//if(TR.visible()) TR << std::endl; //DELETE ME
						psi+=psi_increment;
						continue; //We're not in a sub-bin of the current bin, so do nothing.
					}
					//if(TR.visible()) TR << "\tIN BIN" << std::endl; //DELETE ME

					curbin_cdf.push_back(rama.rama_probability( aatype, phi*phipsi_multiplier, psi*phipsi_multiplier )); //We've found a subbin in the current bin.  Store the rama probability (inverted for D-amino acids).

					//Set the torsion ranges for the current sub-bin
					utility::vector1 < std::pair <core::Real,core::Real> > tors_ranges;
					tors_ranges.resize(3, std::make_pair(0.0,0.0));
					tors_ranges[1].first=phi;
					tors_ranges[1].second=set_in_range(phi+phi_increment);
					tors_ranges[2].first=psi;
					tors_ranges[2].second=set_in_range(psi+psi_increment);
					tors_ranges[3].first=(cis ? -5 : 175); //If this is cis, the sub-bins are from omega = -5 to 5; if this is trans, the sub-bins are from omega = 175 to -175.
					tors_ranges[3].second=(cis ? 5 : -175);

					trim_subbin_edges_and_rescale_subbin( curbin_cdf[curbin_cdf.size()], tors_ranges[1], phimin, phimax );
					trim_subbin_edges_and_rescale_subbin( curbin_cdf[curbin_cdf.size()], tors_ranges[2], psimin, psimax );

					curbin_subbin_defs.push_back(tors_ranges);
					if ( curbin_cdf.size() > 1 ) curbin_cdf[curbin_cdf.size()]+=curbin_cdf[curbin_cdf.size()-1]; //Add the previous bin's probability, to make this a cumulative distribution function.

					psi+=psi_increment;
				}
				phi+=phi_increment;
			}

			//Re-normalize the current cdf:
			for ( core::Size i=1, imax=curbin_cdf.size(); i<=imax; ++i ) {
				curbin_cdf[i] = curbin_cdf[i] / curbin_cdf[imax];
			}

			//Add the current cdf to the cdf list for all bins:
			subbin_cdf_i_.push_back( curbin_cdf );
			//Add the current list of sub-bin ranges to the list for all bins:
			subbin_ranges_i_.push_back( curbin_subbin_defs );
		} //Looping through bins for residue i.

		//Note: Rosetta's cumulative distribution (cdf) functions are strange.
		//They assume a cdf where P(cell i) = cdf[i+1]-cdf[i], rather than
		//P(cell i) = cdf[i]-cdf[i-1].  This means I need to adjust, here:
		for ( core::Size j=1, jmax=subbin_cdf_i_.size(); j<=jmax; ++j ) {
			for ( core::Size i=subbin_cdf_i_[j].size(); i>1; --i ) {
				subbin_cdf_i_[j][i]=subbin_cdf_i_[j][i-1];
			}
			if ( subbin_cdf_i_[j].size()>0 ) subbin_cdf_i_[j][1]=0;
		}
	}
	/////////// END FOR RESIDUE I ///////////

	/////////// FOR RESIDUE I+1 ///////////
	if ( subbin_type_iplus1() == BTSB_L_AA || subbin_type_iplus1() == BTSB_D_AA || subbin_type_iplus1() == BTSB_L_PRO || subbin_type_iplus1() == BTSB_D_PRO || subbin_type_iplus1() == BTSB_GLY ) {
		runtime_assert_string_msg( n_mainchain_torsions_iplus1()==3, "In BinTransitionData::set_up_subbins(): The " + get_subbin_type_name(subbin_type_iplus1()) + " sub-bin type requires that there be 3 mainchain torsions." );

		core::Real phipsi_multiplier(1.0);
		core::chemical::AA aatype = core::chemical::aa_ala;
		if ( subbin_type_iplus1() == BTSB_L_AA ) {
			aatype=core::chemical::aa_ala;
		} else if ( subbin_type_iplus1() == BTSB_D_AA ) {
			aatype=core::chemical::aa_ala;
			phipsi_multiplier=-1.0;
		} else if ( subbin_type_iplus1() == BTSB_L_PRO ) {
			aatype=core::chemical::aa_pro;
		} else if ( subbin_type_iplus1() == BTSB_D_PRO ) {
			aatype=core::chemical::aa_pro;
			phipsi_multiplier=-1.0;
		} else if ( subbin_type_iplus1() == BTSB_GLY ) {
			aatype=core::chemical::aa_gly;
		}

		// A Ramachandran object used for torsion bins.
		core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		core::Size const n_phi_bins( rama.n_phi_bins() ); // The number of phi bins in the rama object
		core::Size const n_psi_bins( rama.n_phi_bins() ); // The number of psi bins in the rama object

		subbin_cdf_iplus1_.clear();
		subbin_ranges_iplus1_.clear();
		for ( core::Size bin_iplus1=1, binimax=n_bins_iplus1(); bin_iplus1<=binimax; ++bin_iplus1 ) { //Loop through bins for the ith residue.
			utility::vector1 /*sub-bins*/ < core::Real > curbin_cdf; //The cumulative distribution function for the sub-bins of the current bin
			utility::vector1 /*sub-bins*/ < utility::vector1 /*mainchain torsions*/ < std::pair /*start and end vals*/ < core::Real, core::Real > > > curbin_subbin_defs; //The definitions of the sub-bins of the current bin.
			core::Real const phimin (binranges_iplus1_[bin_iplus1][1].first);
			core::Real const phimax (binranges_iplus1_[bin_iplus1][1].second);
			core::Real const psimin (binranges_iplus1_[bin_iplus1][2].first);
			core::Real const psimax (binranges_iplus1_[bin_iplus1][2].second);
			core::Real const omegamin (binranges_iplus1_[bin_iplus1][3].first);
			core::Real const omegamax (binranges_iplus1_[bin_iplus1][3].second);

			bool cis ( in_bin( omegamin, omegamax, 0.0 ) );

			core::Real phi=-180.0;
			core::Real psi=-180.0;
			core::Real const phi_increment = 360.0/static_cast<core::Real>(n_phi_bins);
			core::Real const psi_increment = 360.0/static_cast<core::Real>(n_psi_bins);

			for ( core::Size i=1; i<=n_phi_bins; ++i ) { //Loop though all phi bins
				psi=-180.0;
				for ( core::Size j=1; j<=n_psi_bins; ++j ) { //Loop though all psi bins
					//if(TR.visible()) TR << binnames_i_[bin_i] << "\tphi=" << phi << "\tpsi=" << psi << "\tphimin=" << phimin << "\tpsimin=" << psimin; //DELETE ME
					if ( !(
							( in_bin(phimin, phimax, phi) || in_bin(phimin, phimax, phi+phi_increment) )
							&&
							( in_bin(psimin, psimax, psi) || in_bin(psimin, psimax, psi+psi_increment)) )
							) {
						//if(TR.visible()) TR << std::endl; //DELETE ME
						psi+=psi_increment;
						continue; //We're not in a sub-bin of the current bin, so do nothing.
					}
					//if(TR.visible()) TR << "\tIN BIN" << std::endl; //DELETE ME

					curbin_cdf.push_back(rama.rama_probability( aatype, phi*phipsi_multiplier, psi*phipsi_multiplier )); //We've found a subbin in the current bin.  Store the rama probability (inverted for D-amino acids).

					//Set the torsion ranges for the current sub-bin
					utility::vector1 < std::pair <core::Real,core::Real> > tors_ranges;
					tors_ranges.resize(3, std::make_pair(0.0,0.0));
					tors_ranges[1].first=phi;
					tors_ranges[1].second=set_in_range(phi+phi_increment);
					tors_ranges[2].first=psi;
					tors_ranges[2].second=set_in_range(psi+psi_increment);
					tors_ranges[3].first=(cis ? -5 : 175); //If this is cis, the sub-bins are from omega = -5 to 5; if this is trans, the sub-bins are from omega = 175 to -175.
					tors_ranges[3].second=(cis ? 5 : -175);

					trim_subbin_edges_and_rescale_subbin( curbin_cdf[curbin_cdf.size()], tors_ranges[1], phimin, phimax );
					trim_subbin_edges_and_rescale_subbin( curbin_cdf[curbin_cdf.size()], tors_ranges[2], psimin, psimax );

					curbin_subbin_defs.push_back(tors_ranges);
					if ( curbin_cdf.size() > 1 ) curbin_cdf[curbin_cdf.size()]+=curbin_cdf[curbin_cdf.size()-1]; //Add the previous bin's probability, to make this a cumulative distribution function.

					psi+=psi_increment;
				}
				phi+=phi_increment;
			}

			//Re-normalize the current cdf:
			for ( core::Size i=1, imax=curbin_cdf.size(); i<=imax; ++i ) {
				curbin_cdf[i] = curbin_cdf[i] / curbin_cdf[imax];
			}

			//Add the current cdf to the cdf list for all bins:
			subbin_cdf_iplus1_.push_back( curbin_cdf );
			//Add the current list of sub-bin ranges to the list for all bins:
			subbin_ranges_iplus1_.push_back( curbin_subbin_defs );
		} //Looping through bins for residue i.

		//Note: Rosetta's cumulative distribution (cdf) functions are strange.
		//They assume a cdf where P(cell i) = cdf[i+1]-cdf[i], rather than
		//P(cell i) = cdf[i]-cdf[i-1].  This means I need to adjust, here:
		for ( core::Size j=1, jmax=subbin_cdf_iplus1_.size(); j<=jmax; ++j ) {
			for ( core::Size i=subbin_cdf_iplus1_[j].size(); i>1; --i ) {
				subbin_cdf_iplus1_[j][i]=subbin_cdf_iplus1_[j][i-1];
			}
			if ( subbin_cdf_iplus1_[j].size()>0 ) subbin_cdf_iplus1_[j][1]=0;
		}
	}
	/////////// END FOR RESIDUE I+1 ///////////

	return;
} //set_up_subbins

/// @brief Set up the cumulative probability distributions for each bin at position i+1 given a bin at position i, and vice versa.
///
void BinTransitionData::set_up_bin_cdfs() {
	if ( !matrix_initialized() || matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::set_up_bin_cdfs(): The matrix must be initialized but not finalized before calling this function.\n" );

	core::Size const nbins_i( n_bins_i() ); //Number of bins for the ith residue
	core::Size const nbins_iplus1( n_bins_iplus1() ); //Number of bins for the i+1st residue

	//Set up cdf for i+1 bins, given i.
	bin_iplus1_cdf_given_i_.clear(); //Reset the cdfs.
	bin_iplus1_cdf_given_i_.resize(nbins_i);
	for ( core::Size n=1; n<=nbins_i; ++n ) { //Loop through each bin for residue i
		bin_iplus1_cdf_given_i_[n].clear();
		bin_iplus1_cdf_given_i_[n].resize( nbins_iplus1, 0.0 );
		for ( core::Size m=1; m<=nbins_iplus1; ++m ) { //Loop through each bin for residue i+1
			bin_iplus1_cdf_given_i_[n][m] = probability_matrix_[n][m] + ( m>1 ? bin_iplus1_cdf_given_i_[n][m-1] : 0.0 );
		}
		for ( core::Size m=1; m<=nbins_iplus1; ++m ) { //Loop through again to normalize
			bin_iplus1_cdf_given_i_[n][m] /= bin_iplus1_cdf_given_i_[n][nbins_iplus1];
		}
	}

	//Set up cdf for i bins, given i+1.
	bin_i_cdf_given_iplus1_.clear(); //Reset the cdfs.
	bin_i_cdf_given_iplus1_.resize(nbins_iplus1);
	for ( core::Size m=1; m<=nbins_iplus1; ++m ) { //Loop through each bin for residue i
		bin_i_cdf_given_iplus1_[m].clear();
		bin_i_cdf_given_iplus1_[m].resize( nbins_i, 0.0 );
		for ( core::Size n=1; n<=nbins_i; ++n ) { //Loop through each bin for residue i+1
			bin_i_cdf_given_iplus1_[m][n] = probability_matrix_[n][m] + ( n>1 ? bin_i_cdf_given_iplus1_[m][n-1] : 0.0 );
		}
		for ( core::Size n=1; n<=nbins_i; ++n ) { //Loop through again to normalize
			bin_i_cdf_given_iplus1_[m][n] /= bin_i_cdf_given_iplus1_[m][nbins_i];
		}
	}

	//Note: Rosetta's cumulative distribution (cdf) functions are strange.
	//They assume a cdf where P(cell i) = cdf[i+1]-cdf[i], rather than
	//P(cell i) = cdf[i]-cdf[i-1].  This means I need to adjust, here:
	for ( core::Size i=1, imax=nbins_i; i<=imax; ++i ) {
		for ( core::Size j=nbins_iplus1; j>1; --j ) {
			bin_iplus1_cdf_given_i_[i][j] = bin_iplus1_cdf_given_i_[i][j-1];
		}
		if ( nbins_iplus1 > 0 ) bin_iplus1_cdf_given_i_[i][1]=0;
	}
	for ( core::Size i=1, imax=nbins_iplus1; i<=imax; ++i ) {
		for ( core::Size j=nbins_i; j>1; --j ) {
			bin_i_cdf_given_iplus1_[i][j] = bin_i_cdf_given_iplus1_[i][j-1];
		}
		if ( nbins_i > 0 ) bin_i_cdf_given_iplus1_[i][1] = 0;
	}

	return;
} //set_up_bin_cdfs

/// @brief Ensure that a sub-bin doesn't exceed the bounds of a bin.  If it does, scale its back, and reduce its contribution
/// to the cumulative distribution function proportionately.
void BinTransitionData::trim_subbin_edges_and_rescale_subbin(
	core::Real &curbin_cdf_val,
	std::pair <core::Real,core::Real> &tors_range,
	core::Real const &phimin,
	core::Real const &phimax
) const {

	core::Real tors_range_out_start(tors_range.first);
	core::Real tors_range_out_end(tors_range.second);

	{ //Scope 1: adjusting the start of the sub-bin.
		core::Real const binstart(phimin);
		core::Real const binend(phimin < phimax ? phimax : phimax+360.0 );

		core::Real substart(tors_range.first);
		core::Real subend(tors_range.first < tors_range.second ? tors_range.second : tors_range.second + 360.0);
		if ( binend >= 180.0 && subend < binstart ) {
			substart+=360.0; subend+=360.0;
		}

		if ( substart < binstart ) {
			curbin_cdf_val *= (subend - binstart) / (subend - substart);
			tors_range_out_start = phimin;
			//if(TR.visible()) TR << "Trimming start of sub-bin." << std::endl; //DELETE ME
		}

		if ( subend > binend ) {
			curbin_cdf_val *= (binend - substart) / (subend - substart);
			tors_range_out_end = phimax;
			//if(TR.visible()) TR << "Trimming end of sub-bin." << std::endl; //DELETE ME
		}
	}

	tors_range.first = tors_range_out_start;
	tors_range.second = tors_range_out_end;

	return;
} //trim_subbin_edges_and_rescale_subbin

/// @brief Given a sub-bin type enum, return its name as it would appear in
/// a residue params file.
std::string BinTransitionData::get_subbin_type_name (BTSB_SUBBIN_TYPE const type) const
{
	std::string returnstring ("");
	switch(type) {
	case BTSB_NONE :
		returnstring="NONE";
		break;
	case BTSB_L_AA :
		returnstring="L_AA";
		break;
	case BTSB_D_AA :
		returnstring="D_AA";
		break;
	case BTSB_L_PRO :
		returnstring="L_PRO";
		break;
	case BTSB_D_PRO :
		returnstring="D_PRO";
		break;
	case BTSB_GLY :
		returnstring="GLY";
		break;
	default :
		returnstring="UNKNOWN_TYPE";
		break;
	}
	return returnstring;
} //get_subbin_type_name

/// @brief Given a sub-bin type name, return the sub-bin type enum value.
/// @details Returns BTSB_UNKNOWN if not identifiable.
BTSB_SUBBIN_TYPE BinTransitionData::get_subbin_type_from_name( std::string const &name ) const
{
	for ( core::Size i=1; i<static_cast<core::Size>(BTSB_UNKNOWN); ++i ) {
		if ( get_subbin_type_name(static_cast<BTSB_SUBBIN_TYPE>(i) ) == name ) return static_cast< BTSB_SUBBIN_TYPE >(i);
	}
	return BTSB_UNKNOWN;
}

/// @brief Given a property enum value, return the property name as it would appear in a
/// residue params file.
std::string BinTransitionData::get_property_effect_name( BT_PROPERTIES const property) const
{
	std::string returnstring("");

	switch(property){
	case BT_PROTEIN :
		returnstring="PROTEIN";
		break;
	case BT_L_AA :
		returnstring="L_AA";
		break;
	case BT_D_AA :
		returnstring="D_AA";
		break;
	case BT_ALPHA_AA :
		returnstring="ALPHA_AA";
		break;
	case BT_BETA_AA :
		returnstring="BETA_AA";
		break;
	case BT_POLAR :
		returnstring="POLAR";
		break;
	case BT_METALBINDING :
		returnstring="METALBINDING";
		break;
	case BT_CHARGED :
		returnstring="CHARGED";
		break;
	case BT_AROMATIC :
		returnstring="AROMATIC";
		break;
	case BT_DISULFIDE_BONDED :
		returnstring="DISULFIDE_BONDED";
		break;
	case BT_SIDECHAIN_THIOL :
		returnstring="SIDECHAIN_THIOL";
		break;
	case BT_CYCLIC :
		returnstring="CYCLIC";
		break;
	default :
		returnstring="UNKNOWN_PROPERTY";
		break;
	}

	return returnstring;
} //get_property_effect_name()

/// @brief Given a property name, return the property enum value.
///
BT_PROPERTIES BinTransitionData::get_property_from_name( std::string const &name ) const
{
	for ( core::Size i=1; i<static_cast<core::Size>(BT_UNKNOWN_PROPERTY); ++i ) {
		if ( get_property_effect_name(static_cast<BT_PROPERTIES>(i) ) == name ) return static_cast< BT_PROPERTIES >(i);
	}
	return BT_UNKNOWN_PROPERTY;
} //get_property_from_name

/// @brief Given a property, check whether a residue has that property.
/// @details This function provides the link to the suitable property lookup function in ResidueType.
bool BinTransitionData::has_property( BT_PROPERTIES const property, core::conformation::Residue const &rsd ) const {
	switch(property){
	case BT_PROTEIN :
		return rsd.type().is_protein();
		break;
	case BT_L_AA :
		return rsd.type().is_l_aa();
		break;
	case BT_D_AA :
		return rsd.type().is_d_aa();
		break;
	case BT_ALPHA_AA :
		return rsd.type().is_alpha_aa();
		break;
	case BT_BETA_AA :
		return rsd.type().is_beta_aa();
		break;
	case BT_POLAR :
		return rsd.type().is_polar();
		break;
	case BT_METALBINDING :
		return rsd.type().is_metalbinding();
		break;
	case BT_CHARGED :
		return rsd.type().is_charged();
		break;
	case BT_AROMATIC :
		return rsd.type().is_aromatic();
		break;
	case BT_DISULFIDE_BONDED :
		return rsd.type().is_disulfide_bonded();
		break;
	case BT_SIDECHAIN_THIOL :
		return rsd.type().is_sidechain_thiol();
		break;
	case BT_CYCLIC :
		return rsd.type().is_cyclic();
		break;
	default :
		utility_exit_with_message("Unknown property passed to BinTransitionData::has_property()!  Failing!\n");
		break;
	}
	return false; //Should never reach this point.
} //has_property

/// @brief Given a residue, check whether its properties match the required and prohibited properties lists for residue i,
/// and whether its name matches the required and prohibited names list for residue i.
bool BinTransitionData::criteria_match_i( core::conformation::Residue const &rsd ) const
{
	if ( !matrix_initialized() || !matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::criteria_match_i(): The transition probability matrix has not yet been loaded!  Failing.\n" );

	if ( properties_i_.size()>0 ) { //If there are required properties
		for ( core::Size i=1, imax=properties_i_.size(); i<=imax; ++i ) { //Loop through required properties
			if ( !has_property(properties_i_[i], rsd) ) return false; //If the residue lacks this property, return false.
		}
	}

	if ( prohibited_properties_i_.size()>0 ) { //If there are prohibited properties
		for ( core::Size i=1, imax=prohibited_properties_i_.size(); i<=imax; ++i ) { //Loop through the prohibited properties
			if ( has_property(prohibited_properties_i_[i], rsd) ) return false; //If the residue has any prohibited property, return false.
		}
	}

	if ( res_identities_i_.size()>0 ) { //If there are required residue identities
		bool match(false);
		for ( core::Size i=1, imax=res_identities_i_.size(); i<=imax; ++i ) { //Loop through required res_identities
			if ( rsd.name3() == res_identities_i_[i] ) { //If the residue matches ANY of the required residue identities, we're good.
				match=true;
				break;
			}
		}
		if ( !match ) return false; //Otherwise, if NO required residue matched this residue, return false.
	}

	if ( prohibited_res_identities_i_.size()>0 ) { //If there are prohibited res_identities
		for ( core::Size i=1, imax=prohibited_res_identities_i_.size(); i<=imax; ++i ) { //Loop through the prohibited res_identities
			if ( rsd.name3() == prohibited_res_identities_i_[i] ) return false; //If the residue matches any prohibited residue identity, return false.
		}
	}

	//if(TR.visible()) TR << "Residue " << rsd.name3() << rsd.seqpos() << " matches!" << std::endl; //DELETE ME!

	return true; //Nothing has failed at this point, so it's a match!
} //criteria_match_i

/// @brief Given a residue, check whether its properties match the required and prohibited properties lists for residue i+1,
/// and whether its name matches the required and prohibited names list for residue i+1.
bool BinTransitionData::criteria_match_iplus1( core::conformation::Residue const &rsd ) const
{
	if ( !matrix_initialized() || !matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::criteria_match_iplus1(): The transition probability matrix has not yet been loaded!  Failing.\n" );

	if ( properties_iplus1_.size()>0 ) { //If there are required properties
		for ( core::Size i=1, imax=properties_iplus1_.size(); i<=imax; ++i ) { //Loop through required properties
			if ( !has_property(properties_iplus1_[i], rsd) ) return false; //If the residue lacks this property, return false.
		}
	}

	if ( prohibited_properties_iplus1_.size()>0 ) { //If there are prohibited properties
		for ( core::Size i=1, imax=prohibited_properties_iplus1_.size(); i<=imax; ++i ) { //Loop through the prohibited properties
			if ( has_property(prohibited_properties_iplus1_[i], rsd) ) return false; //If the residue has any prohibited property, return false.
		}
	}

	if ( res_identities_iplus1_.size()>0 ) { //If there are required residue identities
		bool match(false);
		for ( core::Size i=1, imax=res_identities_iplus1_.size(); i<=imax; ++i ) { //Loop through required res_identities
			if ( rsd.name3() == res_identities_iplus1_[i] ) { //If the residue matches ANY of the required residue identities, we're good.
				match=true;
				break;
			}
		}
		if ( !match ) return false; //Otherwise, if NO required residue matched this residue, return false.
	}

	if ( prohibited_res_identities_iplus1_.size()>0 ) { //If there are prohibited res_identities
		for ( core::Size i=1, imax=prohibited_res_identities_iplus1_.size(); i<=imax; ++i ) { //Loop through the prohibited res_identities
			if ( rsd.name3() == prohibited_res_identities_iplus1_[i] ) return false; //If the residue matches any prohibited residue identity, return false.
		}
	}

	//if(TR.visible()) TR << "Residue " << rsd.name3() << rsd.seqpos() << " matches!" << std::endl; //DELETE ME!

	return true; //Nothing has failed at this point, so it's a match!
} //criteria_match_iplus1

/// @brief Select a random bin from the bins allowed for residue i, biased by the relative counts for that bin.
///
core::Size BinTransitionData::random_bin_i() const
{
	if ( !matrix_initialized() || !matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::random_bin_i(): Transition probability matrix has not been loaded!\n" );

	return numeric::random::pick_random_index_from_cdf( binsums_i_cdf_, numeric::random::rg() );
} //random_bin_i

/// @brief Select a random bin from the bins allowed for residue i+1, biased by the relative counts for that bin.
///
core::Size BinTransitionData::random_bin_iplus1() const
{
	if ( !matrix_initialized() || !matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::random_bin_iplus1(): Transition probability matrix has not been loaded!\n" );

	return numeric::random::pick_random_index_from_cdf( binsums_iplus1_cdf_, numeric::random::rg() );
} //random_bin_iplus1

/// @brief Select a random bin from the bins allowed for residue i+1, biased by the relative counts for that bin given that residue i is in a particular bin.
///
core::Size BinTransitionData::random_bin_iplus1_given_i( core::Size const bin_i ) const
{
	if ( !matrix_initialized() || !matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::random_bin_iplus1_given_i(): The transition probability matrix has not been loaded!\n" );
	if ( bin_i < 1 || bin_i > bin_iplus1_cdf_given_i_.size() || bin_i > n_bins_i() ) utility_exit_with_message( "In BinTransitionData::random_bin_iplus1_given_i(): The bin index for the ith residue is out of range.\n" );
	return numeric::random::pick_random_index_from_cdf( bin_iplus1_cdf_given_i_[bin_i], numeric::random::rg() );
} //random_bin_iplus1_given_i


/// @brief Given a vector of mainchain torsions for a particular residue, figure out which bin the torsions lie in.
///
core::Size BinTransitionData::which_bin_i( utility::vector1 < core::Real > const &mainchain_torsions ) const
{
	if ( !matrix_initialized() || !matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::which_bin_i(): The transition probability matrix has not been loaded!\n" );
	if ( mainchain_torsions.size() != n_mainchain_torsions_i() ) utility_exit_with_message( "In BinTransitionData::which_bin_i(): The number of mainchain torsions in the test vector does not match the residue type for this transition probability matrix!\n" );
	core::Size bin_index(0);
	bool found(true);
	for ( core::Size i=1, imax=n_bins_i(); i<=imax; ++i ) { //Loop through all bins.
		found=true;
		for ( core::Size j=1, jmax=n_mainchain_torsions_i(); j<=jmax; ++j ) { //Loop through all mainchain torsions.
			if ( !in_bin( binranges_i_[i][j].first, binranges_i_[i][j].second, set_in_range(mainchain_torsions[j])) ) {
				found=false; //If it's not in one of the bin ranges corresponding to one of the mainchain torsions, then it's not in the bin.
				break; //So there's no point checking the other mainchain torsions.
			}
		} //Looping through mainchain torsions
		if ( found ) {
			bin_index=i;
			break; //Take the first bin found.  (There should be only one).
		}
	}

	if ( bin_index==0 || !found ) utility_exit_with_message( "In BinTransitionData::which_bin(): The bin corresponding to the input vector could not be found!\n" );

	return bin_index;
} //which_bin_i

/// @brief Given a vector of mainchain torsions for a particular residue, figure out which bin the torsions lie in.
///
core::Size BinTransitionData::which_bin_iplus1( utility::vector1 < core::Real > const &mainchain_torsions ) const
{
	if ( !matrix_initialized() || !matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::which_bin_iplus1(): The transition probability matrix has not been loaded!\n" );
	if ( mainchain_torsions.size() != n_mainchain_torsions_iplus1() ) utility_exit_with_message( "In BinTransitionData::which_bin_iplus1(): The number of mainchain torsions in the test vector does not match the residue type for this transition probability matrix!\n" );
	core::Size bin_index(0);
	bool found(true);
	for ( core::Size i=1, imax=n_bins_iplus1(); i<=imax; ++i ) { //Loop through all bins.
		found=true;
		for ( core::Size j=1, jmax=n_mainchain_torsions_iplus1(); j<=jmax; ++j ) { //Loop through all mainchain torsions.
			if ( !in_bin( binranges_iplus1_[i][j].first, binranges_iplus1_[i][j].second, set_in_range(mainchain_torsions[j])) ) {
				found=false; //If it's not in one of the bin ranges corresponding to one of the mainchain torsions, then it's not in the bin.
				break; //So there's no point checking the other mainchain torsions.
			}
		} //Looping through mainchain torsions
		if ( found ) {
			bin_index=i;
			break; //Take the first bin found.  (There should be only one).
		}
	}

	if ( bin_index==0 || !found ) utility_exit_with_message( "In BinTransitionData::which_bin(): The bin corresponding to the input vector could not be found!\n" );

	return bin_index;
} //which_bin_iplus1

/// @brief Is a given residue within the bounds of a given bin?
/// @details Uses bin definitions for the ith residue.
bool BinTransitionData::in_bin_i( core::Size const bin_index, core::conformation::Residue const &rsd  ) const
{
	return ( bin_index == which_bin_i( rsd.mainchain_torsions() ) );
} //in_bin_i

/// @brief Is a given residue within the bounds of a given bin?
/// @details Uses bin definitions for the i+1st residue.
bool BinTransitionData::in_bin_iplus1( core::Size const bin_index, core::conformation::Residue const &rsd  ) const
{
	return ( bin_index == which_bin_iplus1( rsd.mainchain_torsions() ) );
} //in_bin_iplus1

/// @brief Do final post-load calculations (e.g. precomputing the sum of the transition probability matrix entries,
/// the sum of the columns, the sum of the rows, etc.)
void BinTransitionData::finalize() {
	if ( !matrix_initialized() ) utility_exit_with_message( "In BinTransitionData::finalize(): Transition probability matrix not yet initialized!\n" );
	if ( matrix_finalized() ) utility_exit_with_message( "In BinTransitionData::finalize(): Transition probability matrix already finalized!\n" );

	core::Size const nrows( n_bins_i() );
	core::Size const ncols( n_bins_iplus1() );
	if ( nrows<1 || ncols<1 ) utility_exit_with_message( "In BinTransitionData::finalize(): The number of bins cannot be zero!\n" );

	//Initialize the binsums vectors
	binsums_i_.clear();
	binsums_i_.resize( nrows, 0.0 );
	binsums_iplus1_.clear();
	binsums_iplus1_.resize( ncols, 0.0 );
	total_counts_=0.0;

	//Sum the columns and rows:
	for ( core::Size i=1; i<=nrows; ++i ) {
		for ( core::Size j=1; j<=ncols; ++j ) {
			binsums_i_[i]+=probability_matrix(i,j);
			binsums_iplus1_[j]+=probability_matrix(i,j);
		}
		total_counts_+=binsums_i_[i]; //The row has now been summed, so the total counts can be derived from the sum of the row sums.
	}

	//Compute the cumulative probability distributions:
	binsums_i_cdf_.clear();
	binsums_i_cdf_.resize( nrows, 0.0 );
	binsums_iplus1_cdf_.clear();
	binsums_iplus1_cdf_.resize( nrows, 0.0 );
	for ( core::Size i=1; i<=nrows; ++i ) {
		if ( i>1 ) binsums_i_cdf_[i]=binsums_i_cdf_[i-1];
		binsums_i_cdf_[i]+=binsums_i_[i]/total_counts_;
	}
	for ( core::Size j=1; j<=ncols; ++j ) {
		if ( j>1 ) binsums_iplus1_cdf_[j]=binsums_iplus1_cdf_[j-1];
		binsums_iplus1_cdf_[j]+=binsums_iplus1_[j]/total_counts_;
	}

	//Note: Rosetta's cumulative distribution (cdf) functions are strange.
	//They assume a cdf where P(cell i) = cdf[i+1]-cdf[i], rather than
	//P(cell i) = cdf[i]-cdf[i-1].  This means I need to adjust, here:
	for ( core::Size i=nrows; i>1; --i ) {
		binsums_i_cdf_[i]=binsums_i_cdf_[i-1];
	}
	if ( nrows>0 ) binsums_i_cdf_[1] = 0;
	for ( core::Size j=ncols; j>1; --j ) {
		binsums_iplus1_cdf_[j]=binsums_iplus1_cdf_[j-1];
	}
	if ( ncols>0 ) binsums_iplus1_cdf_[1] = 0;

	//Set up the sub-bins and their cumulative probability distributions:
	set_up_subbins();

	//Set up cumulative probability distribution functions for probability of bin at i+1 given bin at i, and vice versa.
	set_up_bin_cdfs();

	matrix_finalized_=true;

	return;
} //finalize

/// @brief Does a bin with the specified name exist?
///
bool
BinTransitionData::bin_exists(
	std::string const &name
) const {
	return bin_exists_i(name) || bin_exists_iplus1(name);
}

/// @brief Does a bin with the specified name exist, defined for the ith residue?
///
bool
BinTransitionData::bin_exists_i(
	std::string const &name
) const {
	return is_in_list( name, binnames_i_);
}

/// @brief Does a bin with the specified name exist, defined for the i+1st residue?
///
bool
BinTransitionData::bin_exists_iplus1(
	std::string const &name
) const {
	return is_in_list( name, binnames_iplus1_);
}


/// @brief Writes a report summarizing the data stored in this BinTransitionData object.
/// @details If verbose is true, the full set of sub-bins is written out, too.
std::string BinTransitionData::summarize_data( bool const verbose ) const
{
	std::ostringstream outstream;

	outstream << "\t" << "Mainchain torsions for ith residue: " << n_mainchain_torsions_i() << std::endl;
	outstream << "\t" << "Mainchain torsions for i+1st residue: " << n_mainchain_torsions_iplus1() << std::endl;
	outstream << "\t" << "Bins for ith residue: " << n_bins_i() << std::endl;
	outstream << "\t" << "Bins for i+1st residue: " << n_bins_iplus1() << std::endl << std::endl;

	// Print the bin ranges:
	outstream << "\tBin_i";
	for ( core::Size i=1, imax=n_mainchain_torsions_i(); i<=imax; ++i ) outstream << "\tStart_tors" << i << "\tEnd_tors" << i ;
	outstream << std::endl;
	for ( core::Size i=1, imax=n_bins_i(); i<=imax; ++i ) {
		outstream << "\t" << binnames_i_[i];
		for ( core::Size j=1, jmax=n_mainchain_torsions_i(); j<=jmax; ++j ) {
			outstream << "\t" << bin_boundaries_i(i,j).first << "\t" << bin_boundaries_i(i,j).second;
		}
		outstream << std::endl;
	}
	outstream << std::endl;
	outstream << "\tBin_iplus1";
	for ( core::Size i=1, imax=n_mainchain_torsions_iplus1(); i<=imax; ++i ) outstream << "\tStart_tors" << i << "\tEnd_tors" << i ;
	outstream << std::endl;
	for ( core::Size i=1, imax=n_bins_iplus1(); i<=imax; ++i ) {
		outstream << "\t" << binnames_iplus1_[i];
		for ( core::Size j=1, jmax=n_mainchain_torsions_iplus1(); j<=jmax; ++j ) {
			outstream << "\t" << bin_boundaries_iplus1(i,j).first << "\t" << bin_boundaries_iplus1(i,j).second;
		}
		outstream << std::endl;
	}
	outstream << std::endl;

	// Allowed and prohibited residues:
	outstream << "\t" << "Allowed residues at ith residue: ";
	if ( res_identities_i_.size() > 0 ) {
		for ( core::Size i=1,imax=res_identities_i_.size(); i<=imax; ++i ) {
			outstream << res_identities_i_[i] << " ";
		}
	} else outstream << "(All residues)";
	outstream << std::endl;

	outstream << "\t" << "Prohibited residues at ith residue: ";
	if ( prohibited_res_identities_i_.size() > 0 ) {
		for ( core::Size i=1,imax=prohibited_res_identities_i_.size(); i<=imax; ++i ) {
			outstream << prohibited_res_identities_i_[i] << " ";
		}
	} else outstream << "(No residues)";
	outstream << std::endl;

	outstream << "\t" << "Allowed residues at i+1st residue: ";
	if ( res_identities_iplus1_.size() > 0 ) {
		for ( core::Size i=1,imax=res_identities_iplus1_.size(); i<=imax; ++i ) {
			outstream << res_identities_iplus1_[i] << " ";
		}
	} else outstream << "(All residues)";
	outstream << std::endl;

	outstream << "\t" << "Prohibited residues at i+1st residue: ";
	if ( prohibited_res_identities_iplus1_.size() > 0 ) {
		for ( core::Size i=1,imax=prohibited_res_identities_iplus1_.size(); i<=imax; ++i ) {
			outstream << prohibited_res_identities_iplus1_[i] << " ";
		}
	} else outstream << "(No residues)";
	outstream << std::endl << std::endl;


	// Allowed and prohibited properties:
	outstream << "\t" << "Allowed properties at ith residue: ";
	if ( properties_i_.size() > 0 ) {
		for ( core::Size i=1,imax=properties_i_.size(); i<=imax; ++i ) {
			outstream << get_property_effect_name(properties_i_[i]) << " ";
		}
	} else outstream << "(All properties)";
	outstream << std::endl;

	outstream << "\t" << "Prohibited properties at ith residue: ";
	if ( prohibited_properties_i_.size() > 0 ) {
		for ( core::Size i=1,imax=prohibited_properties_i_.size(); i<=imax; ++i ) {
			outstream << get_property_effect_name(prohibited_properties_i_[i]) << " ";
		}
	} else outstream << "(No properties)";
	outstream << std::endl;

	outstream << "\t" << "Allowed properties at i+1st residue: ";
	if ( properties_iplus1_.size() > 0 ) {
		for ( core::Size i=1,imax=properties_iplus1_.size(); i<=imax; ++i ) {
			outstream << get_property_effect_name(properties_iplus1_[i]) << " ";
		}
	} else outstream << "(All properties)";
	outstream << std::endl;

	outstream << "\t" << "Prohibited properties at i+1st residue: ";
	if ( prohibited_properties_iplus1_.size() > 0 ) {
		for ( core::Size i=1,imax=prohibited_properties_iplus1_.size(); i<=imax; ++i ) {
			outstream << get_property_effect_name(prohibited_properties_iplus1_[i]) << " ";
		}
	} else outstream << "(No properties)";
	outstream << std::endl << std::endl;


	//Print the matrix:
	outstream << "\t\t";
	for ( core::Size i=1, imax=n_bins_iplus1(); i<=imax; ++i ) outstream << binnames_iplus1_[i] << "\t";
	outstream << std::endl;

	for ( core::Size i=1, imax=n_bins_i(); i<=imax; ++i ) {
		outstream << "\t" << binnames_i_[i] << "\t";
		for ( core::Size j=1, jmax=n_bins_iplus1(); j<=jmax; ++j ) outstream << probability_matrix(i,j) << "\t";
		outstream << std::endl;
	}
	outstream << std::endl;

	//Print the total count:
	outstream << "\tTotal count (sum of entries in the probability matrix):\t" << total_counts_ << std::endl << std::endl;

	//Print the sums for each bin and the cumulative probability distribution
	outstream << "\tBin_i\tSum\tCumulative_dist" << std::endl;
	for ( core::Size i=1, imax=n_bins_i(); i<=imax; ++i ) {
		outstream << "\t" << binnames_i_[i] << "\t" << binsums_i_[i] << "\t" << binsums_i_cdf_[i] << std::endl;
	}
	outstream << std::endl;
	outstream << "\tBin_i+1\tSum\tCumulative_dist" << std::endl;
	for ( core::Size i=1, imax=n_bins_iplus1(); i<=imax; ++i ) {
		outstream << "\t" << binnames_iplus1_[i] << "\t" << binsums_iplus1_[i] << "\t" << binsums_iplus1_cdf_[i] << std::endl;
	}
	outstream << std::endl;

	//Print the sub-bin types:
	outstream << "\tSub-bin type for residue i:\t" << get_subbin_type_name(subbin_type_i_) << "." << std::endl;
	if ( verbose ) { //Only print out the full set of sub-bins in verbose mode.
		for ( core::Size i=1,imax=n_bins_i(); i<=imax; ++i ) {
			outstream << "\tSub-bins for " << binnames_i_[i] << " bin:" << std::endl;
			outstream << "\tSub-bins\tcdf";
			for ( core::Size j=1, jmax=n_mainchain_torsions_i(); j<=jmax; ++j ) outstream << "\tstart" << j << "\tend" << j;
			outstream << std::endl;
			for ( core::Size j=1, jmax=subbin_cdf_i_[i].size(); j<=jmax; ++j ) {
				outstream << "\t" << j << "\t" << subbin_cdf_i_[i][j];
				for ( core::Size k=1, kmax=n_mainchain_torsions_i(); k<=kmax; ++k ) {
					outstream << "\t" << subbin_ranges_i_[i][j][k].first << "\t" << subbin_ranges_i_[i][j][k].second;
				}
				outstream << std::endl;
			}
		}
		outstream << std::endl;
	}
	outstream << "\tSub-bin type for residue i+1:\t" << get_subbin_type_name(subbin_type_iplus1_) << "." << std::endl;
	if ( verbose ) { //Only print out the full set of sub-bins in verbose mode.
		for ( core::Size i=1,imax=n_bins_iplus1(); i<=imax; ++i ) {
			outstream << "\tSub-bins for " << binnames_iplus1_[i] << " bin:" << std::endl;
			outstream << "\tSub-bins\tcdf";
			for ( core::Size j=1, jmax=n_mainchain_torsions_iplus1(); j<=jmax; ++j ) outstream << "\tstart" << j << "\tend" << j;
			outstream << std::endl;
			for ( core::Size j=1, jmax=subbin_cdf_iplus1_[i].size(); j<=jmax; ++j ) {
				outstream << "\t" << j << "\t" << subbin_cdf_iplus1_[i][j];
				for ( core::Size k=1, kmax=n_mainchain_torsions_iplus1(); k<=kmax; ++k ) {
					outstream << "\t" << subbin_ranges_iplus1_[i][j][k].first << "\t" << subbin_ranges_iplus1_[i][j][k].second;
				}
				outstream << std::endl;
			}
		}
		outstream << std::endl;
	}
	outstream << std::endl;

	//Print out the cumulative distribution functions for the i+1st residue, given the ith:
	outstream << "\tCumulative distribution functions for i+1st residue's bins, given bin for ith:" << std::endl;
	outstream << "\t\tBin_i+1" << std::endl;
	outstream << "\tBin_i";
	for ( core::Size i=1, imax=n_bins_iplus1(); i<=imax; ++i ) {
		outstream << "\t" << binnames_iplus1_[i];
	}
	outstream << std::endl;
	for ( core::Size i=1, imax=n_bins_i(); i<=imax; ++i ) {
		outstream << "\t" << binnames_i_[i];
		for ( core::Size j=1, jmax=n_bins_iplus1(); j<=jmax; ++j ) {
			outstream << "\t" << bin_iplus1_cdf_given_i_[i][j];
		}
		outstream << std::endl;
	}
	outstream << std::endl;

	//Print out the cumulative distribution functions for the ith residue, given the i+st:
	outstream << "\tCumulative distribution functions for ith residue's bins, given bin for i+1st residue:" << std::endl;
	outstream << "\t\tBin_i+1" << std::endl;
	outstream << "\tBin_i";
	for ( core::Size i=1, imax=n_bins_iplus1(); i<=imax; ++i ) {
		outstream << "\t" << binnames_iplus1_[i];
	}
	outstream << std::endl;
	for ( core::Size i=1, imax=n_bins_i(); i<=imax; ++i ) {
		outstream << "\t" << binnames_i_[i];
		for ( core::Size j=1, jmax=n_bins_iplus1(); j<=jmax; ++j ) {
			outstream << "\t" << bin_i_cdf_given_iplus1_[j][i];
		}
		outstream << std::endl;
	}
	outstream << std::endl;

	outstream << std::endl;

	return outstream.str();
} //summarize_data


/// @brief Set the number of bins for the ith and i+1st residues.
/// @details This also initializes the probability matrix (to all zeros), the binnames vectors
/// (to vectors of empty strings), and the binranges vectors (to all zeros).
void BinTransitionData::set_n_bins( core::Size const n_bins_i, core::Size const n_bins_iplus1 ) {
	if ( matrix_initialized() && TR.Warning.visible() ) TR.Warning << "Warning: re-initializing probability matrix in BinTransitionData object!" << std::endl;

	runtime_assert_string_msg( n_mainchain_torsions_i_ > 0 && n_mainchain_torsions_iplus1_ > 0,
		"In core::scoring::bin_transitions::BinTransitionData::set_n_bins(): the number of mainchain torsions for residues i and i+1 must be set before calling this function." );

	runtime_assert_string_msg( n_bins_i > 0, "In core::scoring::bin_transitions::BinTransitionData::set_n_bins(): The number of bins for the ith residue must be greater than zero!"  );
	runtime_assert_string_msg( n_bins_iplus1 > 0, "In core::scoring::bin_transitions::BinTransitionData::set_n_bins(): The number of bins for the i plus 1st residue must be greater than zero!"  );

	//Set number of bins for residues i and i+1:
	n_bins_i_ = n_bins_i;
	n_bins_iplus1_ = n_bins_iplus1;

	//Initialize the binnames vectors
	binnames_i_.clear();
	binnames_i_.resize(n_bins_i_, "");
	binnames_iplus1_.clear();
	binnames_iplus1_.resize(n_bins_iplus1_, "");

	//Initialize the binranges vectors
	binranges_i_.clear();
	for ( core::Size i=1; i<=n_bins_i_; ++i ) {
		utility::vector1 < std::pair < core::Real, core::Real > > tempvect;
		for ( core::Size j=1; j<=n_mainchain_torsions_i_; ++j ) tempvect.push_back( std::pair< core::Real, core::Real>(0,0) );
		binranges_i_.push_back(tempvect);
	}
	binranges_iplus1_.clear();
	for ( core::Size i=1; i<=n_bins_iplus1_; ++i ) {
		utility::vector1 < std::pair < core::Real, core::Real > > tempvect;
		for ( core::Size j=1; j<=n_mainchain_torsions_iplus1_; ++j ) tempvect.push_back( std::pair< core::Real, core::Real>(0,0) );
		binranges_iplus1_.push_back(tempvect);
	}

	//Initialize the probability matrix:
	probability_matrix_.clear();
	for ( core::Size i=1; i<=n_bins_i_; ++i ) {
		utility::vector1 < core::Real > tempvect;
		tempvect.resize(n_bins_iplus1_, 0.0);
		probability_matrix_.push_back(tempvect);
	}
	matrix_initialized_ = true;
	matrix_finalized_ = false;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Initialized probability matrix to " << n_bins_i_ << " by " << n_bins_iplus1 << " size." << std::endl;
		TR.Debug.flush();
	}

	return;
} //set_n_bins

} //namespace bin_transitions
}//namespace scoring
}//namespace core
