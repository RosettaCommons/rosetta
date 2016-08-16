// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/NMerPSSMEnergy.hh
/// @brief  PSSMerence energy method implementation
/// @author Chris King (dr.chris.king@gmail.com)

// Unit headers
#include <core/scoring/methods/NMerPSSMEnergy.hh>
#include <core/scoring/methods/NMerPSSMEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>

// C++ Headers
#include <string>
#include <vector>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.NMerPSSMEnergy" );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the NMerPSSMEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
NMerPSSMEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new NMerPSSMEnergy );
}

ScoreTypes
NMerPSSMEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( nmer_pssm );
	return sts;
}

void
NMerPSSMEnergy::nmer_length( Size const nmer_length ){
	nmer_length_ = nmer_length;
	//nmer residue energy is attributed to position 1
	nmer_cterm_ = nmer_length_ - 1 ;
}

void
NMerPSSMEnergy::gate_pssm_scores( bool const gate_pssm_scores ){
	gate_pssm_scores_ = gate_pssm_scores;
}

void
NMerPSSMEnergy::nmer_pssm_scorecut( Real const nmer_pssm_scorecut ){
	nmer_pssm_scorecut_ = nmer_pssm_scorecut;
}

void
NMerPSSMEnergy::initialize_from_options()
{
	using namespace basic::options;
	NMerPSSMEnergy::nmer_length( option[ OptionKeys::score::nmer_ref_seq_length ]() );
	NMerPSSMEnergy::gate_pssm_scores( option[ OptionKeys::score::nmer_pssm_scorecut ].user() );
	NMerPSSMEnergy::nmer_pssm_scorecut( option[ OptionKeys::score::nmer_pssm_scorecut ]() );
}

NMerPSSMEnergy::NMerPSSMEnergy() :
	parent( methods::EnergyMethodCreatorOP( new NMerPSSMEnergyCreator ) )
{
	NMerPSSMEnergy::initialize_from_options();
	read_nmer_pssms_from_options();
}

NMerPSSMEnergy::NMerPSSMEnergy( utility::vector1< std::map< chemical::AA, utility::vector1< core::Real > > > const & all_nmer_pssms_in ):
	parent( methods::EnergyMethodCreatorOP( new NMerPSSMEnergyCreator ) )
{
	//TODO: make this an argument of the function call
	NMerPSSMEnergy::initialize_from_options();

	all_nmer_pssms_.clear();
	for ( Size ipssm = 1; ipssm <= all_nmer_pssms_in.size(); ++ipssm ) {
		std::map< chemical::AA, utility::vector1< core::Real > > nmer_pssm;
		std::map< chemical::AA, utility::vector1< core::Real > > const nmer_pssm_in( all_nmer_pssms_in[ ipssm ] );
		//copy contents of input into new copy
		for ( std::map< chemical::AA, utility::vector1< Real > >::const_iterator it = nmer_pssm_in.begin(); it != nmer_pssm_in.end(); ++it ) {
			nmer_pssm.insert( *it );
		}
		//append new copy to our cleared instance
		all_nmer_pssms_.push_back( nmer_pssm );
	}
}

NMerPSSMEnergy::~NMerPSSMEnergy() {}

void NMerPSSMEnergy::read_nmer_pssms_from_options() {

	using namespace basic::options;

	TR << "checking for NMerPSSMEnergy PSSM list" << std::endl;

	//check for pssm list file
	if ( option[ OptionKeys::score::nmer_pssm_list ].user() ) {
		std::string const pssm_list_fname( option[ OptionKeys::score::nmer_pssm_list ] );
		NMerPSSMEnergy::read_nmer_pssm_list( pssm_list_fname );
	}
	//use single pssm file
	if ( option[ OptionKeys::score::nmer_pssm ].user() ) {
		std::string const pssm_fname( option[ OptionKeys::score::nmer_pssm ] );
		NMerPSSMEnergy::read_nmer_pssm( pssm_fname );
	}
}

//read energy table list
void NMerPSSMEnergy::read_nmer_pssm_list( std::string pssm_list_fname ) {
	if ( !utility::file::file_exists( pssm_list_fname ) ) {
		pssm_list_fname = basic::database::full_name( pssm_list_fname, false );
	}
	TR << "reading NMerPSSMEnergy list from " << pssm_list_fname << std::endl;
	utility::io::izstream in_stream( pssm_list_fname );
	if ( !in_stream.good() ) {
		utility_exit_with_message( "[ERROR] Error opening NMerPSSMEnergy list file" );
	}
	//now loop over all names in list
	std::string pssm_fname;
	while ( getline( in_stream, pssm_fname ) ) {
		utility::vector1< std::string > const tokens( utility::split( pssm_fname ) );
		//skip comments
		if ( tokens[ 1 ][ 0 ] == '#' ) continue;
		NMerPSSMEnergy::read_nmer_pssm( pssm_fname );
	}
}

//load PSSM with AA x seqpos scores
// PSSM format is 1 AA per line w/ nmer_length_ score vals
void NMerPSSMEnergy::read_nmer_pssm( std::string pssm_fname ) {

	if ( !utility::file::file_exists( pssm_fname ) ) {
		pssm_fname = basic::database::full_name( pssm_fname, false );
	}
	TR << "reading NMerPSSMEnergy scores from " << pssm_fname << std::endl;
	utility::io::izstream in_stream( pssm_fname );
	if ( !in_stream.good() ) {
		utility_exit_with_message( "[ERROR] Error opening NMerPSSMEnergy file" );
	}

	std::map< chemical::AA, utility::vector1< core::Real > > nmer_pssm;
	std::string line;
	while ( getline( in_stream, line) ) {
		utility::vector1< std::string > const tokens( utility::string_split_multi_delim( line, " \t" ) );
		//skip comments
		if ( tokens[ 1 ][ 0 ] == '#' ) continue;
		char const char_aa( tokens[ 1 ][ 0 ] );
		chemical::AA aa( chemical::aa_from_oneletter_code( char_aa ));
		if ( nmer_pssm.count( aa ) ) {
			utility_exit_with_message( "[ERROR] NMer PSSM energy database file "
				+ pssm_fname + " has double entry for aa " + char_aa );
		}
		if ( tokens.size() != nmer_length_ + 1 ) {
			utility_exit_with_message( "[ERROR] NMer PSSM database file "
				+ pssm_fname + " has wrong number entries at line " + line
				+ "\n\tfound: " + utility::to_string( tokens.size() ) + " expected: " + utility::to_string( Size( nmer_length_ + 1 ) ) + "\nNote: Whitespace delimited!" );
		}
		utility::vector1< Real > seqpos_scores( nmer_length_, 0.0 );
		for ( Size ival = 2; ival <= tokens.size(); ++ival ) {
			Real const score( atof( tokens[ ival ].c_str() ) );
			seqpos_scores[ ival - 1 ] = score;
		}
		nmer_pssm[ aa ] = seqpos_scores;
	}
	all_nmer_pssms_.push_back( nmer_pssm );
}


EnergyMethodOP
NMerPSSMEnergy::clone() const
{
	return EnergyMethodOP( new NMerPSSMEnergy( *this ) );
}

core::Size
NMerPSSMEnergy::n_pssms() const
{
	return all_nmer_pssms_.size();
}

core::Real
NMerPSSMEnergy::pssm_energy_at_frame_seqpos( Size const frame_seqpos, core::chemical::AA const aa, Size const idx_pssm ) const
{
	if ( !all_nmer_pssms_[ idx_pssm ].count( aa ) )  return 0.;
	return all_nmer_pssms_[ idx_pssm ].find( aa )->second[ frame_seqpos ];
}

//retrieves ref energy of NMer centered on seqpos
//we're changing this so energy is computed as sum of all frames that overlap w this seqpos
//that way, res energy is actually reflective
//…unless we recalc the whole pssm for each overlapping frame and reeval gate criterion
void
NMerPSSMEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace chemical;

	if ( all_nmer_pssms_.empty() ) return;
	Size const seqpos( rsd.seqpos() );
	//over each pssm
	for ( Size ipssm = 1; ipssm <= n_pssms(); ++ipssm ) {
		if ( all_nmer_pssms_[ ipssm ].empty() ) continue; //this really shouldn't happen, but just in case
		//calc nmer's score for this pssm
		Real rsd_energy( 0.0 );
		chemical::AA const rsd_aa( pose.residue( seqpos ).aa() );

		//loop effective p1 seqpos over all overlapping positions
		//need chain begin, end so dont run off end of sequence in multi-chain poses
		Size chain_begin( pose.conformation().chain_begin( pose.chain( seqpos ) ) );
		Size chain_end( pose.conformation().chain_end( pose.chain( seqpos ) ) );
		Size p1_seqpos_begin( seqpos - nmer_length_ + 1 < chain_begin ? chain_begin : seqpos - nmer_length_ + 1 );
		//will we run off end if we start p1 at seqpos?
		Size p1_seqpos_end( seqpos + nmer_length_ - 1 > chain_end ? chain_end - nmer_length_ + 1 : seqpos );
		//loop over each frame beginning
		for ( Size p1_seqpos = p1_seqpos_begin; p1_seqpos <= p1_seqpos_end; ++p1_seqpos ) {
			//get pssm index of seqpos in this p1 frame
			Size rsd_iseq_nmer( seqpos - p1_seqpos + 1 );
			//now go ahead and get pssm energy from all_nmer_pssms_[ ipssm ]
			Real rsd_energy_this_nmer( pssm_energy_at_frame_seqpos( rsd_iseq_nmer, rsd_aa, ipssm ) );

			//skip this part if not doing gating
			if ( gate_pssm_scores_ ) {
				Real energy( 0.0 );
				for ( Size iseq_nmer = 1; iseq_nmer <= nmer_length_; ++iseq_nmer ) {
					Size iseq_pose( iseq_nmer + p1_seqpos - 1 );
					//bail if we fall off end of chain
					if ( iseq_pose > chain_end ) break;
					chemical::AA const aa( pose.residue( iseq_pose ).aa() );
					Real this_rsd_energy( pssm_energy_at_frame_seqpos( iseq_nmer, aa, ipssm ) );
					energy += this_rsd_energy;
				}
				//gate energy at pssm_scorecut, thus ignoring low-scoring nmers
				//skip rsd_energy accumulation if total pssm score is low enough
				if ( energy < nmer_pssm_scorecut_ ) continue;
			}
			rsd_energy += rsd_energy_this_nmer;
		}
		//add sum of all frames' rsd energies into emap
		emap[ nmer_pssm ] += rsd_energy;
	}
	//normalize energy by number of pssms used
	//otherwise avg scores would become huge if we use lots of pssms instead of just 1
	emap[ nmer_pssm ] /= n_pssms();
	return;
}


Real
NMerPSSMEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}

/// @brief NMerPSSMEnergy is context independent; indicates that no
/// context graphs are required
void
NMerPSSMEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
NMerPSSMEnergy::version() const
{
	return 1; // Initial versioning
}
} // methods
} // scoring
} // core

