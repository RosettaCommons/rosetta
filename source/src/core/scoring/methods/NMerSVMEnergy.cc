// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/NMerSVMEnergy.cc
/// @brief  SVM sequence profile energy method implementation
/// @author Indigo King (indigo.c.king@gmail.com)
/// @author Vikram K. Mulligan (vmulligan@flatironinstititue.org) -- Implemented lazy, one-time, threadsafe loading of all data read from files, and
/// replaced reads of global options system with reads from local EnergyMethodsOptions class.

// Unit headers
#include <core/scoring/methods/NMerSVMEnergy.hh>
#include <core/scoring/methods/NMerSVMEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/methods/NMerPSSMEnergy.hh>
#include <core/scoring/ScoringManager.hh>
#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>

// C++ Headers
#include <string>
#include <vector>

// Utility Headers
#include <utility/string_util.hh>

#include <utility/libsvm/Svm_rosetta.hh>
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

#include <utility/exit.hh>

static basic::Tracer TR( "core.scoring.methods.NMerSVMEnergy" );

namespace core {
namespace scoring {
namespace methods {

using namespace utility::libsvm;
using utility::vector1;

/// @details This must return a fresh instance of the NMerSVMEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
NMerSVMEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &options
) const {
	return utility::pointer::make_shared< NMerSVMEnergy >( options );
}

ScoreTypes
NMerSVMEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( nmer_svm );
	return sts;
}

void
NMerSVMEnergy::nmer_length( Size const nmer_length ){
	nmer_length_ = nmer_length;
	//nmer residue energy is attributed to position 1
	nmer_cterm_ = nmer_length_ - 1 ;
}

Size
NMerSVMEnergy::n_svms() const
{
	return all_nmer_svms_.size();
}

Size
NMerSVMEnergy::nmer_length() const
{
	return nmer_length_;
}

void
NMerSVMEnergy::term_length( Size const term_length ){
	term_length_ = term_length;
}

Size
NMerSVMEnergy::term_length() const
{
	return term_length_;
}

void
NMerSVMEnergy::use_pssm_features( bool const use_pssm_features ){
	use_pssm_features_ = use_pssm_features;
}

void
NMerSVMEnergy::avg_rank_as_energy( bool const avg_rank_as_energy ){
	avg_rank_as_energy_ = avg_rank_as_energy;
}

void
NMerSVMEnergy::gate_svm_scores( bool const gate_svm_scores ){
	gate_svm_scores_ = gate_svm_scores;
}

void
NMerSVMEnergy::nmer_svm_scorecut( Real const nmer_svm_scorecut ){
	nmer_svm_scorecut_ = nmer_svm_scorecut;
}

void
NMerSVMEnergy::initialize_from_options( EnergyMethodOptions const & options )
{
	NMerSVMEnergy::nmer_length( options.nmer_ref_seq_length() );
	NMerSVMEnergy::gate_svm_scores( options.nmer_svm_scorecut_defined() );
	NMerSVMEnergy::term_length( options.nmer_svm_term_length() );
	NMerSVMEnergy::use_pssm_features( options.nmer_svm_pssm_feat() );
	NMerSVMEnergy::nmer_svm_scorecut( options.nmer_svm_scorecut() );
	NMerSVMEnergy::avg_rank_as_energy( options.nmer_svm_avg_rank_as_energy() );
	//load user-defined encoding?
	if ( options.nmer_svm_aa_matrix_defined() ) {
		NMerSVMEnergy::read_aa_encoding_matrix( options.nmer_svm_aa_matrix() );
	} else {
		NMerSVMEnergy::read_aa_encoding_matrix(
			//or default to database BLOSUM62
			basic::database::full_name( "sequence/substitution_matrix/BLOSUM62.prob.rescale" ) );
	}
}

NMerSVMEnergy::NMerSVMEnergy( EnergyMethodOptions const & options ) :
	parent( methods::EnergyMethodCreatorOP( utility::pointer::make_shared< NMerSVMEnergyCreator >() ) )
{
	NMerSVMEnergy::initialize_from_options( options );
	read_nmer_svms_from_options( options );
}

//init from scratch w/ a vector of libsvm model and rank filenames
NMerSVMEnergy::NMerSVMEnergy( EnergyMethodOptions const & options, utility::vector1< std::string > const & svm_fnames, utility::vector1< std::string > const & svm_rank_fnames ):
	parent( methods::EnergyMethodCreatorOP( utility::pointer::make_shared< NMerSVMEnergyCreator >() ) )
{
	NMerSVMEnergy::initialize_from_options( options );

	all_nmer_svms_.clear();
	all_nmer_svms_ranks_.clear();
	for ( Size isvm = 1; isvm <= svm_fnames.size(); ++isvm ) {
		NMerSVMEnergy::read_nmer_svm( svm_fnames[ isvm ] );
		NMerSVMEnergy::read_nmer_svm_rank( svm_rank_fnames[ isvm ] );
	}
}

// full ctor init with vector of svm filenames with default blosum encoding matrix
NMerSVMEnergy::NMerSVMEnergy(
	core::Size const nmer_length,
	bool const gate_svm_scores,
	core::Size const term_length,
	bool const use_pssm_features,
	bool const avg_rank_as_energy,
	core::Real const nmer_svm_scorecut,
	utility::vector1< std::string > const & svm_fname_vec,
	utility::vector1< std::string > const & svm_rank_fname_vec,
	utility::vector1< std::string > const & pssm_fname_vec
) :
	parent( utility::pointer::make_shared< NMerSVMEnergyCreator >() )
{
	NMerSVMEnergy::nmer_length( nmer_length  );
	NMerSVMEnergy::gate_svm_scores( gate_svm_scores );
	NMerSVMEnergy::term_length( term_length );
	NMerSVMEnergy::use_pssm_features( use_pssm_features );
	NMerSVMEnergy::avg_rank_as_energy( avg_rank_as_energy );
	NMerSVMEnergy::nmer_svm_scorecut( nmer_svm_scorecut );
	NMerSVMEnergy::read_nmer_svm_fname_vector( svm_fname_vec );
	NMerSVMEnergy::read_nmer_svm_rank_fname_vector( svm_rank_fname_vec );
	// set pssm options and pssm file paths
	nmer_pssm_.nmer_length( nmer_length );
	nmer_pssm_.gate_pssm_scores( false );
	nmer_pssm_.nmer_pssm_scorecut( 0.0 );
	nmer_pssm_.read_nmer_pssm_fname_vector( pssm_fname_vec );
	NMerSVMEnergy::read_aa_encoding_matrix(
		//default to database BLOSUM62
		basic::database::full_name( "sequence/substitution_matrix/BLOSUM62.prob.rescale" ) );
}

// full ctor init with vector of svm filenames with custom encoding matrix
NMerSVMEnergy::NMerSVMEnergy(
	Size const nmer_length,
	bool const gate_svm_scores,
	Size const term_length,
	bool const use_pssm_features,
	bool const avg_rank_as_energy,
	core::Real const nmer_svm_scorecut,
	utility::vector1< std::string > const & svm_fname_vec,
	utility::vector1< std::string > const & svm_rank_fname_vec,
	utility::vector1< std::string > const & pssm_fname_vec,
	std::string const & aa_matrix
) :
	parent( utility::pointer::make_shared< NMerSVMEnergyCreator >() )
{
	NMerSVMEnergy::nmer_length( nmer_length  );
	NMerSVMEnergy::gate_svm_scores( gate_svm_scores );
	NMerSVMEnergy::term_length( term_length );
	NMerSVMEnergy::use_pssm_features( use_pssm_features );
	NMerSVMEnergy::avg_rank_as_energy( avg_rank_as_energy );
	NMerSVMEnergy::nmer_svm_scorecut( nmer_svm_scorecut );
	NMerSVMEnergy::read_nmer_svm_fname_vector( svm_fname_vec );
	NMerSVMEnergy::read_nmer_svm_rank_fname_vector( svm_rank_fname_vec );
	// set pssm options and pssm file paths
	nmer_pssm_.nmer_length( nmer_length );
	nmer_pssm_.gate_pssm_scores( false );
	nmer_pssm_.nmer_pssm_scorecut( 0.0 );
	nmer_pssm_.read_nmer_pssm_fname_vector( pssm_fname_vec );
	NMerSVMEnergy::read_aa_encoding_matrix( aa_matrix );
}

NMerSVMEnergy::~NMerSVMEnergy() = default;

void
NMerSVMEnergy::read_nmer_svms_from_options( EnergyMethodOptions const & options ) {

	TR << "checking for NMerSVMEnergy SVM list" << std::endl;

	//check for svm list file
	if ( options.nmer_svm_list_defined() ) {
		NMerSVMEnergy::read_nmer_svm_list( options.nmer_svm_list() );
	}
	//use single svm file
	if ( options.nmer_svm_defined() ) {
		NMerSVMEnergy::read_nmer_svm( options.nmer_svm() );
	}
	//check for svm ranks list file
	if ( options.nmer_svm_rank_list_defined() ) {
		NMerSVMEnergy::read_nmer_svm_rank_list( options.nmer_svm_rank_list() );
	}
	//use single svm ranks file
	if ( options.nmer_svm_rank_defined() ) {
		NMerSVMEnergy::read_nmer_svm_rank( options.nmer_svm_rank() );
	}
}

//load svms from a list file
void
NMerSVMEnergy::read_nmer_svm_list( std::string const & svm_list_fname ) {
	std::string const & svm_list_contents( core::scoring::ScoringManager::get_instance()->get_nmer_svm_list_file_contents( svm_list_fname ) );

	TR << "reading NMerSVMEnergy list from " << svm_list_fname << std::endl;
	TR << "SVM cache is cleared each time a list is loaded!" << std::endl;
	all_nmer_svms_.clear();

	//now loop over all names in list
	utility::vector1< std::string > const lines( utility::split_by_newlines(svm_list_contents) );
	for ( std::string const & line : lines ) {
		utility::vector1< std::string > const tokens( utility::split( line ) );
		//skip comments
		if ( tokens[ 1 ][ 0 ] == '#' ) continue;
		NMerSVMEnergy::read_nmer_svm( line );
	}
}

//load svms from a vector of filenames
void
NMerSVMEnergy::read_nmer_svm_fname_vector( utility::vector1< std::string > const & svm_fname_vec ) {
	//now loop over all names in vector
	for ( Size i = 1; i<= svm_fname_vec.size(); ++i ) {
		NMerSVMEnergy::read_nmer_svm( svm_fname_vec[ i ] );
	}
}

//load svm_ranks from a vector of filenames
void
NMerSVMEnergy::read_nmer_svm_rank_fname_vector( utility::vector1< std::string > const & svm_rank_fname_vec ) {
	//now loop over all names in vector
	for ( Size i = 1; i<= svm_rank_fname_vec.size(); ++i ) {
		NMerSVMEnergy::read_nmer_svm_rank( svm_rank_fname_vec[ i ] );
	}
}

void
NMerSVMEnergy::read_nmer_svm(
	std::string const & svm_fname
) {
	Svm_rosettaCOP nmer_svm( core::scoring::ScoringManager::get_instance()->get_nmer_svm( svm_fname ) ); //Svm_rosetta class has no clone operator, apparently.
	all_nmer_svms_.push_back( nmer_svm );
}

/// @brief load svm_ranks from a list file
void
NMerSVMEnergy::read_nmer_svm_rank_list( std::string const & svm_rank_list_fname ) {
	std::string const & svm_rank_list_contents( core::scoring::ScoringManager::get_instance()->get_nmer_svm_rank_list_file_contents( svm_rank_list_fname ) );

	TR << "reading NMerSVMEnergy rank list from " << svm_rank_list_fname << std::endl;

	//now loop over all names in list
	utility::vector1< std::string > const lines( utility::split_by_newlines(svm_rank_list_contents) );
	for ( std::string const & line : lines ) {
		utility::vector1< std::string > const tokens( utility::split( line ) );
		//skip comments
		if ( tokens[ 1 ][ 0 ] == '#' ) continue; //A bit of an odd way to do this, but okay. --VKM
		NMerSVMEnergy::read_nmer_svm_rank( line );
	}
}

/// @brief this is a sorted (ascending) list of precomputed energies of a bunch (100k) of random human peptides
/// @brief ...used to give a %ile rank score, useful for comparing diff MHC alelle predictions
void
NMerSVMEnergy::read_nmer_svm_rank( std::string const & svm_rank_fname ) {
	utility::vector1< core::Real > const & nmer_svm_rank( core::scoring::ScoringManager::get_instance()->get_nmer_svm_rank( svm_rank_fname ) );
	all_nmer_svms_ranks_.push_back( nmer_svm_rank );
}

//load the map that featurizes each aa in the sequence
//for smooth encoding, load BLOSUM62 or similar
//for exact encoding, load an identity matrix
void
NMerSVMEnergy::read_aa_encoding_matrix( std::string const & fname ){
	//clear the old data in case we're re-init'ing
	aa_encoder_.clear();
	aa_encoder_ = core::scoring::ScoringManager::get_instance()->get_nmer_svm_aa_matrix( fname );
}

//transform score matrix vals to probabilities, rescale to [-1,1]
//NO SILLY! just xform the matrix once and check it into rosetta_database!
EnergyMethodOP
NMerSVMEnergy::clone() const
{
	return utility::pointer::make_shared< NMerSVMEnergy >( *this );
}

//methods called by const methods must be const!
vector1< Svm_node_rosettaOP >
NMerSVMEnergy::get_svm_nodes( vector1< Real > const & feature_vals ) const
{
	vector1< Svm_node_rosettaOP > features;
	Size feature_idx( 1 );
	for ( Size ival = 1; ival <= feature_vals.size(); ++ival ) {
		features.push_back( utility::pointer::make_shared< Svm_node_rosetta >( feature_idx, feature_vals[ ival ] ) );
		++feature_idx;
	}
	return features;
}

//methods called by const methods must be const!
vector1< Real >
NMerSVMEnergy::encode_aa_string( std::string const & seq ) const
{
	vector1< Real > feature_vals;
	for ( Size iseq = 0; iseq <= seq.size() - 1; ++iseq ) {
		//iter through encoding vector for each char of seq
		char aa( seq[ iseq ] );
		//check for aa in encoder, else aa = 'X'
		//dont use map[] for access cuz [] can change map, and this is a const funxn!
		vector1< Real > aa_feature_vals;
		if ( aa_encoder_.find( aa ) != aa_encoder_.end() ) aa_feature_vals = aa_encoder_.find( aa )->second;
		else aa_feature_vals = aa_encoder_.find( 'X' )->second;
		for ( Size ival = 1; ival <= aa_feature_vals.size(); ++ival ) {
			feature_vals.push_back( aa_feature_vals[ ival ] );
		}
	}
	return feature_vals;
}

//return feature vector with features averaged over sequence positions
vector1< Real >
NMerSVMEnergy::encode_wtd_avg_aa_string( std::string const & seq, vector1< Real > const & wts ) const
{
	debug_assert( seq.length() == wts.size() );
	//get scale factors from relative wts, must sum to one!
	Real wtsum = 0;
	for ( Size iwt = 1; iwt <= wts.size(); ++iwt ) wtsum += wts[ iwt ];
	vector1< Real > wts_norm( wts );
	for ( Size iwt = 1; iwt <= wts_norm.size(); ++iwt ) wts_norm[ iwt ] /= wtsum;

	vector1< Real > feature_vals;
	for ( Size iseq = 0; iseq <= seq.size() - 1; ++iseq ) {
		std::string const aa( seq.substr( iseq, 1 ) );
		vector1< Real > aa_feature_vals( encode_aa_string( aa ) );
		for ( Size ival = 1; ival <= aa_feature_vals.size(); ++ival ) {
			Real val_norm( wts_norm[ iseq + 1 ] * aa_feature_vals[ ival ] );
			//if vector element DNE, add wtd feature value
			if ( ival > feature_vals.size() ) feature_vals.push_back( val_norm );
			//else sum it into the correct vector element
			else feature_vals[ ival ] += val_norm;
		}
	}
	return feature_vals;
}

//chain_seqpos is the P1 position of the epitope
vector1< Real >
NMerSVMEnergy::encode_nmer( std::string const & chain_sequence, Size const chain_seqpos, Size const isvm ) const
{
	std::string const nmer_seq( chain_sequence.substr( chain_seqpos - 1, nmer_length_ ) );
	vector1< Real > feature_vals( encode_aa_string( nmer_seq ) );
	if ( term_length_ > 0 ) {
		add_encoded_termini( chain_sequence, chain_seqpos, feature_vals );
	}
	if ( use_pssm_features_ ) {
		add_pssm_features( nmer_seq, isvm, feature_vals );
	}
	return feature_vals;
}

//chain_seqpos is the P1 position of the epitope
//OK to ask for chain seqpos + nmer that goes beyond end of string; just fill X.
void
NMerSVMEnergy::add_encoded_termini( std::string const & chain_sequence, Size const chain_seqpos, vector1< Real > & feature_vals ) const
{
	//get term nmers upstream and downstream of core, use dummy X for nonexistent res
	std::string nterm_seq, cterm_seq;
	for ( core::Size offset(1); offset <= term_length_; ++offset ) {

		//N terminus: if chain_seqpos is not greater than offset, we will underflow the string
		//also check that the index will be valid against the string's length.
		core::Size const nterm_index(( chain_seqpos - 1 ) - offset);
		if ( (chain_seqpos > offset) && (nterm_index < chain_sequence.length()) ) {
			nterm_seq = chain_sequence[ nterm_index ] + nterm_seq;
		} else {
			nterm_seq = "X" + nterm_seq;
		}

		//C terminus: (chain_seqpos - 1) means "chain_seqpos but indexed from 0 instead of 1 because strings".
		//(offset - 1) makes the offset line up: start the c terminus 0 characters after the end of the epitope, not 1 after (which would skip a position)
		//nmer_length_ skips us past the part of the sequence string that's the epitope we are looking at
		//This is not algebraically simplest but it is easiest to understand; the compiler is smart.
		core::Size const cterm_index(( chain_seqpos - 1 ) + nmer_length_ + ( offset - 1 ));

		//index==chain_sequence.length() is valid C++, but then chain_sequence[index] returns a null character
		if ( cterm_index < chain_sequence.length() ) {
			cterm_seq = cterm_seq + chain_sequence[ cterm_index ];
		} else {
			cterm_seq = cterm_seq + "X";
		}
	}

	//useful to debug that the above code actually worked...
	//TR << "chain sequence " << chain_sequence << " nterm " << nterm_seq << " cterm " << cterm_seq << std::endl;

	//construct wts
	vector1< Real > nterm_wts( term_length_, 1.0 ), cterm_wts( term_length_, 1.0 );
	//currently hardcoded to weight linearly down as we go away from central core sequence
	for ( Size iterm = 1; iterm <= term_length_; ++iterm ) {
		nterm_wts[ iterm ] = static_cast< Real >( iterm );
		cterm_wts[ term_length_ - iterm + 1 ] = static_cast< Real >( iterm );
	}
	//get avged term seqs
	vector1< Real > feat_vals_nterm( encode_wtd_avg_aa_string( nterm_seq, nterm_wts ) );
	vector1< Real > feat_vals_cterm( encode_wtd_avg_aa_string( cterm_seq, cterm_wts ) );
	//add to term feature vals to begin and end of feature vals
	for ( Size ival_term = 1; ival_term <= feat_vals_nterm.size(); ++ival_term ) {
		feature_vals.insert( feature_vals.begin() + ival_term - 1, feat_vals_nterm[ ival_term ] );
		feature_vals.push_back( feat_vals_cterm[ ival_term ] );
	}
}

// make sure you prescale these pssms! use the same ones for training and here
// lets prescale them by logoddsing against backgroun distribution
void
NMerSVMEnergy::add_pssm_features( std::string const & seq, Size const isvm, vector1< Real > & feature_vals ) const
{
	if ( isvm > nmer_pssm_.n_pssms() ) {
		utility_exit_with_message(
			"NMerSVMEnergy pssm features require one pssm per svm model. Check your svm and pssm lists for same length! NMerPSSMEnergy has " + std::to_string( nmer_pssm_.n_pssms() ) + " models defined. NMerSVMEnergy has " + std::to_string( n_svms() ) + " models defined\n" );
	}
	for ( Size iseq = 0; iseq <= seq.size() - 1; ++iseq ) {
		//iter through encoding vector for each char of seq
		core::chemical::AA aa( core::chemical::aa_from_oneletter_code( seq[ iseq ] ) );
		Real pssm_energy( nmer_pssm_.pssm_energy_at_frame_seqpos( iseq + 1, aa, isvm ) );
		feature_vals.push_back( pssm_energy );
	}
}

/// @brief the main compute function, rsd_energy_avg is the default return value for residue energy
/// @brief optionally return normalized rank avg value instead
void
NMerSVMEnergy::get_residue_energy_by_svm(
	pose::Pose const & pose,
	Size const & seqpos,
	Real & rsd_energy_avg,
	Real & rsd_rank_avg,
	vector1< Real > & rsd_svm_energies,
	vector1< Real > & rsd_svm_ranks
) const
{
	debug_assert( rsd_svm_energies.size() == n_svms() );
	//for now, just assign all of the p1=seqpos frame's nmer_svm energy to this residue
	//TODO: distribute frame's nmer_svm energy evenly across the nmer
	//TODO: avoid wasting calc time by storing nmer_value, nmer_val_out_of_date in pose cacheable data
	Size const p1_seqpos( seqpos );

	//get the nmer string
	//TODO: how deal w/ sequences shorter than nmer_length_?
	// this matters at both termini… maybe take max of all overlapping frames w/ missing res as 'X'?
	// go ahead and bail if we fall off the end of the chain
	rsd_energy_avg = 0.;
	rsd_rank_avg = 0.;
	if ( p1_seqpos + nmer_length_ - 1 <= pose.conformation().chain_end( pose.chain( p1_seqpos ) ) ) {
		//need p1 position in this chain's sequence, offset with index of first res
		Size const chain_p1_seqpos( p1_seqpos - pose.conformation().chain_begin( pose.chain( p1_seqpos ) ) + 1 );
		std::string const & chain_sequence( pose.chain_sequence( pose.chain( p1_seqpos ) ) );

		//encode nmer string --> feature(Svm_node) vector
		for ( Size isvm = 1; isvm <= n_svms(); ++isvm ) {
			vector1< Svm_node_rosettaOP > p1_seqpos_nmer_features(
				get_svm_nodes( encode_nmer( chain_sequence, chain_p1_seqpos, isvm ) ) );
			//get this svm model
			Svm_rosettaCOP nmer_svm_model( all_nmer_svms_[ isvm ] );
			//Svm_rosetta funxn to wrap svm_predict, see Svm_rosetta::predict_probability for template,
			Real rsd_svm_energy( nmer_svm_model->predict( p1_seqpos_nmer_features ) );
			//gate energy at svm_scorecut, thus ignoring low-scoring nmers
			if ( gate_svm_scores_ && rsd_svm_energy < nmer_svm_scorecut_ ) rsd_svm_energy = nmer_svm_scorecut_;
			//store this svms rsd energy
			rsd_svm_energies[ isvm ] = rsd_svm_energy;
			//normalize energy by number of svms used
			//otherwise avg scores would become huge if we use lots of svms instead of just 1
			rsd_energy_avg += ( rsd_svm_energy / n_svms() );

			// calc svm score rank for each svm
			vector1< core::Real > const & nmer_svm_rank( all_nmer_svms_ranks_[ isvm ] );
			Size rank( 1 );
			for ( Size irank = 1; irank <= nmer_svm_rank.size(); ++irank ) {
				if ( nmer_svm_rank[ irank ] > rsd_svm_energy || irank == nmer_svm_rank.size() ) {
					rank = irank;
					break;
				}
			}
			// calc fraction of rank, low rank is lower scoring subseq
			core::Real const rsd_svm_rank_fraction( core::Real( rank ) / core::Real( nmer_svm_rank.size() ) );
			rsd_svm_ranks[ isvm ] = rsd_svm_rank_fraction;
			rsd_rank_avg += ( rsd_svm_rank_fraction / n_svms() );
		}
	}
}

void
NMerSVMEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	if ( all_nmer_svms_.empty() ) return;
	Size const seqpos( rsd.seqpos() );

	Real rsd_energy_avg( 0. );
	Real rsd_rank_avg( 0. );
	vector1< Real > rsd_svm_energies( n_svms(), Real( 0. ) );
	vector1< Real > rsd_svm_ranks( n_svms(), Real( 0. ) );
	get_residue_energy_by_svm( pose, seqpos, rsd_energy_avg, rsd_rank_avg, rsd_svm_energies, rsd_svm_ranks );
	if ( avg_rank_as_energy_ ) emap[ nmer_svm ] += rsd_rank_avg;
	else emap[ nmer_svm ] += rsd_energy_avg;
}


Real
NMerSVMEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}

/// @brief NMerSVMEnergy is context independent; indicates that no
/// context graphs are required
void
NMerSVMEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
NMerSVMEnergy::version() const
{
	return 1; // Initial versioning
}
} // methods
} // scoring
} // core
