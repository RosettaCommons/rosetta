// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/NMerSVMEnergy.hh
/// @brief  SVMerence energy method implementation
/// @author Chris King (dr.chris.king@gmail.com)

// Unit headers
#include <core/scoring/methods/NMerSVMEnergy.hh>
#include <core/scoring/methods/NMerSVMEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/methods/NMerPSSMEnergy.hh>
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
#include <basic/database/open.hh>

#include <utility/libsvm/Svm_rosetta.hh>
#include <utility/vector1.hh>

static basic::Tracer TR( "core.scoring.methods.NMerSVMEnergy" );

namespace core {
namespace scoring {
namespace methods {
using namespace utility::libsvm;

/// @details This must return a fresh instance of the NMerSVMEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
NMerSVMEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new NMerSVMEnergy;
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
NMerSVMEnergy::gate_svm_scores( bool const gate_svm_scores ){
	gate_svm_scores_ = gate_svm_scores;
}

void
NMerSVMEnergy::nmer_svm_scorecut( Real const nmer_svm_scorecut ){
	nmer_svm_scorecut_ = nmer_svm_scorecut;
}

void
NMerSVMEnergy::initialize_from_options()
{
	using namespace basic::options;
	NMerSVMEnergy::nmer_length( option[ OptionKeys::score::nmer_ref_seq_length ]() ); 
	NMerSVMEnergy::gate_svm_scores( option[ OptionKeys::score::nmer_svm_scorecut ].user() ); 
	NMerSVMEnergy::term_length( option[ OptionKeys::score::nmer_svm_term_length ]() ); 
	NMerSVMEnergy::use_pssm_features( option[ OptionKeys::score::nmer_svm_pssm_feat ]() ); 
	NMerSVMEnergy::nmer_svm_scorecut( option[ OptionKeys::score::nmer_svm_scorecut ]() ); 
	//load user-defined encoding?
	if( option[ OptionKeys::score::nmer_svm_aa_matrix ].user() )
		NMerSVMEnergy::read_aa_encoding_matrix( option[ OptionKeys::score::nmer_svm_aa_matrix ]() ); 
	//or default to database BLOSUM62
	else NMerSVMEnergy::read_aa_encoding_matrix(
			basic::database::full_name( "sequence/substitution_matrix/BLOSUM62.prob.rescale" ) );
}

NMerSVMEnergy::NMerSVMEnergy() :
	parent( new NMerSVMEnergyCreator )
{
	NMerSVMEnergy::initialize_from_options();
	read_nmer_svms_from_options();
}

//init from scratch w/ a vecotr of libsvm model filenames
NMerSVMEnergy::NMerSVMEnergy( utility::vector1< std::string > const & svm_fnames ):
	parent( new NMerSVMEnergyCreator )
{
	NMerSVMEnergy::initialize_from_options();

	all_nmer_svms_.clear();
	for( Size isvm = 1; isvm <= svm_fnames.size(); ++isvm ){
		NMerSVMEnergy::read_nmer_svm( svm_fnames[ isvm ] );
	}
}

NMerSVMEnergy::~NMerSVMEnergy() {}

void
NMerSVMEnergy::read_nmer_svms_from_options() {

	using namespace basic::options;

	TR << "checking for NMerSVMEnergy SVM list" << std::endl;

	//check for svm list file
	if ( option[ OptionKeys::score::nmer_svm_list ].user() ) {
		std::string const svm_list_fname( option[ OptionKeys::score::nmer_svm_list ] );
		NMerSVMEnergy::read_nmer_svm_list( svm_list_fname );
	}
	//use single svm file
	if( option[ OptionKeys::score::nmer_svm ].user() ){
		std::string const svm_fname( option[ OptionKeys::score::nmer_svm ] );
		NMerSVMEnergy::read_nmer_svm( svm_fname );
	}
}

//load svms from a list file
void
NMerSVMEnergy::read_nmer_svm_list( std::string svm_list_fname ) {
  if ( !utility::file::file_exists( svm_list_fname ) ) {
    svm_list_fname = basic::database::full_name( svm_list_fname, false );
  }
	TR << "reading NMerSVMEnergy list from " << svm_list_fname << std::endl;
	utility::io::izstream in_stream( svm_list_fname );
	if (!in_stream.good()) {
		utility_exit_with_message( "Error opening NMerSVMEnergy list file" );
	}
	//now loop over all names in list
	std::string svm_fname;
	while( getline( in_stream, svm_fname ) ){
		utility::vector1< std::string > const tokens( utility::split( svm_fname ) );
		//skip comments
		if( tokens[ 1 ][ 0 ] == '#' ) continue;
		NMerSVMEnergy::read_nmer_svm( svm_fname );
	}
}

void
NMerSVMEnergy::read_nmer_svm( std::string svm_fname ) {

  if ( !utility::file::file_exists( svm_fname ) ) {
    svm_fname = basic::database::full_name( svm_fname, false );
  }
	TR << "reading NMerSVMEnergy scores from " << svm_fname << std::endl;
	const char* svm_fname_ch( svm_fname.c_str() );
	Svm_rosettaOP nmer_svm( new Svm_rosetta( svm_fname_ch ) );
	all_nmer_svms_.push_back( nmer_svm );
}

//load the map that featurizes each aa in the sequence
//for smooth encoding, load BLOSUM62 or similar
//for exact encoding, load an identity matrix
void
NMerSVMEnergy::read_aa_encoding_matrix( std::string const fname ){

	//clear the old data in case we're re-init'ing
	aa_encoder_.clear();

	//and load the data
  TR << "reading NMerSVM encoding matrix " << fname << std::endl;
  utility::io::izstream in_stream( fname );
  if (!in_stream.good()) {
    utility_exit_with_message( "[ERROR] Error opening NMerSVM encoding matrix file" );
  }
  std::string line;
  while( getline( in_stream, line) ) { 
    utility::vector1< std::string > const tokens( utility::string_split_multi_delim( line, " \t" ) );
    //skip comments
    if( tokens[ 1 ][ 0 ] == '#' ) continue;
    char const aa( tokens[ 1 ][ 0 ] );
    if( aa_encoder_.count( aa ) ) utility_exit_with_message( "[ERROR] NMer SVM encoding matrix file "
        + fname + " has double entry for aa " + aa );
    utility::vector1< Real > aa_vals;
    for( Size ival = 2; ival <= tokens.size(); ++ival ){
      Real const val( atof( tokens[ ival ].c_str() ) );
      aa_vals.push_back( val );
    }   
    aa_encoder_[ aa ] = aa_vals; 
  }
}

//transform score matrix vals to probabilities, rescale to [-1,1]
//NO SILLY! just xform the matrix once and check it into rosetta_database!
EnergyMethodOP
NMerSVMEnergy::clone() const
{
	return new NMerSVMEnergy( *this );
}

//methods called by const methods must be const!
vector1< Svm_node_rosettaOP >
NMerSVMEnergy::get_svm_nodes( vector1< Real > const feature_vals ) const
{
	vector1< Svm_node_rosettaOP > features;
	Size feature_idx( 1 );
	for( Size ival = 1; ival <= feature_vals.size(); ++ival ){
		features.push_back( Svm_node_rosettaOP( new Svm_node_rosetta( feature_idx, feature_vals[ ival ] ) ) );
		++feature_idx;
	}
	return features;
}

//methods called by const methods must be const!
vector1< Real >
NMerSVMEnergy::encode_aa_string( std::string const seq ) const
{
	vector1< Real > feature_vals;
	for( Size iseq = 0; iseq <= seq.size() - 1; ++iseq ){
		//iter through encoding vector for each char of seq
		char aa( seq[ iseq ] );
		//check for aa in encoder, else aa = 'X'
		//dont use map[] for access cuz [] can change map, and this is a const funxn!
		vector1< Real > aa_feature_vals;
		if( aa_encoder_.find( aa ) != aa_encoder_.end() ) aa_feature_vals = aa_encoder_.find( aa )->second; 
		else aa_feature_vals = aa_encoder_.find( 'X' )->second; 
		for( Size ival = 1; ival <= aa_feature_vals.size(); ++ival ){
			feature_vals.push_back( aa_feature_vals[ ival ] );
		}
	}
	return feature_vals;
}

//return feature vector with features averaged over sequence positions
vector1< Real >
NMerSVMEnergy::encode_wtd_avg_aa_string( std::string const seq, vector1< Real > const wts ) const
{
	assert( seq.length() == wts.size() );
	//get scale factors from relative wts, must sum to one!
	Real wtsum = 0;
	for( Size iwt = 1; iwt <= wts.size(); ++iwt ) wtsum += wts[ iwt ];
	vector1< Real > wts_norm( wts );
	for( Size iwt = 1; iwt <= wts_norm.size(); ++iwt ) wts_norm[ iwt ] /= wtsum;

	vector1< Real > feature_vals;
	for( Size iseq = 0; iseq <= seq.size() - 1; ++iseq ){
		std::string aa( seq.substr( iseq, 1 ) );
		vector1< Real > aa_feature_vals( encode_aa_string( aa ) );
		for( Size ival = 1; ival <= aa_feature_vals.size(); ++ival ){
			Real val_norm( wts_norm[ iseq + 1 ] * aa_feature_vals[ ival ] );
			//if vector element DNE, add wtd feature value
			if( ival > feature_vals.size() ) feature_vals.push_back( val_norm );
			//else sum it into the correct vector element
			else feature_vals[ ival ] += val_norm;
		}
	}
	return feature_vals;
}

//chain_seqpos is the P1 position of the epitope
vector1< Real >
NMerSVMEnergy::encode_nmer( std::string const chain_sequence, Size const chain_seqpos, Size const isvm ) const
{
	std::string const nmer_seq( chain_sequence.substr( chain_seqpos - 1, nmer_length_ ) );
	vector1< Real > feature_vals( encode_aa_string( nmer_seq ) );
	if( term_length_ > 0 ){
		add_encoded_termini( chain_sequence, chain_seqpos, feature_vals );
	}
	if( use_pssm_features_ ){
		add_pssm_features( nmer_seq, isvm, feature_vals );
	}
	return feature_vals;
}

//chain_seqpos is the P1 position of the epitope
void
NMerSVMEnergy::add_encoded_termini( std::string const chain_sequence, Size const chain_seqpos, vector1< Real > & feature_vals ) const
{
	//get term nmers upstream and downstrem of core, use dummy - for missing res
	std::string nterm_seq, cterm_seq;
	for( Size offset = 1; offset <= term_length_; ++offset ){
		if( chain_seqpos - offset >= 1 )
				nterm_seq = chain_sequence[ ( chain_seqpos - 1 ) - offset ] + nterm_seq;
		//else nterm_seq = '-' + nterm_seq;
		else nterm_seq = "X" + nterm_seq;
		if( chain_seqpos + offset <= chain_sequence.length() )
				cterm_seq = cterm_seq + chain_sequence[ ( chain_seqpos - 1 ) + ( nmer_length_ - 1 ) + offset ];
		//else cterm_seq = cterm_seq + '-';
		else cterm_seq = cterm_seq + "X";
	}
	//construct wts
	vector1< Real > nterm_wts( term_length_, 1.0 ), cterm_wts( term_length_, 1.0 );
	//currently hardcoded to weight linearly down as we go away from central core sequence
	for( Size iterm = 1; iterm <= term_length_; ++iterm ){
		nterm_wts[ iterm ] = static_cast< Real >( iterm );
		cterm_wts[ term_length_ - iterm + 1 ] = static_cast< Real >( iterm );
	}
	//get avged term seqs
	vector1< Real > feat_vals_nterm( encode_wtd_avg_aa_string( nterm_seq, nterm_wts ) );
	vector1< Real > feat_vals_cterm( encode_wtd_avg_aa_string( cterm_seq, cterm_wts ) );
	//add to term feature vals to begin and end of feature vals
	for( Size ival_term = 1; ival_term <= feat_vals_nterm.size(); ++ival_term ){
		feature_vals.insert( feature_vals.begin() + ival_term - 1, feat_vals_nterm[ ival_term ] );
		feature_vals.push_back( feat_vals_cterm[ ival_term ] );
	}
}

// make sure you prescale these pssms! use the same ones for training and here
// lets prescale them by logoddsing against backgroun distribution
void
NMerSVMEnergy::add_pssm_features( std::string const seq, Size const isvm, vector1< Real > & feature_vals ) const
{
	if( isvm > nmer_pssm_.n_pssms() ) utility_exit_with_message(
			"NMerSVMEnergy pssm features require one pssm per svm model. Check your svm and pssm lists for same length!\n" );
	for( Size iseq = 0; iseq <= seq.size() - 1; ++iseq ){
		//iter through encoding vector for each char of seq
		core::chemical::AA aa( core::chemical::aa_from_oneletter_code( seq[ iseq ] ) );
		Real pssm_energy( nmer_pssm_.pssm_energy_at_frame_seqpos( iseq + 1, aa, isvm ) );
		feature_vals.push_back( pssm_energy );
	}
}

Real
NMerSVMEnergy::get_residue_energy_by_svm(
	pose::Pose const & pose,
	Size const & seqpos,
	Real & rsd_energy_avg,
	vector1< Real > & rsd_svm_energies
) const
{
	assert( rsd_svm_energies.size() == n_svms() );
	//for now, just assign all of the p1=seqpos frame's nmer_svm energy to this residue
	//TODO: distribute frame's nmer_svm energy evenly across the nmer
	//TODO: avoid wasting calc time by storing nmer_value, nmer_val_out_of_date in pose cacheable data
	Size p1_seqpos( seqpos );

	//get the nmer string
	//TODO: how deal w/ sequences shorter than nmer_length_?
	// this matters at both termini… maybe take max of all overlapping frames w/ missing res as 'X'?
	// go ahead and bail if we fall off the end of the chain
	if( p1_seqpos + nmer_length_ - 1 <= pose.conformation().chain_end( pose.chain( p1_seqpos ) ) ){
		//need p1 position in this chain's sequence, offset with index of first res
		Size chain_p1_seqpos( p1_seqpos - pose.conformation().chain_begin( pose.chain( p1_seqpos ) ) + 1 );
		std::string chain_sequence( pose.chain_sequence( pose.chain( p1_seqpos ) ) );

		rsd_energy_avg = 0.;
		//encode nmer string --> feature(Svm_node) vector
		for( Size isvm = 1; isvm <= n_svms(); ++isvm ){
			vector1< Svm_node_rosettaOP > p1_seqpos_nmer_features(
					get_svm_nodes( encode_nmer( chain_sequence, chain_p1_seqpos, isvm ) ) );
			//get this svm model
			Svm_rosettaOP nmer_svm_model( all_nmer_svms_[ isvm ] );
			//Svm_rosetta funxn to wrap svm_predict, see Svm_rosetta::predict_probability for template, 
			Real rsd_svm_energy( nmer_svm_model->predict( p1_seqpos_nmer_features ) );
			//gate energy at svm_scorecut, thus ignoring low-scoring nmers
			if( gate_svm_scores_ && rsd_svm_energy < nmer_svm_scorecut_ ) rsd_svm_energy = nmer_svm_scorecut_;
			//store this svms rsd energy
			rsd_svm_energies[ isvm ] = rsd_svm_energy;
			//normalize energy by number of svms used
			//otherwise avg scores would become huge if we use lots of svms instead of just 1
			rsd_energy_avg += ( rsd_svm_energy / n_svms() );
		}
	}
	return rsd_energy_avg;
}

void
NMerSVMEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	if( all_nmer_svms_.empty() ) return;
	Size const seqpos( rsd.seqpos() );

	Real rsd_energy( 0. );
	vector1< Real > rsd_svm_energies( n_svms(), Real( 0. ) );
	get_residue_energy_by_svm( pose, seqpos, rsd_energy, rsd_svm_energies );
	emap[ nmer_svm ] += rsd_energy;
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

