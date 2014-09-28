// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Steven Lewis (smlewi@unc.edu)
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/pack/task/ResfileReader.hh>

// Package Headers
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Project Headers
#include <core/chemical/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/assert.hh> //ASSERT_ONLY makes release build happy
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>

using basic::T;
using basic::Error;
using basic::Warning;

//STL headers
#include <string>
//#include <iostream> //need this for debugging
#include <fstream>
#include <sstream> //automatic checking string to int conversion
#ifdef WIN32
#include <cctype> //for split_lines to handle '\t' tab characters
#endif

#include <algorithm>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>




namespace core {
namespace pack {
namespace task {

/// @details Auto-generated virtual destructor
ResfileCommand::~ResfileCommand() {}

static thread_local basic::Tracer TR( "core.pack.task.ResfileReader" );

using std::string;
using std::endl;
using std::istream;
using std::istringstream;
using std::map;
using std::stringstream;
using utility::vector1;
using core::pose::Pose;

ResfileContents::ResfileContents(
	Pose const & pose,
	istream & resfile ) :
	commands_( pose.total_residue() )
{
	using namespace std;

	map< string, ResfileCommandOP > command_map = create_command_map();
	bool have_read_start_token = false;

	// save the RANGE and CHAIN commands and apply them at the end
	vector1< ResfileCommandOP > default_commands;
	vector1< std::list< ResfileCommandOP > > residue_range_commands(pose.total_residue(), std::list< ResfileCommandOP >());
	vector1< std::list< ResfileCommandOP > > residue_chain_commands(pose.total_residue(), std::list< ResfileCommandOP >());

	uint lineno = 0;
	while ( resfile ) {
		vector1< string > tokens( tokenize_line( resfile ));
		++lineno;

		if (!tokens.size()) continue;
		if (comment_begin(tokens,1)) continue; // ignore the rest of this line

		if (!have_read_start_token) {
			parse_header_line(tokens, command_map, lineno, have_read_start_token);
		} else {
			parse_body_line(pose, tokens, command_map, lineno,
				residue_range_commands, residue_chain_commands);
		}
	}

	if (!have_read_start_token) {
		TR.Warning
			<< "Reached the end of resfile without finding a 'START' token." << endl
			<< "No residue-specific behavior specified in resfile." << endl;
	}

	// apply the RANGE and CHAIN commands
	for(Size i = 1; i <= pose.total_residue(); ++i){
		if(commands_[i].empty()){
			if(!residue_range_commands[i].empty()){
				commands_[i].assign(
					residue_range_commands[i].begin(), residue_range_commands[i].end());
			} else if(!residue_chain_commands[i].empty()){
				commands_[i].assign(
					residue_chain_commands[i].begin(), residue_chain_commands[i].end());
			}
		}
	}
}


ResfileContents::~ResfileContents() {}

std::list< ResfileCommandCOP > const &
ResfileContents::default_commands() const
{
	return default_commands_;
}


bool
ResfileContents::specialized_commands_exist_for_residue( Size resid ) const
{
	if ( resid > commands_.size() || resid == 0 ) {
		utility_exit_with_message( "Residue index out of bounds: resid=" +
			utility::to_string( resid ) + " with only " +
			utility::to_string( commands_.size() ) + " residues total." );
	}
	return ! commands_[ resid ].empty();
}

std::list< ResfileCommandCOP > const &
ResfileContents::commands_for_residue( Size resid ) const
{
	if ( resid > commands_.size() || resid == 0 ) {
		utility_exit_with_message( "Residue index out of bounds: resid=" +
			utility::to_string( resid ) + " with only " +
			utility::to_string( commands_.size() ) + " residues total." );
	}
	return commands_[ resid ];
}

void
ResfileContents::parse_header_line(
	vector1< string > const & tokens,
	map< string, ResfileCommandOP > const & command_map,
	Size const lineno,
	bool & have_read_start_token
) {
	// read in default behaviors, store them, process them later
	if ( get_token( 1, tokens) == "START" ) {
		have_read_start_token = true;
	} else {
		Size which_token = 1, ntokens = tokens.size();
		while(which_token <= ntokens){
			if (comment_begin(tokens, which_token)) break;

			ResfileCommandOP command(
				locate_command(which_token, tokens, command_map, lineno));
			command->initialize_from_tokens(tokens, which_token, 0);
			default_commands_.push_back(command);
		}
	}
}


///@ details Parse body line in resfile
/// expected formats:
///  <residue identifier> <chain identifier> <commands*>
///  <residue identifier> - <residue identifier> <chain identifier> <commands*>
///  * <chain identifer> <commands*>
///
/// The here is how a residue specification is resolved if it is
/// specified multiple times:
///
///  1) If a residue is specified by multiple commands of the same level (SINGLE, RANGE, or CHAIN), then this is an error
///  2) If a residue is specified by multiple commands of different levels then the more restricted level takes precidence, SINGLE over RANGE and CHAIN, and RANGE over CHAIN
void
ResfileContents::parse_body_line(
	pose::Pose const & pose,
	vector1< string > const & tokens,
	map< string, ResfileCommandOP > const & command_map,
	Size const lineno,
	vector1< std::list<ResfileCommandOP > > & residue_range_commands,
	vector1< std::list<ResfileCommandOP > > & residue_chain_commands
) {
	Size which_token = 1, ntokens(tokens.size());
	int PDBnum, PDBnum_end;
	char icode, icode_end, chain;
	residue_identifier_type id_type;

	parse_resid( which_token, tokens, lineno,
		PDBnum, PDBnum_end, icode, icode_end, chain, id_type);

	bool found_commands(false);

	if(id_type == ResfileContents::SINGLE_RESID){
		Size const resid(locate_resid(pose, chain, PDBnum, icode, lineno));

		while (which_token <= ntokens) {
			if (comment_begin(tokens, which_token)) break;

			ResfileCommandOP command(
				locate_command(which_token, tokens, command_map, lineno));
			command->initialize_from_tokens( tokens, which_token, resid );
			commands_[ resid ].push_back( command );
			found_commands = true;
		}
	} else if (id_type == ResfileContents::RANGE_RESID) {
		Size const resid_start(locate_resid(pose, chain, PDBnum, icode, lineno));
		Size const resid_end(locate_resid(pose, chain, PDBnum_end, icode_end, lineno));
		if(resid_start >= resid_end){
			stringstream err_msg;
			err_msg
				<< "On line " << lineno << ", "
				<< "the start residue (PDBnum=" << PDBnum;
			if(icode != ' '){
				err_msg << ", icode=" << icode;
			}
			err_msg
				<< ", chain=" << (chain == ' ' ? '_' : chain) << ") "
				<< "does not come before the end residue "
				<< "(PDBnum=" << PDBnum_end;
			if(icode_end != ' '){
				err_msg << ", icode=" << icode_end;
			}
			err_msg
				<< ", chain=" << (chain == ' ' ? '_' : chain) << ").";
			onError(err_msg.str());
		}

		while ( which_token <= ntokens ) {
			if (comment_begin(tokens, which_token)) break;

			ResfileCommandOP command(
				locate_command(which_token, tokens, command_map, lineno));
			Size const saved_which_token(which_token);
			Size which_token_i;
			// The number in pdb files is not straight
			for(Size i = resid_start; i <= resid_end; ++i){
				which_token_i = saved_which_token;
				ResfileCommandOP command_i(command->clone());
				command_i->initialize_from_tokens(tokens, which_token_i, i);
				residue_range_commands[i].push_back(command_i);
			}
			which_token = which_token_i;
			found_commands = true;
		}
	} else if (id_type == ResfileContents::CHAIN_RESID) {
		if(!pose.pdb_info()) {
			stringstream err_msg;
			err_msg
				<< "On line " << lineno << ", "
				<< "attemting to set task for all residues in a chain, however the pose does not have a PdbInfo object, which is used to look up which residues are part of which chain. Often the PdbInfo object is setup when the pose is intialized. Please check your protocol.";
			utility_exit_with_message(err_msg.str());
		}
		bool found_a_residue_on_chain(false);
		while ( which_token <= ntokens ) {
			if (comment_begin(tokens, which_token)) break;

			ResfileCommandOP command(
				locate_command(which_token, tokens, command_map, lineno));
			Size const saved_which_token(which_token);
			Size which_token_i;
			for(Size i=1; i <= pose.total_residue(); ++i){
				if(pose.pdb_info()->chain(i) == chain){
					found_a_residue_on_chain = true;
					which_token_i = saved_which_token;
					ResfileCommandOP command_i(command->clone());
					command_i->initialize_from_tokens(tokens, which_token_i, i);
					residue_chain_commands[i].push_back(command_i);
				}
			}
			which_token = which_token_i;
			found_commands = true;
		}
		if(!found_a_residue_on_chain){
			stringstream err_msg;
			err_msg
				<< "On line " << lineno << ", "
				<< "there are no residues with chain '"
				<< (chain == ' ' ? '_' : chain) << "'.";
			onError(err_msg.str());
		}

	} else {
		// unrecongized id type
		runtime_assert(false);
	}

	if(!found_commands){
		stringstream err_msg;
		err_msg
			<< "On line " << lineno << ", "
			<< "no commands are specified.";
		onError(err_msg.str());
	}
}


void
ResfileContents::parse_resid(
	Size & which_token,
	vector1< string > const & tokens,
	Size const lineno,
	int & PDBnum,
	int & PDBnum_end, // only defined for RANGE id_type
	char & icode,
	char & icode_end, // only defined for RANGE id_type
	char & chain,
	residue_identifier_type & id_type) const {

	string token(get_token(which_token, tokens));
	++which_token;
	if(*token.begin() == '*'){
		// apply task all residues on the chain
		id_type = ResfileContents::CHAIN_RESID;

		if(token.length() > 1){
			stringstream err_msg;
			err_msg
				<< "On line " << lineno << ", "
				<< "malformed residue identifier specification involving '*'";
			onError(err_msg.str());
		}

	} else {

		id_type = ResfileContents::SINGLE_RESID;

		// this may turn out to be a range specification but we have to
		// wait to the next token to figure that out.

		parse_PDBnum_icode(token, lineno, PDBnum, icode);
	}

	token = get_token(which_token, tokens);
	++which_token;

	if(token == "-"){

		if(id_type == ResfileContents::CHAIN_RESID){
			// In this case we've already encountered a '*' to it can't be a
			// range specification
			stringstream err_msg;
			err_msg
				<< "On line " << lineno << ", "
				<< "malformed residue identifier.";
			onError(err_msg.str());
		}

		id_type = ResfileContents::RANGE_RESID;

		token = get_token( which_token, tokens );
		++which_token;

		parse_PDBnum_icode(token, lineno, PDBnum_end, icode_end);

		// advance token counter in preparation for getting the chain
		// token "again".
		++which_token;
	}

	// The pdb format is case insensitive everwhere except for the chain
	// specifier! Get the token again but this time do not change it to
	// upper case
	token = get_token(which_token-1, tokens, false);
	if (token.length() != 1){
		stringstream err_msg;
		err_msg
			<< "On line " << lineno << ", "
			<< "the chain identifier '" << token << "' "
			<< "must be just a single character in [_A-Za-z] "
			<< "(note the chain identifier is case sensitive).";
		onError(err_msg.str());
	}
	chain = token[0];
	if (chain == '_') chain = ' ';
	if(core::chemical::chr_chains.find(chain) == std::string::npos
		 && chain != ' ' ){
		stringstream err_msg;
		err_msg
			<< "On line " << lineno << ", "
			<< "The chain identifier '" << chain << "' "
			<< "is not in [_A-Za-z] "
			<< "(note the chain identifier is case sensitive).";
		onError(err_msg.str());
	}
}

void
ResfileContents::parse_PDBnum_icode(
	string const & token,
	Size const lineno,
	int & PDBnum,
	char & icode) const {

	istringstream PDBnum_s;
	if ( std::isalpha( *token.rbegin() ) ) {
		PDBnum_s.str(token.substr(0, token.length() - 1));
		icode = *token.rbegin();
	} else {
		PDBnum_s.str(token);
		icode = ' ';
	}

	char remaining;
	if(!(PDBnum_s >> PDBnum) || PDBnum_s.get(remaining)){
		stringstream err_msg;
		err_msg
			<< "On line " << lineno << ", "
			<< "the token '" << token << "' "
			<< "is not a valid <PDBNUM>[<ICODE>] identifier.";
		onError(err_msg.str());
	}

}

Size
ResfileContents::locate_resid(
	Pose const & pose,
	char const chain,
	Size const PDBnum,
	char const icode,
	Size const lineno
) const {

	Size resid(0);
	if(pose.pdb_info() == 0){
		if(1 <= PDBnum && PDBnum <= pose.total_residue()){
			resid = PDBnum;
		}
	} else {
		resid = pose.pdb_info()->pdb2pose().find( chain, PDBnum, icode );
	}

	if(resid == 0){
		std::stringstream err_msg;
		err_msg  << "On line " << lineno << ", the pose does not have residue with chain=" << chain << ", PDBnum=" << PDBnum;
		if(icode != ' '){
			err_msg << ", icode=" << icode;
		}
		err_msg << ".";
		onError(err_msg.str());
	}
	return resid;
}



ResfileCommandOP
ResfileContents::locate_command(
	Size const which_token,
	vector1< string > const & tokens,
	std::map< string, ResfileCommandOP > const & command_map,
	Size const lineno
) const {

	std::map< string, ResfileCommandOP >::const_iterator command(
		command_map.find(get_token(which_token, tokens)));

	if (command == command_map.end()) {
		std::stringstream err_msg;
		err_msg
			<< "On line " << lineno
			<< " command '" << get_token(which_token, tokens) <<"' is not recognized.";
		onError(err_msg.str());
	}
	return(command->second->clone());
}


///////////////////////////////////////////////////////////////////////
///@brief NATRO disables packing and designing at a position, the residue
///will be totally unchanged
void
NATRO::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
NATRO::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).prevent_repacking();
}

///////////////////////////////////////////////////////////////////////
///@brief NATAA allows repacking but no sequence changes (all rotamers are of the original residue)
void
NATAA::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
NATAA::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).restrict_to_repacking();
}

///////////////////////////////////////////////////////////////////////
///@brief ALLAA is deprecated; allows repacking and designing to any canonical residue (default state of PackerTask)
void
ALLAA::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY( tokens ),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
ALLAA::residue_action(
	PackerTask & task,
	Size resid
) const
{
	//warn user not to use ALLAA to mean ALLAAxc.  This repeats for every line on purpose!
	// disable warning with command -nowarn_ALLAA
	// this warning is /really/ not necessary -- people ought to know that cys is an amino acid!
//	T("core.pack.task.ResfileReader") << "RESFILE WARNING: " << ALLAA::name() << " is deprecated.  Use "
//				<< ALLAAxc::name() << " to exclude cysteine, or "
//				<< ALLAAwc::name() << " to include cysteine.  Substituting "
//				<< ALLAAwc::name() << " behavior." << std::endl;

	//as in warning, pass to ALLAAwc
	ALLAAwc allaawc;
	allaawc.residue_action(task,resid);

}


///////////////////////////////////////////////////////////////////////
///@brief ALLAAxc allows repacking and designing to any canonical noncysteine residue
void
ALLAAxc::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
ALLAAxc::residue_action(
	PackerTask & task,
	Size resid
) const
{
	utility::vector1< bool > keep_aas( chemical::num_canonical_aas, true );
	keep_aas[ chemical::aa_cys ] = false;
	std::string const mode( "ALLAAxc" );
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_aas, mode );
}


///////////////////////////////////////////////////////////////////////
///@brief allows repacking and designing to any canonical residue (default state of PackerTask)
void
ALLAAwc::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() || get_token( which_token, tokens ) == ALLAA::name() );
	++which_token;
}

void
ALLAAwc::residue_action(
	PackerTask & /*task*/,
	Size /*resid*/
) const
{
	//do not change anything, default packerTask behavior is to design with everything
}


///////////////////////////////////////////////////////////////////////
///@brief PIKAA allows residues specifed in a following string and packing
///the string should be formatted ALLCAPS with no spaces between residues
///using the standard single letter codes
void
PIKAA::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	using namespace chemical;
	using utility::vector1;

	assert( get_token( which_token, tokens ) == name() );
	keep_canonical_aas_.resize( chemical::num_canonical_aas, false );

	if ( which_token == tokens.size() ) {
		std::stringstream err_msg;
		err_msg << "PIKAA must be followed by a string of allowed amino acids.";
		onError(err_msg.str());
	}
	++which_token;
	std::string const & aas_to_keep = get_token( which_token, tokens );


	// note: stl uses an index-by-0 convention so the for-loop initialization statment
	// and its boundary check are not in the rosetta standard.
	for ( Size ii = 0; ii < aas_to_keep.size(); ++ii ) {
		if ( oneletter_code_specifies_aa( aas_to_keep[ ii ] ) ) {
			AA aa( aa_from_oneletter_code( aas_to_keep[ ii ] ) );
			if ( size_t(aa) <= keep_canonical_aas_.size() ) {
				keep_canonical_aas_[ aa ] = true;
			}
			// allow use of PIKAA for DNA types (letters a,c,g,t)
			else if ( aa == na_ade || aa == na_cyt || aa == na_gua || aa == na_thy ) {
				na_allowed_.push_back( aa );
			}
		} else {
			TR << "Ignoring unknown one-letter amino acid code " << aas_to_keep[ ii ] << " while parsing PIKAA mode for residue " << resid << "." << std::endl;

		}
	}

	++which_token;
}

void
PIKAA::residue_action(
	PackerTask & task,
	Size resid
) const
{
	if ( keep_canonical_aas_.size() != chemical::num_canonical_aas ) {
		utility_exit_with_message( "PIKAA Resfile Command used uninitialized" );
	}
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_canonical_aas_ );

	// nucleic acids
	for ( std::list< chemical::AA >::const_iterator
			iter = na_allowed_.begin(), iter_end = na_allowed_.end();
			iter != iter_end; ++iter ) {
		task.nonconst_residue_task(resid).allow_noncanonical_aa( *iter );
	}
}

///////////////////////////////////////////////////////////////////////
///@brief PIKNA allows nucleic acid residues specifed in a following string
/// uses a string of single letter codes
void
PIKNA::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	using namespace chemical;
	using utility::vector1;

	assert( get_token( which_token, tokens ) == name() );

	if ( which_token == tokens.size() ) {
		std::stringstream err_msg;
		err_msg  << "PIKNA must be followed by nucleic acids in one-letter format.";
		onError(err_msg.str());
	}
	++which_token;
	std::string const & nas_string( get_token( which_token, tokens ) );

	// note: stl uses an index-by-0 convention so the for-loop initialization statment
	// and its boundary check are not in the rosetta standard.
	for ( std::string::const_iterator letter( nas_string.begin() );
		    letter != nas_string.end(); ++letter ) {
		// custom conversion from single letter to aa enum
		AA na( aa_unk );
		if      ( *letter == 'A' || *letter == 'a' ) na = na_ade;
		else if ( *letter == 'C' || *letter == 'c' ) na = na_cyt;
		else if ( *letter == 'G' || *letter == 'g' ) na = na_gua;
		else if ( *letter == 'T' || *letter == 't' ) na = na_thy;
		else {
			std::stringstream err_msg;
			TR << "Ignoring unknown one-letter nucleic acid code. " << *letter <<" while parsing PIKNA option for residue " << resid << ".";
			//onError(err_msg.str());
		}
		keep_nas_.push_back( na );
	}
	++which_token;
}

void
PIKNA::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).restrict_absent_nas( keep_nas_ );
}

///////////////////////////////////////////////////////////////////////
///@brief PIKNA allows nucleic acid residues specifed in a following string
/// uses a string of single letter codes
void
PIKRNA::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	using namespace chemical;
	using utility::vector1;

	assert( tokens[ which_token ] == name() );
	if ( which_token == tokens.size() ) {
		TR.Error << "RESFILE ERROR: PIKRNA must be followed by a string of allowed "
		        << "nucleic acids in single-letter format" << std::endl;
		utility_exit();
	}
	std::string const & nas_string( tokens[ ++which_token ] );

	// note: stl uses an index-by-0 convention so the for-loop initialization statment
	// and its boundary check are not in the rosetta standard.
	for ( std::string::const_iterator letter( nas_string.begin() );
		    letter != nas_string.end(); ++letter ) {
		// custom conversion from single letter to aa enum
		AA na( aa_unk );
		if      ( *letter == 'A' || *letter == 'a' ) na = na_rad;
		else if ( *letter == 'C' || *letter == 'c' ) na = na_rcy;
		else if ( *letter == 'G' || *letter == 'g' ) na = na_rgu;
		else if ( *letter == 'U' || *letter == 'u' ) na = na_ura;
		else {
			TR.Error << "RESFILE ERROR: unknown one-letter nucleic acid code " << *letter
				      << " while parsing PIKRNA option for residue " << resid << std::endl;
			utility_exit();
		}
		keep_rnas_.push_back( na );
	}

	++which_token;
}

void
PIKRNA::residue_action(
	PackerTask & task,
	Size resid
) const
{
	for ( Size ii = 1; ii <= keep_rnas_.size(); ++ii ) {
		task.nonconst_residue_task(resid).allow_aa( keep_rnas_[ ii ] );
	}
}

///////////////////////////////////////////////////////////////////////
///@brief NOTAA disallows residues specified in a following string, and allows packing
///the string should be formatted ALLCAPS with no spaces between residues
///using the standard single letter codes
void
NOTAA::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	using namespace chemical;

	assert( get_token( which_token, tokens ) == name() );
	keep_aas_.resize( chemical::num_canonical_aas, true );

	++which_token;
	std::string const & aas_to_exclude = get_token( which_token, tokens );

	// note: stl uses an index-by-0 convention so the for-loop initialization statment
	// and its boundary check are not in the rosetta standard.
	for ( Size ii = 0; ii < aas_to_exclude.size(); ++ii ) {
		if ( oneletter_code_specifies_aa( aas_to_exclude[ ii ] ) &&
				 aa_from_oneletter_code(aas_to_exclude[ii]) <= chemical::num_canonical_aas  ) {
			keep_aas_[ aa_from_oneletter_code( aas_to_exclude[ ii ] ) ] = false;
		} else {
			std::stringstream err_msg;
			TR << "Ignoring Unknown one-letter amino acid code "<< aas_to_exclude[ ii ] << " while parsing NOTAA option for residue " << resid << ".";
			//onError(err_msg.str());  // keep parsing on error
		}
	}

	++which_token;
}

void
NOTAA::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_aas_ );
}

///////////////////////////////////////////////////////////////////////
///@brief EMPTY disallows canonical residues but leaves packing and designing unchanged
///this is intended for use with noncanonical residues
///it will act like NOTAA QWERTYIPASDFGHKLCVNM (all residues), which essentially prevents repacking; PIKAA with no argument raises error
void
EMPTY::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
EMPTY::residue_action(
	PackerTask & task,
	Size resid
) const
{
	//vector is expected format for PackerTask, but false at all positions
	utility::vector1< bool > keep_aas( chemical::num_canonical_aas, false );
	std::string mode( "EMPTY" );
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_aas, mode );
	task.nonconst_residue_task(resid).disallow_noncanonical_aas();

}

///////////////////////////////////////////////////////////////////////
///@brief RESET disallows noncanonical residues and enables all of the canonical
///this is intended for use when both NC and PIKAA actions are used to allow for noncanonical and canonical residue at the same position
void
RESET::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
RESET::residue_action(
	PackerTask & task,
	Size resid
) const
{
	//vector is expected format for PackerTask, but false at all positions
	utility::vector1< bool > keep_aas( chemical::num_canonical_aas, true );
	std::string mode( "RESET" );
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_aas, mode );
	task.nonconst_residue_task(resid).disallow_noncanonical_aas();
}

///////////////////////////////////////////////////////////////////////
///@brief POLAR allows polar residues and packing
///polar-ness is ultimately determined in residue parameter file
void
POLAR::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
POLAR::residue_action(
	PackerTask & task,
	Size resid
) const
{
	using namespace chemical;

	utility::vector1< bool > keep_aas( chemical::num_canonical_aas, false );

	for ( ResidueLevelTask::ResidueTypeCOPListConstIter
			restype_iter = task.residue_task( resid ).allowed_residue_types_begin(),
			restype_iter_end = task.residue_task( resid ).allowed_residue_types_end();
			restype_iter != restype_iter_end; ++restype_iter ) {
		if ( (*restype_iter)->aa() > num_canonical_aas ) {
			std::stringstream err_msg;
			err_msg  << "POLAR mode read for residue " << resid << " which has been instructed to use non-canonical amino acids.";
			onError(err_msg.str());
			continue;
		}
		if ( (*restype_iter)->is_polar() ) {
			keep_aas[ (*restype_iter)->aa() ] = true;
		}
	}
	std::string mode( "POLAR");
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_aas, mode );

}

///////////////////////////////////////////////////////////////////////
///@brief APOLAR allows nonpolar residues and packing
///apolarity is (ultimately) determined by the lack of a POLAR flag in the residue parameter file
void
APOLAR::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;

}
void
APOLAR::residue_action(
	PackerTask & task,
	Size resid
) const
{
	using namespace chemical;

	utility::vector1< bool > keep_aas( chemical::num_canonical_aas, false );

	for ( ResidueLevelTask::ResidueTypeCOPListConstIter
			restype_iter = task.residue_task( resid ).allowed_residue_types_begin(),
			restype_iter_end = task.residue_task( resid ).allowed_residue_types_end();
			restype_iter != restype_iter_end; ++restype_iter ) {
		if ( (*restype_iter)->aa() > num_canonical_aas ) {
			std::stringstream err_msg;
			err_msg  << "APOLAR mode read for residue " << resid << " which has been instructed to use non-canonical amino acids.";
			onError(err_msg.str());
			continue;
		}
		if ( ! (*restype_iter)->is_polar() ) {
			keep_aas[ (*restype_iter)->aa() ] = true;
		}
	}
	std::string mode( "APOLAR" );
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_aas, mode );

}

///////////////////////////////////////////////////////////////////////
///@brief APOLA is deprecated, it calls APOLAR to allow nonpolar residues and packing
void
APOLA::initialize_from_tokens(
#ifdef NDEBUG
							  utility::vector1< std::string > const & ,
#else
							  utility::vector1< std::string > const & tokens,
#endif
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	TR << "RESFILE NOTE: APOLA command deprecated.  Use APOLAR command instead.  Treating as APOLAR command." << std::endl;
	++which_token;
}

void
APOLA::residue_action(
	PackerTask & task,
	Size resid
) const
{
	APOLAR apolar_command;
	apolar_command.residue_action( task, resid );
}

///////////////////////////////////////////////////////////////////////
///@brief EX handles extrachi options.  one EX command is necessary for each
///chi and sampling level you wish to turn on, so multiple EX commands may
///appear on a line.  EX must be followed by an integer (which chi)
///EX recognizes an optional subcommand LEVEL following the chi integer
///LEVEL must be followed by a second integer for the level you want
void
EX::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	assert( get_token( which_token, tokens ) == name() );

	if ( which_token == tokens.size() ) {
		std::stringstream err_msg;
		err_msg  << "EX command must be followe by a chi-id or \"ARO\" for reside " << resid << ".";
		onError(err_msg.str());
	}

	aro_specified_ = false;
	if ( get_token( which_token + 1, tokens ) == "ARO" ) {
		aro_specified_ = true;
		++which_token;
		if ( which_token == tokens.size() ) {
			std::stringstream err_msg;
			err_msg  << "EX ARO command must be followed by a chi-id for residue " << resid << ".";
			onError(err_msg.str());
		}
	}
	++which_token;
	which_chi_ = atoi( get_token( which_token, tokens ).c_str() );
	Size which_chi_token = which_token;
	chi_sample_level_ = EX_ONE_STDDEV;
	if ( which_token != tokens.size() ) {
		if ( get_token( which_token + 1, tokens ) == "LEVEL" ) {
			++which_token;

			++which_token;
			chi_sample_level_ = static_cast< ExtraRotSample> (atoi( get_token( which_token, tokens ).c_str() ));
			if ( chi_sample_level_ >= ExtraRotSampleCardinality || chi_sample_level_ < 0 ) {
				std::stringstream err_msg;
				err_msg  << "Extra rotamer sample level " << get_token( which_token, tokens ) << " is not in the range [0-" << ExtraRotSampleCardinality <<"] for residue " << resid << ".";
				onError(err_msg.str());
			}
		}
	}

	if ( which_chi_ > 2 && aro_specified_ ) {
		std::stringstream err_msg;
		err_msg  << "The token following EX ARO, " << get_token( which_chi_token, tokens ) << ", should be a valid chi-id level, eg a '1' or a '2' for residue " << resid << ".";
		onError(err_msg.str());
	}
	else if ( which_chi_ > 4 ) {
		std::stringstream err_msg;
		err_msg  << "The given chi-id, '" << get_token( which_chi_token, tokens ) << "' must either be an integer [1-4] or 'ARO' for residue " << resid << ".";
		onError(err_msg.str());
	}
	++which_token;
}

void
EX::residue_action(
	PackerTask & task,
	Size resid
) const
{
	if ( which_chi_ == 1 ) {
		if ( aro_specified_ ) {
			task.nonconst_residue_task( resid ).or_ex1aro_sample_level( chi_sample_level_ );
		} else {
			task.nonconst_residue_task( resid ).or_ex1_sample_level( chi_sample_level_ );
		}
	} else if ( which_chi_ == 2 ) {
		if ( aro_specified_ ) {
			task.nonconst_residue_task( resid ).or_ex2aro_sample_level( chi_sample_level_ );
		} else {
			task.nonconst_residue_task( resid ).or_ex2_sample_level( chi_sample_level_ );
		}
	} else if ( which_chi_ == 3 ) {
		task.nonconst_residue_task( resid ).or_ex3_sample_level( chi_sample_level_ );
	} else if ( which_chi_ == 4 ) {
		task.nonconst_residue_task( resid ).or_ex4_sample_level( chi_sample_level_ );
	}
}

////////////////////////////////////////////////////////////////////
///@brief NC allows a noncanonical residue; use one NC command per noncanonical.
/// The "nc_to_include_" string should match the interchangeability_group of
/// your desired residue type, and the residue type(s) in that group with
/// matching variants will be added to the PackerTask.
void
NC::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
	nc_to_include_ = get_token( which_token, tokens );
	++which_token;
}//end NC

void
NC::residue_action(
	PackerTask & task,
	Size resid
) const
{
	core::chemical::ResidueTypeSet const & residue_set = task.residue_task( resid ).get_original_residue_set();
	if ( residue_set.interchangeability_group_map( nc_to_include_ ).size() != 0 ){
		task.nonconst_residue_task(resid).allow_noncanonical_aa( nc_to_include_ );
	}	else {
		std::stringstream err_msg;
		err_msg  << "Unable to add non-canonical amino acid(s) with interchangeability group " << nc_to_include_ << " because there are no ResidueTypes with that interchangeability group in the ResidueTypeSet for residue " << resid << ".";
		onError(err_msg.str());
	}
}//end NC

/////////////////////////////////////////////////////////////////
///@brief EX_CUTOFF allows setting of the extrachi_cutoff (for determining burial for extra rotamers)
void
EX_CUTOFF::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	assert( get_token( which_token, tokens ) == name() );

	//check that next token is an integer value
	//stringstreams provide type safety!  yay!
	++which_token;
	std::istringstream ex_cutoff_ss( get_token( which_token, tokens ) );

	if (( ex_cutoff_ss >> ex_cutoff_ ).fail() ) { //convert string to Size and check error
		std::stringstream err_msg;
		err_msg << "The leve given for EX_CUTOFF, " << get_token( which_token-1, tokens) << ", must be an integer in the range [0-18] where 18 is the default for residue " << resid << ".";
		onError(err_msg.str());
	}//end error if

	++which_token;
}//end EX_CUTOFF

void
EX_CUTOFF::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).and_extrachi_cutoff(ex_cutoff_);
}//end EX_CUTOFF

///////////////////////////////////////////////////////////
///@brief USE_INPUT_SC turns on inclusion of the current rotamer for the packer
void
USE_INPUT_SC::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}//end USE_INPUT_SC

void
USE_INPUT_SC::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).or_include_current(true);
}//end USE_INPUT_SC

///////////////////////////////////////////////////////////
///@brief AUTO
void
AUTO::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
AUTO::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).add_behavior( name() );
}//end AUTO

///////////////////////////////////////////////////////////
///@brief SCAN
void
SCAN::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
SCAN::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).add_behavior( name() );
}
//end SCAN

///////////////////////////////////////////////////////////
///@brief TARGET
void
TARGET::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	assert( get_token( which_token, tokens ) == name() );

	if ( which_token == tokens.size() ) {
		std::stringstream err_msg;
		err_msg  << "TARGET must be followed by a character or string that identifies a target type.  Use '_' if no target type is desired; for residue " << resid << ".";
		onError(err_msg.str());
	}
	++which_token;
	argstring_ = get_token( which_token, tokens );
	++which_token;
}
void
TARGET::residue_action(
	PackerTask & task,
	Size resid
) const
{
	using namespace chemical;
	using utility::vector1;

	ResidueLevelTask & rtask( task.nonconst_residue_task(resid) );
	task.nonconst_residue_task(resid).add_behavior( name() );

	if ( argstring_.size() == 1 ) {
		char letter( argstring_[0] );
		if ( letter == '_' ) { return; } // dummy character
		else if ( oneletter_code_specifies_aa( letter ) && aa_from_oneletter_code( letter ) < chemical::num_canonical_aas  ) {
			rtask.target_type( aa_from_oneletter_code( letter ) );
		} else {
			std::stringstream err_msg;
			err_msg  << "Unknown one-letter code '"<< letter << "' for TARGET command for residue " << resid << ".";
			onError(err_msg.str());
		}
	}
	else rtask.target_type( argstring_ );
}
//end TARGET

///////////////////////////////////////////////////////////
///@brief NO_ADDUCTS
void
NO_ADDUCTS::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
NO_ADDUCTS::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).or_adducts(false);
}
//end NO_ADDUCTS

///////////////////////////////////////////////////////////
///@brief FIX_HIS_TAUTOMER
void
FIX_HIS_TAUTOMER::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	assert( tokens[ which_token ] == name() );
	++which_token;
}

void
FIX_HIS_TAUTOMER::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task(resid).or_fix_his_tautomer(true); //call is safe against not-histidine
}
//end FIX_HIS_TAUTOMER

///@details this creates a map linking the parsed strings from the resfile
///to the command objects.  NEW COMMANDS MUST BE ADDED HERE, HARD CODED
///note that it uses the command object name() method, not hard coded strings
///(of course, the name() method has hard coded strings...)
std::map< std::string, ResfileCommandOP >
create_command_map()
{
	using namespace std;

	map< string, ResfileCommandOP > command_map;
	command_map[ NATRO::name() ] = ResfileCommandOP( new NATRO );
	command_map[ NATAA::name() ] = ResfileCommandOP( new NATAA );
	command_map[ ALLAA::name() ] = ResfileCommandOP( new ALLAA );
	command_map[ ALLAAxc::name() ] = ResfileCommandOP( new ALLAAxc );
	command_map[ ALLAAwc::name() ] = ResfileCommandOP( new ALLAAwc );
	command_map[ PIKAA::name() ] = ResfileCommandOP( new PIKAA );
	command_map[ PIKNA::name() ] = ResfileCommandOP( new PIKNA );
	command_map[ PIKRNA::name() ] = ResfileCommandOP( new PIKRNA );
	command_map[ NOTAA::name() ] = ResfileCommandOP( new NOTAA );
	command_map[ EMPTY::name() ] = ResfileCommandOP( new EMPTY );
	command_map[ RESET::name() ] = ResfileCommandOP( new RESET );
	command_map[ POLAR::name() ] = ResfileCommandOP( new POLAR );
	command_map[ APOLAR::name() ] = ResfileCommandOP( new APOLAR );
	command_map[ APOLA::name() ] = ResfileCommandOP( new APOLA );
	command_map[ EX::name() ] = ResfileCommandOP( new EX );
	command_map[ EX_CUTOFF::name() ] = ResfileCommandOP( new EX_CUTOFF );
	command_map[ NC::name() ] = ResfileCommandOP( new NC );
	command_map[ USE_INPUT_SC::name() ] = ResfileCommandOP( new USE_INPUT_SC );
	command_map[ AUTO::name() ] = ResfileCommandOP( new AUTO );
	command_map[ SCAN::name() ] = ResfileCommandOP( new SCAN );
	command_map[ TARGET::name() ] = ResfileCommandOP( new TARGET );
	command_map[ NO_ADDUCTS::name() ] = ResfileCommandOP( new NO_ADDUCTS );
	command_map[ FIX_HIS_TAUTOMER::name() ] = ResfileCommandOP( new FIX_HIS_TAUTOMER );

	return command_map;
}
///@brief utility function for resfile reader (checks for a leading # signaling a comment)
bool
comment_begin( utility::vector1< std::string > const & tokens, Size which_token )
{
	return get_token( which_token, tokens )[ 0 ] == '#';
}

///@details resfile parser applies a resfile filename to a PackerTask
///each line of the resfile is broken into whitespace-delimited tokens
///whenever it reads a comment token, it ignores the rest of the line
///commands read before a "start" token are stored for application
///later as defaults lines after the start token alter specific
///ResidueLevelTasks in the PackerTask currently the reader assumes
///residue ID = PDB ID at the end, any residues not explicitly
///specified have the default actions performed on them

void
parse_resfile(
	pose::Pose const & pose,
	PackerTask & the_task)
{
	parse_resfile(pose, the_task, basic::options::option[basic::options::OptionKeys::packing::resfile].value().at(1));
}


void
parse_resfile(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string filename)
{

	std::string resfile;
	utility::io::izstream file( filename );
	if (!file) {
		TR << "File:" << filename << " not found!\n";
		utility_exit_with_message( "Cannot open file " + filename );
	} else {
		//T() << "read file: " << filename << "\n";
	}
	utility::slurp( file, resfile );
	try{
		parse_resfile_string( pose, the_task, resfile );
	} catch (ResfileReaderException &) {
		if (basic::options::option[ basic::options::OptionKeys::run::interactive ].user()){
			throw;
		} else {
			utility_exit();
		}
	}
}




	// question: if no chain is supplied should it be accepted?
	// yes just pass ' ' for the chain
	// how if a symbol is a chain or not?
	// all commands begin with something in the command map, if it's not a command treat it as a chain

void
parse_resfile_string(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string const & resfile_string ) throw(ResfileReaderException)
{
	using namespace std;
	istringstream resfile(resfile_string);
	ResfileContents contents( pose, resfile );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		std::list< ResfileCommandCOP > const & ii_command_list(
			contents.specialized_commands_exist_for_residue( ii ) ?
			contents.commands_for_residue( ii ) : contents.default_commands() );

		for ( std::list< ResfileCommandCOP >::const_iterator
				iter = ii_command_list.begin(), iter_end = ii_command_list.end();
				iter != iter_end; ++iter ) {
			(*iter)->residue_action( the_task, ii );
		}
	}
}


/// @details utility function to increment next token to be parsed
///
/// The PDB specification is case insensitive everywhere except for
/// the chain identifiers. By default get_tokens makes everything
/// upper case. To handle this special case, get_tokens optionally
/// allows it to not convert it to upper case.
std::string
get_token(
	const Size which_token,
	const utility::vector1<std::string> & tokens,
	bool make_upper_case) {

	if (which_token >  tokens.size() ){
		if (which_token == 1){
			std::stringstream err_msg;
			err_msg  << "Resfile does not specify anything.";
			onError(err_msg.str());
		} else {
			std::stringstream err_msg;
			err_msg  << "After token " << tokens[ which_token -1 ] << " another token is expected.";
			onError(err_msg.str());
		}
	}
	std::string token = tokens[ which_token ];
	if(make_upper_case){

		std::transform(token.begin(), token.end(), token.begin(), (int(*)(int)) std::toupper);
	}
	return token;
}

utility::vector1< std::string >
tokenize_line( std::istream & inputstream )
{

	utility::vector1< std::string > tokens;
	std::string input_line;
	std::getline( inputstream, input_line );

	unsigned int llength = input_line.size();
	unsigned int processing = 0;
	unsigned int token_start_char = 0;

	while ( processing < llength ) {
		if ( std::isspace( input_line[ processing ] ) ) { // == ' ') {
			if ( !( std::isspace(input_line[ token_start_char ] ) ) ) { // != ' ' ) {
				std::string token = input_line.substr( token_start_char, processing - token_start_char);
				tokens.push_back( token );
			}
			++processing;
			token_start_char = processing;
		} else {
			++processing;
		}
	}
	if ( processing != token_start_char ) { // last token on the line
		std::string token = input_line.substr( token_start_char, processing - token_start_char + 1 );
		tokens.push_back( token );
	}

	/*
	for ( Size ii = 1; ii <= tokens.size(); ++ii ) {
		std::cout << "token: " << ii << ": \"" << tokens[ ii ] << "\"" << std::endl;
	}
	*/

	return tokens;
}

void onError( std::string message){
	static bool interactive = basic::options::option[ basic::options::OptionKeys::run::interactive ].value();
	TR << message << std::endl;
	if (interactive){
		throw ResfileReaderException(message);
	} else {
		utility_exit();
	}
}


} //namespace task
} //namespace pack
} //namespace core
