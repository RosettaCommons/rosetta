// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Steven Lewis (smlewi@gmail.com)
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/pack/task/ResfileReader.hh>

// Package Headers
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/ResidueIndexDescription.hh>

// Project Headers
#include <core/chemical/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>

// Utility Headers
#include <utility>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/assert.hh> //ASSERT_ONLY makes release build happy
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


//STL headers
#include <string>
//#include <iostream> //need this for debugging
#include <sstream> //automatic checking string to int conversion
#ifdef WIN32
#include <cctype> //for split_lines to handle '\t' tab characters
#endif

#include <algorithm>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>

#include <core/pack/task/ResidueLevelTask.hh> // AUTO IWYU For ResidueLevelTask, ResidueLevelTask::ResidueTy...


namespace core {
namespace pack {
namespace task {

/// @details Auto-generated virtual destructor
ResfileCommand::~ResfileCommand() = default;

static basic::Tracer TR( "core.pack.task.ResfileReader" );

using std::string;
using std::endl;
using std::istream;
using std::istringstream;
using std::map;
using std::stringstream;
using utility::vector1;
using core::pose::Pose;

// Moves fname
ResfileContents::ResfileContents(
	Pose const & pose,
	std::string const & fname,
	istream & resfile
) :
	fname_initialized_from_( fname ),
	commands_( pose.size() )
{
	using namespace std;

	map< string, ResfileCommandOP > command_map = create_command_map();
	bool have_read_start_token = false;

	// save the RANGE and CHAIN commands and apply them at the end
	vector1< ResfileCommandOP > default_commands;
	vector1< std::list< ResfileCommandOP > > residue_range_commands( pose.size(), std::list< ResfileCommandOP >());
	vector1< std::list< ResfileCommandOP > > residue_chain_commands( pose.size(), std::list< ResfileCommandOP >());

	uint lineno = 0;
	while ( resfile ) {
		vector1< string > tokens( tokenize_line( resfile ));

		check_for_deprecated_commands( tokens );

		++lineno;

		if ( !tokens.size() ) continue;
		if ( comment_begin(tokens,1) ) continue; // ignore the rest of this line

		if ( !have_read_start_token ) {
			parse_header_line(tokens, command_map, lineno, have_read_start_token);
		} else {
			parse_body_line(pose, tokens, command_map, lineno,
				residue_range_commands, residue_chain_commands);
		}
	}

	if ( !have_read_start_token ) {
		TR.Warning
			<< "Reached the end of resfile without finding a 'START' token." << endl
			<< "No residue-specific behavior specified in resfile." << endl;
	}

	// apply the RANGE and CHAIN commands
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( commands_[i].empty() ) {
			if ( !residue_range_commands[i].empty() ) {
				commands_[i].assign(
					residue_range_commands[i].begin(), residue_range_commands[i].end());
			} else if ( !residue_chain_commands[i].empty() ) {
				commands_[i].assign(
					residue_chain_commands[i].begin(), residue_chain_commands[i].end());
			}
		}
	}
}


ResfileContents::~ResfileContents() = default;

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

std::string const &
ResfileContents::fname_initialized_from() const {
	return fname_initialized_from_;
}

/// @brief Given a vector of strings corresponding to the words in a whitespace-separated line, look for
/// deprecated commands and issue a suitable warning.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
ResfileContents::check_for_deprecated_commands(
	utility::vector1< std::string > const & tokens
) const {
	static const std::string errmsg( "Error in ResfileContents::check_for_deprecated_commands(): The \"" );
	static const std::string errmsg2( "\" command has been deprecated.  Please refer to the documentation on PackerPalettes (https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/PackerPalette) for details on how to properly use PackerPalettes to create a set of ResidueTypes with which one will design, and TaskOperations to prune away ResidueTypes that are not desired at a particular location.");
	static const utility::vector1< std::string > deprecated_tokens( { "NC", "EMPTY", "PIKRNA", "RESET" } );

	for ( std::string const & token : tokens ) {
		debug_assert( !token.empty() ); //Shouldn't be possible.
		if ( token[0] == '#' ) break; //We're done once we reach a comment.
		for ( std::string const & deprecated_token : deprecated_tokens ) {
			if ( token == deprecated_token ) {
				utility_exit_with_message( errmsg + deprecated_token + errmsg2  );
			}
		}
	}
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
		while ( which_token <= ntokens ) {
			if ( comment_begin(tokens, which_token) ) break;

			ResfileCommandOP command(
				locate_command(which_token, tokens, command_map, lineno));
			command->initialize_from_tokens(tokens, which_token, 0);
			default_commands_.push_back(command);
		}
	}
}


/// @ details Parse body line in resfile
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
	vector1< std::list< ResfileCommandOP > > & residue_range_commands,
	vector1< std::list< ResfileCommandOP > > & residue_chain_commands
) {
	Size which_token = 1, ntokens(tokens.size());
	int PDBnum, PDBnum_end;
	char icode, icode_end;
	std::string chain;
	residue_identifier_type id_type;

	parse_resid( which_token, tokens, lineno,
		PDBnum, PDBnum_end, icode, icode_end, chain, id_type);

	bool found_commands(false);

	if ( id_type == ResfileContents::SINGLE_RESID ) {
		Size const resid(locate_resid(pose, chain, PDBnum, icode, lineno));

		while ( which_token <= ntokens ) {
			if ( comment_begin(tokens, which_token) ) break;

			ResfileCommandOP command(
				locate_command(which_token, tokens, command_map, lineno));
			command->initialize_from_tokens( tokens, which_token, resid );
			commands_[ resid ].push_back( command );
			found_commands = true;
		}
	} else if ( id_type == ResfileContents::RANGE_RESID ) {
		Size const resid_start(locate_resid(pose, chain, PDBnum, icode, lineno));
		Size const resid_end(locate_resid(pose, chain, PDBnum_end, icode_end, lineno));
		if ( resid_start >= resid_end ) {
			stringstream err_msg;
			err_msg
				<< "On line " << lineno << ", "
				<< "the start residue (PDBnum=" << PDBnum;
			if ( icode != ' ' ) {
				err_msg << ", icode=" << icode;
			}
			err_msg
				<< ", chain=" << chain_printable(chain) << ") "
				<< "does not come before the end residue "
				<< "(PDBnum=" << PDBnum_end;
			if ( icode_end != ' ' ) {
				err_msg << ", icode=" << icode_end;
			}
			err_msg
				<< ", chain=" << chain_printable(chain) << ").";
			onError(err_msg.str());
		}

		while ( which_token <= ntokens ) {
			if ( comment_begin(tokens, which_token) ) break;

			ResfileCommandOP command(
				locate_command(which_token, tokens, command_map, lineno));
			Size const saved_which_token(which_token);
			Size which_token_i = -1; // Initialize to an obviously bogus value.
			// The number in pdb files is not straight
			for ( Size i = resid_start; i <= resid_end; ++i ) {
				which_token_i = saved_which_token;
				ResfileCommandOP command_i(command->clone());
				command_i->initialize_from_tokens(tokens, which_token_i, i);
				residue_range_commands[i].push_back(command_i);
			}
			debug_assert( which_token_i != core::Size(-1) );
			which_token = which_token_i;
			found_commands = true;
		}
	} else if ( id_type == ResfileContents::CHAIN_RESID ) {
		if ( !pose.pdb_info() ) {
			stringstream err_msg;
			err_msg
				<< "On line " << lineno << ", "
				<< "attemting to set task for all residues in a chain, however the pose does not have a PdbInfo object, which is used to look up which residues are part of which chain. Often the PdbInfo object is setup when the pose is intialized. Please check your protocol.";
			utility_exit_with_message(err_msg.str());
		}
		bool found_a_residue_on_chain(false);
		while ( which_token <= ntokens ) {
			if ( comment_begin(tokens, which_token) ) break;

			ResfileCommandOP command(
				locate_command(which_token, tokens, command_map, lineno));
			Size const saved_which_token(which_token);
			Size which_token_i;
			for ( Size i=1; i <= pose.size(); ++i ) {
				if ( pose.pdb_info()->chain(i) == chain ) {
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
		if ( !found_a_residue_on_chain ) {
			stringstream err_msg;
			err_msg
				<< "On line " << lineno << ", "
				<< "there are no residues with chain '"
				<< chain_printable(chain) << "'.";
			onError(err_msg.str());
		}

	} else {
		utility_exit_with_message("Unrecognized id type");
	}

	if ( !found_commands ) {
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
	std::string & chain,
	residue_identifier_type & id_type
) const {

	string token(get_token(which_token, tokens));
	++which_token;
	if ( *token.begin() == '*' ) {
		// apply task all residues on the chain
		id_type = ResfileContents::CHAIN_RESID;

		if ( token.length() > 1 ) {
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
		try {
			pose::parse_PDBnum_icode(token, fname_initialized_from_, lineno, PDBnum, icode);
		} catch ( utility::excn::Exception & e ) {
			onError( e.msg() );
		}
	}

	token = get_token(which_token, tokens);
	++which_token;

	if ( token == "-" ) {

		if ( id_type == ResfileContents::CHAIN_RESID ) {
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

		try {
			pose::parse_PDBnum_icode(token, fname_initialized_from_, lineno, PDBnum_end, icode_end);
		} catch ( utility::excn::Exception & e ) {
			onError( e.msg() );
		}

		// advance token counter in preparation for getting the chain
		// token "again".
		++which_token;
	}

	// The pdb format is case insensitive everwhere except for the chain
	// specifier! Get the token again but this time do not change it to
	// upper case
	token = get_token(which_token-1, tokens, false);
	chain = token;
	if ( chain == "_" ) chain = ' ';
}


Size
ResfileContents::locate_resid(
	Pose const & pose,
	std::string const & chain,
	Size const PDBnum,
	char const icode,
	Size const lineno
) const {

	Size resid(0);
	if ( pose.pdb_info() == nullptr ) {
		if ( 1 <= PDBnum && PDBnum <= pose.size() ) {
			resid = PDBnum;
		}
	} else {
		resid = pose.pdb_info()->pdb2pose().find( chain, PDBnum, icode );
	}

	if ( resid == 0 ) {
		std::stringstream err_msg;
		err_msg  << "On line " << lineno << ", the pose does not have residue with chain=" << chain << ", PDBnum=" << PDBnum;
		if ( icode != ' ' ) {
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

	auto command( command_map.find( get_token( which_token, tokens ) ) );

	if ( command == command_map.end() ) {
		std::stringstream err_msg;
		err_msg
			<< "On line " << lineno
			<< " command '" << get_token(which_token, tokens) <<"' is not recognized.";
		onError(err_msg.str());
	}
	return(command->second->clone());
}

std::string
ResfileContents::chain_printable(std::string const & chain) const {
	if ( chain.empty() || chain == " " ) {
		return "_";
	} else {
		return chain;
	}
}

///////////////////////////////////////////////////////////////////////
/// @brief NATRO disables packing and designing at a position, the residue
///will be totally unchanged
void
NATRO::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
NATRO::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task( resid ).prevent_repacking();
}

///////////////////////////////////////////////////////////////////////
/// @brief NATAA allows repacking but no sequence changes (all rotamers are of the original residue)
void
NATAA::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
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
/// @brief ALLAA is deprecated; allows repacking and designing to any canonical residue (default state of PackerTask)
void
ALLAA::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY( tokens ),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
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
	// T("core.pack.task.ResfileReader") << "RESFILE WARNING: " << ALLAA::name() << " is deprecated.  Use "
	//    << ALLAAxc::name() << " to exclude cysteine, or "
	//    << ALLAAwc::name() << " to include cysteine.  Substituting "
	//    << ALLAAwc::name() << " behavior." << std::endl;

	//as in warning, pass to ALLAAwc
	ALLAAwc allaawc;
	allaawc.residue_action(task,resid);
}


///////////////////////////////////////////////////////////////////////
/// @brief ALLAAxc allows repacking and designing to any canonical noncysteine residue
void
ALLAAxc::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
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
/// @brief allows repacking and designing to any canonical residue (default state of PackerTask)
void
ALLAAwc::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() || get_token( which_token, tokens ) == ALLAA::name() );
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
/// @brief PIKAA allows residues specifed in a following string.
/// @details In actuality, it is PROHIBITING any residue that is NOT in the
/// following string.  The string should be formatted as an all-caps string of
/// one-letter codes.  Noncanonical amino acids can be included using X[<full base name>].
/// For example, to allow tyrosine, threonine, tryptophan, and 2-aminoisobutyric acid,
/// you would use "PIKAA YTWX[AIB]".
/// @author Original author unknown.
/// @author Noncanonical pruning support added by Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PIKAA::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	using namespace chemical;
	using utility::vector1;

	static const std::string errmsg( "Error in core::pack::task::PIKAA::initialize_from_tokens(): " );

	debug_assert( get_token( which_token, tokens ) == name() );

	// Clear stored data:
	basenames_to_keep_.clear();

	runtime_assert_string_msg( which_token != tokens.size(), errmsg + "PIKAA must be followed by a string of allowed amino acids." );

	++which_token;
	std::string const & aas_to_keep = get_token( which_token, tokens );


	for ( core::Size i(0), imax(aas_to_keep.length()); i<imax; ++i ) {
		char const aa_to_keep( aas_to_keep[i] );
		if ( aa_to_keep == 'X' ) {
			// NCAAs are specified with "X[<full base name>]".
			std::string basename;
			++i;
			runtime_assert_string_msg( i != imax && aas_to_keep[i] == '[', errmsg + "An \"X\" in a \"PIKAA\" string must be followed by square brackets containing the full base name of a residue type." );
			std::string base_name;
			do {
				++i;
				runtime_assert_string_msg( i != imax , errmsg + "An opening bracket in a \"PIKAA\" string must be closed with a closing bracket." );
				char const curlet( aas_to_keep[i] );
				if ( curlet == ']' ) break;
				runtime_assert_string_msg( !( curlet == ' ' || curlet == '\t' || curlet == '\n' || curlet == '\0' ), errmsg + "A \"PIKAA\" string cannot contain whitespace." );
				base_name += curlet;
			} while(true);
			add_base_name_to_keep( base_name );
		} else if ( oneletter_code_specifies_aa( aa_to_keep ) ) {
			AA aa( aa_from_oneletter_code( aa_to_keep ) );
			if ( static_cast<Size>(aa) <= static_cast<Size>(core::chemical::num_canonical_aas ) ) {
				add_base_name_to_keep( core::chemical::name_from_aa( aa ) ); //Canonical base names are the same as name3s.
			} else if ( aa == na_ade || aa == na_cyt || aa == na_gua || aa == na_thy ) {
				// allow use of PIKAA for DNA types (letters a,c,g,t)
				add_base_name_to_keep( core::chemical::name_from_aa( aa ) );
			}
		} else {
			TR << "Ignoring unknown one-letter amino acid code " << aa_to_keep << " while parsing PIKAA mode for residue " << resid << "." << std::endl;
		}
	}

	if ( !initialized_ ) {
		TR.Warning << "Warning!  No sensible information was parsed from \"PIKAA " + tokens[which_token] + "\"." << std::endl;
	}

	++which_token;
}

/// @brief Add a base name to the list of base names to keep.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PIKAA::add_base_name_to_keep(
	std::string const & basename
) {
	static const std::string errmsg( "Error in core::pack::task::PIKAA::add_base_name_to_keep(): " );
	runtime_assert_string_msg( !basename.empty(), errmsg + "A base name cannot be an empty string." );
	runtime_assert_string_msg( !basenames_to_keep_.has_value( basename ), errmsg + "The base name \"" + basename + "\" was added more than once." );

	basenames_to_keep_.push_back(basename);
	initialized_ = true;
}

void
PIKAA::residue_action(
	PackerTask & task,
	Size resid
) const
{
	runtime_assert_string_msg( initialized_, "Error in core::pack::task::PIKAA::residue_action(): PIKAA resfile command used uninitialized." );
	task.nonconst_residue_task( resid ).restrict_restypes( basenames_to_keep_ );
}

///////////////////////////////////////////////////////////////////////
/// @brief PIKNA allows nucleic acid residues specifed in a following string
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

	debug_assert( get_token( which_token, tokens ) == name() );

	if ( which_token == tokens.size() ) {
		std::stringstream err_msg;
		err_msg  << "PIKNA must be followed by nucleic acids in one-letter format.";
		onError(err_msg.str());
	}
	++which_token;
	std::string const & nas_string( get_token( which_token, tokens ) );

	for ( char const letter : nas_string ) {
		// custom conversion from single letter to aa enum
		AA na( aa_unk );
		if      ( letter == 'A' || letter == 'a' ) na = na_ade;
		else if ( letter == 'C' || letter == 'c' ) na = na_cyt;
		else if ( letter == 'G' || letter == 'g' ) na = na_gua;
		else if ( letter == 'T' || letter == 't' ) na = na_thy;
		else {
			TR << "Ignoring unknown one-letter nucleic acid code. " << letter <<" while parsing PIKNA option for residue " << resid << ".";
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
/// @brief NOTAA disallows residues specified in a following string, and allows packing
///the string should be formatted ALLCAPS with no spaces between residues
///using the standard single letter codes
void
NOTAA::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	static const std::string errmsg( "Error in core::pack::task::NOTAA::initialize_from_tokens(): " );
	using namespace chemical;

	debug_assert( get_token( which_token, tokens ) == name() );
	basenames_to_exclude_.clear();

	++which_token;

	std::string const & aas_to_exclude = get_token( which_token, tokens );

	for ( core::Size i(0), imax(aas_to_exclude.length()); i<imax; ++i ) {
		char const aa_to_exclude( aas_to_exclude[i] );
		if ( aa_to_exclude == 'X' ) {
			// NCAAs are specified with "X[<full base name>]".
			std::string basename;
			++i;
			runtime_assert_string_msg( i != imax && aas_to_exclude[i] == '[', errmsg + "An \"X\" in a \"NOTAA\" string must be followed by square brackets containing the full base name of a residue type." );
			std::string base_name;
			do {
				++i;
				runtime_assert_string_msg( i != imax , errmsg + "An opening bracket in a \"NOTAA\" string must be closed with a closing bracket." );
				char const curlet( aas_to_exclude[i] );
				if ( curlet == ']' ) break;
				runtime_assert_string_msg( !( curlet == ' ' || curlet == '\t' || curlet == '\n' || curlet == '\0' ), errmsg + "A \"NOTAA\" string cannot contain whitespace." );
				base_name += curlet;
			} while(true);
			add_base_name_to_exclude( base_name );
		} else if ( oneletter_code_specifies_aa( aa_to_exclude ) ) {
			AA aa( aa_from_oneletter_code( aa_to_exclude ) );
			if ( static_cast<Size>(aa) <= static_cast<Size>(core::chemical::num_canonical_aas ) ) {
				add_base_name_to_exclude( core::chemical::name_from_aa( aa ) ); //Canonical base names are the same as name3s.
			} else if ( aa == na_ade || aa == na_cyt || aa == na_gua || aa == na_thy ) {
				// allow use of PIKAA for DNA types (letters a,c,g,t)
				add_base_name_to_exclude( core::chemical::name_from_aa( aa ) );
			}
		} else {
			TR << "Ignoring unknown one-letter amino acid code " << aa_to_exclude << " while parsing NOTAA mode for residue " << resid << "." << std::endl;
		}
	}

	if ( !initialized_ ) {
		TR.Warning << "Warning!  No sensible information was parsed from \"NOTAA " + tokens[which_token] + "\"." << std::endl;
	}

	++which_token;
}

void
NOTAA::residue_action(
	PackerTask & task,
	Size resid
) const {
	runtime_assert_string_msg( initialized_, "Error in core::pack::task::NOTAA::residue_action(): NOTAA resfile command used uninitialized." );
	task.nonconst_residue_task( resid ).disable_restypes( basenames_to_exclude_ );
}

/// @brief Add a base name to the list of base names to exclude.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
NOTAA::add_base_name_to_exclude(
	std::string const & basename
) {
	static const std::string errmsg( "Error in core::pack::task::NOTAA::add_base_name_to_exclude(): " );
	runtime_assert_string_msg( !basename.empty(), errmsg + "A base name cannot be an empty string." );
	runtime_assert_string_msg( !basenames_to_exclude_.has_value( basename ), errmsg + "The base name \"" + basename + "\" was added more than once." );

	basenames_to_exclude_.push_back(basename);
	initialized_ = true;
}

///////////////////////////////////////////////////////////////////////
/// @brief POLAR allows polar residues and packing
///polar-ness is ultimately determined in residue parameter file
void
PROPERTY::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
	++which_token;
	property_ = core::chemical::ResidueProperties::get_property_from_string( get_token( which_token, tokens ) );
	runtime_assert_string_msg( property_ != core::chemical::NO_PROPERTY, "Error in core::pack::task::PROPERTY::initialize_from_tokens(): The string \"" + get_token( which_token, tokens ) + "\" corresponds to no known residue property." );
	++which_token;
}

void
PROPERTY::residue_action(
	PackerTask & task,
	Size resid
) const
{
	using namespace chemical;

	runtime_assert_string_msg( property_ != core::chemical::NO_PROPERTY, "Error in core::pack::task::PROPERTY::residue_action(): The PROPERTY object cannot be used uninitialized." );

	task.nonconst_residue_task(resid).restrict_to_restypes_with_all_properties( utility::vector1< core::chemical::ResidueProperty >( { property_ } ) );
}



///////////////////////////////////////////////////////////////////////
/// @brief POLAR allows polar residues and packing
///polar-ness is ultimately determined in residue parameter file
void
POLAR::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
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

void
CHARGED::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
CHARGED::residue_action(
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
			err_msg  << "CHARGED mode read for residue " << resid << " which has been instructed to use non-canonical amino acids.";
			onError(err_msg.str());
			continue;
		}
		if ( (*restype_iter)->is_charged() ) {
			keep_aas[ (*restype_iter)->aa() ] = true;
		}
	}
	std::string mode( "CHARGED");
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_aas, mode );

}

void
AROMATIC::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
AROMATIC::residue_action(
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
			err_msg  << "AROMATIC mode read for residue " << resid << " which has been instructed to use non-canonical amino acids.";
			onError(err_msg.str());
			continue;
		}
		if ( (*restype_iter)->is_aromatic() ) {
			keep_aas[ (*restype_iter)->aa() ] = true;
		}
	}
	std::string mode( "AROMATIC");
	task.nonconst_residue_task(resid).restrict_absent_canonical_aas( keep_aas, mode );

}


///////////////////////////////////////////////////////////////////////
/// @brief APOLAR allows nonpolar residues and packing
///apolarity is (ultimately) determined by the lack of a POLAR flag in the residue parameter file
void
APOLAR::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
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
/// @brief APOLA is deprecated, it calls APOLAR to allow nonpolar residues and packing
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
	debug_assert( get_token( which_token, tokens ) == name() );
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
/// @brief EX handles extrachi options.  one EX command is necessary for each
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
	debug_assert( get_token( which_token, tokens ) == name() );

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
			chi_sample_level_ = static_cast< ExtraRotSample >( atoi( get_token( which_token, tokens ).c_str() ) );
			if ( chi_sample_level_ >= ExtraRotSampleCardinality ) { //|| chi_sample_level_ < 0 ) {
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
	} else if ( which_chi_ > 4 ) {
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

/////////////////////////////////////////////////////////////////
/// @brief EX_CUTOFF allows setting of the extrachi_cutoff (for determining burial for extra rotamers)
void
EX_CUTOFF::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	debug_assert( get_token( which_token, tokens ) == name() );

	//check that next token is an integer value
	//stringstreams provide type safety!  yay!
	++which_token;
	std::istringstream ex_cutoff_ss( get_token( which_token, tokens ) );

	if ( ( ex_cutoff_ss >> ex_cutoff_ ).fail() ) { //convert string to Size and check error
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
/// @brief USE_INPUT_SC turns on inclusion of the current rotamer for the packer
void
USE_INPUT_SC::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
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
/// @brief AUTO
void
AUTO::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
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
/// @brief SCAN
void
SCAN::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
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
/// @brief TARGET
void
TARGET::initialize_from_tokens(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	Size resid
)
{
	debug_assert( get_token( which_token, tokens ) == name() );

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
	} else rtask.target_type( argstring_ );
}
//end TARGET

///////////////////////////////////////////////////////////
/// @brief NO_ADDUCTS
void
NO_ADDUCTS::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( get_token( which_token, tokens ) == name() );
	++which_token;
}

void
NO_ADDUCTS::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task( resid ).or_adducts( false );
}
//end NO_ADDUCTS

///////////////////////////////////////////////////////////
/// @brief FIX_HIS_TAUTOMER
void
FIX_HIS_TAUTOMER::initialize_from_tokens(
	utility::vector1< std::string > const & ASSERT_ONLY(tokens),
	Size & which_token,
	Size /*resid*/
)
{
	debug_assert( tokens[ which_token ] == name() );
	++which_token;
}

void
FIX_HIS_TAUTOMER::residue_action(
	PackerTask & task,
	Size resid
) const
{
	task.nonconst_residue_task( resid ).or_fix_his_tautomer( true ); //call is safe against not-histidine
}
//end FIX_HIS_TAUTOMER

/// @details this creates a map linking the parsed strings from the resfile
///to the command objects.  NEW COMMANDS MUST BE ADDED HERE, HARD CODED
///note that it uses the command object name() method, not hard coded strings
///(of course, the name() method has hard coded strings...)
std::map< std::string, ResfileCommandOP >
create_command_map()
{
	using namespace std;

	map< string, ResfileCommandOP > command_map;
	command_map[ NATRO::name() ] = utility::pointer::make_shared< NATRO >();
	command_map[ NATAA::name() ] = utility::pointer::make_shared< NATAA >();
	command_map[ ALLAA::name() ] = utility::pointer::make_shared< ALLAA >();
	command_map[ ALLAAxc::name() ] = utility::pointer::make_shared< ALLAAxc >();
	command_map[ ALLAAwc::name() ] = utility::pointer::make_shared< ALLAAwc >();
	command_map[ PIKAA::name() ] = utility::pointer::make_shared< PIKAA >();
	command_map[ PIKNA::name() ] = utility::pointer::make_shared< PIKNA >();
	command_map[ NOTAA::name() ] = utility::pointer::make_shared< NOTAA >();
	command_map[ PROPERTY::name() ] = utility::pointer::make_shared< PROPERTY >();
	command_map[ POLAR::name() ] = utility::pointer::make_shared< POLAR >();
	command_map[ APOLAR::name() ] = utility::pointer::make_shared< APOLAR >();
	command_map[ APOLA::name() ] = utility::pointer::make_shared< APOLA >();
	command_map[ AROMATIC::name() ] = utility::pointer::make_shared< AROMATIC >();
	command_map[ CHARGED::name() ] = utility::pointer::make_shared< CHARGED >();
	command_map[ EX::name() ] = utility::pointer::make_shared< EX >();
	command_map[ EX_CUTOFF::name() ] = utility::pointer::make_shared< EX_CUTOFF >();
	command_map[ USE_INPUT_SC::name() ] = utility::pointer::make_shared< USE_INPUT_SC >();
	command_map[ AUTO::name() ] = utility::pointer::make_shared< AUTO >();
	command_map[ SCAN::name() ] = utility::pointer::make_shared< SCAN >();
	command_map[ TARGET::name() ] = utility::pointer::make_shared< TARGET >();
	command_map[ NO_ADDUCTS::name() ] = utility::pointer::make_shared< NO_ADDUCTS >();
	command_map[ FIX_HIS_TAUTOMER::name() ] = utility::pointer::make_shared< FIX_HIS_TAUTOMER >();

	return command_map;
}
/// @brief utility function for resfile reader (checks for a leading # signaling a comment)
bool
comment_begin( utility::vector1< std::string > const & tokens, Size which_token )
{
	return get_token( which_token, tokens )[ 0 ] == '#';
}

/// @details resfile parser applies a resfile filename to a PackerTask
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
	parse_resfile( pose, the_task, basic::options::option[ basic::options::OptionKeys::packing::resfile ].value().at( 1 ) );
}


void
parse_resfile(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string const & filename)
{
	std::string resfile;
	utility::io::izstream file( filename );
	if ( !file ) {
		TR << "File:" << filename << " not found!\n";
		utility_exit_with_message( "Cannot open file " + filename );
	}
	utility::slurp( file, resfile );
	file.close();
	try {
		parse_resfile_string( pose, the_task, filename, resfile );
	} catch ( ResfileReaderException & ) {
		if ( basic::options::option[ basic::options::OptionKeys::run::interactive ].value() ) {
			throw;
		} else {
			utility_exit();
		}
	}
}


/// @brief changes the state of the given PackerTask according to the commands in the resfile.
/// @details This version calls the overloaded parse_resfile_string and just passes it a ResidueSubset
/// that's set to "true" for every residue of the pose.
void
parse_resfile_string(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string const & resfile_fname,
	std::string const & resfile_string
)
{
	core::select::residue_selector::ResidueSubset const mask( pose.size(), true );
	parse_resfile_string( pose, the_task, resfile_fname, resfile_string, mask );
	return;
}

/// @brief changes the state of the given PackerTask according to the commands in the resfile.
/// @details This version accepts a ResidueSubset (a utility::vector1<bool>) that serves as a
/// mask.  Residues set to "false" don't have their packer behaviour altered by this TaskOperation.
// question: if no chain is supplied should it be accepted?
// yes just pass ' ' for the chain
// how if a symbol is a chain or not?
// all commands begin with something in the command map, if it's not a command treat it as a chain
void
parse_resfile_string(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string const & resfile_fname,
	std::string const & resfile_string,
	core::select::residue_selector::ResidueSubset const & mask
)
{
	using namespace std;
	istringstream resfile( resfile_string );
	ResfileContents contents( pose, resfile_fname, resfile );

	runtime_assert_string_msg( mask.size() == pose.size(), "Error in core::pack::task::parse_resfile_string(): The mask passed to this function is not the same size as the pose.  The mask must have one entry for every residue." );

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {

		if ( !mask[ii] ) continue; //Skip masked residues.

		std::list< ResfileCommandCOP > const & ii_command_list(
			contents.specialized_commands_exist_for_residue( ii ) ?
			contents.commands_for_residue( ii ) : contents.default_commands() );

		for ( auto const & iter : ii_command_list ) {
			iter->residue_action( the_task, ii );
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
	bool make_upper_case
) {

	if ( which_token > tokens.size() ) {
		if ( which_token == 1 ) {
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
	if ( make_upper_case ) {
		std::transform( token.begin(), token.end(), token.begin(), (int(*)(int)) std::toupper);
	}
	return token;
}

utility::vector1< std::string >
tokenize_line( std::istream & inputstream )
{

	std::string input_line;
	std::getline( inputstream, input_line );

	return tokenize_line( input_line );
}

utility::vector1< std::string >
tokenize_line( std::string const & input_line){

	utility::vector1< std::string > tokens;

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

void onError( std::string message ){
	TR << message << std::endl;
	throw CREATE_EXCEPTION(ResfileReaderException, message);
}


} //namespace task
} //namespace pack
} //namespace core
