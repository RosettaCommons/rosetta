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
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>

// Project Headers
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



namespace core {
namespace pack {
namespace task {


ResfileContents::ResfileContents( PackerTask const & the_task, std::istream & resfile ) :
	commands_( the_task.total_residue() )
{
	using namespace std;

	map< string, ResfileCommandOP > command_map = create_command_map();
	bool have_read_start_token = false;

	utility::vector1< bool > non_default_lines( the_task.total_residue(), false );
	utility::vector1< std::string > default_tokens;
	utility::vector1< Size > origin_lines_of_default_tokens;

	uint lineno = 0;
	while ( resfile ) {
		utility::vector1< string > tokens( tokenize_line( resfile ));
		++lineno;

		// for debug
		//std::cout << "line->";
		//for( Size i=1; i <= tokens.size(); i++){
		//	std::cout << tokens[ i ] << ", ";
		//}
		//std::cout << std::endl;

		Size ntokens( tokens.size() );
		if ( ntokens == 0 ) continue;
		if ( comment_begin( tokens, 1 ) ) continue; // ignore the rest of this line

		if ( have_read_start_token ) {
			Size which_token = 1;

			// expected format: <residue identifier> <chain identifier> <commands*>
			//( the res/chain combo is used to get the pose's resid)

			// PDB numbering can be negative
			std::string const PDBnum_token = get_token( which_token, tokens );
			int PDBnum;
			char icode = ' ';
			if ( std::isalpha( *PDBnum_token.rbegin() ) ) {
				PDBnum = atoi( PDBnum_token.substr( 0, PDBnum_token.length() - 1 ).c_str() );
				icode = *PDBnum_token.rbegin();
			} else { // no insertion code
				PDBnum = atoi( PDBnum_token.c_str() );
			}
			++which_token;
			char chain;
			chain = get_token( which_token, tokens )[ 0 ];
			if (chain == '_') chain = ' ';
			++which_token;

			Size resid;
			resid = the_task.pdb_pose_map()->find( chain, PDBnum, icode );
			if (resid == 0){
				std::stringstream err_msg;
				err_msg  << "On line " << lineno << ", the pose does not have residue (" << chain << ", " << PDBnum << ").";
				onError( err_msg.str());
			}
			non_default_lines[ resid ] = true;

			while ( which_token <= ntokens ) {
				if ( comment_begin( tokens, which_token ) ) break; // ignore the rest of this line
				if ( command_map.find( get_token( which_token, tokens ) ) == command_map.end() ) {
					std::stringstream err_msg;
					err_msg  << "On line " << lineno << " command '" << get_token( which_token, tokens) <<"' is not recognized.";
					onError(err_msg.str());
					which_token++;
					continue;
				}

				ResfileCommandOP command = command_map[ get_token( which_token, tokens ) ]->clone();

				try{
					command->initialize_from_tokens( tokens, which_token, resid );
				} catch ( ResfileReaderException() ){
					// there was a problem with this command.  If we're doing error recovery skip to next command.
					while( which_token <= ntokens && command_map.find( get_token(which_token, tokens ) ) == command_map.end() )
						which_token++;
					continue;
				}
				commands_[ resid ].push_back( command );
			}

		} else { // the start token has not been read
			// read in default behaviors, store them, process them later
			if ( get_token( 1, tokens) == "START" ) {
				have_read_start_token = true;
			} else {
				for ( Size ii = 1; ii <= ntokens; ++ii ) {
					if ( comment_begin( tokens, ii ) ) break; // ignore the rest of this line
					default_tokens.push_back( get_token( ii, tokens ) );
					origin_lines_of_default_tokens.push_back( lineno );
				}
			}

		}
	}

	if ( ! have_read_start_token ) {
		T("core.pack.task.ResfileReader") << "RESFILE WARNING: reached the end of resfile without finding a 'start' token." << std::endl;
		T("core.pack.task.ResfileReader") << "RESFILE WARNING: No residue-specific behavior specified in resfile" << std::endl;
	}

	// now process default behaviors

	for ( Size ii = 1; ii <= non_default_lines.size(); ++ii ) {
		if ( ! non_default_lines[ ii ] ) {
			Size which_token = 1, ntokens = default_tokens.size();

			while( which_token <= ntokens ){
				if ( command_map.find( get_token( which_token, default_tokens ) ) == command_map.end() ) {
					std::stringstream err_msg;
					err_msg  << "The default  command '" << get_token( which_token, default_tokens) <<"' is not recognized.";
					onError(err_msg.str());
					which_token++;
					continue;
				}

				ResfileCommandOP command = command_map[ get_token( which_token, default_tokens ) ]->clone();

				try{
					command->initialize_from_tokens( default_tokens, which_token, 0 );
				} catch ( ResfileReaderException() ){
					// there was a problem with this command.  If we're doing error recovery skip to next command.
					while( which_token <= ntokens && command_map.find( get_token(which_token, default_tokens ) ) == command_map.end() )
						which_token++;
					continue;
				}
				default_commands_.push_back( command );
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
			T("core.pack.task.ResfileReader") << "Ignoring unknown one-letter amino acid code " << aas_to_keep[ ii ] << " while parsing PIKAA mode for residue " << resid << "." << std::endl;

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
			T("core.pack.task.ResfileReader")  << "Ignoring unknown one-letter nucleic acid code. " << *letter <<" while parsing PIKNA option for residue " << resid << ".";
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
		Error() << "RESFILE ERROR: PIKRNA must be followed by a string of allowed "
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
			Error() << "RESFILE ERROR: unknown one-letter nucleic acid code " << *letter
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
			T("core.pack.task.ResfileReader")  << "Ignoring Unknown one-letter amino acid code "<< aas_to_exclude[ ii ] << " while parsing NOTAA option for residue " << resid << ".";
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

	for ( ResidueLevelTask::ResidueTypeCAPListConstIter
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

	for ( ResidueLevelTask::ResidueTypeCAPListConstIter
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
	T("core.pack.task.ResfileReader") << "RESFILE NOTE: APOLA command deprecated.  Use APOLAR command instead.  Treating as APOLAR command." << std::endl;
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
///@brief NC allows a noncanonical residue; use one NC command per noncanonical
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
	if ( residue_set.has_name( nc_to_include_ )){
		task.nonconst_residue_task(resid).allow_noncanonical_aa( nc_to_include_ );
	}	else {
		std::stringstream err_msg;
		err_msg  << "Unable to add non-canonical amino acid " << nc_to_include_ << " because it is not in the ResidueTypeSet for residue " << resid << ".";
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
	command_map[ NATRO::name() ]        = new NATRO;
	command_map[ NATAA::name() ]        = new NATAA;
	command_map[ ALLAA::name() ]        = new ALLAA;
	command_map[ ALLAAxc::name() ]      = new ALLAAxc;
	command_map[ ALLAAwc::name() ]      = new ALLAAwc;
	command_map[ PIKAA::name() ]        = new PIKAA;
	command_map[ PIKNA::name() ]        = new PIKNA;
	command_map[ PIKRNA::name() ]        = new PIKRNA;
	command_map[ NOTAA::name() ]        = new NOTAA;
	command_map[ EMPTY::name() ]        = new EMPTY;
	command_map[ POLAR::name() ]        = new POLAR;
	command_map[ APOLAR::name() ]       = new APOLAR;
	command_map[ APOLA::name() ]        = new APOLA;
	command_map[ EX::name() ]           = new EX;
	command_map[ EX_CUTOFF::name() ]    = new EX_CUTOFF;
	command_map[ NC::name() ]           = new NC;
	command_map[ USE_INPUT_SC::name() ] = new USE_INPUT_SC;
	command_map[ AUTO::name() ]         = new AUTO;
	command_map[ SCAN::name() ]         = new SCAN;
	command_map[ TARGET::name() ]       = new TARGET;
	command_map[ NO_ADDUCTS::name() ]   = new NO_ADDUCTS;
	command_map[ FIX_HIS_TAUTOMER::name() ]   = new FIX_HIS_TAUTOMER;

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
parse_resfile( PackerTask & the_task, std::string filename)
{

	std::string resfile;
	utility::io::izstream file( filename );
	if (!file) {
		T("core.pack.task.ResfileReader") << "File:" << filename << " not found!\n";
		utility_exit_with_message( "Cannot open file " + filename );
	} else {
		//T() << "read file: " << filename << "\n";
	}
	utility::slurp( file, resfile );
	try{
		parse_resfile_string( the_task, resfile );
	} catch (ResfileReaderException) {
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
parse_resfile_string( PackerTask & the_task, std::string const & resfile_string ) throw(ResfileReaderException)
{
	using namespace std;
	istringstream resfile(resfile_string);
	ResfileContents contents( the_task, resfile );

	for ( Size ii = 1; ii <= the_task.total_residue(); ++ii ) {

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

std::string
get_token( const Size which_token, const utility::vector1<std::string> & tokens ) {

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
	std::transform( token.begin(), token.end(), token.begin(), (int(*)(int)) std::toupper);
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
	static bool interactive = basic::options::option[ basic::options::OptionKeys::run::interactive ].user();
	T("core.pack.task.ResfileReader") << message << std::endl;
	if (interactive){
		throw ResfileReaderException(message);
	} else {
		utility_exit();
	}
}


} //namespace task
} //namespace pack
} //namespace core
