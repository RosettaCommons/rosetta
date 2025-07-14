// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/LoopsFileIO.cc
/// @brief
/// @author Brian D. Weitzner

// Unit header
#include <protocols/loops/LoopsFileIO.hh>

// Package headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_definers/LoopsDefiner.hh>

// Project headers

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh> // needed to parse old style loop files
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>

#include <fstream>

#include <core/pose/ResidueIndexDescription.hh> // AUTO IWYU For ResidueIndexDescriptionPoseNum, RID_FileSource
#include <utility/json_spirit/json_spirit_reader.h> // AUTO IWYU For read

namespace protocols {
namespace loops {

using core::pose::ResidueIndexDescriptionCOP;

/// @details Auto-generated virtual destructor
LoopsFileData::~LoopsFileData() = default;

static basic::Tracer tr( "protocols.loops.LoopsFileIO" );

bool JSONFormattedLoopsFileReader::initialized_( false );
utility::vector1<std::string> JSONFormattedLoopsFileReader::valid_loop_file_keys_;

void
validate_loop_start_stop(
	bool prohibit_single_residue_loops,
	core::Size start,
	core::Size stop,
	std::string const & filename,
	core::Size linecount
)
{
	if ( start > stop || ( start == stop && prohibit_single_residue_loops ) ) {
		utility_exit_with_message(
			"[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) +
			" ): " + " Invalid loop definition (start residue " +
			( prohibit_single_residue_loops ? ">=" : ">" )  + " end residue) - ERROR"  );
	}
}


LoopFromFileData::LoopFromFileData() :
	skip_rate_( 0.0 ),
	extended_( false ),
	prohibit_single_residue_loops_( true )
{}

LoopFromFileData::LoopFromFileData(
	ResidueIndexDescriptionCOP start_res,
	ResidueIndexDescriptionCOP cutpoint_res,
	ResidueIndexDescriptionCOP end_res,
	core::Real skip_rate,
	bool extended,
	bool prohibit_single_residue_loops
) :
	start_res_( start_res ),
	cutpoint_res_( cutpoint_res ),
	end_res_( end_res ),
	skip_rate_( skip_rate ),
	extended_( extended ),
	prohibit_single_residue_loops_( prohibit_single_residue_loops )
{}

core::pose::ResidueIndexDescriptionCOP
make_rid_posenum( std::string const & fname, core::Size lineno, core::Size resnum ) {
	using namespace core::pose;
	RID_SourceCOP source( new RID_FileSource( fname, lineno ) );
	ResidueIndexDescriptionCOP rid( new ResidueIndexDescriptionPoseNum( source, resnum ) );
	return rid;
}

/// constructed the other way around (for the the PoseNumberedLoopReader)
LoopFromFileData::LoopFromFileData(
	SerializedLoop const & loop,
	std::string const & fname,
	bool prohibit_single_residue_loops
) :
	start_res_( make_rid_posenum( fname, 0, loop.start ) ),
	cutpoint_res_( make_rid_posenum( fname, 0, loop.cut ) ),
	end_res_( make_rid_posenum( fname, 0, loop.stop ) ),
	skip_rate_( loop.skip_rate ),
	extended_( loop.extended ),
	prohibit_single_residue_loops_( prohibit_single_residue_loops )
{}

/// @brief Convert from the information read from the loop-definition file into
/// residue indices that make sense in the context of this pose, and validate
/// the input loop definition (i.e. that the index of the start residue is
/// less than the index of the end residue ).
SerializedLoop
LoopFromFileData::resolve_as_serialized_loop_from_pose( core::pose::Pose const & pose ) const
{

	debug_assert( start_res_ );
	debug_assert( cutpoint_res_ );
	debug_assert( end_res_ );

	core::Size const start_res_index    = start_res_->resolve_index(    pose, true );
	core::Size const cutpoint_res_index = cutpoint_res_->resolve_index( pose, true );
	core::Size const end_res_index      = end_res_->resolve_index(      pose, true );

	std::string fname("");
	core::Size lineno(0);
	auto rid_source( utility::pointer::dynamic_pointer_cast< core::pose::RID_FileSource const >( start_res_->get_source() ) );
	if ( rid_source != nullptr ) {
		fname = rid_source->fname();
		lineno = rid_source->linenum();
	}

	validate_loop_start_stop(
		prohibit_single_residue_loops_, start_res_index, end_res_index,
		fname, lineno );

	SerializedLoop loop;
	loop.start = start_res_index;
	loop.cut = cutpoint_res_index;
	loop.stop = end_res_index;
	loop.skip_rate = skip_rate_;
	loop.extended = extended_;

	return loop;
}

LoopsOP
LoopsFileData::resolve_loops(
	core::pose::Pose const & pose
) const
{
	SerializedLoopList sloops = resolve_as_serialized_loops( pose );
	return utility::pointer::make_shared< Loops >( sloops );
}

SerializedLoopList
LoopsFileData::resolve_as_serialized_loops(
	core::pose::Pose const & pose
) const
{
	SerializedLoopList sloops; sloops.reserve( size() );
	for ( core::Size ii = 1; ii <= size(); ++ii ) {
		sloops.push_back( (*this)[ ii ].resolve_as_serialized_loop_from_pose( pose ) );
	}
	return sloops;
}

core::Size LoopsFileData::size() const
{
	return loops_file_data_.size();
}

void LoopsFileData::resize( core::Size new_size )
{
	loops_file_data_.resize( new_size );
}

void LoopsFileData::push_back( LoopFromFileData const & loop )
{
	loops_file_data_.push_back( loop );
}

void LoopsFileData::insert_loop_at_index( LoopFromFileData const & loop, core::Size i )
{
	loops_file_data_[ i ] = loop;
}

LoopFromFileData const & LoopsFileData::operator [] ( core::Size const i ) const
{
	return loops_file_data_[ i ];
}


/// @details default constructor: set state to not expect a pose, but to return an empty Loops object
GuardedLoopsFromFile::GuardedLoopsFromFile() :
	in_charge_( true ),
	pose_has_resolved_loop_indices_( true ),
	rely_on_loopfile_indices_( false ),
	loops_( utility::pointer::make_shared< Loops >() )
{}

/// @details construct from LoopsFileData: set state to expect a Pose
GuardedLoopsFromFile::GuardedLoopsFromFile( LoopsFileData const & lfd ) :
	in_charge_( true ),
	pose_has_resolved_loop_indices_( false ),
	rely_on_loopfile_indices_( true ),
	loops_file_data_( utility::pointer::make_shared< LoopsFileData >(lfd) ),
	loops_( utility::pointer::make_shared< Loops >() )
{}

/// @details construct from LoopsFileData: set state to expect a Pose
GuardedLoopsFromFile::GuardedLoopsFromFile( loops_definers::LoopsDefinerOP ld ) :
	in_charge_( true ),
	pose_has_resolved_loop_indices_( false ),
	rely_on_loopfile_indices_( true ),
	loops_definer_( ld ),
	loops_( utility::pointer::make_shared< Loops >() )
{}

/// @details construct from a LoopsOP
GuardedLoopsFromFile::GuardedLoopsFromFile( LoopsOP loops ) :
	in_charge_( false ),
	pose_has_resolved_loop_indices_( true ),
	rely_on_loopfile_indices_( false ),
	loops_(std::move( loops )) // copy the pointer -- assume this pointer comes from some other GuardedLoopsFromFile object and that I am not in charge
{}


/// @details construct from a LoopsOP
GuardedLoopsFromFile::GuardedLoopsFromFile( Loops const & loops ) :
	in_charge_( false ),
	pose_has_resolved_loop_indices_( true ),
	rely_on_loopfile_indices_( false ),
	loops_( utility::pointer::make_shared< Loops >( loops ) ) // copy the contents into a new loops object -- assume this pointer comes from some other GuardedLoopsFromFile object and that I am not in charge
{}

/// @details Shallow copy of the LoopsOP data so that it can be shared between
/// multiple objects; take the "in_charge_" bit from the source in the event that
/// this is a called as part of a clone.
GuardedLoopsFromFile::GuardedLoopsFromFile( GuardedLoopsFromFile const & ) = default;

/// @details Shallow copy of the LoopsOP data so that it can be shared between
/// multiple objects -- also assume that the GuardedLoopsFromFile object is the one that is in charge.
GuardedLoopsFromFile::GuardedLoopsFromFile( GuardedLoopsFromFile const & src, bool ) :
	GuardedLoopsFromFile(src) // Delegate copy for most items
{
	in_charge_ = false;
}

GuardedLoopsFromFile::~GuardedLoopsFromFile() = default;

GuardedLoopsFromFile &
GuardedLoopsFromFile::operator = ( GuardedLoopsFromFile const & ) = default;

void
GuardedLoopsFromFile::in_charge( bool setting )
{
	in_charge_ = setting;
}

bool GuardedLoopsFromFile::in_charge() const { return in_charge_; }

/// @details This will re-resolve the loop indices with the input pose, even if they had been resolved in the past.
void
GuardedLoopsFromFile::resolve_loop_indices( core::pose::Pose const & pose )
{
	if ( ! in_charge_ ) return;
	if ( ! rely_on_loopfile_indices_ ) { pose_has_resolved_loop_indices_ = true; return; }
	if ( loops_file_data_ != nullptr ) {
		(*loops_) = *( loops_file_data_->resolve_loops( pose ));
	} else if ( loops_definer_ != nullptr ) {
		Loops loops( loops_definer_->apply(pose) );
		(*loops_) = loops;
	} else {
		utility_exit_with_message("GuardedLoopsFromFile was supposed to resolve loops, but it has nothing to resolve.");
	}
	pose_has_resolved_loop_indices_ = true;
}

/// @details This will only resolve the loop indices once, so repeated calls to this will not alter the Loops
/// object with (possibly) new loop indices.
void
GuardedLoopsFromFile::resolve_loop_indices_once( core::pose::Pose const & pose )
{
	if ( ! in_charge_ ) return;
	if ( pose_has_resolved_loop_indices_ ) return;
	resolve_loop_indices( pose );
}

/// @brief request the LoopsCOP pointer; asserts that the loop indices
/// have been resolved or that "I am not in charge".
LoopsCOP
GuardedLoopsFromFile::loops() const
{
	debug_assert( !in_charge_ || pose_has_resolved_loop_indices_ );
	return loops_;
}

/// @brief request the LoopsOP pointer; asserts that the loop indices
/// have been resolved or that "I am not in charge".
LoopsOP
GuardedLoopsFromFile::loops()
{
	debug_assert( !in_charge_ || pose_has_resolved_loop_indices_ );
	return loops_;
}

/// @brief Deep copy into the loops data; this updates the single Loops object that is / can be shared
/// among multiple objects.
void
GuardedLoopsFromFile::loops( Loops const & setting )
{
	pose_has_resolved_loop_indices_ = true; // the assumption is that the user has resolved the indices themselves.
	rely_on_loopfile_indices_ = false; // the user has indicated the Loops data should be used and not the LoopsFileData
	*loops_ = setting;

}

/// @brief Shallow copy of the loops data
void
GuardedLoopsFromFile::set_loops_pointer( LoopsOP setting )
{
	debug_assert( !in_charge_ );
	pose_has_resolved_loop_indices_ = true; // the assumption is that the user has resolved the indices themselves.
	rely_on_loopfile_indices_ = false; // the user has indicated the Loops data should be used and not the LoopsFileData
	loops_ = setting;

}

/// @brief set the LoopsFileData object directly
void
GuardedLoopsFromFile::loops( LoopsFileData const & setting )
{
	debug_assert( in_charge_ ); // if not in_charge, then this data will never be used.
	pose_has_resolved_loop_indices_ = false;
	rely_on_loopfile_indices_ = true;
	loops_file_data_ = utility::pointer::make_shared< LoopsFileData >(setting);
}

/// @brief read access to the LoopsFileData
LoopsFileDataCOP
GuardedLoopsFromFile::loops_file_data() const
{
	return loops_file_data_;
}

/// @brief set the LoopsFileData object directly
void
GuardedLoopsFromFile::loops( loops_definers::LoopsDefinerOP const & setting )
{
	debug_assert( in_charge_ ); // if not in_charge, then this data will never be used.
	pose_has_resolved_loop_indices_ = false;
	rely_on_loopfile_indices_ = true;
	loops_definer_ = setting;
}

/// @brief read access to the LoopsFileData
loops_definers::LoopsDefinerCOP
GuardedLoopsFromFile::loops_definer() const
{
	return loops_definer_;
}


LoopsFileIO::LoopsFileIO() : utility::VirtualBase()
{

}

LoopsFileIO::LoopsFileIO( const LoopsFileIO & ) : utility::VirtualBase()
{
}

// destructor
LoopsFileIO::~LoopsFileIO()= default;

//////////////////////////////////////////////////////////////////////
std::ostream & operator<< ( std::ostream & os, const LoopsFileIO & /*loops*/ ) {
	/*
	os << "LOOP  begin  end  cut  skip_rate  extended" << std::endl;
	for ( Loops::const_iterator it = loops.begin(), it_end = loops.end();
	it != it_end; ++it ) {
	os << *it << std::endl;
	}
	*/
	return os;
}


LoopsFileDataOP
LoopsFileIO::read_loop_file(
	std::string const & filename,
	bool prohibit_single_residue_loops
)
{
	std::ifstream infile( filename.c_str() );

	if ( !infile.good() ) {
		utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + filename + "'" );
	}
	return read_loop_file_stream( infile, filename, prohibit_single_residue_loops );
}

LoopsFileDataOP
LoopsFileIO::read_loop_file_stream(
	std::istream & loopfstream,
	std::string const & filename,
	bool prohibit_single_residue_loops
)
{
	// if the first line is a comment, consume this line and determine if the line is
	// specifying the new file format
	// "# FORMAT JSON" is the specification for the new format.
	// There may or may not be a newline immediately after the "JSON" token.
	core::Size lines_pre_read( 0 );
	if ( loopfstream.peek() == '#' ) {
		std::string line;
		getline( loopfstream, line );
		++lines_pre_read;
		utility::vector1< std::string > tokens ( utility::split( line ) );

		if ( tokens.size() >= 3 && tokens[2] == "FORMAT" ) {
			if ( tokens[3] == "JSON" || tokens[3] == "JSON\n" || tokens[3] == "JSON\r" ) {
				JSONFormattedLoopsFileReader reader;
				reader.set_linecount_offset( lines_pre_read );
				LoopsFileDataOP loops = reader.read_loop_file( loopfstream, filename, prohibit_single_residue_loops );
				return loops;
			}
		}

	}

	PoseNumberedLoopFileReader reader;
	reader.set_linecount_offset( lines_pre_read );
	SerializedLoopList sloops = reader.read_pose_numbered_loops_file( loopfstream, filename, prohibit_single_residue_loops );
	LoopsFileDataOP loops( new LoopsFileData );
	loops->resize( sloops.size() );
	for ( core::Size ii = 1; ii <= sloops.size(); ++ii ) {
		loops->insert_loop_at_index( LoopFromFileData( sloops[ ii ], filename, prohibit_single_residue_loops ), ii );
	}
	return loops;
}

PoseNumberedLoopFileReader::PoseNumberedLoopFileReader() :
	loop_line_begin_token_( "LOOP" ),
	linecount_offset_( 0 )
{}

SerializedLoopList
PoseNumberedLoopFileReader::read_pose_numbered_loops_file(
	std::istream & is,
	std::string const & filename,
	bool strict_looprelax_checks
)
{
	std::string line;
	core::Size linecount = linecount_offset_;
	int errcount=50; //if we reach 0 we bail!

	SerializedLoopList loops;

	while ( getline( is, line) ) {
		linecount++;
		utility::vector1< std::string > tokens ( utility::split( line ) );

		SerializedLoop current_loop;
		if ( tokens.size() > 0 ) {
			if ( tokens[1].substr(0,3) == "END" ) break;
			if ( tokens[1] == loop_line_begin_token_ ) {
				if ( tokens.size() < 3 ) {
					utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Minimum of 3 tokens necessary (begin, end, cutpoint)"  );
				}
				if ( tokens.size() > 6 ) {
					utility_exit_with_message( "[ERROR] Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + " Maximum of 6 tokens allowed (LOOP begin end cutpoint skiprate extended)"  );
				}
				current_loop.start = (core::Size) atoi(tokens[2].c_str());
				current_loop.stop = (core::Size) atoi(tokens[3].c_str());
				current_loop.cut = 0;        // default - let LoopRebuild choose cutpoint
				current_loop.skip_rate = 0.0;  // default - never skip
				//std::string extend_loop_str;
				bool extend_loop = false;

				if ( tokens.size() > 3 ) {
					current_loop.cut = (core::Size) atoi(tokens[4].c_str());
				}
				if ( tokens.size() > 4 ) {
					current_loop.skip_rate = atof(tokens[5].c_str());
				}
				if ( tokens.size() > 5 ) {
					if ( tokens[6] == "X" ) {
						tr.Error << "Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "[WARNING] DEPRECATED old style extended marker X is used" << std::endl;
						extend_loop = true;
						if ( errcount > 0 ) errcount--;
						else {
							utility_exit_with_message( "too many errors in loop-file " + filename );
						}
					} else {
						int extended_token = atoi(tokens[6].c_str());
						if  ( extended_token == 0 ) extend_loop = false;
						else                        extend_loop = true;
					}
				}

				current_loop.extended = extend_loop;
				validate_loop_start_stop( strict_looprelax_checks, current_loop.start, current_loop.stop, filename, linecount );
				loops.push_back( current_loop );
			} else if ( tokens[1][0] != '#' ) {
				if ( tokens.size() >= 2 ) {
					tr.Error << "Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "DEPRECATED r++ style loopfile" << std::endl;

					if ( errcount>0 ) {
						errcount--;
					} else {
						utility_exit_with_message( "too many errors in loop-file " + filename );
					}

					current_loop.start = (core::Size) atoi(tokens[1].c_str());
					current_loop.stop   = (core::Size) atoi(tokens[2].c_str());
					current_loop.cut = 0;        // default - let LoopRebuild choose cutpoint
					current_loop.skip_rate = 0.0;  // default - never skip
					bool extend_loop = false;    // default - not extended
					if ( tokens.size() > 2 ) {
						current_loop.cut = (core::Size) atoi(tokens[3].c_str());
					}
					if ( tokens.size() > 3 ) {
						current_loop.skip_rate = atof(tokens[4].c_str());
					}
					if ( tokens.size() > 4 ) {
						if ( tokens[5] == "X" ) {
							tr.Error << "Error parsing " + filename + " ( line " + ObjexxFCL::string_of( linecount ) + " ): " + "[WARNING] DEPRECATED old style extended marker X is used" << std::endl;
							extend_loop = true;
						} else {
							int extended_token = atoi(tokens[5].c_str());
							if ( extended_token == 0 ) extend_loop = false;
							else                extend_loop = true;
						}
					}
					current_loop.extended = extend_loop;

					validate_loop_start_stop( strict_looprelax_checks, current_loop.start, current_loop.stop, filename, linecount );
					loops.push_back( current_loop );

				} else {
					tr.Warning << "Skipping line '" << line << "'" << std::endl;
				}
			}
		}
	} //while

	return loops;
}

void PoseNumberedLoopFileReader::set_linecount_offset( core::Size setting ) {
	linecount_offset_ = setting;
}

/// @brief For code that relys on reading loop-file-formatted ranges if residues
/// but which really ought to use
void
PoseNumberedLoopFileReader::hijack_loop_reading_code_set_loop_line_begin_token(
	std::string const & token
)
{
	loop_line_begin_token_ = token;
}


LoopsFileDataOP
JSONFormattedLoopsFileReader::read_loop_file(
	std::istream & is,
	std::string const & filename,
	bool prohibit_single_residue_loops
){
	utility::json_spirit::mValue mapped_json;

	// sanity check
	if ( utility::json_spirit::read(is, mapped_json) ) utility_exit_with_message( "Unable to read loops file '" + filename + "'." );
	if ( mapped_json.type() != utility::json_spirit::obj_type ) utility_exit_with_message( "The loop file '" + filename + "' is not formatted correctly.  Please see the documenation.");

	return parse_json_formatted_data( mapped_json, prohibit_single_residue_loops, filename );
}

LoopsFileDataOP
JSONFormattedLoopsFileReader::parse_json_formatted_data(
	utility::json_spirit::mValue & json_data,
	bool prohibit_single_residue_loops,
	std::string const & filename
){
	std::string LoopSetKey = "LoopSet";
	if ( ! json_data.get_obj().count( LoopSetKey ) ) utility_exit_with_message( "The loop file '" + filename + "' is not formatted correctly.  Please see the documenation.");
	utility::json_spirit::mArray & array = json_data.get_obj()[LoopSetKey].get_array();
	if ( ! array.size() ) utility_exit_with_message( "The LoopList appears to be empty.  Please check your input file, '" + filename + "'." );

	LoopsFileDataOP loops( new LoopsFileData );
	core::Size count_lines_approximate = linecount_offset_;
	for ( auto & i : array ) {
		ensure_all_fields_are_valid( i, filename );
		LoopFromFileData current_loop = LoopFromFileData();

		current_loop.start_res(    parse_json_residue_info( i, start,     filename, count_lines_approximate ) );
		current_loop.end_res(      parse_json_residue_info( i, stop,      filename, count_lines_approximate ) );
		current_loop.cutpoint_res( parse_json_residue_info( i, cut_point, filename, count_lines_approximate ) );
		current_loop.prohibit_single_residue_loops( prohibit_single_residue_loops );

		parse_configuration_options( i, current_loop );
		loops->push_back( current_loop );
	}
	return loops;
}

void
JSONFormattedLoopsFileReader::ensure_all_fields_are_valid( utility::json_spirit::mValue & json_data, std::string const & fname )
{
	if ( json_data.type() != utility::json_spirit::obj_type ) {
		utility_exit_with_message( "The loop file '" + fname + "' is not formatted correctly.  Please see the documenation.");
	}

	setup_residue_type_map();
	for ( auto it = json_data.get_obj().begin(); it != json_data.get_obj().end(); ++it ) {

		if ( std::find( valid_loop_file_keys_.begin(), valid_loop_file_keys_.end(), it->first ) == valid_loop_file_keys_.end() ) {
			utility_exit_with_message( "Unknown key \"" + it->first + ".\"  Please check your input file." );
		}
	}
}

void JSONFormattedLoopsFileReader::set_linecount_offset( core::Size setting )
{
	linecount_offset_ = setting;
}


ResidueIndexDescriptionCOP
JSONFormattedLoopsFileReader::parse_json_residue_info(
	utility::json_spirit::mValue & json_loop_data,
	ResidueIdentifier residue_identifier,
	std::string const & filename,
	core::Size & approximate_linenumber
) {
	core::Size resNo = 0;
	std::string chain_id;
	char insert_code = ' ';
	bool usesPDBNumbering = false;

	std::string res_identity =  name_from_residue_identifier( residue_identifier );
	if ( json_loop_data.get_obj().count( res_identity ) ) {
		++approximate_linenumber; // we found the token, ergo, increment the approximate linenumber

		ensure_all_fields_are_valid( json_loop_data.get_obj()[ res_identity ], filename );
		utility::json_spirit::mObject json_representation_of_residue_data = json_loop_data.get_obj()[ res_identity ].get_obj();

		// resSeq - this is the residue number that will be used.
		std::string residue_number = name_from_loop_configuration( resSeq );
		if ( json_representation_of_residue_data.count( residue_number ) ) {
			if ( json_representation_of_residue_data[ residue_number ].type() != utility::json_spirit::int_type ) {
				utility_exit_with_message( "The \"resSeq\" field must be an integer.  Please check your input file, '" + filename + "'." );
			}
			resNo = json_representation_of_residue_data[ residue_number ].get_int();
		} else {
			utility_exit_with_message( "The \"resSeq\" field for " + res_identity + " must be defined.  Please check your input file, '" + filename + "'." );
		}

		// chainID - One letter chainID for the residue number.  This is necessary when using PDB numbering.
		std::string chain_identifier = name_from_loop_configuration( chainID );
		if ( json_representation_of_residue_data.count( chain_identifier ) ) {
			if ( json_representation_of_residue_data[ chain_identifier ].type() != utility::json_spirit::str_type ) {
				utility_exit_with_message( "The \"chainID\" field must be a one character string.  Please check your input file, '" +filename + "'." );
			}
			chain_id = json_representation_of_residue_data[ chain_identifier ].get_str();
			usesPDBNumbering = true;
		}

		// iCode - insertion codes are sometimes used in PDBs.  This allows any PDB residue to be used to define the loop.
		std::string insertion_code = name_from_loop_configuration( iCode );
		if ( json_representation_of_residue_data.count( insertion_code ) ) {
			if ( ! usesPDBNumbering ) {
				utility_exit_with_message( "Using an insertion code requires specifying the residue's chainID.  Please check your input file, '" + filename + "'." );
			}
			if ( json_representation_of_residue_data[ insertion_code ].type() != utility::json_spirit::str_type ) {
				utility_exit_with_message( "The \"iCode\" field must be a one character string or omitted.  Please check your input file, '" + filename + "'." );
			}
			std::string tmp_iCode = json_representation_of_residue_data[ insertion_code ].get_str();
			if ( tmp_iCode.length() != 1 ) {
				utility_exit_with_message( "Insertion codes must be exactly one character long.  Please check your input file, '" + filename + "'." );
			}
			insert_code = char( tmp_iCode[0] );
		}
	} else if ( res_identity != name_from_residue_identifier( cut_point ) ) {
		utility_exit_with_message( "The \"" + res_identity + "\" residue must be specified.  Please check your input file '" + filename + "'." );
	}

	core::pose::RID_SourceCOP rid_source( new core::pose::RID_FileSource( filename, approximate_linenumber ) );
	if ( usesPDBNumbering ) {
		return utility::pointer::make_shared< core::pose::ResidueIndexDescriptionPDB >( rid_source, chain_id, resNo, insert_code );
	} else {
		return utility::pointer::make_shared< core::pose::ResidueIndexDescriptionPoseNum >( rid_source, resNo );
	}
}

void
JSONFormattedLoopsFileReader::parse_configuration_options(
	utility::json_spirit::mValue & json_loop_data,
	LoopFromFileData & loop
){
	if ( json_loop_data.get_obj().count( name_from_loop_configuration( extras ) ) ) {

		utility::json_spirit::mObject json_representation_of_residue_data = json_loop_data.get_obj()[ name_from_loop_configuration( extras ) ].get_obj();

		// skip rate
		std::string skip_rate_name = name_from_loop_configuration( skip_rate );
		if ( json_representation_of_residue_data.count( skip_rate_name ) ) {
			if ( json_representation_of_residue_data[ skip_rate_name ].type() != utility::json_spirit::real_type && json_representation_of_residue_data[ skip_rate_name ].type() != utility::json_spirit::int_type ) {
				utility_exit_with_message( "Skip rates in loop files must be a floating point number in [0, 1).  Please check your input file." );
			}
			loop.skip_rate( json_representation_of_residue_data[ skip_rate_name ].get_real() );

			if ( !(loop.skip_rate() >= 0.0 && loop.skip_rate() < 1.0) ) {
				utility_exit_with_message( "Skip rates in loop files must be in [0, 1).  Please check your input file." );
			}
		}

		// extend
		std::string extend_name = name_from_loop_configuration( extend );
		if ( json_representation_of_residue_data.count( extend_name ) ) {
			if ( json_representation_of_residue_data[ extend_name ].type() != utility::json_spirit::bool_type ) {
				utility_exit_with_message( "The \"extend\" field must be a boolean.  Use the string literals \"true\" or \"false\" to toggle the behavior.  If this field is omitted, it will default to \"false.\" Please check your input file." );
			}
			loop.extended( json_representation_of_residue_data[ extend_name ].get_bool() );
		}

		// use pose numbering
		std::string pose_numbering_name = name_from_loop_configuration( use_pose_numbering );
		if ( json_representation_of_residue_data.count( pose_numbering_name ) ) {
			if ( json_representation_of_residue_data[ pose_numbering_name ].type() != utility::json_spirit::bool_type ) {
				utility_exit_with_message( "The \"use_pose_numbering\" field must be a boolean.  Use the string literals \"true\" or \"false\" to toggle the behavior.  If this field is omitted, it will default to \"false.\" Please check your input file.  Note that this field exists only to provide compatibility with the old style loops files." );
			}
			if ( json_representation_of_residue_data[ pose_numbering_name ].get_bool() ) {

				/// Talk to Andrew about this stuff.
				/*if ( ! loop.start_res().pose_numbered() || loop.start_res().pose_numbered() != loop.end_res().pose_numbered() || loop.start_res().pose_numbered() != loop.cutpoint_res().pose_numbered() ) {
				utility_exit_with_message( "Unable to reconcile the mixture of Rosetta and PDB numbering.  Please check your input file to ensure the correct numbering scheme is being used." );
				}
				return;
				*/
			}
		}
	}

	/// Talk to Andrew about this stuff.


	// make sure the numbering scheme is consistent
	/*if ( loop.start_res().pose_numbered() || loop.start_res().pose_numbered() != loop.end_res().pose_numbered() || loop.start_res().pose_numbered() != loop.cutpoint_res().pose_numbered() ) ) {
	utility_exit_with_message( "Unable to reconcile the mixture of Rosetta and PDB numbering.  Please check your input file to ensure the correct numbering scheme is being used." );
	} */
	return;
}


void JSONFormattedLoopsFileReader::setup_residue_type_map()
{
	if ( initialized_ ) return;
	initialized_ = true;

	valid_loop_file_keys_.resize( number_of_configuration_keywords );
	valid_loop_file_keys_[start]              = "start";
	valid_loop_file_keys_[stop]               = "stop";
	valid_loop_file_keys_[cut_point]          = "cut";
	valid_loop_file_keys_[extras]             = "extras";
	valid_loop_file_keys_[resSeq]             = "resSeq";
	valid_loop_file_keys_[iCode]              = "iCode";
	valid_loop_file_keys_[chainID]            = "chainID";
	valid_loop_file_keys_[skip_rate]          = "skip_rate";
	valid_loop_file_keys_[extend]             = "extend";
	valid_loop_file_keys_[use_pose_numbering] = "use_pose_numbering";
}

std::string JSONFormattedLoopsFileReader::name_from_residue_identifier( ResidueIdentifier residue_identifier )
{
	setup_residue_type_map();
	return valid_loop_file_keys_[ residue_identifier ];
}

std::string JSONFormattedLoopsFileReader::name_from_loop_configuration( LoopConfiguration loop_configuration )
{
	setup_residue_type_map();
	return valid_loop_file_keys_[ loop_configuration ];
}

} // namespace loops
} // namespace protocols
