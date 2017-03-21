// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/SilentFileData.cc
///
/// @brief silent input file reader for mini.
/// @author James Thompson

// C++ Headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <set>
#include <sstream>

// mini headers
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/chemical/util.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>

#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/stream_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/prof.hh>

#include <numeric/random/random.hh>
#include <numeric/random/reservoir_sample.hh>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {

static THREAD_LOCAL basic::Tracer tr( "core.io.silent.SilentFileData" );

//SilentFileData::SilentFileData() :
// SilentFileData( SilentFileOptions() )
//{}

SilentFileData::SilentFileData( SilentFileOptions const & options ) :
	filename_(),
	store_argv_in_file_( false ),
	strict_column_mode_( false ),
	record_source_( false ),
	silent_struct_type_(""), // by default its option controlled.
	verbose_( true ),
	options_( options )
{}

//SilentFileData:: SilentFileData( std::string const& filename ) :
// SilentFileData( filename, SilentFileOptions() )
//{}

SilentFileData::SilentFileData( std::string const& filename, SilentFileOptions const & options ) :
	filename_( filename ),
	store_argv_in_file_( false ),
	strict_column_mode_( false ),
	record_source_( false ),
	silent_struct_type_(""), // by default its option controlled.
	verbose_( true ),
	options_( options )
{}

//SilentFileData::SilentFileData(
// const std::string &filename,
// bool  store_argv_in_file,
// bool  strict_column_mode,
// const std::string & silent_struct_type
//) :
// SilentFileData( filename, store_argv_in_file, strict_column_mode, silent_struct_type, SilentFileOptions() )
//{}

SilentFileData::SilentFileData(
	const std::string &filename,
	bool  store_argv_in_file,
	bool  strict_column_mode,
	const std::string & silent_struct_type,
	SilentFileOptions const & options
) :
	filename_( filename ),
	store_argv_in_file_( store_argv_in_file ),
	strict_column_mode_( strict_column_mode ),
	record_source_( false ),
	silent_struct_type_( silent_struct_type ),
	verbose_( true ),
	options_( options )
{}


SilentFileData::~SilentFileData() {
	clear_structure_map();
}

utility::vector1< std::string >
SilentFileData::read_tags_fast( std::string const & filename ) const {
	utility::vector1< std::string > tags_in_file;
	read_tags_fast( filename, tags_in_file );
	return tags_in_file;
} // read_tags_fast

utility::vector1< std::string > SilentFileData::tags() const {
	utility::vector1< std::string > tag_list;

	for ( const_iterator it=structure_map_.begin(),
			it_end = structure_map_.end(); it != it_end; ++it ) {
		tag_list.push_back( it->decoy_tag() );
	}
	return tag_list;
} // tags


std::string
SilentFileData::get_sequence( std::string const & filename )
{
	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message(
			"ERROR: Unable to open silent_input file: '" + filename + "'"
		);
	}

	std::string sequence_line;
	getline( data, sequence_line ); // sequence line

	if ( sequence_line.substr(0,9) != "SEQUENCE:" ) {
		tr.Error << "bad format in first line of silent file " << filename << " (function read_stream):"
			<< std::endl;
		tr.Error << sequence_line << std::endl;
		return "";
	}
	std::string sequence = sequence_line.substr( 9 );

	// might be better to continue reading into FULL_MODEL_PARAMETERS line, if it exists.
	// That would have more detailed information on the modeled sequence -- rhiju

	return ObjexxFCL::strip_whitespace( sequence );
}


std::pair< utility::vector1< int >, utility::vector1< char > >
SilentFileData::get_resnum( std::string const & filename )
{
	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message(
			"ERROR: Unable to open silent_input file: '" + filename + "'"
		);
	}

	utility::vector1< int > resnum;
	utility::vector1< char > chain;

	std::string line;
	getline( data, line ); // sequence line
	getline( data, line ); // score line
	while ( getline(data,line) ) {
		if ( line.substr(0,8) == "RES_NUM "  ) {
			std::stringstream line_stream( line );
			figure_out_residue_numbers_from_line( line_stream, resnum, chain );
			break;
		}
	}
	if ( resnum.size() == 0 ) {
		std::string sequence  = get_sequence( filename );
		for ( Size n = 1; n <= sequence.size(); n++ ) {
			resnum.push_back( n ); chain.push_back( 'A' );
		}
	}
	return std::make_pair( resnum, chain );
}

bool SilentFileData::read_tags_fast(
	std::string const & filename,
	utility::vector1< std::string > & tags_in_file
) const {

	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message(
			"ERROR: Unable to open silent_input file: '" + filename + "'"
		);
	}

	std::string line;
	getline( data, line ); // sequence line
	getline( data, line ); // score line

	while ( getline(data,line) ) {
		if ( line.substr(0,7) == "SCORE: " && line.substr(8,20).find( "score" ) == std::string::npos ) {
			std::istringstream l( line );

			std::string current_word;
			while ( !l.fail() ) {
				l >> current_word;
			}
			tags_in_file.push_back( current_word );
		}
	} // while( getline(data,line) )
	return true;
}

bool
SilentFileData::matched_tags(
	std::string const & expression,
	std::string const & mode,
	utility::vector1< std::string >& tags_in_file//,
	// utility::vector1< SilentStructOP >& decoys_in_file,
	// bool ignore_decoys
) const {

	//open file
	utility::io::izstream data( filename().c_str() );
	if ( !data.good() ) {
		//   utility_exit_with_message(
		//    "ERROR: Unable to open silent_input file: '" + filename() + "'"
		//   );
		tr.Info << "file: " << filename() << " not found" << std::endl;
		return false;
	}

	//figure out search mode
	bool return_first=false;
	bool return_last=false;
	if ( mode == "first" ) {
		return_first=true;
	} else if ( mode == "last" ) {
		return_last=true;
	} else if ( mode != "all" ) {
		utility_exit_with_message( "illegal mode selection in SilentFileData::matched_tag(). available modes are first, last, and all " );
	}

	//search the tags
	std::string line;
	std::string current_tag;
	std::string final_tag;
	//bool ignore_this_decoy = true;
	utility::vector1< std::string > the_lines;
	while ( getline(data,line) ) {
		//find a line with SCORE: but without "score"
		if ( line.substr(0,7) == "SCORE: " && line.substr(8,20).find( "score" ) == std::string::npos ) {
			tr.Info << "reading the lines in file:\n" << line << std::endl;
			//find last word in line -> current_tag
			std::istringstream l( line );
			while ( !l.fail() ) {
				l >> current_tag;
			}
			tr.Info << "current_tag:" << current_tag << std::endl;
			//does it match ?
			if ( current_tag.find(expression)!=std::string::npos ) {
				//    ignore_this_decoy = ignore_decoys;
				if ( !return_last ) {
					tags_in_file.push_back( current_tag );
				} else final_tag = current_tag;
				if ( return_first ) return true;  //if we want decoys this will be wrong
			}
		} // SCORE: header of decoy
		//  if ( !ignore_this_decoy ) {
		//   utility_exit_with_message( "SilentFileData::matched_tag() function stubbed out -- cannot return decoys yet" );
		//  }
	} // while( getline(data,line) )
	//  if ( return_last ) {
	//   if ( current_tag.find(expression)!=std::string::npos ) {
	//    tags_in_file.push_back( current_tag );
	//    return true;
	//   }
	//  }
	if ( return_last ) {
		if ( final_tag.find( expression ) != std::string::npos ) {
			tags_in_file.push_back( final_tag );
			return true;
		}
	}
	return tags_in_file.size();
}


//////////////////////////////////////////////////////
void SilentFileData::add_structure_replace_tag_if_necessary(
	SilentStructOP & new_struct
) {

	std::string old_tag = new_struct->decoy_tag();
	std::string new_tag = old_tag;

	using ObjexxFCL::string_of;

	int n_tries = 0;
	while ( has_tag( new_tag ) ) {
		// assign a new tag. Stupid algorithm for the moment, there's nothing to
		// keep the length of tags from exploding into long strings. If these are
		// getting to be too long, use the renumber_all_decoys() method to assign
		// more sensible tags.
		++n_tries;
		new_tag = old_tag + "_" + string_of( n_tries );
	}

	if ( n_tries > 0 ) {
		tr.Warning << "renamed tag " << old_tag << " to " << new_tag;
		tr.Warning << "    (SilentStruct with " << old_tag;
		tr.Warning << " already exists!)" << std::endl;
	}
	new_struct->decoy_tag( new_tag );

	add_structure( new_struct );
}

//////////////////////////////////////////////////////
void SilentFileData::add_structure(
	SilentStructOP const & new_struct
) {

	if ( add_as_other_struct_if_relevant( new_struct ) ) return;

	std::string const & new_tag( new_struct->decoy_tag() );

	if ( !has_tag( new_tag ) ) {
		structure_map_[ new_tag ] = new_struct;
		structure_list_.push_back(new_struct);
		tr.Debug << "added structure with tag " << new_tag << std::endl;
		tr.Debug << "now have " << structure_map_.size() << " structures."
			<< std::endl;
		tr.flush();
	} else {
		SilentStructOP new_struct_clone = new_struct->clone();
		add_structure_replace_tag_if_necessary( new_struct_clone );
		structure_list_.push_back(new_struct_clone);
	}

} // add_structure

void SilentFileData::add_structure(
	SilentStruct const & new_struct
) {
	add_structure( new_struct.clone() );
} // add_structure

/// @details stepwise modeling can involve holding references to 'other poses' which, taken together
///    correspond to the full model of the pose.
bool
SilentFileData::add_as_other_struct_if_relevant( SilentStructOP const & new_struct ) {

	std::string const & new_tag( new_struct->decoy_tag() );

	std::map< std::string, std::string > const & comments = new_struct->get_all_comments();
	if ( comments.find( "OTHER_POSE" ) == comments.end() ) return false;
	runtime_assert( new_struct->scoreline_prefix() == "OTHER:" );
	runtime_assert( has_tag( new_tag ) );
	SilentStructOP parent_silent_struct = structure_map_.find( new_tag )->second;
	parent_silent_struct->add_other_struct( new_struct );
	return true;
}


void SilentFileData::renumber_all_decoys() {
	utility::vector1< SilentStructOP > silent_structs;
	silent_structs.reserve( structure_map_.size() );

	for ( auto const & map_elem : structure_map_ ) {
		silent_structs.push_back( map_elem.second );
	}
	clear_structure_map();

	using ObjexxFCL::string_of;
	int count = 1;
	for ( auto const & ss : silent_structs ) {
		std::string new_tag = ss->decoy_tag().substr(0,2) + string_of( count );
		ss->decoy_tag( new_tag );
		structure_map_[ new_tag ] = ss;
		structure_list_.push_back( ss );
		++count;
	}
} // renumber_all_decoys

bool SilentFileData::write_silent_struct(
	SilentStruct & s,
	std::string const & filename,
	bool bWriteScoreOnly
) const {
	bool success = false;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::io::ozstream output;

	std::stringstream header;
	if ( bWriteScoreOnly ) s.print_score_header( header );
	else s.print_header( header ); //this is SEQUENCE line + print_score_header

#ifndef USEMPI
	//this seems to segfault in mpi mode...
	//moreover it is a rather specialized feature which I'd rather see in a more highlevel code OL 6/2/09
	// which is also not consistently implemented, since it is missing in the write_all method.
	if ( store_argv_in_file() ) {
		std::string cmd_executed = option.get_argv();
		write_comment( header, cmd_executed );
	}
#endif

	// std::string header_string( header.str() );
	// std::cout << header_string << std::endl;

	//atomic version of opening procedure... writes header only if file is new -- appends if file already exists.
	output.open_append_if_existed( filename, header );

	// to prevent jobs from proceeding blithely when a directory doesn't exist.
	if ( !output.good() ) {
		utility_exit_with_message( "Could not make "+filename );
	}

	write_silent_struct( s, output, bWriteScoreOnly );

	output.close();

	success = true;
	return success;
}

/// @details Unlike write_silent_struct(), this method writes both the header
/// and the body of the silent file to the given stream.  Therefore this is
/// probably the more appropriate method to use in high-level codes.
bool SilentFileData::_write_silent_struct(
	SilentStruct & s,
	std::ostream & out,
	bool bWriteScoreOnly) const {
	s.print_header(out);
	return write_silent_struct(s,out,bWriteScoreOnly);
}

/// @details Note that this method does not write the header to the given
/// stream.  As such, I would guess that method is meant for internal use only,
/// even though it isn't declared as private.  If you want to write a complete
/// silent file to a stream, use the _write_silent_struct() method instead.
bool SilentFileData::write_silent_struct(
	SilentStruct & s,
	std::ostream & out,
	bool bWriteScoreOnly
) const {
	bool success = false;

	EnergyNamesOP enames;
	// set minimal set of EnergyNames if we know about them ...
	if ( has_shared_silent_data( energynames ) ) {
		enames = EnergyNamesOP(
			utility::pointer::static_pointer_cast< core::io::silent::EnergyNames > ( get_shared_silent_data( energynames ) )
		);
	} else {
		enames = EnergyNamesOP( new EnergyNames() );
		enames->energy_names( s.energy_names().energy_names() );
		set_shared_silent_data( energynames, enames );
	}

	if ( strict_column_mode() ) {
		s.set_valid_energies( enames->energy_names() );
	} else {
		// check to make sure the scores in this silent-struct are the same as the
		// last one that I printed. If not, re-print the header! Also, if the
		// user says to print a SCORE header for each decoy, print that header.
		if ( s.energy_names().energy_names() != enames->energy_names() ||
				//option[ out::file::silent_print_all_score_headers ]()
				options_.print_all_score_headers()
				) {
			s.print_header( out );
			enames = EnergyNamesOP( new EnergyNames() );
			enames->energy_names( s.energy_names().energy_names() );
			set_shared_silent_data( energynames, enames );
		}
	}

	s.print_scores( out );
	if ( !bWriteScoreOnly ) {
		s.print_residue_numbers( out );
		s.print_submotif_info(out);
		s.print_conformation( out );
	}

	if ( !bWriteScoreOnly ) { // other silent structs (used in stepwise modeling)
		for ( Size n = 1; n <= s.other_struct_list().size(); n++ ) {
			SilentStructOP other_struct = s.nonconst_other_struct_list()[n];
			other_struct->set_valid_energies( enames->energy_names() ); // prevents rewrite of score header, etc.
			write_silent_struct( *other_struct, out, false );
		}
	}

	// the following code was written to satisfy the American No Child
	// Left Behind act, so that the SilentFileData::write_silent_struct
	// method will always be successful.
	success = true;
	return success;
}

void SilentFileData::write_comment(
	std::ostream & out,
	std::string const & line
) const {
	out << "# " << line << '\n';
}

int SilentFileData::nres() const {
	return static_cast< int > ( begin()->nres() );
}

SilentStructOP SilentFileData::operator[] (std::string tag) {
	runtime_assert( has_tag( tag ) );
	return structure_map_[ tag ];
}

void SilentFileData::write_all(
	std::string const & filename, bool bWriteScoreOnly
) const {
	if ( begin() == end() ) return;

	utility::io::ozstream output;

	if ( begin() != end() ) {
		std::stringstream header;
		begin()->print_header( header );
		output.open_append_if_existed( filename, header );
	}

	for ( core::io::silent::SilentFileData::const_iterator iter = begin(),
			it_end = end(); iter != it_end; ++iter ) {
		write_silent_struct( **iter, output, bWriteScoreOnly );
		//    iter->print_scores      ( output );
		//    if ( !bWriteScoreOnly ) {
		//     iter->print_conformation( output );
		//    }
	} // for ( iter )
} // write_all


bool
SilentFileData::read_file(
	std::string const & filename
) {

	bool success = _read_file(
		filename,
		//!option[ OptionKeys::in::file::silent_read_through_errors ]() /*throw_exception_on_bad_structs*/
		!options_.read_through_errors()
	);
	//if ( !success && !option[ OptionKeys::in::file::silent_read_through_errors ]() ) {
	if ( !success && !options_.read_through_errors() ) {
		throw utility::excn::EXCN_BadInput("no success reading silent file "+filename);
	}
	return success;
} // read_file


bool
SilentFileData::_read_file(
	std::string const & filename,
	bool throw_exception_on_bad_structs /*default false*/
) {

	using std::string;
	using utility::vector1;

	//empty tag-vector to signal read ALL
	utility::vector1< std::string > tags_wanted;

	utility::vector1< std::string > all_tags = read_tags_fast( filename );
	//if ( option[ in::file::silent_select_random ].user() ) {
	if ( options_.select_random_set() ) {
		//core::Size const n_wanted( option[ in::file::silent_select_random ]() );
		core::Size const n_wanted( options_.select_random() );
		runtime_assert( n_wanted > 0.0 );
		runtime_assert( n_wanted <= all_tags.size() );

		tags_wanted = numeric::random::reservoir_sample< string >(
			all_tags, n_wanted, numeric::random::rg()
		);
		tr.Debug << "selected " << tags_wanted.size() << " tags from " <<
			all_tags.size() << std::endl;
	}

	//if ( option[ in::file::silent_select_range_start ]() >= 0 ) {
	if ( options_.select_range_start() >= 0 ) {

		// get start index and limit to size of silent file (all_tags.size() ) by wrapping around (rather then by throwing an error)
		//core::Size range_start = ( option[ in::file::silent_select_range_start ]() *
		// option[ in::file::silent_select_range_mul ] ) % all_tags.size();
		core::Size range_start = ( options_.select_range_start() * options_.select_range_mul() ) % all_tags.size();

		core::Size range_end   = range_start + options_.select_range_len();

		// setting in::file::silent_select_range_len to -1 means run to the end.
		//if ( option[ in::file::silent_select_range_len ]() < 0 ) range_end = all_tags.size();
		if ( options_.select_range_len() < 0 ) range_end = all_tags.size();

		// do a range check just in case.
		range_end = std::min( (int)range_end, (int)all_tags.size() );

		tr << "Reading range: "  << range_start << " <= i < " << range_end << std::endl;
		for ( core::Size position = range_start; position < range_end; position++ ) {
			tags_wanted.push_back( all_tags[position+1] ); // the +1 converts to 1 based counting ov the utility::vector1
		}

	}

	if ( all_tags.size() == 0 ) return true;

	bool success = _read_file( filename, tags_wanted, throw_exception_on_bad_structs  );
	return success;
} // read_file


bool
SilentFileData::read_file(
	std::string const & filename,
	utility::vector1< std::string > const & tags
) {

	bool success = _read_file(
		filename,
		tags,
		//!option[ OptionKeys::in::file::silent_read_through_errors ]() /*throw_exception_on_bad_structs*/
		!options_.read_through_errors()
	);
	//if ( !success && !option[ OptionKeys::in::file::silent_read_through_errors ]() ) {
	if ( !success && !options_.read_through_errors() ) {
		throw utility::excn::EXCN_BadInput("no success reading silent file "+filename);
	}
	return success;
}


bool
SilentFileData::_read_file(
	std::string const & filename,
	utility::vector1< std::string > const & tags,
	bool throw_exception_on_bad_structs /*default false*/
) {
	bool success = false; // default value for the early return statements
	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		if ( throw_exception_on_bad_structs ) {
			throw utility::excn::EXCN_FileNotFound( filename );
		} else {
			utility_exit_with_message(
				"ERROR:: Unable to open silent_input file: '" +
				filename + "'"
			);
		}
	}
	success = read_stream(data,tags,throw_exception_on_bad_structs,filename);
	return success;
}


bool
SilentFileData::read_stream(
	std::istream & data,
	utility::vector1< std::string > const & tags,
	bool throw_exception_on_bad_structs, /*default false*/
	std::string filename /** for error reporting and suppressing bitflip**/
) {

	// put the tags on a set to avoid repetition and fast search
	std::set<std::string> tagset(tags.begin(), tags.end());
	std::string sequence_line, score_type_line, line;
	bool success = false;
	// read SEQUENCE line:
	std::getline( data, sequence_line );
	if ( sequence_line.substr(0,9) != "SEQUENCE:" ) {
		tr.Error << "bad format in first line of silent file " << filename << " (function read_stream):"
			<< std::endl;
		tr.Error << sequence_line << std::endl;
		return success;
	}
	check_if_rna_from_sequence_line( sequence_line );

	// read header SCORE line
	getline( data, score_type_line );
	if ( !line_starts_with_score( score_type_line ) || !is_score_type_line( score_type_line ) ) {
		tr.Error << "bad format in second line of silent file " << filename << " (function read_stream):" << std::endl;
		tr.Error << score_type_line << std::endl;
		return success;
	}


	getline( data, line );

	utility::vector1< std::string > mylines; // used for initialization of particular Silent-Structs
	std::map < SilentFileHeaderLine, std::string > header_lines; //Used for initialization of particular Silent-Structs from header lines.  These persist from silent structure to silent structure.
	header_lines[ SFHEADER_SEQUENCE_LINE ] = sequence_line;
	header_lines[ SFHEADER_SCORE_TYPE_LINE ] = score_type_line;

	if ( verbose_ ) {
		if ( !tagset.empty()/*size()*/ ) {
			tr.Info << "Reading " << tagset.size()
				<< " structures from " << filename
				<< std::endl;
		} else  tr.Info << "Reading all structures from " << filename << std::endl;
	}
	bool all_tags = tagset.empty();//tagset.size() == 0; //if no tags selected we read all decoys

	// start looping over the structures
	bool line_ok( true );
	bool reading_header( true ); //True if we're reading header lines, false if we're reading data lines.
	bool premature_stop( false ); // True if we were reading header lines when we reached the end of a file.

	//tr.Debug << "Entering loop..." << std::endl; //DELETE ME

	while ( line_ok ) {
		//tr.Debug << line << std::endl; //DELETE ME

		if ( line.substr(0,1) == "#" ) {
			tr.Info << "Read a comment line from " << filename << std::endl;
			tr.Info << line << std::endl;
			comment_lines_.push_back( line );
			getline(data,line);
			line_ok = data.good();
			if ( !line_ok ) { premature_stop = true; break; }
			continue;
		}

		if ( reading_header ) {
			//tr.Debug << "Trying to read a header line..." << std::endl; //DELETE ME

			if ( line_starts_with_score( line ) ) { //This is a score type or score line:
				header_lines[ is_score_type_line(line) ? SFHEADER_SCORE_TYPE_LINE : SFHEADER_SCORE_LINE ] = line;
				getline(data,line);
				line_ok = data.good();
				if ( !line_ok ) { premature_stop = true; break; }
				continue;
			} else if ( line_starts_with_sequence( line ) ) {
				header_lines[ SFHEADER_SEQUENCE_LINE ] = line;
				getline(data,line);
				line_ok = data.good();
				if ( !line_ok ) { premature_stop = true; break; }
				continue;
			} else if ( line_starts_with_remark( line ) ) {
				read_silent_struct_type_from_remark( line, true /* header */);
				read_full_model_parameters_from_remark( line, true /* header */ );
				mylines.push_back(line);
				getline(data,line);
				line_ok = data.good();
				if ( !line_ok ) { premature_stop = true; break; }
				continue;
			} else if (
					line_starts_with_other(line) ||
					line_starts_with_other_header_string(line)
					) {
				mylines.push_back(line);
				getline(data,line);
				line_ok = data.good();
				if ( !line_ok ) { premature_stop = true; break; }
				continue;
			} else {
				if ( !data.good() || data.eof() ) { premature_stop = true; break; }
				reading_header=false;
			}
		}

		if ( !reading_header ) {
			// If we're not reading header lines...
			//tr.Debug << "Trying to read a non-header line..." << std::endl; //DELETE ME

			if ( line_starts_with_score( line ) || line_starts_with_sequence( line ) || line_starts_with_other( line ) || line_starts_with_remark( line ) || line_starts_with_other_header_string( line ) ) {
				//If another header has started after we've been reading data lines, then we're done reading data lines for the current silent struct, and it's time to
				//parse what we've read.

				reading_header = true;
				lines_from_header_line_collection( header_lines, mylines ); //Copy the header lines to the lines collection.

				SilentStructOP tmp_struct( create_SilentStructOP() );
				if ( options_.force_silent_bitflip_on_read() && filename != "suppress_bitflip" ) {
					tmp_struct->set_force_bitflip(true); //Option to force flipping from big-endian to little-endian or the converse.
				}

				bool init_good = tmp_struct->init_from_lines( mylines, *this );

				if ( !init_good && throw_exception_on_bad_structs ) {
					throw utility::excn::EXCN_BadInput(
						"failure to read decoy "+tmp_struct->decoy_tag()+
						" from silent-file " + filename
					);
				}

				if ( init_good ) {
					//tr.Debug << "candidate structure " << tmp_struct->decoy_tag()
					// << std::endl;
					PROF_START( basic::SILENT_READ_TAG_TEST );

					// Look for the tag in the tagset, if it is present set
					// add the structure and remove the tag from the tagset.
					// if afer removal of the tag tagset becames empty break
					// because that means that all the structures that we want
					// have been read.
					bool good_tag = false;
					if ( !all_tags ) {
						auto tag_it = tagset.find(tmp_struct->decoy_tag());
						if ( tag_it != tagset.end() ) {
							good_tag = true;
							tagset.erase(tag_it); //expensive restructering of set -- how about saving a bunch of bools to do this book-keeping.
						}
					}
					bool add_struct( init_good && ( all_tags || good_tag ));
					PROF_STOP( basic::SILENT_READ_TAG_TEST );

					if ( record_source_ ) tmp_struct->add_comment("SOURCE",filename);

					if ( add_struct ) {
						add_structure( tmp_struct );
						// check if there are any tags left in the tagset, if not break
						// since there is no need to read the rest of the file.
					}
				}

				mylines.clear();
				mylines.reserve( 10 ); //+10 since there are at least the SCORE lines extra, maybe some REMARK lines
				premature_stop = true; // Temporarily consider this to be a premature stop, since we've read a header line.  This is switched back in a few lines.
			} //if line starts with "SCORE:", "SEQUENCE:", or "OTHER:"

			mylines.push_back( line );

			getline(data,line);
			line_ok = data.good();
			if ( line_ok ) premature_stop = false;

		} //if we're not reading the header

		// If the loop is interrupted prematurely, we must set premature_stop = true to avoid segfault.

	} // while( getline(data,line) )

	// don't forget to initialize last structure!
	if ( !premature_stop ) {
		lines_from_header_line_collection( header_lines, mylines ); //Copy the header lines to the lines collection.
		SilentStructOP tmp_struct( create_SilentStructOP() );
		if ( options_.force_silent_bitflip_on_read() && filename != "suppress_bitflip" ) {
			tmp_struct->set_force_bitflip(true); //Option to force flipping from big-endian to little-endian or the converse.
		}
		bool init_good = tmp_struct->init_from_lines( mylines, *this );
		if ( !init_good && throw_exception_on_bad_structs ) throw utility::excn::EXCN_BadInput( "failure to read decoy "+tmp_struct->decoy_tag()+" from silent-file " +filename);

		bool good_tag = false;
		std::set<std::string>::iterator tag_it = tagset.find(tmp_struct->decoy_tag());
		if ( tag_it != tagset.end() ) {
			good_tag = true;
		}

		bool add_struct( init_good && ( all_tags || good_tag ));

		if ( add_struct ) {
			add_structure( tmp_struct );
		}
	}

	if ( verbose_ ) {
		tr.Info << "Finished reading " << structure_map_.size() << " structures from "
			<< filename << std::endl;
		tr.flush();
	}

	success = true;

	return success;
}


/// @detail The first remarks line in a silent file block
///described the type of silent file that is comming. For example,
/// REMARK BINARY SILENTFILE
bool
SilentFileData::read_silent_struct_type_from_remark(
	std::string const& line,
	bool const header
) {

	bool changed( false );

	if ( line.substr(0,6) != "SCORE:" ) {
		if ( ( line.find( "BINARY_SILENTFILE" ) != std::string::npos ) ||
				( line.find ("BINARY SILENTFILE" ) != std::string::npos ) ) {
			silent_struct_type_ = "binary";
			changed = true;
		} else if ( header && ( line.find( "RNA" ) != std::string::npos ) &&
				silent_struct_type_ != "binary_rna" && silent_struct_type_ != "rna" ) {
			silent_struct_type_ = "rna";
			changed = true;
		} else if ( ( line.find( "PROTEIN_SILENTFILE" ) != std::string::npos ) || ( line.find( "PROTEIN SILENTFILE" ) != std::string::npos ) ) {
			if ( line.find("SINGLE") != std::string::npos ) {
				silent_struct_type_ ="protein_float";
			} else {
				silent_struct_type_ = "protein";
			}
			changed = true;
		}
	} else if ( ( line.find ("SOCREJUMP SILENTFILE" ) != std::string::npos ) ) {
		silent_struct_type_ = "score_jump";
		changed = true;
	}

	if ( changed ) tr.Trace << "found new silent_struct_type_ " << silent_struct_type_ << " from line " << line << std::endl;

	return changed;
}

/// @brief Look for FULL_MODEL_PARAMETERS in REMARK line from header.
bool
SilentFileData::read_full_model_parameters_from_remark(
	std::string const& line,
	bool const header
) {
	runtime_assert( header );
	Size pos = line.find( "FULL_MODEL_PARAMETERS" );
	if ( pos == std::string::npos ) return false;

	using namespace core::pose::full_model_info;
	std::istringstream ss( line.substr( pos ) );
	FullModelParametersOP full_model_parameters_from_line( new FullModelParameters );
	ss >> *full_model_parameters_from_line;
	full_model_parameters_ = full_model_parameters_from_line;
	return true;
}

/// @brief This is somewhat redundant if the silent file has a "REMARK RNA" line, but some older silent files didn't do that.
bool SilentFileData::check_if_rna_from_sequence_line( std::string const& line ) {
	std::istringstream l( line );
	std::string dummy, sequence;
	l >> dummy;
	l >> sequence;
	bool is_rna( true );
	for ( Size n = 0; n < sequence.size(); n++ ) {
		if ( sequence[n] != 'a' && sequence[n] != 'g' && sequence[n] != 'c' && sequence[n] != 'u' ) {
			is_rna = false; break;
		}
	}

	if ( is_rna ) silent_struct_type_ = "rna";
	// could also be binary_rna type, but historically we always made sure to include a "REMARK BINARY_RNA" line in those silent files, and
	// this will be picked up later.

	return false;
}

/// @brief Are the first six characters of a line "SCORE:"?
/// @details Both score type lines and score lines start with "SCORE:".
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
SilentFileData::line_starts_with_score(
	std::string const &line
) const {
	return line.substr(0,6) == "SCORE:";
}

/// @brief Are the first nine characters of a line "SEQUENCE:"?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
SilentFileData::line_starts_with_sequence(
	std::string const &line
) const {
	return line.substr(0,9) == "SEQUENCE:";
}

/// @brief Are the first six characters of a line "OTHER:"?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
SilentFileData::line_starts_with_other(
	std::string const &line
) const {
	return line.substr(0,6) == "OTHER:";
}

/// @brief Are the first six characters of a line "REMARK "?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
SilentFileData::line_starts_with_remark(
	std::string const &line
) const {
	return line.substr(0,7) == "REMARK ";
}

/// @brief Does a line start with another header string?
/// @details Includes strings like "ANNOTATED_SEQUENCE", "FOLD_TREE", etc.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
SilentFileData::line_starts_with_other_header_string(
	std::string const &line
) const {
	return (
		line.substr(0,10) == "FOLD_TREE " ||
		line.substr(0,2) == "RT" ||
		line.substr(0,24) == "NONCANONICAL_CONNECTION:" ||
		line.substr(0,19) == "ANNOTATED_SEQUENCE:" ||
		line.substr(0,4) == "JUMP" ||
		line.substr(0,13) == "SYMMETRY_INFO" ||
		line.substr(0,13) == "CHAIN_ENDINGS" ||
		line.substr(0,7) == "RES_NUM" ||
		line.substr(0,11) == "SEGMENT_IDS" ||
		line.substr(0,13) == "SUBMOTIF_INFO"
	);
}


/// @brief Is a line (that is assumed to start with "SCORE:") a score type line or a score line?
/// @details Considered a score type line if and only if the second whitespace-separated chunk is "score" or "total_score"; otherwise, considered a score line.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
SilentFileData::is_score_type_line(
	std::string const &line
) const {
	std::istringstream linestream( line );
	std::string buffer;
	linestream >> buffer;
	if ( !linestream.good() || linestream.eof() ) return false;
	linestream >> buffer;
	if ( linestream.bad() || linestream.fail() ) return false;
	return buffer == "score" || buffer == "total_score" ;
}

/// @brief Given a list of non-header lines read from a silent file and a map of certain types of single-occurance header lines to
/// strings, copy the single-occurance header lines to the top of the list of lines.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
SilentFileData::lines_from_header_line_collection(
	std::map< SilentFileHeaderLine, std::string> const & header_lines,
	utility::vector1< std::string> & all_lines
) const {
	runtime_assert_string_msg( header_lines.count( SFHEADER_SEQUENCE_LINE ), "Error in core::io::silent::SilentFileData::lines_from_header_line_collection(): The silent structure must have a \"SEQUENCE:\" line." );
	runtime_assert_string_msg( header_lines.count( SFHEADER_SCORE_TYPE_LINE ), "Error in core::io::silent::SilentFileData::lines_from_header_line_collection(): The silent structure must have a \"SCORE:\" line listing score types." );
	//runtime_assert_string_msg( header_lines.count( SFHEADER_SCORE_LINE ), "Error in core::io::silent::SilentFileData::lines_from_header_line_collection(): The silent structure must have a \"SCORE:\" line listing scores." );

	//Must be inserted in reverse order (each to the start of the list):
	if ( header_lines.count( SFHEADER_SCORE_LINE ) ) all_lines.insert( all_lines.begin(), header_lines.at(SFHEADER_SCORE_LINE) );
	all_lines.insert( all_lines.begin(), header_lines.at(SFHEADER_SCORE_TYPE_LINE) );
	all_lines.insert( all_lines.begin(), header_lines.at(SFHEADER_SEQUENCE_LINE) );

	/*if(tr.Debug.visible()) {
	tr.Debug << "\n";
	for(core::Size i=1, imax=all_lines.size(); i<=imax; ++i) {
	tr.Debug << all_lines[i] << "\n";
	}
	tr.Debug << std::endl;
	tr.Debug.flush();
	}*/

}


/// @brief creates a SilentStructOP using command-line options. Just a wrapper
/// around SilentStructFactory::get_instance()->get_silent_struct_in().
SilentStructOP SilentFileData::create_SilentStructOP() {

	//static SilentStructFactory ssf;
	SilentStructOP new_ss_op;

	if ( silent_struct_type_ == "" ) new_ss_op = SilentStructFactory::get_instance()->get_silent_struct_in( options_ );
	else                           new_ss_op = SilentStructFactory::get_instance()->get_silent_struct( silent_struct_type_, options_ );

	if ( full_model_parameters_ ) new_ss_op->set_full_model_parameters( full_model_parameters_ );

	return new_ss_op;
}

/// @brief Returns the SharedSilentDataOP associated with the given type.
SharedSilentDataOP SilentFileData::get_shared_silent_data(
	SharedSilentDataType ssdt
) const {
	return shared_silent_data_[ ssdt ];
}

void SilentFileData::set_shared_silent_data(
	SharedSilentDataType ssdt,
	SharedSilentDataOP ssd_op
) const {
	shared_silent_data_[ ssdt ] = ssd_op;
}

bool SilentFileData::has_shared_silent_data( SharedSilentDataType ssdt ) const {
	return ( shared_silent_data_.find( ssdt ) != shared_silent_data_.end() );
}

void
SilentFileData::reverse_score_filter(
	Real const score_fraction
) {
	// make sure fraction
	runtime_assert( score_fraction <=  0 );
	runtime_assert( score_fraction >= -1 );

	// go through all of the scores, and calculate the threshold for accepting a
	// decoy by score.
	utility::vector1< Real > scores;
	for ( iterator iter = begin(), it_end = end(); iter != it_end; ++iter ) {
		scores.push_back( iter->get_energy( "score" ) );
	}

	//Real local_score_fraction( score_fraction * -1 );
	std::sort( scores.begin(), scores.end() );
	std::reverse( scores.begin(), scores.end() );
	Size const idx( static_cast< Size > ( -1 * score_fraction * scores.size() ) );
	Real const boundary( *( scores.begin() + idx ) );
	tr.Debug << "reverse_score_filter: " << std::endl;
	tr.Debug << "filtering for decoys with score worse than " << boundary << std::endl;

	Structure_Map new_structure_map_;
	utility::vector1 < SilentStructOP > new_structure_list_;
	for ( iterator iter = begin(), it_end = end(); iter != it_end; ++iter ) {
		if ( !iter->has_energy( "score" ) ) {
			std::string msg(
				"Error: can't find score in SilentStuct!"
			);
			msg += "\nSilentStruct object has scores:\n";

			// little bit of indirection here to get the SilentStruct scores.
			std::ostringstream mystream;
			iter->print_score_header( mystream );
			mystream << std::endl;
			iter->print_scores( mystream );
			mystream << std::endl;
			msg += mystream.str();

			utility_exit_with_message( msg );
		}

		if ( iter->get_energy( "score" ) > boundary ) {
			new_structure_map_[ iter->decoy_tag() ] = *iter;
			new_structure_list_.push_back(*iter);
		}
	}

	structure_map_ = new_structure_map_;
	structure_list_ = new_structure_list_;
} // score_filter

void
SilentFileData::score_filter(
	Real const score_fraction
) {
	if ( score_fraction < 0 ) {
		return reverse_score_filter( score_fraction );
	}
	// make sure fraction is in bounds
	runtime_assert( score_fraction >= 0 );
	runtime_assert( score_fraction <= 1 );

	// go through all of the scores, and calculate the threshold for accepting a
	// decoy by score.
	utility::vector1< Real > scores;
	for ( iterator iter = begin(), it_end = end(); iter != it_end; ++iter ) {
		scores.push_back( iter->get_energy( "score" ) );
	}

	std::sort( scores.begin(), scores.end() );
	Size const idx( static_cast< Size > (score_fraction * scores.size() ) );
	Real const boundary( *( scores.begin() + idx ) );
	tr.Debug << "score_filter: " << std::endl;
	tr.Debug << "filtering for decoys with score worse than " << boundary << std::endl;

	Structure_Map new_structure_map_;
	utility::vector1 < SilentStructOP > new_structure_list_;
	for ( iterator iter = begin(), it_end = end(); iter != it_end; ++iter ) {
		if ( !iter->has_energy( "score" ) ) {
			std::string msg(
				"Error: can't find score in SilentStuct!"
			);
			msg += "\nSilentStruct object has scores:\n";

			// little bit of indirection here to get the SilentStruct scores.
			std::ostringstream mystream;
			iter->print_score_header( mystream );
			mystream << std::endl;
			iter->print_scores( mystream );
			mystream << std::endl;
			msg += mystream.str();

			utility_exit_with_message( msg );
		}

		if ( iter->get_energy( "score" ) < boundary ) {
			new_structure_map_[ iter->decoy_tag() ] = *iter;
			new_structure_list_.push_back(*iter);
		}
	}

	structure_map_ = new_structure_map_;
	structure_list_ = new_structure_list_;
} // score_filter

////////////////////////////////////////////////////////////////////////////////////////
// This is a little funky -- the "tags" by which the silent structs are
// keyed in the SilentFileData need to be renamed. Weird because the silent
// struct's own "decoy_tag" will no longer match their key within the silent file data.
////////////////////////////////////////////////////////////////////////////////////////
void
SilentFileData::order_by_energy()
{

	using namespace core::io::silent;

	// go through all of the scores, and then order.
	typedef std::list< std::pair< std::pair< Real, Size >, SilentStructOP > > ScoreTagList;
	ScoreTagList score_tag_list;
	Size count( 0 ); // to break ties in energies.
	for ( iterator iter = begin(), it_end = end(); iter != it_end; ++iter ) {
		Real const & silent_score = (*iter)->get_energy( "score" );
		score_tag_list.push_back( std::make_pair( std::make_pair( silent_score, ++count) , *iter /*SilentStructOP*/ ) );
	}

	score_tag_list.sort();

	Structure_Map new_structure_map_;
	utility::vector1 < SilentStructOP > new_structure_list_;
	count = 0;
	for ( ScoreTagList::const_iterator iter = score_tag_list.begin();
			iter != score_tag_list.end(); ++iter ) {
		SilentStructOP const & silent_struct_op( iter->second );

		///////////////////////////////////////////////
		// Note that we have to supply a new tag when
		// we reorder the silent structs by energy...
		//  maps are ordered by key in C++.
		///////////////////////////////////////////////
		std::string const new_tag = "S_" +
			ObjexxFCL::lead_zero_string_of( count++,6 )+
			silent_struct_op->decoy_tag();

		new_structure_map_[ new_tag ] = silent_struct_op;
		new_structure_list_.push_back( silent_struct_op );
	}

	structure_map_ = new_structure_map_;
	structure_list_ = new_structure_list_;

} // score_filter


} // namespace silent
} // namespace io
} // namespace core
