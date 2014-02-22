// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SharedSilentData.hh>

#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/prof.hh>

#include <numeric/random/random.hh>
#include <numeric/random/reservoir_sample.hh>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {

static basic::Tracer tr("core.io.silent.SilentFileData");
static numeric::random::RandomGenerator RG(21458);

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
			it_end = structure_map_.end(); it != it_end; ++it ){
		tag_list.push_back( it->decoy_tag() );
	}
	return tag_list;
} // tags


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

	std::string line,tag;
	getline( data, line ); // sequence line
	getline( data, line ); // score line

	while( getline(data,line) ) {
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
	//	utility::vector1< SilentStructOP >& decoys_in_file,
	//	bool ignore_decoys
) const {

	//open file
	utility::io::izstream data( filename().c_str() );
	if ( !data.good() ) {
// 		utility_exit_with_message(
// 			"ERROR: Unable to open silent_input file: '" + filename() + "'"
// 		);
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
	while( getline(data,line) ) {
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
				//				ignore_this_decoy = ignore_decoys;
 				if ( !return_last ) {
					tags_in_file.push_back( current_tag );
				} else final_tag = current_tag;
				if ( return_first ) return true;  //if we want decoys this will be wrong
			}
		} // SCORE: header of decoy
		//		if ( !ignore_this_decoy ) {
		//			utility_exit_with_message( "SilentFileData::matched_tag() function stubbed out -- cannot return decoys yet" );
		//		}
	} // while( getline(data,line) )
// 	if ( return_last ) {
// 		if ( current_tag.find(expression)!=std::string::npos ) {
// 			tags_in_file.push_back( current_tag );
// 			return true;
// 		}
// 	}
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

	std::string const & new_tag( new_struct->decoy_tag() );

	if ( !has_tag( new_tag ) ) {
		structure_map_[ new_tag ] = new_struct;
		tr.Debug << "added structure with tag " << new_tag << std::endl;
		tr.Debug << "now have " << structure_map_.size() << " structures."
						 << std::endl;
		tr.flush();
	} else {
		SilentStructOP new_struct_clone = new_struct->clone();
		add_structure_replace_tag_if_necessary( new_struct_clone );
	}

} // add_structure

void SilentFileData::add_structure(
	SilentStruct const & new_struct
) {
	add_structure( new_struct.clone() );
} // add_structure

void SilentFileData::renumber_all_decoys() {
	utility::vector1< SilentStructOP > silent_structs;
	silent_structs.reserve( structure_map_.size() );

	Structure_Map::iterator map_iter;
	for ( map_iter = structure_map_.begin(); map_iter != structure_map_.end(); ++map_iter ) {
		silent_structs.push_back( map_iter->second );
	}
	clear_structure_map();

	using ObjexxFCL::string_of;
	int count = 1;
	utility::vector1< SilentStructOP >::const_iterator ss_iter;
	for ( ss_iter = silent_structs.begin(); ss_iter != silent_structs.end(); ++ss_iter ) {
		std::string new_tag = (*ss_iter)->decoy_tag().substr(0,2) + string_of( count );
		(*ss_iter)->decoy_tag( new_tag );
		structure_map_[ new_tag ] = *ss_iter;
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
	if ( bWriteScoreOnly )	s.print_score_header( header );
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

	//	std::string header_string( header.str() );
	//	std::cout << header_string << std::endl;

	//atomic version of opening procedure... writes header only if file is new -- appends if file already exists.
	output.open_append_if_existed( filename, header );

	// to prevent jobs from proceeding blithely when a directory doesn't exist.
	if ( !output.good() ){
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
			static_cast< EnergyNames * > ( get_shared_silent_data( energynames )() )
		);
	} else {
		enames = new EnergyNames();
		enames->energy_names( s.energy_names().energy_names() );
		set_shared_silent_data( energynames, enames );
	}

	if ( strict_column_mode() ) {
		s.set_valid_energies( enames->energy_names() );
	} else {
		// check to make sure the scores in this silent-struct are the same as the
		// last one that I printed. If not, re-print the header! Also, if the
		// user says to print a SCORE header for each decoy, print that header.
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if ( s.energy_names().energy_names() != enames->energy_names() ||
				option[ out::file::silent_print_all_score_headers ]()
		) {
			s.print_header( out );
			enames = new EnergyNames();
			enames->energy_names( s.energy_names().energy_names() );
			set_shared_silent_data( energynames, enames );
		}
	}

	s.print_scores( out );
	if ( !bWriteScoreOnly ) {
		s.print_residue_numbers( out );
		s.print_conformation( out );
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

	if ( begin() != end() )	{
		std::stringstream header;
		begin()->print_header( header );
		output.open_append_if_existed( filename, header );
	}

	for ( core::io::silent::SilentFileData::const_iterator iter = begin(),
					it_end = end(); iter != it_end; ++iter ) {
		write_silent_struct( **iter, output, bWriteScoreOnly );
		// 			iter->print_scores      ( output );
		// 			if ( !bWriteScoreOnly ) {
		// 				iter->print_conformation( output );
		// 			}
	} // for ( iter )
} // write_all


bool
SilentFileData::read_file(
	std::string const & filename
) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool success = _read_file(
		filename,
		!option[ OptionKeys::in::file::silent_read_through_errors ]() /*throw_exception_on_bad_structs*/
	);
	if ( !success && !option[ OptionKeys::in::file::silent_read_through_errors ]() ) {
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
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//empty tag-vector to signal read ALL
	utility::vector1< std::string > tags_wanted;

	utility::vector1< std::string > all_tags = read_tags_fast( filename );
	if ( option[ in::file::silent_select_random ].user() ) {
		core::Size const n_wanted( option[ in::file::silent_select_random ]() );
		runtime_assert( n_wanted > 0.0 );
		runtime_assert( n_wanted <= all_tags.size() );

		tags_wanted = numeric::random::reservoir_sample< string >(
			all_tags, n_wanted, RG
		);
		tr.Debug << "selected " << tags_wanted.size() << " tags from " <<
			all_tags.size() << std::endl;
	}

	if ( option[ in::file::silent_select_range_start ]() >= 0 ) {

	// get start index and limit to size of silent file (all_tags.size() ) by wrapping around (rather then by throwing an error)
	core::Size range_start = ( option[ in::file::silent_select_range_start ]() *
			option[ in::file::silent_select_range_mul ] ) % all_tags.size();
	core::Size range_end   = range_start + option[ in::file::silent_select_range_len ]();

	// setting in::file::silent_select_range_len to -1 means run to the end.
	if( option[ in::file::silent_select_range_len ]() < 0 ) range_end = all_tags.size();

	// do a range check just in case.
	range_end = std::min( (int)range_end, (int)all_tags.size() );

	tr << "Reading range: "  << range_start << " <= i < " << range_end << std::endl;
	for( core::Size position = range_start; position < range_end; position++){
		tags_wanted.push_back( all_tags[position+1] ); // the +1 converts to 1 based counting ov the utility::vector1
	}

  }

  if( all_tags.size() == 0 ) return true;

	bool success = _read_file( filename, tags_wanted, throw_exception_on_bad_structs  );
	return success;
} // read_file


bool
SilentFileData::read_file(
	std::string const & filename,
	utility::vector1< std::string > const & tags
) {
using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool success = _read_file(
    filename,
		tags,
		!option[ OptionKeys::in::file::silent_read_through_errors ]()	/*throw_exception_on_bad_structs*/
	);
	if ( !success && !option[ OptionKeys::in::file::silent_read_through_errors ]() ) {
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
														std::string filename /** for error reporting **/
){

	// put the tags on a set to avoid repetition and fast search
  std::set<std::string> tagset(tags.begin(), tags.end());
	std::string sequence_line, score_line, line;
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
	getline( data, score_line );
	if ( score_line.substr(0,7) != "SCORE: " ) {
		tr.Error << "bad format in second line of silent file " << filename << " (function read_stream):" << std::endl;
		tr.Error << score_line << std::endl;
		return success;
	}

	// read potential REMARK line
	getline( data, line );

	//////////////////////////////
	// Checking for remarks...
	// NOTE THAT THIS IS COMPLETELY TERRIBLE -- THESE TAGS SHOULD BE UNIFIED
	// WITH THE PRINT_HEADER() REMARKS THAT ARE IN THE INDIVIDUAL SILENT STRUCTS, ETC.
	// Probably should all be handled by the SilentStructFactory, but I don't have
	// time to refactor this now -- rhiju, Dec. 2009.

	// fixed a bug when reading silent-files with lines like this:
	// SCORE: score template description
	// SCORE: -2.00 1ERNA.pdb S_000000001


	read_silent_struct_type_from_remark( line, true );
	SilentStructOP tmp_struct = create_SilentStructOP();

	utility::vector1< std::string > mylines; // used for initialization of particular Silent-Structs
	mylines.push_back( sequence_line );
	mylines.push_back( score_line    );

	if ( tagset.size() )  tr.Info << "Reading " << tagset.size()
		<< " structures from " << filename
		<< std::endl;
	else  tr.Info << "Reading all structures from " << filename << std::endl;
	bool all_tags = tagset.size() == 0; //if no tags selected we read all decoys

	// start looping over the structures
	bool line_ok( true );
	while ( line_ok ) {
		if ( line.substr(0,1) == "#" ) {
			tr.Info << "Read a comment line from " << filename << std::endl;
			tr.Info << line << std::endl;
			comment_lines_.push_back( line );
			line_ok = getline(data,line);
			continue;
		}
		//is the next decoy coming up ? then stop collecting lines and initialize the silent-struct from the read lines
		if ( line.substr(0,10) == "SEQUENCE: " && mylines.size() == 2 && mylines[ 1 ].substr(0,10) == "SEQUENCE: " ) {
			//we have just read the two SEQUENCE: and SCORE: header lines and a new header starts... ignore this
			mylines.clear();
		}
		if ( ( line.substr(0,7) == "SCORE: " || line.substr(0,10) == "SEQUENCE: " ) && mylines.size() > 3 ) {

			bool init_good = tmp_struct->init_from_lines( mylines, *this );

			if ( !init_good && throw_exception_on_bad_structs ) {
				 	throw utility::excn::EXCN_BadInput(
						"failure to read decoy "+tmp_struct->decoy_tag()+
						" from silent-file " + filename
					);
			}

			if ( init_good ) {
				//tr.Debug << "candidate structure " << tmp_struct->decoy_tag()
				//	<< std::endl;
				PROF_START( basic::SILENT_READ_TAG_TEST );

				// Look for the tag in the tagset, if it is present set
				// add the structure and remove the tag from the tagset.
				// if afer removal of the tag tagset becames empty break
				// because that means that all the structures that we want
				// have been read.
				bool good_tag = false;
				if ( !all_tags ) {
					std::set<std::string>::iterator tag_it = tagset.find(tmp_struct->decoy_tag());
					if ( tag_it != tagset.end() ) {
						good_tag = true;
						tagset.erase(tag_it); //expensive restructering of set -- how about saving a bunch of bools to do this book-keeping.
			   	}
				}
				bool add_struct( init_good && ( all_tags || good_tag ));
				PROF_STOP( basic::SILENT_READ_TAG_TEST );

				if (record_source_) tmp_struct->add_comment("SOURCE",filename);

				if ( add_struct ) {
					add_structure( tmp_struct );
					// check if there are any tags left in the tagset, if not break
					// since there is no need to read the rest of the file.
				}

				tmp_struct = create_SilentStructOP();
			}

			mylines.clear();
			mylines.reserve( tmp_struct->sequence().size() + 10 ); //+10 since there are at least the SCORE lines extra, maybe some REMARK lines
		}

		if ( read_silent_struct_type_from_remark( line ) ) tmp_struct = create_SilentStructOP();

		mylines.push_back( line );

		line_ok = getline(data,line);
		/// in no case interrupt loop here, because then mylines will have incomplete silent-struct causing
		/// seqfault in the "init_from_lines" for last decoy
		/* CAUSES SEGFAULT  OL 12/10/11
 if(!all_tags && tagset.empty())
			 break;  */
	} // while( getline(data,line) )

	// don't forget to initialize last structure!
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

	tr.Info << "Finished reading " << structure_map_.size() << " structures from "
		<< filename << std::endl;
	tr.flush();

	success = true;

	return success;
}


///@detail The first remarks line in a silent file block
///described the type of silent file that is comming. For example,
///
bool
SilentFileData::read_silent_struct_type_from_remark(
	std::string const& line,
	bool header
) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool changed( false );

	if ( line.substr(0,6) != "SCORE:" ) {
		if (( line.find( "BINARY_SILENTFILE" ) != std::string::npos ) || ( line.find ("BINARY SILENTFILE" ) != std::string::npos )) {
			if ( line.find( "RNA" ) != std::string::npos || option[ in::file::residue_type_set ]() == "rna" ) {
				silent_struct_type_ = "binary_rna";
			} else {
				silent_struct_type_ = "binary";
			}
			changed = true;
		} else if ( header && ( line.find( "RNA" ) != std::string::npos || option[ in::file::residue_type_set ]() == "rna") &&
								silent_struct_type_ != "binary_rna" && silent_struct_type_ != "rna" ) {
			silent_struct_type_ = "rna";
			changed = true;
		} else if (( line.find( "PROTEIN_SILENTFILE" ) != std::string::npos ) || ( line.find( "PROTEIN SILENTFILE" ) != std::string::npos ) ) {
			if ( line.find("SINGLE") != std::string::npos ) {
				silent_struct_type_ ="protein_float";
			} else {
				silent_struct_type_ = "protein";
			}
			changed = true;
		}
	} else if (( line.find ("SOCREJUMP SILENTFILE" ) != std::string::npos )) {
		silent_struct_type_ = "score_jump";
		changed = true;
 }

	if ( changed ) tr.Trace << "found new silent_struct_type_ " << silent_struct_type_ << " from line " << line << std::endl;

	return changed;
}


/// @brief This is somewhat redundant if the silent file has a "REMARK RNA" line, but some older silent files didn't do that.
bool SilentFileData::check_if_rna_from_sequence_line( std::string const& line ){
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

/// @brief creates a SilentStructOP using command-line options. Just a wrapper
/// around SilentStructFactory::get_instance()->get_silent_struct_in().
SilentStructOP SilentFileData::create_SilentStructOP() {

	//static SilentStructFactory ssf;
	SilentStructOP new_ss_op;

	//if( silent_struct_type_ == "" )	new_ss_op	= ssf.get_silent_struct_in();
	//else                          	new_ss_op	= ssf.get_silent_struct( silent_struct_type_ );
	if( silent_struct_type_ == "" )	new_ss_op	= SilentStructFactory::get_instance()->get_silent_struct_in();
	else                          	new_ss_op	= SilentStructFactory::get_instance()->get_silent_struct( silent_struct_type_ );

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
	for ( iterator iter = begin(), it_end = end(); iter != it_end; ++iter )
		scores.push_back( iter->get_energy( "score" ) );

	//Real local_score_fraction( score_fraction * -1 );
	std::sort( scores.begin(), scores.end() );
	std::reverse( scores.begin(), scores.end() );
	Size const idx( static_cast< Size > ( -1 * score_fraction * scores.size() ) );
	Real const boundary( *( scores.begin() + idx ) );
	tr.Debug << "reverse_score_filter: " << std::endl;
	tr.Debug << "filtering for decoys with score worse than " << boundary << std::endl;

	Structure_Map new_structure_map_;
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
		}
	}

	structure_map_ = new_structure_map_;
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
	for ( iterator iter = begin(), it_end = end(); iter != it_end; ++iter )
		scores.push_back( iter->get_energy( "score" ) );

	std::sort( scores.begin(), scores.end() );
	Size const idx( static_cast< Size > (score_fraction * scores.size() ) );
	Real const boundary( *( scores.begin() + idx ) );
	tr.Debug << "score_filter: " << std::endl;
	tr.Debug << "filtering for decoys with score worse than " << boundary << std::endl;

	Structure_Map new_structure_map_;
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
		}
	}

	structure_map_ = new_structure_map_;
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
	typedef std::list< std::pair< Real, SilentStructOP > > 	ScoreTagList;
	ScoreTagList score_tag_list;
	for ( iterator iter = begin(), it_end = end(); iter != it_end; ++iter ) {
		Real const & silent_score = (*iter)->get_energy( "score" );
		score_tag_list.push_back( std::make_pair( silent_score, *iter /*SilentStructOP*/ ) );
	}

	score_tag_list.sort();

	Structure_Map new_structure_map_;
	Size count( 0 );
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
	}

	structure_map_ = new_structure_map_;

} // score_filter


} // namespace silent
} // namespace io
} // namespace core
