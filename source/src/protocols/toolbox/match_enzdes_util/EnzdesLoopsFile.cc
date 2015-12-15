// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009


#include <protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.hh>

#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

#include <basic/Tracer.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/LexicographicalIterator.hh>

#include <sstream>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.enzdes.EnzdesLoopsFile" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

/// @details Auto-generated virtual destructor
EnzdesLoopsFile::~EnzdesLoopsFile() {}

/// @details Auto-generated virtual destructor
EnzdesLoopInfo::~EnzdesLoopInfo() {}


ResInteractions::ResInteractions() :
	targ_res_(0), num_interactions_(1),
	dis_(/* NULL */), loop_ang_(NULL), targ_ang_(NULL),
	loop_dih_(/* NULL */), targ_dih_(NULL), lt_dih_(NULL)
{
	targ_atom_names_.push_back( "CA" );
	loopres_atom_names_.push_back( "CA" );
}

ResInteractions::~ResInteractions(){}

/// @details returns false if the file fails to specify the target residue number
bool
ResInteractions::read_data( utility::io::izstream & data )
{

	utility::vector1< std::string > tokens;
	tokens.push_back("");

	bool block_end(false);

	while ( !block_end ) {

		if ( data.eof() ) {
			tr << "Error: end of file reached before RES_CONTACT_END tag was found." << std::endl;
			return false;
		}

		std::string line;
		getline( data, line );

		tokens.clear(); tokens.push_back(""); //weird utilvec1 copy behaviour makes this necessary

		tokens = utility::split( line );

		if ( tokens.size() < 1 ) continue;

		if ( tokens[1] == "RES_CONTACT_END" ) {
			block_end = true;
		} else if ( !process_input_line_tokens( tokens ) ) return false;

	} //file reading

	if ( targ_res_ == 0 ) return false;

	return true;
}

void
ResInteractions::write_data() const {

	tr << "target_res " << targ_res_ << std::endl;
	tr << "num_contacts " << num_interactions_ << std::endl;
	//tr << "distance " << dist_ << std::endl;
	//tr << "tolerance " << tolerance_ << std::endl;

	tr << "targ_atom_names ";
	for ( core::Size i = 1; i <= targ_atom_names_.size(); ++i ) tr << targ_atom_names_[i] << " ";
	tr << std::endl;

	tr << "loopres_atom_names ";
	for ( core::Size i = 1; i <= loopres_atom_names_.size(); ++i ) tr << loopres_atom_names_[i] << " ";
	tr << std::endl;

}

toolbox::match_enzdes_util::GeomSampleInfoCOP
ResInteractions::dis() const {
	return dis_; }

toolbox::match_enzdes_util::GeomSampleInfoCOP
ResInteractions::loop_ang() const {
	return loop_ang_; }

toolbox::match_enzdes_util::GeomSampleInfoCOP
ResInteractions::targ_ang() const {
	return targ_ang_; }

toolbox::match_enzdes_util::GeomSampleInfoCOP
ResInteractions::loop_dih() const {
	return loop_dih_; }

toolbox::match_enzdes_util::GeomSampleInfoCOP
ResInteractions::targ_dih() const {
	return targ_dih_; }

toolbox::match_enzdes_util::GeomSampleInfoCOP
ResInteractions::lt_dih() const {
	return lt_dih_; }

bool
ResInteractions::process_input_line_tokens(
	utility::vector1< std::string > const & tokens )
{
	using namespace toolbox::match_enzdes_util;

	bool to_return(false);
	core::Real const generic_force_K(100.0);
	core::Real const generic_period(360.0);

	if ( ( tokens.size() < 2 ) || (tokens[2] == "" ) ) {
		tr << "Error when processing res_interactions block. Line containing " << tokens[1] << " seems to have no useable data." << std::endl;

		return false;
	}

	if ( tokens[1] == "target_res" ) {
		targ_res_ = (core::Size) atoi( tokens[2].c_str() );
		to_return = true;
	} else if ( tokens[1] == "num_contacts" ) {
		num_interactions_ = (core::Size) atoi( tokens[2].c_str() );
		to_return = true;
	} else if ( tokens[1] == "targ_atom_names" ) {
		//else if( tokens[1] == "distance" ){
		// dist_ = (core::Real) atof( tokens[2].c_str() );
		// to_return = true;
		//}
		//else if( tokens[1] == "tolerance" ) {
		// tolerance_ = (core::Real) atof( tokens[2].c_str() );
		// to_return = true;
		//}
		targ_atom_names_.clear();
		for ( core::Size i = 2; i <= tokens.size(); ++i ) targ_atom_names_.push_back( tokens[ i ] );
		to_return = true;
	} else if ( tokens[1] == "targ_base_atom_names" ) {
		targ_base_atom_names_.clear();
		for ( core::Size i = 2; i <= tokens.size(); ++i ) targ_base_atom_names_.push_back( tokens[ i ] );
		to_return = true;
	} else if ( tokens[1] == "targ_base2_atom_names" ) {
		targ_base2_atom_names_.clear();
		for ( core::Size i = 2; i <= tokens.size(); ++i ) targ_base2_atom_names_.push_back( tokens[ i ] );
		to_return = true;
	} else if ( tokens[1] == "loopres_atom_names" ) {
		loopres_atom_names_.clear();
		for ( core::Size i = 2; i <= tokens.size(); ++i ) loopres_atom_names_.push_back( tokens[ i ] );
		to_return = true;
	} else if ( tokens[1] == "loopres_base_atom_names" ) {
		loopres_base_atom_names_.clear();
		for ( core::Size i = 2; i <= tokens.size(); ++i ) loopres_base_atom_names_.push_back( tokens[ i ] );
		to_return = true;
	} else if ( tokens[1] == "loopres_base2_atom_names" ) {
		loopres_base2_atom_names_.clear();
		for ( core::Size i = 2; i <= tokens.size(); ++i ) loopres_base2_atom_names_.push_back( tokens[ i ] );
		to_return = true;
	} else if ( tokens[1] == "distance" ) {
		if ( tokens.size() < 3 ) { tr << "too little information given for distance_LT." << std::endl;
			return false;
		}
		if ( tokens.size() == 3 ) dis_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), generic_force_K, 0.0 ) );
		else dis_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), 0.0 ) );
		to_return = true;
	} else if ( tokens[1] == "angle_loop" ) {
		if ( tokens.size() < 3 ) { tr << "too little information given for angle_loop." << std::endl;
			return false;
		}
		if ( tokens.size() == 3 ) loop_ang_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), generic_force_K, generic_period ) );
		else loop_ang_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), generic_period ) );
		to_return = true;
	} else if ( tokens[1] == "angle_targ" ) {
		if ( tokens.size() < 3 ) { tr << "too little information given for angle_targ." << std::endl;
			return false;
		}
		if ( tokens.size() == 3 ) targ_ang_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), generic_force_K, generic_period ) );
		else targ_ang_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), generic_period ) );
		to_return = true;
	} else if ( tokens[1] == "dih_loop" ) {
		if ( tokens.size() < 3 ) { tr << "too little information given for dih_loop." << std::endl;
			return false;
		}
		if ( tokens.size() == 3 ) loop_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), generic_force_K, generic_period ) );
		else if ( tokens.size() == 4 ) loop_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), generic_period ) );
		else loop_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), atof( tokens[5].c_str() ) ) );
		to_return = true;
	} else if ( tokens[1] == "dih_targ" ) {
		if ( tokens.size() < 3 ) { tr << "too little information given for dih_loop." << std::endl;
			return false;
		}
		if ( tokens.size() == 3 ) targ_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), generic_force_K, generic_period ) );
		else if ( tokens.size() == 4 ) targ_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), generic_period ) );
		else targ_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), atof( tokens[5].c_str() ) ) );
		to_return = true;
	} else if ( tokens[1] == "dih_LT" ) {
		if ( tokens.size() < 3 ) { tr << "too little information given for dih_loop." << std::endl;
			return false;
		}
		if ( tokens.size() == 3 ) lt_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), generic_force_K, generic_period ) );
		else if ( tokens.size() == 4 ) lt_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), generic_period ) );
		else lt_dih_ = toolbox::match_enzdes_util::GeomSampleInfoOP( new GeomSampleInfo( atof( tokens[2].c_str() ), atof( tokens[3].c_str() ), atof( tokens[4].c_str() ), atof( tokens[5].c_str() ) ) );
		to_return = true;
	}


	if ( !to_return ) tr << "Unspecified error when processing res_interactions block line containing " << tokens[1] << std::endl;

	return to_return;
}

CstResInteractions::CstResInteractions() :
	ResInteractions(),
	resA_(false), cst_block_(0)
{}


bool
CstResInteractions::read_data( utility::io::izstream & data )
{

	utility::vector1< std::string > tokens;
	tokens.push_back("");

	bool cst_line_encountered(false), block_end(false);

	while ( !block_end ) {

		//std::cerr << "reading CstResInp, token[1] is: '" << tokens[1] << "'" << std::endl;
		if ( data.eof() ) {
			tr << "Error: end of file reached before CST_TARGET_END tag was found." << std::endl;
			this->write_data();
			return false;
		}

		std::string line;
		getline( data, line );

		//std::cerr << "reading CstResInteractions, line: " << line << std::endl;
		tokens.clear(); tokens.push_back(""); //weird utilvec1 copy behaviour makes this necessary

		tokens = utility::split( line );

		if ( tokens.size() < 1 ) continue;

		if ( tokens[1] == "CST_TARGET_END" ) {
			block_end = true;
		} else if ( tokens[1] == "cst_target" ) {

			cst_block_ = (core::Size ) atoi( tokens[2].c_str() );

			if ( tokens[3] == "A" ) resA_ = true;
			else if ( tokens[3] == "B" ) resA_ = false;
			else {
				tr << "When specifing a cst_target line, the 3rd element has to be either A or B, corresponding to which of the two residues in the cstfile is desired." << std::endl;
				return false;
			}

			cst_line_encountered = true;

		} else if ( !process_input_line_tokens( tokens ) ) return false; //if cst_target

	} //file reading
	//std::cerr << "reading CstResInp after while, token[1] is: '" << tokens[1] << "'" << std::endl;

	if ( !cst_line_encountered ) {
		tr << "Did not find 'cst_target' line in CST block." << std::endl;
		return false;
	}

	if ( this->targ_res() != 0 ) {

		tr << "An explicit residue number has been specified in a CST_TARGET block" << std::endl;
		return false;
	}

	return true;
}

void
CstResInteractions::write_data() const {

	ResInteractions::write_data();

	tr << "cst_target " << cst_block_ << " ";

	if ( resA_ ) tr << "A";
	else tr << "B";

	tr << std::endl;
}

EnzdesLoopInfo::EnzdesLoopInfo() :
	loop_start_(0), loop_end_(0),
	loop_start_pdb_(0), loop_end_pdb_(0),
	loop_start_pdb_chain_(' '), loop_end_pdb_chain_(' '),
	pose_numb_(false), pdb_numb_(false),
	min_length_(0), max_length_(0),
	preserve_buried_contacts_(false), contact_buried_problematic_res_(false)
{
	ss_strings_.clear();
	res_interactions_.clear();
	cstres_interactions_.clear();
}

bool
EnzdesLoopInfo::read_loops_file_block(
	utility::io::izstream & data
){

	ss_strings_.clear();
	res_interactions_.clear();
	cstres_interactions_.clear();
	core::Size linenum(0);

	utility::vector1< std::string > tokens;
	tokens.push_back("");

	bool loop_end( false );

	bool min_length_tag_found(false), max_length_tag_found(false);

	utility::vector1< std::string > ss_blueprints;

	while ( !loop_end ) {

		if ( data.eof() ) {
			tr << "Error: end of file reached before LOOP_END tag was found." << std::endl;
			return false;
		}
		linenum++;

		std::string line;
		getline( data, line );

		tokens.clear(); tokens.push_back(""); //weird utilvec1 copy behaviour makes this necessary

		tokens = utility::split( line );

		if ( tokens.size() < 1 ) continue;

		//std::cerr << "eli line is '" << line << "'  , tokens[1] is '" << tokens[1] << "'" << std::endl;

		if ( tokens[1] == "LOOP_END" ) loop_end = true;

		else if ( tokens[1] == "start" ) {

			loop_start_ = (core::Size) atoi( tokens[2].c_str() );
			pose_numb_ = true;

		} else if ( tokens[1] == "stop" ) {

			loop_end_ = (core::Size) atoi( tokens[2].c_str() );
			pose_numb_ = true;

		} else if ( tokens[1] == "pdb_start" ) {

			loop_start_pdb_chain_ = tokens[2][0];
			loop_start_pdb_ = (core::Size) atoi( tokens[3].c_str() );
			pdb_numb_ = true;

		} else if ( tokens[1] == "pdb_stop" ) {

			loop_end_pdb_chain_ = tokens[2][0];
			loop_end_pdb_ = (core::Size) atoi( tokens[3].c_str() );
			pdb_numb_ = true;

		} else if ( tokens[1] == "min_length" ) {
			min_length_ = (core::Size) atoi( tokens[2].c_str() );
			min_length_tag_found = true;
		} else if ( tokens[1] == "max_length" ) {
			max_length_ = (core::Size) atoi( tokens[2].c_str() );
			max_length_tag_found = true;
		} else if ( tokens[1] == "ss_string" ) {
			for ( core::Size i = 2; i <= tokens.size(); ++i ) ss_strings_.push_back( tokens[ i ] );
		} else if ( tokens[1] == "ss_blueprint" ) {
			for ( core::Size i = 2; i <= tokens.size(); ++i ) ss_blueprints.push_back( tokens[i] );
		} else if ( tokens[1] == "CST_TARGET_BEGIN" ) {

			//std::cerr << "detected CST_TARG block on line " << linenum << ", tokens[1] is " << tokens[1] << std::endl;
			cstres_interactions_.push_back( CstResInteractions() );

			if ( ! cstres_interactions_[ cstres_interactions_.size() ].read_data( data ) ) {
				tr << "Error occured when processing a CST_TARGET_BEGIN block." << std::endl;
				return false;
			}
			//std::cerr << "successfully read CST_TARGET block " << std::endl;
		} else if ( tokens[1] == "RES_CONTACT_BEGIN" ) {

			res_interactions_.push_back( ResInteractions() );

			if ( ! res_interactions_[ res_interactions_.size() ].read_data( data ) ) {
				tr << "Error occured when processing a RES_CONTACT_BEGIN block." << std::endl;
				return false;
			}
		} else if ( tokens[1] == "instruction" ) {
			//unimplemented for now
		}


	} //file reading loop

	tr << res_interactions_.size() << " residue contact blocks and " << cstres_interactions_.size() << " cst residue contact blocks were read in. " << std::endl;

	for ( core::Size i = 1; i <= ss_blueprints.size(); ++i ) generate_ss_strings_from_blueprint( ss_blueprints[i] );

	//in case the user forgot these tags, we take this to mean that the loop length isn't changing
	if ( !min_length_tag_found ) {
		min_length_ = loop_end_ - loop_start_ + 1;
		tr << "No min_length tag found, min_length_ set to " << min_length_ << "." << std::endl;
	}
	if ( !max_length_tag_found ) {
		max_length_ = loop_end_ - loop_start_ + 1;
		tr << "No max_length tag found, max_length_ set to " << max_length_ << "." << std::endl;
	}

	return check_data_consistency( true );

} //read_loops_file_block

bool
EnzdesLoopInfo::check_data_consistency(
	bool report
) const
{
	if ( loop_start_ == 0 && loop_start_pdb_ == 0 ) {
		if ( report ) tr << "file didn't seem to specify loop start" << std::endl;
		return false;
	}

	if ( loop_end_ == 0 && loop_start_pdb_ == 0 ) {
		if ( report ) tr << "file didn't seem to specify loop stop" << std::endl;
		return false;
	}

	if ( pose_numb_ && pdb_numb_ ) {
		if ( report ) tr << "specify loop start and stop by PDB or pose numbering" << std::endl;
		return false;
	}

	if ( min_length_ == 0 ) {
		if ( report ) tr << "file specified illegal minimum loop length (min_length tag) of 0" << std::endl;
		return false;
	}

	if ( max_length_ == 0 ) {
		if ( report ) tr << "file specified illegal maximum loop length (max_length tag) of 0" << std::endl;
		return false;
	}

	if ( loop_end_ <= loop_start_ && loop_end_pdb_ <= loop_start_pdb_ ) {
		if ( report ) tr << "specified start of loop is after stop of loop. go work at mcdonalds." << std::endl;
		return false;
	}

	if ( min_length_ > max_length_ ) {
		tr << "min_length of loop is larger than max_length of loop. go work at mcdonalds" << std::endl;
		return false;
	}

	if ( min_length_ < 4 ) {
		if ( report ) tr << "min length of loop is too short, has to be at least 4." << std::endl;
		return false;
	}

	for ( core::Size i = 1; i <= ss_strings_.size(); ++i ) {

		if ( ss_strings_[i].length() < min_length_ ) tr << "WARNING: secondary structure string " << i << ", " << ss_strings_[i] << ", specified in the enzdes loops file is shorter than min_length for the loop. This overrides the specified min_length." << std::endl;

		if ( ss_strings_[i].length() > max_length_ ) tr << "WARNING: secondary structure string " << i <<  ", " << ss_strings_[i] << ", specified in the enzdes loops file is longer than max_length for the loop. This overrides the specified max_length." << std::endl;
	}

	return true;
} //check_data_consistency


void
EnzdesLoopInfo::generate_ss_strings_from_blueprint(
	std::string const & ss_blueprint )
{

	utility::vector1< std::string > blueprint_elements = utility::string_split( ss_blueprint, ')' );

	utility::vector1< char > blueprint_element_ss_chars;
	utility::vector1< core::Size > blueprint_element_num_lengths;
	utility::vector1< core::Size > blueprint_element_min_lengths;
	core::Size num_previous_ss_strings( ss_strings_.size() );
	core::Size num_combos(1);

	for ( utility::vector1< std::string >::iterator ble_it = blueprint_elements.begin(); ble_it != blueprint_elements.end(); ++ble_it ) {

		if ( *ble_it == "" ) continue;

		std::string ble_length_info = ble_it->substr(2);
		utility::vector1< std::string > min_max_strings( utility::string_split( ble_length_info, '-' ) );

		if ( ( (min_max_strings.size() != 1 ) && (min_max_strings.size() != 2 ) ) || ( (*ble_it)[1] != '(' ) ) {
			utility_exit_with_message("SS_blueprint "+ *ble_it + " could not be understood when trying to generate ss_strings from it.");
		}

		blueprint_element_ss_chars.push_back( (*ble_it)[0] );

		if ( min_max_strings.size() == 1 ) {
			blueprint_element_num_lengths.push_back(1);
			blueprint_element_min_lengths.push_back( (core::Size) atoi( min_max_strings[1].c_str() ) );
		} else {  //implies min_max_strings.size() == 2, see above check
			core::Size ble_min_length( (core::Size) atoi( min_max_strings[1].c_str() ) );
			core::Size ble_max_length( (core::Size) atoi(min_max_strings[2].c_str() ) );
			//runtime_assert( ble_min_length <= ble_max_length );
			blueprint_element_num_lengths.push_back( ble_max_length - ble_min_length + 1);
			blueprint_element_min_lengths.push_back( ble_min_length );
		}

		num_combos = num_combos * blueprint_element_num_lengths[ blueprint_element_num_lengths.size() ];
	}
	core::Size num_elements( blueprint_element_num_lengths.size() );
	//using LexicographicalIterator to go through all the length combinations
	core::Size too_long_strings(0), too_short_strings(0);

	utility::LexicographicalIterator lex( blueprint_element_num_lengths );
	std::set< std::string > observed_ss_strings; //safeguard to prevent the same string appearing more than once
	core::Size redundant_strings(0);
	while ( !lex.at_end() ) {

		core::Size length_this_string(0);

		for ( core::Size i = 1; i <= num_elements; ++i ) {
			length_this_string += lex[i]+blueprint_element_min_lengths[i] - 1;
		}
		if ( length_this_string < min_length_ ) {
			too_short_strings++;
			++lex;
			continue;
		}
		if ( length_this_string > max_length_ ) {
			too_long_strings++;
			++lex;
			continue;
		}

		//ok, this particular string is within the desired range. now let's build it
		std::string ss_string("");
		for ( core::Size i = 1; i <= num_elements; ++i ) {

			core::Size length_this_element( lex[i]+blueprint_element_min_lengths[i] - 1 );
			for ( core::Size j = 1; j <= length_this_element; ++j ) {
				ss_string.push_back( blueprint_element_ss_chars[i] );
			}
		}
		//runtime_assert( ss_string.size() == length_this_string );
		if ( observed_ss_strings.find( ss_string ) == observed_ss_strings.end() ) {
			ss_strings_.push_back( ss_string );
			observed_ss_strings.insert( ss_string );
		} else redundant_strings++;

		++lex;
	}

	tr << "SS blueprint string " << ss_blueprint << " led to the following " << ss_strings_.size() - num_previous_ss_strings << " secondary structure strings out of a total of " << num_combos << " possible combinations: " << std::endl;

	for ( core::Size i = num_previous_ss_strings + 1; i <= ss_strings_.size(); ++i ) tr << ss_strings_[i] << ", ";

	tr << std::endl << "A total of " << too_long_strings << " ss_strings were ignored because they were too long, a total of " << too_short_strings << " ss_strings were ignored because they were too short, and a total of " << redundant_strings << " ss_strings were ignored because they were redundant." << std::endl;

}


EnzdesLoopsFile::EnzdesLoopsFile()
: file_read_(false )
{
	enzloops_.clear();
}

bool
EnzdesLoopsFile::read_loops_file(
	std::string filename )
{

	utility::io::izstream data( filename.c_str() );
	std::istringstream line_stream;

	enzloops_.clear();

	if ( !data ) {
		std::cerr << "ERROR:: Unable to open enzdes loops file: "
			<< filename << std::endl;
		std::exit( 1 );
	}
	tr.Info << "reading enzdes loops from  " << filename << " ..." << std::endl;

	core::Size counted_loops(0);

	while ( !data.eof() ) {

		std::string key, line;
		getline(data,line);
		line_stream.clear();
		line_stream.str(line);
		line_stream >> key;

		if ( key == "LOOP_BEGIN" ) {

			counted_loops++;
			tr << "reading loop block " << counted_loops << "... " << std::endl;

			EnzdesLoopInfoOP el( new EnzdesLoopInfo() );

			if ( el->read_loops_file_block( data ) ) {

				enzloops_.push_back( el );

				tr << "Data read: start(" << el->start() << "), stop(" << el->stop() << "), min_length(" << el->min_length() << "),";
				tr << " max_length(" << el->max_length() << ") " << std::endl;
				tr << el->ss_strings().size() << " secstruct strings, " << el->res_interactions().size() << " res interactions" << std::endl;
				tr << " ... done reading block " << counted_loops << "." << std::endl;
			} else {
				tr << "Error when reading file " << filename << ". Block " << counted_loops << " or the info therein was corrupted." << std::endl;
				return false;
			}
		}
	}

	file_read_ = true;

	return true;

} //read_loops_file


EnzdesLoopInfoCOP
EnzdesLoopsFile::loop_info( core::Size l ) const
{
	assert( l <= enzloops_.size() );

	return enzloops_[l];

}

void
EnzdesLoopsFile::clear()
{
	file_read_ = false;
	enzloops_.clear();
}

} //match_enzdes_util
} //toolbox
} //protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::EnzdesLoopInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( loop_start_ ) ); // core::Size
	arc( CEREAL_NVP( loop_end_ ) ); // core::Size
	arc( CEREAL_NVP( loop_start_pdb_ ) ); // core::Size
	arc( CEREAL_NVP( loop_end_pdb_ ) ); // core::Size
	arc( CEREAL_NVP( loop_start_pdb_chain_ ) ); // char
	arc( CEREAL_NVP( loop_end_pdb_chain_ ) ); // char
	arc( CEREAL_NVP( pose_numb_ ) ); // _Bool
	arc( CEREAL_NVP( pdb_numb_ ) ); // _Bool
	arc( CEREAL_NVP( min_length_ ) ); // core::Size
	arc( CEREAL_NVP( max_length_ ) ); // core::Size
	arc( CEREAL_NVP( ss_strings_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( preserve_buried_contacts_ ) ); // _Bool
	arc( CEREAL_NVP( contact_buried_problematic_res_ ) ); // _Bool
	arc( CEREAL_NVP( res_interactions_ ) ); // utility::vector1<ResInteractions>
	arc( CEREAL_NVP( cstres_interactions_ ) ); // utility::vector1<CstResInteractions>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::EnzdesLoopInfo::load( Archive & arc ) {
	arc( loop_start_ ); // core::Size
	arc( loop_end_ ); // core::Size
	arc( loop_start_pdb_ ); // core::Size
	arc( loop_end_pdb_ ); // core::Size
	arc( loop_start_pdb_chain_ ); // char
	arc( loop_end_pdb_chain_ ); // char
	arc( pose_numb_ ); // _Bool
	arc( pdb_numb_ ); // _Bool
	arc( min_length_ ); // core::Size
	arc( max_length_ ); // core::Size
	arc( ss_strings_ ); // utility::vector1<std::string>
	arc( preserve_buried_contacts_ ); // _Bool
	arc( contact_buried_problematic_res_ ); // _Bool
	arc( res_interactions_ ); // utility::vector1<ResInteractions>
	arc( cstres_interactions_ ); // utility::vector1<CstResInteractions>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::EnzdesLoopInfo );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::EnzdesLoopInfo )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::EnzdesLoopsFile::save( Archive & arc ) const {
	arc( CEREAL_NVP( file_read_ ) ); // _Bool
	arc( CEREAL_NVP( enzloops_ ) ); // utility::vector1<EnzdesLoopInfoOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::EnzdesLoopsFile::load( Archive & arc ) {
	arc( file_read_ ); // _Bool
	arc( enzloops_ ); // utility::vector1<EnzdesLoopInfoOP>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::EnzdesLoopsFile );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::EnzdesLoopsFile )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::CstResInteractions::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::toolbox::match_enzdes_util::ResInteractions >( this ) );
	arc( CEREAL_NVP( resA_ ) ); // _Bool
	arc( CEREAL_NVP( cst_block_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::CstResInteractions::load( Archive & arc ) {
	arc( cereal::base_class< protocols::toolbox::match_enzdes_util::ResInteractions >( this ) );
	arc( resA_ ); // _Bool
	arc( cst_block_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::CstResInteractions );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::CstResInteractions )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::ResInteractions::save( Archive & arc ) const {
	arc( CEREAL_NVP( targ_res_ ) ); // core::Size
	arc( CEREAL_NVP( targ_atom_names_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( targ_base_atom_names_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( targ_base2_atom_names_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( num_interactions_ ) ); // core::Size
	arc( CEREAL_NVP( dis_ ) ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( CEREAL_NVP( loop_ang_ ) ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( CEREAL_NVP( targ_ang_ ) ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( CEREAL_NVP( loop_dih_ ) ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( CEREAL_NVP( targ_dih_ ) ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( CEREAL_NVP( lt_dih_ ) ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( CEREAL_NVP( loopres_atom_names_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( loopres_base_atom_names_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( loopres_base2_atom_names_ ) ); // utility::vector1<std::string>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::ResInteractions::load( Archive & arc ) {
	arc( targ_res_ ); // core::Size
	arc( targ_atom_names_ ); // utility::vector1<std::string>
	arc( targ_base_atom_names_ ); // utility::vector1<std::string>
	arc( targ_base2_atom_names_ ); // utility::vector1<std::string>
	arc( num_interactions_ ); // core::Size
	arc( dis_ ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( loop_ang_ ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( targ_ang_ ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( loop_dih_ ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( targ_dih_ ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( lt_dih_ ); // toolbox::match_enzdes_util::GeomSampleInfoOP
	arc( loopres_atom_names_ ); // utility::vector1<std::string>
	arc( loopres_base_atom_names_ ); // utility::vector1<std::string>
	arc( loopres_base2_atom_names_ ); // utility::vector1<std::string>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::ResInteractions );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::ResInteractions )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_match_enzdes_util_EnzdesLoopsFile )
#endif // SERIALIZATION
