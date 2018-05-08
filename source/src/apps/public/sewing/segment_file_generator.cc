// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/segment_file_generator.cc
/// @brief a generator of sewing segment files
/// @author frankdt (frankdt@email.unc.edu)

#include <apps/public/sewing/segment_file_generator.hh>
#include <basic/Tracer.hh>
#include <algorithm>
#include <regex>

#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <boost/algorithm/string.hpp>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/dssp/Dssp.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <numeric/random/random.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <core/conformation/Atom.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

static basic::Tracer TR( "apps.pilot.frankdt.segment_file_generator" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

FileOptionKey const motif_file( "motif_file");
FileOptionKey const pdb_list_file( "pdb_list_file" );
BooleanOptionKey const strict_dssp_changes( "strict_dssp_changes" );

namespace apps {
namespace pilot {
namespace frankdt {

segment_file_generator::segment_file_generator():
	utility::pointer::ReferenceCount()
{

}

segment_file_generator::~segment_file_generator(){}

segment_file_generator::segment_file_generator( segment_file_generator const & ) {

}



segment_file_generatorOP
segment_file_generator::clone() const {
	return segment_file_generatorOP( new segment_file_generator( *this ) );
}



utility::vector1< Motif >
import_motifs( std::string motif_filename ){
	TR << "Parsing motif file" << std::endl;
	utility::vector1< Motif > motifs;
	utility::io::izstream motif_file( motif_filename );
	std::string line;
	while ( getline( motif_file, line) ) {
		utility::vector1< std::string > parsed_motif;
		boost::split( parsed_motif, line, boost::is_any_of(",") );
		Motif stored_motif;
		for ( std::string sec_struc : parsed_motif ) {
			boost::trim( sec_struc );
			utility::vector1< std::string > parsed_sec_struc;
			boost::split( parsed_sec_struc, sec_struc, boost::is_any_of("\t ") );
			if ( parsed_sec_struc.size() == 3 ) {
				SecondaryStruct sec_struct = SecondaryStruct( parsed_sec_struc[1],
					utility::string2Size( parsed_sec_struc[2] ),
					utility::string2Size( parsed_sec_struc[3] ) );
				stored_motif.abv_motif_string_ += sec_struct.dssp_;
				stored_motif.motif_string_ = stored_motif.motif_string_ +  "_" + sec_struct.dssp_+ "_" + utility::to_string( sec_struct.min_ ) + "_" + utility::to_string( sec_struct.max_ );
				stored_motif.secondary_structs_.push_back( sec_struct );
			} else {
				TR << parsed_sec_struc << std::endl;
				utility_exit_with_message( "Invalid motif file! Error: \n" + line );
			}
		}
		motifs.push_back( stored_motif );
	}
	return motifs;
}


void
invalid_motif( std::string bad_motif, std::string error_msg ){
	utility_exit_with_message( "Bad motif format: " + bad_motif + ". " + error_msg );
}

Motif
reverse_motif( Motif const motif ){
	Motif rev_motif;
	for ( SecondaryStruct sec_struct : motif.secondary_structs_ ) {
		rev_motif.secondary_structs_.insert( rev_motif.secondary_structs_.begin(), sec_struct );
	}
	for ( SecondaryStruct rev_sec_struct : rev_motif.secondary_structs_ ) {
		rev_motif.motif_string_ = rev_motif.motif_string_ +  "_" + rev_sec_struct.dssp_+ "_" + utility::to_string( rev_sec_struct.min_ ) + "_" + utility::to_string( rev_sec_struct.max_ );
	}
	std::string rev_abv_motif_string = motif.abv_motif_string_;
	std::reverse( rev_abv_motif_string.begin(), rev_abv_motif_string.end() );
	rev_motif.abv_motif_string_ = rev_abv_motif_string;
	return rev_motif;
}

bool
lt ( Motif i, Motif j ){ return i.motif_string_ < j.motif_string_; }

std::set< Motif, bool(*)( Motif, Motif ) >
check_motifs( utility::vector1< Motif > motifs ){
	bool ( *fn_ptr )( Motif, Motif ) = lt; //create function pointer
	std::set< Motif, bool(*)( Motif, Motif ) > motifs_to_use( fn_ptr ); //pass function pointer as comparor for Motif set
	for ( Motif motif : motifs ) {
		std::string motif_string = motif.abv_motif_string_;
		if ( motif_string.size() >= 1 ) {
			//check for invalid chars or two of the same char in a row
			std::regex valid_dssp("[HLEUYRN]+");
			std::regex invalid_doubles("([HLEUYRN])\1");
			if ( !regex_match( motif_string, valid_dssp ) || regex_match( motif_string, invalid_doubles ) ) {
				invalid_motif( motif_string, "Invalid characters, or two of the same dssp in a row." );
			}
			//check for terminal loops
			if ( ( motif_string[0] != 'H' && motif_string[0] != 'E' && motif_string[0] != 'Y') ||
					( motif_string[ motif_string.size() - 1 ] != 'H' && motif_string[ motif_string.size() - 1 ] != 'E' && motif_string[ motif_string.size() - 1 ] != 'Y' ) ) {
				invalid_motif( motif_string, "Terminal loop segments are not allowed." );
			}
			//if all checks pass, motif is valid.
			motifs_to_use.insert( motif );
			//if the first and last char aren't equal, add reverse motif to set, too.
			if ( motif_string[1] != motif_string[ motif_string.size() ] ) {
				motifs_to_use.insert( reverse_motif( motif ) );
			}
		} else {
			invalid_motif( motif_string, "Must have at least one segment" );
		}
	}
	return motifs_to_use;
}

void
check_for_chain_breaks( std::map< core::Size, protocols::sewing::data_storage::SmartSegmentOP > & pdb_segments ){
	//Make sure there aren't any chain breaks between residues
	TR << "Checking for chain breaks:" << std::endl;
	bool chain_break_found = false;
	for ( std::pair< core::Size, protocols::sewing::data_storage::SmartSegmentOP > current_seg : pdb_segments ) {
		for ( core::Size res_id = 1; res_id < current_seg.second->get_length(); ++res_id ) { //check for any chain breaks between residues in segment
			if ( current_seg.second->get_residue( res_id )->get_atom( 2 ).xyz().distance(
					current_seg.second->get_residue( res_id + 1 )->get_atom( 2 ).xyz() ) > 5.0 // typical caplha distance is ~4, but have seen as high as 4.3
					) {
				TR << "Found chain break inside segment." << std::endl;
				chain_break_found = true;
				break;
			}
		}
		if ( chain_break_found ) { break; }
		if ( current_seg.second->get_c_terminal_neighbor() ) { //if current segment has a neighbor, check the distance between the calpha's of segment with c-term neighboring segment
			if ( current_seg.second->get_residue( current_seg.second->get_length() )->get_atom( 2 ).xyz().distance(
					current_seg.second->get_c_terminal_neighbor()->get_residue( 1 )->get_atom( 2 ).xyz() ) > 5.0// typical caplha distance is ~4, but have seen as high as 4.3
					) {
				TR << "Found chain break between segments." << std::endl;
				chain_break_found = true;
			}
		}
		if ( chain_break_found ) { break; }
	}

	/*
	while( pdb_segments.size() > 0 && !chain_break_found && !all_segments_checked){
	if( pdb_segments.count( current_segid )  > 0 ){ // for the current segment, check for any chainbreaks between residues
	for( core::Size res_id = 1; res_id < pdb_segments[ current_segid ]->get_length(); ++res_id ){
	if( pdb_segments[ current_segid ]->get_residue( res_id )->get_atom( 2 ).xyz().distance(
	pdb_segments[ current_segid ]->get_residue( res_id + 1 )->get_atom( 2 ).xyz() ) > 5.0 // typical caplha distance is ~4, but have seen as high as 4.3
	){
	TR << "Found chain break inside segment." << std::endl;
	chain_break_found = true;
	break;
	}
	}
	if( pdb_segments[ current_segid ]->get_c_terminal_neighbor() ){ //if current segment has a neighbor, check the distance between the calpha's of segment with c-term neighboring segment
	if( pdb_segments[ current_segid ]->get_residue( pdb_segments[ current_segid ]->get_length() )->get_atom( 2 ).xyz().distance(
	pdb_segments[ current_segid ]->get_c_terminal_neighbor()->get_residue( 1 )->get_atom( 2 ).xyz() ) > 5.0// typical caplha distance is ~4, but have seen as high as 4.3
	){
	TR << "Found chain break between segments." << std::endl;
	chain_break_found = true;
	}
	}
	++current_segid;
	}
	else{
	TR << "No chain breaks found." << std::endl;
	all_segments_checked = true;
	}
	}
	*/
	if ( chain_break_found ) { // if there is a chain break, added segments are removed
		pdb_segments.clear();
	}
}



utility::vector1< Motif >
copy_set_to_vector( std::set< Motif, bool(*)( Motif, Motif ) > const copy_set ){
	utility::vector1< Motif > copy_vec;
	for ( Motif motif : copy_set ) {
		copy_vec.push_back( motif );
	}
	return copy_vec;
}

bool
dssp_code_matches( std::string motif_dssp, std::string segment_dssp ){
	if ( motif_dssp == "N" ) { // if anything is okay, don't check
		return true;
	}
	if ( motif_dssp == segment_dssp ) { //check for exact match
		return true;
	}
	if ( motif_dssp == "H" || motif_dssp == "E" || motif_dssp == "L" ) { //these matches should be exact so return false
		return false;
	}
	// check for incorrect "not" matches
	if ( motif_dssp == "U" && segment_dssp == "H" ) {
		return false;
	}
	if ( motif_dssp == "Y" && segment_dssp == "L" ) {
		return false;
	}
	if ( motif_dssp == "R" && segment_dssp == "E" ) {
		return false;
	}
	// if those are okay, it's a match.
	return true;
}



utility::vector1< Motif >
compare_segment_to_motif( utility::vector1< Motif > & motifs_to_match, core::Size position, protocols::sewing::data_storage::SmartSegmentOP current_segment ){
	utility::vector1< Motif > matched_motifs;
	utility::vector1< core::Size > motifs_dont_match;
	for ( core::Size i = motifs_to_match.size(); i > 0; --i ) {
		Motif motif = motifs_to_match[ i ];
		if ( position > motif.secondary_structs_.size() ) {
			motifs_dont_match.push_back( i );
			continue;
		}
		SecondaryStruct sec_struct = motif.secondary_structs_[ position ];
		//check for segment / motif match
		if ( dssp_code_matches( sec_struct.dssp_, utility::to_string( current_segment->get_dssp_code() ) ) && // check dssp codes
				!( current_segment->get_length() < sec_struct.min_ ) && // segment isn't too short
				!( current_segment->get_length() > sec_struct.max_ ) // segment isn't too long
				) {
			if ( position == motif.secondary_structs_.size() ) { // if it does match, check to see if it's the last segment to check in the motif
				matched_motifs.push_back( motif );
			}
		} else { // motif doesn't match the segment, so mark it to be deleted
			motifs_dont_match.push_back( i );
		}
	}
	// delete motifs that don't match the segment from the motifs_to_match
	for ( core::Size pos : motifs_dont_match ) {
		motifs_to_match.erase( motifs_to_match.begin() + ( pos - 1 ) );
	}
	return matched_motifs;
}

protocols::sewing::hashing::SegmentVectorOP
store_segment_motif_match( protocols::sewing::data_storage::SmartSegmentOP segment, core::Size last_segment_position ){
	protocols::sewing::hashing::SegmentVectorOP seg_vec = protocols::sewing::hashing::SegmentVectorOP( new protocols::sewing::hashing::SegmentVector() );
	core::Size position = 1;
	bool all_segments_added = false;
	while ( !all_segments_added ) {
		seg_vec->push_back( protocols::sewing::data_storage::SmartSegmentOP( new protocols::sewing::data_storage::SmartSegment( *segment ) ) );
		if ( position > 1 ) { //if we've already added other segments, link them up.
			protocols::sewing::data_storage::SmartSegment::link_to( seg_vec->at( position - 1 ), seg_vec->at( position ) );
		}
		if ( segment->get_c_terminal_neighbor() && position < last_segment_position ) {
			segment = segment->get_c_terminal_neighbor();
			++position;
		} else {
			if ( position != last_segment_position ) {
				utility_exit_with_message( "Matched Segment is shorter than expected. Size: " + utility::to_string( position ) + " Expected: " + utility::to_string( last_segment_position ) );
			}
			all_segments_added = true;
		}
	}
	return seg_vec;
}

void
write_segments_to_file( protocols::sewing::hashing::SegmentVectorOP seg_vec, utility::pointer::shared_ptr< utility::io::ozstream > file ){
	for ( protocols::sewing::data_storage::SmartSegmentOP out_segment : *seg_vec ) {
		*file << "SEGMENT " << out_segment->get_dssp_code();
		if ( !out_segment->is_n_terminus_fixed() ) {
			*file << " N_TERM";
		} else {
			*file << " N_FIXD";
		}
		if ( !out_segment->is_c_terminus_fixed() ) {
			*file << " C_TERM";
		} else {
			*file << " C_FIXD";
		}
		*file << std::endl;
		for ( protocols::sewing::data_storage::SmartSewingResidueOP out_residue : out_segment->get_residue_vector() ) {
			*file << "RESIDUE " << out_residue->get_amino_acid_type();
			*file << " " << out_residue->get_full_type_name();
			for ( core::Real out_chi : out_residue->get_chis() ) {
				*file << " " << out_chi;
			}
			*file << std::endl;
			for ( core::conformation::Atom out_atom : out_residue->get_atom_vector() ) {
				*file << "ATOM " << out_atom.xyz().at(0) << " " << out_atom.xyz().at(1) << " " << out_atom.xyz().at(2) << std::endl;
			}
		}
	}
	// write out here
	TR << seg_vec->size() << "Segments stored." << std::endl;
}


int
main( int argc, char * argv [] )
{
	option.add( motif_file, "File containing a description of what segment types to include including min and max length. Accepted dssp codes: HLEUYRN");
	option.add( pdb_list_file, "File containing filenames of all pdbs to be used separated by newline characters.");
	option.add( strict_dssp_changes, "Should we require secondary structure elements to be more than 1 residue long?" ).def( true );
	devel::init(argc,argv);

	if ( argc == 1 || !option[ pdb_list_file ].user() || !option[ motif_file ].user() ) {
		utility_exit_with_message("Too few args passed. A motif_file and a pdb_list_file must be provided.");
	}

	utility::io::izstream pdblist( option[ pdb_list_file ].value() );
	utility::vector1< Motif > imported_motifs = import_motifs( option[ motif_file ].value() );
	std::set< Motif, bool(*)( Motif, Motif ) > motif_set = check_motifs( imported_motifs );

	std::map<std::string, utility::pointer::shared_ptr< utility::io::ozstream > > motif_segment_files;
	for ( Motif motif : motif_set ) {
		motif_segment_files[ motif.motif_string_ ] = utility::pointer::shared_ptr< utility::io::ozstream >( new utility::io::ozstream );
		motif_segment_files[ motif.motif_string_ ]->open( "smotifs" + motif.motif_string_ + ".segments" );
		core::Size version = numeric::random::random_range(1,INT_MAX); // create a random integer within the range of possible core sizes.
		*( motif_segment_files[ motif.motif_string_ ] ) << "VERSION " << version << motif.motif_string_ << std::endl; // and write it into the segment file with the motif string
	}

	std::string line;
	protocols::sewing::hashing::SegmentVectorCOP empty_const_segvec = protocols::sewing::hashing::SegmentVectorCOP( new protocols::sewing::hashing::SegmentVector() );
	//right now generator doesn't handle ligands....
	utility::vector1< protocols::sewing::data_storage::LigandDescription > ligands;
	utility::vector1< protocols::sewing::data_storage::LigandDescription > expanded_ligands;
	std::map< core::Size, protocols::sewing::data_storage::LigandResidueCOP > partner_ligands;
	//utility::vector1< std::pair< core::Size, core::Size> > pose_segments;
	std::string pose_segment_starts_string;
	std::string pose_segment_ends_string;


	while ( getline (pdblist, line) ) {
		std::map< core::Size, protocols::sewing::data_storage::SmartSegmentOP > pdb_segments;
		TR << "Reading pose from file: " << line << std::endl;
		core::pose::PoseOP pose = core::import_pose::pose_from_file(line);
		bool strict_ss = basic::options::option[ strict_dssp_changes ].value();
		std::string pose_dssp_string;
		std::string required_resnums;
		core::select::residue_selector::ResidueSelectorCOP required_selector;
		core::Size num_added_segments = protocols::sewing::hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector( *pose, nullptr, empty_const_segvec, pdb_segments, pose_segment_starts_string, pose_segment_ends_string, pose_dssp_string,ligands, partner_ligands, expanded_ligands, required_resnums, required_selector, strict_ss );
		check_for_chain_breaks( pdb_segments );
		if ( num_added_segments < 1 || pdb_segments.size() == 0 ) {
			TR << "No segments added from pose" << std::endl;
			continue;
		}
		TR << "pdb_segments has " << pdb_segments.size() << " raw segments." << std::endl;

		//extract all proper motifs
		for ( std::pair< core::Size, protocols::sewing::data_storage::SmartSegmentOP > current_seg_pair : pdb_segments ) { // since we'll be checking every segment, checking them in order doesn't matter.
			utility::vector1< Motif > motifs_to_match = copy_set_to_vector( motif_set );
			bool viable_motif = true;
			bool segments_remain = true;
			core::Size position = 1;
			protocols::sewing::data_storage::SmartSegmentOP out_segment = current_seg_pair.second;
			//find all possible matching motifs, starting with this current segment
			while ( viable_motif && segments_remain ) {
				utility::vector1< Motif > matched_motifs = compare_segment_to_motif( motifs_to_match, position, out_segment );
				//store matched motifs - we don't need to worry about deleting these from motifs_to_match since it will be done for us
				for ( Motif motif : matched_motifs ) {
					protocols::sewing::hashing::SegmentVectorOP seg_vec_to_write = store_segment_motif_match( current_seg_pair.second, position );
					write_segments_to_file( seg_vec_to_write, motif_segment_files[ motif.motif_string_ ] );
				}
				//update variables for next loop
				viable_motif = !( motifs_to_match.size() == 0 );
				if ( out_segment->get_c_terminal_neighbor() ) {
					out_segment = out_segment->get_c_terminal_neighbor();
					++position;
				} else {
					segments_remain = false;
				}
			}
		}
		pdb_segments.clear();
	}
	return 0;
}

} //apps
} //pilot
} //frankdt

int
main( int argc, char * argv [] )
{
	apps::pilot::frankdt::main( argc, argv );
}
