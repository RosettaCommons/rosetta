// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/BasisMapGenerator.cc
/// @brief Given a model file, edge file ,and one or more input structures, generates alignment files for use with AppendAssemblyMover.
/// @author guffysl (guffy@email.unc.edu)
/// @author minniel (minnie@email.unc.edu)

#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/hashing/BasisMapGenerator.hh>
#include <protocols/sewing/hashing/AlignmentGenerator.hh>
#include <protocols/sewing/hashing/Hasher.hh>
#include <protocols/sewing/hashing/EdgeMapGenerator.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>
//Temporary, functions will move later
//#include <protocols/sewing/util/io.hh >

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR( "protocols.sewing.BasisMapGenerator" );


namespace protocols {
namespace sewing {
namespace hashing {

BasisMapGenerator::BasisMapGenerator()
{
	// initialize_from_options();
	// alignment_settings_from_options();
}


BasisMapGenerator::BasisMapGenerator(
	std::string model_file_name,
	std::string edge_file_name,
	//std::string alignment_file_name,
	utility::vector1< core::Size > match_segments,
	utility::vector1< core::Size > pose_segment_starts,
	utility::vector1< core::Size > pose_segment_ends,
	core::Size recursive_depth )
{
	// alignment_settings_ = alignment_settings;
	//alignment_file_name_ = alignment_file_name;
	set_private_data( model_file_name, edge_file_name, match_segments, pose_segment_starts, pose_segment_ends, recursive_depth );

}

BasisMapGenerator::BasisMapGenerator(
	EdgeMapGeneratorOP edge_file_reader,
	std::string model_file_name,
	std::string alignment_file_name )
{
	basis_map_ = utility::pointer::shared_ptr< BasisMap > ( new BasisMap() );
	set_model_file( model_file_name );
	set_edge_file_reader( edge_file_reader );
	read_alignments_from_file( alignment_file_name );
	alignment_settings_.recursive_depth_ = 1; //override the recursive depth used to create the alignment file, we'll only need to calculate one layer out.
}

BasisMapGenerator::BasisMapGenerator(
	EdgeMapGeneratorOP edge_file_reader,
	std::string model_file_name )
{
	basis_map_ = BasisMapOP ( new BasisMap );
	// alignment_settings_from_options(); // initialize to (placeholder) default values. Won't be used since there will be no pose starting structure.
	alignment_settings_.recursive_depth_ = 1; //since we aren't creating an alignment file, we only need to calculate one layer out.
	set_model_file( model_file_name );
	set_edge_file_reader( edge_file_reader );
}

BasisMapGenerator::BasisMapGenerator(
	EdgeMapGeneratorOP edge_file_reader,
	SegmentVectorCOP segments )
{
	basis_map_ = BasisMapOP ( new BasisMap );
	// alignment_settings_from_options(); // initialize to (placeholder) default values. Won't be used since there will be no pose starting structure.
	alignment_settings_.recursive_depth_ = 1; //since we aren't creating an alignment file, we only need to calculate one layer out.
	segment_vector_ = segments;
	version_ = segment_vector_->get_version();
	set_edge_file_reader( edge_file_reader );
}




void
BasisMapGenerator::edge_file_from_options(){
	if ( !basic::options::option[ basic::options::OptionKeys::sewing::edge_file_name ].user() ) {
		utility_exit_with_message("Error: Please provide a SEWING edge file to use for alignment generation." );
	}
	std::string edge_file_name = basic::options::option[ basic::options::OptionKeys::sewing::edge_file_name ].value();


	set_edge_file(edge_file_name);
}

void
BasisMapGenerator::model_file_from_options(){
	if ( !basic::options::option[ basic::options::OptionKeys::sewing::model_file_name ].user() ) {
		utility_exit_with_message("Error: Please provide a SEWING model file to use for alignment generation." );
	}
	std::string model_file_name = basic::options::option[ basic::options::OptionKeys::sewing::model_file_name ].value();
	set_model_file( model_file_name );
}

void
BasisMapGenerator::alignment_settings_from_options(){
	utility::vector1< core::Size > match_segments = basic::options::option[ basic::options::OptionKeys::sewing::match_segments ].value();
	utility::vector1< core::Size > pose_segment_starts = basic::options::option[ basic::options::OptionKeys::sewing::pose_segment_starts ].value();
	utility::vector1< core::Size > pose_segment_ends = basic::options::option[ basic::options::OptionKeys::sewing::pose_segment_ends ].value();
	core::Size recursive_depth = basic::options::option[ basic::options::OptionKeys::sewing::recursive_depth ].value();
	set_alignment_settings( match_segments, pose_segment_starts, pose_segment_ends, recursive_depth );
}

void
BasisMapGenerator::set_model_file( std::string model_file_name ){
	//Read in model file
	segment_vector_ = ModelFileReader::read_model_file( model_file_name );
	//Assign data from model file
	version_ = segment_vector_->get_version();

}

void
BasisMapGenerator::set_edge_file( std::string edge_file_name ){
	EdgeMapGeneratorOP edge_file_reader = EdgeMapGeneratorOP( new EdgeMapGenerator( edge_file_name ) );
	edge_file_reader->segment_vector( segment_vector_ );
	set_edge_file_reader( edge_file_reader );

}


//Edge file readers no longer exist!! Now we'll need an EdgeMapGenerator
void
BasisMapGenerator::set_edge_file_reader( EdgeMapGeneratorOP edge_file_reader ){
	edge_file_reader_ = edge_file_reader;
	edge_file_reader_->read_edge_file();
	//Check that the versions are the same. If not, throw an error and exit.
	if ( version_ != edge_file_reader_->version() ) {
		utility_exit_with_message( " Error! The provided edge file does not correspond to the provided model file. Please ensure that your model file and edge file have matching versions. " );
	}
	//Assign data from edge file
	hasher_settings_ = edge_file_reader_->hasher_settings();
	alignment_generator_ = AlignmentGeneratorOP( new AlignmentGenerator( hasher_settings_, segment_vector_, pdb_segments_ ) );
}

void
BasisMapGenerator::set_alignment_settings(
	utility::vector1< core::Size > match_segments,
	utility::vector1< core::Size > pose_segment_starts,
	utility::vector1< core::Size > pose_segment_ends,
	core::Size recursive_depth){

	//Make sure the number of starts and ends are the same
	if ( pose_segment_starts.size() != pose_segment_ends.size() ) {
		utility_exit_with_message( "You must supply the same number of pose segment starts and ends!" );
	}
	//Construct AlignmentSettings object from arguments

	//The two core::Size are the numbers of the first and last residue in each segment
	utility::vector1< std::pair< core::Size, core::Size > > segments;
	segments.clear();
	for ( core::Size i = 1; i <= pose_segment_starts.size(); ++i ) {
		segments.push_back( std::make_pair( pose_segment_starts[ i ], pose_segment_ends[ i ] ) );
	}
	alignment_settings_ = AlignmentSettings( match_segments, segments, recursive_depth );
}

//NOTE: Add recursive_depth and edge_file_name and alignment_file_name as options
//(Alignment_file_name will actually be a prefix for now--later add in functionality to give each one a unique name with JobQueen system)
void BasisMapGenerator::initialize_from_options() {
	//  //Get model file name
	model_file_from_options();
	//  //Get edge file name
	edge_file_from_options();
	//  //Get options required for AlignmentSettings
	alignment_settings_from_options();
}

void BasisMapGenerator::set_private_data( std::string model_file_name,
	std::string edge_file_name,
	utility::vector1< core::Size > match_segments,
	utility::vector1< core::Size > pose_segment_starts,
	utility::vector1< core::Size > pose_segment_ends,
	core::Size recursive_depth ){

	// Store Alignment Settings
	set_alignment_settings( match_segments, pose_segment_starts, pose_segment_ends, recursive_depth );

	//Read in model file
	set_model_file( model_file_name );
	//Read in edge file
	set_edge_file( edge_file_name );
	// //Create BasisMap;
	basis_map_ = utility::pointer::shared_ptr< BasisMap > ( new BasisMap() );
}


BasisMapGenerator::~BasisMapGenerator(){}

BasisMapGenerator::BasisMapGenerator( BasisMapGenerator const & src )
{
	hasher_settings_ = src.hasher_settings();
	version_ = src.version();
	segment_vector_ = src.segment_vector(); //This is safe now since it's a COP
	pdb_segments_ = src.const_pdb_segments();

	// edge_map_ = src.edge_map();
	edge_file_reader_ = src.edge_file_reader();
	alignment_settings_ = src.alignment_settings();
	alignment_generator_ = src.alignment_generator();
	basis_map_ = src.basis_map();
}

BasisMapGeneratorOP
BasisMapGenerator::clone() const{
	return BasisMapGeneratorOP( new BasisMapGenerator( *this ) );
}


void
BasisMapGenerator::write_alignments_to_file(std::string new_file_name) const{
	//All we really need for the file name is the PDB info and the version
	utility::io::ozstream file;
	file.open_append( new_file_name );
	//First write version
	file << edge_file_reader_->version() << std::endl;
	//Then write alignment settings
	//We need to convert the alignment settings to strings
	for ( core::Size i = 1; i <= alignment_settings_.segments_.size(); ++i ) {
		file << utility::to_string( alignment_settings_.segments_[ i ].first )  << " "
			<< utility::to_string( alignment_settings_.segments_[ i ].second );
		if ( i != alignment_settings_.segments_.size() ) {
			file << " ";
		}
	}
	//Then match segments on a new line
	file << std::endl;
	utility::vector1< core::Size >::const_iterator end_iter = alignment_settings_.match_segments_.end();
	--end_iter; //End on the second to last element for spacing purposes
	for ( utility::vector1< core::Size >::const_iterator iter = alignment_settings_.match_segments_.begin();
			iter != end_iter;
			++iter ) {
		file << utility::to_string( *iter ) << " ";
	}
	file << utility::to_string( *end_iter ) << std::endl;
	//Finally recursive depth
	file << utility::to_string( alignment_settings_.recursive_depth_ ) << std::endl;
	//Now we'll write the actual alignments as efficiently as we can.

	//Iterate through BasisMap
	BasisMap::const_iterator bmiterator = basis_map_->begin();
	for ( ; bmiterator != basis_map_->end(); ++bmiterator ) {
		//Tag outer map S1
		file << "S1 " << utility::to_string( bmiterator->first ) << std::endl;
		std::map< core::Size, utility::vector1< std::pair< core::Size, core::Size > > >::const_iterator inner_map_iter = bmiterator->second.begin();
		//Iterate through inner map of BasisMap
		for ( ; inner_map_iter != bmiterator->second.end(); ++inner_map_iter ) {
			//Tag inner map S2
			file << "S2 " << utility::to_string( inner_map_iter->first ) << std::endl;
			utility::vector1< std::pair< core::Size, core::Size > >::const_iterator pair_iter = inner_map_iter->second.begin();
			//Then just list basis pairs, one per line
			//Iterate through vector of pairs
			for ( ; pair_iter != inner_map_iter->second.end(); ++pair_iter ) {
				file << utility::to_string( pair_iter->first )  << " "
					<< utility::to_string( pair_iter->second ) << std::endl;
			}

		}
	}
	file.close();
}


void
BasisMapGenerator::read_alignments_from_file(std::string alignment_file_name) {
	//attempt to open alignment file
	utility::io::izstream alignment_file( alignment_file_name );
	if ( !alignment_file.good() ) {
		utility_exit_with_message( "Could not find file with name " + alignment_file_name );
	}

	std::string line;

	//The first line is the version
	alignment_file.getline( line );
	version_ = utility::string2Size( line );

	//The second line contains the AlignmentSettings.segments_ data
	alignment_file.getline( line );
	utility::vector1< std::string > segment_line = utility::string_split( line );

	utility::vector1< std::pair< core::Size, core::Size > > segments;
	segments.clear();

	//Get values for segments_
	for ( core::Size i = 1; i <= segment_line.size() - 1; ++i ) {
		segments.push_back( std::make_pair( utility::string2Size( segment_line[ i ] ), utility::string2Size( segment_line[ i+1 ] ) ) );
	}

	//The third line contains the AlignmentSettings.match_segments_ data
	alignment_file.getline( line );
	utility::vector1< std::string > match_segment_line = utility::string_split( line );

	utility::vector1< core::Size > match_segments;
	match_segments.clear();

	//Get values for match_segments_
	for ( core::Size i = 1; i <= match_segment_line.size(); ++i ) {
		match_segments.push_back(  utility::string2Size( match_segment_line[ i ] ) );
	}

	//The fourth line contains the AlignmentSettings.recursive_depth_
	alignment_file.getline( line );
	core::Size recursive_depth = utility::string2Size( line );

	//initialize AllignmentSettings data member from values
	alignment_settings_ = AlignmentSettings( match_segments, segments, recursive_depth );

	//The rest of the file contains basis pairs
	utility::vector1< std::string > current_line;
	//int mod1;
	core::Size current_outer_key;
	core::Size current_inner_key;

	while ( alignment_file.getline( line ) ) {
		current_line = utility::string_split( line );
		if ( current_line[ 1 ].compare( "S1" ) ) { // line stores outer key segment ( segment_id)
			//   mod1 = utility::string2int( current_line[ 2 ] );
			//      seg1 = utility::string2Size( current_line[ 2 ] );
			current_outer_key = utility::string2Size( current_line[ 2 ] );
		} else if ( current_line[ 1 ].compare( "S2" ) ) { // line stores inner key segment (model_id, segment_id)
			//mod1 = utility::string2int( current_line[ 2 ] );
			//seg1 = utility::string2Size( current_line[ 3 ] );
			current_inner_key = utility::string2Size( current_line[ 2 ] );
		} else { //line stores the basis pair info
			for ( core::Size i = 1; i <= current_line.size() - 1; ++i ) { //read in pairs & store in basis map
				( * basis_map_ )[ current_outer_key ][ current_inner_key ].push_back( std::make_pair( utility::string2Size( current_line[ i ] ), utility::string2Size( current_line[ i+1 ] ) ) );
			}
		}
	}
	//Now we've read in all the basis pairs
	alignment_file.close();
}

BasisMap::iterator
BasisMapGenerator::get_alignments( data_storage::SmartSegmentOP query_segment )
{
	core::Size basis_map_key =  query_segment->get_segment_id();
	if ( basis_map_->find( basis_map_key ) == basis_map_->end() ) { //if alignments haven't been stored already, add them to the basis_map_
		std::set< data_storage::SmartSegmentCOP > query_segments;
		query_segments.insert( query_segment );
		core::Size start_depth = 0;
		recurse_over_segments( start_depth, query_segments );
	}
	if ( basis_map_->find( basis_map_key ) == basis_map_->end() ) { //if alignments still haven't been stored -- THIS SHOULDN'T HAPPEN!
		utility_exit_with_message("Error: recurse_over_segments did not add alignments. How did this happen?" );
	}
	return basis_map_->find( basis_map_key );
}

void
BasisMapGenerator::delete_segment_basis_pairs( core::Size segment_to_delete )
{
	BasisMap::iterator seg_it = basis_map_->find( segment_to_delete );
	if ( seg_it != basis_map_->end() ) { //check to make sure segment is in the basis map
		basis_map_->erase( seg_it ); //delete it from the map
	}
}



//NOTE: This function now assumes that all segments are in the EdgeMap (new behavior!!)
void
BasisMapGenerator::recurse_over_segments( core::Size current_depth, std::set< data_storage::SmartSegmentCOP > segments_to_check )
{
	//Check if we've finished recursing
	if ( current_depth >= alignment_settings_.recursive_depth_ ) {
		return;
	}

	++current_depth;
	std::set< data_storage::SmartSegmentCOP > next_set_to_check;
	next_set_to_check.clear();

	std::set< data_storage::SmartSegmentCOP >::iterator seg_it = segments_to_check.begin();
	for ( ; seg_it != segments_to_check.end(); ++seg_it ) {
		//Check this segment vs all segments it has edges for
		core::Size current_segment = ( *seg_it )->get_segment_id();
		//POSSIBLE CASES
		//Case 1: Segment is already a key in the BasisMap-> do nothing as we've clearly already looked it up
		if ( basis_map_->find( current_segment ) != basis_map_->end() ) {
			continue;
		} else {
			//Case 2: Segment is not a key in BasisMap but is in the EdgeMap -> use the EdgeMap to figure out which other segments we need to check for basis pairs.
			//NOTE: Under the new behavior everything should already be in EdgeMap, so this will sort of be the only reasonable case.
			//Edges_per_segment has first=current_segment, second=set of matching segment IDs
			EdgeMap::iterator edges_per_segment = edge_file_reader_->get_edges(current_segment);
			//edge_iter iterates over the set of segment IDs that match to current_segment (value is a core::Size )
			std::set< core::Size >::iterator edge_iter = edges_per_segment->second.begin();
			for ( ; edge_iter != edges_per_segment->second.end(); ++edge_iter ) {
				data_storage::SmartSegmentCOP seg2;
				if ( pdb_segments_.count( *edge_iter ) == 0 ) {
					seg2 = segment_vector_->at(*edge_iter );
				} else {
					seg2 = pdb_segments_.at( *edge_iter );
				}
				//Get alignments
				utility::vector1< std::pair< core::Size, core::Size > > alignments = alignment_generator_->get_all_alignments( *seg_it, seg2 );
				utility::vector1< std::pair< core::Size, core::Size > > reverse_alignments = alignment_generator_->get_all_alignments( seg2, *seg_it );
				//Add to basis map
				( ( *basis_map_ )[ current_segment ] )[ seg2->get_segment_id() ] = alignments;
				( ( *basis_map_ )[ seg2->get_segment_id() ] )[ current_segment ] = reverse_alignments;


				//Add OTHER SEGMENTS FROM THIS MODEL to next_set_to_check

				if ( !( seg2->is_n_terminus_fixed() ) ) { //The N-terminus is not fixed, meaning the segment is N-terminal
					//Add the C-terminal-most segment to the set
					next_set_to_check.insert( data_storage::SmartSegment::get_c_most_segment(seg2 , false ));
				} else if ( !( seg2->is_c_terminus_fixed() ) ) { //The C-terminus is not fixed, so the segment is C-terminal.
					//I use an else if because if the segment is both N- and C-terminal, it's already added, so no need to check again
					//Add the N-terminal-most segment to the set
					next_set_to_check.insert( data_storage::SmartSegment::get_c_most_segment(seg2 , false ));
				}
			}
		}

		//Finally perform our recursive call using our new set of segments
		recurse_over_segments( current_depth, next_set_to_check );
	}
}

//NOTE This function assumes you're handing it the N-terminal segment
data_storage::BasisPair
BasisMapGenerator::get_basis_pair_for_segment( core::Size segment_id )
{
	//If this segment is not in the BasisMap, add it
	if ( basis_map_->find( segment_id ) == basis_map_->end() ) {
		std::set< data_storage::SmartSegmentCOP > segs_to_check;

		if ( pdb_segments_.count( segment_id ) == 0 ) {
			segs_to_check.insert( ( *segment_vector_ )[ segment_id ] );
		} else {
			segs_to_check.insert( pdb_segments_.at( segment_id ) );
		}

		recurse_over_segments( 0, segs_to_check );
	}
	//Now we need to randomly select one basis pair for this segment
	utility::vector1< core::Size > second_segs;
	//Get all the keys from the map for this segment
	for ( std::map< core::Size, utility::vector1< std::pair< core::Size, core::Size > > >::iterator it = ( *basis_map_ )[ segment_id ].begin();
			it != ( *basis_map_ )[ segment_id ].end();
			++it ) {
		second_segs.push_back( it->first );
	}
	//Randomly select second segment id
	core::Size second_seg_id = 0;


	//We should probably just call get_basis_pair_for_segment_pair here to avoid code duplication
	if ( second_segs.size() > 0 ) {
		second_seg_id = numeric::random::rg().random_element( second_segs );
	}

	return get_basis_pair_for_segment_pair( segment_id, second_seg_id, true );

}


//If the returned bool is false, the first segment can never chimerize
std::pair< bool, data_storage::BasisPair >
BasisMapGenerator::get_basis_pair_for_local_segments( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2, bool first_is_N ){
	data_storage::BasisPair basis_pair = std::make_pair( data_storage::Basis( 0, 0 ), data_storage::Basis( 0, 0 ) );
	std::pair< bool, bool > can_chimerize = segs_can_chimerize( seg1, seg2, first_is_N );
	if ( can_chimerize.first == false ) {
		return std::make_pair( false, basis_pair );
	} else if ( can_chimerize.second == false ) {
		return std::make_pair( true, basis_pair );
	}
	//Now we know there's at least some possibility of these chimerizing

	//Check on vital residues
	std::set< core::Size > seg1_vital_res = seg1->get_vital_residues();
	std::set< core::Size > seg2_vital_res = seg2->get_vital_residues();
	core::Size i_first_vital_residue = 0;
	core::Size j_first_vital_residue = 0;
	if ( seg1_vital_res.size() != 0 ) {
		i_first_vital_residue = *(seg1_vital_res.begin() );
	}
	if ( seg2_vital_res.size() != 0 ) {
		j_first_vital_residue = *(seg2_vital_res.begin() );
	}

	std::pair< std::pair< core::Size, core::Size >, std::pair< core::Size, core::Size > > res_bounds = get_basis_res_bounds( seg1, seg2, first_is_N );
	//Note that these are bounds relative to the LOCAL segments.

	//check for chimaerae
	core::Size segID1 = 0;
	core::Size segID2 = 0;
	int seg1_correction_factor = 0;
	int seg2_correction_factor = 0;
	//CORRECTION = Cbasis - Nbasis
	if ( seg1->is_chimaeric() ) {
		//TODO
		if ( first_is_N ) {
			//If seg1 is N-terminal, we'll want to chimerize off its C-terminal parent
			//NOTE this assumes we aren't chimerizing onto something that's already a double chimaera
			segID1 = seg1->get_c_terminal_parent()->get_segment_id();

			//Residue numbering will need to be adjusted in later checks
			if ( seg1->get_basis_pair().first.segment_id() == segID1 ) { //C term is first, N term is second
				seg1_correction_factor = seg1->get_basis_pair().first.resnum() - seg1->get_basis_pair().second.resnum();
			} else { //Nterm is first, Cterm is second
				seg1_correction_factor = seg1->get_basis_pair().second.resnum() - seg1->get_basis_pair().first.resnum();
			}
		} else {
			//Otherwise chimerize off the N-terminal parent
			segID1 = seg1->get_n_terminal_parent()->get_segment_id();
			//Residue numbering will stay the same
		}
	} else {
		segID1 = seg1->get_segment_id();
	}
	if ( seg2->is_chimaeric() ) {
		if ( first_is_N ) {
			//If seg2 is C-terminal, we'll want to chimerize off its N-terminal parent
			segID2 = seg2->get_n_terminal_parent()->get_segment_id();
			//Residue numbering will stay the same
		} else {
			//Otherwise chimerize off the N-terminal parent
			segID2 = seg2->get_c_terminal_parent()->get_segment_id();
			//Residue numbering will need to be adjusted in later checks
			if ( seg2->get_basis_pair().first.segment_id() == segID2 ) { //C term is first, N term is second
				seg2_correction_factor = seg2->get_basis_pair().first.resnum() - seg2->get_basis_pair().second.resnum();
			} else { //Nterm is first, Cterm is second
				seg2_correction_factor = seg2->get_basis_pair().second.resnum() - seg2->get_basis_pair().first.resnum();
			}
		}
	} else {
		segID2 = seg2->get_segment_id();
	}

	bool found_valid_basis_pair = false;
	core::Size ntries = 0;
	while ( !found_valid_basis_pair && ntries < 1000 ) {
		++ntries;
		basis_pair = get_basis_pair_for_segment_pair( segID1, segID2, first_is_N );
		//See if the resnums fall within the chimaera bounds
		if ( seg1_correction_factor + int(basis_pair.first.resnum() ) <= 0 || seg2_correction_factor + int( basis_pair.second.resnum() ) <= 0 ) {
			continue;
		}
		//See if the resnums fall within res_bounds
		core::Size new_i = basis_pair.first.resnum() + seg1_correction_factor;
		core::Size new_j = basis_pair.second.resnum() + seg2_correction_factor;
		if ( new_i < res_bounds.first.first || new_i > res_bounds.first.second || new_j < res_bounds.second.first || new_j > res_bounds.second.second ) {
			continue;
		}
		//Final check for vital residues
		if ( ( j_first_vital_residue != 0 &&  new_i <= ( new_j - j_first_vital_residue ) ) || ( i_first_vital_residue != 0 && new_j <= ( new_i - i_first_vital_residue )  ) ) {
			//The basis pair is not valid
			continue;
		}

		//Convert back to local segment numbering
		basis_pair.first.segment_id( seg1->get_segment_id() );
		basis_pair.second.segment_id( seg2->get_segment_id() );
		basis_pair.first.resnum( new_i );
		basis_pair.second.resnum( new_j );
		found_valid_basis_pair = true;
	}

	if ( !found_valid_basis_pair ) { //If we never found a valid basis pair, return zeros
		basis_pair = std::make_pair( data_storage::Basis( 0, 0 ), data_storage::Basis( 0, 0 ) );
	}
	return std::make_pair( true, basis_pair );
}



std::pair< std::pair< core::Size, core::Size >, std::pair< core::Size, core::Size > >
BasisMapGenerator::get_basis_res_bounds( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2, bool first_is_N ){
	core::Size seg1_size = seg1->get_length();
	core::Size seg2_size = seg2->get_length();
	core::Size window_width = this->hasher_settings().min_hash_score / 4;
	//i is index within the first segment (already in the assembly, might be N or C terminal)
	core::Size i_start = 1;
	core::Size i_end = ( seg1_size + 1 - window_width );
	if ( seg1->is_chimaeric() ) { //if we are adding to a chimaera it will always be segment_1 since it'll already be in the assembly
		TR << "Adding to chimaera" << std::endl;
		core::Size basis_res;
		//we have to find which Basis is n_terminal
		//If the first basis is N-terminal, the numbering will be fine
		if ( seg1->get_n_terminal_parent()->get_segment_id() == seg1->get_basis_pair().first.segment_id() ) {
			basis_res = seg1->get_basis_pair().first.resnum();
		} else {
			//If the first basis is C-terminal, we should go by the N-terminal numbering
			basis_res = seg1->get_basis_pair().second.resnum();
		}
		TR << "Chimaera has basis res of: " << basis_res << std::endl;
		if ( first_is_N ) {
			//we are adding to the n-term
			//   TR << "adding to n_term. Setting i_end: ";
			i_end = basis_res - 1;
			TR << i_end << std::endl;
		} else {
			//we are adding to c term
			//   TR << "adding to c_term. Setting i: ";
			i_start = basis_res + 1;
			//TR << i << std::endl;
		}
	}
	std::set< core::Size > seg1_vital_res = seg1->get_vital_residues();
	std::set< core::Size > seg2_vital_res = seg2->get_vital_residues();
	if ( seg1_vital_res.size() != 0 ) {
		if ( first_is_N ) {
			//std::set is a sorted container, so this should work
			core::Size last_vital_residue = *(seg1_vital_res.rbegin() );
			i_start = std::max( i_start, last_vital_residue - window_width + 1 );
		} else {
			core::Size first_vital_residue = *(seg1_vital_res.begin() );
			//If seg1 is C-teminal, then the new basis residue should be <= ( first vital residue + window_width - 1 )
			i_end = std::min( i_end, first_vital_residue + window_width - 1 );
		}
	}
	core::Size j_start = 1;
	core::Size j_end = (seg2_size + 1 - window_width );
	//If seg2 is N-terminal, then the new basis residue should be >= ( last vital residue - window_width_ + 1 )
	//If seg2 is C-teminal, then the new basis residue should be <= ( first vital residue + window_width - 1 )
	if ( seg2_vital_res.size() != 0 ) {
		if ( first_is_N ) {
			core::Size last_vital_residue = *(seg2_vital_res.rbegin() );
			j_start = std::max( j_start, last_vital_residue - window_width + 1 );
		} else {
			core::Size first_vital_residue = *(seg2_vital_res.begin() );
			j_end = std::min( j_end, first_vital_residue + window_width - 1 );
		}
	}
	return std::make_pair( std::make_pair( i_start, i_end ), std::make_pair( j_start, j_end ) );
}





data_storage::BasisPair
BasisMapGenerator::get_basis_pair_for_segment_pair( core::Size seg1, core::Size seg2, bool first_is_N ){
	//Checking which is N-terminal is necessary because all pairs are stored in the basis map N first then C
	core::Size segID1 = 0;
	core::Size segID2 = 0;
	if ( first_is_N ) {
		segID1 = seg1;
		segID2 = seg2;
	} else {
		segID1 = seg2;
		segID2 = seg1;
	}
	//If this segment is not in the BasisMap, add it
	if ( basis_map_->find( segID1 ) == basis_map_->end() ) {
		std::set< data_storage::SmartSegmentCOP > segs_to_check;

		if ( pdb_segments_.count( segID1 ) == 0 ) {
			segs_to_check.insert( ( *segment_vector_ )[ segID1 ] );
		} else {
			segs_to_check.insert( pdb_segments_.at( segID1 ) );
		}
		recurse_over_segments( 0, segs_to_check );
	}

	std::pair< core::Size, core::Size > resnums = std::make_pair( 0, 0 );
	if ( basis_map_->at( segID1 ).at( segID2 ).size() > 0 ) {
		resnums = numeric::random::rg().random_element( ( *basis_map_ )[ segID1 ][ segID2 ] );
	}
	if ( resnums.first == 0 || segID2 == 0 ) {
		return std::make_pair( data_storage::Basis( 0, 0), data_storage::Basis( 0, 0) );
	}
	//Now make the basis pair from this info
	data_storage::Basis first_basis( segID1, resnums.first );
	data_storage::Basis second_basis( segID2, resnums.second );
	//now return the bases in the order they were given (not necessarily N first)
	if ( first_is_N ) {
		return std::make_pair( first_basis, second_basis );
	} else {
		return std::make_pair( second_basis, first_basis );
	}
}


std::pair< bool, bool>
BasisMapGenerator::segs_can_chimerize( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2, bool first_is_N ){
	//Note that this DOES allow chimaeric segments to be considered
	core::Size window_width = this->hasher_settings().min_hash_score / 4;
	if ( !( seg1->is_hashable() && seg2->is_hashable()
			&& ( seg1->get_dssp_code() == seg2->get_dssp_code() ) )
			) {
		//TR << "Cannot hash segments!" << std::endl;
		if ( !seg1->is_hashable() ) { //seg1 can never be hashed, first bool is false
			return std::make_pair( false, false );
		} else {
			return std::make_pair(true, false);
		}
	}//endif

	// TR << "Producing basis pairs within " << utility::to_string( window_width_ ) << " residues of termini" << std::endl;
	//Find number of residues in each segment
	core::Size seg1_size = seg1->get_length();
	core::Size seg2_size = seg2->get_length();
	//If either segment is shorter than required_num_residues, no basis pairs possible.
	if ( seg1_size < window_width || seg2_size < window_width ) {
		TR << "One of the segments is too short!" << std::endl;
		if ( seg1_size < window_width ) {
			//Seg1 is bad
			return std::make_pair( false, false );
		} else {
			//Seg1 is OK
			return std::make_pair( true, false );
		}
	}
	std::pair< std::pair< core::Size, core::Size>, std::pair< core::Size, core::Size > > res_bounds = get_basis_res_bounds( seg1, seg2, first_is_N );

	if ( res_bounds.first.first > res_bounds.first.second + window_width || seg1->get_dssp_code() == 'L' ) {
		return std::make_pair( false, false);
	}
	if ( res_bounds.second.first > res_bounds.second.second + window_width || seg2->get_dssp_code() == 'L' ) {
		return std::make_pair( true, false);
	}
	return std::make_pair( true, true );
}


bool
BasisMapGenerator::basis_pair_is_valid( data_storage::BasisPair basis_pair, bool first_is_N ){
	core::Size segID1;
	core::Size segID2;
	core::Size res1;
	core::Size res2;
	if ( first_is_N ) {
		segID1 = basis_pair.first.segment_id();
		segID2 = basis_pair.second.segment_id();
		res1 = basis_pair.first.resnum();
		res2 = basis_pair.second.resnum();
	} else {
		segID1 = basis_pair.second.segment_id();
		segID2 = basis_pair.first.segment_id();
		res1 = basis_pair.second.resnum();
		res2 = basis_pair.first.resnum();
	}
	if ( basis_map_->count( segID1 ) == 0 ) {
		//We need to recurse for this segID
		std::set< data_storage::SmartSegmentCOP > segs_to_check;
		if ( pdb_segments_.count( segID1 ) == 0 ) {
			segs_to_check.insert( ( *segment_vector_ )[ segID1 ] );
		} else {
			segs_to_check.insert( pdb_segments_.at( segID1 ) );
		}
		recurse_over_segments( 0, segs_to_check );
	}
	if ( basis_map_->count( segID1 ) == 0 || basis_map_->at( segID1 ).count( segID2 ) == 0 ) {
		return false;
	}
	// if( basis_map_->at( segID1 ).at( segID2 ).find( std::make_pair( res1, res2 ) ) != basis_map_->at( segID1 ).at( segID2 ).end() ){
	if ( std::find( basis_map_->at( segID1 ).at( segID2 ).begin(), basis_map_->at( segID1 ).at( segID2 ).end(), std::make_pair( res1, res2 ) ) != basis_map_->at( segID1 ).at( segID2 ).end() ) {
		return true;
	}
	return false;
}



//Getters
core::Size
BasisMapGenerator::version() const {
	return version_;
}

HasherSettings
BasisMapGenerator::hasher_settings() const{
	return hasher_settings_;
}

SegmentVectorCOP
BasisMapGenerator::segment_vector() const{
	return segment_vector_;
}
EdgeMapGeneratorOP
BasisMapGenerator::edge_file_reader() const{
	return edge_file_reader_;
}

AlignmentSettings
BasisMapGenerator::alignment_settings() const{
	return alignment_settings_;
}

AlignmentGeneratorOP
BasisMapGenerator::alignment_generator() const{
	return alignment_generator_;
}

utility::pointer::shared_ptr< BasisMap >
BasisMapGenerator::basis_map() const{
	return basis_map_;
}

std::map< core::Size, data_storage::SmartSegmentOP > &
BasisMapGenerator::pdb_segments(){
	return pdb_segments_;
}


std::map< core::Size, data_storage::SmartSegmentOP >
BasisMapGenerator::const_pdb_segments() const{
	return pdb_segments_;
}

//Setters

void
BasisMapGenerator::version( core::Size ver ){
	version_ = ver;
}

void
BasisMapGenerator::pdb_segments( std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs ){
	pdb_segments_ = pdbsegs;
}

void
BasisMapGenerator::segment_vector( SegmentVectorCOP segvec ){
	version_ = segvec->get_version();
	segment_vector_ = segvec;
}

void
BasisMapGenerator::alignment_settings( AlignmentSettings & as ){
	alignment_settings_ = as;
}
void
BasisMapGenerator::alignment_generator( AlignmentGeneratorOP ag ){
	alignment_generator_ = ag;
}

void
BasisMapGenerator::basis_map( utility::pointer::shared_ptr< BasisMap > bm ){
	basis_map_ = bm;
}

void
BasisMapGenerator::hasher_settings( HasherSettings & hs ){
	hasher_settings_ = hs;
}



} //hashing
} //sewing
} //protocols


