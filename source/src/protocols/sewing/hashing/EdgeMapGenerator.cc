// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/EdgeMapGenerator.cc
/// @brief Class used to generate edge files for use with SEWING protocols.
/// @author guffysl (guffy@email.unc.edu)

//Sewing headers
#include <protocols/sewing/hashing/EdgeMapGenerator.hh>
#include <protocols/sewing/hashing/Hasher.hh>
//#include <protocols/sewing/conformation/Model.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
//#include <protocols/sewing/sampling/SewGraph.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/data_storage/Basis.hh>
//Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

//Utility headers
#include <utility/mpi_util.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
static basic::Tracer TR( "protocols.sewing.hashing.EdgeMapGenerator" );


namespace protocols {
namespace sewing {
namespace hashing {

EdgeMapGenerator::EdgeMapGenerator():
	utility::pointer::ReferenceCount(),
	edge_file_name_(""),
	is_edge_file_loaded_( false )
{
	pdb_segments_.clear();
}

EdgeMapGenerator::EdgeMapGenerator( std::string edge_file_name ):
	edge_file_name_( edge_file_name ),
	is_edge_file_loaded_( false )
{
	//edges_ = nullptr;
	//We'll try to read the file with the given name.
	//Since this class can read or write edge files, maybe we shouldn't read it in a constructor?
	//initialize_from_options();
	//just override the edge file name
	pdb_segments_.clear();
}

EdgeMapGenerator::EdgeMapGenerator(std::string model_file_name, core::Size boxes_per_dimension, core::Size min_hash_score, core::Size max_clash_score, std::string edge_file_name, bool hash_opposite_termini ):
	edge_file_name_(edge_file_name),
	is_edge_file_loaded_( false )
{
	//initialize_from_options();
	//edges_ = nullptr;
	pdb_segments_.clear();
	set_private_data( model_file_name, boxes_per_dimension, min_hash_score, max_clash_score , hash_opposite_termini);
}

EdgeMapGenerator::EdgeMapGenerator( EdgeMapGenerator const & other ):
	hasher_settings_( other.hasher_settings_ ),
	version_(other.version_ ),
	edges_( other.edges_ ),
	segment_vector_( other.segment_vector_ ),
	edge_file_name_(other.edge_file_name_),
	is_edge_file_loaded_( other.is_edge_file_loaded_ ),
	pdb_segments_( other.pdb_segments_ )
{}






EdgeMapGenerator::~EdgeMapGenerator(){}

EdgeMapGeneratorOP
EdgeMapGenerator::clone(){
	return EdgeMapGeneratorOP( new EdgeMapGenerator( *this ) );
}


void EdgeMapGenerator::initialize_from_options(){
	core::Size max_clash_score = basic::options::option[ basic::options::OptionKeys::sewing::max_clash_score ].value();
	core::Size min_hash_score = basic::options::option[ basic::options::OptionKeys::sewing::min_hash_score ].value();
	int boxes_per_dimension = basic::options::option[ basic::options::OptionKeys::sewing::box_length ].value();
	if ( !basic::options::option[ basic::options::OptionKeys::sewing::model_file_name ].user() ) {
		//utility_exit_with_message("Error: Please provide a SEWING model file to use for edge generation." );
		//We no longer always need a model file, so just exit later if it isn't there
	}
	std::string model_file_name = basic::options::option[ basic::options::OptionKeys::sewing::model_file_name ].value();

	if ( !basic::options::option[ basic::options::OptionKeys::sewing::edge_file_name ].user() ) {
		utility_exit_with_message("Error: Please provide a name for the edge file." );
	}
	bool hash_opposite_termini = basic::options::option[ basic::options::OptionKeys::sewing::score_between_opposite_terminal_segments].value();

	edge_file_name_ = basic::options::option[ basic::options::OptionKeys::sewing::edge_file_name ].value();
	set_private_data( model_file_name, core::Size( boxes_per_dimension ), min_hash_score, max_clash_score, hash_opposite_termini );
}

void EdgeMapGenerator::set_private_data( std::string model_file_name, core::Size boxes_per_dimension, core::Size min_hash_score, core::Size max_clash_score, bool hash_opposite_termini ){

	if ( model_file_name == "" ) {
		segment_vector_ = SegmentVectorCOP( new SegmentVector );
	} else {
		//Not actually sure what the name of this class/method will be, but let's go with this one
		segment_vector_ = ModelFileReader::read_model_file( model_file_name );
		version_ = segment_vector_->get_version();
	}
	//model_map_ = model_output.second;
	if ( boxes_per_dimension != 3 && boxes_per_dimension != 5 ) {
		TR << "Invalid box dimensions provided. Reverting to default (3)." << std::endl;
		boxes_per_dimension = 3;
	}
	hasher_settings_ = HasherSettings(max_clash_score, min_hash_score, boxes_per_dimension, hash_opposite_termini);

}

//Getters

core::Size EdgeMapGenerator::version() const { return version_;}
std::string EdgeMapGenerator::edge_file_name() const { return edge_file_name_; }
HasherSettings EdgeMapGenerator::hasher_settings() const { return hasher_settings_; }
//utility::pointer::shared_ptr< std::map< int, Model > > EdgeMapGenerator::model_map() { return model_map_; }
SegmentVectorCOP EdgeMapGenerator::segment_vector() const{ return segment_vector_; }

EdgeMapOP EdgeMapGenerator::edge_map() const { return edges_; }
bool EdgeMapGenerator::is_edge_file_loaded() const { return is_edge_file_loaded_; }
std::map< core::Size, data_storage::SmartSegmentOP > &
EdgeMapGenerator::pdb_segments(){
	return pdb_segments_;
}

std::map< core::Size, data_storage::SmartSegmentOP > const
EdgeMapGenerator::const_pdb_segments() const{
	return pdb_segments_;
}


//Setters

void EdgeMapGenerator::version( core::Size ver){ version_ = ver;}
void EdgeMapGenerator::edge_file_name( std::string name){ edge_file_name_ = name; }
//void EdgeMapGenerator::model_map( utility::pointer::shared_ptr< std::map< int, Model > > modelmap){ model_map_ = modelmap; }
void EdgeMapGenerator::segment_vector( SegmentVectorCOP segvector )
{
	segment_vector_ = segvector;
}
void
EdgeMapGenerator::pdb_segments( std::map< core::Size, data_storage::SmartSegmentOP > & pdbsegs ){
	pdb_segments_ = pdbsegs;
}
void EdgeMapGenerator::hasher_settings( HasherSettings & hs ){ hasher_settings_ = hs; }



bool EdgeMapGenerator::edge_exists( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2, HasherOP hasher ){
	//We'll be getting two iterators from this function that indicate the beginning and end of a list of basis pairs.
	//This list is stored in the hasher and only contains residue numbers since we already have the segment here.
	std::pair< utility::vector1< std::pair<core::Size, core::Size > >::const_iterator, utility::vector1< std::pair< core::Size, core::Size > >::const_iterator > iters = hasher->iterate_over_basis_pairs(seg1, seg2 );
	if ( iters.first == iters.second ) {
		TR << "No basis pairs detected" << std::endl;
		return false;
	}
	//Should work, even though it looks a little weird
	//This is apparently never being called :(
	for ( ; iters.first != iters.second; ++iters.first ) {
		//Create the basis pair (nothing more than a pair of pairs of
		//iters.first points to the current pair of residue numbers
		TR.Debug << "New basis pair: Segment " << seg1->get_segment_id() << " residue " << iters.first->first << ", segment " << seg2->get_segment_id() << " residue " << iters.first->second << std::endl;
		data_storage::BasisPair current_basis = std::make_pair( data_storage::Basis(seg1->get_segment_id(), iters.first->first), data_storage::Basis( seg2->get_segment_id(), iters.first->second ) );
		if ( hasher->score_basis_pair( current_basis ) ) {
			return true;
		}
	}
	return false;

}

///@details This function is intended to be used with AppendAssemblyMover
///( and possibly with some later applications) to add a set of connected segments (formerly known as a model) to an existing EdgeMap

void EdgeMapGenerator::append_to_edge_map( core::Size base_seg_id ){
	//First we'll want to find the N-terminal and C-terminal segments for this segment
	//Then hash those segments into the EdgeMap
	data_storage::SmartSegmentCOP n_seg;
	data_storage::SmartSegmentCOP c_seg;

	//N-terminal segment
	if ( pdb_segments_.count( base_seg_id ) == 0 ) {
		n_seg = data_storage::SmartSegment::get_n_most_segment( segment_vector_->at( base_seg_id ) , false );
	} else {
		n_seg = data_storage::SmartSegment::get_n_most_segment( pdb_segments_.at( base_seg_id ), false );
	}
	//C-terminal segment
	if ( pdb_segments_.count( base_seg_id ) == 0 ) {
		c_seg = data_storage::SmartSegment::get_c_most_segment( segment_vector_->at( base_seg_id ) , false );
	} else {
		n_seg = data_storage::SmartSegment::get_c_most_segment( pdb_segments_.at( base_seg_id ), false );
	}

	//So this function will also need to have a hasher & to generate lists of N- and C-terminal segments just like when making the edge file
	//Hasher
	HasherOP append_hasher( new Hasher( hasher_settings_, segment_vector_, pdb_segments_ ) );
	//Lists of terminal segments
	std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > term_segs = get_terminal_segments();
	utility::vector1< core::Size > n_segs = term_segs.first;
	utility::vector1< core::Size > c_segs = term_segs.second;
	//Now hash N-terminal segment against all C-terminal segments and C-terminal segment against all N-terminal segments
	for ( core::Size c_index = 1; c_index <= c_segs.size(); ++c_index ) {
		//for readability
		core::Size c_seg_id = c_segs[ c_index ];
		data_storage::SmartSegmentCOP test_c_seg;
		if ( pdb_segments_.count( c_seg_id ) == 0 ) {
			test_c_seg = segment_vector_->at( c_seg_id );
		} else {
			test_c_seg = pdb_segments_.at( c_seg_id );
		}
		if ( edge_exists( n_seg, test_c_seg , append_hasher ) ) {
			//Add in both directions
			(*edges_)[ n_seg->get_segment_id() ].insert( c_seg_id );
			(*edges_)[ c_seg_id ].insert( n_seg->get_segment_id() );
		}
	}
	//Now the same for the C-terminal segment
	//for( core::Size n_index = 1; n_index <= n_segs.size(); ++n_index ) {
	for ( core::Size n_seg_id: n_segs ) {
		//core::Size n_seg_id = n_segs[ n_index ];

		data_storage::SmartSegmentCOP test_n_seg;
		if ( pdb_segments_.count( n_seg_id ) == 0 ) {
			test_n_seg = segment_vector_->at( n_seg_id );
		} else {
			test_n_seg = pdb_segments_.at( n_seg_id );
		}
		if ( edge_exists( test_n_seg, c_seg, append_hasher ) ) {
			(*edges_)[ n_seg_id ].insert( c_seg->get_segment_id() );
			(*edges_)[ c_seg->get_segment_id() ].insert( n_seg_id );
		}

	}
}



void EdgeMapGenerator::write_edges_for_segment(core::Size current_n_segment, HasherOP hasher, utility::io::ozstream & output_file, utility::vector1< core::Size > const & c_segments ) {

	//New function is simple--just score current_n_segment against all c_segments
	//Iterate over all segments in c_segments
	//for( core::Size c_seg_id = 1; c_seg_id <= c_segments.size(); ++c_seg_id ) {
	for ( core::Size c_seg_id: c_segments ) {
		TR.Debug << "Searching for an edge between segment " << current_n_segment << " and segment " << c_seg_id << std::endl;


		data_storage::SmartSegmentCOP n_seg;
		data_storage::SmartSegmentCOP c_seg;
		//if( pdb_segments_.count( current_n_segment ) == 0 ){
		if ( current_n_segment <= segment_vector_->size() ) {
			n_seg = segment_vector_->at( current_n_segment );
		} else {
			TR << "Reading segment " << current_n_segment << " from pdb_segments_" << std::endl;
			n_seg = pdb_segments_.at( current_n_segment );
		}

		if ( c_seg_id <= segment_vector_->size() ) {
			c_seg = segment_vector_->at( c_seg_id );
		} else {
			TR << "Reading segment " << c_seg_id << " from pdb_segments_" << std::endl;
			c_seg = pdb_segments_.at( c_seg_id );
		}

		if ( edge_exists( n_seg, c_seg, hasher ) ) {

			TR.Debug << "Found an edge between " << current_n_segment << " and " << c_seg_id << "." << std::endl;
			output_file << utility::to_string( current_n_segment ) << " " << utility::to_string( c_seg_id ) << std::endl;

		}
	}
}



std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > EdgeMapGenerator::get_terminal_segments()
{
	//Initialize the two vectors
	utility::vector1< core::Size > n_segs;
	utility::vector1< core::Size > c_segs;
	//First loop through the segment_vector_
	for ( core::Size i = 1; i <= segment_vector_->size(); ++i ) {
		//For each segment, if it is N-terminal, append to n_segs, and if it is C-terminal, append to c_segs
		if ( ! ( (* segment_vector_ )[ i ]->is_n_terminus_fixed() ) ) { //This segment is N-terminal (its N terminus is free)
			n_segs.push_back( i );
		}
		if ( !( ( *segment_vector_ )[ i ]->is_c_terminus_fixed() ) ) { //This segment is C-terminal (note it is possible for it to be both)
			c_segs.push_back( i );
		}
	}
	//Now loop through pdb_segments_
	for ( std::pair< core::Size, data_storage::SmartSegmentOP > pdbseg: pdb_segments_ ) {
		if ( !pdbseg.second->is_n_terminus_fixed() ) {
			n_segs.push_back( pdbseg.first );
		}
		if ( !pdbseg.second->is_c_terminus_fixed() ) {
			c_segs.push_back( pdbseg.first );
		}
	}
	return std::make_pair( n_segs, c_segs );
}


///@details This function handles the distribution of labor for edge file
///generation. The actual calls to the hasher & writing to the file (really
///everything directly dealing with SEWING elements) happens elsewhere.
void EdgeMapGenerator::generate_edge_file(){
	//Different for MPI vs non-MPI (we'll almost always be using MPI)
	TR << "Begin generate_edge_file" << std::endl;
	//Since we're definitely generating an edge file here, we'll need a model file, so we can go ahead and exit if one is not provided to us.
	if ( segment_vector_->size() == 0 && pdb_segments_.empty() ) {
		utility_exit_with_message( "Cannot generate edge file without an input model/segment file!" );
	}
	//We should also double check that a name was provided for the edge file.
	if ( edge_file_name_ == "" ) {
		utility_exit_with_message( "No edge file name provided for output edge file!" );
	}



#ifdef USEMPI
	//In MPI mode, we'll need to handle the master node and worker nodes separately
	core::Size rank = utility::mpi_rank();
	core::Size num_procs = utility::mpi_nprocs();
	TR << " MPI mode. Rank " << rank << " with " << num_procs << " processors." << std::endl;


	// ***********BEGIN MASTER NODE ********************
	if(rank == 0) {

		//This is the master node. The only thing it writes to the edge file is the version.
		//The file name will be <edge_file_name_>.0
		std::string master_edge_file_name = edge_file_name_ + "." + utility::to_string(rank);
		utility::io::ozstream file;
		file.open_append(master_edge_file_name);
		file << version_ << std::endl;
		//We'll also want to write the hasher settings in this file for AppendAssemblyMover, etc.
		file << utility::to_string( hasher_settings_.max_clash_score ) << " "
				 << utility::to_string( hasher_settings_.min_hash_score ) << " "
				 << utility::to_string( hasher_settings_.boxes_per_dimension ) << " "
				 << utility::to_string( hasher_settings_.hash_between_termini ) << std::endl;
		file.close();


		//Function to get vectors of N-terminal and C-terminal segments
		std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > terminal_seg_ids = get_terminal_segments();
		utility::vector1< core::Size > n_terminal_segs = terminal_seg_ids.first;

		//Now the master node will hand out the first set of models (1 to each processor)

		//Make sure the user isn't getting overzealous with the processors (this seems like an unlikely scenario . . . )
		if(num_procs > n_terminal_segs.size() ) {
			utility_exit_with_message("You have more processors than number of segments to hash. Reduce your number of processors for efficiency");
		}

		TR << "MPI mode with " << num_procs << " processors" << std::endl;
		//Start at processor #1
		core::Size curr_proc = 1;
		//		for(; it != model_map_->end(); ++it) {
		core::Size i = 1;
		for( ; i <= n_terminal_segs.size(); ++i ) {
			//Send the segment ID of the desired segment to the current processor
			utility::send_integer_to_node(curr_proc, n_terminal_segs[ i ]  );
			TR << "Master node sent a job for segment " << n_terminal_segs[ i ] << " to processor " << curr_proc << std::endl;
			//Move on to the next processor
			++curr_proc;

			//Stop when we've finished the last processor (the last processor # will be num_procs - 1)
			if(curr_proc == num_procs) { break; }
		}
		TR << "Finished sending jobs for each processor. Final value of curr_proc: " << curr_proc << std::endl;
		++i;//Increment counter once more to account for the last segment sent in the above loop

		while( i <= n_terminal_segs.size() ){
			//When a processor has completed the segment it was sent, send it a new one
			core::Size received_node = utility::receive_integer_from_anyone();
			TR << "Master node received a message from processor " << received_node << std::endl;
			utility::send_integer_to_node(received_node, n_terminal_segs[ i ] );
			TR << "Master node sent a new job for segment " << n_terminal_segs[ i ] << " to processor " << received_node << std::endl;
			++i;
		}

		//Once we're done with all the work, tell each processor to spin down
		core::Size counter = 1;
		while(counter != num_procs){
			core::Size received_node = utility::receive_integer_from_anyone();
			utility::send_integer_to_node(received_node,0);
			++counter;
		}
		TR << "Master node finished sending jobs" << std::endl;
  }//End master node

	// ***********BEGIN WORKER NODE ********************
	else{
		TR << "Tracer file for worker node " << rank << std::endl;
		utility::io::ozstream node_file;
		std::string worker_edge_file_name = edge_file_name_ + "." + utility::to_string(rank);
		TR << "Node will output to file " << worker_edge_file_name << std::endl;
		node_file.open_append(worker_edge_file_name);

		//Each worker node should have its own Hasher object. Go ahead and create those now so we only have to do it once.
 HasherOP worker_node_hasher( new Hasher( hasher_settings_, segment_vector_, pdb_segments_ ) );
		TR << "Hasher created" << std::endl;

		//Get vectors of n- and c-terminal segments
		std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > terminal_seg_ids = get_terminal_segments();
		utility::vector1< core::Size > c_segs = terminal_seg_ids.second;

		while( true ){ //Loop until master tells you to stop!

			//Receive message from master node
			core::Size n_seg_id = utility::receive_integer_from_node( 0 );
			TR << "Preparing to find edges for segment " << n_seg_id << std::endl;
			//Master node will send 0 when work is finished
			if( n_seg_id == 0 ){ break; }

			//Now we'll score this individual model against all previous models
			TR << "Processor " << rank << " received job for segment number " << n_seg_id << std::endl;

			//Provide the file for this node to write to so we don't have to keep opening and closing it
			write_edges_for_segment( n_seg_id, worker_node_hasher, node_file, c_segs );

			//Let the master node know when we're done with this model
			utility::send_integer_to_node(0, rank);
		}
		node_file.close();
		TR << "Processor " << rank << " was told to stop working" << std::endl;
	}
	MPI_Barrier( MPI_COMM_WORLD ); // make all nodes reach this point together.
	MPI_Finalize();


#else
	//Not in MPI mode


	//Create the hasher
	HasherOP hasher( new Hasher( hasher_settings_, segment_vector_, pdb_segments_ ) );


	//First we need to open the file and write the version at the top
	utility::io::ozstream file;
	file.open_append(edge_file_name_);
	file << version_ << std::endl;
	file << utility::to_string( hasher_settings_.max_clash_score ) << " "
		<< utility::to_string( hasher_settings_.min_hash_score ) << " "
		<< utility::to_string( hasher_settings_.boxes_per_dimension ) << " "
		<< utility::to_string( hasher_settings_.hash_between_termini ) << std::endl;


	//Get lists of n- and c-terminal segments

	std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > terminal_seg_ids = get_terminal_segments();
	utility::vector1< core::Size > n_segs = terminal_seg_ids.first;
	TR << "Found " << n_segs.size() << "n_terminal segments." << std::endl;
	utility::vector1< core::Size > c_segs = terminal_seg_ids.second;
	TR << "Found " << c_segs.size() << "c_terminal segments." << std::endl;
	//Now simply loop through the segments one at a time and write to the given edge file name
	for ( core::Size ii = 1; ii <= n_segs.size(); ++ii ) {
		write_edges_for_segment( n_segs[ ii ] , hasher, file, c_segs );
	}

#endif

}//generate_edge_file



void
EdgeMapGenerator::load_edge_file( std::string edge_file_name ) {
	edge_file_name_ = edge_file_name;
	read_edge_file();
}//load_edge_file

void
EdgeMapGenerator::read_edge_file() {
	//attempt to open edge file
	utility::io::izstream edge_file( edge_file_name_ );
	if ( !edge_file.good() ) {
		utility_exit_with_message( "Could not find file with name " + edge_file_name_ );
	}

	std::string line;

	//The first line is the version
	edge_file.getline( line );
	version_ = utility::string2Size( line );

	//The second line contains the HasherSettings data
	edge_file.getline( line );
	utility::vector1< std::string > settings = utility::string_split( line );

	//Get values for HasherSettings
	core::Size max_clash_score = utility::string2Size( settings[ 1 ] );
	core::Size min_hash_score = utility::string2Size( settings[ 2 ] );
	core::Size boxes_per_dimension = utility::string2Size( settings[ 3 ] );
	bool hash_opposite_termini = bool( utility::string2int( settings[ 4 ] ) );
	hasher_settings_ = HasherSettings( max_clash_score, min_hash_score, boxes_per_dimension, hash_opposite_termini );

	//The rest of the file contains edges
	utility::vector1< std::string > current_edge;
	core::Size seg1;
	core::Size seg2;
	//std::pair < core::Size, core::Size > current_key;
	while ( edge_file.getline( line ) ) {
		//Each line has seg1 seg2
		current_edge = utility::string_split( line );
		seg1 = utility::string2Size( current_edge[ 1 ] );
		seg2 = utility::string2Size( current_edge[ 2 ] );
		//current_key = std::make_pair( mod1, seg1 );
		//If this segment doesn't have edges yet, add it to the edge map
		//Actually not sure if this is necessary anymore . . .
		/*
		if ( edges_.find( seg1 ) == edges_.end() ) { //If the N-terminal edge is not yet in the map
		//   edges_[ current_key ] = utility::vector1< std::pair< core::Size, core::Size > >;
		edges_[ seg1 ].clear();
		}
		*/
		//In any case, add this edge
		( (*edges_)[ seg1 ] ).insert( seg2 );
		( (*edges_)[ seg2 ] ).insert( seg1 );
	}
	//Now we've read in all the edges

	edge_file.close();
	is_edge_file_loaded_ = true;
}//read_edge_file()






EdgeMap::iterator
EdgeMapGenerator::get_edges( core::Size sew_segment ) {
	return edges_->find( sew_segment );
}

} //hashing
} //sewing
} //protocols






