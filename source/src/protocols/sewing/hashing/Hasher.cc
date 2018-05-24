// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/hashing/Hasher.cc
/// @brief A geometric hashing class used by the SEWING protocol
/// @author Minnie Langlois (minnie@email.unc.edu)


//External headers
#include <boost/math/special_functions/round.hpp>

//Project headers
#include <protocols/sewing/hashing/Hasher.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
//Protocols headers

//Core headers
#include <core/id/AtomID.hh>
#include <core/conformation/Atom.hh>

//Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

//Utility headers
#include <utility/vector1.hh>
#include <numeric/HomogeneousTransform.hh>
#include <utility/utility.functions.hh>

static basic::Tracer TR( "protocols.sewing.hashing.Hasher" );


namespace protocols {
namespace sewing {
namespace hashing {


// Constructors

void
Hasher::allocate_segments(){
	core::Size const larger_than_possible_segment_size = 100;
	seg1_ = data_storage::SmartSegmentOP( new data_storage::SmartSegment( false, larger_than_possible_segment_size ) );
	seg2_ = data_storage::SmartSegmentOP( new data_storage::SmartSegment( false, larger_than_possible_segment_size ) );
	int const max_substructure_length = 5; //max number of seg in a substructure = 5, max num neighbors = 4
	seg1_substructure_ = utility::vector1< data_storage::SmartSegmentOP >( max_substructure_length );
	seg2_substructure_ = utility::vector1< data_storage::SmartSegmentOP >( max_substructure_length );
	for ( int i=1; i <= max_substructure_length; i++ ) {
		seg1_substructure_[i] = data_storage::SmartSegmentOP( new data_storage::SmartSegment( false, larger_than_possible_segment_size ) );
		seg2_substructure_[i] = data_storage::SmartSegmentOP( new data_storage::SmartSegment( false, larger_than_possible_segment_size ) );
	}
}

Hasher::Hasher():
	utility::pointer::ReferenceCount()
{
	all_basis_pairs_.clear();
	//initialize_from_options();
	allocate_segments();
}

Hasher::Hasher( Hasher const & other) {
	//Make a deep copy of the segment vector
	/*
	segment_vector_.clear();
	for( core::Size i = 1; i <= other.segment_vector().size(); ++i ){
	data_storage::SmartSegmentOP seg = other.segment_vector().at( i );
	segment_vector_.push_back( data_storage::SmartSegmentOP( new data_storage::SmartSegment( *seg ) ) );
	if( i > 1 && seg->is_n_terminus_fixed() ){
	data_storage::SmartSegment::link_to( segment_vector_[ i - 1 ], segment_vector_[ i ] );
	}
	}
	*/
	segment_vector_ = other.segment_vector();
	all_basis_pairs_.clear();
	hasher_settings_ = other.hasher_settings();
	pdb_segments_ = other.const_pdb_segments();
	allocate_segments();
}

void
Hasher::initialize_from_options() {
	core::Size max_clash_score = basic::options::option[ basic::options::OptionKeys::sewing::max_clash_score ].value();
	core::Size min_hash_score = basic::options::option[ basic::options::OptionKeys::sewing::min_hash_score ].value();
	core::Size boxes_per_dimension = basic::options::option[ basic::options::OptionKeys::sewing::box_length ].value();
	bool hash_between_termini = basic::options::option[ basic::options::OptionKeys::sewing::score_between_opposite_terminal_segments].value();
	hasher_settings_ = HasherSettings(max_clash_score, min_hash_score, boxes_per_dimension, hash_between_termini);
}

Hasher::Hasher( HasherSettings & hasher_settings, SegmentVectorCOP segment_vector ) {
	if ( hasher_settings.boxes_per_dimension != 3 &&
			hasher_settings.boxes_per_dimension != 5 ) {
		TR << "boxes_per_dimension must be either 3 or 5!!" << std::endl;
		exit(1);
	}
	hasher_settings_ = hasher_settings;
	all_basis_pairs_.clear();
	//Make a deep copy of the segment vector
	//segment_vector_.clear();
	//for( data_storage::SmartSegmentOP seg: segment_vector ){
	/*
	for( core::Size i = 1; i <= segment_vector.size(); ++i ){
	data_storage::SmartSegmentOP seg = segment_vector[ i ];
	segment_vector_.push_back( data_storage::SmartSegmentOP( new data_storage::SmartSegment( *seg ) ) );
	if( i > 1 && seg->is_n_terminus_fixed() ){
	data_storage::SmartSegment::link_to( segment_vector_[ i - 1 ], segment_vector_[ i ] );
	}
	}*/
	segment_vector_ = segment_vector;
	allocate_segments();
}


Hasher::Hasher( HasherSettings & hasher_settings, SegmentVectorCOP segment_vector, std::map< core::Size, data_storage::SmartSegmentOP > pdb_segments )
{
	if ( hasher_settings.boxes_per_dimension != 3 &&
			hasher_settings.boxes_per_dimension != 5 ) {
		TR << "boxes_per_dimension must be either 3 or 5!!" << std::endl;
		exit(1);
	}
	hasher_settings_ = hasher_settings;
	all_basis_pairs_.clear();
	segment_vector_ = segment_vector;
	pdb_segments_ = pdb_segments;
	allocate_segments();

}


HasherOP
Hasher::clone() const {
	return HasherOP( new Hasher( *this ) );
}

Hasher::~Hasher(){}

// Getters/Setters
SegmentVectorCOP
Hasher::segment_vector() const {
	return segment_vector_;
}

void
Hasher::segment_vector( SegmentVectorCOP segvec ) {
	/*
	//Make a deep copy of the segment vector
	segment_vector_.clear();
	for( core::Size i = 1; i <= segvec.size(); ++i ){
	data_storage::SmartSegmentOP seg = segvec[ i ];
	segment_vector_.push_back( data_storage::SmartSegmentOP( new data_storage::SmartSegment( *seg ) ) );
	if( i > 1 && seg->is_n_terminus_fixed() ){
	data_storage::SmartSegment::link_to( segment_vector_[ i - 1 ], segment_vector_[ i ] );
	}
	}
	*/
	segment_vector_ = segvec;
}

HasherSettings
Hasher::hasher_settings() const{
	return hasher_settings_;
}

void
Hasher::hasher_settings( HasherSettings const & hasher_settings){
	hasher_settings_ = hasher_settings;
}

std::map< core::Size, data_storage::SmartSegmentOP > &
Hasher::pdb_segments(){
	return pdb_segments_;
}

void
Hasher::pdb_segments( std::map< core::Size, data_storage::SmartSegmentOP > & pdbsegs ){
	pdb_segments_ = pdbsegs;
}
std::map< core::Size, data_storage::SmartSegmentOP >
Hasher::const_pdb_segments() const{
	return pdb_segments_;
}


///@details Construct HomogenousTransform using the 3 points in the BasisSet. Use this
///HomogenousTransform to transform all features into the local coordinate frame
void
//Hasher::transform_segments( data_storage::SmartSegmentOP n_term_segment, data_storage::SmartSegmentOP basis_segment, core::Size basis_residue ) {
Hasher::transform_segments( data_storage::Basis current_basis, bool first ){
	TR.Debug << "Beginning transform" << std::endl;
	//Copy the corresponding segments
	core::Size sub_index = 1;
	//If this is already in pdb_segments, just copy the owning pointers
	data_storage::SmartSegmentCOP current_segment;
	if ( pdb_segments_.count( current_basis.segment_id() ) == 0 ) {
		current_segment = data_storage::SmartSegment::get_n_most_segment( segment_vector_->at( current_basis.segment_id() ), false );
	} else {
		current_segment = pdb_segments_.at( current_basis.segment_id() );
	}

	while ( current_segment != nullptr ) {
		//Set all the appropriate values
		if ( first ) {
			data_storage::SmartSegment::become( seg1_substructure_[ sub_index ], current_segment );
			if ( current_segment->get_segment_id() == current_basis.segment_id() ) {
				TR.Debug << "SETTING UP SEG1" << std::endl;
				seg1_ = seg1_substructure_[ sub_index ];
			}
		} else {
			data_storage::SmartSegment::become( seg2_substructure_[ sub_index ], current_segment );
			if ( current_segment->get_segment_id() == current_basis.segment_id() ) {
				TR.Debug << "SETTING UP SEG2" << std::endl;
				seg2_ = seg2_substructure_[ sub_index ];
			}
		}
		if ( sub_index > 1 ) {
			if ( first ) {
				data_storage::SmartSegment::link_to( seg1_substructure_[ sub_index - 1 ], seg1_substructure_[ sub_index ] );
			} else {
				data_storage::SmartSegment::link_to( seg2_substructure_[ sub_index - 1 ], seg2_substructure_[ sub_index ] );
			}
		}
		current_segment = current_segment->get_c_terminal_neighbor();
		++sub_index;
	}


	numeric::HomogeneousTransform<core::Real> ht;
	if ( first ) {
		utility::vector1< core::conformation::Atom > & basis_atoms = seg1_->get_residue( current_basis.resnum() )->get_atom_vector();
		ht = numeric::HomogeneousTransform< core::Real > (
			numeric::xyzVector< core::Real >( basis_atoms[1].xyz() ),
			numeric::xyzVector< core::Real > ( basis_atoms[2].xyz() ),
			numeric::xyzVector< core::Real > (basis_atoms[3].xyz() ) );
	} else {
		utility::vector1< core::conformation::Atom > & basis_atoms = seg2_->get_residue( current_basis.resnum() )->get_atom_vector();
		ht = numeric::HomogeneousTransform< core::Real > (
			numeric::xyzVector< core::Real >( basis_atoms[1].xyz() ),
			numeric::xyzVector< core::Real > ( basis_atoms[2].xyz() ),
			numeric::xyzVector< core::Real > (basis_atoms[3].xyz() ) );
	}

	bool traverse_segments = true;
	core::Size curr_res = 1;
	core::Size curr_atom = 1;

	data_storage::SmartSegmentOP seg;
	if ( first ) {
		seg = seg1_substructure_[ 1 ];
	} else {
		seg = seg2_substructure_[ 1 ];
	}


	while ( traverse_segments ) {
		//= segment_vector_->at( current_basis.segment_id() )->get_residue( current_basis.resnum() )->get_const_atom_vector();
		numeric::xyzVector < core::Real > new_coords;

		for ( data_storage::SmartSewingResidueOP sew_res : seg->get_residue_vector() ) {
			curr_atom = 1;
			for ( core::conformation::Atom & atom : sew_res->get_atom_vector() ) {
				//TESTING maybe this was ok to start with
				new_coords = ht.to_local_coordinate( atom.xyz() );
				atom.xyz( numeric::xyzVector< core::Real > (new_coords ) );
				//TR << "After transformation in transform_segments, atom " << curr_atom << " of Residue " << curr_res << " of segment " << current_segment->get_segment_id() << " has coordinates " << new_coords.to_string() << std::endl;
				curr_atom++;
			}
			curr_res++;
		}
		if ( seg->is_c_terminus_fixed() ) {
			seg = seg->get_c_terminal_neighbor();
			++sub_index;
		} else {
			traverse_segments = false;
		}
	}
}


///@details Method to generate a HashKey for a SewAtom based on its coordinates.
///Essentially, this discretizes the coordinates into 0.25 Angstrom bins.
HashKey
Hasher::generate_key(
	core::conformation::Atom const & atom
) const {
	HashKey key;
	/*
	key[1] = boost::math::iround( atom.xyz().x()*4 );
	key[2] = boost::math::iround( atom.xyz().y()*4 );
	key[3] = boost::math::iround( atom.xyz().z()*4 );
	*/
	key[1] = std::round( atom.xyz().x()*4 );
	key[2] = std::round( atom.xyz().y()*4 );
	key[3] = std::round( atom.xyz().z()*4 );
	//TR << "Generated key " << key[ 1 ] << " " << key[ 2 ] << " " << key[ 3 ] << std::endl;
	return key;
}



///@Details For each of the transformed features, insert the appropriate
///key-value pair into the hash table. The key is the 3D-voxel (or bin)
///corresponding to the basis set. The value is the residues number for the
///residue that generated basis set that the transformed atom coordinates are in frame of
void
Hasher::hash_segments( data_storage::SmartSegmentOP const transformed_n_term_segment, data_storage::Basis const & basis_residue ) {
	//Mostly copied from old Hasher
	//For each set of transformed query_residues, count up and save hits to other models
	bool traverse_segments = true;
	data_storage::SmartSegmentOP current_segment = transformed_n_term_segment;
	while ( traverse_segments ) {
		TR.Debug << "Hashing segment " << current_segment->get_segment_id() << std::endl;

		core::Size segid = current_segment->get_segment_id();
		for ( core::Size res_num = 1; res_num < segment_vector_->at( segid )->get_length(); ++res_num ) {
			if ( current_segment->get_segment_id() == basis_residue.segment_id() && res_num == basis_residue.resnum() ) {
				continue; // don't hash the basis residue
			}
			utility::vector1< core::conformation::Atom > & atom_vector = current_segment->get_residue( res_num )->get_atom_vector();
			for ( core::Size atom_num = 1; atom_num <= 4; ++atom_num ) { //Not all of the atoms were being hashed; maybe this was the issue?
				HashValue value;
				value.basis_resnum = basis_residue.resnum();
				value.basis_segment_id = basis_residue.segment_id();
				value.segment_id = current_segment->get_segment_id();
				value.resnum = res_num;
				value.atomno = atom_num;
				HashKey const key = generate_key( atom_vector[ atom_num ] );
				hash_map_[key] = value;
			}
		}
		if ( current_segment->is_c_terminus_fixed() ) {
			current_segment = current_segment->get_c_terminal_neighbor();
		} else {
			traverse_segments = false;
		}

	}
}
/*
void
Hasher::copy_segment_neighbors( data_storage::SmartSegmentOP term_segment, utility::vector1< data_storage::SmartSegmentOP > & neighbor_vector ){
data_storage::SmartSegmentOP current_segment = term_segment;
//neighbor_vector.clear(); //no memory leaks for you, vector.
TR << "Cleared vector" << std::endl;
//AHA! We don't know for sure that "n_term_segment" is N-terminal. It might be C-terminal!
//int segment_index = 1;
int segment_index = term_segment->get_segment_id()+1;
int neighbor_vector_index = 1;
while( current_segment->is_c_terminus_fixed() ){
TR << "iterating through copy loop" << std::endl;
neighbor_vector[ neighbor_vector_index ]->become( segment_vector_[segment_index] );// store copy in vector
segment_index++;
neighbor_vector_index++;
current_segment = segment_vector_[segment_index];
}
neighbor_vector[ neighbor_vector_index ]->become( segment_vector_[segment_index] );// store copy in vector
}
*/

data_storage::SmartSegmentOP
Hasher::initialize_hashmap( data_storage::BasisPair const & basis_pair ) {
	core::Size starttime = time( NULL );
	//This uses the first model and first basis residue in the transform_model function
	//Get the two models

	// data_storage::SmartSegmentOP seg1 = segment_vector_[ basis_pair.first.segment_id() ];
	//data_storage::SmartSegmentOP seg2 = segment_vector_[ basis_pair.second.segment_id() ];

	//BEGIN INIT SECTION

	core::Size current_segment_id = 0;
	core::Size current_basis = 0;
	data_storage::SmartSegmentOP  segments_to_hash;
	data_storage::SmartSegmentOP  segments_to_score;
	if ( hash_map_.size() > 0 ) { //If the hash_map_ is not empty
		current_segment_id = (hash_map_.begin())->second.basis_segment_id;
		current_basis = (hash_map_.begin())->second.basis_resnum;
	}
	//If the first basis pair's segment model is already in the hashmap
	if ( current_segment_id == basis_pair.first.segment_id() && current_basis == basis_pair.first.resnum() ) {
		//basis 1 is already in the hash map, so we'll score model 2 against it.
		TR.Debug << "Basis 1 already in hash map" << std::endl;
		transform_segments( basis_pair.second, false );
		segments_to_score =  seg2_substructure_[ 1 ]; //First in vector will always be N-terminal
		//data_storage::SmartSegment::get_n_most_segment( segment_vector_[ basis_pair.second.segment_id() ] , false ); //basis pair might not be in the n_most segment
	} else if ( current_segment_id == basis_pair.second.segment_id() && current_basis == basis_pair.second.resnum() ) {
		//basis 2 is in the hash map, so we'll score model 1 against it.
		TR.Debug << "Model 2 already in hash map" << std::endl;
		transform_segments( basis_pair.first, true );
		segments_to_score =  seg1_substructure_[ 1 ];
		//data_storage::SmartSegment::get_n_most_segment( segment_vector_[ basis_pair.first.segment_id() ] , false ); //basis pair might not be in the n_most segment
	} else {
		hash_map_.clear();
		TR.Debug << "Transform both segments " << std::endl;
		TR.Debug << "Transforming segment " << basis_pair.first.segment_id() << " at residue " << basis_pair.first.resnum() << std::endl;
		transform_segments( basis_pair.first, true );
		TR.Debug << "Transformed segment. Preparing to hash" << std::endl;
		segments_to_hash = seg1_substructure_[ 1 ];
		//data_storage::SmartSegment::get_n_most_segment( segment_vector_[ basis_pair.first.segment_id() ] , false ); //basis pair might not be in the n_most segment
		hash_segments( segments_to_hash, basis_pair.first );
		TR.Debug << "Hashed basis 1. Scoring basis 2 against it" << std::endl;
		TR.Debug << "Transforming segment " << basis_pair.second.segment_id() << " at residue " << basis_pair.second.resnum() << std::endl;
		transform_segments( basis_pair.second, false );
		segments_to_score =  seg2_substructure_[ 1 ];
		//data_storage::SmartSegment::get_n_most_segment( segment_vector_[ basis_pair.second.segment_id() ] , false ); //basis pair might not be in the n_most segment
		TR.Debug << "Transformed segment." << std::endl;
	}
	core::Size endtime = time( NULL );
	TR << "Initialized hashmap in " << endtime - starttime << " seconds!" << std::endl;
	return segments_to_score;
}




///@details what this function does
bool
Hasher::score_basis_pair(
	data_storage::BasisPair const & basis_pair
) {
	core::Size starttime = time(NULL);
	//this will score the basis_pair
	//Assemblies assume each segment will occur at most once (if this changes later, comment out this statement!)
	// if( data_storage::SmartSegment::get_n_most_segment( segment_vector_[ basis_pair.first.segment_id() ], false ) == data_storage::SmartSegment::get_n_most_segment( segment_vector_[ basis_pair.second.segment_id() ], false ) ) {
	/*
	if( seg1_substructure_[ 1 ] == seg2_substructure_[ 1 ] ){
	TR << "These segments are the same segment apparently" << std::endl;
	return false;
	}
	*/
	//Initialize variables for clash score and hash score
	core::Size clash_score = 0;
	core::Size hash_score = 0;
	TR.Debug << "scoring " << basis_pair.first.segment_id() << " and " << basis_pair.second.segment_id() << std::endl;
	//Initialize the HashMap
	data_storage::SmartSegmentOP current_segment = initialize_hashmap( basis_pair );
	//The hash_map_ now contains one of the models, and the other must be scored against it.
	//TR << "HashMap initialized!" << std::endl;
	//TR << "Current size of HashMap: " << hash_map_.size() << std::endl;
	//Now we need to find the scores.
	//bool traverse_segments = true;
	//data_storage::SmartSegmentOP current_segment = segments_to_score;

	while ( current_segment != nullptr ) {
		core::Size segid = current_segment->get_segment_id();
		for ( core::Size res_num = 1; res_num <= segment_vector_->at( segid )->get_length(); ++res_num ) {
			//if ( current_segment->get_segment_id() == basis_pair.first.segment_id() && res_num == basis_pair.first.resnum() ) {
			//  continue; // don't score the basis residue
			//}
			//else if ( current_segment->get_segment_id() == basis_pair.second.segment_id() && res_num == basis_pair.second.resnum() ) {
			//  continue; // don't score the basis residue
			//}
			//TR << "Success! This isn't the basis residue" << std::endl;
			utility::vector1< core::conformation::Atom > & atom_vector = current_segment->get_residue( res_num )->get_atom_vector();
			//TR << "Got atom vector!" << std::endl;
			for ( core::Size atom_num = 1; atom_num <= 4; ++atom_num ) {
				//TR << "During alignment in score_basis_pair, atom " << atom_num << " of residue " << res_num << " of segment " << current_segment->get_segment_id() << " has coordinates " << atom_vector[atom_num].xyz().to_string() << std::endl;
				//TR << "Checking score of atom " << atom_num << std::endl;
				std::set< HashValue > aligned_atoms = find_aligned_atoms( atom_vector[ atom_num ] );
				//Iterate through the aligned atoms
				for ( std::set< HashValue >::iterator hashit = aligned_atoms.begin(); hashit != aligned_atoms.end(); ++hashit ) {
					//If these are the same atom type, they could potentially superimpose
					if ( hashit->atomno == atom_num ) {
						//TR << "Atom numbers " << hashit->atomno <<" of residue " << res_num << " align" << std::endl;
						//Only count as a good hash if:
						//1) the user has specified that they don't want to restrict scoring to
						//opposite terminal segments, or
						//2) The atoms are in opposite terminal segments (can_hash( seg1, seg2) )
						data_storage::SmartSegmentCOP hashmap_segment;
						if ( pdb_segments_.count( hashit->segment_id ) == 0 ) {
							hashmap_segment = segment_vector_->at( hashit->segment_id  );
						} else {
							hashmap_segment = pdb_segments_.at( hashit->segment_id );
						}
						//All we actually do with this segment is check if the two are hashable, so it's ok
						if ( can_hash( current_segment, hashmap_segment ) || !hasher_settings_.hash_between_termini ) {
							hash_score = hash_score + 1;
						}//end if hash
					} else { //end if these are the same atom //If these are different atoms
						//TR << "Atom numbers " << hashit->atomno << " and " << atom_num << " of residue " << res_num << " clash" << std::endl;
						clash_score = clash_score + 1;
					}//end if these are different atoms
				}//End for loop
				//TR << "Done with this atom" << std::endl;
			}//End iterator over model atoms
		}//End iterator over residues
		current_segment = current_segment->get_c_terminal_neighbor();
	}
	core::Size endtime = time( NULL );
	TR << "Scored basis pair " << basis_pair.first.segment_id() << " " << basis_pair.first.resnum() << ", " << basis_pair.second.segment_id() << " " << basis_pair.second.resnum() << " in " << endtime - starttime << " seconds!" << std::endl;
	TR << "Clash score " << clash_score << " and hash score " << hash_score << std::endl;
	return ( clash_score <= hasher_settings_.max_clash_score && hash_score >= hasher_settings_.min_hash_score );

}


std::set< HashValue >
Hasher::find_aligned_atoms( core::conformation::Atom const & start_atom ) const{
	//TR << "Current HashMap size: " << hash_map_.size() << std::endl;
	//Generate a key for that atom
	//core::conformation::Atom cur_atom = start_atom;
	HashKey const cur_key = generate_key ( start_atom );
	HashKey modified_key;
	//Determine size of box from hasher settings
	int limit = ( hasher_settings_.boxes_per_dimension  - 1 )/ 2;
	int lower_limit = -1*limit;
	//TR << "Checking between limits " << lower_limit << " and " << limit << std::endl;
	//Initialize the aligned atoms set
	std::set< HashValue > aligned_atoms;
	//aligned_atoms.clear();

	//We'll actually have (boxes_per_dimension)**3 keys to look up for each atom.
	//In this case we could just use a vector1 of the hash keys and look up the hash values from the map
	//directly since we changed it to only have one HashValue per HashKey
	for ( int i=lower_limit; i <= limit; ++i ) {
		for ( int j=lower_limit; j <= limit; ++j ) {
			for ( int k=lower_limit; k <= limit; ++k ) {
				modified_key[ 1 ] = cur_key[ 1 ] + i;
				modified_key[ 2 ] = cur_key[ 2 ] + j;
				modified_key[ 3 ] = cur_key[ 3 ] + k;
				if ( hash_map_.find( modified_key ) != hash_map_.end() ) {
					aligned_atoms.insert( hash_map_.find( modified_key )->second );
					//TR << "Found aligned atom!" << std::endl;
				}//End if
			}//end k loop
		}//end j loop
	}//end i loop
	return aligned_atoms;
} //end find_aligned_atoms



///@details what this function does
std::pair<  utility::vector1< std::pair< core::Size, core::Size > >::const_iterator,
utility::vector1< std::pair< core::Size, core::Size > >::const_iterator >
Hasher::iterate_over_basis_pairs(
	data_storage::SmartSegmentCOP segment_1,
	data_storage::SmartSegmentCOP segment_2
)
{
	//Fill in our vector
	all_basis_pairs_.clear();

	//If either of the segments is not hashable, or if they are different SS types, leave the vector empty.
	if ( !( segment_1->is_hashable() && segment_2->is_hashable()
			&& ( segment_1->get_dssp_code() == segment_2->get_dssp_code() ) )
			) {
		//in this case, we will leave all_basis_pairs_ empty.
		TR << "Cannot hash segments!" << std::endl;
		return std::make_pair( all_basis_pairs_.end(), all_basis_pairs_.end() );
	}//endif

	//Go ahead and calculate how many residues forward and back we'll have to go to get the desired overlap.
	//Divide the # atoms by 4 and round up. Atoms in the basis residues don't count.
	//This is a core::Real for division purposes (no integer division)
	const core::Real atoms_per_residue = 4.0;
	core::Size required_num_residues = core::Size ( std::ceil( hasher_settings_.min_hash_score / atoms_per_residue ) );
	if ( required_num_residues == 0 ) {
		required_num_residues = 1; //1 residue must always be available as basis even if nothing is required to align
	}
	TR.Debug << "Producing basis pairs within " << utility::to_string( required_num_residues ) << " residues of termini" << std::endl;

	//Find number of residues in each segment
	core::Size seg1_size = segment_1->get_length();
	core::Size seg2_size = segment_2->get_length();
	TR.Debug << "seg1_size " << seg1_size << " seg2_size " << seg2_size << std::endl;
	//If either segment is shorter than required_num_residues, no basis pairs possible. Since basis residues don't count, use <=.
	if ( seg1_size <= required_num_residues || seg2_size <= required_num_residues ) {
		TR << "One of the segments is too short!" << std::endl;
		return std::make_pair( all_basis_pairs_.end(), all_basis_pairs_.end() );
	}
	//We need to check that one segment is N-terminal and one segment is C-terminal, but we don't really need to know which is which.
	core::Size i_start = 1;
	core::Size i_end = ( seg1_size + 1 - required_num_residues );
	core::Size j_start = 1;
	core::Size j_end = (seg2_size - required_num_residues + 1 );
	if ( can_hash( segment_1, segment_2 ) ) {
		TR.Debug <<" Hashing segments" << std::endl;
		for ( core::Size i = i_start; i <= i_end; ++i ) {
			for ( core::Size j = j_start; j <= j_end; ++j ) {
				//Only insert basis pairs that could possibly produce an overlap of the size we want
				//For instance, if we're using the last residue in segment 1, then we must be at least at required_num_residues + 1
				//(the basis residues don't count toward the overlap) from the end of segment 2

				//Compare the number of possible overlap residues C-terminal of i with the number N-terminal of j.
				//To test this, figure out min (how many residues in seg1 are N-terminal of i, how many residues in seg2 are N-terminal of j)
				//+ min(how many in seg2 are C-terminal of j, how many in seg1 are C-terminal of i )
				//This will be the size of the largest possible overlap for a given alignment.
				//Compare this to required_num_residues.
				TR.Debug << "Basis residues " << i << " and " << j << std::endl;
				core::Size max_overlap_size = utility::min( i - 1, j - 1 ) + utility::min ( seg1_size - i, seg2_size - j );
				TR.Debug << "maximum detected overlap " << max_overlap_size << std::endl;
				if ( max_overlap_size >= required_num_residues ) {
					TR.Debug << "Adding basis pair" << std::endl;
					all_basis_pairs_.push_back( std::make_pair( i, j) );
				}//end if
			}//end for j
		}//end for i
	}//end if
	//this will give you some iterators, yay!
	return std::make_pair( all_basis_pairs_.begin(), all_basis_pairs_.end() );
}//end iterate_over_basis_pairs function



bool
Hasher::can_hash( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2 ){
	//Minnie's updated SmartSewing version (9/16/2016), still only hashing terminal segments
	//if a terminus is not fixed, that means the segment IS that terminus
	//If seg1 is the c terminus and seg2 is the N terminus
	TR.Debug << "Checking if we can hash these segments . . . " << std::endl;
	return ( !(seg1->is_c_terminus_fixed() || seg2->is_n_terminus_fixed()) // seg1 cterm and seg2 nterm are free
		// or seg2 is the c terminus and seg1 is the n terminus
		||(!seg2->is_c_terminus_fixed() && !seg1->is_n_terminus_fixed()) // seg2 cterm and seg1 nterm are free
	);
}//can_hash


//Debugging purposes only
core::Size
Hasher::get_hash_map_size() const{
	return hash_map_.size();
}

data_storage::SmartSegmentCOP
Hasher::get_seg1() const{
	return seg1_;
}

data_storage::SmartSegmentCOP
Hasher::get_seg2() const{
	return seg2_;
}




} //hashing
} //protocols
} //sewing

