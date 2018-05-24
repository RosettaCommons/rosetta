// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/hashing/Hasher.hh
/// @brief A geometric hashing class used by the SEWING protocol
/// @author Minnie Langlois (minnie@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_hashing_Hasher_hh
#define INCLUDED_protocols_sewing_hashing_Hasher_hh


//Protocols headers
#include <protocols/sewing/hashing/Hasher.fwd.hh>
#include <protocols/sewing/data_storage/Basis.hh>
#include <protocols/sewing/data_storage/SmartSegment.fwd.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.fwd.hh>
#include <protocols/sewing/hashing/hasher_data.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

//Core headers
#include <core/types.hh>
#include <core/conformation/Atom.hh>


namespace protocols {
namespace sewing {
namespace hashing {

///@brief A geometric hashing class used by the SEWING protocol
class Hasher : public utility::pointer::ReferenceCount {

public:

	// constructors
	Hasher();

	virtual ~Hasher();

	Hasher( HasherSettings & hasher_settings, SegmentVectorCOP segment_vector );

	Hasher( HasherSettings & hasher_settings, SegmentVectorCOP segment_vector, std::map< core::Size, data_storage::SmartSegmentOP > pdb_segments );

	// copy constructor
	Hasher(Hasher const & src);

	HasherOP clone() const;


	bool can_hash( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2 );

	void initialize_from_options();

	///@brief returns pair of iterators to beginning and end of all_basis_pairs
	///between the given SewSegments
	std::pair< utility::vector1< std::pair< core::Size, core::Size > >::const_iterator,
	utility::vector1< std::pair< core::Size, core::Size > >::const_iterator >
	iterate_over_basis_pairs( data_storage::SmartSegmentCOP segment_1, data_storage::SmartSegmentCOP segment_2 );


	//Getters & Setters
	SegmentVectorCOP segment_vector() const;

	void segment_vector( SegmentVectorCOP segment_vector );

	HasherSettings hasher_settings() const;

	void hasher_settings( HasherSettings const & hasher_settings);

	std::map< core::Size, data_storage::SmartSegmentOP > &
	pdb_segments();
	std::map< core::Size, data_storage::SmartSegmentOP >
	const_pdb_segments() const;
	void
	pdb_segments(  std::map< core::Size, data_storage::SmartSegmentOP > & );

	//AtomMap get_atom_map( data_storage::BasisPair & basis_pair );


	//Hashing and Transformations
	///@brief iterates through neighbors of n_terminal segments transforms xyz coordninates into basis_residue coordinate frame.
	void transform_segments( data_storage::Basis current_basis, bool first );

	///@brief adds segments to the hash_map, excluding the basis residue from being scored.
	void hash_segments ( data_storage::SmartSegmentOP const transformed_segment, data_storage::Basis const & basis_residue );

	///@brief creates key for HashMap
	HashKey generate_key( core::conformation::Atom const & atom ) const;

	///@brief Given a ModelIterator to a particular atom, returns a set of HashValues that correspond to it.
	std::set< HashValue > find_aligned_atoms( core::conformation::Atom const & start_atom ) const;

	///@brief Returns a transformed segments that will be scored against the model currently in the HashMap.
	data_storage::SmartSegmentOP initialize_hashmap( data_storage::BasisPair const & basis_pair );

	///@brief returns whether the basis pair meets edge requirements for max_clash
	///and min_hash_score
	bool score_basis_pair( data_storage::BasisPair const & basis_pair);

	core::Size get_hash_map_size() const;
	data_storage::SmartSegmentCOP get_seg1() const;
	data_storage::SmartSegmentCOP get_seg2() const;




private:
	///@brief iterates through all neighbors of the n_term segment and copies them into preallocated vector via the data_storage::SmartSegment become() method.
	//void copy_segment_neighbors( data_storage::SmartSegmentOP n_term_segment, utility::vector1< data_storage::SmartSegmentOP > & neighbor_vector );
	///@brief initializes preallocated segments and respective neighbor vectors
	void allocate_segments();

	HasherSettings hasher_settings_;
	utility::vector1< std::pair< core::Size, core::Size > > all_basis_pairs_;
	SegmentVectorCOP segment_vector_;
	HashMap hash_map_;
	//pre-allocated smartsegments for hashmap
	data_storage::SmartSegmentOP seg1_;
	data_storage::SmartSegmentOP seg2_;
	utility::vector1< data_storage::SmartSegmentOP > seg1_substructure_;
	utility::vector1< data_storage::SmartSegmentOP > seg2_substructure_;
	std::map< core::Size, data_storage::SmartSegmentOP > pdb_segments_;
};


} //hashing
} //protocols
} //sewing



#endif //INCLUDED_protocols_sewing_hashing_Hasher_hh





