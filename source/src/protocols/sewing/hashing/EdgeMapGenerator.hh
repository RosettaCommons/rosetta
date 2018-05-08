// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/EdgeMapGenerator.hh
/// @brief Class used to generate edge files for use with SEWING protocols.
/// @author guffysl (guffy@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_hashing_EdgeMapGenerator_hh
#define INCLUDED_protocols_sewing_hashing_EdgeMapGenerator_hh

#include <protocols/sewing/hashing/EdgeMapGenerator.fwd.hh>

//#include <protocols/sewing/conformation/Model.fwd.hh>
#include <protocols/sewing/data_storage/SmartSegment.fwd.hh>
#include <protocols/sewing/hashing/hasher_data.hh> //Need this for HasherSettings
#include <protocols/sewing/data_storage/SmartSewingResidue.fwd.hh>

//Pretty sure this one doesn't actually exist yet . . .
#include <protocols/sewing/hashing/Hasher.fwd.hh>
//Core headers
#include <core/types.hh>
#include <map>
#include <string>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.fwd.hh>
namespace protocols {
namespace sewing {
namespace hashing {

///@brief Class used to generate and load in edge files for use with SEWING protocols.
class EdgeMapGenerator : public utility::pointer::ReferenceCount {

public:

	EdgeMapGenerator();
	EdgeMapGenerator(std::string model_file_name, core::Size boxes_per_dimension, core::Size min_hash_score, core::Size max_clash_score, std::string edge_file_name, bool hash_opposite_termini );
	EdgeMapGenerator(EdgeMapGenerator const & other);
	EdgeMapGenerator(std::string edge_file_name );

	virtual ~EdgeMapGenerator();

	EdgeMapGeneratorOP
	clone();

	/// @brief Check whether two SewSegments can form an edge.
	bool edge_exists(data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2, HasherOP hasher);

	/// @brief Method called by each individual processor to write to its edge file
	void write_edges_for_segment(core::Size current_n_segment, HasherOP hasher, utility::io::ozstream & output_file, utility::vector1< core::Size > const & c_segs );

	/// @brief Generate vectors of IDs of N-terminal and C-terminal segments
	std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > get_terminal_segments();
	///@brief Used to add a given segment's n-most and c-most terminal segments to EdgeMap
	void append_to_edge_map( core::Size base_seg_id );

	/// @brief Main method of this class that is responsible for dividing work between nodes and iterating over models
	//Note--since there now ARE NO MODELS, this logic will also need to change.
	void generate_edge_file();
	/// @brief Get values of command-line options
	void initialize_from_options();
	/// @brief Shared setup of private data members
	void set_private_data( std::string model_file_name, core::Size boxes_per_dimension, core::Size min_hash_score, core::Size max_clash_score, bool hash_opposite_termini );

	//Getters
	EdgeMapOP edge_map() const;
	core::Size version() const;
	std::string edge_file_name() const;
	HasherSettings hasher_settings() const;
	std::map< core::Size, data_storage::SmartSegmentOP > &
	pdb_segments();
	std::map< core::Size, data_storage::SmartSegmentOP > const
	const_pdb_segments() const;


	// utility::pointer::shared_ptr< std::map< int, Model > > model_map();
	SegmentVectorCOP segment_vector() const; //Vector of pointers, so not too bad to throw around
	bool is_edge_file_loaded() const;

	//Setters

	void version( core::Size ver);
	void edge_file_name( std::string name);
	//void model_map( utility::pointer::shared_ptr< std::map< int, Model > >  modelmap);
	void segment_vector( SegmentVectorCOP segvector );
	void hasher_settings( HasherSettings & hs );
	void
	pdb_segments(  std::map< core::Size, data_storage::SmartSegmentOP > & );

	///@brief Main function for reading edge files (sets values of other member variables)
	void load_edge_file( std::string edge_file_name ); // all other data members come from the edge file


	///@brief Alternate function to read edge file when name is already known
	void read_edge_file();

	///@brief function returns a set of segments that the provided segment can connect to.
	EdgeMap::iterator get_edges( core::Size );

private:
	HasherSettings hasher_settings_;
	//Each thread will need its own Hasher object with shared HasherSettings
	core::Size version_=0;
	//utility::pointer::shared_ptr< std::map< int, Model > > model_map_;
	EdgeMapOP edges_ = EdgeMapOP( new EdgeMap );

	SegmentVectorCOP segment_vector_;
	std::string edge_file_name_;
	bool is_edge_file_loaded_;
	std::map< core::Size, data_storage::SmartSegmentOP > pdb_segments_;
};

} //hashing
} //sewing
} //protocols



#endif //INCLUDED_protocols_sewing_hashing_EdgeMapGenerator_hh






