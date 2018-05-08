// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/AlignmentGenerator.hh
/// @brief Determines all possible alignments for an edge given the two SewSegments involved. Also contains the mechanisms for reading an edge file.
/// @author guffysl (guffy@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_hashing_AlignmentGenerator_hh
#define INCLUDED_protocols_sewing_hashing_AlignmentGenerator_hh

#include <protocols/sewing/hashing/AlignmentGenerator.fwd.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/sewing/data_storage/SmartSegment.fwd.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sewing {
namespace hashing {

/// @brief This class contains the infrastructure for parsing an edge file and for generating all possible alignments for a given edge.
class AlignmentGenerator : public utility::pointer::ReferenceCount {

public:

	AlignmentGenerator();
	AlignmentGenerator( HasherSettings const & hasher_settings, SegmentVectorCOP segment_vector );
	AlignmentGenerator( HasherSettings const & hasher_settings, SegmentVectorCOP segment_vector, std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs );
	AlignmentGenerator(AlignmentGenerator const & src);

	virtual ~AlignmentGenerator();

	///@brief Determines all possible basis pairs for an edge given the two SewSegments involved.
	utility::vector1< std::pair< core::Size, core::Size > > get_all_alignments( data_storage::SmartSegmentCOP s1, data_storage::SmartSegmentCOP s2 );

	//What else will AlignmentGenerator need to be able to do?
	//Somewhere we'll need an atom map for each basis pair, but we won't really use that until Assembly--maybe Assembly could
	//just get it directly from the Hasher?

	AlignmentGeneratorOP
	clone() const;


	// static EdgeContainer read_edge_file( std::string edge_file_name );


	//Getters
	SegmentVectorCOP segment_vector() const;
	HasherSettings hasher_settings() const;
	std::map< core::Size, data_storage::SmartSegmentOP > & pdb_segments();
	std::map< core::Size, data_storage::SmartSegmentOP > const_pdb_segments() const;
	//Setters
	void segment_vector( SegmentVectorCOP segvec );
	void hasher_settings( HasherSettings & hs );
	void pdb_segments( std::map< core::Size, data_storage::SmartSegmentOP > & );
private:
	//We need the hasher settings for the Hasher that will generate the alignments
	HasherSettings hasher_settings_;
	//The Hasher takes a basis pair and needs the actual segments for its initialization
	SegmentVectorCOP segment_vector_;
	std::map< core::Size, data_storage::SmartSegmentOP > pdb_segments_;

};

} //hashing
} //sewing
} //protocols



#endif //INCLUDED_protocols_sewing_hashing_AlignmentGenerator_hh





