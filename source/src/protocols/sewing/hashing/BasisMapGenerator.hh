// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/BasisMapGenerator.hh
/// @brief Given a model file, edge file ,and one or more input structures, generates alignment files for use with AppendAssemblyMover.
/// @author guffysl (guffy@email.unc.edu)
/// @author minniel (minnie@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_hashing_BasisMapGenerator_hh
#define INCLUDED_protocols_sewing_hashing_BasisMapGenerator_hh

// Unit headers
#include <protocols/sewing/hashing/BasisMapGenerator.fwd.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/sewing/data_storage/Basis.hh>
// Project Headers
#include <protocols/sewing/hashing/AlignmentGenerator.fwd.hh>
#include <protocols/sewing/hashing/EdgeMapGenerator.fwd.hh>

namespace protocols {
namespace sewing {
namespace hashing {



///@brief Given a model file, edge file ,and one or more input structures, generates alignment files for use with AppendAssemblyMover.
class BasisMapGenerator {

public:

	BasisMapGenerator();

	// copy constructor
	BasisMapGenerator( BasisMapGenerator const & src );


	BasisMapGenerator(
		EdgeMapGeneratorOP edge_file_reader,
		SegmentVectorCOP segments );
	BasisMapGenerator(
		EdgeMapGeneratorOP edge_file_reader,
		std::string model_file_name,
		std::string alignment_file_name
	);

	BasisMapGenerator(
		EdgeMapGeneratorOP edge_file_reader,
		std::string model_file_name
	);

	BasisMapGenerator(
		std::string model_file_name,
		std::string edge_file_name,
		utility::vector1< core::Size > match_segments,
		utility::vector1< core::Size > pose_segment_starts,
		utility::vector1< core::Size > pose_segment_ends,
		core::Size recursive_depth);

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~BasisMapGenerator();

	//Helper methods for constructors
	void initialize_from_options();
	void set_private_data( std::string model_file_name, std::string edge_file_name, utility::vector1< core::Size >  match_segments, utility::vector1< core::Size >pose_segment_starts, utility::vector1< core::Size > pose_segment_ends, core::Size recursive_depth );
	void write_alignments_to_file( std::string file_name ) const;
	void edge_file_from_options();
	void model_file_from_options();
	void alignment_settings_from_options();
	void set_model_file(std::string model_file_name);
	void set_edge_file(std::string edge_file);
	void set_edge_file_reader(EdgeMapGeneratorOP edge_file_reader);
	void set_alignment_settings(
		utility::vector1< core::Size > match_segments,
		utility::vector1< core::Size > pose_segment_starts,
		utility::vector1< core::Size > pose_segment_ends,
		core::Size recursive_depth);

	//Getters
	core::Size version() const;
	// std::string alignment_file_name() const;
	HasherSettings hasher_settings() const;
	SegmentVectorCOP segment_vector() const;
	EdgeMapGeneratorOP edge_file_reader() const;
	AlignmentSettings alignment_settings() const;
	AlignmentGeneratorOP alignment_generator() const;
	BasisMapOP basis_map() const;
	std::map< core::Size, data_storage::SmartSegmentOP > &
	pdb_segments();

	std::map< core::Size, data_storage::SmartSegmentOP >
	const_pdb_segments() const;


	BasisMap::iterator get_alignments( data_storage::SmartSegmentOP query_segment );
	void delete_segment_basis_pairs( core::Size );

	//Setters
	void version( core::Size ver);
	void segment_vector( SegmentVectorCOP segvec);
	void hasher_settings( HasherSettings & hs );
	void alignment_generator( AlignmentGeneratorOP ag );
	void alignment_settings( AlignmentSettings & as );
	void basis_map( BasisMapOP bm );
	void
	pdb_segments(  std::map< core::Size, data_storage::SmartSegmentOP >);
	/// @brief required in the context of the parser/scripting scheme
	BasisMapGeneratorOP
	clone() const;

	void recurse_over_segments( core::Size current_depth, std::set< data_storage::SmartSegmentCOP > segments_to_check);
	void read_alignments_from_file( std::string file_name );
	///@brief Given a segmentID, return a randomly chosen basis pair from the BasisMap for this segment
	data_storage::BasisPair
	get_basis_pair_for_segment( core::Size segment_id );


	///@brief This function wraps several of the functions below to retrieve a basis pair for local (possibly chimaeric)
	///input segments.
	///@details This function first calls segs_can_chimerize to determine whether the segments can ever chimerize
	///It then retrieves the basis residue bounds for the two segmnets
	///If seg1 (or seg2) is chimaeric, it gets the appropriate parent to use for finding a basis pair. It also handles
	///all conversions between parent and chimaeric residue numbering
	///It has 1000 chances to choose a valid basis pair for the chimaera that fits within the residue bounds
	///If it hasn't found something by then, it just gives up
	std::pair< bool, data_storage::BasisPair >
	get_basis_pair_for_local_segments( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2, bool first_is_N );

	///@brief Given two (possibly chimaeric) segments, return ranges of valid basis residues for each
	///@details Accounds for both chimaeric basis residues and vital residues. Returns pairs (resnum ranges) in
	///the same order that the segments were provided.
	std::pair< std::pair< core::Size, core::Size >, std::pair< core::Size, core::Size > >
	get_basis_res_bounds( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2, bool first_is_N );

	///@brief Given two (original, non-chimaeric) segments, return a possible basis pair for those two segments.
	///@details Does not take chimaeric basis residue bounds into account; currently we handle those by first checking
	///if the chimaeric segment *can* chimerize with the partner and then call this repeatedly until we get something that works
	data_storage::BasisPair
	get_basis_pair_for_segment_pair( core::Size seg1, core::Size seg2, bool first_is_N );

	///@brief Checks that a pair of segments (including chimaerae!) can chimerize.
	///@details First bool indicates whether seg1 can ever be chimerized with ANYTHING;
	///second indicates if it can chimerize with seg2 specifically.
	///Calls get_basis_res_bounds to handle chimaera( only checks on the correct terminal-side of the old basis res for chimaerae
	std::pair< bool, bool>
	segs_can_chimerize( data_storage::SmartSegmentCOP seg1, data_storage::SmartSegmentCOP seg2, bool first_is_N );

	///@brief Checks that the given basis pair exists in the BasisMap after calling recurse_over_segments on seg1
	///@details Must be called for original segment IDs (no chimaerae or segment duplicates)
	bool
	basis_pair_is_valid( data_storage::BasisPair basis_pair, bool first_is_N );























private:
	HasherSettings hasher_settings_;
	core::Size version_;
	SegmentVectorCOP segment_vector_;
	EdgeMapGeneratorOP edge_file_reader_;
	AlignmentSettings alignment_settings_;
	BasisMapOP basis_map_;
	AlignmentGeneratorOP alignment_generator_;
	std::map< core::Size, data_storage::SmartSegmentOP > pdb_segments_;

};

//std::ostream &operator<< (std::ostream &os, BasisMapGenerator const &mover);


} //hashing
} //sewing
} //protocols


#endif //protocols_sewing_hashing_BasisMapGenerator_hh







