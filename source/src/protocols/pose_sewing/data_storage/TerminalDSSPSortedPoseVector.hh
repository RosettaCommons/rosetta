// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/TerminalDSSPSortedPoseVector.hh
/// @brief a vector of Poses sorted by the secondary structure of the termini
/// @author frankdt (frankdt@email.unc.edu)



#ifndef INCLUDED_protocols_pose_sewing_data_storage_TerminalDSSPSortedPoseVector_hh
#define INCLUDED_protocols_pose_sewing_data_storage_TerminalDSSPSortedPoseVector_hh

#include <protocols/pose_sewing/data_storage/TerminalDSSPSortedPoseVector.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

#include <core/pose/Pose.hh>
#include <core/simple_metrics/SimpleMetric.fwd.hh>

#include <protocols/pose_sewing/data_storage/SegmentEnvelope.fwd.hh>
#include <protocols/pose_sewing/data_storage/PoseSegment.fwd.hh>
#include <protocols/pose_sewing/data_storage/PoseWithTerminalSegmentsOfKnownDSSP.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <map>

namespace protocols {
namespace pose_sewing {
namespace data_storage {

/// @brief a vector of Poses sorted by the secondary structure of the termini
class TerminalDSSPSortedPoseVector : public utility::VirtualBase {

public:

	TerminalDSSPSortedPoseVector();
	TerminalDSSPSortedPoseVector(TerminalDSSPSortedPoseVector const & src);

	virtual ~TerminalDSSPSortedPoseVector();

	TerminalDSSPSortedPoseVectorOP
	clone() const;

	void
	parse_motif_file(std::string motif_file);

	void
	add_segment_envelopes_from_string(std::string motifs);

	void
	populate_from_pdb_list(std::string pdb_list);

	void
	populate_from_pdb_list(std::string motif_file, std::string pdb_list);

	void
	store_pose(PoseSegmentOP,PoseSegmentOP);

	void
	store_pose(PoseWithTerminalSegmentsOfKnownDSSPOP);

	PoseWithTerminalSegmentsOfKnownDSSPOP
	get_random_pose(std::string,std::string);

	utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>
	get_vector(std::string nterm_secstructs,std::string cterm_secstructs);

	void
	set_segment_envelopes(utility::vector1<utility::vector1<SegmentEnvelopeOP>>);

	void
	add_to_segment_envelopes(utility::vector1<SegmentEnvelopeOP> new_envelope);

	void
	clear_segment_envelopes();

	void
	populate_from_segment_file(std::string, bool clear = false);

	///@brief Read the segment file, return a vector of segments of max size.
	///@details This function also returns a random permutation of the total set
	///   n_ter and c_ter control which SS can be possible for N and C term.
	///   This function is to reduce the memory footprint associated with whole PDB segment reads.
	utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>
	populate_from_segment_file_and_get_random_vector_set(
		utility::vector1< std::string > const & paths,
		std::string const & allowed_n_ter_ss,
		std::string const & allowed_c_ter_ss,
		core::Size max_vector_size = 5000
	) const ;

	///@brief Read the segment file, return a vector of segments of max size.
	///@details This function also returns a random permutation of the total set
	///   n_ter and c_ter control which SS can be possible for N and C term.
	///   This function is to reduce the memory footprint associated with whole PDB segment reads.
	utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>
	populate_from_segment_file_and_get_random_vector_set(
		std::string const & path,
		std::string const & allowed_n_ter_ss,
		std::string const & allowed_c_ter_ss,
		core::Size max_vector_size = 5000
	) const ;

	void
	create_segment_file(std::string );

	void
	populate_all_OP_vector();

	utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>
	get_all_poseOPs();

	void
	clear_all_vectors();

	void
	clear_terminal_vectors();

	void
	purge_missing_density();

	///@brief Write segment file using a pdb list.  Append segment file if needed.
	void
	simultaneously_populate_and_write_segment_file (
		std::string const & motif_filename,
		std::string const & pdb_list,
		std::string const & file_path,
		bool store_reference_pdbs,
		bool use_pdbs_only,
		utility::vector1< protocols::filters::FilterOP > const & filters,
		utility::vector1< core::simple_metrics::SimpleMetricOP > const & metrics,
		bool output_elements = false);

	///@brief Write segment file on an individual pose. Append segment file if needed.
	void
	simultaneously_populate_and_write_segment_file_pdb (
		core::pose::Pose const & pose,
		std::string const & motif_filename,
		std::string const & file_path,
		bool store_reference_pdbs,
		bool use_pdbs_only,
		utility::vector1< protocols::filters::FilterOP > const & filters,
		utility::vector1< core::simple_metrics::SimpleMetricOP > const & metrics,
		bool output_elements = false);

	void
	simultaneously_populate_and_write_segment_file_pdb_using_stored_motifs (
		core::pose::Pose const & pose,
		std::string const & file_path,
		bool store_reference_pdbs,
		bool use_pdbs_only,
		utility::vector1< protocols::filters::FilterOP > const & filters,
		utility::vector1< core::simple_metrics::SimpleMetricOP > const & metrics,
		bool output_elements);

	core::Size
	get_total_size();

	utility::vector1<core::pose::PoseCOP>
	convert_pose_into_segment_vector(core::pose::PoseOP current_pose);

private:

	utility::vector1<utility::vector1<SegmentEnvelopeOP>> segment_envelopes_;

	utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP> all_poseOPs_;

	std::map< char, std::map<char,utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>>> poseOP_map_;

	std::string nterm_secstructs_;
	std::string cterm_secstructs_;

};


} //protocols
} //pose_sewing
} //data_storage



#endif //INCLUDED_protocols_pose_sewing_data_storage_TerminalDSSPSortedPoseVector_hh





