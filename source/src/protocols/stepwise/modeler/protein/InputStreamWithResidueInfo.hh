// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file InputStreamWithResidueInfo.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_InputStreamWithResidueInfo_HH
#define INCLUDED_protocols_stepwise_InputStreamWithResidueInfo_HH

#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <map>

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {

	void
	initialize_input_streams(  	utility::vector1< protocols::stepwise::modeler::protein::InputStreamWithResidueInfoOP > & input_streams );

	void
	initialize_input_streams_with_residue_info( utility::vector1< InputStreamWithResidueInfoOP > & input_streams_with_residue_info,
																							utility::vector1< std::string > const & pdb_tags,
																							utility::vector1< std::string > const & silent_files_in,
																							utility::vector1< Size > const & input_res,
																							utility::vector1< Size > const & input_res2
																							);

	import_pose::pose_stream::PoseInputStreamOP
	setup_pose_input_stream(
													utility::options::StringVectorOption const & option_s1,
													utility::options::StringVectorOption const & option_silent1,
													utility::options::StringVectorOption const & option_tags1	);


	//////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////
	class InputStreamWithResidueInfo:public utility::pointer::ReferenceCount {
	public:

		InputStreamWithResidueInfo( import_pose::pose_stream::PoseInputStreamOP pose_input_stream,
																utility::vector1< Size > const & input_res,
																utility::vector1< Size > const & slice_res );

		InputStreamWithResidueInfo( import_pose::pose_stream::PoseInputStreamOP pose_input_stream,
																utility::vector1< Size > const & input_res);

		~InputStreamWithResidueInfo();

		import_pose::pose_stream::PoseInputStreamOP & pose_input_stream();
		utility::vector1< Size > const & input_res();
		utility::vector1< Size > const & slice_res();
		std::map< Size, Size > & full_to_sub();

		void set_slice_res( 	utility::vector1< Size > const & slice_res );
		void set_full_to_sub( std::map< Size, Size > const & full_to_sub );
		void set_rsd_set( chemical::ResidueTypeSetCAP & rsd_set );

		void reset();

		bool has_another_pose() const;

		void
		advance_to_next_pose_segment();

		void copy_next_pose_segment( pose::Pose & pose );

		void copy_next_pose_segment( pose::Pose & pose,
																 pose::Pose & import_pose,
																 bool const check_sequence_matches,
																 bool const align_pose_to_import_pose = false );

		void
		apply_current_pose_segment( pose::Pose & pose );

		void
		apply_current_pose_segment( pose::Pose & pose,
																pose::Pose & import_pose,
																bool const check_sequence_matches,
																bool const align_pose_to_import_pose = false );

		void set_backbone_only( bool const setting );

		Size compute_size();

	private:

		void
		initialize_defaults();

		void cleanup_pose( pose::Pose & import_pose ) const;

		void check_sequence( pose::Pose const & pose, pose::Pose const & import_pose );

	private:

		import_pose::pose_stream::PoseInputStreamOP pose_input_stream_;
		utility::vector1< Size > input_res_;
		utility::vector1< Size > slice_res_;
		std::map< Size, Size > full_to_sub_;
		chemical::ResidueTypeSetCAP rsd_set_;
		bool backbone_only_;

		pose::PoseOP import_pose_;
	};


} //protein
} //modeler
} //stepwise
} //protocols

#endif
