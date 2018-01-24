// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/multistage_rosetta_scripts/TagManager.hh
/// @brief This class was designed to help the MRSJobQueen handle tags for many input jobs
/// @author Jack Maguire, jack@med.unc.edu


#ifndef INCLUDED_protocols_multistage_rosetta_scripts_TagManager_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_TagManager_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/multistage_rosetta_scripts/TagManager.fwd.hh>

#include <protocols/jd3/CompletedJobOutput.hh>
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobResult.fwd.hh>
#include <protocols/jd3/JobSummary.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>

#include <core/pose/Pose.hh>

#include <basic/datacache/DataMap.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>


namespace protocols {
namespace multistage_rosetta_scripts {

///@brief This class does not fail when you try to add an element that already exists. Instead, it just keeps the older element.
class NoFailDataMap : public basic::datacache::DataMap {
public:

	NoFailDataMap() :
		basic::datacache::DataMap()
	{}

	bool add(
		std::string const & type,
		std::string const & name,
		utility::pointer::ReferenceCountOP const op
	) override;

};

struct ParsedTagCache{
	ParsedTagCache(
		core::Size input_pose_id_in,
		NoFailDataMapOP data_map_in,
		TagMapOP mover_tags_in,
		TagMapOP filter_tags_in,
		utility::pointer::shared_ptr< moves::Movers_map > movers_map_in,
		utility::pointer::shared_ptr< filters::Filters_map > filters_map_in
	) :
		input_pose_id( input_pose_id_in ),
		data_map( data_map_in ),
		mover_tags( mover_tags_in ),
		filter_tags( filter_tags_in ),
		movers_map( movers_map_in ),
		filters_map( filters_map_in )
	{}

	core::Size input_pose_id;
	NoFailDataMapOP data_map;
	TagMapOP mover_tags;
	TagMapOP filter_tags;
	utility::pointer::shared_ptr< moves::Movers_map > movers_map;
	utility::pointer::shared_ptr< filters::Filters_map > filters_map;
};

class TagManager : public utility::pointer::ReferenceCount {
public:
	TagManager();
	TagManager( core::Size num_input_pose_ids );

	virtual ~TagManager();

	void set_num_input_pose_ids( core::Size num_input_pose_ids );

	inline void set_common_data_tags( std::list< utility::tag::TagCOP > && tags_in_order ){
		tag_list_for_input_pose_id_[ 0 ]->splice(
			tag_list_for_input_pose_id_[ 0 ]->begin(),
			tags_in_order
		);
	}

	inline void set_common_data_tags( std::list< utility::tag::TagCOP > const & tags_in_order ){
		std::list< utility::tag::TagCOP > temp_copy = tags_in_order;
		set_common_data_tags( std::move( temp_copy ) );
	}

	///@brief register all tags that are not present in common tags
	void register_data_tags_for_input_pose_id(
		core::Size input_pose_id,
		std::list< utility::tag::TagCOP > const & tags_in_order
	);

	///@brief overload of register_data_tags_for_input_pose_id that uses move semantics
	void register_data_tags_for_input_pose_id(
		core::Size input_pose_id,
		std::list< utility::tag::TagCOP > && tags_in_order
	);

	ParsedTagCacheOP generate_data_for_input_pose_id(
		core::Size input_pose_id,
		core::pose::Pose const &
	);

	inline ParsedTagCacheOP generate_data_for_input_pose_id( core::Size input_pose_id ){
		return generate_data_for_input_pose_id( input_pose_id, dummy_pose_ );
	}

protected:
	std::vector< TagListOP > const & tag_list_for_input_pose_id() const {
		return tag_list_for_input_pose_id_;
	}

	std::vector< TagListOP > & tag_list_for_input_pose_id() {
		return tag_list_for_input_pose_id_;
	}

private:
	std::vector< TagListOP > tag_list_for_input_pose_id_;//[0] only holds common data

	///@brief cache the most recent request because it may get called again soon
	ParsedTagCacheOP most_recent_request_;

	core::pose::Pose dummy_pose_;
};

} //multistage_rosetta_scripts
} //protocols

#endif
