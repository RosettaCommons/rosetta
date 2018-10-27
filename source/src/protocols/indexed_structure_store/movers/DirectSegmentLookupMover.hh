// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/indexed_structure_store/movers/DirectSegmentLookupMover.hh

#pragma once

// C++ Headers
#include <string>
#include <map>

// Utility Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <protocols/indexed_structure_store/StructureStore.fwd.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/indexed_structure_store/DirectSegmentLookup.hh>
#include <protocols/indexed_structure_store/movers/DirectSegmentLookupMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>


namespace protocols { namespace indexed_structure_store { namespace  movers {

class DirectSegmentLookupMover : public protocols::moves::Mover {
public:
	DirectSegmentLookupMover() = default;

	void
	structure_store_path(std::string path);

	std::string
	structure_store_path();

	void lookup_config(DirectSegmentLookupConfig config);

	DirectSegmentLookupConfig lookup_config();

	std::string stored_subset_name() { return stored_subset_name_; }
	void stored_subset_name(std::string name) {  stored_subset_name_ = name; }

	bool overwrite_stored_subset() { return overwrite_stored_subset_; }
	void overwrite_stored_subset(bool overwrite) {  overwrite_stored_subset_ = overwrite; }

	std::string label_insertion() { return label_insertion_; }
	void label_insertion(std::string name) {  label_insertion_ = name; }

	std::string output_lookup_count() { return output_lookup_count_; }
	void output_lookup_count(std::string name) {  output_lookup_count_ = name; }

	std::string output_lookup_length() { return output_lookup_length_; }
	void output_lookup_length(std::string name) {  output_lookup_length_ = name; }

	Size from_chain() { return from_chain_; }
	void from_chain(Size chain) { from_chain_ = chain; }

	Size to_chain() { return to_chain_; }
	void to_chain(Size chain) { to_chain_ = chain; }

	StructureStoreOP
	structure_store();

	search::StructureDatabaseOP
	structure_database();

	moves::MoverOP
	clone() const override;

	void apply( Pose & pose ) override;

	core::pose::PoseOP
	get_additional_output() override;

	std::string
	get_name() const override;

	static std::string
	mover_name();

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

private:
	void init_structure_store();

	core::pose::PoseOP
	generate_result_pose(core::pose::Pose & source_pose, DirectSegmentLookupResult lookup_result);

	std::string structure_store_path_;

	DirectSegmentLookupConfig lookup_config_;

	core::Size lookup_context_ = 2;

	std::string stored_subset_name_;
	bool overwrite_stored_subset_ = false;

	std::string label_insertion_;
	std::string output_lookup_count_;
	std::string output_lookup_length_;

	Size from_chain_ = 1;
	Size to_chain_ = 2;

	StructureStoreOP structure_store_;
	search::StructureDatabaseOP structure_database_;

	core::pose::PoseOP source_pose_;
	std::vector<DirectSegmentLookupResult> lookup_results_;
	Size lookup_result_index_ = 0;
	Size max_num_results_ = 0;
};


} } }
