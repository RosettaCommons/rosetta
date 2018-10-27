// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/indexed_structure_store/movers/SegmentSequenceProfileMover.hh

#pragma once

// C++ Headers
#include <string>
#include <map>

// Utility Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <core/sequence/SequenceProfile.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

#include <protocols/indexed_structure_store/StructureStore.fwd.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/indexed_structure_store/SegmentSequenceProfile.hh>
#include <protocols/indexed_structure_store/movers/SegmentSequenceProfileMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>


namespace protocols { namespace indexed_structure_store { namespace  movers {

class SegmentSequenceProfileMover : public protocols::moves::Mover {
public:
	static std::string
		pssm_output_order;

public:
	SegmentSequenceProfileMover() = default;

	void
	structure_store_path(std::string path);

	std::string
	structure_store_path();

	core::select::residue_selector::ResidueSelectorCOP selector() { return selector_; }
	void selector(core::select::residue_selector::ResidueSelectorCOP new_selector) { selector_ = new_selector; }

	SegmentSequenceProfileConfig profile_config() { return profile_config_; }
	void profile_config(SegmentSequenceProfileConfig new_profile_config) { profile_config_ = new_profile_config; }

	std::string output_lookup_count() { return output_lookup_count_; }
	void output_lookup_count(std::string new_output_lookup_count) { output_lookup_count_ = new_output_lookup_count; }

	std::string output_pssm_filename() { return output_pssm_filename_; }
	void output_pssm_filename(std::string new_output_pssm_filename) { output_pssm_filename_ = new_output_pssm_filename; }

	std::string output_pssm_inline() { return output_pssm_inline_; }
	void output_pssm_inline(core::SSize new_output_pssm_inline) { output_pssm_inline_ = new_output_pssm_inline; }

	bool apply_profile() { return apply_profile_; }
	void apply_profile(bool new_apply_profile) { apply_profile_ = new_apply_profile; }

	std::string scaling() { return scaling_; }
	void scaling(std::string new_scaling) { scaling_ = new_scaling; }

	core::Real weight() { return weight_; }
	void weight(core::Real new_weight) { weight_ = new_weight; }

	indexed_structure_store::StructureStoreOP
	structure_store();

	search::StructureDatabaseOP
	structure_database();

	moves::MoverOP
	clone() const override;

	void apply( Pose & pose ) override;

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

	core::sequence::SequenceProfile to_sequence_profile( SegmentSequenceProfileResult::ArrayXaa profile);

	void write_profile(core::sequence::SequenceProfile profile, std::ostream & outstream);

	std::string structure_store_path_;

	core::select::residue_selector::ResidueSelectorCOP selector_;

	SegmentSequenceProfileConfig profile_config_;

	std::string output_lookup_count_ = "";
	std::string output_pssm_filename_ = "";
	std::string output_pssm_inline_ = "";

	bool apply_profile_ = true;

	std::string scaling_ = "prob";
	core::Real weight_ = 1;

	StructureStoreOP structure_store_;
	search::StructureDatabaseOP structure_database_;
};


} } }
