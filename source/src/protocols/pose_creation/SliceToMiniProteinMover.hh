// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/SliceToMiniProteinMover.hh
/// @brief This mover chops a big protein into miniProteins
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_protocols_pose_creation_SliceToMiniProteinMover_hh
#define INCLUDED_protocols_pose_creation_SliceToMiniProteinMover_hh

#include <basic/datacache/DataMap.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/pose_creation/SliceToMiniProteinMover.fwd.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <core/pack/task/TaskFactory.hh>

namespace protocols {
namespace pose_creation {

class SliceToMiniProteinMover : public moves::Mover {
public:
	struct SSElement{
		Size start_res;
		Size end_res;
		std::string type;
		SSElement(Size start_res_i, Size end_res_i, Size type_res_i){
			start_res = start_res_i;
			end_res = end_res_i;
			type = type_res_i;
		}
	};
	struct Chunk{
		core::pose::Pose pose;
		Size start_res;
		Size end_res;
		core::Real ddg;
		bool outputed;
		Chunk(core::pose::Pose pose_i, Size start_res_i, Size end_res_i){
			pose=pose_i;
			start_res=start_res_i;
			end_res=end_res_i;
			outputed=false;
		}
	};
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSelectorCOP ResidueSelectorCOP;

	/// @brief parse pdb into ssecondary structure elements based on input selector
	utility::vector1<SliceToMiniProteinMover::SSElement> parse_ss(core::pose::Pose const & pose) const;

	/// @brief parse pdb into initial chunks
	utility::vector1<SliceToMiniProteinMover::Chunk> parse_chunks(core::pose::Pose const & pose,utility::vector1<SliceToMiniProteinMover::SSElement> ss_elements);

	/// @brief filter chunks based on ddg
	utility::vector1<SliceToMiniProteinMover::Chunk> filter_chunks(utility::vector1<SliceToMiniProteinMover::Chunk> chunks);

	/// @brief trim any overhanging helices and helices that are slightly too long
	utility::vector1<SliceToMiniProteinMover::Chunk> trim_chunks(utility::vector1<SliceToMiniProteinMover::Chunk> chunks);

	/// @brief setups chunks for output. all, longest or best_ddg
	void setup_chunks_for_output(utility::vector1<SliceToMiniProteinMover::Chunk> & chunks);

	/// @brief  Constructor
	SliceToMiniProteinMover();

	core::pose::PoseOP get_additional_output() override;

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::scoring::ScoreFunctionOP scorefxn_; //dflt NULL
	protocols::moves::MoverOP relax_mover_; //dflt NULL; in the unbound state, prior to taking the energy, should we do any relaxation

	core::Size max_length_;
	core::Size min_sse_count_;
	core::Size min_sse_length_;
	core::Real ddg_ala_slice_score_;
	core::select::residue_selector::ResidueSelectorCOP selector_;
	utility::vector1<SliceToMiniProteinMover::Chunk> final_chunks_;
	std::string output_mode_;
};

} // pose_creation
} // protocols

#endif

