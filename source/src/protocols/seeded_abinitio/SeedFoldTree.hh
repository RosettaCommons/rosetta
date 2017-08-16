// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/seeded_abinitio/SeedFoldTree.cc
/// @brief
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_SeedFoldTree_hh
#define INCLUDED_protocols_seeded_abinitio_SeedFoldTree_hh

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/string_util.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <set>

namespace protocols {
namespace seeded_abinitio {

class SeedFoldTree : public protocols::moves::Mover
{
public:
	SeedFoldTree();
	SeedFoldTree( core::kinematics::FoldTreeOP ft );
	~SeedFoldTree() override;
	void fold_tree(core::kinematics::FoldTreeOP ft );
	core::kinematics::FoldTreeOP fold_tree() const;
	void apply( core::pose::Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	utility::vector1 < core::Size > get_cutpoints();
	core::Size best_by_ala_scan( core::Size start, core::Size end, core::pose::PoseOP & ts_pose );

	std::set< core::Size > get_folding_vertices();

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;


	core::kinematics::FoldTreeOP set_foldtree(
		core::pose::PoseOP & seed_target_chain,
		std::string secstr,
		protocols::loops::Loops & loops,
		bool protein_not_folded_yet );

	bool ddg_based();
	void ddg_based( bool ddgb );
	void set_anchor_res( utility::vector1< core::Size > anchor );
	void anchor_specified( bool anchor_specified );
	bool anchor_specified();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// newly generated foldtree OP
	core::kinematics::FoldTreeOP fold_tree_;

	/// is there a second target chain to add
	bool add_target_; //= false; not implemented

	///check whether there is a second target chain loaded
	bool pdb_contains_target_;// = false;

	///is the input pose submitted with the target
	bool twochains_;

	///vector containing all cutpoints
	utility::vector1<Size> cut_points_;

	///seed info
	protocols::loops::Loops all_seeds_;

	/// for manual set up
	bool set_jumps_manually;

	/// anchor residues
	utility::vector1< core::Size > anchors_;
	bool anchor_specified_;

	/// should the jump atoms be computed based on ddG
	bool ddg_based_;

	/// scorefunction for ala scan
	core::scoring::ScoreFunctionOP scorefxn_;

	//utility::vector1<core::Size> manual_jumps;
	utility::vector1< std::pair< Size, Size > > manual_jump_pairs_;
	core::pose::PoseOP template_pdb_;
	core::pose::PoseOP target_chain_;
	core::pose::PoseOP seeds_only_;
	core::pose::PoseOP only_seeds_chain_;
	std::set< core::Size > folding_vertices_;
};
}//end seeded_abinitio
}//end protocols

#endif
