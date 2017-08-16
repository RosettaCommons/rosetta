// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
///
/// @file protocols/seeded_abinitio/SwapSegment.hh
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_SwapSegment_hh
#define INCLUDED_protocols_seeded_abinitio_SwapSegment_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace seeded_abinitio {

class SwapSegment : public protocols::moves::Mover {
public:
	typedef core::pose::Pose Pose;

	SwapSegment();

	void apply( core::pose::Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;


	~SwapSegment() override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	void copying_side_chains(
		core::pose::Pose & pose,
		core::pose::PoseOP & swap_segment,
		protocols::loops::Loops & seeds);

	void swap_segment(
		core::pose::Pose & pose,
		core::pose::PoseOP & swap_segment,
		protocols::loops::Loops & seeds);

	void swap_chain(
		core::pose::Pose & pose,
		core::pose::PoseOP & target_chain,
		core::Size chain_to_swap);

	bool copy_sidechains_;

	bool swap_segment_;

	core::Size swap_chain_;

	///input pdb that contains the segments that should be swapped
	core::pose::PoseOP seeds_pdb_;

	///check for the segments
	bool seeds_presence_;

	protocols::loops::Loops all_seeds_;

	///chain that contains the seed in the seed_pdb
	core::Size from_chain_;

	///chain in which the segments should be swapped/side chain replaced
	core::Size to_chain_;

	core::scoring::ScoreFunctionOP scorefxn_;

	///switch to determine what numbering needs to be used since
	///parse time is different from computing time and if the pose has changed, numbering will be off
	bool previously_grown_;
};

}
}

#endif
