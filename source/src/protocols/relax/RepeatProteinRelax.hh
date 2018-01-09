// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RepeatProteinRelax.cc
/// @brief This code is to be used instead of relax/minimize for repeat proteins. It maintains perfect symmetry while running relax
/// @author TJ Brunette

// Unit headers

#ifndef INCLUDED_protocols_simple_moves_RepeatProteinRelax_hh
#define INCLUDED_protocols_simple_moves_RepeatProteinRelax_hh

#include <protocols/relax/RepeatProteinRelax.fwd.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// C++ Headers
#include <string>
#include <utility/vector1.hh>

//unit header
#include <basic/Tracer.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/moves/Mover.hh>


namespace protocols {
namespace relax {

class RepeatProteinRelax : public RelaxProtocolBase {
public:
	RepeatProteinRelax();
	Size get_midpoint_longest_loop(core::pose::Pose const & pose,Size const repeat_length);
	void setup_repeat_symminfo(Size const repeatlen, Size const nrepeat, Size const base_repeat, bool const with_jumping, core::conformation::symmetry::SymmetryInfo & symminfo);
	void setup_repeat_pose(core::pose::Pose & pose);
	void setup_repeat_pose_jumping(core::pose::Pose & pose);
	void minimize_pose(core::pose::Pose & pose);
	void relax_pose(core::pose::Pose & pose);
	void setup_movemap(core::pose::Pose & pose);
	void seal_jumps(core::pose::Pose & pose);
	void add_residue_labels_back(core::pose::Pose & pose, std::map<core::Size, utility::vector1<std::string> > res_label_map,int symmetry_resid_offset);
	void apply(core::pose::Pose & pose ) override;
	moves::MoverOP clone() const override { return moves::MoverOP( new RepeatProteinRelax( *this ) ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &) override;
	std::string get_name() const override;
	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
private:
	Size numb_repeats_;
	bool minimize_;
	bool loop_cutpoint_mode_;
	Size relax_iterations_;
	bool cartesian_;
	bool ramp_down_constraints_;
	bool remove_symm_;
	bool modify_symmetry_and_exit_;
	std::string min_type_;
	core::scoring::ScoreFunctionOP sfxn_;

};
} //protocols
} //moves
#endif
