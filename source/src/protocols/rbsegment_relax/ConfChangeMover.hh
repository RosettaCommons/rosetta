// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rbsegment_relax/ConfChangeMover.hh
/// @author Diego del Alamo (diego.delalamo@gmail.com) and Davide Sala (d.sala1388@gmail.com)

#ifndef INCLUDED_protocols_rbsegment_relax_ConfChangeMover_hh
#define INCLUDED_protocols_rbsegment_relax_ConfChangeMover_hh

#include <protocols/rbsegment_relax/ConfChangeMover.fwd.hh>
#include <protocols/rbsegment_relax/ConfChangeMoverCreator.hh>

#include <basic/datacache/DataMap.fwd.hh>

#include <boost/unordered/unordered_map.hpp>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.fwd.hh>

#include <core/select/residue_selector/TrueResidueSelector.fwd.hh>
#include <core/select/residue_selector/util.hh>

#include <core/types.hh>

#include <numeric/xyzTransform.fwd.hh>

#include <protocols/hybridization/CartesianHybridize.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rbsegment_relax/RBSegmentMover.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace rbsegment_relax {

class ConfChangeMover : public hybridization::CartesianHybridize {
public:

	ConfChangeMover();

	~ConfChangeMover() = default;

	protocols::moves::MoverOP clone() const override;

	void apply(core::pose::Pose &pose) override;

	void stage1_pose_setup(core::pose::Pose &pose);

	void stage1(core::pose::Pose &pose);

	bool check_gaps(core::pose::Pose &pose);

	rbsegment_relax::RBSegment multi_rb(core::pose::Pose &pose, bool const &multiple_sses);

	utility::vector1<utility::vector1<core::Real>> sse_contact_strength(core::pose::Pose &pose) const;

	std::set<core::Size> get_residues_in_rbsegment(rbsegment_relax::RBSegment const &rbseg) const;

	void recursive_residues_from_rbsegs(rbsegment_relax::RBSegment const &rbseg, std::set<core::Size> &residues) const;

	void check_or_create_fragments(core::pose::Pose &pose);

	void stage2(core::pose::Pose &pose, core::pose::Pose const &original_pose);

	void add_dihedral_csts_to_rb(core::pose::Pose &pose, core::pose::Pose const &original_pose);

	void add_dihedral_csts(core::pose::Pose &pose, core::pose::Pose const &original_pose);

	core::pose::PoseOP get_additional_output() override;

	static std::string mover_name();

	static void provide_xml_schema(utility::tag::XMLSchemaDefinition &xsd);

private:

	void parse_my_tag(utility::tag::TagCOP tag, basic::datacache::DataMap &data) override;

	std::string get_name() const override;

private:

	//////////////////
	// XML VARS

	utility::vector1<core::Size> newstart_, addstart_;
	utility::vector1<core::Size> newend_, addend_;
	utility::vector1<core::Size> seg_;
	utility::vector1<core::pose::PoseOP> store_poseOPs_;

	core::Size stage1_moves_ = 10000;
	core::Real stage1_multi_freq_ = 0.5;
	core::Real stage1_frag_freq_ = 0.5;
	core::Real stage1_twist_freq_ = 0.2;
	core::Real stage1_temp_ = 1.0;
	core::Size stage1_move_restore_ = 0;

	core::Size stage2_moves_ = 1000;
	core::Real stage2_segment_freq_ = 0.25;
	core::Real stage2_targgaps_freq_ = 0.25;
	core::Real stage2_temp_ = 1.0;
	core::Size stage2_models_ = 1;

	core::Real rot_ = 10.0;
	core::Real transl_ = 1.0;

	core::scoring::ScoreFunctionOP stage1_scorefxn_;
	core::scoring::ScoreFunctionOP stage2_scorefxn_;

	//////////////////
	// MOVER VARS
	rbsegment_relax::GaussianRBSegmentMover rbmover_;
	rbsegment_relax::HelicalGaussianMover helixmover_;
	core::fragment::FragSetOP frags3_;
	core::fragment::FragSetOP frags9_;

	std::map<core::Size, core::fragment::Frame> frames_;
	utility::vector1<rbsegment_relax::RBSegment> rbsegs_, rbsegs_2_,  stored_rebsegs_ ;
	loops::Loops loops_;

	core::pose::PoseOP template_;

	std::string rb_file_ = "AUTO";

	std::string const s1_rb_single_ = "1_RB1";
	std::string const s1_rb_multi_ = "1_RB2+";
	std::string const s1_twist_ = "1_TWIST";
	std::string const s1_frag_move_ = "1_FRAG3";

	std::string const s2_homolog_move_ = "2_TEMPLATE";
	std::string const s2_frame_select_ = "2_FRAG9T";
	std::string const s2_frame_random_ = "2_FRAG9R";

	utility::vector1<core::Size> res_adj_ = {3, 4, 5};

	bool stage1_minimization_ = false;

	core::select::residue_selector::ResidueSelectorCOP s1_res_ = nullptr;
	core::select::residue_selector::ResidueSelectorCOP s2_res_ = nullptr;
	core::select::residue_selector::ResidueSelectorCOP dihed2_res_ = nullptr;
	core::optimization::CartesianMinimizer minimizer_;
	utility::vector1<core::select::residue_selector::ResidueSelectorCOP> s1_rigid_res_;
};

}
}

#endif
