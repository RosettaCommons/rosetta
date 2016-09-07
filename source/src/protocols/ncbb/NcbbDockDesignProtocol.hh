// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file   protocols/ncbb/NcbbDockDesignProtocol.hh
///
/// @brief
/// @author Andrew Watkins


#ifndef INCLUDED_protocols_ncbb_NcbbDockDesignProtocl_hh
#define INCLUDED_protocols_ncbb_NcbbDockDesignProtocl_hh

#include <protocols/ncbb/NcbbDockDesignProtocol.fwd.hh>
#include <protocols/ncbb/NcbbDockDesignProtocolCreator.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>

namespace protocols {
namespace ncbb {

class NcbbDockDesignProtocol : public moves::Mover
{

public:

	//default ctor
	//NcbbDockDesignProtocol(): Mover("NcbbDockDesignProtocol"){}
	NcbbDockDesignProtocol();

	NcbbDockDesignProtocol(
		core::scoring::ScoreFunctionOP score_function,
		core::Real const mc_temp,
		core::Real const pert_mc_temp,
		core::Real const pert_dock_rot_mag,
		core::Real const pert_dock_trans_mag,
		core::Real const pert_pep_small_temp,
		core::Real const pert_pep_small_H,
		core::Real const pert_pep_small_L,
		core::Real const pert_pep_small_E,
		core::Real const pert_pep_shear_temp,
		core::Real const pert_pep_shear_H,
		core::Real const pert_pep_shear_L,
		core::Real const pert_pep_shear_E,

		core::Size const pert_pep_num_rep,
		core::Size const pert_num,
		core::Size const dock_design_loop_num,

		bool const no_design,
		bool const final_design_min,
		bool const use_soft_rep,
		bool const mc_initial_pose,
		bool const ncbb_design_first,

		bool const pymol,
		bool const keep_history

	);

	NcbbDockDesignProtocol(
		core::scoring::ScoreFunctionOP score_function,
		core::Real const mc_temp,
		core::Real const pert_dock_rot_mag,
		core::Real const pert_dock_trans_mag,
		core::Size const dock_design_loop_num,
		bool const no_design,
		bool const final_design_min,
		bool const pymol,
		bool const keep_history
	);


	//default dtor
	~NcbbDockDesignProtocol() override= default;

	//methods
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override { return "NcbbDockDesignProtocol"; }
	protocols::moves::MoverOP fresh_instance() const override { return NcbbDockDesignProtocolOP( new NcbbDockDesignProtocol ); }
	protocols::moves::MoverOP clone() const override;
	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	core::Real mc_temp_;
	core::Real pert_mc_temp_;
	core::Real pert_dock_rot_mag_;
	core::Real pert_dock_trans_mag_;
	core::Real pert_pep_small_temp_;
	core::Real pert_pep_small_H_;
	core::Real pert_pep_small_L_;
	core::Real pert_pep_small_E_;
	core::Real pert_pep_shear_temp_;
	core::Real pert_pep_shear_H_;
	core::Real pert_pep_shear_L_;
	core::Real pert_pep_shear_E_;

	core::Size pert_pep_num_rep_;
	core::Size pert_num_;
	core::Size dock_design_loop_num_;

	bool no_design_;
	bool final_design_min_;
	bool use_soft_rep_;
	bool mc_initial_pose_;
	bool ncbb_design_first_;

	bool pymol_;
	bool keep_history_;

};


} // namespace ncbb
} // namespace protocols

#endif // INCLUDED_protocols_ncbb_ncbb_NcbbDockDesignProtocl_hh
