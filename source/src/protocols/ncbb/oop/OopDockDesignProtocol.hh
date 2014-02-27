// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file   protocols/ncbb/oop/OopDockDesignProtocol.hh
///
/// @brief
/// @author Kevin Drew


#ifndef INCLUDED_protocols_ncbb_oop_OopDockDesignProtocl_hh
#define INCLUDED_protocols_ncbb_oop_OopDockDesignProtocl_hh

#include <protocols/ncbb/oop/OopDockDesignProtocol.fwd.hh>
#include <protocols/ncbb/oop/OopDockDesignProtocolCreator.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>

namespace protocols {
namespace ncbb {
namespace oop {

class OopDockDesignProtocol : public moves::Mover 
{

	public:

		//default ctor
		//OopDockDesignProtocol(): Mover("OopDockDesignProtocol"){}
		OopDockDesignProtocol();

		OopDockDesignProtocol(
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
			bool const oop_design_first,

			bool const pymol,
			bool const keep_history

		);

		OopDockDesignProtocol(
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
		virtual ~OopDockDesignProtocol(){}

		//methods
		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "OopDockDesignProtocol"; }
		void setup_filter_stats();
        protocols::moves::MoverOP fresh_instance() const { return OopDockDesignProtocolOP( new OopDockDesignProtocol ); }
        protocols::moves::MoverOP clone() const;
        void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	private:
		void setup_pert_foldtree( core::pose::Pose & pose);

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
		bool oop_design_first_;

		bool pymol_;
		bool keep_history_;

};


} // namespace oop
} // namespace ncbb
} // namespace protocols

#endif // INCLUDED_protocols_ncbb_oop_OopDockDesignProtocl_hh
