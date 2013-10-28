// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/symmetric_docking/SymDockProtocol.hh
///
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_symmetric_docking_SymDockProtocol_hh
#define INCLUDED_protocols_symmetric_docking_SymDockProtocol_hh

#include <protocols/moves/Mover.hh>
#include <protocols/symmetric_docking/SymDockProtocol.fwd.hh>
// AUTO-REMOVED #include <protocols/symmetric_docking/SymDockBaseProtocol.hh>

#include <core/types.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.fwd.hh>
#include <protocols/symmetric_docking/SymDockingHiRes.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/vector1.hh>

#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace symmetric_docking {

void SymDock_main();

class SymDockProtocol : public moves::Mover
{
public:

typedef core::conformation::symmetry::SymmetricConformation SymmetricConformation;
typedef core::conformation::symmetry::SymmetryInfo SymmetryInfo;

public:

  SymDockProtocol();

	SymDockProtocol(
		bool const fullatom,
		bool const local_refine,
		bool const view=false
	);

	SymDockProtocol(
		bool const fullatom,
		bool const local_refine,
		bool const view,
		core::scoring::ScoreFunctionOP docking_score_low,
		core::scoring::ScoreFunctionOP docking_score_high
	);

	virtual ~SymDockProtocol();

	/// @brief setup that is called from constructor
	void set_default();

	void register_options();

	/// @brief setter

	void set_dock_rtmin( bool dock_rtmin_in );

	void set_sc_min( bool sc_min_in );
	void set_max_repeats( Size const max_repeats_in );
	void set_dock_ppk( bool dock_ppk_in );

	void set_fullatom( bool const fullatom_in );

	void set_local_refine( bool const local_refine_in );

	void set_view( bool view_in );

	void set_lowres_scorefxn( core::scoring::ScoreFunctionOP docking_score_low_in );

	void set_highres_scorefxn( core::scoring::ScoreFunctionOP docking_score_high_in );

	void set_highres_scorefxn(
		core::scoring::ScoreFunctionOP docking_score_high_in,
		core::scoring::ScoreFunctionOP docking_score_pack_in );

	bool docking_lowres_filter( core::pose::Pose & pose );
	bool docking_highres_filter( core::pose::Pose & pose );

	core::Real
	calc_interaction_energy( core::pose::Pose & pose );

	core::Real
	calc_rms( core::pose::Pose & pose );

	/// @brief recovers the side-chains from the native-pose
	void recover_sidechains( core::pose::Pose & pose, const core::pose::Pose & native_pose );

	void task_factory( core::pack::task::TaskFactoryOP task_factory );

	// turn on design of partner2 during docking. Not thoroughly tested!
	void design( bool const des );
	bool design() const;

	// skip population of the score map
	void hurry( bool const hurry );

	core::pack::task::TaskFactoryOP task_factory() const;
	core::pack::task::TaskFactoryOP & task_factory();

	void score_only( core::pose::Pose & pose );

	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );


	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;

private:

 void
  classic_mcm_protocol(
    core::pose::Pose & pose,
    core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn,
    protocols::moves::MonteCarloOP monteCarlo,
    core::Size num_cycles,
    core::Size repack_every_Nth
  ) const;

  protocols::moves::MoverOP
  make_dockmcm_mover(
    core::pose::Pose const & pose,
    protocols::moves::MoverOP repack_mover,
    protocols::moves::MoverOP rigbod_mover,
    core::kinematics::MoveMapOP movemap, //< would be COP but MinMover wants OP
    core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn,
    protocols::moves::MonteCarloOP monteCarlo
  ) const;

private:

/// information about the mode
	bool fullatom_;
	bool local_refine_;
	bool rtmin_;
	bool sc_min_;
	Size max_repeats_;
	bool dock_ppk_;

	/// the jump number across which to do rigid_body transformations
	utility::vector1<core::Size> movable_jumps_;
	//core::Size rb_jump_;
	/// should docking change the foldtree?
	bool autofoldtree_;
	/// whether or not to initialize the viewer (for opengl)
	bool view_;
	bool design_;
	bool passed_lowres_filter_;
	bool passed_highres_filter_;
	bool hurry_; // skip populating the score map
	// for outputting scorefiles
	std::map < std::string, core::Real > score_map_;

	// for scoring
	core::scoring::ScoreFunctionOP docking_score_low_;
	core::scoring::ScoreFunctionOP docking_score_high_;
	core::scoring::ScoreFunctionOP docking_score_high_min_;
	core::scoring::ScoreFunctionOP docking_score_pack_;

	moves::MonteCarloOP mc_;

	//protocols
	protocols::symmetric_docking::SymDockingLowResOP docking_low_;
	protocols::symmetric_docking::SymDockingHiResOP docking_high_;
	core::pack::task::TaskFactoryOP init_task_factory_; // use this to restrict the packer task for docking protocol

};

} // symmetric_docking
} // protocols
#endif //INCLUDED_protocols_symmetric_docking_SymDockProtocol_HH
