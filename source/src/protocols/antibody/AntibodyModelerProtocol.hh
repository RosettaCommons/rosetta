// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/AntibodyModelerProtocol.hh
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody_AntibodyModelerProtocol_hh
#define INCLUDED_protocols_antibody_AntibodyModelerProtocol_hh

#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyModelerProtocol.fwd.hh>
#include <protocols/constraint_movers/ConstraintSetMover.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>


namespace protocols {
namespace antibody {

class AntibodyModelerProtocol: public moves::Mover {
public:

	// default constructor
	AntibodyModelerProtocol();

	// default destructor
	~AntibodyModelerProtocol() override;

	protocols::moves::MoverOP clone() const override;

	/// @brief Assigns default values to primitive members
	void set_default();

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

	/// @brief Associates relevant options with the AntibodyModeler class
	static void register_options();

	// simple inline setters
	void set_BenchMark(bool setting) {
		benchmark_ = setting;
	}
	void set_ModelH3(bool setting) {
		model_h3_ = setting;
	}
	void set_SnugFit(bool setting) {
		snugfit_ = setting;
	}
	void set_refine_h3(bool setting) {
		refine_h3_ = setting;
	}
	void set_sc_min(bool setting) {
		sc_min_ = setting ;
	}
	void set_rt_min(bool setting) {
		rt_min_ = setting ;
	}
	void set_packonly_after_graft (bool setting) {
		packonly_after_graft_ = setting;
	}
	void set_cst_weight ( core::Real const cst_weight) {
		cst_weight_ = cst_weight;
	}
	void set_use_constraints(bool const use_csts) {
		use_csts_ = use_csts;
	}
	void set_constrain_cter(bool const setting) {
		constrain_cter_ = setting;
	}
	void set_auto_constraint(bool const setting) {
		auto_constraint_ = setting;
	}
	void set_all_atom_mode_kink_constraint(bool const setting) {
		all_atom_mode_constraint_ = setting;
	}
	void set_constrain_vlvh_qq(bool const setting) {
		constrain_vlvh_qq_ = setting;
	}
	void set_middle_pack_min( bool middle_pack_min) {
		middle_pack_min_ = middle_pack_min;
	}
	void set_extend_h3_before_modeling(bool setting) {
		extend_h3_before_modeling_ = setting;
	}

	void display_constraint_residues( core::pose::Pose & pose );

	void show( std::ostream & out=std::cout ) const override;
	friend std::ostream & operator<<(std::ostream& out, const AntibodyModelerProtocol & ab_m );


private:
	bool model_h3_;
	bool snugfit_;
	bool refine_h3_;
	bool LH_repulsive_ramp_;
	bool sc_min_;
	bool rt_min_;
	bool camelid_constraints_;
	bool middle_pack_min_;
	bool packonly_after_graft_;
	core::Real cen_cst_, high_cst_;
	bool extend_h3_before_modeling_;
	bool use_csts_;
	bool constrain_vlvh_qq_, constrain_cter_;
	bool auto_constraint_;
	bool all_atom_mode_constraint_;
	bool run_snugdock_ = false;

	// Benchmark mode for shorter_cycles
	bool benchmark_;

	bool user_defined_; // for constructor options passed to init

	bool flags_and_objects_are_in_sync_;

	// used as a flag to enable reading in of cst files
	core::Real cst_weight_;

	// score functions
	core::scoring::ScoreFunctionOP loop_scorefxn_highres_;
	core::scoring::ScoreFunctionOP loop_scorefxn_centroid_;
	core::scoring::ScoreFunctionOP dock_scorefxn_highres_;
	core::scoring::ScoreFunctionOP pack_scorefxn_;

	// constraint set mover
	protocols::constraint_movers::ConstraintSetMoverOP cdr_constraint_;

	// external objects
	AntibodyInfoOP ab_info_;


	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	void finalize_setup( core::pose::Pose & pose );

	/// @brief Sets up the instance of AntibodyModeler and initializes all members based on values passed in at construction
	///  or via the command line.
	void init();

	void setup_objects();

	/// @brief Output of various metrics of final model
	void echo_metrics_to_output( core::pose::Pose & pose );

	/// @brief Access AntibodyInfo to return loop for modeling, current 93-102 in Chothia numbering
	protocols::loops::Loop get_cdr_h3_loop();

}; // class


} // namespace antibody
} // namespace protocols

#endif

