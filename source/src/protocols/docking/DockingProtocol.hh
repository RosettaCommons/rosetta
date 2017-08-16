// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   DockingProtocol.hh
///
/// @brief
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_docking_DockingProtocol_hh
#define INCLUDED_protocols_docking_DockingProtocol_hh

// Unit Headers
#include <protocols/docking/DockingProtocol.fwd.hh>

// Package Headers
#include <protocols/docking/types.hh>
#include <protocols/docking/DockingEnsemble.fwd.hh>
#include <protocols/docking/DockingLowResEnsemble.fwd.hh>
#include <protocols/docking/DockFilters.fwd.hh>
#include <protocols/docking/DockingLowRes.fwd.hh>
#include <protocols/docking/DockingHighRes.fwd.hh>
#include <protocols/docking/DockingInitialPerturbation.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <protocols/toolbox/task_operations/InterfaceTaskOperation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>  // Needs to be the full header so the scorefxn can default to NULL
#include <core/kinematics/FoldTree.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/simple_moves/ConstraintSetMover.fwd.hh>
#include <protocols/simple_moves/ReturnSidechainMover.fwd.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers

// ObjexxFCL Headers

// C++ headers

namespace protocols {
namespace docking {

/// @brief This is the standard RosettaDock protocol
/// @details RosettaDock protocol based on [refs...Gray2003, Wang2005, Chaudhury 2007 ... ]
class DockingProtocol : public moves::Mover
{
public:
	/// @brief
	///  empty constructor fills values with the expected defaults
	///  rb_jump will be assigned as 1 (meaning the first jump will
	///  be used as the jump across which all rigid-body perturbations
	///  will occur
	DockingProtocol();

	// the logic for "low_res_protocol_only" and "docking_local_refine" is
	// counterintuitive: what it means:
	// "low-res-only" "local-refine" => lowres / highres will be run
	//  false          false            true     true   (counterintuitive)
	//  true           false            true     false  (good)
	//  false          true             false    true   (good)
	//  true           true             false    false  (counterintuitive)

	DockingProtocol(
		Size const rb_jump_in,
		bool const low_res_protocol_only = false, // if true: skip high resolution docking
		bool const docking_local_refine = false, // if true: skip low resolution docking
		bool const autofoldtree = true,
		core::scoring::ScoreFunctionOP docking_score_low = core::scoring::ScoreFunctionOP( nullptr ),
		core::scoring::ScoreFunctionOP docking_score_high = core::scoring::ScoreFunctionOP( nullptr )
	);

	DockingProtocol(
		DockJumps const movable_jumps,
		bool const low_res_protocol_only = false, // if true: skip high resolution docking
		bool const docking_local_refine = false, // if true: skip low resolution docking
		bool const autofoldtree = true,
		core::scoring::ScoreFunctionOP docking_score_low = core::scoring::ScoreFunctionOP( nullptr ),
		core::scoring::ScoreFunctionOP docking_score_high = core::scoring::ScoreFunctionOP( nullptr )
	);

	/// @brief Assigns default values to primitive members
	void set_default();

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	//destructor
	~DockingProtocol() override;

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	bool reinitialize_for_new_input() const override;// {return true;};

	/// @brief copy ctor
	DockingProtocol( DockingProtocol const & rhs );

	/// @brief assignment operator
	DockingProtocol & operator=( DockingProtocol const & rhs );

	/// @brief Associates relevant options with the DockingProtocol class
	static void register_options();

	/// @brief Sets the score function that will be used in the low-resolution phase
	void set_lowres_scorefxn( core::scoring::ScoreFunctionOP docking_scorefxn_low );

	/// @brief Sets the score function that will be used in the high-resolution phase.
	/// The same score function will be used for evaluating moves, packing and discriminating
	void set_highres_scorefxn( core::scoring::ScoreFunctionOP docking_scorefxn_high );

	/// @brief Sets the score function that will be used in the high-resolution phase.
	/// The first scorefunction will be used for evaluating moves and discriminating, the second will be used for packing
	void set_highres_scorefxn(
		core::scoring::ScoreFunctionOP docking_scorefxn_high,
		core::scoring::ScoreFunctionOP docking_scorefxn_pack
	);

	/// @brief Sets the score function that will be used in the high-resolution phase.
	/// The first scorefunction will be used for evaluating moves, the second will be used for packing and the third for
	/// discriminating
	void set_highres_scorefxn(
		core::scoring::ScoreFunctionOP docking_scorefxn_high,
		core::scoring::ScoreFunctionOP docking_scorefxn_pack,
		core::scoring::ScoreFunctionOP docking_scorefxn_output
	);

	void set_sc_min( bool sc_min );
	void set_rt_min( bool rt_min );
	void set_dock_min( bool const dock_min );

	void set_no_filters( bool no_filters );
	void set_low_res_protocol_only( bool const low_res_protocol_only );
	void set_docking_local_refine( bool const docking_local_refine );
	void set_use_legacy_protocol( bool const use_legacy_protocol );
	void set_cst_weight( core::Real const cst_weight );
	void set_use_constraints( bool const use_csts );
	void set_interface_definition_task_operation( protocols::toolbox::task_operations::InterfaceTaskOperationOP interface_definition );

	void set_additional_task_operarations( utility::vector1< core::pack::task::operation::TaskOperationOP > additional_task_operations );
	void add_additional_task_operaration( core::pack::task::operation::TaskOperationOP task_operation );
	utility::vector1< core::pack::task::operation::TaskOperationOP > get_additional_task_operarations();

	void apply( core::pose::Pose & pose ) override;

	// score_only is no longer implemented.  It remains here until a decision is made about what to do with it.
	// void score_only( core::pose::Pose & pose );

	// inline getters
	/// @brief returns the docking partners chain identifiers
	std::string partners() const { return partners_;}


	/// @brief returns ref to the jumps vector for docking
	DockJumps & movable_jumps(){ return movable_jumps_;}

	/// @brief returns a const ref to the jumps vector for docking
	DockJumps const & movable_jumps() const { return movable_jumps_; }
	core::pack::task::TaskFactory const & task_factory() { return *init_task_factory_; }

	//getters for const access to movers and data of docking protocol
	protocols::simple_moves::SwitchResidueTypeSetMoverCOP to_centroid() const;
	protocols::moves::MoverCOP to_all_atom() const;
	protocols::docking::DockingLowResCOP docking_lowres_mover() const;
	protocols::docking::DockingHighResCOP docking_highres_mover() const;
	protocols::docking::DockingInitialPerturbationCOP perturber() const;

	// Allow a developer to set a custom high resolution mover
	void set_docking_highres_mover( protocols::docking::DockingHighResOP docking_highres_mover );

	// Allow a developer to add some additional steps to the end of the low resolution phase.
	// Each mover that's added in this way is added to a SequenceMover that is applied after low resolution
	// docking concludes.
	void add_additional_low_resolution_step( protocols::moves::MoverOP additional_low_resolution_mover );

	// inline setters
	void set_autofoldtree( bool const autofoldtree ){ autofoldtree_ = autofoldtree; }
	void set_partners( std::string const& partners ){ partners_=partners; }
	void set_inner_cycles( core::Size inner_cycles ) { lowres_inner_cycles_=inner_cycles; }
	void set_outer_cycles( core::Size outer_cycles ) { lowres_outer_cycles_=outer_cycles; }

	// For RosettaScripts
	// TODO: Remove when legacy high res disappears
	void set_design( bool const design ) { design_ = design; }
	void set_task_factory( core::pack::task::TaskFactoryOP task_factory ){ init_task_factory_ = task_factory; }
	void set_ignore_default_docking_task(bool const ignore_default_docking_task){ignore_default_docking_task_ = ignore_default_docking_task;}
	void set_movable_jumps( DockJumps const jump_numbers ){ movable_jumps_ = jump_numbers; }
	void set_reporting( bool report ) { reporting_ = report; }
	void set_ensemble1( std::string  const& ensemble1 ) { ensemble1_filename_ = ensemble1; }
	void set_ensemble2( std::string  const& ensemble2 ) { ensemble2_filename_ = ensemble2; }
	void set_recover_sidechains_filename( std::string const& file ) { recover_sidechains_filename_ = file; }

	// Other member functions
	void add_jump( core::Size const jump_number ){ movable_jumps_.push_back( int( jump_number ) ); }

	void show( std::ostream & out=std::cout ) const override;

	// function for the parser with lots of accessors
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// information about the mode
	bool user_defined_; // for constructor options passed to init
	bool low_res_protocol_only_;
	bool reporting_;
	bool autofoldtree_;

	bool flags_and_objects_are_in_sync_;
	bool first_apply_with_current_setup_;

	bool sc_min_;
	bool rt_min_;
	bool dock_min_;

	bool if_ensemble_;  // Add by DK

	bool no_filters_;
	bool use_legacy_protocol_;
	bool docking_local_refine_;

	bool use_csts_;
	core::Real cst_weight_;
	core::Real cst_fa_weight_;

	core::Real score_cutoff_;
	core::kinematics::FoldTree fold_tree_;
	std::string partners_;

	std::string previous_sequence_;

	// low res options
	core::SSize lowres_inner_cycles_, lowres_outer_cycles_;

	/// jumps that rigid_body transformations can occur over
	DockJumps movable_jumps_;

	// score functions
	core::scoring::ScoreFunctionOP docking_scorefxn_low_;
	core::scoring::ScoreFunctionOP docking_scorefxn_high_;
	core::scoring::ScoreFunctionOP docking_scorefxn_pack_;
	core::scoring::ScoreFunctionOP docking_scorefxn_output_;

	// success criteria enforcers
	protocols::moves::MonteCarloOP mc_; //not used currently
	protocols::docking::DockingLowResFilterOP lowres_filter_;
	protocols::docking::DockingHighResFilterOP highres_filter_;


	//protocols
	protocols::docking::DockingLowResOP docking_lowres_mover_;
	protocols::docking::DockingHighResOP docking_highres_mover_;
	protocols::docking::DockingInitialPerturbationOP perturber_;

	// atom set switch movers
	protocols::simple_moves::SwitchResidueTypeSetMoverOP to_centroid_;
	protocols::moves::MoverOP to_all_atom_;

	// ensemble objects
	protocols::docking::DockingEnsembleOP ensemble1_;
	protocols::docking::DockingEnsembleOP ensemble2_;
	std::string ensemble1_filename_, ensemble2_filename_;

	// constraint set mover
	protocols::simple_moves::ConstraintSetMoverOP docking_constraint_;

	protocols::simple_moves::ReturnSidechainMoverOP recover_sidechains_;

	protocols::moves::SequenceMoverOP additional_low_resolution_steps_;

	//if side-chains are to be taken from specified pdb file... it is set here...
	std::string recover_sidechains_filename_;

	core::pack::task::TaskFactoryOP init_task_factory_; // use this to restrict the packer task for docking protocol
	bool design_; // for RosettaScripts.  to be deprecated when legacy high res disappears
	bool ignore_default_docking_task_; //passed down to DockingHighRes, prevents the default DockingTaskFactory from being built

	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	void finalize_setup( core::pose::Pose & pose );

	/// @brief Sets up the instance of DockingProtocol and initializes all members based on values passed in at construction
	///  or via the command line.
	void init(
		DockJumps const movable_jumps,
		bool const low_res_protocol_only_,
		bool const docking_local_refine,
		bool const autofoldtree,
		core::scoring::ScoreFunctionOP docking_score_low,
		core::scoring::ScoreFunctionOP docking_score_high
	);

	void initForEqualOperatorAndCopyConstructor(DockingProtocol & lhs, DockingProtocol const & rhs);

	void setup_objects();

	void setup_constraints( core::pose::Pose & pose );
	void add_constraints_to_scorefunction();
	void check_high_res_protocol();

};

std::ostream & operator<<(std::ostream & out, DockingProtocol const & dp );

} // docking
} // protocols

#endif

