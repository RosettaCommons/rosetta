// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   TemperedDocking.hh
///
/// @brief allows low-resolution docking using simulated or parallel tempering
/// @author Oliver Lange

#ifndef INCLUDED_protocols_docking_TemperedDocking_hh
#define INCLUDED_protocols_docking_TemperedDocking_hh

// Unit Headers
#include <protocols/docking/TemperedDocking.fwd.hh>

// Package Headers
#include <protocols/docking/types.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/canonical_sampling/TemperingBase.fwd.hh>

#include <utility/tag/Tag.fwd.hh>


// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers

// ObjexxFCL Headers

// C++ headers

namespace protocols {
namespace docking {

/// @brief allows docking using simulated or parallel tempering
/// @detailed
class TemperedDocking : public moves::Mover
{
public:
	/// @brief Associates relevant options with the TemperedDocking class
	static void register_options();

	/// @brief default constructor fills values with the expected defaults
	TemperedDocking();

	/// @brief clone
	virtual protocols::moves::MoverOP clone() const;

	///@brief copy ctor
	TemperedDocking( TemperedDocking const & rhs );

	///@brief assignment operator
	TemperedDocking & operator=( TemperedDocking const & rhs );

	/// @brief Assigns default values to primitive members
	void set_defaults();

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	/// @brief Sets the score function that will be used in the low-resolution phase
	void set_lowres_scorefxn( core::scoring::ScoreFunctionOP docking_scorefxn_low );

	/// @brief Sets the score function that will be used in the high-resolution phase.
	///		The same score function will be used for evaluating moves, packing and discriminating
	void set_highres_scorefxn( core::scoring::ScoreFunctionOP docking_scorefxn_high );

	/// @brief Sets the score function that will be used in the high-resolution phase.
	///		The first scorefunction will be used for evaluating moves and discriminating, the second will be used for packing
	void set_highres_scorefxn(
		core::scoring::ScoreFunctionCOP docking_scorefxn_high,
		core::scoring::ScoreFunctionCOP docking_scorefxn_pack );

	/// @brief Sets the score function that will be used in the high-resolution phase.
	///		The first scorefunction will be used for evaluating moves, the second will be used for packing and the third for discriminating
	void set_highres_scorefxn(
		core::scoring::ScoreFunctionCOP docking_scorefxn_high,
		core::scoring::ScoreFunctionCOP docking_scorefxn_pack,
		core::scoring::ScoreFunctionCOP docking_scorefxn_output);

	void set_sc_min( bool sc_min );

	virtual void apply( core::pose::Pose & pose );


	//getters for const access to movers and data of docking protocol
	protocols::moves::MoverCOP to_centroid() const;

	std::string partners() const { return partners_;} /// @brief returns the docking partners chain identifiers

	virtual std::string get_name() const { return "TemperedDocking"; }

	DockJumps & movable_jumps(){ return movable_jumps_;} ///@brief returns ref to the jumps vector for docking

	DockJumps const & movable_jumps() const { return movable_jumps_; } ///@ return const ref to the jumps vector for docking


	//setters
	void set_autofoldtree( bool setting ){ autofoldtree_ = setting; }

	void set_partners( std::string const& setting ){ partners_=setting; }

	void set_movable_jumps( DockJumps const& setting ){ movable_jumps_ = setting; }

	void set_use_constraints( bool setting );


	void add_jump( core::SSize const jump_number ){ movable_jumps_.push_back( int( jump_number ) ); }
	void show( std::ostream & out=std::cout ) const;
	friend std::ostream & operator<<(std::ostream& out, const TemperedDocking & dp );

	// function for the parser with lots of accessors
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

protected:
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	void finalize_setup( core::pose::Pose & pose );

	/// @brief Sets up the instance of TemperedDocking and initializes all members based on values passed in at construction
	///		or via the command line.
  void init(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionCOP docking_score_low,
		core::scoring::ScoreFunctionCOP docking_score_high
	);

	void copy( TemperedDocking & lhs, TemperedDocking const & rhs);

	void setup_objects();

private:
	/// information about the mode
	bool autofoldtree_;

	bool flags_and_objects_are_in_sync_;
	bool first_apply_with_current_setup_;

	bool sc_min_;

	/// how many cycles of monte-carlo per trajectory in the sampler ?
	core::Size n_cycles_;

	core::kinematics::FoldTree fold_tree_;
	std::string partners_;

	std::string previous_sequence_;

	/// jumps that rigid_body transformations can occur over
	DockJumps movable_jumps_;

	// score functions
	core::scoring::ScoreFunctionCOP docking_scorefxn_low_;
	core::scoring::ScoreFunctionCOP docking_scorefxn_high_;
	core::scoring::ScoreFunctionCOP docking_scorefxn_pack_;
	core::scoring::ScoreFunctionCOP docking_scorefxn_output_;

	// constraint set mover
	protocols::moves::MoverOP docking_constraint_;

	protocols::rigid::RigidBodyPerturbNoCenterMoverOP rb_mover_;
	protocols::canonical_sampling::TemperingBaseOP tempering_;
	protocols::canonical_sampling::MetropolisHastingsMoverOP sampler_;
	//if side-chains are to be taken from specified pdb file... it is set here...
	std::string recover_sidechains_filename_;

	protocols::moves::MoverOP to_centroid_;

	bool use_csts_;

	core::Real rigid_rot_mag_;
	core::Real rigid_trans_mag_;

	static bool options_registered_;

};
} // docking
} // protocols

#endif

