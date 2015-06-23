// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file docking_initialization_protocols
/// @brief initialization protocols for docking
/// @details
///		This contains the functions that create initial positions for docking
///		You can either randomize partner 1 or partner 2, spin partner 2, or
///		perform a simple perturbation.
/// @author Monica Berrondo


#ifndef INCLUDED_protocols_docking_DockingInitialPerturbation_hh
#define INCLUDED_protocols_docking_DockingInitialPerturbation_hh

#include <protocols/docking/types.hh>
#include <protocols/docking/DockingInitialPerturbation.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/docking/RigidBodyInfo.fwd.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <string>

//Numeric Headers
#include <numeric/xyzVector.hh>

namespace protocols {
namespace docking {

/// @brief this mover carries out the initial perturbation phase of the RosettaDock algorithm
/// based on user-inputted command line options
class DockingInitialPerturbation : public moves::Mover
{
public:
	/// @brief Default constructor
	DockingInitialPerturbation();

	/// @brief Constructor with two arguments. The first is the jump number to dock over, the second is a boolean (true
	///		will use slide into contact, false will not).
	DockingInitialPerturbation(
		core::Size const rb_jump,
		bool const slide=true
	);

	/// @brief Constructor with two arguments. The first is the DockJumps, the second is a boolean (true
	///		will use slide into contact, false will not).
	DockingInitialPerturbation(
		DockJumps const movable_jumps,
		bool const slide=true
	);

	//destructor
	~DockingInitialPerturbation();

	// protocol functions
	virtual void apply( core::pose::Pose & pose );
	void apply_body(core::pose::Pose & pose, core::Size jump_number );
	virtual std::string get_name() const;

	/// @brief Calls set_dault, register_from_options and init_from_options
	void init();

	/// @brief Sets members to default values
	void set_default();

	/// @brief Associates relevant options with the DockingInitialPerturbation class
	static void register_options();

	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

	/// @brief set functions
	void set_randomize1(bool randomize1){ randomize1_ = randomize1; }
	void set_randomize2(bool randomize2){ randomize2_ = randomize2; }
	void set_use_ellipsoidal_randomization(bool use_ellipsoidal_randomization){
		use_ellipsoidal_randomization_ = use_ellipsoidal_randomization;
	}
	void set_dock_pert(utility::vector1< core::Real > dock_pert){
		dock_pert_ = dock_pert;
		if_dock_pert_ = true;
	}
	void set_tilt(utility::vector1< core::Real > tilt){
		tilt_ = tilt;
	}
	void set_tilt1_center(std::string tilt_center){
		tilt1_center_ = tilt_center;
	}
	void set_tilt2_center(std::string tilt_center){
		tilt2_center_ = tilt_center;
	}
	void set_uniform_trans(core::Real uniform_trans){
		uniform_trans_ = uniform_trans;
		if_uniform_trans_ = true;
	}
	void set_spin( bool spin){ spin_ = spin;}
	void set_center( bool center) { center_at_interface_ = center;}

	/// zhe for rosetta_scripts
	void parse_my_tag(
										utility::tag::TagCOP tag,
										basic::datacache::DataMap &,
										protocols::filters::Filters_map const &,
										protocols::moves::Movers_map const &,
										core::pose::Pose const &
	);
	virtual protocols::moves::MoverOP clone() const; //zhe for scripts

private:
	/// do slide into context?
	bool slide_;

	// docking
	DockJumps movable_jumps_;

	bool randomize1_;
	bool randomize2_;
	bool use_ellipsoidal_randomization_;
	bool spin_;
	bool if_dock_pert_;
	bool if_uniform_trans_;
	bool multiple_jumps_;
	bool center_at_interface_;

	utility::vector1< core::Real > dock_pert_;
	core::Real uniform_trans_;

	utility::vector1< core::Real > tilt_;
	std::string tilt1_center_;
	std::string tilt2_center_;

	core::Vector slide_axis_;
	core::Vector spin_center_;

	protocols::docking::RigidBodyInfoOP rigid_body_info_;
};  // class DockingInitialPerturbation


/// @brief Contrary to the name, slides things apart first, then together.
/// OK for proteins, bad for ligands (because they may escape the pocket permanently).
class DockingSlideIntoContact : public moves::Mover
{
public:
	DockingSlideIntoContact();

	// constructor with arguments
	DockingSlideIntoContact( core::Size const rb_jump );

	DockingSlideIntoContact( core::Size const rb_jump, core::Vector const slide_axis );

	//destructor
	~DockingSlideIntoContact();

	// protocol functions
	virtual void apply( core::pose::Pose & pose );
	virtual void show(std::ostream & output=std::cout) const;

	// getters
	virtual std::string get_name() const;
	core::Size get_jump_num() const { return rb_jump_; }

private:
	core::scoring::ScoreFunctionOP scorefxn_;

	// which jump to use for docking
	core::Size rb_jump_;
	core::Vector slide_axis_; //used if a specific slide axis is specified in the constructor
	
};  // class DockingSlideIntoContact

std::ostream &operator<< ( std::ostream &os, DockingSlideIntoContact const &mover );


/// @brief Slides docking partners together by monitoring fa_rep.
/// @details
///		If partners are already touching, no change is made.
///		Separation will be 1A or less after calling this function.
class FaDockingSlideIntoContact : public moves::Mover
{
public:
	FaDockingSlideIntoContact();
	FaDockingSlideIntoContact( core::Size const rb_jump);
	FaDockingSlideIntoContact( utility::vector1<core::Size> rb_jumps);
	FaDockingSlideIntoContact( core::Size const rb_jump, core::Vector const slide_axis );

	//destructor
	~FaDockingSlideIntoContact();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual void show(std::ostream & output=std::cout) const;
	core::Size get_jump_num() const { return rb_jump_; }
	core::Real get_tolerance() const { return tolerance_; }

private:
	core::Size rb_jump_; // use this or rb_jumps_, not both
	utility::vector1<core::Size> rb_jumps_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real tolerance_; ///< how accurate do you want to be?
	core::Vector slide_axis_; //used if a specific slide axis is specified in the constructor

};  // class FaDockingSlideIntoContact

std::ostream &operator<< ( std::ostream &os, FaDockingSlideIntoContact const &fadock );

} // docking
} // protocols

#endif
