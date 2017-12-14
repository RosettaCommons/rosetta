// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Monica Berrondo
/// @author Ian Davis
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_rigid_RigidBodyMover_hh
#define INCLUDED_protocols_rigid_RigidBodyMover_hh

// Unit headers
#include <protocols/rigid/RigidBodyMover.fwd.hh>

// Package headers
#include <protocols/canonical_sampling/ThermodynamicMover.hh>


// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/conformation/symmetry/SymDof.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>

namespace protocols {
namespace rigid {

/// @brief Partner, which partner gets moved
enum Partner {
	partner_upstream = 1,
	partner_downstream //=2
} ;

/// @brief Direction, which direction
enum Direction {
	c2n=-1,
	random=0,
	n2c=1
} ;

///////////////////////////////////////////////////////////////////////////////
/// @brief Rigid-body random translate/rotate around centroid of downstream side of a jump.
/// @details We operate on a single jump rather than e.g. a randomly selected
/// jump from the mobile jumps in a MoveMap.
/// If you want a random choice among jumps, put multiple RigidBodyMovers
/// into a RandomMover (which will give you more control, anyway).
class RigidBodyMover : public protocols::canonical_sampling::ThermodynamicMover {
public:
	typedef ThermodynamicMover parent;

public:

	// default constructor
	RigidBodyMover();

	// constructor with arguments
	RigidBodyMover(
		int const rb_jump_in,
		Direction dir_in=n2c
	);

	RigidBodyMover( RigidBodyMover const & src );

	~RigidBodyMover() override;

	/// @brief Manual override of rotation center.
	void rot_center( core::Vector const & rot_center_in ) { rot_center_ = rot_center_in; }

	void apply( core::pose::Pose & pose ) override = 0;
	std::string get_name() const override;
	void show(std::ostream & output=std::cout) const override;


	bool
	preserve_detailed_balance() const override { return true; }

	/// @brief set whether detailed balance is preserved (i.e. no branch angle optimization during moves)

	void
	set_preserve_detailed_balance(
		bool
	) override {};

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges

	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges( core::pose::Pose & pose ) override;

	// data
protected:
	int rb_jump_;
	/// direction of folding (n-term to c-term)
	Direction dir_;

	/// center of rotation
	core::Vector rot_center_;
	bool freeze_; // use the same movement as before (if one is set) during apply

public:
	// Can't change jump number after creation b/c it determines center of rotation!
	//void rb_jump_( int const rb_jump_in ) { rb_jump_ = rb_jump_in; }
	int rb_jump() const { return rb_jump_; }
	void rb_jump(int jump_id){rb_jump_ = jump_id;} // set jump

	void unfreeze(){freeze_=false;}  // create a new trans_axis during apply
	void freeze(){freeze_=true;} // use the same trans_axis as before (if one is set) during apply
};  // class RigidBodyMover


/// @brief This Mover does a perturbation defined by the rotational and translational magnitudes.
class RigidBodyPerturbMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:

	// default constructor
	RigidBodyPerturbMover();

	// constructor with arguments (rb_jump defined)
	RigidBodyPerturbMover(
		int const rb_jump_in,
		core::Real const rot_mag_in,
		core::Real const trans_mag_in,
		Partner const partner_in=partner_downstream,
		bool interface_in=false //rot_center calculated at interface
	);

	// constructor with arguments (rb_jump defined)
	RigidBodyPerturbMover(
		int const rb_jump_in,
		core::Real const rot_mag_in,
		core::Real const trans_mag_in,
		Partner const partner_in,
		utility::vector1< bool > const & ok_for_centroid_calculation
	);

	// constructor with arguments (movable jumps defined by a movemap)
	RigidBodyPerturbMover(
		core::pose::Pose const & pose_in,
		core::kinematics::MoveMap const & mm,
		core::Real const rot_mag_in,
		core::Real const trans_mag_in,
		Partner const partner_in=partner_downstream,
		bool interface_in=false //rot_center calculated at interface
	);

	// constructor with arguments (overloaded for default jump, but defined rot/mag)
	RigidBodyPerturbMover(
		core::Real const rot_mag_in,
		core::Real const trans_mag_in,
		Partner const partner_in=partner_downstream,
		bool interface_in=false //rot_center calculated at interface
	);

	RigidBodyPerturbMover( RigidBodyPerturbMover const & );
	~RigidBodyPerturbMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	void show(std::ostream & output=std::cout) const override;

	void rot_magnitude( core::Real const magnitude ) { rot_mag_ = magnitude; }

	void trans_magnitude( core::Real const magnitude ) { trans_mag_ = magnitude; }
	core::Distance get_trans_mag() const;
	core::Angle get_rot_mag() const;

	/// @brief Manual override of rotation center.
	void rot_center( core::Vector const & /*rot_center_in*/ ); // recreate unless freeze is specified.

protected:
	// perturbation magnitudes (rotational and translational)
	core::Real rot_mag_;
	core::Real trans_mag_;
	utility::vector1<core::Real>  rb_delta_; // size 6: translateXYZ, rotateXYZ

private:
	Partner partner_;
	bool interface_;
	utility::vector1<core::Size> movable_jumps_;
	utility::vector1< bool > ok_for_centroid_calculation_;
};  // class RigidBodyPerturbMover

std::ostream &operator<< ( std::ostream &os, RigidBodyPerturbMover const &mover );


class RigidBodyPerturbRandomJumpMover : public RigidBodyMover{
public:
	RigidBodyPerturbRandomJumpMover();

	RigidBodyPerturbRandomJumpMover(
		core::Real const& rot_mag_in,
		core::Real const& trans_mag_in,
		core::Size const& num_jump_in);
	void apply(core::pose::Pose& pose) override;
	std::string get_name() const override;

	~RigidBodyPerturbRandomJumpMover() override;

private:
	core::Real rot_mag_in_;
	core::Real trans_mag_in_;
	core::Size num_jump_;
};

/// @brief does a perturbation defined by the rotational and translational magnitudes
///  without setting up the center
///  Can be defined through a move map or with rb_jump
///  Defining through a movemap with multiple jumps leads to a random jump being
///  chosen at apply time, NOT at construction time! This is done to simplify
///  docking with more than one active jump.
class RigidBodyPerturbNoCenterMover : public RigidBodyMover{
	typedef RigidBodyMover Parent;
public:
	// default constructor
	RigidBodyPerturbNoCenterMover();

	// constructor with arguments (rb_jump defined)
	RigidBodyPerturbNoCenterMover(
		int const rb_jump_in,
		core::Real const rot_mag_in,
		core::Real const trans_mag_in
	);

	// constructor with arguments (rb_jump not defined)
	// movemap used instead, to allow multiple jumps
	RigidBodyPerturbNoCenterMover(
		core::pose::Pose const & pose_in,
		core::kinematics::MoveMap const & mm,
		core::Real const rot_mag_in,
		core::Real const trans_mag_in,
		Direction dir_in
	);

	// function for the parser with lots of accessors
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	RigidBodyPerturbNoCenterMover( RigidBodyPerturbNoCenterMover const & src );
	~RigidBodyPerturbNoCenterMover() override;

	void apply( core::pose::Pose & pose ) override;
	void rot_magnitude( core::Real const magnitude ) { rot_mag_ = magnitude; }
	void trans_magnitude( core::Real const magnitude ) { trans_mag_ = magnitude; }

	void add_jump( core::Size );
	void clear_jumps();

	moves::MoverOP clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	/// perturbation magnitudes (rotational and translational)
	core::Real rot_mag_;
	core::Real trans_mag_;
private:
	typedef utility::vector1< core::Size > JumpList;
	JumpList movable_jumps_;
};  // class RigidBodyPerturbNoCenterMover


class RigidBodyRandomizeMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyRandomizeMover();

	// constructor with arguments
	RigidBodyRandomizeMover(
		core::pose::Pose const & pose_in,
		int const rb_jump_in=1,
		Partner const partner_in=partner_downstream,
		int phi_angle=360,
		int psi_angle=360,
		bool update_center_after_move=true
	);

	RigidBodyRandomizeMover( RigidBodyRandomizeMover const & );
	~RigidBodyRandomizeMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	void show(std::ostream & output=std::cout) const override;

	core::Size get_phi() const;
	core::Size get_psi() const;

private:
	Partner partner_; // which partner gets randomized
	core::Size phi_angle_;
	core::Size psi_angle_;
	numeric::xyzMatrix_double rotation_matrix_;
	bool update_center_after_move_;
};  // class RigidBodyRandomizeMover

std::ostream &operator<< ( std::ostream &os, RigidBodyRandomizeMover const &randommover );


/// @brief A Mover that spins about a random axis.
class RigidBodySpinMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodySpinMover();

	/// @brief constructor with arguments
	///       spin axis is initialized to 0 and then calculated during apply()
	RigidBodySpinMover( int const rb_jump_in );

	RigidBodySpinMover( RigidBodySpinMover const & src );
	~RigidBodySpinMover() override;

	void spin_axis( core::Vector spin_axis_in );
	void rot_center( core::Vector const & rot_center_in );

	void spin_mag( core::Real const & spin_mag );

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

protected:
	core::Vector spin_axis_;
	bool update_spin_axis_;

	core::Real spin_mag_;
	bool default_spin_mag_;

};  // class RigidBodySpinMover


class RigidBodyDeterministicSpinMover : public RigidBodySpinMover {
public:
	typedef RigidBodySpinMover parent;

public:
	//default ctor
	RigidBodyDeterministicSpinMover();

	/// @brief constructor with arguments
	/// spin axis is initialized to 0 then calculated during apply()
	/// if spin_axis is not already set
	RigidBodyDeterministicSpinMover( int const rb_jump_in, core::Vector spin_axis, core::Vector rotation_center, float angle_magnitude );

	//copy ctor
	RigidBodyDeterministicSpinMover( RigidBodyDeterministicSpinMover const & src );

	//dtor
	~RigidBodyDeterministicSpinMover() override;
	void angle_magnitude( float angle_magnitude );
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	float angle_magnitude_;
};  // RigidBodyDeterministicSpinMover

////////////////////////////////////////////////////////////////////////////////

/// @brief A Mover that tilts around the spin axis.
class RigidBodyTiltMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyTiltMover();

	/// @brief constructor with arguments
	///       spin axis is initialized to 0 and then calculated during apply()
	RigidBodyTiltMover(
		int const rb_jump_in,
		core::Real const tilt1_mag_in, //max tilt of partner1 in degrees
		core::Real const tilt2_mag_in, //max tilt of partner2 in degrees
		core::Size const tilt1_center_in, //resID for residue tilt of partner1 is centered at
		core::Size const tilt2_center_in //resID for residue tilt of partner2 is centered at
	);

	RigidBodyTiltMover( RigidBodyTiltMover const & src );
	~RigidBodyTiltMover() override;

	void spin_axis( core::Vector const & spin_axis_in );
	void tilt1_mag(core::Real const tilt1_mag_in ){tilt1_mag_ =tilt1_mag_in;}
	void tilt2_mag(core::Real const tilt2_mag_in ){tilt2_mag_ =tilt2_mag_in;}
	void tilt1_center(core::Size const tilt1_center_in ){tilt1_center_ =tilt1_center_in;}
	void tilt2_center(core::Size const tilt2_center_in ){tilt2_center_ =tilt2_center_in;}

	void apply( core::pose::Pose & pose ) override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:

	void tilt(
		core::pose::Pose & pose,
		std::string const & which,
		core::Real tilt_mag,
		core::Vector const & tilt_center,
		core::Vector const & spin_axis,
		core::Vector const & partner_center,
		Direction dir,
		core::kinematics::Stub const & upstream_stub,
		core::kinematics::Stub const & downstream_stub ) const;

	core::Vector find_tilt_center(
		core::pose::Pose const & pose,
		core::Size tilt_center_res,
		core::Vector const & partner_center) const;

protected:
	core::Real tilt1_mag_;
	core::Real tilt2_mag_;
	core::Size tilt1_center_;
	core::Size tilt2_center_;
	core::Vector spin_axis_;

};  // class RigidBodyTiltMover


/// @brief This Mover translate down an axis.
class RigidBodyTransMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyTransMover( bool vary_stepsize=false );

	// constructor with arguments
	RigidBodyTransMover(
		core::pose::Pose const & pose_in,
		int const rb_jump_in=1,
		bool vary_stepsize=false
	);

	// constructor with arguments that specify the trans axis
	RigidBodyTransMover( core::Vector const & trans_axis, int const rb_jump_in=1, bool vary_stepsize=false );


	~RigidBodyTransMover() override;

	RigidBodyTransMover( RigidBodyTransMover const & src );
	core::Vector & trans_axis() { return trans_axis_; }
	core::Vector trans_axis() const { return trans_axis_; }
	void trans_axis( core::Vector trans_axis_in ) { trans_axis_ = trans_axis_in; }

	void step_size( core::Real step_size_in ) { step_size_ = step_size_in; }
	core::Real step_size() { return step_size_ ;}
	void vary_stepsize( bool vary ) { vary_stepsize_ = vary; }

	void apply( core::pose::Pose & pose ) override;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Vector centroid_axis(core::pose::Pose const & pose_in) const;

	core::Real step_size_;
	core::Vector trans_axis_;

	bool vary_stepsize_;

};  // class RigidBodyTransMover

////////////////////////////////////////////////////////////////////////////////

/// @brief Rigid-body move that evenly samples the space within a sphere
class UniformSphereTransMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	UniformSphereTransMover();

	// constructor with arguments
	UniformSphereTransMover(
		int const rb_jump_in,
		core::Real step_size_in
	);

	UniformSphereTransMover( UniformSphereTransMover const & );
	~UniformSphereTransMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	void reset_trans_axis();

private:
	core::Real step_size_;
	core::Real random_step_; // by saving these we can apply the same random step to other things after we freeze
	core::Vector trans_axis_;// by saving these we can apply the same random step to other things after we freeze
};  // class UniformSphereTransMover


/// @brief A Mover that initializes all DOFs in the system randomly. It starts with rotation angles only.
class RigidBodyDofRandomizeMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyDofRandomizeMover();

	// constructor with arguments
	RigidBodyDofRandomizeMover(
		int const rb_jump_in,
		core::conformation::symmetry::SymDof
	);

	RigidBodyDofRandomizeMover( RigidBodyDofRandomizeMover const & src );
	~RigidBodyDofRandomizeMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	core::conformation::symmetry::SymDof dof_;
}; // class RigidBodyDofRandomizeMover


/// @brief A Mover that initializes all DOFs in the system randomly. It starts with rotation angles only.
class RigidBodyDofSeqRandomizeMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyDofSeqRandomizeMover();

	// constructor with arguments
	RigidBodyDofSeqRandomizeMover(
		std::map< Size, core::conformation::symmetry::SymDof > const & dofs
	);

	RigidBodyDofSeqRandomizeMover( RigidBodyDofSeqRandomizeMover const & );
	~RigidBodyDofSeqRandomizeMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	std::map< Size, core::conformation::symmetry::SymDof > dofs_;

};


/// @brief A Mover that translates down an axis determined by the available DOFs.
/// Translations are made along all allowed directions (x,y or z) for a selected jump.
class RigidBodyDofTransMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyDofTransMover();

	// constructor with arguments
	RigidBodyDofTransMover(
		core::conformation::symmetry::SymDof dof,
		int const rb_jump_in,
		core::Real step_size
	);

	// constructor with arguments
	RigidBodyDofTransMover(
		std::map< Size, core::conformation::symmetry::SymDof > dofs
	);

	RigidBodyDofTransMover( RigidBodyDofTransMover const & );
	~RigidBodyDofTransMover() override;

	core::Vector & trans_axis() { return trans_axis_; }
	core::Vector trans_axis() const { return trans_axis_; }
	void trans_axis( core::Vector trans_axis_in ) { trans_axis_ = trans_axis_in; }

	void step_size( core::Real step_size_in ) { step_size_ = step_size_in; }

	bool last_slide_good( ) { return last_slide_good_; }

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	bool last_slide_good_;
	int jump_dir_;
	core::Real step_size_;
	core::Vector trans_axis_;
	core::conformation::symmetry::SymDof dof_;
}; // class RigidBodyDofTransMover


/// @brief A Mover that translates down an axis determined by the available DOFs.
/// Translations are made along all allowed directions (x,y or z) for a selected jump.
/// Jumps are visited in random order
class RigidBodyDofSeqTransMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyDofSeqTransMover();

	// constructor with arguments
	RigidBodyDofSeqTransMover(
		std::map< Size, core::conformation::symmetry::SymDof > dofs
	);

	RigidBodyDofSeqTransMover( RigidBodyDofSeqTransMover const & );
	~RigidBodyDofSeqTransMover() override;

	core::Vector & trans_axis() { return trans_axis_; }
	core::Vector trans_axis() const { return trans_axis_; }
	void trans_axis( core::Vector trans_axis_in ) { trans_axis_ = trans_axis_in; }
	void step_size( core::Real step_size_in ) { step_size_ = step_size_in; }

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	// allowed dofs
	std::map< Size, core::conformation::symmetry::SymDof > dofs_;
	// allowed jumps
	utility::vector1 < int > rb_jumps_;
	core::Real step_size_;
	// silly, stores the direction only!!!
	core::Vector trans_axis_;
}; // class RigidBodyDofSeqTransMover


/// @brief A Mover that translates down an axis determined by the available DOFs.
/// Translations are made along all allowed directions (x,y or z) for a randomly selected jump.
class RigidBodyDofRandomTransMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyDofRandomTransMover();

	// constructor with arguments
	RigidBodyDofRandomTransMover(
		std::map< Size, core::conformation::symmetry::SymDof > dofs
	);

	RigidBodyDofRandomTransMover( RigidBodyDofRandomTransMover const & src );
	~RigidBodyDofRandomTransMover() override;

	core::Vector & trans_axis() { return trans_axis_; }
	core::Vector trans_axis() const { return trans_axis_; }
	void trans_axis( core::Vector trans_axis_in ) { trans_axis_ = trans_axis_in; }
	void step_size( core::Real step_size_in ) { step_size_ = step_size_in; }

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	/// allowed dofs
	std::map< Size, core::conformation::symmetry::SymDof > dofs_;
	// allowed jumps
	utility::vector1 < int > rb_jumps_;
	core::Real step_size_;
	// silly, stores the direction only!!!
	core::Vector trans_axis_;
}; // class RigidBodyDofRandomTransMover


/// @brief This Mover does a perturbation defined by the rotational and translational magnitudes.
///   Allowed dofs are specified by a map.
///   Can be defined through a move map or with rb_jump. A single jump is selected.
class RigidBodyDofPerturbMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor makes no sense at all!!!
	//  RigidBodyDofPerturbMover() : RigidBodyMover(), rot_mag_( 1.0 ), trans_mag_( 3.0 )
	//  {
	//    moves::Mover::type( "RigidBodyDofPerturbMover" );
	//  }

	// constructor with arguments (rb_jump not defined)
	// movemap used instead
	RigidBodyDofPerturbMover(
		std::map< Size, core::conformation::symmetry::SymDof > dofs,
		core::Real const rot_mag_in = 1.0,
		core::Real const trans_mag_in = 3.0
	);

	RigidBodyDofPerturbMover(
		int const rb_jump_in,
		core::conformation::symmetry::SymDof dof,
		core::Real const rot_mag_in = 1.0,
		core::Real const trans_mag_in = 3.0
	);

	RigidBodyDofPerturbMover( RigidBodyDofPerturbMover const & );
	~RigidBodyDofPerturbMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void rot_magnitude( core::Real const magnitude ) { rot_mag_ = magnitude; }
	void trans_magnitude( core::Real const magnitude ) { trans_mag_ = magnitude; }
	void dof( core::conformation::symmetry::SymDof dof ) { dof_ = dof; }

private:
	/// allowed dofs
	core::conformation::symmetry::SymDof dof_;
	/// perturbation magnitudes (rotational and translational)
	core::Real rot_mag_;
	core::Real trans_mag_;
}; // class RigidBodyDofPerturbMover


/// @brief This Mover does a perturbation defined by the rotational and translational magnitudes.
///   Allowed dofs are specified by a map.
///   Can be defined through a move map or with rb_jump. All jumps are selected in random order.
class RigidBodyDofSeqPerturbMover : public RigidBodyMover{
public:
	typedef RigidBodyMover parent;

public:
	// constructor with arguments (rb_jump not defined)
	// movemap used instead
	RigidBodyDofSeqPerturbMover(
		std::map< Size, core::conformation::symmetry::SymDof > dofs,
		core::Real const rot_mag_in = 1.0,
		core::Real const trans_mag_in = 3.0
	);
	RigidBodyDofSeqPerturbMover( RigidBodyDofSeqPerturbMover const & );
	~RigidBodyDofSeqPerturbMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void rot_magnitude( core::Real const magnitude ) { rot_mag_ = magnitude; }
	void trans_magnitude( core::Real const magnitude ) { trans_mag_ = magnitude; }
	void dofs( std::map< Size, core::conformation::symmetry::SymDof > dofs ) { dofs_ = dofs; }

private:
	// allowed dofs
	std::map< Size, core::conformation::symmetry::SymDof > dofs_;
	// allowed jumps
	utility::vector1 < int > rb_jumps_;
	/// perturbation magnitudes (rotational and translational)
	core::Real rot_mag_;
	core::Real trans_mag_;
}; // class RigidBodyDofSeqPerturbMover

} // rigid
} // protocols

#endif //INCLUDED_protocols_rigid_RigidBodyMover_HH
