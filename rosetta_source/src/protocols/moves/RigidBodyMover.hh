// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Monica Berrondo
/// @author Ian Davis
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_moves_RigidBodyMover_hh
#define INCLUDED_protocols_moves_RigidBodyMover_hh

// Unit headers
#include <protocols/moves/RigidBodyMover.fwd.hh>

// Package headers
#include <protocols/moves/ThermodynamicMover.hh>
//#include <protocols/geometry/RB_geometry.hh>

#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/conformation/symmetry/SymDof.hh>

// ObjexxFCL Headers

// C++ Headers
#include <map>

// Utility Headers
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace protocols {
namespace moves {

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
class RigidBodyMover : public ThermodynamicMover {
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

	virtual ~RigidBodyMover();

	/// @brief Manual override of rotation center.
	void rot_center( core::Vector const rot_center_in )	{	rot_center_ = rot_center_in; }

	// Can't change jump number after creation b/c it determines center of rotation!
	//void rb_jump_( int const rb_jump_in ) { rb_jump_ = rb_jump_in; }
	int rb_jump() const { return rb_jump_; }

	virtual void apply( core::pose::Pose & pose ) = 0;
	virtual std::string get_name() const;

	virtual
	bool
	preserve_detailed_balance() const { return true; }

	/// @brief set whether detailed balance is preserved (i.e. no branch angle optimization during moves)
	virtual
	void
	set_preserve_detailed_balance(
		bool preserve_detailed_balance
	) {};

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges( core::pose::Pose & pose );


	// data
protected:
	int rb_jump_;
	/// direction of folding (n-term to c-term)
	Direction dir_;

	/// center of rotation
	core::Vector rot_center_;
};

// does a perturbation defined by the rotational and translational magnitudes
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
		utility::vector1< bool > ok_for_centroid_calculation
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
	virtual ~RigidBodyPerturbMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;




	void rot_magnitude( core::Real const magnitude ) { rot_mag_ = magnitude; }

	void trans_magnitude( core::Real const magnitude ) { trans_mag_ = magnitude; }

	/// @brief Manual override of rotation center.
	void rot_center( core::Vector const /*rot_center_in*/ );

protected:
	/// perturbation magnitudes (rotational and translational)
	core::Real rot_mag_;
	core::Real trans_mag_;
private:
	Partner partner_;
	bool interface_;
	utility::vector1<core::Size> movable_jumps_;
	utility::vector1< bool > ok_for_centroid_calculation_;
};

///@brief does a perturbation defined by the rotational and translational magnitudes
/// 	without setting up the center
///		Can be defined through a move map or with rb_jump
///		Defining through a movemap with multiple jumps leads to a random jump being
///		chosen at apply time, NOT at construction time! This is done to simplify
///		docking with more than one active jump.
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
			 utility::tag::TagPtr const tag,
			 protocols::moves::DataMap &,
			 protocols::filters::Filters_map const &,
			 protocols::moves::Movers_map const &,
			 core::pose::Pose const &
	);

	RigidBodyPerturbNoCenterMover( RigidBodyPerturbNoCenterMover const & src );
	virtual ~RigidBodyPerturbNoCenterMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void rot_magnitude( core::Real const magnitude ) { rot_mag_ = magnitude; }
	void trans_magnitude( core::Real const magnitude ) { trans_mag_ = magnitude; }

	void add_jump( core::Size );
	void clear_jumps();
protected:
	/// perturbation magnitudes (rotational and translational)
	core::Real rot_mag_;
	core::Real trans_mag_;
private:
	utility::vector1< core::Size > movable_jumps_;
};

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
		int psi_angle=360
	);

	RigidBodyRandomizeMover( RigidBodyRandomizeMover const & );
	virtual ~RigidBodyRandomizeMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	Partner partner_; // which partner gets randomized
	core::Size phi_angle_;
	core::Size psi_angle_;

};

// spin about a random axis
class RigidBodySpinMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodySpinMover();

	///@brief constructor with arguments
	///       spin axis is initialized to 0 and then calculated during apply()
	RigidBodySpinMover( int const rb_jump_in );

	RigidBodySpinMover( RigidBodySpinMover const & src );
	~RigidBodySpinMover();

	void spin_axis( core::Vector spin_axis_in );
	void rot_center( core::Vector const rot_center_in );

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::Vector spin_axis_;
	bool update_spin_axis_;

};

// translate down an axis
class RigidBodyTransMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor
	RigidBodyTransMover();

	// constructor with arguments
	RigidBodyTransMover(
		core::pose::Pose const & pose_in,
		int const rb_jump_in=1
	);

	RigidBodyTransMover( RigidBodyTransMover const & src );
	virtual ~RigidBodyTransMover();

	core::Vector & trans_axis() { return trans_axis_; }
	core::Vector trans_axis() const { return trans_axis_; }
	void trans_axis( core::Vector trans_axis_in ) { trans_axis_ = trans_axis_in; }

	void step_size( core::Real step_size_in ) { step_size_ = step_size_in; }

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::Real step_size_;
	core::Vector trans_axis_;

};


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
	~UniformSphereTransMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::Real step_size_;

};

// Initialize all dofs in the system randomly. Start by rotation angles
// only.
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
	virtual ~RigidBodyDofRandomizeMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::conformation::symmetry::SymDof dof_;
};

// Initialize all dofs in the system randomly. Start by rotation angles
// only.
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
	virtual ~RigidBodyDofSeqRandomizeMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	std::map< Size, core::conformation::symmetry::SymDof > dofs_;

};


// Translate down axis determined by the available dofs.
// Translations are made along all allowed directions (x,y or z)
// for a selected jump
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
	~RigidBodyDofTransMover();

	core::Vector & trans_axis() { return trans_axis_; }
	core::Vector trans_axis() const { return trans_axis_; }
	void trans_axis( core::Vector trans_axis_in ) { trans_axis_ = trans_axis_in; }

	void step_size( core::Real step_size_in ) { step_size_ = step_size_in; }

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	int jump_dir_;
	core::Real step_size_;
	core::Vector trans_axis_;

};

// Translate down axis determined by the available dofs.
// Translation are made along all allowed directions (x,y or z)
// for all available jumps. Jumps are visited in random order
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
	~RigidBodyDofSeqTransMover();

	core::Vector & trans_axis() { return trans_axis_; }
	core::Vector trans_axis() const { return trans_axis_; }
	void trans_axis( core::Vector trans_axis_in ) { trans_axis_ = trans_axis_in; }
	void step_size( core::Real step_size_in ) { step_size_ = step_size_in; }

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	 /// allowed dofs
  std::map< Size, core::conformation::symmetry::SymDof > dofs_;
	// allowed jumps
  utility::vector1 < int > rb_jumps_;
	core::Real step_size_;
	// silly, stores the direction only!!!
	core::Vector trans_axis_;
};

// Translate down axis determined by the available dofs.
// Translation are made along all allowed directions (x,y or z)
// for a randomly selected jump.
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
	~RigidBodyDofRandomTransMover();

	core::Vector & trans_axis() { return trans_axis_; }
	core::Vector trans_axis() const { return trans_axis_; }
	void trans_axis( core::Vector trans_axis_in ) { trans_axis_ = trans_axis_in; }
	void step_size( core::Real step_size_in ) { step_size_ = step_size_in; }

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	 /// allowed dofs
	std::map< Size, core::conformation::symmetry::SymDof > dofs_;
	// allowed jumps
	utility::vector1 < int > rb_jumps_;
	core::Real step_size_;
	// silly, stores the direction only!!!
	core::Vector trans_axis_;
};




///@brief does a perturbation defined by the rotational and translational magnitudes
///   Allowed dofs are specified by a map
///   Can be defined through a move map or with rb_jump. A single jump is selected
class RigidBodyDofPerturbMover : public RigidBodyMover {
public:
	typedef RigidBodyMover parent;

public:
	// default constructor makes no sense at all!!!
	//  RigidBodyDofPerturbMover() : RigidBodyMover(), rot_mag_( 1.0 ), trans_mag_( 3.0 )
	//  {
	//    Mover::type( "RigidBodyDofPerturbMover" );
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
	virtual ~RigidBodyDofPerturbMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void rot_magnitude( core::Real const magnitude ) { rot_mag_ = magnitude; }
	void trans_magnitude( core::Real const magnitude ) { trans_mag_ = magnitude; }
	void dof( core::conformation::symmetry::SymDof dof ) { dof_ = dof; }

private:
	/// allowed dofs
	core::conformation::symmetry::SymDof dof_;
  /// perturbation magnitudes (rotational and translational)
  core::Real rot_mag_;
  core::Real trans_mag_;
};

///@brief does a perturbation defined by the rotational and translational magnitudes
///   Allowed dofs are specified by a map
///   Can be defined through a move map or with rb_jump. All jumps are selected in random order
class RigidBodyDofSeqPerturbMover : public RigidBodyMover{
public:
	typedef RigidBodyMover parent;

public:
  // default constructor makes no sense at all!!!
//  RigidBodyDofSeqPerturbMover() : RigidBodyMover(), rot_mag_( 3.0 ), trans_mag_( 8.0 )
//  {
//    Mover::type( "RigidBodyDofSeqPerturbMover" );
//  }

	// constructor with arguments (rb_jump not defined)
	// movemap used instead
	RigidBodyDofSeqPerturbMover(
		std::map< Size, core::conformation::symmetry::SymDof > dofs,
		core::Real const rot_mag_in = 1.0,
		core::Real const trans_mag_in = 3.0
	);
	RigidBodyDofSeqPerturbMover( RigidBodyDofSeqPerturbMover const & );
	virtual ~RigidBodyDofSeqPerturbMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void rot_magnitude( core::Real const magnitude ) { rot_mag_ = magnitude; }
	void trans_magnitude( core::Real const magnitude ) { trans_mag_ = magnitude; }
	void dofs( std::map< Size, core::conformation::symmetry::SymDof > dofs ) { dofs_ = dofs; }

private:
	/// allowed dofs
	std::map< Size, core::conformation::symmetry::SymDof > dofs_;
	// allowed jumps
	utility::vector1 < int > rb_jumps_;
  /// perturbation magnitudes (rotational and translational)
  core::Real rot_mag_;
  core::Real trans_mag_;
};

} // moves
} // protocols


#endif //INCLUDED_protocols_moves_RigidBodyMover_HH
