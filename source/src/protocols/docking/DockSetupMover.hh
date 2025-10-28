// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief allows low-resolution docking using simulated or parallel tempering
/// @author Oliver Lange (oliver.lange@tum.de)

#ifndef INCLUDED_protocols_docking_DockSetupMover_hh
#define INCLUDED_protocols_docking_DockSetupMover_hh

// Unit Headers
#include <protocols/docking/DockSetupMover.fwd.hh>

// Package Headers
#include <protocols/docking/types.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/DockingPartners.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/docking/RigidBodyInfo.fwd.hh>

#include <utility/tag/Tag.fwd.hh>


// Utility Headers

// Numeric Headers

// ObjexxFCL Headers

// C++ headers

namespace protocols {
namespace docking {

/// @brief allows docking using simulated or parallel tempering
/// @details
class DockSetupMover : public moves::Mover
{
public:
	/// @brief default constructor fills values with the expected defaults
	DockSetupMover();

	/// @brief clone
	protocols::moves::MoverOP clone() const override;

	/// @brief copy ctor
	DockSetupMover( DockSetupMover const & rhs );

	/// @brief assignment operator
	DockSetupMover & operator=( DockSetupMover const & rhs );

	/// @brief Assigns default values to primitive members
	void set_defaults();

	void apply( core::pose::Pose & pose ) override;

	core::pose::DockingPartners const& partners() const { return partners_;} /// @brief returns the docking partners chain identifiers


	// DockJumps & movable_jumps(){ return movable_jumps_;} /// @brief returns ref to the jumps vector for docking
	// DockJumps const & movable_jumps() const { return movable_jumps_; } /// @ return const ref to the jumps vector for docking

	void set_partners( core::pose::DockingPartners const& setting ){ partners_=setting; }
	// void set_movable_jumps( DockJumps const& setting ){ movable_jumps_ = setting; }
	// void add_jump( core::SSize const jump_number ){ movable_jumps_.push_back( int( jump_number ) ); }

	void show( std::ostream & out=std::cout ) const override;
	friend std::ostream & operator<<(std::ostream& out, const DockSetupMover & dp );

	// function for the parser with lots of accessors
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	void copy( DockSetupMover & lhs, DockSetupMover const & rhs);

private:
	/// --- configurables -----
	core::pose::DockingPartners partners_;
	protocols::rigid::RigidBodyPerturbNoCenterMoverOP rb_mover_;
	DockJumps movable_jumps_; //vector1_int
	protocols::docking::RigidBodyInfoOP rigid_body_info_;

	/// --- state ----
	// core::kinematics::FoldTree fold_tree_;
	// std::string previous_sequence_;

};

} // docking
} // protocols

#endif

