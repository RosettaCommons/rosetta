// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rigid/UniformRigidBodyMover.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_rigid_UniformRigidBodyMover_hh
#define INCLUDED_protocols_rigid_UniformRigidBodyMover_hh

// Unit Headers
#include <protocols/rigid/UniformRigidBodyMover.fwd.hh>
#include <protocols/environment/ClientMover.hh>

// Package headers
#include <protocols/environment/claims/EnvClaim.hh>

#include <protocols/canonical_sampling/ThermodynamicMover.hh>

// Project headers
#include <core/pose/Pose.hh>

#include <basic/datacache/WriteableCacheableMap.fwd.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace rigid {

class UniformRigidBodyMover : public canonical_sampling::ThermodynamicMover {
	typedef int JumpNumber;

public:
	UniformRigidBodyMover();

	UniformRigidBodyMover( JumpNumber target_jump,
		core::Real rotation_mag_ = 3.0,
		core::Real translation_mag_ = 8.0 );


	~UniformRigidBodyMover() override = default;


	void apply( core::pose::Pose& ) override;

	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	moves::MoverOP fresh_instance() const override;


	moves::MoverOP clone() const override;

	void jump_number( JumpNumber );

	JumpNumber jump_number() const;


	void
	set_preserve_detailed_balance( bool ) override {};


	bool
	preserve_detailed_balance() const override { return true; }


	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges( core::pose::Pose & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	utility::tag::XMLSchemaComplexTypeGeneratorOP
	complex_type_gen();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	JumpNumber target_jump_;
	core::Real rotation_mag_, translation_mag_;

}; // end UniformRigidBodyMover base class

} // abinitio
} // protocols

#endif //INCLUDED_protocols_rigid_UniformRigidBodyMover_hh
