// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rigid/UniformRigidBodyCM.hh
/// @author Justin R. Porter
/// @author Brian D. Weitzner
/// @author Oliver F. Lange

#ifndef INCLUDED_protocols_rigid_UniformRigidBodyCM_hh
#define INCLUDED_protocols_rigid_UniformRigidBodyCM_hh

// Unit Headers
#include <protocols/rigid/UniformRigidBodyCM.fwd.hh>
#include <protocols/environment/ClientMover.hh>

// Package headers
#include <protocols/rigid/UniformRigidBodyMover.hh>

#include <protocols/environment/claims/EnvClaim.hh>

// Project headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.hh>


#include <basic/datacache/WriteableCacheableMap.fwd.hh>


// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace rigid {

class UniformRigidBodyCM : public environment::ClientMover {
	typedef core::environment::LocalPosition LocalPosition;
	typedef environment::claims::EnvClaims EnvClaims;
	typedef int JumpNumber;

public:
	UniformRigidBodyCM();

	UniformRigidBodyCM( std::string const& name,
		LocalPosition const& mobile,
		LocalPosition const& stationary,
		core::Real rotation_magnitude = 3.0,
		core::Real translation_magnitude = 8.0 );


	~UniformRigidBodyCM() override = default;


	EnvClaims yield_claims( core::pose::Pose const&,
		basic::datacache::WriteableCacheableMapOP ) override;


	void initialize( core::pose::Pose& pose ) override;

	void apply( core::pose::Pose& ) override;

	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string const& name() const { return name_; }

	void name( std::string const& name ) { name_ = name; }


	moves::MoverOP fresh_instance() const override;


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
	void passport_updated() override;

private:

	std::string name_;
	LocalPosition mobile_label_, stationary_label_;
	UniformRigidBodyMoverOP mover_;

}; // end UniformRigidBodyCM base class

} // abinitio
} // protocols

#endif //INCLUDED_protocols_rigid_UniformRigidBodyCM_hh
