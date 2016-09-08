// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/StructPerturberCMCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/StructPerturberCM.hh>
#include <protocols/abinitio/abscript/StructPerturberCMCreator.hh>

// Package headers
#include <core/environment/DofPassport.hh>

#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>

#include <core/kinematics/MoveMap.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>


//Utility Headers
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

#include <boost/foreach.hpp>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

//Req'd on WIN32
#include <basic/datacache/WriteableCacheableMap.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.abinitio.abscript.StructPerturberCM", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;

// creator
std::string
StructPerturberCMCreator::keyname() const {
	return StructPerturberCMCreator::mover_name();
}

protocols::moves::MoverOP
StructPerturberCMCreator::create_mover() const {
	return protocols::moves::MoverOP( new StructPerturberCM );
}

std::string
StructPerturberCMCreator::mover_name() {
	return "StructPerturberCM";
}

StructPerturberCM::StructPerturberCM():
	Parent()
{}

StructPerturberCM::StructPerturberCM( std::string const& label,
	core::Real magnitude ):
	magnitude_( magnitude ),
	label_( label )
{}

claims::EnvClaims StructPerturberCM::yield_claims( core::pose::Pose const& in_pose,
	basic::datacache::WriteableCacheableMapOP ){
	claims::EnvClaims claims;

	claims::TorsionClaimOP claim( new claims::TorsionClaim(
		utility::pointer::static_pointer_cast< ClientMover > ( get_self_ptr() ),
		label(), std::make_pair( 1, in_pose.size() ) ) );
	claim->strength( claims::CAN_CONTROL, claims::DOES_NOT_CONTROL );

	claims.push_back( claim );

	return claims;
}

void StructPerturberCM::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap&,
	protocols::filters::Filters_map const&,
	protocols::moves::Movers_map const&,
	core::pose::Pose const& ){
	magnitude( tag->getOption< core::Real >( "magnitude", 2.0 ) );
	label( tag->getOption< std::string >( "label", "BASE" ) );
}

void StructPerturberCM::apply( core::pose::Pose& pose ){
	if ( passport() ) {
		DofUnlock unlock( pose.conformation(), passport() );
		core::kinematics::MoveMapOP mm = passport()->render_movemap();

		for ( Size i = 1; i <= pose.size(); ++i ) {
			for ( DofPassport::const_iterator it = passport()->begin();
					it != passport()->end(); ++it ) {
				pose.set_dof( *it, pose.dof( *it) + ( numeric::random::rg().gaussian() * magnitude_ ) );
			}
		}
	} else {
		for ( Size i = 1; i <= pose.size(); ++i ) {
			pose.set_phi( i, pose.phi( i ) + numeric::random::rg().gaussian() * magnitude_ );
			pose.set_psi( i, pose.psi( i ) + numeric::random::rg().gaussian() * magnitude_ );
		}
	}
}

std::string StructPerturberCM::get_name() const {
	return "StructPerturberCM";
}

protocols::moves::MoverOP
StructPerturberCM::clone() const {
	return protocols::moves::MoverOP( new StructPerturberCM( *this ) );
}

} // abscript
} // abinitio
} // protocols
