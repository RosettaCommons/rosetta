// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rigid/CoMTrackerCM.cc
/// @author Justin R. Porter
/// @author Brian D. Weitzner
/// @author Oliver F. Lange

// Unit Headers
#include <protocols/environment/CoMTrackerCM.hh>
#include <protocols/environment/CoMTrackerCMCreator.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <core/environment/LocalPosition.hh>

#include <protocols/environment/DofUnlock.hh>

#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/VirtResClaim.hh>

// Project headers
#include <core/id/NamedStubID.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/util.hh>

//Utility Headers
#include <utility>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

#include <utility/string_util.hh>
#include <utility/vector1.hh>

#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <basic/datacache/DataMap.hh>

// tracer
#include <basic/Tracer.hh>


#ifdef WIN32
#include <basic/datacache/WriteableCacheableMap.hh>
#endif


// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

std::string const GENERATE_STATIONARY_ATTACHMENT_POINT = "[NOT_SET]";

static THREAD_LOCAL basic::Tracer tr( "protocols.rigid.CoMTrackerCM", basic::t_info );

using namespace core::environment;
using namespace protocols::environment;

// creator
std::string
CoMTrackerCMCreator::keyname() const {
	return CoMTrackerCMCreator::mover_name();
}

protocols::moves::MoverOP
CoMTrackerCMCreator::create_mover() const {
	return ClientMoverOP( new CoMTrackerCM );
}

std::string
CoMTrackerCMCreator::mover_name() {
	return "CoMTrackerCM";
}

CoMTrackerCM::CoMTrackerCM():
	ClientMover()
{}

CoMTrackerCM::CoMTrackerCM( std::string  name,
	core::select::residue_selector::ResidueSelectorCOP mobile_selector,
	std::string  stationary_label ):
	ClientMover(),
	name_(std::move( name )),
	stationary_label_(std::move( stationary_label )),
	mobile_selector_(std::move( mobile_selector ))
{}

CoMTrackerCM::CoMTrackerCM( std::string  name,
	core::select::residue_selector::ResidueSelectorCOP mobile_selector ):
	ClientMover(),
	name_(std::move( name )),
	stationary_label_( GENERATE_STATIONARY_ATTACHMENT_POINT ),
	mobile_selector_(std::move( mobile_selector ))
{}

void CoMTrackerCM::passport_updated(){
	if ( this->has_passport() ) {
		EnvironmentCOP env( active_environment() );
		SequenceAnnotationCOP ann = env->annotations();

		// utility::vector1< core::Size > mobile_residues = ann->resolve_seq( mobile_label_ );
		// com_residues_ = std::set< core::Size >( mobile_residues.begin(), mobile_residues.end() );
	} else {
		//configure a null moveset.
	}
}

numeric::xyzVector< core::Real > com_calc( core::pose::Pose const& pose,
	utility::vector1_bool const& rsds ) {
	using namespace numeric;
	using namespace core::conformation;

	utility::vector1< xyzVector< core::Real > > coords;

	for ( core::Size i = 1; i <= rsds.size(); ++i ) {
		if ( rsds[i] ) {
			for ( auto it = pose.residue( i ).atom_begin();
					it != pose.residue( i ).heavyAtoms_end(); ++it ) {
				coords.push_back( it->xyz() );
			}
		}
	}

	assert( coords.size() > 0 );

	return center_of_mass( coords );
}

void CoMTrackerCM::update_tracking_residue( core::kinematics::RT::Vector new_position,
	core::Size tracking_residue_id,
	core::pose::Pose & pose ) const {
	using core::Size;
	using core::kinematics::Jump;
	using core::kinematics::RT;
	using core::kinematics::Stub;
	using utility::vector1;

	// TODO: Update this to guarantee test_point is downstream from the tracking_residue.
	// Get the initial position of a (presumably) downstream atom to compare after the tracking resiude is updated.
	if ( mobile_residues_.index( false ) == 0 ||
			mobile_residues_.index( false ) >= (int) pose.total_residue() ) {
		std::ostringstream ss;
		ss << "The CoMTrackerCM '" << name() << "' was configured to move all the residues in the pose.  "
			<< "This probably isn't what you meant, check your selectors and input pose." << std::endl;
		throw utility::excn::EXCN_BadInput( ss.str() );
	}
	assert( (Size) mobile_residues_.index( false ) <= pose.total_residue() );

	RT::Vector test_point = pose.residue( mobile_residues_.index( false ) ).xyz( 1 );

	// TODO: throw an excpetion here instead of using an assert
	// Make sure only virutal residues are being used as tracking residues
	assert( pose.residue( tracking_residue_id ).name() == "VRT" );
	assert( pose.fold_tree().root() != (int)tracking_residue_id );

	// By definition, the virtual residue we are placing must be the downstream partner in exactly one jump.
	// We refer to this jump as the "positioning_jump".
	// It can be the upstream partner in zero or more jumps, so the first thing we are going to do is list the jumps
	// in which the virtual residue is participating.
	// We refer to these jumps as "positioned_jumps".
	vector1< int > positioned_jump_ids;
	int positioning_jump_id = 0;

	for ( int i = 1; i <= (int)pose.fold_tree().num_jump(); ++i ) {
		if ( pose.fold_tree().downstream_jump_residue( i ) == (int)tracking_residue_id ) {
			assert( positioning_jump_id == 0 );
			positioning_jump_id = i;
		} else if ( pose.fold_tree().upstream_jump_residue( i ) == (int)tracking_residue_id ) {
			positioned_jump_ids.push_back( i );
		}
	}

	// TODO: throw an excpetion instead of using an assert
	// Make sure the tracking residue's jump was detected properly, then make the jump number const
	assert( positioning_jump_id != 0 );

	// Get the stubs from the jump that is responsible for positioning the tracking residue and set the center of the
	// tracking residue's stub to the new location.
	Stub tracking_res_stub = pose.conformation().downstream_jump_stub( positioning_jump_id );
	Stub upstream_postioning_stub = pose.conformation().upstream_jump_stub( positioning_jump_id );

	tracking_res_stub.v = new_position;

	// Update all of the jumps the tracking residue positions.
	// It is critical that this is done prior to actually adjusting the position of the tracking residue so the
	// downstream stubs are still in the correct positions.
	for ( vector1< int >::const_iterator it = positioned_jump_ids.begin(); it != positioned_jump_ids.end(); ++it ) {
		if ( passport()->has_jump_access( *it ) ) {
			//We don't nessecarily always have access to all jumps that are built by this guy, and sometimes that's ok.
			pose.set_jump( *it, Jump( RT( tracking_res_stub, pose.conformation().downstream_jump_stub( *it ) ) ) );
		} else {
			core::Size non_vrt_residue = pose.fold_tree().downstream_jump_residue( *it );
			if ( mobile_residues_.size() >= non_vrt_residue &&
					mobile_residues_[ non_vrt_residue ] == true ) {
				// In this case, it's not ok. We MUST adjust the jump if it builds
				assert( false );
			}
		}
	}

	// TODO: throw an excpetion here instead of using an assert
	// Update the position of the tracking residue and ensure that it ends up where we want it
	pose.set_jump( positioning_jump_id, Jump( RT( upstream_postioning_stub, tracking_res_stub ) ) );
	assert( pose.residue( tracking_residue_id ).xyz( "ORIG" ).distance_squared( new_position ) < 1e-10 );

	assert( pose.residue( mobile_residues_.index( false ) ).xyz( 1 ).distance_squared( test_point ) < 1e-10 );
}

void CoMTrackerCM::update_com( core::pose::Pose& pose ) const {
	using namespace numeric;
	using core::kinematics::Jump;
	using core::kinematics::RT;

	// active_environment() returns NULL without a valid Passport!
	EnvironmentCOP env( active_environment() );
	SequenceAnnotationCOP ann = env->annotations();

	// get the sequence number of the virtual residue that should track the CoM
	core::Size const vrt_resid = ann->resolve_seq( LocalPosition( com_name_, 1 ) );

	//calculate the center of mass in laboratory coordinate frame
	RT::Vector com = com_calc( pose, mobile_residues_ );

	update_tracking_residue( com, vrt_resid, pose );
}

void CoMTrackerCM::initialize( core::pose::Pose& pose ){
	DofUnlock activeation( pose.conformation(), passport() );
	core::environment::DofPassportCOP pp = passport();
	update_com( pose );
}

void CoMTrackerCM::apply( core::pose::Pose& pose ){
	DofUnlock activation( pose.conformation(), passport() );
	update_com( pose );
}

void CoMTrackerCM::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const&,
	protocols::moves::Movers_map const&,
	core::pose::Pose const& ) {

	//  mobile_label_ = tag->getOption< std::string >( "mobile_label" );
	stationary_label_ = tag->getOption< std::string >( "stationary_label", GENERATE_STATIONARY_ATTACHMENT_POINT );
	name_ = tag->getOption< std::string >( "name" );

	// the mobile selector is inverted during brokering to determine which residues are stationary
	using namespace core::select::residue_selector;
	mobile_selector_ = datamap.get_ptr< ResidueSelector >( "ResidueSelector", tag->getOption<std::string>( "mobile_selector" ) );


}

claims::EnvClaims CoMTrackerCM::yield_claims( core::pose::Pose const& pose,
	basic::datacache::WriteableCacheableMapOP ){
	using core::Size;
	using core::select::residue_selector::ResidueSubset;
	claims::EnvClaims claim_list;

	// com_name_ = mobile_label_ + "CoM";
	com_name_ = name();
	com_jump_name_ = name() + "_jump";

	// Get the position of the first residue in the "Mobile Selection"
	mobile_residues_ = mobile_selector_->apply( pose );
	Size mobile_connection_point = mobile_residues_.index( true )+1;

	if ( std::find( mobile_residues_.begin(), mobile_residues_.end(), true ) == mobile_residues_.end() ) {
		std::ostringstream ss;
		ss << "The mobile_selector for '" << this->get_name() << "' made an empty selection. This is not allowed.";
		throw utility::excn::EXCN_BadInput( ss.str() );
	}
	assert( std::find( mobile_residues_.begin(), mobile_residues_.end(), true ) != mobile_residues_.end() );

	moves::MoverOP this_ptr = get_self_ptr();
	claims::VirtResClaimOP vclaim( new claims::VirtResClaim( utility::pointer::static_pointer_cast< ClientMover >(this_ptr),
		LocalPosition( "BASE", mobile_connection_point ),
		com_jump_name_,
		com_name_ ) );

	vclaim->jump().strength( claims::MUST_CONTROL, claims::MUST_CONTROL );
	claim_list.push_back( vclaim );

	LocalPosition stationary_attchmnt_pt;
	if ( stationary_label_ == GENERATE_STATIONARY_ATTACHMENT_POINT ) {
		ResidueSubset stationary_residues = mobile_residues_.invert();
		stationary_attchmnt_pt = LocalPosition( "BASE", stationary_residues.index( true ) );
	} else {
		stationary_attchmnt_pt = LocalPosition( stationary_label_, 1 );
	}

	return claim_list;
}

std::string CoMTrackerCM::get_name() const {
	return "CoMTrackerCM('"+name_+"')";
}

moves::MoverOP CoMTrackerCM::fresh_instance() const {
	return ClientMoverOP( new CoMTrackerCM() );
}

moves::MoverOP CoMTrackerCM::clone() const{
	return ClientMoverOP( new CoMTrackerCM( *this ) );
}

} // rigid
} // protocols
