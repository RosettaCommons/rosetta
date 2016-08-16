// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  inheried from RigidBodyPerturbNoCenterMover, 1. use Rigid_Body_Info to parse the movable_jump to
/// get independant from DockSetupMover; 2. use unbiased rotation sampling and gaussing translation instead
/// of gaussian_move() in which rotation is not unbiased. 3. restrict the search space in RT level
/// @author Zhe Zhang


#include <devel/replica_docking/UnbiasedRigidBodyMover.hh>
#include <devel/replica_docking/UnbiasedRigidBodyMoverCreator.hh>

#include <protocols/docking/RigidBodyInfo.hh>
#include <basic/datacache/DataMap.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/tag/Tag.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>

//#include <numeric/rotation.functions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.io.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/trig.functions.hh>
#include <numeric/conversions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer tr( "devel.replica_docking.UnbiasedRigidBodyMover" );

namespace devel {
namespace replica_docking {

std::string
UnbiasedRigidBodyPerturbNoCenterMoverCreator::keyname() const {
	return UnbiasedRigidBodyPerturbNoCenterMoverCreator::mover_name();
}

protocols::moves::MoverOP
UnbiasedRigidBodyPerturbNoCenterMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new UnbiasedRigidBodyPerturbNoCenterMover );
}

std::string
UnbiasedRigidBodyPerturbNoCenterMoverCreator::mover_name() {
	return "UnbiasedRigidBodyPerturbNoCenter";
}

UnbiasedRigidBodyPerturbNoCenterMover::UnbiasedRigidBodyPerturbNoCenterMover() :
	rigid_body_info_( /* NULL */ ),
	initialized_( false ),
	restrict_( false ),
	max_move_( false )
{}

UnbiasedRigidBodyPerturbNoCenterMover::UnbiasedRigidBodyPerturbNoCenterMover( UnbiasedRigidBodyPerturbNoCenterMover const& other ) : RigidBodyPerturbNoCenterMover( other ) {
	///copy value of every private variables
	rigid_body_info_ = other.rigid_body_info_;
	movable_jumps_ = other.movable_jumps_;
	initialized_ = other.initialized_;
	restrict_ = other.restrict_;
	max_move_ = other.max_move_;
	ref_file_ = other.ref_file_;
	max_trans_dist_ = other.max_trans_dist_;
	max_rot_angle_ = other.max_rot_angle_;
}

UnbiasedRigidBodyPerturbNoCenterMover::~UnbiasedRigidBodyPerturbNoCenterMover() {}

std::string
UnbiasedRigidBodyPerturbNoCenterMover::get_name() const
{
	return "UnbiasedRigidBodyPerturbNoCenter";
}

protocols::moves::MoverOP
UnbiasedRigidBodyPerturbNoCenterMover::clone() const
{
	return protocols::moves::MoverOP( new UnbiasedRigidBodyPerturbNoCenterMover(*this) );
}

protocols::moves::MoverOP
UnbiasedRigidBodyPerturbNoCenterMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new UnbiasedRigidBodyPerturbNoCenterMover );
}

void
UnbiasedRigidBodyPerturbNoCenterMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {

	if ( !data.has( "RigidBodyInfo", "docking_setup" ) ) {
		tr << "RigidBodyInfo not found in basic::datacache::DataMap" << std::endl;
		rigid_body_info_ = protocols::docking::RigidBodyInfoOP( new protocols::docking::RigidBodyInfo );
		data.add( "RigidBodyInfo", "docking_setup", rigid_body_info_ );
	} else {
		rigid_body_info_ = data.get_ptr<protocols::docking::RigidBodyInfo>( "RigidBodyInfo", "docking_setup" );
		tr.Debug << "get RigidBodyInfo pointer from basic::datacache::DataMap" << std::endl;
	}

	if ( tag->hasOption( "restrict" ) && tag->getOption< bool >("restrict") ) {
		restrict_ = true;
		max_move_ = tag->getOption< bool >("max_move", false);
		max_trans_dist_ = tag->getOption<core::Real>( "restrict_trans", 9999.0);
		max_rot_angle_ = numeric::conversions::radians( tag->getOption< core::Real >("restrict_rot",180.0) );
		if ( tag->hasOption( "ref" ) ) {
			ref_file_ = tag->getOption< std::string >( "ref" );
		} else { /// if not specified, restrict the search space from the starting input structure
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			ref_file_ = option[ in::file::s ](1);
		}
		tr.Debug << "Searching space is restricted relative to " << ref_file_ << " by max_trans_dist: " << max_trans_dist_ << " and max_rot_angle (radian): " << max_rot_angle_ << std::endl;
	}

	Parent::parse_my_tag( tag, data, filters, movers, pose );
}


void UnbiasedRigidBodyPerturbNoCenterMover::initialize ( core::pose::Pose const & pose )
{
	tr.Debug << "initialize: " << std::endl;
	if ( rigid_body_info_ ) {
		movable_jumps_ = rigid_body_info_->movable_jumps();
		tr.Debug << "finished reading movable_jumps_ from RigidBodyInfo" << std::endl;
		if ( movable_jumps_.empty() ) {
			utility_exit_with_message( "DockSetupMover has to be applied before this !" );
		}
	}
	runtime_assert( !movable_jumps_.empty() );
	for ( protocols::docking::DockJumps::const_iterator it=movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {
		Parent::add_jump( *it );
		tr.Debug << "movable_jump " << *it << " was added!" << std::endl;
	}

	// set baseclass rb_jump_ randomly from list of movable jumps
	if ( movable_jumps_.size() > 1 ) {
		rb_jump_ = numeric::random::rg().random_element( movable_jumps_ );

	} else if ( movable_jumps_.size() == 1 ) {
		rb_jump_ = movable_jumps_[1];
		tr.Debug <<"set rb_jump_ movable_jumps_[1]" << std::endl;
	} else {
		rb_jump_ = 1;
		tr.Debug <<"set rb_jump_ 1 " << std::endl;
	}
	tr.Debug << "Set movable jump# " << rb_jump_ << std::endl;
	tr.Debug <<"rot_mag "<< rot_mag_ <<" ; trans_mag " << trans_mag_ << std::endl;

	if ( restrict_ ) {
		core::pose::PoseOP ref_pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *ref_pose, ref_file_ , core::import_pose::PDB_file);
		ref_pose->fold_tree( pose.fold_tree() ); // set its fold_tree same as the docking pose
		core::kinematics::Jump ref_jump = ref_pose->jump( rb_jump_ );
		ref_R_ = ref_jump.get_rotation();
		ref_T_ = ref_jump.get_translation();
		tr.Debug << "ref_T_: " << ref_T_.x() << " "<< ref_T_.y() << " "<< ref_T_.z() << std::endl;
	}
	initialized_ = true;
}

void UnbiasedRigidBodyPerturbNoCenterMover::apply( core::pose::Pose& pose ) {
	// Parent::apply( pose );
	// overload to use random unit quaternion to generate uniformly distributed rotation.
	using namespace numeric;
	using namespace numeric::random;

	if ( !initialized_ ) initialize( pose );

	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	xyzMatrix< core::Real> const rot=flexible_jump.get_rotation();
	xyzVector< core::Real> const trans=flexible_jump.get_translation();
	xyzVector<core::Real> delta_trans;
	xyzMatrix<core::Real> delta_rot;
	//  tr.Debug << "pose.trans: " << trans << std::endl;
	if ( !max_move_ ) {
		bool flag( true );
		while ( flag ) {
			delta_trans = random_translation( trans_mag_, numeric::random::rg() );
			if ( restrict_ ) {
				xyzVector<core::Real> diff_T = delta_trans + trans - ref_T_;
				//   tr.Debug << "diff_T: " << diff_T << std::endl;
				if ( diff_T.length() <= max_trans_dist_ ) flag = false;
			} else {
				flag = false;
			}
		}

		flag = true;
		while ( flag ) {
			core::Real theta = random_rotation_angle( rot_mag_, numeric::random::rg() ); // non-negative, [0 pi]
			xyzVector<core::Real> axis = random_point_on_unit_sphere< core::Real >( numeric::random::rg() );
			delta_rot = numeric::rotation_matrix_radians( axis, theta );
			if ( restrict_ ) {
				xyzMatrix<core::Real> diff_R = delta_rot * rot * ref_R_.transposed();
				core::Real diff_theta;
				rotation_axis( diff_R, diff_theta );
				//   tr.Debug << "diff_theta: " << diff_theta << std::endl;
				if ( diff_theta <= max_rot_angle_ ) flag = false;
			} else {
				flag = false;
			}
		}
	} else { // end multiple search moves
		/// single move with the maximum value of the restricted parameter
		/// used to generate initial input structures with the same Lrmsd from the native structure
		delta_trans = max_trans_dist_ * random_point_on_unit_sphere< core::Real >( numeric::random::rg() );
		xyzVector<core::Real> axis = random_point_on_unit_sphere< core::Real >( numeric::random::rg() );
		delta_rot = numeric::rotation_matrix_radians( axis, max_rot_angle_ );
	} // end max_move

	tr.Debug << "delta_trans: " << delta_trans << std::endl;
	tr.Debug << "delta_rot: " << delta_rot << std::endl;

	flexible_jump.set_translation( delta_trans + trans );
	flexible_jump.set_rotation( delta_rot*rot );
	//  flexible_jump.gaussian_trans_random_rotation( dir_, trans_mag_, rot_mag_);
	pose.set_jump( rb_jump_, flexible_jump );


}

}
}
