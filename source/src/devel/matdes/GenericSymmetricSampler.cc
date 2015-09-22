// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file
/// @brief
/// @author Neil King ( neilking@uw.edu )
/// @author Javier Castellanos ( javiercv@uw.edu )
/// @author Frank DiMaio

// Unit headers
#include <devel/matdes/GenericSymmetricSampler.hh>
#include <devel/matdes/GenericSymmetricSamplerCreator.hh>

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <numeric/random/random.hh>

static THREAD_LOCAL basic::Tracer TR( "devel.matdes.GenericSymmetricSampler" );


namespace devel {
namespace matdes {

using namespace core;
using namespace utility;

//////////////////////////////////////////////
/// creator
std::string
GenericSymmetricSamplerCreator::keyname() const
{
	return GenericSymmetricSamplerCreator::mover_name();
}

protocols::moves::MoverOP
GenericSymmetricSamplerCreator::create_mover() const {
	return protocols::moves::MoverOP( new GenericSymmetricSampler );
}

std::string
GenericSymmetricSamplerCreator::mover_name()
{
	return "GenericSymmetricSampler";
}

//////////////////////////////////////////////
/// mover
GenericSymmetricSampler::GenericSymmetricSampler() :
	dof_id_(""),
	angle_min_(-1.0),
	angle_max_(1.0),
	angle_step_(1.0),
	radial_disp_min_(-1.0),
	radial_disp_max_(1.0),
	radial_disp_step_(1.0),
	usescore_(false),
	maxscore_(false)
{ }

protocols::moves::MoverOP
GenericSymmetricSampler::clone() const {
	return protocols::moves::MoverOP( new GenericSymmetricSampler( *this ) );
}

protocols::moves::MoverOP
GenericSymmetricSampler::fresh_instance() const {
	return protocols::moves::MoverOP( new GenericSymmetricSampler() );
}


void
GenericSymmetricSampler::apply(Pose & pose) {
	using core::conformation::symmetry::SymmetryInfoCOP;
	using core::conformation::symmetry::SymDof;
	using core::pose::Pose;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_BAD_INPUT;
	using protocols::moves::FAIL_RETRY;


	// for each movable dof sample
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);

	std::map<Size,SymDof> dofs = sym_info->get_dofs();
	int sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, dof_id_ );

	// grid search of possibly 6 DOFs
	utility::vector1 <Real> dof_min(6,0.0);
	utility::vector1 <Real> dof_max(6,0.0);

	for ( int dofid=1; dofid<=6; ++dofid ) {
		if ( dofs[sym_aware_jump_id].allow_dof(dofid) ) {
			TR << "move jump " << dof_id_ << " along RB " << dofid << std::endl;
			dof_min[dofid] = (dofid>3) ? angle_min_:radial_disp_min_;
			dof_max[dofid] = (dofid>3) ? angle_max_:radial_disp_max_;
		}
	}

	core::kinematics::Jump j = pose.jump(sym_aware_jump_id);
	j.fold_in_rb_deltas();

	Pose best_pose = pose;
	Real best_score = 1e20;
	protocols::moves::MoverStatus ms( FAIL_RETRY );

	for ( core::Real x=dof_min[1]; x<=dof_max[1]; x+=radial_disp_step_ ) {
		for ( core::Real y=dof_min[2]; y<=dof_max[2]; y+=radial_disp_step_ ) {
			for ( core::Real z=dof_min[3]; z<=dof_max[3]; z+=radial_disp_step_ ) {
				for ( core::Real a=dof_min[4]; a<=dof_max[4]; a+=radial_disp_step_ ) {
					for ( core::Real b=dof_min[5]; b<=dof_max[5]; b+=radial_disp_step_ ) {
						for ( core::Real g=dof_min[6]; g<=dof_max[6]; g+=radial_disp_step_ ) {
							core::kinematics::Jump j0 = j;
							j0.set_rb_delta(1,1,x); j0.set_rb_delta(2,1,y); j0.set_rb_delta(3,1,z);
							j0.set_rb_delta(4,1,a); j0.set_rb_delta(5,1,b); j0.set_rb_delta(6,1,g);
							TR << "perturb " << x << "," << y << "," << z << " ; " << a << "," << b << "," << g << std::endl;
							j0.fold_in_rb_deltas();
							pose.set_jump( sym_aware_jump_id, j0 );

							// call mover on perturbed pose
							mover_->apply( pose );
							ms = mover_->get_last_move_status();
							if ( ms == FAIL_RETRY ) {
								TR.Warning << "Mover failed. Retry. " << std::endl;
								continue;
							} else if ( ms == FAIL_DO_NOT_RETRY || ms == FAIL_BAD_INPUT ) {
								TR.Error << "Mover failed. Exit." << std::endl;
								break;
							}
							core::Real myscore;
							if ( usescore_ ) {
								myscore = (*scorefxn_)(pose);
							} else {
								myscore = filter_->report_sm( pose );
							}

							if ( maxscore_ )  myscore*=-1;
							if ( myscore < best_score ) {
								best_score = myscore;
								best_pose = pose;
							}
						}
					}
				}
			}
		}
	}

	pose = best_pose;
}

void
GenericSymmetricSampler::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &filters,
	Movers_map const &movers,
	Pose const & ) {
	dof_id_ = tag->getOption<std::string>("dof_id", "");
	angle_min_ = tag->getOption<Real>("angle_min", -1.0);
	angle_max_ = tag->getOption<Real>("angle_max",  1.0);
	angle_step_ = tag->getOption<Real>("angle_step", 1.0);
	radial_disp_min_ = tag->getOption<Real>("radial_disp_min", -1.0);
	radial_disp_max_ = tag->getOption<Real>("radial_disp_max", 1.0);
	radial_disp_step_ = tag->getOption<Real>("radial_disp_step", 1.0);

	maxscore_ = tag->getOption<bool>("select_max", false);

	// manditory: mover
	std::string mover_name( tag->getOption< std::string >( "mover_name"));
	Movers_map::const_iterator  find_mover ( movers.find( mover_name ));
	runtime_assert( find_mover!=movers.end() );
	mover_ = find_mover->second;

	// optional: filter
	std::string filter_name( tag->getOption< std::string >( "filter_name", "true_filter" ) );
	Filters_map::const_iterator find_filter( filters.find( filter_name ));
	runtime_assert( find_filter!=filters.end() );
	filter_ = find_filter->second;
	runtime_assert( filter_ != 0 );


	// optional scorefunction ... overrides filter
	std::string sfxn ( tag->getOption< std::string >( "scorefxn", "" ) );
	if ( sfxn != "" ) {
		scorefxn_ = data.get< core::scoring::ScoreFunction * >( "scorefxns", sfxn )->clone();   //fpd use clone
		TR << "Evaluation during is done by SCOREFUNCTION" << sfxn << std::endl;
		usescore_ = true;
	} else {
		TR << "Evaluation during is done by FILTER" << filter_name << std::endl;
		usescore_ = false;
	}
}

} // matdes
} // devel
