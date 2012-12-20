// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CCDLoopCloser.cc
///
/// @brief
/// @author Tim Jacobs

//unit
#include <devel/loop_creation/CCDLoopCloser.hh>
#include <devel/loop_creation/CCDLoopCloserCreator.hh>

//core
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Edge.hh>

//protocols
#include <protocols/loops/loop_closure/ccd/ccd_closure.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/toolbox/pose_metric_calculators/ClashCountCalculator.hh>

//basic
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>

//numeric
#include <numeric/random/random.hh>

//utility
#include <utility/tag/Tag.hh>

namespace devel {
namespace loop_creation {

static basic::Tracer TR( "devel.loop_creation.CCDLoopCloser" );

//****CREATOR METHODS****//
std::string
CCDLoopCloserCreator::keyname() const
{
	return CCDLoopCloserCreator::mover_name();
}

protocols::moves::MoverOP
CCDLoopCloserCreator::create_mover() const {
	return new CCDLoopCloser;
}

std::string
CCDLoopCloserCreator::mover_name()
{
	return "CCDLoopCloser";
}
//****END CREATOR METHODS****//

	
CCDLoopCloser::CCDLoopCloser():
max_closure_attempts_(10),
max_ccd_moves_per_closure_attempt_(10000)
{
	init();
}

protocols::moves::MoverOP
CCDLoopCloser::clone() const {
	return( protocols::moves::MoverOP( new CCDLoopCloser( *this ) ) );
}
protocols::moves::MoverOP
CCDLoopCloser::fresh_instance() const {
	return protocols::moves::MoverOP( new CCDLoopCloser );
}
	
std::string
CCDLoopCloser::get_name() const {
	return "CCDLoopCloser";
}
	
void
CCDLoopCloser::init(){
	prevent_nonloop_modifications(true);
	clash_calculator_ = new protocols::toolbox::pose_metric_calculators::ClashCountCalculator(2.0);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "clash_calculator", clash_calculator_ );
}

void
CCDLoopCloser::apply(
	core::pose::Pose & pose
){

	using namespace core;

	core::kinematics::FoldTree saved_ft = pose.fold_tree();
	if(prevent_nonloop_modifications())
	{
		prepare_fold_tree(pose);
	}
	
	
	basic::MetricValue<core::Size> bb_clash_metric;
	clash_calculator_->get("bb", bb_clash_metric, pose);
	core::Size starting_clashes=bb_clash_metric.value();
	TR << "Number of BB clashes in input pose: " << starting_clashes << std::endl;
		
	protocols::loops::add_single_cutpoint_variant(pose, loop());
	
	// setup movemap
	kinematics::MoveMap mm;
	for ( Size ii=loop().start(); ii<=loop().stop(); ++ii )
	{
		mm.set_bb( ii, true );
		
		//don't change phi for prolines
		if ( pose.residue(ii).aa() == chemical::aa_pro )
		{
			mm.set( id::TorsionID( id::phi_torsion, id::BB, ii ), false );
		}
	}
	TR.Debug << "CCD residues: " << loop().start() << " " << loop().stop() << std::endl;
	TR.Debug << "CCD cutpoint: " << loop().cut() << std::endl;
	
	//MAKE THESE PARAMETERS
	core::Real max_rama_score_increase_( 2.0 ),
		max_total_delta_helix_( 15 ),
		max_total_delta_strand_( 15 ),
		max_total_delta_loop( 15 ),
		tolerance_( 0.01 );
	
	//output variable for fast_ccd_closure
	core::Real forward_deviation_, backward_deviation_, torsion_delta_, rama_delta_;
	core::pose::Pose saved_pose = pose;
	for(core::Size i=1; i<=max_closure_attempts_; ++i)
	{
		core::Size num_moves_needed = protocols::loops::loop_closure::ccd::fast_ccd_loop_closure(
			pose,
			mm,
			loop().start(),
			loop().stop(),
			loop().cut(),
			max_ccd_moves_per_closure_attempt_,
			tolerance_,
			true/*rama_check*/,
			max_rama_score_increase_,
			max_total_delta_helix_,
			max_total_delta_strand_,
			max_total_delta_loop,
			forward_deviation_,
			backward_deviation_,
			torsion_delta_,
			rama_delta_
		);
		
		TR.Debug << "Loop closure attempt " << i << " returned in " << num_moves_needed << " iterations (max of "
			<< max_ccd_moves_per_closure_attempt_ << ")" << std::endl;
			
		TR.Debug << "torsion rmsd: " << torsion_delta_ << std::endl;
		TR.Debug << "rama delta: " << rama_delta_ << std::endl;
		TR.Debug << "forward deviation: " << forward_deviation_ << std::endl;
		TR.Debug << "backward deviation: " << backward_deviation_ << std::endl;
	
		if(forward_deviation_ < 0.1 || backward_deviation_ < 0.1)
		{
			//check for clashes
			basic::MetricValue<core::Size> loop_clash_metric;
			clash_calculator_->get("bb", loop_clash_metric, pose);
			core::Size loop_clashes = loop_clash_metric.value() - starting_clashes;
			TR << "Number of BB clashes from loop: " << loop_clashes << std::endl;
			if(loop_clashes == 0)
			{
				pose.fold_tree(saved_ft);
				TR << "Loop " << loop().start() << " " << loop().stop() << " successfully closed in " << i
					<< " attempts (max of " << max_closure_attempts_ << ")" << std::endl;
				success_=true;
				break;
			}
		}
		//restore pose to original and try closure again
		pose=saved_pose;
	}
	pose.fold_tree(saved_ft);
}

///@brief Create a fold tree that prevents downstream propogation
/// of loop insertions
void
CCDLoopCloser::prepare_fold_tree(
	core::pose::Pose & pose
){
	using namespace core;

	//prepare special foldtree
	core::kinematics::FoldTree new_ft;
	new_ft.add_edge(1, loop().start()-1, kinematics::Edge::PEPTIDE);
	new_ft.add_edge(loop().start()-1, loop().cut(), kinematics::Edge::PEPTIDE);
	new_ft.add_edge(loop().start()-1, loop().stop()+1, 1);
	new_ft.add_edge(loop().stop()+1, loop().cut()+1, kinematics::Edge::PEPTIDE);
	new_ft.add_edge(loop().stop()+1, pose.total_residue(), kinematics::Edge::PEPTIDE);
			
	pose.fold_tree(new_ft);
}
	
///@brief parse tag for use in RosettaScripts
void
CCDLoopCloser::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){
	prevent_nonloop_modifications(tag->getOption< bool >("prevent_nonloop_modifications", true));
	
	max_ccd_moves_per_closure_attempt_ = tag->getOption< core::Size >("max_ccd_moves_per_closure_attempt", 10000);
	max_closure_attempts_ = tag->getOption< core::Size >("max_closure_attempts", 10);
}

} //loop creation
} //devel
