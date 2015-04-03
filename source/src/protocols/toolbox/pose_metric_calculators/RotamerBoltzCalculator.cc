// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/PoseMetricCalculators/RotamerBoltzCalculator.cc
/// @brief Calculates Rotamer occupancy of each rotameric state in a given set of residues.
/// @author Hetu Kamisetty
#include <protocols/toolbox/pose_metric_calculators/RotamerBoltzCalculator.hh>
#include <utility/stream_util.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/graph/Graph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMoverLazy.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>

#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <cmath>
#include <boost/foreach.hpp>

//Auto Headers
#include <protocols/simple_filters/ScoreTypeFilter.hh>
namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {
  RotamerBoltzCalculator::RotamerBoltzCalculator(core::scoring::ScoreFunctionOP scorefxn, core::Real temp, core::Real repack_rad){
    //pose_ = pose;
    scorefxn_ = scorefxn;
    //mm_ = new core::kinematics::MoveMap;
    temperature(temp);
    rb_jump_ = 1;
    repacking_radius(repack_rad);

  }

  utility::vector1<core::Real> RotamerBoltzCalculator::computeAllBoltz(core::pose::Pose& pose){
    utility::vector1<core::Real> allboltz;
    for(Size i=1;i<pose.total_residue(); i++){
      allboltz.push_back(computeBoltzWeight(pose, i));
    }
    all_boltz_ = allboltz;
    return allboltz;
  }
  core::Real RotamerBoltzCalculator::computeBoltzWeight(core::pose::Pose& pose, Size resi){
    core::pack::task::PackerTaskOP task = init_task(pose,resi);
    protocols::simple_moves::MinMoverOP  mm = init_minmover(pose, resi, true, task);
    return computeBoltzWeight(pose, resi, mm, task);//assume unbound
  }
  core::Real RotamerBoltzCalculator::computeBoltzWeight(core::pose::Pose& pose, Size resi,  protocols::simple_moves::MinMoverOP min_mover, core::pack::task::PackerTaskOP task){
    ///user doesn't have movemap, so can't initialize min_mover himself right now.

    //core::pose::Pose& pose = pose();
    //std::cout<<" in rotamer boltz calculator with resi "<<resi<<std::endl;
    using namespace protocols::moves;
    using namespace core::pack::rotamer_set;
    using namespace core::pack::task;
    using namespace core::conformation;
    protocols::simple_moves::PackRotamersMoverLazy pmover(scorefxn());
    task->set_bump_check(false);
    pmover.task(task);
    pmover.call_setup(pose);
    protocols::simple_filters::ScoreTypeFilter stf( scorefxn_, core::scoring::total_score, 0 );
    //std::cout<<" in rotamer boltz done with call_setup "<<resi<<std::endl;
    pmover.apply(pose);//what happens if setup and apply are called with different poses?
	  min_mover->apply( pose );
	  core::pose::Pose const const_min_pose( pose );
	  core::Real const init_score( stf.compute( const_min_pose ) );

    RotamerSetsCOP rotsets = pmover.rotamer_sets();
    Size moltenResid = rotsets->resid_2_moltenres(resi);
	  utility::vector1< core::Real > scores;
    utility::vector0<int> rot_to_pack;


    /*
    std::cout<<" num rotamers acc to rotsets "<<rotsets->nrotamers_for_moltenres(moltenResid)<<std::endl;
    std::cout<<" num rotamers acc to this "<<rotset_->num_rotamers()<<std::endl;
    RotamerSetCOP rset = rotsets->rotamer_set_for_residue(resi);
    for(Size i=1;i<=rotset_->num_rotamers(); i++){
      ResidueCOP r = rotset_->rotamer(i);
      //std::cout<<i<<"th rotamer acc to this "<<r->type()<<std::endl;
      utility::vector1<core::Real> chi = r->chi();
      for(Size j=1;j<=chi.size(); j++){
        std::cout<<j<<"th chi acc to "<<i<<"th rotamer is "<<chi.at(j)<<std::endl;
      }
    }
    for(Size i=1;i<=rset->num_rotamers(); i++){
      ResidueCOP r = rset->rotamer(i);
      utility::vector1<core::Real> chi = r->chi();
      for(Size j=1;j<=chi.size(); j++){
        std::cout<<j<<"th chi acc to "<<i<<"th rotamer is "<<chi.at(j)<<std::endl;
      }
    }
    */

    for(Size i=1;i<rotset_->num_rotamers(); i++){
      core::pose::Pose pose = const_min_pose;
    PROF_START(basic::TEST3);
      rot_to_pack = init_rot_to_pack(rotsets, moltenResid, i);
      pmover.apply_to_rotpack(pose, rot_to_pack);
		  min_mover->apply( pose );
    PROF_STOP(basic::TEST3);
		  core::Real const score( stf.compute( pose ) );
		  scores.push_back( score );
    }
    return computeBoltzSum(init_score, scores);
  //gives rot_to_pack where entry is rotameric id into rotamersets
  }
  /*RotamerBoltzCalculator::RotamerBoltzCalculator(Pose& pose, utility::vector0<int> moltres_to_pack, utility::vector0<int> rotid_in_moltres, core::Real temp){//gives rotamer id into rotamerset for each position of interest.
  }
  */

  core::pack::task::PackerTaskOP RotamerBoltzCalculator::init_task(core::pose::Pose& pose, core::Size resi){
	  using namespace core::pack::task;
	  using namespace core::pack::task::operation;
	  using namespace core::pack::rotamer_set;
	  using namespace core::conformation;

    Residue const & res = pose.residue( resi );
    RotamerSetFactory rsf;
    rotset_ = rsf.create_rotamer_set( res );
    rotset_->set_resid( resi );

    TaskFactoryOP tf( new core::pack::task::TaskFactory() );
    tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ));
    PackerTaskOP ptask( tf->create_task_and_apply_taskoperations( pose ) );
    ResidueLevelTask & restask( ptask->nonconst_residue_task( resi ) );
    restask.restrict_to_repacking();///hmk: what is this doing?
    core::graph::GraphOP packer_graph( new core::graph::Graph( pose.total_residue() ) );
    ptask->set_bump_check(false);
    rotset_->build_rotamers( pose, *scorefxn_, *ptask, packer_graph, false );

    /// des_around will mark all residues around resi for design, the rest for packing.
    protocols::toolbox::task_operations::DesignAroundOperationOP des_around( new protocols::toolbox::task_operations::DesignAroundOperation );
    des_around->design_shell( repacking_radius() );
    des_around->include_residue( resi );
    tf->push_back( des_around );
    PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );
    return task;
  }
  protocols::simple_moves::MinMoverOP RotamerBoltzCalculator::init_minmover(core::pose::Pose& pose, core::Size, bool unbound, core::pack::task::PackerTaskOP  task){
	  using namespace core::conformation;

    core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );

    mm->set_bb( false );
    if( unbound ) // the complex was split, don't minimize rb dof
      mm->set_jump( false );
    else // minimize rb if bound
      mm->set_jump( (int)(rb_jump()), true );
    for( core::Size i=1; i<=pose.total_residue(); ++i ){
      if( task->being_designed( i ) ){
        task->nonconst_residue_task( i ).restrict_to_repacking(); // mark all des around to repacking only
        mm->set_chi( i, true );
      }
      else{
        task->nonconst_residue_task( i ).prevent_repacking(); /// mark all non-desaround positions for no repacking
        //task->nonconst_residue_task( i ).restrict_to_repacking();
        mm->set_chi( i, false );
      }
    }
	  //task->nonconst_residue_task( resi ).prevent_repacking();

    protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover() );
    min_mover->score_function( scorefxn() );
    min_mover->movemap( mm );
    min_mover->min_options()->min_type( "dfpmin_armijo_nonmonotone" );
    return min_mover;
  }


  core::Real RotamerBoltzCalculator::computeBoltzSum(core::Real init_score, utility::vector1<core::Real> scores){
    core::Real boltz_sum ( 0.0 );
    BOOST_FOREACH( core::Real const score, scores )
      boltz_sum += exp(( init_score - score )/temperature());

    return( 1/boltz_sum );
  }

  void RotamerBoltzCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const{
    if ( key == "boltz" ) {
      basic::check_cast( valptr, &all_boltz_, "boltz expects to return a util:vector1<real>" );
		  (static_cast<basic::MetricValue<utility::vector1< core::Real > > *>(valptr))->set( all_boltz_ );

    } else {
      basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
      utility_exit();
    }

  }
  std::string RotamerBoltzCalculator::print( std::string const & key ) const{
    if ( key == "boltz" ) {
      return utility::to_string(all_boltz_);
    }
    else{
      return "error";
    }
  }
	void RotamerBoltzCalculator::recompute( core::pose::Pose const & this_pose ){

    core::pose::Pose cpose = this_pose;
    computeAllBoltz(cpose);
  }
  utility::vector0 <int> RotamerBoltzCalculator::init_rot_to_pack(core::pack::rotamer_set::RotamerSetsCOP rotamer_sets, core::Size moltenres, core::Size rot_to_fix){
	  utility::vector0 < int > rot_to_pack;
		rot_to_pack.reserve( rotamer_sets->nrotamers() );
    ///for ( Size rot = 1; rot <= num_rots_to_pack(); ++rot ) {
    for ( Size rot = 1; rot <= rotamer_sets->nrotamers(); ++rot ) {
      core::Size this_moltenres = rotamer_sets->moltenres_for_rotamer(rot);
      core::Size this_rotid = rotamer_sets->rotid_on_moltenresidue(rot);
      if(this_moltenres ==moltenres && !(this_rotid==rot_to_fix)){
        continue;
      }else{
        rot_to_pack.push_back( rot );
      }
		}
    return rot_to_pack;
  }


}//pose_metrics
}//toolbox
}//protocols
