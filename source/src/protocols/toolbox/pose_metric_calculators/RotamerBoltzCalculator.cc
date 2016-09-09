// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/PoseMetricCalculators/RotamerBoltzCalculator.cc
/// @brief Calculates Rotamer occupancy of each rotameric state in a given set of residues.
/// @author Hetu Kamisetty
/// @author Tom Linsky (tlinsky at uw dot edu)

// Unit Headers
#include <protocols/toolbox/pose_metric_calculators/RotamerBoltzCalculator.hh>

// Protocol Headers
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/toolbox/EnergyLandscapeEvaluator.hh>

// Core Headers
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/rms_util.tmpl.hh>

// Basic/Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/stream_util.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <utility/graph/Graph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/types.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMoverLazy.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>

#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <cmath>
#include <boost/foreach.hpp>

//Auto Headers
#include <protocols/simple_filters/ScoreTypeFilter.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
#endif // SERIALIZATION

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.pose_metric_calculators.RotamerBoltzCalculator" );

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

RotamerBoltzCalculator::RotamerBoltzCalculator(
	core::scoring::ScoreFunctionOP const scorefxn,
	core::Real const temp,
	core::Real const repack_rad ):
	core::pose::metrics::StructureDependentCalculator(),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	evaluator_( new protocols::toolbox::RotamerBoltzmannWeightEvaluator( temp, true ) ),
	lazy_( false ),
	rb_jump_( 1 ),
	repacking_radius_( repack_rad ),
	scorefxn_( scorefxn ),
	temperature_( temp ),
	probabilities_(),
	rotset_()
{}

core::pose::metrics::PoseMetricCalculatorOP
RotamerBoltzCalculator::clone() const
{
	return core::pose::metrics::PoseMetricCalculatorOP( new RotamerBoltzCalculator(*this) );
}

void
RotamerBoltzCalculator::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
	notify_structure_change();
}

void
RotamerBoltzCalculator::set_energy_landscape_evaluator( EnergyLandscapeEvaluatorCOP evaluator )
{
	evaluator_ = evaluator;
	notify_structure_change();
}

void
RotamerBoltzCalculator::set_lazy( bool const lazy )
{
	lazy_ = lazy;
	notify_structure_change();
}

void
RotamerBoltzCalculator::set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn )
{
	scorefxn_ = scorefxn->clone();
}

void
RotamerBoltzCalculator::computeAllBoltz( core::pose::Pose const & pose )
{
	if ( !selector_ ) utility_exit_with_message( "RotamerBoltzCalculator: residue selector not set!" );

	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );
	debug_assert( subset.size() == pose.size() );

	probabilities_.clear();
	for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
		if ( !subset[ resid ] ) continue;
		TR << "Computing Boltz Weight for residue " << resid << std::endl;
		core::pose::Pose work_pose = pose;
		probabilities_[ resid ] = computeBoltzWeight( work_pose, resid );
	}

	all_boltz_.clear();
	all_boltz_.assign( pose.size(), 0.0 );
	for ( RotamerProbabilities::const_iterator rp=probabilities_.begin(); rp!=probabilities_.end(); ++rp ) {
		all_boltz_[ rp->first ] = rp->second;
	}
}

core::Real RotamerBoltzCalculator::computeBoltzWeight(core::pose::Pose& pose, Size resi){
	core::pack::task::PackerTaskOP task = init_task(pose,resi);
	protocols::simple_moves::MinMoverOP  mm = init_minmover(pose, resi, true, task);
	return computeBoltzWeight(pose, resi, mm, task);//assume unbound
}

core::Real
RotamerBoltzCalculator::computeBoltzWeight(
	core::pose::Pose & pose,
	core::Size resi,
	protocols::simple_moves::MinMoverOP min_mover,
	core::pack::task::PackerTaskOP task )
{
	if ( lazy_ ) {
		return computeBoltzWeight_lazy( pose, resi, min_mover, task );
	} else {
		return compute_boltz_weight_packrotamers( pose, resi, min_mover, task );
	}
}

core::Real
compute_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2, ObjexxFCL::FArray1D< bool > const & subset )
{
	return core::scoring::rmsd_no_super_subset( pose1, pose2, subset, core::scoring::is_scatom );
}

core::Real
RotamerBoltzCalculator::compute_boltz_weight_packrotamers(
	core::pose::Pose & pose,
	core::Size const resi,
	protocols::simple_moves::MinMoverOP min_mover,
	core::pack::task::PackerTaskOP task ) const
{
	// filter used for scoring
	protocols::simple_filters::ScoreTypeFilter stf( scorefxn_, core::scoring::total_score, 0 );

	// freeze input residue
	task->nonconst_residue_task( resi ).prevent_repacking();

	// repack/min the input pose
	core::pack::pack_rotamers( pose, *scorefxn_, task );
	min_mover->apply( pose );
	core::pose::Pose const const_min_pose( pose );

	// build set of test rotamers
	core::pack::rotamer_set::RotamerSetCOP rotset = create_rotamer_set( pose, resi );

	// for RMSD "subset"
	ObjexxFCL::FArray1D< bool > subset( pose.size(), false );
	subset( resi ) = true;

	ScoreRmsPoints scores( ScoreRmsPoint( stf.compute( pose ), 0.0 ) );
	TR << "Initial input pose has score " << scores.bg().score() << std::endl;

	for ( core::pack::rotamer_set::Rotamers::const_iterator rotamer=rotset->begin(); rotamer!=rotset->end(); ++rotamer ) {
		pose = const_min_pose;
		pose.replace_residue( resi, **rotamer, false/*orient bb*/ );
		core::pack::pack_rotamers( pose, *scorefxn_, task );
		min_mover->apply( pose );
		core::Real const score = stf.compute( pose );
		core::Real rmsd;
		if ( pose.residue(resi).name1() == 'G' ) {
			rmsd = 0.0;
		} else {
			rmsd = compute_rmsd( const_min_pose, pose, subset );
		}
		TR << "This rotamer has score " << score << " and sidechain atom RMSD " << rmsd << std::endl;
		scores.push_back( ScoreRmsPoint( score, rmsd ) );
	}

	return evaluator_->compute( scores );
}


core::Real RotamerBoltzCalculator::computeBoltzWeight_lazy(core::pose::Pose& pose, Size resi,  protocols::simple_moves::MinMoverOP min_mover, core::pack::task::PackerTaskOP task){
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

	for ( Size i=1; i<rotset_->num_rotamers(); i++ ) {
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

core::pack::task::PackerTaskOP
RotamerBoltzCalculator::init_task( core::pose::Pose const & pose, core::Size const resi )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::pack::rotamer_set;
	using namespace core::conformation;

	// build rotamer set for current residue
	Residue const & res = pose.residue( resi );
	RotamerSetFactory rsf;
	rotset_ = rsf.create_rotamer_set( res );
	rotset_->set_resid( resi );

	TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ));
	PackerTaskOP ptask( tf->create_task_and_apply_taskoperations( pose ) );
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pack::make_symmetric_PackerTask_by_truncation( pose, ptask );
	}
	ResidueLevelTask & restask( ptask->nonconst_residue_task( resi ) );
	restask.restrict_to_repacking();///hmk: what is this doing?
	utility::graph::GraphOP packer_graph( new utility::graph::Graph( pose.size() ) );
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
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );
	}

	mm->set_bb( false );
	if ( unbound ) { // the complex was split, don't minimize rb dof
		mm->set_jump( false );
	} else { // minimize rb if bound
		mm->set_jump( (int)(rb_jump()), true );
	}
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( task->being_designed( i ) ) {
			task->nonconst_residue_task( i ).restrict_to_repacking(); // mark all des around to repacking only
			mm->set_chi( i, true );
		} else {
			task->nonconst_residue_task( i ).prevent_repacking(); /// mark all non-desaround positions for no repacking
			//task->nonconst_residue_task( i ).restrict_to_repacking();
			mm->set_chi( i, false );
		}
	}

	protocols::simple_moves::MinMoverOP min_mover;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		min_mover = protocols::simple_moves::MinMoverOP( new protocols::simple_moves::symmetry::SymMinMover(
			mm, scorefxn(), "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );
	} else {
		min_mover = protocols::simple_moves::MinMoverOP( new protocols::simple_moves::MinMover(
			mm, scorefxn(), "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );
	}
	return min_mover;
}

core::pack::rotamer_set::RotamerSetOP
RotamerBoltzCalculator::create_rotamer_set( core::pose::Pose const & pose, core::Size const resi ) const
{
	core::conformation::Residue const & res = pose.residue( resi );
	core::pack::rotamer_set::RotamerSetFactory rsf;
	core::pack::rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( res );
	rotset->set_resid( resi );

	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	core::pack::task::PackerTaskOP ptask( tf->create_task_and_apply_taskoperations( pose ) );
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pack::make_symmetric_PackerTask_by_truncation( pose, ptask );
	}
	ptask->or_include_current( false );
	ptask->nonconst_residue_task( resi ).restrict_to_repacking();
	ptask->set_bump_check( false );
	utility::graph::GraphOP packer_graph( new utility::graph::Graph( pose.size() ) );
	rotset->build_rotamers( pose, *scorefxn_, *ptask, packer_graph, false );
	return rotset;
}

core::Real RotamerBoltzCalculator::computeBoltzSum(core::Real init_score, utility::vector1<core::Real> scores){
	core::Real boltz_sum ( 0.0 );
	BOOST_FOREACH ( core::Real const score, scores ) {
		boltz_sum += exp(( init_score - score )/temperature());
	}

	return( 1/boltz_sum );
}

void
RotamerBoltzCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const
{
	if ( key == "boltz" ) {
		basic::check_cast( valptr, &all_boltz_, "boltz expects to return a util:vector1<real>" );
		(static_cast<basic::MetricValue<utility::vector1< core::Real > > *>(valptr))->set( all_boltz_ );
	} else if ( key == "probabilities" ) {
		basic::check_cast( valptr, &probabilities_, "RotamerBoltzCalculator::lookup(): probabilties expects to return a RotamerProbabilities object" );
		static_cast< basic::MetricValue< RotamerProbabilities > *>( valptr )->set( probabilities_ );
	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}
}

std::string
RotamerBoltzCalculator::print( std::string const & key ) const
{
	if ( key == "boltz" ) {
		return utility::to_string(all_boltz_);
	} else if ( key == "probabilities" ) {
		std::stringstream ss;
		for ( RotamerProbabilities::const_iterator rp=probabilities_.begin(); rp!=probabilities_.end(); ++rp ) {
			ss << rp->first << " " << rp->second << std::endl;
		}
		return ss.str();
	} else {
		return "error";
	}
}

void RotamerBoltzCalculator::recompute( core::pose::Pose const & this_pose )
{
	computeAllBoltz( this_pose );
}

utility::vector0 <int> RotamerBoltzCalculator::init_rot_to_pack(core::pack::rotamer_set::RotamerSetsCOP rotamer_sets, core::Size moltenres, core::Size rot_to_fix){
	utility::vector0 < int > rot_to_pack;
	rot_to_pack.reserve( rotamer_sets->nrotamers() );
	///for ( Size rot = 1; rot <= num_rots_to_pack(); ++rot ) {
	for ( Size rot = 1; rot <= rotamer_sets->nrotamers(); ++rot ) {
		core::Size this_moltenres = rotamer_sets->moltenres_for_rotamer(rot);
		core::Size this_rotid = rotamer_sets->rotid_on_moltenresidue(rot);
		if ( this_moltenres ==moltenres && !(this_rotid==rot_to_fix) ) {
			continue;
		} else {
			rot_to_pack.push_back( rot );
		}
	}
	return rot_to_pack;
}


}//pose_metrics
}//toolbox
}//protocols

#ifdef    SERIALIZATION
/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculator::RotamerBoltzCalculator() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( selector_ ) ); // ResidueSelectorCOP
	arc( CEREAL_NVP( evaluator_ ) ); // EnergyLandscapeEvaluatorCOP
	arc( CEREAL_NVP( lazy_ ) ); // bool
	arc( CEREAL_NVP( rb_jump_ ) ); // core::Real
	arc( CEREAL_NVP( repacking_radius_ ) ); // core::Real
	arc( CEREAL_NVP( scorefxn_ ) ); // core::scoring::ScoreFunctionOP
	arc( CEREAL_NVP( temperature_ ) ); // core::Real
	arc( CEREAL_NVP( probabilities_ ) ); // std::map< core::Size, core::Real >
	arc( CEREAL_NVP( all_boltz_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( rotset_ ) ); // core::pack::rotamer_set::RotamerSetOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculator::load( Archive & arc ) {
	arc( selector_ );
	arc( evaluator_ );
	arc( lazy_ );
	arc( rb_jump_ ); // core::Real
	arc( repacking_radius_ ); // core::Real
	arc( scorefxn_ ); // core::scoring::ScoreFunctionOP
	arc( temperature_ ); // core::Real
	arc( probabilities_ ); // std::map< core::Size, core::Real >
	arc( all_boltz_ ); // utility::vector1<core::Real>
	arc( rotset_ ); // core::pack::rotamer_set::RotamerSetOP
}
SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculator );
CEREAL_REGISTER_TYPE( protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_RotamerBoltzCalculator )
#endif // SERIALIZATION
