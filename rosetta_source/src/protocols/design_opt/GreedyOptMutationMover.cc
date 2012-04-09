// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Chris King (chrisk1@uw.edu)
//#include <algorithm >
#include <protocols/design_opt/GreedyOptMutationMover.hh>
#include <protocols/design_opt/GreedyOptMutationMoverCreator.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/pose/PDBInfo.hh>
#include <fstream>
// AUTO-REMOVED #include <utility/file/FileName.hh>
#include <iostream>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

//Auto Headers
#include <basic/options/keys/OptionKeys.hh>

namespace protocols {
namespace design_opt {

static basic::Tracer TR( "protocols.design_opt.GreedyOptMutationMover" );
using namespace core;
using namespace chemical;

///@brief default ctor
GreedyOptMutationMover::GreedyOptMutationMover() :
	Mover( GreedyOptMutationMoverCreator::mover_name() ),
	task_factory_( NULL ),
	filter_( NULL ),
	relax_mover_( NULL ),
	scorefxn_( NULL ),
	sample_type_( "low" ),
	dump_pdb_( false ),
	stopping_condition_( NULL )
{
	if( sample_type_ == "high" ){
		flip_sign_ = Real( -1 );
	}else if( sample_type_ == "low" ){
		flip_sign_ = Real( 1 );
	}else{
		TR << "WARNING: the sample type, " << sample_type_ << ", is not defined. Use \'high\' or \'low\'." << std::endl;
		runtime_assert( false );
  }
}
//TODO: non-default ctor

protocols::moves::MoverOP
GreedyOptMutationMover::relax_mover() const{
	return relax_mover_;
}

void
GreedyOptMutationMover::relax_mover( protocols::moves::MoverOP mover ){
	relax_mover_ = mover;
}

protocols::filters::FilterOP GreedyOptMutationMover::filter() const{
	return filter_;
}

void
GreedyOptMutationMover::filter( protocols::filters::FilterOP filter ){
	filter_ = filter;
}

core::pack::task::TaskFactoryOP
GreedyOptMutationMover::task_factory() const
{
	return task_factory_;
}

void
GreedyOptMutationMover::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

bool
GreedyOptMutationMover::report_all() const{
  return report_all_;
}

void
GreedyOptMutationMover::report_all( bool const report_all ){
  report_all_ = report_all;
}

void
GreedyOptMutationMover::dump_pdb( bool const dump_pdb ){
  dump_pdb_ = dump_pdb;
}

bool
GreedyOptMutationMover::dump_pdb() const{
  return dump_pdb_;
}

void
GreedyOptMutationMover::sample_type( std::string const sample_type ){
  sample_type_ = sample_type;
}

std::string
GreedyOptMutationMover::sample_type() const{
  return sample_type_;
}

//for comparing values in sort
bool
cmp_pair_by_second(
	std::pair< AA, Real > const pair1,
	std::pair< AA, Real > const pair2 )
{
	return pair1.second < pair2.second;
}

//for comparing values in sort
bool
cmp_triple_by_third(
	std::pair< Size, std::pair< AA, Real > > const pair1,
	std::pair< Size, std::pair< AA, Real > > const pair2 )
{
	return pair1.second.second < pair2.second.second;
}

core::scoring::ScoreFunctionOP
GreedyOptMutationMover::scorefxn() const{
	return scorefxn_;
}

void
GreedyOptMutationMover::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

GreedyOptMutationMover::~GreedyOptMutationMover(){}

protocols::moves::MoverOP
GreedyOptMutationMover::clone() const{
	return new GreedyOptMutationMover( *this );
}

std::string
GreedyOptMutationMoverCreator::keyname() const
{
	return GreedyOptMutationMoverCreator::mover_name();
}

protocols::moves::MoverOP
GreedyOptMutationMoverCreator::create_mover() const {
	return new GreedyOptMutationMover;
}

std::string
GreedyOptMutationMoverCreator::mover_name()
{
	return "GreedyOptMutationMover";
}

std::string
GreedyOptMutationMover::get_name() const {
  return GreedyOptMutationMoverCreator::mover_name();
}

void
GreedyOptMutationMover::apply(core::pose::Pose & pose )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	//store input pose
	core::pose::Pose pose_orig( pose );

	//apply input task factory to pose
	PackerTaskCOP task = task_factory()->create_task_and_apply_taskoperations( pose_orig );

	//get vector< Size > of protein seqpos being designed
	utility::vector1< core::Size > being_designed;
	being_designed.clear();
	for( core::Size resi = 1; resi <= pose_orig.total_residue(); ++resi ){
		if( task->residue_task( resi ).being_designed() && pose_orig.residue(resi).is_protein() )
			being_designed.push_back( resi );
	}
	if( being_designed.empty() ) {
		TR.Warning << "WARNING: No residues are listed as designable." << std::endl;
	}

	//create map of seqpos -> vector of AA/val pairs that pass input filter
	std::map< core::Size, utility::vector1< std::pair< AA, Real > > > residue_id_map;
	//create map of seqpos/AA pair -> input filter score/pass pair
	std::map< std::pair< core::Size, AA >, std::pair< core::Real, bool > > residue_id_val_map;
	residue_id_map.clear(); residue_id_val_map.clear();

	//for each seqpos in being_designed vector
	foreach( core::Size const resi, being_designed ){
		//create vector< AA > of allowed residue types at seqpos
		typedef std::list< ResidueTypeCAP > ResidueTypeCAPList;
		ResidueTypeCAPList const & allowed( task->residue_task( resi ).allowed_residue_types() );
		utility::vector1< AA > allow_temp;
		allow_temp.clear();
		foreach( ResidueTypeCAP const t, allowed ){
			if(std::find(allow_temp.begin(),allow_temp.end(),t->aa())!=allow_temp.end()) continue;
			allow_temp.push_back( t->aa() );
		}
		//for each allowed AA
		foreach( AA const target_aa, allow_temp ){
			//dont need to test if target_aa is same as original aa
			if( target_aa == pose_orig.residue( resi ).type().aa() ) continue;
			//reset pose to original
			pose = pose_orig;
			//make bool vector of allowed aa's [1,20], all false except for target_aa
			utility::vector1< bool > allowed_aas;
			allowed_aas.clear();
			allowed_aas.assign( num_canonical_aas, false );
			allowed_aas[ target_aa ] = true;
			//make mut_res task factory, mutates resi to target_aa  repack 8A shell
			core::pack::task::TaskFactoryOP mut_res = new core::pack::task::TaskFactory( *task_factory() );
			protocols::toolbox::task_operations::DesignAroundOperationOP repack_around_op =
					new protocols::toolbox::task_operations::DesignAroundOperation;
			repack_around_op->design_shell( -1.0 ); //neg radius insures no designing nbrs
			repack_around_op->repack_shell( 8.0 );
			repack_around_op->allow_design( true ); //because we still want to design resi
			repack_around_op->include_residue( resi );
			mut_res->push_back( repack_around_op );
			//make mutate_residue packer task, mutate target res, repack all others
			PackerTaskOP mutate_residue = mut_res->create_task_and_apply_taskoperations( pose );
			mutate_residue->initialize_from_command_line().or_include_current( true );
			mutate_residue->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
/*
			for( core::Size resj = 1; resj <= pose.total_residue(); ++resj ){
				//TODO: retrict to nbr residues?
				if( resi != resj )
					mutate_residue->nonconst_residue_task( resj ).restrict_to_repacking();
				else
					mutate_residue->nonconst_residue_task( resj ).restrict_absent_canonical_aas( allowed_aas );
			}
*/
			TR<<"Mutating residue "<<pose.residue( resi ).name3()<<resi<<" to ";
			//run PackRotamers with mutate_residue task
			protocols::simple_moves::PackRotamersMoverOP pack;
			if( core::pose::symmetry::is_symmetric( pose ) )
				pack =  new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn(), mutate_residue );
			else
				pack = new protocols::simple_moves::PackRotamersMover( scorefxn(), mutate_residue );
			pack->apply( pose );
			TR<<pose.residue( resi ).name3()<<". Now relaxing..."<<std::endl;
			//then run input relax mover
			relax_mover()->apply( pose );
			//then check if passes input filter, bail out if it doesn't
			//store filter value and pass/fail in residue_is_val_map, store aa/val pair in residue_id_map
			//TODO: if no filter defined, just use total_score
			bool const filter_pass( filter()->apply( pose ) );
			//val sign is switched if type is high
			Real const val( flip_sign_ * filter()->report_sm( pose ) );
			if( !filter_pass ){
				TR<<"Filter fails with value "<< val << std::endl;
				if(report_all_) {
					residue_id_val_map[ std::pair< core::Size, AA >( resi, target_aa ) ] = std::pair< core::Real, bool >( val, false);
				}
				continue;
			}
			TR<<"Filter succeeds with value "<< val << std::endl;
			residue_id_val_map[ std::pair< core::Size, AA >( resi, target_aa ) ] =
					std::pair< core::Real, bool >( val, true);
			residue_id_map[ resi ].push_back( std::pair< AA, Real >( target_aa, val ) );
			//dump pdb?
			if( dump_pdb() ){
				using namespace protocols::jd2;
				JobOP job( JobDistributor::get_instance()->current_job() );
				std::stringstream fname;
				fname << job->input_tag() << pose_orig.residue( resi ).name3() << resi << pose.residue( resi ).name3()<<".pdb";
				TR<<"Saving pose "<<fname.str() << std::endl;
				pose.dump_scored_pdb( fname.str(), *scorefxn() );
			}
			TR.flush();
		}//foreach target_aa
	}//foreach resi

	//iterate over residue_id_map
	//some elements of residue_id_map may be empty vectors because didn't pass filter
	//sort each seqpos AA/val vector by value
	//append best scoring restype, val to seqpos_resid_vals
	utility::vector1< std::pair< Size, std::pair< AA, Real > > > seqpos_resid_vals;
	for( std::map< core::Size, utility::vector1< std::pair< AA, Real > > >::iterator resid_vals = residue_id_map.begin();
			resid_vals != residue_id_map.end(); ++resid_vals ){
		//skip if resid_vals vector is empty
		if( resid_vals->second.empty() ) continue;
		//sort resid_vals in incr val order
		std::sort( resid_vals->second.begin(), resid_vals->second.end(), cmp_pair_by_second );
		//best val is lowest, store it in seqpos_resid_vals
		//TODO: option to include best N restypes?
		std::pair< AA, Real > best_resid_val( resid_vals->second[ 1 ] );
		std::pair< Size, std::pair< AA, Real > > seqpos_resid_val( resid_vals->first, best_resid_val );
		seqpos_resid_vals.push_back( seqpos_resid_val );
	}

	//now sort seqpos_resid_vals by value, low to high
	std::sort( seqpos_resid_vals.begin(), seqpos_resid_vals.end(), cmp_triple_by_third );

	TR<<"Combining sorted independently optimal muations… " << std::endl;
	//reset pose to original, init starting filter val
	//must use same relax mover before scoring so comparison is fair!
	pose = pose_orig;
	relax_mover()->apply( pose );
	Real best_val( flip_sign_ * filter()->report_sm( pose ) );
	pose::Pose best_pose( pose );
	//now try each AA at each position
	for( Size iseq = 1; iseq <= seqpos_resid_vals.size(); ++iseq ){
		std::pair< Size, std::pair< AA, Real > > seqpos_resid_val( seqpos_resid_vals[ iseq ] );
		pose = best_pose;
		Size resi( seqpos_resid_val.first );
		AA target_aa( seqpos_resid_val.second.first );

		//make bool vector of allowed aa's [1,20], all false except for target_aa
		utility::vector1< bool > allowed_aas;
		allowed_aas.clear();
		allowed_aas.assign( num_canonical_aas, false );
		allowed_aas[ target_aa ] = true;
		//make mut_res task factory, mutates resi to target_aa
		core::pack::task::TaskFactoryOP mut_res = new core::pack::task::TaskFactory();
		//make mutate_residue packer task, mutate target res, repack all others
		protocols::toolbox::task_operations::DesignAroundOperationOP repack_around_op =
			new protocols::toolbox::task_operations::DesignAroundOperation;
		repack_around_op->design_shell( -1.0 ); //neg radius insures no designing nbrs
		repack_around_op->repack_shell( 8.0 );
		repack_around_op->allow_design( true ); //because we still want to design resi
		repack_around_op->include_residue( resi );
		mut_res->push_back( repack_around_op );
		PackerTaskOP mutate_residue = mut_res->create_task_and_apply_taskoperations( pose );
		mutate_residue->initialize_from_command_line().or_include_current( true );
		mutate_residue->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
	/*
		for( core::Size resj = 1; resj <= pose.total_residue(); ++resj ){
			//TODO: retrict to nbr residues?
			if( resi != resj )
				mutate_residue->nonconst_residue_task( resj ).restrict_to_repacking();
			else
				mutate_residue->nonconst_residue_task( resj ).restrict_absent_canonical_aas( allowed_aas );
		}
	*/
		TR<<"Mutating residue "<<pose.residue( resi ).name3()<<resi<<" to ";
		//run PackRotamers with mutate_residue task
		protocols::simple_moves::PackRotamersMoverOP pack;
		if( core::pose::symmetry::is_symmetric( pose ) )
			pack =  new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn(), mutate_residue );
		else
			pack = new protocols::simple_moves::PackRotamersMover( scorefxn(), mutate_residue );
		pack->apply( pose );
		TR<<pose.residue( resi ).name3()<<". Now relaxing..."<<std::endl;
		//then run input relax mover
		relax_mover()->apply( pose );
		//then check if passes input filter, bail out if it doesn't
		//TODO: if no filter defined, just use total_score
		bool const filter_pass( filter()->apply( pose ) );
		Real const val( flip_sign_ * filter()->report_sm( pose ) );
		if( !filter_pass ){
			TR<<"Filter fails with value "<< val << std::endl;
			continue;
		}
		TR<<"Filter succeeds with value "<< val << std::endl;
		//score mutation, reset best_pose, best val if is lower
		if( val > best_val ){
			TR<<"Mutation rejected. Current best value is "<< best_val << std::endl;
			continue;
		}
		TR<<"Mutation accepted. New best value is "<< val << std::endl;
		best_val = val;
		best_pose = pose;
		if( stopping_condition() && stopping_condition()->apply( pose ) ){
			TR<<"Stopping condition evaluates to true. Stopping early."<<std::endl;
			return;
		}
	}

//reset to best pose after last step
pose = best_pose;

}

//parse rosetta scripts tags
void
GreedyOptMutationMover::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & )
{
	TR << "GreedyOptMutationMover"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	//load filter
	std::string const filter_name( tag->getOption< std::string >( "filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator filter_it( filters.find( filter_name ) );
	if( filter_it == filters.end() )
		utility_exit_with_message( "Filter "+filter_name+" not found" );
	filter( filter_it->second );
	//load relax mover
	std::string const relax_mover_name( tag->getOption< std::string >( "relax_mover", "null" ) );
	protocols::moves::Movers_map::const_iterator mover_it( movers.find( relax_mover_name ) );
	if( mover_it == movers.end() )
		utility_exit_with_message( "Relax mover "+relax_mover_name+" not found" );
	relax_mover( mover_it->second );
	//load sample_type
	sample_type( tag->getOption< std::string >( "sample_type", "low" ) );
	//load report_all
	report_all( tag->getOption< bool >( "report_all", false ) );
	//load scorefxn
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	//load dump_pdb
	dump_pdb( tag->getOption< bool >( "dump_pdb", false ) );
	if( tag->hasOption( "stopping_condition" ) ){
		std::string const stopping_filter_name( tag->getOption< std::string >( "stopping_condition" ) );
		stopping_condition( protocols::rosetta_scripts::parse_filter( stopping_filter_name, filters ) );
		TR<<"Defined stopping condition "<<stopping_filter_name<<std::endl;
	}
}


} // moves
} // protocols
