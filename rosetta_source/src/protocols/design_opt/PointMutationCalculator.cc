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
#include <protocols/design_opt/PointMutationCalculator.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/pose/PDBInfo.hh>
#include <fstream>
// AUTO-REMOVED #include <utility/file/FileName.hh>
#include <iostream>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <numeric/random/random.hh>
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
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <core/pose/symmetry/util.hh>

//Auto Headers
#include <basic/options/keys/OptionKeys.hh>

namespace protocols {
namespace design_opt {

static basic::Tracer TR( "protocols.design_opt.PointMutationCalculator" );
static numeric::random::RandomGenerator RG( 54 );
using namespace core;
using namespace chemical;
using utility::vector1;
using std::pair;

///@brief default ctor
PointMutationCalculator::PointMutationCalculator() :
	task_factory_( NULL ),
//	filters_( NULL ), /*TODO: this throws a warning!*/
	relax_mover_( NULL ),
	scorefxn_( NULL ),
	dump_pdb_( false ),
	sample_type_( "low" )
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

//full ctor
PointMutationCalculator::PointMutationCalculator(
	core::pack::task::TaskFactoryOP task_factory,
	core::scoring::ScoreFunctionOP scorefxn,
	vector1< protocols::filters::FilterOP > filters,
	protocols::moves::MoverOP relax_mover,
	bool dump_pdb,
	std::string sample_type
)
{
	task_factory_ = task_factory;
	filters_ = filters;
	relax_mover_ = relax_mover;
	scorefxn_ = scorefxn;
	dump_pdb_ = dump_pdb;
	sample_type_ = sample_type;
	
	if( sample_type_ == "high" ){
		flip_sign_ = Real( -1 );
	}else if( sample_type_ == "low" ){
		flip_sign_ = Real( 1 );
	}else{
		TR << "WARNING: the sample type, " << sample_type_ << ", is not defined. Use \'high\' or \'low\'." << std::endl;
		runtime_assert( false );
  }
}

//backcompatible; single-filter convenience ctor
PointMutationCalculator::PointMutationCalculator(
	core::pack::task::TaskFactoryOP task_factory,
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::filters::FilterOP filter,
	protocols::moves::MoverOP relax_mover,
	bool dump_pdb,
	std::string sample_type
)
{
	vector1< protocols::filters::FilterOP > filters;
	filters.push_back( filter );
	task_factory_ = task_factory;
	filters_ = filters;
	relax_mover_ = relax_mover;
	scorefxn_ = scorefxn;
	dump_pdb_ = dump_pdb;
	sample_type_ = sample_type;
	
	if( sample_type_ == "high" ){
		flip_sign_ = Real( -1 );
	}else if( sample_type_ == "low" ){
		flip_sign_ = Real( 1 );
	}else{
		TR << "WARNING: the sample type, " << sample_type_ << ", is not defined. Use \'high\' or \'low\'." << std::endl;
		runtime_assert( false );
  }
}

//destruction!
PointMutationCalculator::~PointMutationCalculator(){}

//creators
protocols::design_opt::PointMutationCalculatorOP
PointMutationCalculator::clone() const{
	return new PointMutationCalculator( *this );
}

// setter - getter pairs
void
PointMutationCalculator::relax_mover( protocols::moves::MoverOP mover ){
	relax_mover_ = mover;
}

protocols::moves::MoverOP
PointMutationCalculator::relax_mover() const{
	return relax_mover_;
}

void
PointMutationCalculator::filters( vector1< protocols::filters::FilterOP > filters ){
	filters_ = filters;
}

vector1< protocols::filters::FilterOP > PointMutationCalculator::filters() const{
	return filters_;
}

void
PointMutationCalculator::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

core::pack::task::TaskFactoryOP
PointMutationCalculator::task_factory() const
{
	return task_factory_;
}

void
PointMutationCalculator::dump_pdb( bool const dump_pdb ){
  dump_pdb_ = dump_pdb;
}

bool
PointMutationCalculator::dump_pdb() const{
  return dump_pdb_;
}

void
PointMutationCalculator::sample_type( std::string const sample_type ){
  sample_type_ = sample_type;
}

std::string
PointMutationCalculator::sample_type() const{
  return sample_type_;
}

void
PointMutationCalculator::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionOP
PointMutationCalculator::scorefxn() const{
	return scorefxn_;
}

/*
//utility funxns for comparing values in sort
bool
cmp_pair_by_second(
	pair< AA, Real > const pair1,
	pair< AA, Real > const pair2 )
{
	return pair1.second < pair2.second;
}

bool
cmp_triple_by_third(
	pair< Size, pair< AA, Real > > const pair1,
	pair< Size, pair< AA, Real > > const pair2 )
{
	return pair1.second.second < pair2.second.second;
}

bool
cmp_pair_vec_by_first_vec_val(
  pair< Size, vector1< pair< AA, Real > > > const pair1,
  pair< Size, vector1< pair< AA, Real > > > const pair2 )
{
  return pair1.second[ 1 ].second < pair2.second[ 1 ].second;
}
*/

//backcompatibility; overloaded interface that allows the same data struct but wth one val/aa instead of a vector
void
PointMutationCalculator::calc_point_mut_filters(
	pose::Pose const & pose,
	vector1< pair< Size, vector1< pair< AA, Real > > > > & seqpos_aa_val_vec
)
{
	//call the default with a new tmp container
	vector1< pair< Size, vector1< pair< AA, vector1< Real > > > > > seqpos_aa_vals_vec;
	calc_point_mut_filters( pose, seqpos_aa_vals_vec );
	//iter thru tmp container and use vals to populate input container
	for( vector1< pair< core::Size, vector1< pair< AA, vector1< Real > > > > >::iterator seqpos_aa_vals = seqpos_aa_vals_vec.begin();
			seqpos_aa_vals != seqpos_aa_vals_vec.end(); ++seqpos_aa_vals ){
		Size seqpos( seqpos_aa_vals->first );
		//seqpos_aa_vals->second is a seqpos' vector of aa/vals pairs
		assert( !seqpos_aa_vals->second.empty() );
		vector1< pair< AA, Real > > aa_val;
		for( vector1< pair< AA, vector1< Real > > >::iterator aa_vals = seqpos_aa_vals->second.begin();
				aa_vals != seqpos_aa_vals->second.end(); ++aa_vals ){
			//aa_vals->second is an aa's vector of vals
			assert( !aa_vals->second.empty() );
			aa_val.push_back( pair< AA, Real >( aa_vals->first, aa_vals->second[ 1 ] ) );
		}
		seqpos_aa_val_vec.push_back( pair< Size, vector1< pair< AA, Real > > >( seqpos, aa_val ) );
	}
}

void
PointMutationCalculator::calc_point_mut_filters(
	pose::Pose const & start_pose,
	vector1< pair< Size, vector1< pair< AA, vector1< Real > > > > > & seqpos_aa_vals_vec
)
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	//create vec of pairs of seqpos, vector of AA/vals pairs that pass input filter
	//clear the input first!
	seqpos_aa_vals_vec.clear();
	//apply input task factory to pose
//	PackerTaskCOP task = task_factory()->create_task_and_apply_taskoperations( start_pose );
	PackerTaskCOP task = task_factory_->create_task_and_apply_taskoperations( start_pose );

	//get vector< Size > of protein seqpos being designed
	vector1< core::Size > being_designed;
	being_designed.clear();
	for( core::Size resi = 1; resi <= start_pose.total_residue(); ++resi ){
		if( task->residue_task( resi ).being_designed() && start_pose.residue(resi).is_protein() )
			being_designed.push_back( resi );
	}
	if( being_designed.empty() ) {
		TR.Warning << "WARNING: No residues are listed as designable." << std::endl;
	}

	//for each seqpos in being_designed vector
	foreach( core::Size const resi, being_designed ){
		//create vector< AA > of allowed residue types at seqpos
		typedef std::list< ResidueTypeCAP > ResidueTypeCAPList;
		ResidueTypeCAPList const & allowed( task->residue_task( resi ).allowed_residue_types() );
		vector1< AA > allow_temp;
		allow_temp.clear();
		foreach( ResidueTypeCAP const t, allowed ){
			if(std::find(allow_temp.begin(),allow_temp.end(),t->aa())!=allow_temp.end()) continue;
			allow_temp.push_back( t->aa() );
		}
		//temp store vector of aa/val pairs
		vector1< pair< AA, vector1< Real > > > aa_vals;
		//for each allowed AA
		foreach( AA const target_aa, allow_temp ){
			//make copy of original
			pose::Pose pose( start_pose );
			//make bool vector of allowed aa's [1,20], all false except for target_aa
			vector1< bool > allowed_aas;
			allowed_aas.clear();
			allowed_aas.assign( num_canonical_aas, false );
			allowed_aas[ target_aa ] = true;
			//make mut_res task factory by copying input task_factory,
			//then restrict to mutates resi to target_aa and repack 8A shell
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
			//store aa/val pair in seqpos_aa_vals_vec
			//TODO: if no filter defined, just use total_score
			bool filter_pass( true );
			vector1< Real > vals;
			for( vector1< protocols::filters::FilterOP >::iterator filter = filters_.begin(); filter != filters_.end(); ++filter ){
				//check if this filter passes, AND it with current value of pass/fail
				filter_pass = filter_pass && ( *filter )->apply( pose );
				//val sign is switched if type is high
				Real const val( flip_sign_ * ( *filter )->report_sm( pose ) );
				//bail at first fail
				if( !filter_pass ){
					TR<<"Filter fails with value "<< val << std::endl;
					break;
				}
				TR<<"Filter succeeds with value "<< val << std::endl;
				vals.push_back( val );
			}
			//don't store this aa/val if any filter failed
			if( !filter_pass ) continue;

			//filter passed, store this aa/vals pair
			assert( !vals.empty() );
			aa_vals.push_back( pair< AA, vector1< Real > >( target_aa, vals ) );
			//dump pdb?
			if( dump_pdb() ){
				using namespace protocols::jd2;
				JobOP job( JobDistributor::get_instance()->current_job() );
				std::stringstream fname;
				fname << job->input_tag() << pose.residue( resi ).name3() << resi << pose.residue( resi ).name3()<<".pdb";
				TR<<"Saving pose "<<fname.str() << std::endl;
				pose.dump_scored_pdb( fname.str(), *scorefxn() );
			}
			TR.flush();
		}//foreach target_aa
		//store the aa/vals for this seqpos in the big struct if there are any
		if( !aa_vals.empty() ){
			seqpos_aa_vals_vec.push_back( pair< Size, vector1< pair< AA, vector1< Real > > > >( resi, aa_vals ) );
		}
	}//foreach resi

/*
	//this part sorts the seqpos/aa/val data
	//first over each seqpos by aa val, then over all seqpos by best aa val
	for( vector1< pair< core::Size, vector1< pair< AA, Real > > > >::iterator aa_vals = seqpos_aa_vals_vec.begin();
			aa_vals != seqpos_aa_vals_vec.end(); ++aa_vals ){
		//skip if aa_vals vector is empty
		if( aa_vals->second.empty() ) continue;
		//sort aa_vals in incr val order
		std::sort( aa_vals->second.begin(), aa_vals->second.end(), cmp_pair_by_second );
		//best val is lowest, store all in sorted_seqpos_aa_vals_vec
//		pair< AA, Real > best_resid_val( aa_vals->second[ 1 ] );
		//create the pair of seqpos and sorted AA/val pairs
		pair< Size, vector1< pair< AA, Real > > > sorted_aa_vals( aa_vals->first, aa_vals->second );
		seqpos_aa_vals_vec.push_back( sorted_aa_vals );
	}

	//now sort seqpos_aa_vals_vec by *first* (lowest) val in each seqpos vector, low to high
	//uses cmp_pair_vec_by_first_vec_val to sort based on second val in
	//first pair element of pair( size, vec( pair ) )
	std::sort( seqpos_aa_vals_vec.begin(), seqpos_aa_vals_vec.end(), cmp_pair_vec_by_first_vec_val );
*/

}

} // design_opt
} // protocols
