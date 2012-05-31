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
#include <numeric/random/random.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/GreenPacker.hh>
#include <utility/vector0.hh>
#include <core/pose/symmetry/util.hh>

//Auto Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

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
	scorefxn_( NULL ),
	relax_mover_( NULL ),
//	filters_( NULL ), /*TODO: this throws a warning!*/
//	sample_type_( "low" )
	dump_pdb_( false ),
	rtmin_( false )
{}

//full ctor
PointMutationCalculator::PointMutationCalculator(
	core::pack::task::TaskFactoryOP task_factory,
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::moves::MoverOP relax_mover,
	vector1< protocols::filters::FilterOP > filters,
	utility::vector1< std::string > sample_types,
	bool dump_pdb
)
{
	task_factory_ = task_factory;
	filters_ = filters;
	relax_mover_ = relax_mover;
	scorefxn_ = scorefxn;
	dump_pdb_ = dump_pdb;
	sample_types_ = sample_types;

	for( Size isamp = 1; isamp <= sample_types.size(); ++isamp ){
		if( sample_types_[ isamp ] != "high" && sample_types_[ isamp ] != "low" ){
			TR << "WARNING: the sample type, " << sample_types_[ isamp ] << ", is not defined. Use \'high\' or \'low\'." << std::endl;
			runtime_assert( false );
		}
	}
}

//backcompatible; single-filter convenience ctor
PointMutationCalculator::PointMutationCalculator(
	core::pack::task::TaskFactoryOP task_factory,
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::moves::MoverOP relax_mover,
	protocols::filters::FilterOP filter,
	std::string sample_type,
	bool dump_pdb
)
{
	vector1< protocols::filters::FilterOP > filters;
	filters.push_back( filter );
	vector1< std::string > sample_types;
	sample_types.push_back( sample_type );
	task_factory_ = task_factory;
	filters_ = filters;
	relax_mover_ = relax_mover;
	scorefxn_ = scorefxn;
	dump_pdb_ = dump_pdb;
	sample_types_ = sample_types;

	for( Size isamp = 1; isamp <= sample_types.size(); ++isamp ){
		if( sample_types_[ isamp ] != "high" && sample_types_[ isamp ] != "low" ){
			TR << "WARNING: the sample type, " << sample_types_[ isamp ] << ", is not defined. Use \'high\' or \'low\'." << std::endl;
			runtime_assert( false );
		}
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
PointMutationCalculator::sample_types( vector1< std::string > const sample_types ){
  sample_types_ = sample_types;
}

vector1< std::string >
PointMutationCalculator::sample_types() const{
  return sample_types_;
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

void
PointMutationCalculator::mutate_and_relax(
	pose::Pose & pose,
	Size const & resi,
	AA const & target_aa
){
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

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
	protocols::simple_moves::RotamerTrialsMinMoverOP rtmin;
	if( core::pose::symmetry::is_symmetric( pose ) ) {
		mutate_residue->request_symmetrize_by_union();
		pack = new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn(), mutate_residue );
	} else {
		pack = new protocols::simple_moves::PackRotamersMover( scorefxn(), mutate_residue );
	}
	pack->apply( pose );
	if( rtmin() ){
		// definition/allocation of RTmin mover must flag dependant, as some scoreterms are incompatable with RTmin initilization
		if( core::pose::symmetry::is_symmetric( pose ) ) {
			utility_exit_with_message("Cannot currently use PointMutationCalculator (GreedyOptMutation/ParetoOptMutation) with rtmin on a symmetric pose!");
		}
		rtmin = new protocols::simple_moves::RotamerTrialsMinMover( scorefxn(), *mutate_residue );
		rtmin->apply( pose );
		TR<<"Finished rtmin"<<std::endl;
	}
	TR<<pose.residue( resi ).name3()<<". Now relaxing..."<<std::endl;
	//then run input relax mover
	if( relax_mover() ) {
		relax_mover()->apply( pose );
	}
}

void
PointMutationCalculator::mutate_and_relax(
	pose::Pose & pose,
	Size const & resi,
	AA const & target_aa,
	protocols::simple_moves::GreenPackerOP green_packer
){
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

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
	core::pack::task::operation::RestrictAbsentCanonicalAASOP restrict_to_aa_op =
			new core::pack::task::operation::RestrictAbsentCanonicalAAS;
	restrict_to_aa_op->include_residue( resi );
	restrict_to_aa_op->keep_aas( allowed_aas );
	core::pack::task::operation::InitializeFromCommandlineOP init_from_cmd_op =
			new core::pack::task::operation::InitializeFromCommandline;
	core::pack::task::operation::IncludeCurrentOP incl_curr_op =
		new core::pack::task::operation::IncludeCurrent;
	mut_res->push_back( repack_around_op );
	mut_res->push_back( restrict_to_aa_op );
	mut_res->push_back( init_from_cmd_op );
	mut_res->push_back( incl_curr_op );
//	TR<<"Mutating residue "<<pose.residue( resi ).name3()<<resi<<" to " << target_aa << std::endl;
//	TR<<"Mutating residue "<<pose.residue( resi ).name3()<<resi<<" to ";
	//only use green packer if not symmetric!
	assert( !core::pose::symmetry::is_symmetric( pose ) );
	green_packer->set_task_factory( mut_res );
	green_packer->apply( pose );
//	TR<<"Now relaxing "<<pose.residue( resi ).name3()<<std::endl;
//	TR<<pose.residue( resi ).name3()<<". Now relaxing..."<<std::endl;
	//then run input relax mover
	relax_mover()->apply( pose );
}

void
PointMutationCalculator::eval_filters(
	pose::Pose & pose,
	bool & filter_pass,
	vector1< Real > & vals
){
	//now run filters
	filter_pass = true;
	vals.clear();
	for( Size ifilt = 1; ifilt <= filters_.size(); ++ifilt ){
		//check if this filter passes, AND it with current value of pass/fail
		filter_pass = filter_pass && ( filters_[ ifilt ] )->apply( pose );
		//val sign is switched if type is high
		Real const flip_sign( sample_types_[ ifilt ] == "high" ? -1 : 1 );
		Real const val( flip_sign * ( filters_[ ifilt ] )->report_sm( pose ) );
		//TODO: option to bail at first fail??
		if( !filter_pass )
			TR<< "Filter " << ifilt << " fails with value "<< val << std::endl;
		else
			TR<< "Filter " << ifilt << " succeeds with value "<< val<< std::endl;
		vals.push_back( val );
	}
}

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
	PackerTaskOP tmptask = task_factory_->create_task_and_apply_taskoperations( start_pose );
	if(core::pose::symmetry::is_symmetric(start_pose)) {
		tmptask->request_symmetrize_by_union();
		tmptask = core::pack::make_new_symmetric_PackerTask_by_requested_method(start_pose,tmptask);
	}
	PackerTaskCOP task = tmptask;
	//get vector< Size > of protein seqpos being designed and group id mask for green packer (0 is designed, 1 is rp only)
	vector1< core::Size > being_designed;
	utility::vector1< Size > group_ids;
	being_designed.clear();
	for( core::Size resi = 1; resi <= start_pose.total_residue(); ++resi ){
		if(core::pose::symmetry::is_symmetric(start_pose)) {
			if( resi > core::pose::symmetry::symmetry_info(start_pose)->num_independent_residues() ) {
				break;
			}
		}
		if( task->residue_task( resi ).being_designed() && start_pose.residue(resi).is_protein() ){
			being_designed.push_back( resi );
			group_ids.push_back( 0 );
		} else{
			group_ids.push_back( 1 );
		}
	}
	if( being_designed.empty() ) {
		TR.Warning << "WARNING: No residues are listed as designable." << std::endl;
	}

	//GreenPacker stuff, precompute non-designable residues' interaxn graph
	//dont use green packer if symmetric (symm not supported for green packer)
	//dont use green packer if user specifies linmem interaxn graph (is calc on the fly)
	bool use_precomp_rot_pair_nrgs( true );
	if( basic::options::option[ basic::options::OptionKeys::packing::linmem_ig ].user() ){
		TR << "Note: you are using linmem_ig in your options: " <<
				"packing will be slower because GreedyOpt can't use GreenPacker precomputed rotamer pair energies" << std::endl;
		use_precomp_rot_pair_nrgs = false;
	}
	if( core::pose::symmetry::is_symmetric( start_pose ) ){
		TR << "Note: you are using symmetry: " <<
				"packing will be slower because GreedyOpt can't use GreenPacker precomputed rotamer pair energies" << std::endl;
		use_precomp_rot_pair_nrgs = false;
	}
	protocols::simple_moves::UserDefinedGroupDiscriminatorOP user_defined_group_discriminator(
			new protocols::simple_moves::UserDefinedGroupDiscriminator );
	user_defined_group_discriminator->set_group_ids( group_ids );
	protocols::simple_moves::GreenPackerOP green_packer( new protocols::simple_moves::GreenPacker );
	green_packer->set_group_discriminator( user_defined_group_discriminator );
	green_packer->set_scorefunction( *scorefxn() );
	green_packer->set_reference_round_task_factory( task_factory() );

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
			//make the mutation and relax
			//then check if passes input filter, bail out if it doesn't
			//TODO: if no filter defined, just use total_score
			bool filter_pass;
			vector1< Real > vals;
			if( use_precomp_rot_pair_nrgs ) mutate_and_relax( pose, resi, target_aa, green_packer );
			else mutate_and_relax( pose, resi, target_aa );
//			mutate_and_relax( pose, resi, target_aa );
			eval_filters( pose, filter_pass, vals );

			//don't store this aa/val if any filter failed
			if( !filter_pass ) continue;
			assert( !vals.empty() );
			//store aa/val pair in seqpos_aa_vals_vec
			aa_vals.push_back( pair< AA, vector1< Real > >( target_aa, vals ) );
			//dump pdb?
			if( dump_pdb() ){
				std::stringstream fname;
				fname << protocols::jd2::current_output_name() << start_pose.residue( resi ).name3() << resi << pose.residue( resi ).name3()<<".pdb";
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

bool
PointMutationCalculator::rtmin() const{
	return rtmin_;
}

void
PointMutationCalculator::rtmin( bool const r ){ rtmin_ = r;}



} // design_opt
} // protocols
