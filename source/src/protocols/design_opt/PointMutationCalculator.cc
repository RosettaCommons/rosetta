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
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/string_util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/util.hh>
#include <boost/foreach.hpp>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/GreenPacker.hh>
#include <core/pose/symmetry/util.hh>

//Auto Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#ifdef USEMPI
#include <mpi.h>
#include <utility/mpi_util.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#endif

namespace protocols {
namespace design_opt {

static thread_local basic::Tracer TR( "protocols.design_opt.PointMutationCalculator" );
using namespace core;
using namespace chemical;
using utility::vector1;
using utility::vector0;
using std::pair;

///@brief default ctor
PointMutationCalculator::PointMutationCalculator() :
	task_factory_( NULL ),
	scorefxn_( NULL ),
	relax_mover_( NULL ),
//	filters_( NULL ), /*TODO: this throws a warning!*/
//	sample_type_( "low" )
	dump_pdb_( false ),
	rtmin_( false ),
	parallel_( false ),
	design_shell_( -1.0),
	repack_shell_( 8.0 )
{}

//full ctor
PointMutationCalculator::PointMutationCalculator(
	core::pack::task::TaskFactoryOP task_factory,
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::moves::MoverOP relax_mover,
	vector1< protocols::filters::FilterOP > filters,
	utility::vector1< std::string > sample_types,
	bool dump_pdb,
	bool rtmin,
	bool parallel,
	core::Real design_shell,
	core::Real repack_shell
)
{
	task_factory_ = task_factory;
	filters_ = filters;
	relax_mover_ = relax_mover;
	scorefxn_ = scorefxn;
	dump_pdb_ = dump_pdb;
	rtmin_ = rtmin;
	parallel_ = parallel;
	sample_types_ = sample_types;
	design_shell_ = design_shell;
	repack_shell_ = repack_shell;

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
	bool dump_pdb,
	bool rtmin,
	bool parallel,
	core::Real design_shell,
	core::Real repack_shell

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
	rtmin_ = rtmin;
	parallel_ = parallel;
	sample_types_ = sample_types;
	design_shell_ = design_shell;
	repack_shell_ = repack_shell;

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

void
PointMutationCalculator::set_design_shell( core::Real dz_shell ){
	design_shell_ = dz_shell;
}

void
PointMutationCalculator::set_repack_shell( core::Real rp_shell ){
	repack_shell_ = rp_shell;
}

bool
PointMutationCalculator::rtmin() const{
	return rtmin_;
}

void
PointMutationCalculator::rtmin( bool const r ){ rtmin_ = r;}

bool
PointMutationCalculator::parallel() const{
	return parallel_;
}

void
PointMutationCalculator::parallel( bool const r ){ parallel_ = r;}

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

#ifdef USEMPI
Size
get_nstruct(){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  if ( option[ run::shuffle ]() ) {
    return option[ out::shuffle_nstruct ]();
  } else {
    return option[ out::nstruct ]();
  }
}
#endif

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
	repack_around_op->design_shell( design_shell_ ); //neg radius insures no designing nbrs, positive will do so!
	repack_around_op->repack_shell( repack_shell_ );
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
	if( rtmin ){
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
	repack_around_op->design_shell( design_shell_ ); //neg radius insures no designing nbrs
	repack_around_op->repack_shell( repack_shell_ );
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
	TR << "Mutation " << pose.residue( resi ).name1() << "_" << resi;
	//only use green packer if not symmetric!
	assert( !core::pose::symmetry::is_symmetric( pose ) );
	green_packer->set_task_factory( mut_res );
	green_packer->apply( pose );
	TR << "_" << pose.residue( resi ).name1();
	//then run input relax mover
	relax_mover()->apply( pose );
}

//TODO HEY!~ I think we're evalling the filter TWICE here, once for pass, once for value
// can't we just get both with one calculation?
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
		bool this_filter_pass( ( filters_[ ifilt ] )->apply( pose ) );
		filter_pass = filter_pass && this_filter_pass;
		//val sign is switched if type is high
		Real const flip_sign( sample_types_[ ifilt ] == "high" ? -1 : 1 );
		Real const val( flip_sign * ( filters_[ ifilt ] )->report_sm( pose ) );
		//TODO: option to bail at first fail??
		TR<< " :: Filter " << ifilt;
		if( !this_filter_pass ) TR << " fail, ";
		else TR << " pass, ";
		TR << " value "<< val << " ::";
		vals.push_back( val );
	}
	TR << std::endl;
}

void
insert_point_mut_filter_vals(
	Size const seqpos,
	chemical::AA const aa,
	vector1< Real > const vals,
	vector1< pair< Size, vector1< pair< AA, vector1< Real > > > > > & seqpos_aa_vals_vec
)
{
	using namespace core::chemical;
	//create the aa,vals pair
	pair< AA, vector1< Real > > aa_vals_pair( pair< AA, vector1< Real > >( aa, vals ) );
	//first check if we've assigned anything for seqpos
	//if we have, just append this aa,vals pair onto that seqpos' data
	bool new_pos( true ); // It's a new position unless we find it
	for( Size iseq = 1; iseq <= seqpos_aa_vals_vec.size(); ++iseq ){
		if( seqpos == seqpos_aa_vals_vec[ iseq ].first ){
			new_pos = false;
			bool replaced( false );
			//we need to check if we already have vals for this seqpos,aa in our data
			for( core::Size iaa = 1; iaa <= seqpos_aa_vals_vec[ iseq ].second.size(); ++iaa ){
				char this_aa_char( chemical::oneletter_code_from_aa( seqpos_aa_vals_vec[ iseq ].second[ iaa ].first ) );
				if( this_aa_char == chemical::oneletter_code_from_aa( aa ) ){
					seqpos_aa_vals_vec[ iseq ].second[ iaa ].second = vals;
					replaced = true;
				}
			}
			//dont append new data if we're just replacing
			if( replaced ) break;
			seqpos_aa_vals_vec[ iseq ].second.push_back( aa_vals_pair );
			break;
		}
	}
	//if this is the first instance of data at seqpos,
	//create a 1-element vector and add to the ptmut data
	if( new_pos ){
		vector1< pair< AA, vector1< Real > > > aa_vals_vec( 1, aa_vals_pair );
		seqpos_aa_vals_vec.push_back( pair< Size, vector1< pair< AA, vector1< Real > > > >( seqpos, aa_vals_vec ) );
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

#ifdef USEMPI
	int mpi_rank( 0 ), mpi_nprocs( 1 ), mpi_rank_low( 0 );
	if( parallel() ){
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_nprocs);
		//Get the lowest rank proc that's running this mover
		if( dynamic_cast< protocols::jd2::MPIWorkPoolJobDistributor* >( protocols::jd2::JobDistributor::get_instance() ) ){
			//!!WARNING!! We're assuming we have one head node (0) and nprocs-1 workers !!WARNING!!
			TR << "Detected jd2::MPIWorkPoolJobDistributor... excluding proc 0 from calculations" << std::endl;
			mpi_rank_low = 1;
			//We must have one job( nstruct ) for each worker in pool or we'll freeze later because nodes w/ no job will get killed by jd2
			if( (int) get_nstruct() < mpi_nprocs - mpi_rank_low ) utility_exit_with_message(
					"You must specify nstruct >= " + utility::to_string( mpi_nprocs - mpi_rank_low ) +
					" when using " + utility::to_string( mpi_nprocs ) + " processors for MPI PointMutationCalculator" +
					" when called from rosetta_scripts or any other app using jd2::MPIWorkPoolJobDistributor!" );
		}
		else if( dynamic_cast< protocols::jd2::MPIFileBufJobDistributor* >(protocols::jd2::JobDistributor::get_instance() ) ){
			protocols::jd2::MPIFileBufJobDistributor* jd2 =
					dynamic_cast< protocols::jd2::MPIFileBufJobDistributor* >( protocols::jd2::JobDistributor::get_instance() );
			mpi_rank_low = jd2->min_client_rank();
			TR << "Detected jd2::MPIFileBufJobDistributor... excluding procs 0-" << ( mpi_rank_low - 1 ) << " from calculations" << std::endl;
			//We must have one job( nstruct ) for each worker in pool or we'll freeze later because nodes w/ no job will get killed by jd2
			if( (int) get_nstruct() < mpi_nprocs - mpi_rank_low ) utility_exit_with_message(
					"You must specify nstruct >= " + utility::to_string( mpi_nprocs - mpi_rank_low ) +
					" when using " + utility::to_string( mpi_nprocs ) + " processors for MPI PointMutationCalculator" +
					" when called from rosetta_scripts or any other app using jd2::MPIFileBufJobDistributor!" );
		}
	/*
		//create a group of worker nodes and then a communicator
		MPI_Group mpi_pool_group, mpi_all_group;
		MPI_Comm MPI_COMM_POOL;
		int returnval;
		//create mpi_all_group group
		returnval = MPI_Comm_group( MPI_COMM_WORLD, &mpi_all_group);
		if ( returnval != MPI_SUCCESS ) utility_exit_with_message("failed in creating a new communicator!");
		//create the pool group
		// ranks is node ranks to include in your new group
		int const mpi_pool_nprocs( mpi_nprocs - mpi_rank_low );
		int ranks[ mpi_pool_nprocs ];
		for( int irank = 0; irank < mpi_pool_nprocs; ++irank ){
			ranks[ irank ] = mpi_rank_low + irank;
			TR << "MPI group incl Proc " << mpi_rank_low + irank << std::endl;
		}
		TR << "Creating pool group..." << std::endl;
		returnval = MPI_Group_incl( mpi_all_group, mpi_pool_nprocs, ranks, &mpi_pool_group );
		if ( returnval != MPI_SUCCESS ) utility_exit_with_message("failed in creating a new communicator!");
		TR << "Creating comm group..." << std::endl;
		returnval = MPI_Comm_create( MPI_COMM_WORLD, mpi_pool_group, &MPI_COMM_POOL );
		if ( returnval != MPI_SUCCESS ) utility_exit_with_message("failed in creating a new communicator!");
		TR << "MPI Comm created!" << std::endl;
	*/
	}
#endif

	//make a single list of seqpos,aa pairs
	vector1< pair< Size, AA > > all_muts;
	for( Size iresi = 1; iresi <= being_designed.size(); ++iresi ){
		Size const resi( being_designed[ iresi ] );
		//create vector< AA > of allowed residue types at seqpos
		typedef std::list< ResidueTypeCOP > ResidueTypeCOPList;
		ResidueTypeCOPList const & allowed( task->residue_task( resi ).allowed_residue_types() );
		vector1< AA > allow_temp;
		BOOST_FOREACH( ResidueTypeCOP const t, allowed ){
			if(std::find(allow_temp.begin(),allow_temp.end(),t->aa())!=allow_temp.end()) continue;
			allow_temp.push_back( t->aa() );
		}
		//for each allowed AA
		BOOST_FOREACH( AA const target_aa, allow_temp ){
			all_muts.push_back( pair< Size, AA >( resi, target_aa ) );
		}
	}

	vector1< pair< Size, AA > > my_muts( all_muts );
	if( parallel() ){
#ifdef USEMPI
		//split up my_muts into smaller sublists for diff procs
		my_muts.clear();
		//asign muts to each proc
		for( int imut = 1; imut <= (int) all_muts.size(); ++imut ){
			//e.g. for 4 procs, hand out muts like 1,2,3,1,2,3,etc (nothing given to proc 0)
			Size this_mpi_rank( ( imut - 1 ) % ( mpi_nprocs - mpi_rank_low ) + mpi_rank_low );
			if( (int) this_mpi_rank == mpi_rank ){
				my_muts.push_back( all_muts[ imut ] );
			}
		}
		//TR << "Proc " << mpi_rank << " takes mutations: ";
		//for( Size imut = 1; imut <= my_muts.size(); ++imut ) TR << my_muts[ imut ].first << my_muts[ imut ].second << " ";
		//TR << std::endl;
#endif
	}

	for( Size imut = 1; imut <= my_muts.size(); ++imut ){
		Size seqpos( my_muts[ imut ].first );
		AA target_aa( my_muts[ imut ].second );
		//make copy of original
		pose::Pose pose( start_pose );
		//make the mutation and relax
		//then check if passes input filter, bail out if it doesn't
		bool filter_pass;
		vector1< Real > vals;
		if( use_precomp_rot_pair_nrgs ) mutate_and_relax( pose, seqpos, target_aa, green_packer );
		else mutate_and_relax( pose, seqpos, target_aa );
//			mutate_and_relax( pose, seqpos, target_aa );
		eval_filters( pose, filter_pass, vals );

		//don't store this aa/val if any filter failed
		if( !filter_pass ) continue;
		assert( !vals.empty() );
		//dump pdb? (only if filter passes)
		if( dump_pdb() ){
			std::stringstream fname;
			fname << protocols::jd2::current_output_name() << start_pose.residue( seqpos ).name3() << seqpos << pose.residue( seqpos ).name3()<<".pdb";
			TR<<"Saving pose "<<fname.str() << std::endl;
			pose.dump_scored_pdb( fname.str(), *scorefxn() );
		}
		TR.flush();
		insert_point_mut_filter_vals( seqpos, target_aa, vals, seqpos_aa_vals_vec );
	}//for mut

	if( parallel() ){
#ifdef USEMPI
		//MPI_Barrier( MPI_COMM_POOL );
		//sync everybody's mutation filter data
		//worker sends ptmut data to pool leader
		if( mpi_rank > mpi_rank_low ){
			utility::send_integer_to_node( mpi_rank_low, seqpos_aa_vals_vec.size() );	//send int
			for( Size iseq = 1; iseq <= seqpos_aa_vals_vec.size(); ++iseq ){
				utility::send_integer_to_node( mpi_rank_low, seqpos_aa_vals_vec[ iseq ].first );	//send int
				utility::vector1< std::pair< core::chemical::AA, vector1< core::Real > > > const & aa_pairs( seqpos_aa_vals_vec[ iseq ].second );
				utility::send_integer_to_node( mpi_rank_low, aa_pairs.size() );	//send int
				for( core::Size iaa = 1; iaa <= aa_pairs.size(); ++iaa ){
					utility::send_char_to_node( mpi_rank_low, chemical::oneletter_code_from_aa( aa_pairs[ iaa ].first ) );	//send char
					//TR << "Proc " << mpi_rank << " sending seqpos,aa: " << seqpos_aa_vals_vec[ iseq ].first << aa_pairs[ iaa ].first
					//		<< ": " << aa_pairs[ iaa ].second[ 1 ] << " to Proc 0" << std::endl;
					for( Size ival = 1; ival <= ( filters() ).size(); ++ival ){
						utility::send_double_to_node( mpi_rank_low, aa_pairs[ iaa ].second[ ival ] );	//send Real
					}
				}
			}
		}
		//pool leader receives ptmut data from workers and combines with its own
		else if( mpi_rank == mpi_rank_low ){
			for( int iproc = mpi_rank_low + 1; iproc < mpi_nprocs; ++iproc ){
				//get data for one mut (seqpos, AA, and filter vals )
				//need to know how many seqpos
				Size n_seqpos( utility::receive_integer_from_node( iproc ) ); //rec int
				for( Size imut = 1; imut <= n_seqpos; ++imut ){
					Size seqpos( utility::receive_integer_from_node( iproc ) ); //rec int
					//need to know how many muts at this seqpos
					Size n_aas( utility::receive_integer_from_node( iproc ) );	//rec int
					for( Size iaa = 1; iaa <= n_aas; ++iaa ){
						char aa_char( utility::receive_char_from_node( iproc ) );	//rec char
						chemical::AA aa( chemical::aa_from_oneletter_code( aa_char ) );
						vector1< Real > vals( ( filters() ).size(), 0. );
						for( Size ival = 1; ival <= vals.size(); ++ival ){
							vals[ ival ] = utility::receive_double_from_node( iproc );	//rec Real
						}
						//TR << "Proc " << mpi_rank << " received seqpos,aa: " << seqpos << aa << ": "
						//		<< vals[ 1 ] << " from Proc " << iproc << std::endl;
						insert_point_mut_filter_vals( seqpos, aa, vals, seqpos_aa_vals_vec );
					}
				}
			}
		}
		//MPI_Barrier( MPI_COMM_POOL );
		//then pool leader sends combined ptmut data back to workers
		if( mpi_rank == mpi_rank_low ){
			for( int iproc = mpi_rank_low + 1; iproc < mpi_nprocs; ++iproc ){
				utility::send_integer_to_node( iproc, seqpos_aa_vals_vec.size() );	//send int
				for( Size iseq = 1; iseq <= seqpos_aa_vals_vec.size(); ++iseq ){
					utility::send_integer_to_node( iproc, seqpos_aa_vals_vec[ iseq ].first );	//send int
					utility::vector1< std::pair< core::chemical::AA, vector1< core::Real > > > const & aa_pairs( seqpos_aa_vals_vec[ iseq ].second );
					utility::send_integer_to_node( iproc, aa_pairs.size() );	//send int
					for( core::Size iaa = 1; iaa <= aa_pairs.size(); ++iaa ){
						utility::send_char_to_node( iproc, chemical::oneletter_code_from_aa( aa_pairs[ iaa ].first ) );	//send char
						//TR << "Proc " << mpi_rank << " sending seqpos,aa: " << seqpos_aa_vals_vec[ iseq ].first << aa_pairs[ iaa ].first << std::endl;
						for( Size ival = 1; ival <= ( filters() ).size(); ++ival ){
							utility::send_double_to_node( iproc, aa_pairs[ iaa ].second[ ival ] );	//send Real
						}
					}
				}
			}
		}
		//workers receive combined ptmut data from pool leader
		else if( mpi_rank > mpi_rank_low ){
			//need to know how many seqpos
			Size n_seqpos( utility::receive_integer_from_node( mpi_rank_low ) ); //rec int
			for( Size imut = 1; imut <= n_seqpos; ++imut ){
				Size seqpos( utility::receive_integer_from_node( mpi_rank_low ) ); //rec int
				//need to know how many muts at this seqpos
				Size n_aas( utility::receive_integer_from_node( mpi_rank_low ) );	//rec int
				for( Size iaa = 1; iaa <= n_aas; ++iaa ){
					char aa_char( utility::receive_char_from_node( mpi_rank_low ) );	//rec char
					chemical::AA aa( chemical::aa_from_oneletter_code( aa_char ) );
					//TR << "Proc " << mpi_rank << " received seqpos,aa: " << seqpos << aa << std::endl;
					vector1< Real > vals( ( filters() ).size(), 0. );
					for( Size ival = 1; ival <= vals.size(); ++ival ){
						vals[ ival ] = utility::receive_double_from_node( mpi_rank_low );	//rec Real
					}
					insert_point_mut_filter_vals( seqpos, aa, vals, seqpos_aa_vals_vec );
				}
			}
		}
		//	MPI_Finalize();
#endif
	}

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
