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
#include <protocols/design_opt/GreedyOptMutationMover.hh>
#include <protocols/design_opt/GreedyOptMutationMoverCreator.hh>
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
static numeric::random::RandomGenerator RG( 2718 );
using namespace core;
using namespace chemical;
using utility::vector1;
using std::pair;

///@brief default ctor
GreedyOptMutationMover::GreedyOptMutationMover() :
	Mover( GreedyOptMutationMoverCreator::mover_name() ),
	task_factory_( NULL ),
	scorefxn_( NULL ),
	relax_mover_( NULL ),
	filter_( NULL ),
	filter_delta_( 0 ),
	sample_type_( "low" ),
	dump_pdb_( false ),
	diversify_lvl_( 1 ),
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
	//filter_delta should always be a scalar!
	if( filter_delta_ < Real( 0 ) ) filter_delta_ = ( Real( -1 ) * filter_delta_ );
	//default diversity to all 20 aa's if specified filter_delta but did not spec diversify_lvl
	if( filter_delta_ != Real( 0 ) && diversify_lvl_ == Size( 1 ) ) diversify_lvl_ = 20;
}

//full ctor
GreedyOptMutationMover::GreedyOptMutationMover(
	core::pack::task::TaskFactoryOP task_factory,
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::moves::MoverOP relax_mover,
	protocols::filters::FilterOP filter,
	core::Real filter_delta,
	std::string sample_type,
	bool dump_pdb,
	core::Size diversify_lvl,
	protocols::filters::FilterOP stopping_condition
) :
	Mover( GreedyOptMutationMoverCreator::mover_name() )
{
	task_factory_ = task_factory;
	relax_mover_ = relax_mover;
	filter_ = filter;
	filter_delta_ = filter_delta;
	scorefxn_ = scorefxn;
	sample_type_ = sample_type;
	diversify_lvl_ = diversify_lvl;
	dump_pdb_ = dump_pdb;
	stopping_condition_ = stopping_condition;
	
	if( sample_type_ == "high" ){
		flip_sign_ = Real( -1 );
	}else if( sample_type_ == "low" ){
		flip_sign_ = Real( 1 );
	}else{
		TR << "WARNING: the sample type, " << sample_type_ << ", is not defined. Use \'high\' or \'low\'." << std::endl;
		runtime_assert( false );
  }
	//filter_delta should always be a scalar!
	if( filter_delta_ < Real( 0 ) ) filter_delta_ = ( Real( -1 ) * filter_delta_ );
}

//destruction!
GreedyOptMutationMover::~GreedyOptMutationMover(){}

//creators
protocols::moves::MoverOP
GreedyOptMutationMoverCreator::create_mover() const {
	return new GreedyOptMutationMover;
}

protocols::moves::MoverOP
GreedyOptMutationMover::clone() const{
	return new GreedyOptMutationMover( *this );
}

//name getters
std::string
GreedyOptMutationMoverCreator::keyname() const
{
	return GreedyOptMutationMoverCreator::mover_name();
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

// setter - getter pairs
void
GreedyOptMutationMover::relax_mover( protocols::moves::MoverOP mover ){
	relax_mover_ = mover;
	clear_cached_data();
}

protocols::moves::MoverOP
GreedyOptMutationMover::relax_mover() const{
	return relax_mover_;
}

void
GreedyOptMutationMover::filter( protocols::filters::FilterOP filter ){
	filter_ = filter;
	clear_cached_data();
}

protocols::filters::FilterOP GreedyOptMutationMover::filter() const{
	return filter_;
}

void
GreedyOptMutationMover::filter_delta( Real filter_delta ){
	filter_delta_ = filter_delta;
}

Real GreedyOptMutationMover::filter_delta() const{
	return filter_delta_;
}

void
GreedyOptMutationMover::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
	clear_cached_data();
}

core::pack::task::TaskFactoryOP
GreedyOptMutationMover::task_factory() const
{
	return task_factory_;
}

void
GreedyOptMutationMover::dump_pdb( bool const dump_pdb ){
  dump_pdb_ = dump_pdb;
	clear_cached_data();
}

bool
GreedyOptMutationMover::dump_pdb() const{
  return dump_pdb_;
}

void
GreedyOptMutationMover::sample_type( std::string const sample_type ){
  sample_type_ = sample_type;
	clear_cached_data();
}

std::string
GreedyOptMutationMover::sample_type() const{
  return sample_type_;
}

void
GreedyOptMutationMover::diversify_lvl( core::Size const diversify_lvl ){
  diversify_lvl_ = diversify_lvl;
}

core::Size
GreedyOptMutationMover::diversify_lvl() const{
  return diversify_lvl_;
}

void
GreedyOptMutationMover::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
	clear_cached_data();
}

core::scoring::ScoreFunctionOP
GreedyOptMutationMover::scorefxn() const{
	return scorefxn_;
}

//utility funxns for comparing values in sort
bool
cmp_pair_by_second(
	pair< AA, Real > const pair1,
	pair< AA, Real > const pair2 )
{
	return pair1.second < pair2.second;
}

bool
cmp_pair_vec_by_first_vec_val(
  pair< Size, vector1< pair< AA, Real > > > const pair1,
  pair< Size, vector1< pair< AA, Real > > > const pair2 )
{
  return pair1.second[ 1 ].second < pair2.second[ 1 ].second;
}

void
GreedyOptMutationMover::clear_cached_data(){
  seqpos_aa_val_vec_.clear();
  ref_pose_.clear();
}

//TODO: this should also compare fold trees
bool
GreedyOptMutationMover::pose_coords_are_same( core::pose::Pose const & pose1, core::pose::Pose const & pose2 )
{
	//first check for all restype match, also checks same number res
	if( !pose1.conformation().sequence_matches( pose2.conformation() ) ) return false;
	//then check for all coords identical
	for ( Size i = 1; i <= pose1.total_residue(); ++i ) {
		core::conformation::Residue const & rsd1( pose1.residue( i ) );
		core::conformation::Residue const & rsd2( pose2.residue( i ) );
		//check same n atoms
		if( rsd1.natoms() != rsd2.natoms() ) return false;
		//and coords
		for( Size ii = 1; ii <= rsd1.natoms(); ++ii ) {
			if( rsd1.xyz( ii ).x() != rsd2.xyz( ii ).x() ) return false;
			if( rsd1.xyz( ii ).y() != rsd2.xyz( ii ).y() ) return false;
			if( rsd1.xyz( ii ).z() != rsd2.xyz( ii ).z() ) return false;
		}    
	}    
	return true;
}

void
GreedyOptMutationMover::apply(core::pose::Pose & pose )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	//store input pose
	core::pose::Pose start_pose( pose );
	design_opt::PointMutationCalculatorOP ptmut_calc( new design_opt::PointMutationCalculator(
				task_factory(), scorefxn(), relax_mover(), filter(), sample_type(), dump_pdb() ) );

	//create vec of pairs of seqpos, vector of AA/val pairs that pass input filter
	//only calc the ptmut data and sort once per pose, not at every nstruct iteration
	//but how will we know if that data is still valid? what if the backbone has moved?
	//the best answer is to store the pose passed to apply in a private variable (ref_pose_)
	//and only calc and sort if ref_pose_ is still undef or doesnt match apply pose
	if( ref_pose_.empty() || !pose_coords_are_same( start_pose, ref_pose_ ) ){
		//clear cached data
		clear_cached_data();
		//and (re)set ref_pose_ to this pose
		ref_pose_ = start_pose;
		//get the point mut values
		ptmut_calc->calc_point_mut_filters( start_pose, seqpos_aa_val_vec_ );

		//this part sorts the seqpos/aa/val data
		//first over each seqpos by aa val, then over all seqpos by best aa val
		for( Size ivec = 1; ivec <= seqpos_aa_val_vec_.size(); ++ivec ){
			//skip if aa/vals vector is empty
			if( seqpos_aa_val_vec_[ ivec ].second.empty() ) continue;
			//sort aa/vals in incr val order
			std::sort( seqpos_aa_val_vec_[ ivec ].second.begin(),
					seqpos_aa_val_vec_[ ivec ].second.end(), cmp_pair_by_second );
		}
		//now sort seqpos_aa_val_vec_ by *first* (lowest) val in each seqpos vector, low to high
		//uses cmp_pair_vec_by_first_vec_val to sort based on second val in
		//first pair element of vector in pair( size, vec( pair ) )
		std::sort( seqpos_aa_val_vec_.begin(), seqpos_aa_val_vec_.end(), cmp_pair_vec_by_first_vec_val );

	}

	TR<<"Combining sorted independently optimal mutations… " << std::endl;
	//reset pose to original, init starting filter val
	//must use same relax mover before scoring so comparison is fair!
	pose = start_pose;
	relax_mover()->apply( pose );
	Real best_val( flip_sign_ * filter()->report_sm( pose ) );
	pose::Pose best_pose( pose );
	//now try the best AA at each position, in order
	for( Size iseq = 1; iseq <= seqpos_aa_val_vec_.size(); ++iseq ){
		pose = best_pose;
		//the seqpos index is the first part of the pair
		Size seqpos( seqpos_aa_val_vec_[ iseq ].first );
		//the best aa is the first part of the first element of the aa/val vector
		AA target_aa( seqpos_aa_val_vec_[ iseq ].second[ 1 ].first );

		//allow stochastic sampling of suboptimal restypes
		if( diversify_lvl() > 1 ){
			//smaller of user-def lvl and actual size of vector
			Size max_diversify_lvl( std::min( diversify_lvl_, seqpos_aa_val_vec_[ iseq ].second.size() ) );
			//ifdef filter_delta, redef max div lvl for this seqpos
			if( filter_delta() != Real( 0 ) ){
				Real best_val( seqpos_aa_val_vec_[ iseq ].second[ 1 ].second );
				for( Size iaa = 2; iaa <= max_diversify_lvl; ++iaa ){
					Real val( seqpos_aa_val_vec_[ iseq ].second[ iaa ].second );
					//if this aa is worse than filter_delta break out
					if( val - best_val > filter_delta() ) break;
					//else set max to this one
					else max_diversify_lvl = iaa;
				}
			}
			TR << "Randomly choosing 1 of " << max_diversify_lvl << " allowed mutations at residue " << seqpos << std::endl;
			Size aa_rank( static_cast< Size >( RG.uniform() * max_diversify_lvl + 1 ) );
			target_aa = seqpos_aa_val_vec_[ iseq ].second[ aa_rank ].first;
		}

		//dont need to make a "mutation" if target_aa is same as original aa
		if( target_aa == pose.residue( seqpos ).type().aa() ) continue;

		//then check if passes input filter, bail out if it doesn't
		bool filter_pass;
		vector1< Real > vals;
		ptmut_calc->mutate_and_relax( pose, seqpos, target_aa );
		ptmut_calc->eval_filters( pose, filter_pass, vals );
		Real const val( vals[ 1 ] );

		if( !filter_pass ) continue;
		//score mutation, reset best_pose, best val if is lower
		if( val > best_val ){
			TR << "Mutation " << start_pose.residue( seqpos ).name1() << seqpos << pose.residue( seqpos ).name1() << " rejected. Current best value is "<< best_val << std::endl;
			continue;
		}
		TR << "Mutation " << start_pose.residue( seqpos ).name1() << seqpos << pose.residue( seqpos ).name1() << " accepted. New best value is "<< val << std::endl;
		best_val = val;
		best_pose = pose;
		if( stopping_condition() && stopping_condition()->apply( pose ) ){
			TR<<"Stopping condition evaluates to true. Stopping early."<<std::endl;
			return;
		}
	}
//recover best pose after last step
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
	//load diversify_lvl
	diversify_lvl( tag->getOption< core::Size >( "diversify_lvl", core::Size( 1 ) ) );
	//load filter_delta
	filter_delta( tag->getOption< core::Size >( "filter_delta", core::Real( 0 ) ) );
	//filter_delta should always be a scalar!
	if( filter_delta() < Real( 0 ) ) filter_delta( -1 * filter_delta() );
	//default diversity to all 20 aa's if specified filter_delta but did not spec diversify_lvl
	if( filter_delta() != Real( 0 ) && diversify_lvl() == Size( 1 ) ) diversify_lvl( 20 );
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
