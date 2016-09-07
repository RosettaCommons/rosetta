// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Chris King (chrisk1@uw.edu)
//#include <algorithm >
#include <protocols/design_opt/PointMutationCalculator.hh>
#include <protocols/design_opt/GreedyOptMutationMover.hh>
#include <protocols/design_opt/GreedyOptMutationMoverCreator.hh>
#include <protocols/simple_filters/DeltaFilter.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/pose/PDBInfo.hh>
#include <fstream>
#include <iostream>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/util.hh>
#include <boost/foreach.hpp>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

namespace protocols {
namespace design_opt {

static THREAD_LOCAL basic::Tracer TR( "protocols.design_opt.GreedyOptMutationMover" );
using namespace core;
using namespace chemical;
using utility::vector1;
using std::pair;

/// @brief default ctor
GreedyOptMutationMover::GreedyOptMutationMover() :
	Mover( GreedyOptMutationMoverCreator::mover_name() ),
	task_factory_( /* NULL */ ),
	// filters_( NULL ), /* how set default vecgtor of NULLs? */
	// sample_type_( "low" ),
	scorefxn_( /* NULL */ ),
	relax_mover_( /* NULL */ ),
	dump_pdb_( false ),
	dump_table_( false ),
	parallel_( false ),
	stopping_condition_( /* NULL */ ),
	nstruct_iter_( 1 ),
	stop_before_condition_( false ),
	skip_best_check_( false ),
	rtmin_( false ),
	shuffle_order_( false ),
	diversify_( false ),
	incl_nonopt_( false )
{}

//full ctor
GreedyOptMutationMover::GreedyOptMutationMover(
	core::pack::task::TaskFactoryOP task_factory,
	core::scoring::ScoreFunctionOP scorefxn,
	protocols::moves::MoverOP relax_mover,
	vector1< protocols::filters::FilterOP > filters,
	vector1< std::string > sample_types,
	vector1< core::Real > /*filter_deltas*/,
	bool dump_pdb,
	bool dump_table,
	bool parallel,
	bool stop_before_condition,
	bool skip_best_check,
	bool rtmin,
	bool shuffle_order,
	bool diversify,
	bool incl_nonopt,
	protocols::filters::FilterOP stopping_condition
) :
	Mover( GreedyOptMutationMoverCreator::mover_name() )
{
	task_factory_ = task_factory;
	filters_ = filters;
	relax_mover_ = relax_mover;
	scorefxn_ = scorefxn;
	sample_types_ = sample_types;
	dump_pdb_ = dump_pdb;
	dump_table_ = dump_table;
	parallel_ = parallel;
	stop_before_condition_ = stop_before_condition;
	skip_best_check_ = skip_best_check;
	rtmin_ = rtmin;
	shuffle_order_ = shuffle_order;
	diversify_ = diversify;
	incl_nonopt_ = incl_nonopt;
	stopping_condition_ = stopping_condition;
	nstruct_iter_ = 1;
}

//destruction!
GreedyOptMutationMover::~GreedyOptMutationMover()= default;

//creators
protocols::moves::MoverOP
GreedyOptMutationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new GreedyOptMutationMover );
}

protocols::moves::MoverOP
GreedyOptMutationMover::clone() const{
	return protocols::moves::MoverOP( new GreedyOptMutationMover( *this ) );
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
GreedyOptMutationMover::filters( vector1< protocols::filters::FilterOP > filters ){
	filters_ = filters;
	clear_cached_data();
}

vector1< protocols::filters::FilterOP > GreedyOptMutationMover::filters() const{
	return filters_;
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
GreedyOptMutationMover::dump_table( bool const dump_table ){
	dump_table_ = dump_table;
	clear_cached_data();
}

bool
GreedyOptMutationMover::dump_table() const{
	return dump_table_;
}

void
GreedyOptMutationMover::stop_before_condition( bool const stop_before_condition ){
	stop_before_condition_ = stop_before_condition;
}

bool
GreedyOptMutationMover::stop_before_condition() const{
	return stop_before_condition_;
}

void
GreedyOptMutationMover::skip_best_check( bool const skip_best_check ){
	skip_best_check_ = skip_best_check;
}

bool
GreedyOptMutationMover::skip_best_check() const{
	return skip_best_check_;
}

void
GreedyOptMutationMover::rtmin( bool const b ){
	rtmin_ = b;
}

bool
GreedyOptMutationMover::rtmin() const{ return rtmin_; }

void
GreedyOptMutationMover::shuffle_order( bool const b ){
	shuffle_order_ = b;
}

bool
GreedyOptMutationMover::shuffle_order() const{ return shuffle_order_; }

void
GreedyOptMutationMover::diversify( bool const b ){
	diversify_ = b;
}

bool
GreedyOptMutationMover::diversify() const{ return diversify_; }

void
GreedyOptMutationMover::incl_nonopt( bool const b ){
	incl_nonopt_ = b;
}

bool
GreedyOptMutationMover::incl_nonopt() const{ return incl_nonopt_; }

utility::vector1< protocols::simple_filters::DeltaFilterOP >
GreedyOptMutationMover::delta_filters() const {
	return reset_delta_filters_;
}

void
GreedyOptMutationMover::delta_filters( utility::vector1< protocols::simple_filters::DeltaFilterOP > const d ){
	reset_delta_filters_ = d;
}

void
GreedyOptMutationMover::sample_types( vector1< std::string > const sample_types ){
	sample_types_ = sample_types;
	clear_cached_data();
}

vector1< std::string >
GreedyOptMutationMover::sample_types() const{
	return sample_types_;
}

void
GreedyOptMutationMover::filter_deltas( vector1< core::Real > const filter_deltas ){
	filter_deltas_ = filter_deltas;
	clear_cached_data();
}

vector1< core::Real >
GreedyOptMutationMover::filter_deltas() const{
	return filter_deltas_;
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

void
GreedyOptMutationMover::parallel( bool const parallel ){
	parallel_ = parallel;
}

bool
GreedyOptMutationMover::parallel() const{
	return parallel_;
}

//utility funxns for comparing values in sort
bool
cmp_pair_by_second(
	pair< Size, Real > const pair1,
	pair< Size, Real > const pair2 )
{
	return pair1.second < pair2.second;
}

bool
cmp_pair_by_first_vec_val(
	pair< AA, vector1< Real > > const pair1,
	pair< AA, vector1< Real > > const pair2 )
{
	return pair1.second[ 1 ] < pair2.second[ 1 ];
}

bool
cmp_pair_vec_by_first_vec_val(
	pair< Size, vector1< pair< AA, vector1< Real > > > > const pair1,
	pair< Size, vector1< pair< AA, vector1< Real > > > > const pair2 )
{
	return pair1.second[ 1 ].second[ 1 ] < pair2.second[ 1 ].second[ 1 ];
}

void
GreedyOptMutationMover::downsample_pfront_poses(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//sample nstruct numbers from *poses_.size() member data vectors and reassign
	Size nstruct( 1 );
	if ( option[ run::shuffle ]() ) {
		nstruct = option[ out::shuffle_nstruct ]();
	} else {
		nstruct= option[ out::nstruct ]();
	}
	// skip if want more structs than exist poses
	if ( nstruct < pfront_poses_.size() ) {
		vector1< Size > ds_idxs;
		for ( Size i=1; i <= pfront_poses_.size(); ++i ) ds_idxs.push_back( i );
		numeric::random::random_permutation( ds_idxs, numeric::random::rg() );
		ds_idxs.resize( nstruct );
		vector1< pose::Pose > ds_poses;
		vector1< vector1< core::Real > > ds_filter_vals;
		for ( Size i = 1; i <= ds_idxs.size(); ++i ) {
			Size ds_idx( ds_idxs[ i ] );
			ds_poses.push_back( pfront_poses_[ ds_idx ] );
			ds_filter_vals.push_back( pfront_poses_filter_vals_[ ds_idx ] );
		}
		pfront_poses_ = ds_poses;
		pfront_poses_filter_vals_ = ds_filter_vals;
	}
}

void
GreedyOptMutationMover::clear_cached_data(){
	seqpos_aa_vals_vec_.clear();
	pfront_poses_.clear();
	pfront_poses_filter_vals_.clear();
	pfront_poses_filter_ranks_.clear();
	ref_pose_.clear();
}

//TODO: this should also compare fold trees
bool
GreedyOptMutationMover::pose_coords_are_same(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2
){
	//first check for all restype match, also checks same number res
	if ( !pose1.conformation().sequence_matches( pose2.conformation() ) ) return false;
	//then check for all coords identical
	for ( Size i = 1; i <= pose1.total_residue(); ++i ) {
		core::conformation::Residue const & rsd1( pose1.residue( i ) );
		core::conformation::Residue const & rsd2( pose2.residue( i ) );
		//check same n atoms
		if ( rsd1.natoms() != rsd2.natoms() ) return false;
		//and coords
		for ( Size ii = 1; ii <= rsd1.natoms(); ++ii ) {
			if ( rsd1.xyz( ii ).x() != rsd2.xyz( ii ).x() ) return false;
			if ( rsd1.xyz( ii ).y() != rsd2.xyz( ii ).y() ) return false;
			if ( rsd1.xyz( ii ).z() != rsd2.xyz( ii ).z() ) return false;
		}
	}
	return true;
}

void
GreedyOptMutationMover::calc_pfront_poses_filter_ranks(){
	//(re)init ranks w/ bogus zero data
	pfront_poses_filter_ranks_ = vector1< vector1< Size > >( pfront_poses_filter_vals_.size(),
		vector1< Size >( pfront_poses_filter_vals_[ 1 ].size(), Size( 0 ) ) );
	//for each filter type
	for ( Size ifilt = 1; ifilt <= pfront_poses_filter_vals_[ 1 ].size(); ++ifilt ) {
		//copy all this filter vals into a vector of pair< index, val >
		vector1< pair< Size, Real > > filter_index_vals;
		for ( Size ipose = 1; ipose <= pfront_poses_filter_vals_.size(); ++ipose ) {
			filter_index_vals.push_back( pair< Size, Real >( ipose, pfront_poses_filter_vals_[ ipose ][ ifilt ] ) );
		}
		//and sort by value
		std::sort( filter_index_vals.begin(), filter_index_vals.end(), cmp_pair_by_second );
		//now can get rank and index
		for ( Size rank = 1; rank <= filter_index_vals.size(); ++rank ) {
			Size ipose( filter_index_vals[ rank ].first );
			pfront_poses_filter_ranks_[ ipose ][ ifilt ] = rank;
		}
	}
}

void
calc_pareto_front_nbrs(
	vector1< vector1< Real > > coords,
	vector1< bool > const & is_pfront,
	vector1< bool > & is_pfront_nbr,
	vector1< Real > const & nbr_dist
){
	assert( coords.size() == is_pfront.size() );
	//n coords, d dimensions
	Size n( coords.size() );
	Size d( coords[ 1 ].size() );
	assert( nbr_dist.size() == d );
	//for non pfront points
	for ( Size inp = 1; inp <= n; ++inp ) {
		//init possible nbr to false
		is_pfront_nbr[ inp ] = false;
		//skip front pts
		if ( is_pfront[ inp ] ) continue;
		//over pfront pts
		for ( Size ip = 1; ip <= n; ++ip ) {
			if ( ip == inp ) continue;
			//skip non-front pts
			if ( !is_pfront[ ip ] ) continue;
			//now see if ip pfront pt is in inp's nbr box
			//init true, reset to false if fails any dimension
			bool this_is_nbr( true );
			for ( Size k = 1; k <= d; ++k ) {
				if ( std::abs( coords[ ip ][ k ] - coords[ inp ][ k ] ) > nbr_dist[ k ] ) {
					this_is_nbr = false;
					break;
				}
			}
			if ( this_is_nbr ) {
				is_pfront_nbr[ inp ] = true;
				break;
			}
		}
	}
}

//Yes, I'm aware this is the slowest and simplest implementation of pareto front
// identification, but this needs to be finished 2 days ago
// and this calc is hardly the rate limiting step with a bunch of repacking and filter evals
//calc_single_pos_pareto_front
//takes a set (vec1) of coords
//populates boolean vector w/ is_pareto?
/*
instead of just jiggling all the points, we really want to return is_pareto? and is_pareto_nbr( filter_delta )
maybe just another loop over each pareto point to find all pts w/in delta that re not already pareto
then later, have an option for incl_nbrs at each call to clc_pareto_front; e.g. for pfront poses can say no nbrs to still get one lowest pose after trying multiple muts per position
*/
void
calc_pareto_front(
	vector1< vector1< Real > > coords,
	vector1< bool > & is_pfront,
	vector1< Real > const & coord_perts,
	bool const div,
	bool const incl_nbrs
){
	assert( coords.size() == is_pfront.size() );
	//n coords, d dimensions
	Size n( coords.size() );
	Size d( coords[ 1 ].size() );
	assert( coord_perts.size() == d );
	//randomly offset values by jiggle-factors defined in filter_deltas?
	if ( div ) {
		for ( Size i = 1; i <= n; ++i ) {
			for ( Size k = 1; k <= d; ++k ) {
				coords[ i ][ k ] += ( ( numeric::random::rg().uniform() - 0.5 ) * coord_perts[ k ] );
			}
		}
	}
	//for each point
	for ( Size i = 1; i <= n; ++i ) {
		bool is_dom( false );
		//check if strictly dominated by another point
		for ( Size j = 1; j <= n; ++j ) {
			if ( j == i ) continue;
			is_dom = true;
			bool is_equal = true;
			for ( Size k = 1; k <= d; ++k ) {
				//check for same coords
				is_equal = is_equal && ( coords[ i ][ k ] == coords[ j ][ k ] );
				//i not dominated by j if less in any dimension
				if ( coords[ i ][ k ] < coords[ j ][ k ] ) {
					is_dom = false;
					break;
				}
			}
			if ( is_equal ) is_dom = false;
			if ( is_dom ) break;
		}
		if ( is_dom ) is_pfront[ i ] = false;
		else is_pfront[ i ] = true;
	}
	if ( incl_nbrs ) {
		//now get pfront nbrs based on filter_deltas
		vector1< bool > is_pfront_nbr( is_pfront.size(), false );
		calc_pareto_front_nbrs( coords, is_pfront, is_pfront_nbr, coord_perts );
		//reset is_pfront to incl nbrs
		for ( Size ip = 1; ip <= is_pfront.size(); ip++ ) {
			if ( !is_pfront[ ip ] && is_pfront_nbr[ ip ] ) is_pfront[ ip ]  = true;
		}
	}
}

//removes seqpos aa/vals from set that are not pareto optimal
void
GreedyOptMutationMover::filter_seqpos_pareto_opt_ptmuts(){
	for ( Size iseq = 1; iseq <= seqpos_aa_vals_vec_.size(); ++iseq ) {
		vector1< bool > is_pfront( seqpos_aa_vals_vec_[ iseq ].second.size(), false );
		vector1< vector1< Real > > vals;
		//get a vec of vecs from the data
		for ( Size i = 1; i <= seqpos_aa_vals_vec_[ iseq ].second.size(); ++i ) {
			vals.push_back( seqpos_aa_vals_vec_[ iseq ].second[ i ].second );
		}
		//get is_pareto bool vec
		calc_pareto_front( vals, is_pfront, filter_deltas_, diversify_, incl_nonopt_ );
		//replace aa/vals vector with pareto opt only
		vector1< pair< AA, vector1< Real > > > pfront_aa_vals;
		for ( Size iaa = 1; iaa <= seqpos_aa_vals_vec_[ iseq ].second.size(); ++iaa ) {
			if ( is_pfront[ iaa ] ) pfront_aa_vals.push_back( seqpos_aa_vals_vec_[ iseq ].second[ iaa ] );
			//   else if( is_pfront_nbr[ iaa ] ) pfront_aa_vals.push_back( seqpos_aa_vals_vec_[ iseq ].second[ iaa ] );
		}
		seqpos_aa_vals_vec_[ iseq ].second = pfront_aa_vals;
	}
}

//filters a pose vec for pareto opt only
void
GreedyOptMutationMover::filter_pareto_opt_poses(){
	//get is_pareto bool vec
	vector1< bool > is_pfront( pfront_poses_.size(), false );
	//TODO: this implementation reverts to near-deterministic behavior even when filter_delta is set >0,
	// add bool option deversify/perturb vals?
	vector1< Real > null_filter_deltas( filter_deltas_.size(), Real( 0. ) );
	calc_pareto_front( pfront_poses_filter_vals_, is_pfront, filter_deltas_, diversify_, false );
	//remove entries that are not pareto
	vector1< pose::Pose > new_pfront_poses;
	vector1< vector1< Real > > new_pfront_poses_filter_vals;
	for ( Size ipose = 1; ipose <= pfront_poses_.size(); ++ipose ) {
		if ( is_pfront[ ipose ] ) {
			new_pfront_poses.push_back( pfront_poses_[ ipose ] );
			new_pfront_poses_filter_vals.push_back( pfront_poses_filter_vals_[ ipose ] );
		}
	}
	pfront_poses_ = new_pfront_poses;
	pfront_poses_filter_vals_ = new_pfront_poses_filter_vals;
}

//this should be in PointMutCalc
void
GreedyOptMutationMover::dump_scoring_table( std::string filename, core::pose::Pose const & ref_pose ) const{
	utility::io::ozstream outtable(filename, std::ios::out | std::ios::app ); // Append if logfile already exists.
	if ( outtable ) {
		for ( core::Size ii(1); ii <= seqpos_aa_vals_vec_.size(); ++ii ) {
			core::Size pos( seqpos_aa_vals_vec_[ii].first );
			utility::vector1< std::pair< core::chemical::AA, utility::vector1< core::Real > > > const & aa_pairs( seqpos_aa_vals_vec_[ii].second );
			outtable << pos ;
			if ( ref_pose.pdb_info() ) {
				outtable << " (" << ref_pose.pdb_info()->pose2pdb(pos) << ")";
			}
			outtable << '\t';
			for ( core::Size jj(1); jj <= aa_pairs.size(); ++jj ) {
				outtable << aa_pairs[jj].first << ((ref_pose.aa(pos) == aa_pairs[jj].first)?"*:":":");
				for ( core::Size kk( 1 ); kk <= aa_pairs[ jj ].second.size(); ++kk ) {
					outtable << aa_pairs[jj].second[ kk ] << ":";
				}
				outtable << " ";
			}
			outtable << std::endl;
		}
		outtable << std::endl; // Blank line at end to seperate.
	} else {
		TR.Warning << "WARNING: Unable to open file " << filename << " for writing GreedyOptMutationMover table output." << std::endl;
	}
	outtable.close();
}


void
GreedyOptMutationMover::apply( core::pose::Pose & pose )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	//store input pose
	core::pose::Pose start_pose( pose );
	design_opt::PointMutationCalculatorOP ptmut_calc( new design_opt::PointMutationCalculator(
		task_factory(), scorefxn(), relax_mover(), filters(), sample_types(), dump_pdb(), false, parallel() ) );

	//create vec of pairs of seqpos, vector of AA/val pairs that pass input filter
	//then combine them into a pareto opt pose set
	//only calc the ptmut data and pareto set once per pose, not at every nstruct iteration
	//but how will we know if that data is still valid? what if the pose has chnged?
	//the best answer is to store the pose passed to apply in a private variable (ref_pose_)
	//and only calc if ref_pose_ is still undef or doesnt match apply pose
	//also recalc if pareto pose set is empty
	if ( pfront_poses_.empty() || ref_pose_.empty() || !pose_coords_are_same( start_pose, ref_pose_ ) ) {
		//reset our private data
		clear_cached_data();
		//and (re)set ref_pose_ to this pose
		ref_pose_ = start_pose;

		//get the point mut values
		ptmut_calc->calc_point_mut_filters( start_pose, seqpos_aa_vals_vec_ );
		if ( seqpos_aa_vals_vec_.size() < 1 ) {
			utility_exit_with_message( "ERROR: No acceptable mutations found. All possible mutations failed at least one filter!" );
		}
		//this part sorts the seqpos/aa/val data so that we init with something good (1st)
		//first over each seqpos by aa val, then over all seqpos by best aa val
		for ( Size ivec = 1; ivec <= seqpos_aa_vals_vec_.size(); ++ivec ) {
			//skip if aa/vals vector is empty
			if ( seqpos_aa_vals_vec_[ ivec ].second.empty() ) continue;
			//sort aa/vals in incr val order
			std::sort( seqpos_aa_vals_vec_[ ivec ].second.begin(),
				seqpos_aa_vals_vec_[ ivec ].second.end(), cmp_pair_by_first_vec_val );
		}
		//now sort seqpos_aa_vals_vec_ by *first* (lowest) val in each seqpos vector, low to high
		//uses cmp_pair_vec_by_first_vec_val to sort based on second val in
		//first pair element of vector in pair( size, vec( pair ) )
		std::sort( seqpos_aa_vals_vec_.begin(), seqpos_aa_vals_vec_.end(), cmp_pair_vec_by_first_vec_val );

		//finally, dump table to file, if requested.
		if ( dump_table() ) {
			std::string fname( "GreedyOptTable" );
			if ( protocols::jd2::jd2_used() ) {
				fname += "_" + protocols::jd2::current_output_name();
			}
			fname += ".tab";
			dump_scoring_table( fname, start_pose );
		}

		//this part gets rid of ptmuts that are not pareto opt
		filter_seqpos_pareto_opt_ptmuts();

		//now randomize the sequence position order?
		if ( shuffle_order() ) {
			numeric::random::random_permutation( seqpos_aa_vals_vec_.begin(), seqpos_aa_vals_vec_.end(), numeric::random::rg() );
			TR<<"Combining shuffled mutations… " << std::endl;
		} else TR<<"Combining sorted mutations… " << std::endl;

		//init pareto opt poses with first mutations for now
		//TODO: is there a better way to init?
		Size iseq_init( 1 );
		//the resi index is the first part of the pair
		Size resi_init( seqpos_aa_vals_vec_[ iseq_init ].first );
		for ( Size iaa = 1; iaa <= seqpos_aa_vals_vec_[ iseq_init ].second.size(); ++iaa ) {
			AA target_aa( seqpos_aa_vals_vec_[ iseq_init ].second[ iaa ].first );
			pose::Pose new_pose( start_pose );
			ptmut_calc->mutate_and_relax( new_pose, resi_init, target_aa );
			//this is where the pose gets saved
			pfront_poses_.push_back( new_pose );
			vector1< Real > vals;
			bool filter_pass;
			ptmut_calc->eval_filters( new_pose, filter_pass, vals );
			pfront_poses_filter_vals_.push_back( vals );
		}

		//now try to combine pareto opt mutations
		for ( Size iseq = 2; iseq <= seqpos_aa_vals_vec_.size(); ++iseq ) {
			//create new pose vec to hold all combinations
			vector1< pose::Pose > new_poses;
			vector1< vector1< Real > > new_poses_filter_vals;
			//the resi index is the first part of the pair
			Size seqpos( seqpos_aa_vals_vec_[ iseq ].first );
			TR << "Combining " << pfront_poses_.size() << " structures with mutations at residue " << seqpos << std::endl;
			//over each current pfront pose
			//pfront_poses_ contains all the current poses
			for ( Size ipose = 1; ipose <= pfront_poses_.size(); ++ipose ) {
				//over all aa's at seqpos
				for ( Size iaa = 1; iaa <= seqpos_aa_vals_vec_[ iseq ].second.size(); ++iaa ) {
					//inside this double loop, all pfront_poses_ are combined with all muts at this position
					AA target_aa( seqpos_aa_vals_vec_[ iseq ].second[ iaa ].first );
					pose::Pose new_pose( pfront_poses_[ ipose ] );

					bool filter_pass;
					vector1< Real > vals;
					ptmut_calc->mutate_and_relax( new_pose, seqpos, target_aa );
					ptmut_calc->eval_filters( new_pose, filter_pass, vals );
					//only check this guy for pareto if passes
					if ( !filter_pass ) continue;

					//this is where the pose gets saved
					new_poses.push_back( new_pose );
					new_poses_filter_vals.push_back( vals );
				}
			}
			//   TR << "Generated " << new_poses.size() << " new poses from mutations at residue "
			//     << seqpos << ". Filtering... " << std::endl;

			//only update pfront poses if we found any new ones, else just skip this position,
			// because we know our current pose does pass all the filters
			if ( new_poses.size() < 1 ) {
				TR << "Unable to generate any new poses that pass all filters at position " << iseq << std::endl;
				continue;
			}

			//end optimization if *any* of the new pareto poses trips the stopping_condition
			bool stop( false );
			for ( Size ipose = 1; ipose <= new_poses.size(); ++ipose ) {
				pose::Pose new_pose( new_poses[ ipose ] );
				if ( stopping_condition() && stopping_condition()->apply( new_pose ) ) {
					stop = true;
					if ( !stop_before_condition() ) {
						TR<<"Stopping condition evaluates to true. Stopping early and acceptingn the last mutation: "<<
							start_pose.residue( seqpos ).name1() << seqpos << new_pose.residue( seqpos ).name1() << std::endl;
					} else {
						TR<<"Stopping condition evaluates to true. Stopping early and rejecting the last mutation."
							<< start_pose.residue( seqpos ).name1() << seqpos << new_pose.residue( seqpos ).name1() << std::endl;
					}
					break;
				}
			}

			if ( stop && stop_before_condition() ) break;
			//if forcing new mutations, then reset pose cache and filter just the new guys
			if ( skip_best_check() ) {
				pfront_poses_ = new_poses;
				pfront_poses_filter_vals_ = new_poses_filter_vals;
			} else {
				//default: add the new poses to the pose cache then pareto filter
				for ( Size ipose = 1; ipose <= new_poses.size(); ++ipose ) {
					pfront_poses_.push_back( new_poses[ ipose ] );
					pfront_poses_filter_vals_.push_back( new_poses_filter_vals[ ipose ] );
				}
			}

			//filter new_poses for the pareto opt set
			assert( pfront_poses_.size() == pfront_poses_filter_vals_.size() );
			filter_pareto_opt_poses();
			assert( pfront_poses_.size() == pfront_poses_filter_vals_.size() );

			//break out if we've reached our stopping condition
			if ( stop ) break;

			//Optionally reset baseline for Delta Filters (useful so that the mutations are still evaluated on an individual basis, in the context of the current best pose).
			BOOST_FOREACH ( protocols::simple_filters::DeltaFilterOP const delta_filter, reset_delta_filters_ ) {
				std::string const fname( delta_filter->get_user_defined_name() );
				core::Real const fbaseline( delta_filter->filter()->report_sm( pose ) );
				delta_filter->baseline( fbaseline );
				delta_filter->ref_baseline( fbaseline );
				TR<<"Reset baseline for DeltaFilter "<<fname<<" to "<<fbaseline<<std::endl;
			}
			TR << "Generated  " << pfront_poses_.size() << " new structures at residue " << seqpos << " :";
			downsample_pfront_poses();
			TR << " Downsampled pareto front to  " << pfront_poses_.size() << " structures" << std::endl;
		}
		//randomly downsample cached parreto front data to -nstruct structures
		downsample_pfront_poses();
		TR << "Generated  " << pfront_poses_.size() << " final structures." << std::endl;
		calc_pfront_poses_filter_ranks();
	}
	//assign the apply pose to the next pfront_pose
	Size pfront_pose_iter( ( nstruct_iter_ - 1 ) % pfront_poses_.size() + 1 );
	pose = pfront_poses_[ pfront_pose_iter ];
	//print out filter vals for this pose
	TR << "Structure " << nstruct_iter_ << " filter values: ";
	for ( Size ival = 1; ival <= pfront_poses_filter_vals_[ pfront_pose_iter ].size(); ++ival ) {
		TR << " " << filters()[ ival ]->get_user_defined_name() << ": "
			<< pfront_poses_filter_vals_[ pfront_pose_iter ][ ival ];
	}
	TR << std::endl;
	TR << "Structure " << nstruct_iter_ << " filter ranks: ";
	for ( Size ival = 1; ival <= pfront_poses_filter_ranks_[ pfront_pose_iter ].size(); ++ival ) {
		TR << " " << filters()[ ival ]->get_user_defined_name() << ": "
			<< pfront_poses_filter_ranks_[ pfront_pose_iter ][ ival ];
	}
	TR << std::endl;
	TR.flush();

	//increment pose iterator to get new pfront_pose at nstruct+1
	nstruct_iter_ += 1;
}


void
GreedyOptMutationMover::add_filter( protocols::filters::FilterOP filter, std::string const sample_type, core::Real filter_delta )
{
	//filter_delta should always be a scalar!
	if ( filter_delta < Real( 0. ) ) filter_delta = -1 * filter_delta;
	filters_.push_back( filter );
	sample_types_.push_back( sample_type );
	filter_deltas_.push_back( filter_delta );
}

//parse rosetta scripts tags
void
GreedyOptMutationMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & )
{
	TR << "GreedyOptMutationMover"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	//load relax mover
	std::string const relax_mover_name( tag->getOption< std::string >( "relax_mover", "null" ) );
	auto mover_it( movers.find( relax_mover_name ) );
	if ( mover_it == movers.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Relax mover "+relax_mover_name+" not found" );
	}
	relax_mover( mover_it->second );
	//load scorefxn
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	//load dump_pdb
	dump_pdb( tag->getOption< bool >( "dump_pdb", false ) );
	//load dump_table
	dump_table( tag->getOption< bool >( "dump_table", false ) );
	parallel( tag->getOption< bool >( "parallel", false ) );
	if ( tag->hasOption( "stopping_condition" ) ) {
		std::string const stopping_filter_name( tag->getOption< std::string >( "stopping_condition" ) );
		stopping_condition( protocols::rosetta_scripts::parse_filter( stopping_filter_name, filters ) );
		TR<<"Defined stopping condition "<<stopping_filter_name<<std::endl;
	}

	//load multiple filters from branch tags
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	BOOST_FOREACH ( utility::tag::TagCOP const btag, branch_tags ) {
		if ( btag->getName() == "Filters" ) {
			utility::vector1< utility::tag::TagCOP > const filters_tags( btag->getTags() );
			BOOST_FOREACH ( utility::tag::TagCOP const ftag, filters_tags ) {
				std::string const filter_name( ftag->getOption< std::string >( "filter_name" ) );
				auto find_filt( filters.find( filter_name ));
				if ( find_filt == filters.end() ) {
					TR.Error << "Error !! filter not found in map: \n" << tag << std::endl;
					runtime_assert( find_filt != filters.end() );
				}
				std::string const samp_type( ftag->getOption< std::string >( "sample_type", "low" ));
				core::Real filter_delta( tag->getOption< core::Real >( "filter_delta", core::Real( 0. ) ) );
				add_filter( find_filt->second->clone(), samp_type, filter_delta );
			} //foreach ftag
		} else { // fi Filters
			throw utility::excn::EXCN_RosettaScriptsOption( "tag name " + btag->getName() + " unrecognized." );
		}
	}//foreach btag
	//load single filter
	{
		std::string const filter_name( tag->getOption< std::string >( "filter", "true_filter" ) );
		if ( filter_name != "true_filter" || filters_.size() < 1 ) {
			auto find_filt( filters.find( filter_name ) );
			if ( find_filt == filters.end() ) {
				throw utility::excn::EXCN_RosettaScriptsOption( "Filter "+filter_name+" not found" );
			}
			std::string const samp_type( tag->getOption< std::string >( "sample_type", "low" ) );
			core::Real filter_delta( tag->getOption< core::Real >( "filter_delta", core::Real( 0. ) ) );
			//only add the default dummy filter if we dont have any others, allows user to define filters in branch tags only
			add_filter( find_filt->second->clone(), samp_type, filter_delta );
		}
	}

	//get filters to reset each time a mutation is accepted. For instance, reset the baseline value of delta filters to be the best pose.
	utility::vector1< std::string > delta_filter_names;
	delta_filter_names.clear();
	if ( tag->hasOption( "reset_delta_filters" ) ) {
		delta_filter_names = utility::string_split( tag->getOption< std::string >( "reset_delta_filters" ), ',' );
		BOOST_FOREACH ( std::string const fname, delta_filter_names ) {
			reset_delta_filters_.push_back( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::DeltaFilter > ( protocols::rosetta_scripts::parse_filter( fname, filters ) ) );
			TR<<"The baseline for Delta Filter "<<fname<<" will be reset upon each accepted mutation"<<std::endl;
		}
	}
	//should mutations be allowed around the tested/introduced point mutation, if so, what shell radius
	design_shell_ = tag->getOption< core::Real >( "design_shell", -1.0 );
	//repack which radius after mutating
	repack_shell_ = tag->getOption< core::Real >( "repack_shell", 8.0 );
	//stop mover once the stopping_condition is reached and do not accept the last mutation (ie, reject the mutation that set the stopping_condition to true)
	stop_before_condition( tag->getOption< bool >( "stop_before_condition", false ) );
	//accept mutations during the combining stage as long as they pass the filter(s), regardless of whether or not the value is the best so far.
	skip_best_check( tag->getOption< bool >( "skip_best_check", false ) );
	rtmin( tag->getOption< bool >( "rtmin", false ) );
	//load shuffle_order
	shuffle_order( tag->getOption< bool >( "shuffle_order", false ) );
	diversify( tag->getOption< bool >( "diversify", true ) );
	incl_nonopt( tag->getOption< bool >( "incl_nonopt", false ) );

}


} // moves
} // protocols
