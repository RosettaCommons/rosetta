// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Chris King (chrisk1@uw.edu)
// THIS CODE IS DEPRECATED AND WILL SOON DISAPPEAR
// currently this is just a wrapper for GreedyOptMutationMover
//#include <algorithm >
#include <protocols/design_opt/PointMutationCalculator.hh>
#include <protocols/design_opt/GreedyOptMutationMover.hh>
#include <protocols/design_opt/ParetoOptMutationMover.hh>
#include <protocols/design_opt/ParetoOptMutationMoverCreator.hh>
#include <protocols/simple_filters/DeltaFilter.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/pose/PDBInfo.hh>
#include <fstream>
// AUTO-REMOVED #include <utility/file/FileName.hh>
#include <iostream>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
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

namespace protocols {
namespace design_opt {

static thread_local basic::Tracer TR( "protocols.design_opt.ParetoOptMutationMover" );
using namespace core;
using namespace chemical;
using utility::vector1;
using std::pair;

///@brief default ctor
ParetoOptMutationMover::ParetoOptMutationMover() :
	Mover( ParetoOptMutationMoverCreator::mover_name() ),
	task_factory_( /* NULL */ ),
//	filters_( NULL ), /* how set default vecgtor of NULLs? */
//	sample_type_( "low" ),
	scorefxn_( /* NULL */ ),
	relax_mover_( /* NULL */ ),
	dump_pdb_( false ),
	dump_table_( false ),
	parallel_( false ),
	stopping_condition_( /* NULL */ ),
  stop_before_condition_( false ),
  skip_best_check_( false ),
  rtmin_( false ),
  shuffle_order_( false )
{}

//full ctor
ParetoOptMutationMover::ParetoOptMutationMover(
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
	protocols::filters::FilterOP stopping_condition
) :
	Mover( ParetoOptMutationMoverCreator::mover_name() )
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
	stopping_condition_ = stopping_condition;
}

//destruction!
ParetoOptMutationMover::~ParetoOptMutationMover(){}

//creators
protocols::moves::MoverOP
ParetoOptMutationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ParetoOptMutationMover );
}

protocols::moves::MoverOP
ParetoOptMutationMover::clone() const{
	return protocols::moves::MoverOP( new ParetoOptMutationMover( *this ) );
}

//name getters
std::string
ParetoOptMutationMoverCreator::keyname() const
{
	return ParetoOptMutationMoverCreator::mover_name();
}

std::string
ParetoOptMutationMoverCreator::mover_name()
{
	return "ParetoOptMutationMover";
}

std::string
ParetoOptMutationMover::get_name() const {
  return ParetoOptMutationMoverCreator::mover_name();
}

// setter - getter pairs
void
ParetoOptMutationMover::relax_mover( protocols::moves::MoverOP mover ){
	relax_mover_ = mover;
}

protocols::moves::MoverOP
ParetoOptMutationMover::relax_mover() const{
	return relax_mover_;
}

void
ParetoOptMutationMover::filters( vector1< protocols::filters::FilterOP > filters ){
	filters_ = filters;
}

vector1< protocols::filters::FilterOP > ParetoOptMutationMover::filters() const{
	return filters_;
}

void
ParetoOptMutationMover::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

core::pack::task::TaskFactoryOP
ParetoOptMutationMover::task_factory() const
{
	return task_factory_;
}

void
ParetoOptMutationMover::dump_pdb( bool const dump_pdb ){
  dump_pdb_ = dump_pdb;
}

bool
ParetoOptMutationMover::dump_pdb() const{
  return dump_pdb_;
}

void
ParetoOptMutationMover::dump_table( bool const dump_table ){
  dump_table_ = dump_table;
}

bool
ParetoOptMutationMover::dump_table() const{
  return dump_table_;
}

void
ParetoOptMutationMover::stop_before_condition( bool const stop_before_condition ){
  stop_before_condition_ = stop_before_condition;
}

bool
ParetoOptMutationMover::stop_before_condition() const{
  return stop_before_condition_;
}

void
ParetoOptMutationMover::skip_best_check( bool const skip_best_check ){
  skip_best_check_ = skip_best_check;
}

bool
ParetoOptMutationMover::skip_best_check() const{
  return skip_best_check_;
}

void
ParetoOptMutationMover::rtmin( bool const b ){
  rtmin_ = b;
}

bool
ParetoOptMutationMover::rtmin() const{ return rtmin_; }

void
ParetoOptMutationMover::shuffle_order( bool const b ){
  shuffle_order_ = b;
}

bool
ParetoOptMutationMover::shuffle_order() const{ return shuffle_order_; }

utility::vector1< protocols::simple_filters::DeltaFilterOP > 
ParetoOptMutationMover::delta_filters() const { 
  return reset_delta_filters_; 
}

void
ParetoOptMutationMover::delta_filters( utility::vector1< protocols::simple_filters::DeltaFilterOP > const d ){
  reset_delta_filters_ = d;
}

void
ParetoOptMutationMover::sample_types( vector1< std::string > const sample_types ){
  sample_types_ = sample_types;
}

vector1< std::string >
ParetoOptMutationMover::sample_types() const{
  return sample_types_;
}

void
ParetoOptMutationMover::filter_deltas( vector1< core::Real > const filter_deltas ){
  filter_deltas_ = filter_deltas;
}

vector1< core::Real >
ParetoOptMutationMover::filter_deltas() const{
  return filter_deltas_;
}

void
ParetoOptMutationMover::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionOP
ParetoOptMutationMover::scorefxn() const{
	return scorefxn_;
}

void
ParetoOptMutationMover::parallel( bool const parallel ){
  parallel_ = parallel;
}

bool
ParetoOptMutationMover::parallel() const{
  return parallel_;
}


void
ParetoOptMutationMover::apply( core::pose::Pose & pose )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	//this is now but a hollow shell of a mover
	design_opt::GreedyOptMutationMoverOP greedy_opt( new design_opt::GreedyOptMutationMover(
				task_factory(), scorefxn(), relax_mover(), filters(), sample_types(), filter_deltas(), dump_pdb(), dump_table(),
						parallel(), stop_before_condition(), skip_best_check(), rtmin(), shuffle_order(), false /* diversify */, false /* incl_nonopt */, stopping_condition() ) );
	greedy_opt->apply( pose );
}


void
ParetoOptMutationMover::add_filter( protocols::filters::FilterOP filter, std::string const sample_type, core::Real filter_delta )
{
	//filter_delta should always be a scalar!
	if( filter_delta < Real( 0. ) ) filter_delta = -1 * filter_delta;
  filters_.push_back( filter );
  sample_types_.push_back( sample_type );
  filter_deltas_.push_back( filter_delta );
}

//parse rosetta scripts tags
void
ParetoOptMutationMover::parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & )
{
	TR << "ParetoOptMutationMover"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	//load relax mover
	std::string const relax_mover_name( tag->getOption< std::string >( "relax_mover", "null" ) );
	protocols::moves::Movers_map::const_iterator mover_it( movers.find( relax_mover_name ) );
	if( mover_it == movers.end() )
		throw utility::excn::EXCN_RosettaScriptsOption( "Relax mover "+relax_mover_name+" not found" );
	relax_mover( mover_it->second );
	//load scorefxn
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	//load dump_pdb
	dump_pdb( tag->getOption< bool >( "dump_pdb", false ) );
	//load dump_table
	dump_table( tag->getOption< bool >( "dump_table", false ) );
	parallel( tag->getOption< bool >( "parallel", false ) );
	if( tag->hasOption( "stopping_condition" ) ){
		std::string const stopping_filter_name( tag->getOption< std::string >( "stopping_condition" ) );
		stopping_condition( protocols::rosetta_scripts::parse_filter( stopping_filter_name, filters ) );
		TR<<"Defined stopping condition "<<stopping_filter_name<<std::endl;
	}

	//load multiple filters from branch tags
  utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
  BOOST_FOREACH( utility::tag::TagCOP const btag, branch_tags ){
    if( btag->getName() == "Filters" ){
      utility::vector1< utility::tag::TagCOP > const filters_tags( btag->getTags() );
      BOOST_FOREACH( utility::tag::TagCOP const ftag, filters_tags ){
        std::string const filter_name( ftag->getOption< std::string >( "filter_name" ) );
        Filters_map::const_iterator find_filt( filters.find( filter_name ));
        if( find_filt == filters.end() ) {
          TR.Error << "Error !! filter not found in map: \n" << tag << std::endl;
          runtime_assert( find_filt != filters.end() );
        }
        std::string const samp_type( ftag->getOption< std::string >( "sample_type", "low" ));
				core::Real filter_delta( tag->getOption< core::Real >( "filter_delta", core::Real( 0. ) ) );
				add_filter( find_filt->second->clone(), samp_type, filter_delta );
      } //foreach ftag
    }// fi Filters
    else
      throw utility::excn::EXCN_RosettaScriptsOption( "tag name " + btag->getName() + " unrecognized." );
  }//foreach btag
	//load single filter
	{
		std::string const filter_name( tag->getOption< std::string >( "filter", "true_filter" ) );
		if( filter_name != "true_filter" || filters_.size() < 1 ){
			protocols::filters::Filters_map::const_iterator find_filt( filters.find( filter_name ) );
			if( find_filt == filters.end() )
				throw utility::excn::EXCN_RosettaScriptsOption( "Filter "+filter_name+" not found" );
			std::string const samp_type( tag->getOption< std::string >( "sample_type", "low" ) );
			core::Real filter_delta( tag->getOption< core::Real >( "filter_delta", core::Real( 0. ) ) );
			//only add the default dummy filter if we dont have any others, allows user to define filters in branch tags only
			add_filter( find_filt->second->clone(), samp_type, filter_delta );
		}
	}

  //get filters to reset each time a mutation is accepted. For instance, reset the baseline value of delta filters to be the best pose.
  utility::vector1< std::string > delta_filter_names;
  delta_filter_names.clear();
  if( tag->hasOption( "reset_delta_filters" ) ){
    delta_filter_names = utility::string_split( tag->getOption< std::string >( "reset_delta_filters" ), ',' );
    BOOST_FOREACH( std::string const fname, delta_filter_names ){
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


}


} // moves
} // protocols
