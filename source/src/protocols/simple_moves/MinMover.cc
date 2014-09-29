// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/MinMover.cc
/// @brief  method definitions for MinMover
/// @author ashworth

// Unit headers
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/MinMoverCreator.hh>

// Package headers
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh> // get_score_function
#include <core/pose/Pose.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>

// Basic headers
#include <basic/prof.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>

// Boost Headers
#include <boost/foreach.hpp>
#include <basic/Tracer.hh>


namespace protocols {
namespace simple_moves {

using namespace core;
using namespace kinematics;
using namespace optimization;
using namespace scoring;
using core::pack::task::PackerTaskOP;


static thread_local basic::Tracer TR( "protocols.simple_moves.MinMover" );

std::string
MinMoverCreator::keyname() const
{
	return MinMoverCreator::mover_name();
}

protocols::moves::MoverOP
MinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MinMover );
}

std::string
MinMoverCreator::mover_name()
{
	return "MinMover";
}

// default constructor
// proper lightweight default constructor
MinMover::MinMover() :
    Parent("MinMover"),
		movemap_(/* 0 */),
		scorefxn_(/* 0 */),
		min_options_(/* 0 */),
		cartesian_(false),
		dof_tasks_()
{
	omega_ = true;
	min_options_ = MinimizerOptionsOP( new MinimizerOptions( "linmin", 0.01, true, false, false ) );
}

MinMover::MinMover( std::string const & name ) :
    Parent( name ),
		movemap_(/* 0 */),
		scorefxn_(/* 0 */),
		min_options_(/* 0 */),
		cartesian_(false),
		dof_tasks_()
{
	omega_ = true;
	min_options_ = MinimizerOptionsOP( new MinimizerOptions( "linmin", 0.01, true, false, false ) );
}

MinMover::~MinMover(){}

// constructor with arguments
MinMover::MinMover(
	MoveMapOP movemap_in,
	ScoreFunctionCOP scorefxn_in,
	std::string const & min_type_in,
	Real tolerance_in,
	bool use_nb_list_in,
	bool deriv_check_in /* = false */,
	bool deriv_check_verbose_in /* = false */
) : Parent("MinMover"),
		movemap_( movemap_in ),
		scorefxn_( scorefxn_in ),
		min_options_(/* 0 */),
		threshold_(1000000.0), // TODO: line can be deleted?
		cartesian_(false),
		dof_tasks_()
{
	omega_ = true;
	min_options_ = MinimizerOptionsOP( new MinimizerOptions(
		min_type_in, tolerance_in, use_nb_list_in, deriv_check_in, deriv_check_verbose_in ) );
}

/// @brief allow non-const access to the internal minimizer options object
MinimizerOptionsOP
MinMover::min_options() {
	return min_options_;
}

/// @brief allow const access to the internal minimizer options object
MinimizerOptionsCOP
MinMover::min_options() const {
	return min_options_;
}

void
MinMover::min_options( MinimizerOptionsOP min_options) {
	min_options_ = min_options;
}

void
MinMover::movemap( MoveMapCOP movemap_in )
{
	runtime_assert( movemap_in != 0 );
	movemap_ = core::kinematics::MoveMapOP( new MoveMap( *movemap_in ) );
}

void MinMover::set_movemap( MoveMapCOP movemap_in ){
  movemap( movemap_in );
}

MoveMapCOP
MinMover::movemap() const
{
	return movemap_;
}

void
MinMover::score_function( ScoreFunctionCOP scorefxn_in )
{
	runtime_assert( scorefxn_in != 0 );
	scorefxn_ = scorefxn_in;
}

// Why would you want a deep copy of the ScoreFunction here?
// In any case, commenting this out, because PyRosetta cannot distinguish between these two function signature, always
// chosing this one and leading to unexpected behavior. ~Labonte
/*void
MinMover::score_function( ScoreFunction const & scorefxn_in )
{
	scorefxn_ = scorefxn_in.clone();
}*/

ScoreFunctionCOP
MinMover::score_function() const
{
	return scorefxn_;
}

void MinMover::min_type( std::string min_type_in ) { min_options_->min_type( min_type_in ); }
std::string MinMover::min_type() const { return min_options_->min_type(); }


void MinMover::tolerance( Real tolerance_in ) { min_options_->minimize_tolerance( tolerance_in ); }
Real MinMover::tolerance() const { return min_options_->minimize_tolerance(); }


void MinMover::nb_list( bool nb_list_in ) { min_options_->use_nblist( nb_list_in ); }
bool MinMover::nb_list() const { return min_options_->use_nblist(); }


void MinMover::deriv_check( bool deriv_check_in ) { min_options_->deriv_check( deriv_check_in ); }
bool MinMover::deriv_check() const { return min_options_->deriv_check(); }


///@details restrict the move map by the packer task:
///If a residue is not designable, the backbone is fixes
///If a residue is not packable, the sidechain is fixed
///
///WARNING: This is extending the abuse of using task operations for
///general ResidueSubsets
///
///When TaskOperations replaced with ResidueSetOperations please
///change this too!
void
MinMover::apply_dof_tasks_to_movemap(
	core::pose::Pose const & pose,
	MoveMap & movemap
) const {

	for(
		DOF_TaskMap::const_iterator t=dof_tasks_.begin(), te=dof_tasks_.end();
		t != te; ++t){
		//generate task
		PackerTaskOP task( t->second->create_task_and_apply_taskoperations( pose ) );

		//modify movemap by task
		Size const nres( task->total_residue() );

		for ( Size i(1); i <= nres; ++i ) {
			if ( !task->pack_residue( i ) ){
				if( t->first.first == core::id::PHI ){
					movemap.set( t->first.second, false );
				} else {
					movemap.set( t->first.first, false );
				}
			}
		}
	}
}

void
MinMover::max_iter( Size max_iter_in ) { min_options_->max_iter( max_iter_in ); }

void
MinMover::apply(pose::Pose & pose) {
	// lazy default initialization
	MoveMapOP active_movemap;
	if ( ! movemap() ) movemap() = MoveMapOP( new MoveMap );
	else active_movemap = movemap()->clone();

	apply_dof_tasks_to_movemap(pose, *active_movemap);


	if ( ! scorefxn_ ) scorefxn_ = get_score_function(); // get a default (INITIALIZED!) ScoreFunction

	PROF_START( basic::MINMOVER_APPLY );
	if (!cartesian( )) {
		AtomTreeMinimizer minimizer;
		if( !omega_ ){
			TR<<"shutting off omega dihedral angle minimization"<<std::endl;
			for( core::Size i = 1; i <= pose.total_residue(); ++i ){
			  if( !pose.residue( i ).is_protein() ) continue;
				active_movemap->set( core::id::TorsionID( core::id::omega_torsion, BB, i), false );
		  }
    }
		(*scorefxn_)(pose);
		minimizer.run( pose, *active_movemap, *scorefxn_, *min_options_ );
	} else {
		CartesianMinimizer minimizer;
		(*scorefxn_)(pose);
		minimizer.run( pose, *active_movemap, *scorefxn_, *min_options_ );
	}
	PROF_STOP( basic::MINMOVER_APPLY );

  // emit statistics
  scorefxn_->show(TR.Debug, pose);
  TR.Debug << std::endl;
}

std::string
MinMover::get_name() const {
	return MinMoverCreator::mover_name();
}

void
MinMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Minimization type:\t" << min_type() << "\nScorefunction:\t\t";
	if ( score_function() != 0 ) {
		output  << score_function()->get_name() << std::endl;
	}
	else { output << "none" << std::endl; }
	output << "Score tolerance:\t" << tolerance() << "\nNb list:\t\t" << (nb_list() ? "True" : "False") <<
			"\nDeriv check:\t\t" << (deriv_check() ? "True" : "False") << std::endl << "Movemap:" << std::endl;
	if (movemap() != 0) {
		movemap()->show(output);
	}
}

protocols::moves::MoverOP MinMover::clone() const { return protocols::moves::MoverOP( new protocols::simple_moves::MinMover( *this ) ); }
protocols::moves::MoverOP MinMover::fresh_instance() const { return protocols::moves::MoverOP( new MinMover ); }

void MinMover::parse_def_opts( utility::lua::LuaObject const & def,
	utility::lua::LuaObject const & score_fxns,
	utility::lua::LuaObject const & /*tasks*/,
	protocols::moves::MoverCacheSP /*cache*/ ) {
	if( def["scorefxn"] ) {
		score_function( protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns ) );
	} else {
		score_function( score_fxns["score12"].to<ScoreFunctionSP>()->clone()  );
	}

	if ( ! movemap_ ) movemap_ = core::kinematics::MoveMapOP( new MoveMap );
	if( def["jump"] ) {
		for (utility::lua::LuaIterator i=def["jump"].begin(), end; i != end; ++i) {
			if( (*i).to<int>() == -1 ) {
				movemap_->set_jump( true );
				break;
			} else if( (*i).to<core::Size>() == 0 ) {
				movemap_->set_jump( false );
				break;
			} else {
				TR << "Setting min on jump " << (*i).to<core::Size>() << std::endl;
				movemap_->set_jump( (*i).to<core::Size>(), true );
			}
		}
	}

	min_type( def["type"] ? def[ "type" ].to<std::string>() : "dfpmin_armijo_nonmonotone" );
	tolerance( def["tolerance"] ? def[ "tolerance" ].to<core::Real>() : 0.01 );
	cartesian( def["cartesian"] ? def[ "cartesian" ].to<bool>() : false );
	//fpd if cartesian default to lbfgs minimization
	if ( cartesian() && ! def["type"] ) {
		min_type( "lbfgs_armijo_nonmonotone" );
	}

	// TODO: Add parsing of task operations?
}

void MinMover::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & score_fxns,
				utility::lua::LuaObject const & tasks,
				protocols::moves::MoverCacheSP cache ) {
	if ( ! movemap_ ) movemap_ = core::kinematics::MoveMapOP( new MoveMap );

	parse_def_opts( def, score_fxns, tasks, cache );

	if( def["chi"] )
		movemap_->set_chi(def["chi"].to<bool>());
	if( def["bb"] )
		movemap_->set_chi(def["bb"].to<bool>());


	if( def["movemap"] )
		protocols::elscripts::parse_movemapdef( def["movemap"], movemap_ );
}

void MinMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose )
{
	if ( ! movemap_ ) movemap_ = core::kinematics::MoveMapOP( new MoveMap );
	parse_opts( tag, data, filters, movers, pose );
	parse_chi_and_bb( tag );
	parse_dof_tasks( tag, data );

/// parse_movemap will reset the movemap to minimize all if nothing is stated on it. The following section ensures that the protocol specifies a MoveMap before going that way
//	utility::vector1< TagCOP > const branch_tags( tag->getTags() );
//	utility::vector1< TagCOP >::const_iterator tag_it;
	protocols::rosetta_scripts::parse_movemap( tag, pose, movemap_, data, false );

}

void MinMover::parse_opts(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{
	score_function( protocols::rosetta_scripts::parse_score_function( tag, data ) );

	if ( tag->hasOption("jump") ) {
		if ( ! movemap_ ) movemap_ = core::kinematics::MoveMapOP( new MoveMap );
		utility::vector1<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );
		// string 'ALL' makes all jumps movable
		if (jumps.size() == 1 && (jumps[1] == "ALL" || jumps[1] == "All" || jumps[1] == "all" || jumps[1] == "*") ) {
			movemap_->set_jump( true );
		} else if( tag->getOption< core::Size > ( "jump" ) == 0 ) {
			movemap_->set_jump( false );
		} else {
			BOOST_FOREACH(std::string jump, jumps){
				Size const value = std::atoi( jump.c_str() ); // convert to C string, then convert to integer, then set a Size (phew!)
				TR << "Setting min on jump " << value << std::endl;
				movemap_->set_jump( value, true );
			}
		}
	}
	max_iter( tag->getOption< int >( "max_iter", 200 ) );
	min_type( tag->getOption< std::string >( "type", "dfpmin_armijo_nonmonotone" ) );
	tolerance( tag->getOption< core::Real >( "tolerance", 0.01 ) );
	cartesian( tag->getOption< core::Real >( "cartesian", false ) );

	//fpd  if cartesian or nonideal default to lbfgs minimization otherwise the runtime is horrible
	if ( cartesian() && !tag->hasOption("type") ) {
		min_type( "lbfgs_armijo_nonmonotone" );
	}
	if ( ( tag->getOption<bool>("bondangle",0) ||tag->getOption<bool>("bondlength",0) )
	    && !tag->hasOption("type") ) {
		min_type( "lbfgs_armijo_nonmonotone" );
	}
}

void MinMover::parse_chi_and_bb( TagCOP const tag )
{
	if ( ! movemap_ ) movemap_ = core::kinematics::MoveMapOP( new MoveMap );
	bool const chi( tag->getOption< bool >( "chi" ) ), bb( tag->getOption< bool >( "bb" ) );
  omega_ = tag->getOption< bool >( "omega", true );
	movemap_->set_chi( chi );
	movemap_->set_bb( bb );
	TR<<"Options chi, bb: "<<chi<<", "<<bb<<" omega: "<<omega_<<std::endl;
	if ( tag->hasOption("bondangle") ) {
		bool const value( tag->getOption<bool>("bondangle") );
		movemap_->set( core::id::THETA, value );
	}
	if ( tag->hasOption("bondlength") ) {
		bool const value( tag->getOption<bool>("bondlength") );
		movemap_->set( core::id::D, value );
	}
}

///@detail helper function for parse_of_tasks
void
MinMover::parse_dof_task_type(
	std::string const & tag_name,
	core::id::DOF_Type dof_type,
	core::id::TorsionType torsion_type,
	TagCOP const tag,
	basic::datacache::DataMap & data
) {

	if(!tag->hasOption(tag_name)){
		return;
	}

  using namespace core::pack::task;
  using namespace core::pack::task::operation;
	using std::string;
  typedef utility::vector1< std::string > StringVec;

	TaskFactoryOP task_factory( new TaskFactory() );
	string const t_o_val( tag->getOption<string>(tag_name) );
	StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );

	BOOST_FOREACH( string t_o_key, t_o_keys ){
		if ( data.has( "task_operations", t_o_key ) ) {
			task_factory->push_back( data.get_ptr< TaskOperation >( "task_operations", t_o_key ) );
			TR << "\t " << tag_name << ": " << t_o_key << std::endl;
		} else {
			utility_exit_with_message("TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
		}
	}
	dof_tasks_[	std::make_pair(dof_type, torsion_type)] =	task_factory;
}

void
MinMover::parse_dof_tasks(
	TagCOP const tag,
	basic::datacache::DataMap & data
) {

	if(
		tag->hasOption("bb_task_operations") ||
		tag->hasOption("chi_task_operations") ||
		tag->hasOption("bondangle_task_operations") ||
		tag->hasOption("bondlength_task_operations")){
		TR
			<< "Adding the following task operations to mover " << tag->getName() << " "
			<< "called " << tag->getOption<std::string>( "name", "no_name" ) << ":" << std::endl;
	}

	parse_dof_task_type( "bb_task_operations", core::id::PHI, core::id::BB, tag, data );
	parse_dof_task_type( "chi_task_operations", core::id::PHI, core::id::CHI, tag, data );
	parse_dof_task_type(
		"bondangle_task_operations",
		core::id::THETA,
		core::id::BB, // (dummy parameter)
		tag, data );

	parse_dof_task_type(
		"bondlength_task_operations",
		core::id::D,
		core::id::BB, // (dummy parameter)
		tag, data );
}

std::ostream &operator<< (std::ostream &os, MinMover const &mover)
{
	mover.show(os);
	return os;
}

}  // simple_moves
}  // protocols
