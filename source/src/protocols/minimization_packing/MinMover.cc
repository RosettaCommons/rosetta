// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/minimization_packing/MinMover.cc
/// @brief  method definitions for MinMover
/// @author ashworth

// Unit headers
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/MinMoverCreator.hh>

// Package headers
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/jump_selector/JumpIndexSelector.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh> // get_score_function
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/TaskOperation.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

// Basic headers
#include <basic/prof.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>

// Boost Headers
#include <basic/Tracer.hh>


namespace protocols {
namespace minimization_packing {

using namespace core;
using namespace kinematics;
using namespace optimization;
using namespace scoring;
using core::pack::task::PackerTaskOP;


static basic::Tracer TR( "protocols.minimization_packing.MinMover" );

std::string
MinMoverCreator::keyname() const
{
	return MinMoverCreator::mover_name();
}

protocols::moves::MoverOP
MinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MinMover );
}

void MinMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	MinMover::provide_xml_schema( xsd );
}

std::string
MinMoverCreator::mover_name()
{
	return "MinMover";
}

std::string
MinMover::mover_name()
{
	return MinMoverCreator::mover_name();
}

// default constructor
// proper lightweight default constructor
MinMover::MinMover() :
	Parent("MinMover"),
	movemap_(/* 0 */),
	scorefxn_(/* 0 */),
	min_options_(/* 0 */),
	abs_score_convergence_threshold_(0.0),
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
	abs_score_convergence_threshold_(0.0),
	cartesian_(false),
	dof_tasks_()
{
	omega_ = true;
	min_options_ = MinimizerOptionsOP( new MinimizerOptions( "linmin", 0.01, true, false, false ) );
}

MinMover::~MinMover()= default;

// constructor with arguments
MinMover::MinMover(
	MoveMapOP movemap_in,
	ScoreFunctionCOP scorefxn_in,
	std::string const & min_type_in,
	Real tolerance_in,
	bool use_nb_list_in,
	bool deriv_check_in /* = false */,
	bool deriv_check_verbose_in /* = false */
) :
	Parent("MinMover"),
	movemap_(std::move( movemap_in )),
	scorefxn_(std::move( scorefxn_in )),
	min_options_(/* 0 */),
	abs_score_convergence_threshold_(0.0),
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
	movemap_ = core::kinematics::MoveMapOP( new MoveMap( *movemap_in ) );
}

void MinMover::set_movemap( MoveMapCOP movemap_in ){
	movemap( movemap_in );
}

void
MinMover::movemap_factory( core::select::movemap::MoveMapFactoryCOP mmf ) {
	movemap_factory_ = mmf;
}

MoveMapCOP
MinMover::movemap( core::pose::Pose const & pose ) const {
	if ( movemap_ ) {
		return movemap_;
	} else if ( movemap_factory_ ) {
		return movemap_factory_->create_movemap_from_pose( pose );
	} else {
		return core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	}
}

/// @brief Get the explicitly set movemap object, if present, nullptr if not
/// Will not create a default or look at the MoveMapFactory
core::kinematics::MoveMapCOP
MinMover::explicitly_set_movemap() const {
	return movemap_;
}

void
MinMover::score_function( ScoreFunctionCOP scorefxn_in )
{
	runtime_assert( scorefxn_in != nullptr );
	scorefxn_ = scorefxn_in;
}

// Why would you want a deep copy of the ScoreFunction here?
// In any case, commenting this out, because PyRosetta cannot distinguish between these two function signatures, always
// choosing this one and leading to unexpected behavior. ~Labonte
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

/// @details restrict the move map by the packer task:
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

	for ( auto const & dof_task : dof_tasks_ ) {
		//generate task
		PackerTaskOP task( dof_task.second->create_task_and_apply_taskoperations( pose ) );

		//modify movemap by task
		Size const nres( task->total_residue() );

		//Set DOFs by task operations:
		for ( Size i(1); i <= nres; ++i ) {
			if ( dof_task.first.first == core::id::PHI ) {
				core::kinematics::MoveMap::MoveMapTorsionID mmid;
				mmid.first = i;
				mmid.second = dof_task.first.second;
				movemap.set( mmid, task->pack_residue( i ) );
			} else {
				for ( Size ia=1,iamax=pose.residue(i).natoms(); ia<=iamax; ++ia ) {
					movemap.set( core::id::DOF_ID( core::id::AtomID(ia,i), dof_task.first.first  ), task->pack_residue(i) );
				}
				//movemap.set( t->first.first, task->pack_residue( i ) );
			}
		}
	}
}

void
MinMover::max_iter( Size max_iter_in ) { min_options_->max_iter( max_iter_in ); }

void
MinMover::minimize(pose::Pose & pose, core::kinematics::MoveMap & active_movemap){
	PROF_START( basic::MINMOVER_APPLY );
	inner_run_minimizer( pose, active_movemap );
	if ( abs_score_convergence_threshold() > 0.0 ) {
		while ( abs_score_diff_after_minimization() > abs_score_convergence_threshold() ) {
			// Make sure further minimizations lead to the same answer
			TR << "running another iteration of minimization. difference is: " << abs_score_diff_after_minimization() << std::endl;
			inner_run_minimizer( pose, active_movemap );
		}
	}
	PROF_STOP( basic::MINMOVER_APPLY );

	// emit statistics
	scorefxn_->show(TR.Debug, pose);
	TR.Debug << std::endl;
}

void
MinMover::inner_run_minimizer( core::pose::Pose & pose, core::kinematics::MoveMap & active_movemap ) {
	if ( !cartesian( ) ) {
		AtomTreeMinimizer minimizer;
		if ( !omega_ ) {
			TR<<"shutting off omega dihedral angle minimization"<<std::endl;
			for ( core::Size i = 1; i <= pose.size(); ++i ) {
				if ( !pose.residue( i ).is_protein() ) continue;
				active_movemap.set( core::id::TorsionID( core::id::omega_torsion, BB, i), false );
			}
		}
		score_before_minimization_ = (*scorefxn_)(pose);
		score_after_minimization_ = minimizer.run( pose, active_movemap, *scorefxn_, *min_options_ );
	} else {
		CartesianMinimizer minimizer;
		score_before_minimization_ = (*scorefxn_)(pose);
		score_after_minimization_ = minimizer.run( pose, active_movemap, *scorefxn_, *min_options_ );
	}

}

void
MinMover::apply(pose::Pose & pose) {
	// lazy default initialization
	MoveMapOP active_movemap( movemap(pose)->clone() );

	apply_dof_tasks_to_movemap(pose, *active_movemap);

	if ( ! scorefxn_ ) scorefxn_ = get_score_function(); // get a default (INITIALIZED!) ScoreFunction

	minimize( pose, *active_movemap );
}

void
MinMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Minimization type:\t" << min_type() << "\nScorefunction:\t\t";
	if ( score_function() != nullptr ) {
		output  << score_function()->get_name() << std::endl;
	} else { output << "none" << std::endl; }
	output << "Score tolerance:\t" << tolerance() << "\nNb list:\t\t" << (nb_list() ? "True" : "False") <<
		"\nDeriv check:\t\t" << (deriv_check() ? "True" : "False") << std::endl << "Movemap:" << std::endl;
	if ( movemap_ != nullptr ) {
		movemap_->show(output);
	}
}

protocols::moves::MoverOP MinMover::clone() const { return protocols::moves::MoverOP( new protocols::minimization_packing::MinMover( *this ) ); }
protocols::moves::MoverOP MinMover::fresh_instance() const { return protocols::moves::MoverOP( new MinMover ); }

void MinMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & )
{
	parse_opts( tag, data, filters, movers );
	parse_movemap_factory( tag, data );
	parse_dof_tasks( tag, data );
}

void MinMover::parse_opts(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const & )
{
	score_function( protocols::rosetta_scripts::parse_score_function( tag, data ) );

	max_iter( tag->getOption< int >( "max_iter", 200 ) );
	min_type( tag->getOption< std::string >( "type", "lbfgs_armijo_nonmonotone" ) );
	tolerance( tag->getOption< core::Real >( "tolerance", 0.01 ) );
	abs_score_convergence_threshold( tag->getOption< core::Real >( "abs_score_convergence_threshold", 0.0 ) );
	cartesian( tag->getOption< bool >( "cartesian", false ) );

	//fpd  if cartesian or nonideal default to lbfgs minimization otherwise the runtime is horrible
	if ( cartesian() && !tag->hasOption("type") ) {
		min_type( "lbfgs_armijo_nonmonotone" );
	}
	if ( ( tag->getOption<bool>("bondangle",false) ||tag->getOption<bool>("bondlength",false) )
			&& !tag->hasOption("type") ) {
		min_type( "lbfgs_armijo_nonmonotone" );
	}
}

void MinMover::parse_movemap_factory( TagCOP const tag, basic::datacache::DataMap & data )
{
	using namespace core::select::movemap;

	MoveMapFactoryOP mmf( new MoveMapFactory );

	if ( tag->hasOption("jump") ) {
		utility::vector1<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );
		// string 'ALL' makes all jumps movable
		if ( jumps.size() == 1 && (jumps[1] == "ALL" || jumps[1] == "All" || jumps[1] == "all" || jumps[1] == "*") ) {
			mmf->all_jumps( true );
		} else if ( tag->getOption< core::Size > ( "jump" ) == 0 ) {
			mmf->all_jumps( false );
		} else {
			for ( std::string const & jump : jumps ) {
				int const value = std::atoi( jump.c_str() ); // convert to C string, then convert to integer (phew!)
				core::select::jump_selector::JumpSelectorCOP jump_select( new core::select::jump_selector::JumpIndexSelector( value ) );
				TR << "Setting min on jump " << value << std::endl;
				mmf->add_jump_action( mm_enable, jump_select );
			}
		}
	}

	bool const chi( tag->getOption< bool >( "chi" ) ), bb( tag->getOption< bool >( "bb" ) );
	omega_ = tag->getOption< bool >( "omega", true );
	mmf->all_chi( chi );
	mmf->all_bb( bb );
	TR<<"Options chi, bb: "<<chi<<", "<<bb<<" omega: "<<omega_<<std::endl;
	if ( tag->hasOption("bondangle") ) {
		bool const value( tag->getOption<bool>("bondangle") );
		mmf->all_bondangles( value );
	}
	if ( tag->hasOption("bondlength") ) {
		bool const value( tag->getOption<bool>("bondlength") );
		mmf->all_bondlengths( value );
	}

	movemap_factory( protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data, false, mmf ) );
}

/// @detail helper function for parse_of_tasks
void
MinMover::parse_dof_task_type(
	std::string const & tag_name,
	core::id::DOF_Type dof_type,
	core::id::TorsionType torsion_type,
	TagCOP const tag,
	basic::datacache::DataMap & data
) {

	if ( !tag->hasOption(tag_name) ) {
		return;
	}

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using std::string;
	using StringVec = utility::vector1<std::string>;

	TaskFactoryOP task_factory( new TaskFactory() );
	string const t_o_val( tag->getOption<string>(tag_name) );
	StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );

	for ( string const & t_o_key : t_o_keys ) {
		if ( data.has( "task_operations", t_o_key ) ) {
			task_factory->push_back( data.get_ptr< TaskOperation >( "task_operations", t_o_key ) );
			TR << "\t " << tag_name << ": " << t_o_key << std::endl;
		} else {
			utility_exit_with_message("TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
		}
	}
	dof_tasks_[ std::make_pair(dof_type, torsion_type)] = task_factory;
}

void
MinMover::parse_dof_tasks(
	TagCOP const tag,
	basic::datacache::DataMap & data
) {

	if (
			tag->hasOption("bb_task_operations") ||
			tag->hasOption("chi_task_operations") ||
			tag->hasOption("bondangle_task_operations") ||
			tag->hasOption("bondlength_task_operations") ) {
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


utility::tag::XMLSchemaComplexTypeGeneratorOP
MinMover::complex_type_generator_for_min_mover( utility::tag::XMLSchemaDefinition & xsd ){

	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute(
		"jump", xs_string,
		"Comma-separated list of jumps to minimize over (be sure this jump exists!). "
		"If set to \"ALL\", all jumps will be set to minimize. "
		"If set to \"0\", jumps will be set not to minimize" )
		+ XMLSchemaAttribute( "abs_score_convergence_threshold", xsct_real , "Keep minimizing until difference before minimization is less than this threshold" )
		+ XMLSchemaAttribute::attribute_w_default(
		"max_iter",  xsct_non_negative_integer,
		"maximum number of iterations allowed. This default is also very loose. "
		"This and the tolerance setting both affect if you will reach convergence",
		"200" )
		+ XMLSchemaAttribute::attribute_w_default(
		"type", xsct_minimizer_type,
		"Minimizer type. linmin, dfpmin, dfpmin_armijo, dfpmin_armijo_nonmonotone. "
		"dfpmin minimzers can also be used with absolute tolerance (add \"atol\" to the minimizer type).",
		"lbfgs_armijo_nonmonotone" )
		+ XMLSchemaAttribute::attribute_w_default(
		"tolerance", xsct_real ,
		"Criteria for convergence of minimization. The default is very loose, "
		"it's recommended to specify something less than 0.01. "
		"max_iter also affects convergence",
		"0.01" )
		+ XMLSchemaAttribute::attribute_w_default( "cartesian", xsct_rosetta_bool , "Perform cartesian minimization?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "bondangle",  xsct_rosetta_bool , "Minimize bond angles?", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "bondlength", xsct_rosetta_bool , "Minimize bond lengths?", "0" );
	attributes
		+ XMLSchemaAttribute::required_attribute( "chi", xsct_rosetta_bool , "Minimize chi angles?" )
		+ XMLSchemaAttribute::required_attribute( "bb",  xsct_rosetta_bool , "Minimize backbone torsion angles?" )
		+ XMLSchemaAttribute::attribute_w_default( "omega", xsct_rosetta_bool , "Minimize omega torsions?", "true" );
	//All of these are lists of task operations, but none use parse_task_operations
	attributes
		+ XMLSchemaAttribute( "bb_task_operations", xsct_task_operation_comma_separated_list , "Task operations specifying residues for backbone minimization" )
		+ XMLSchemaAttribute( "chi_task_operations", xsct_task_operation_comma_separated_list , "Task operations specifying residues for sidechain minimization" )
		+ XMLSchemaAttribute( "bondangle_task_operations", xsct_task_operation_comma_separated_list , "Task operations specifying residues for bond angle minimization" )
		+ XMLSchemaAttribute( "bondlength_task_operations", xsct_task_operation_comma_separated_list , "Task operation specifying residues for bond length minimization" );
	rosetta_scripts::attributes_for_parse_score_function( attributes );
	XMLSchemaSimpleSubelementList subelements;
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, subelements );

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen->complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attributes )
		.set_subelements_repeatable( subelements )
		.add_optional_name_attribute();
	return ct_gen;
}

void MinMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Check parse_movemap

	using namespace utility::tag;

	auto ct_gen = complex_type_generator_for_min_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Does minimization over sidechain and/or backbone" )
		.write_complex_type_to_schema( xsd );
}

Real MinMover::score_diff_after_minimization() const {
	return score_after_minimization_ - score_before_minimization_;
}

Real MinMover::abs_score_diff_after_minimization() const {
	return std::abs( score_diff_after_minimization() );
}

}  // minimization_packing
}  // protocols
