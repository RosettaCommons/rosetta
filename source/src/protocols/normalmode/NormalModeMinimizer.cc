// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/normalmode/NormalModeMinimizer.cc
/// @brief  High-level atom tree minimizer class
/// @author Phil Bradley


// Unit headers
#include <protocols/normalmode/NormalMode.hh>
#include <protocols/normalmode/NormalModeMinimizer.hh>
#include <protocols/normalmode/NormalModeMinimizerCreator.hh>
#include <protocols/normalmode/NormalModeMultiFunc.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/jump_selector/JumpIndexSelector.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <iostream>
#include <protocols/rosetta_scripts/util.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using namespace ObjexxFCL::format;
using namespace core;
using namespace kinematics;
using namespace optimization;
using namespace scoring;

namespace protocols {
namespace normalmode {


////////////

static THREAD_LOCAL basic::Tracer nma_tr( "protocols.normalmode", basic::t_debug );


NormalModeMinimizer::NormalModeMinimizer() {
	options_ = MinimizerOptionsOP( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );
	options_->max_iter(200);
}

NormalModeMinimizer::~NormalModeMinimizer() = default;

void
NormalModeMinimizer::apply( pose::Pose & pose ) {
	core::kinematics::MoveMapOP movemap;
	if ( movemap_factory_ ) {
		movemap = movemap_factory_->create_movemap_from_pose( pose );
	} else {
		movemap = core::kinematics::MoveMapOP( new MoveMap );
	}
	if ( ! scorefxn_ ) scorefxn_ = get_score_function(); // get a default (INITIALIZED!) ScoreFunction

	if ( options_->deriv_check() ) {
		deriv_check_result_ = NumericalDerivCheckResultOP( new NumericalDerivCheckResult );
		deriv_check_result_->send_to_stdout( options_->deriv_check_to_stdout() );
	}

	bool const use_nblist( options_->use_nblist() );

	//fpd << fix this!
	runtime_assert( !core::pose::symmetry::is_symmetric( pose ) );

	Real const start_score( (*scorefxn_)( pose ) );
	MinimizerMap min_map;
	min_map.setup( pose, *movemap );

	if ( use_nblist ) {
		pose.energies().set_use_nblist( pose, min_map.domain_map(), options_->nblist_auto_update() );
	}

	scorefxn_->setup_for_minimizing( pose, min_map );

	// Normalmode setup: solve for starting pose
	// We can put in more options for NormalMode here, e.g. distance cutoffs
	// but let's just use default options right now
	protocols::normalmode::NormalMode normalmode;
	normalmode.torsion( true );
	normalmode.solve( pose );

	// setup the function that we will pass to the low-level minimizer
	protocols::normalmode::NormalModeMultifunc
		f( pose, min_map, *scorefxn_, normalmode,
		false, //omega
		options_->deriv_check(), options_->deriv_check_verbose() );

	// Reset modes to be minimized in MultiFunc when defined by user
	if ( modes_using_.size() > 0 ) {
		f.set_modes( modes_using_ );
	}

	if ( deriv_check_result_ ) deriv_check_local( pose, f );

	// starting position -- "dofs" = Degrees Of Freedom
	Multivec dofs( min_map.nangles() );
	min_map.copy_dofs_from_pose( pose, dofs );

	scorefxn_->show( pose );

	// Convert into vars
	Multivec vars( f.dofs_to_vars( dofs ) );

	Real const start_func( f( vars ) );

	nma_tr << "start_func: " << start_func <<  std::endl;

	// now do the optimization with the low-level minimizer function
	Minimizer minimizer( f, *options_ );
	minimizer.run( vars );

	Real const end_func( f( vars ) );

	nma_tr << "end_func: " << end_func << std::endl;

	// turn off nblist
	if ( use_nblist ) pose.energies().reset_nblist();

	dofs = f.vars_to_dofs( vars );
	min_map.reset_jump_rb_deltas( pose, dofs );

	Real const end_score( (*scorefxn_)( pose ) );

	nma_tr << "NormalModeMinimizer::run: nangles= " << min_map.nangles() <<
		" start_score: " << F(12,3,start_score) <<
		" start_func: "  << F(12,3,start_func ) <<
		" end_score: "   << F(12,3,end_score  ) <<
		" end_func: "    << F(12,3,end_func   ) << std::endl;
}


void NormalModeMinimizer::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{

	std::string const scorefxn_name( tag->getOption< std::string >( "scorefxn", "score12" ) );
	scorefxn_ = data.get_ptr<ScoreFunction>( "scorefxns", scorefxn_name );

	// minimizer
	options_->max_iter( tag->getOption< int >( "max_iter", 200 ) );
	options_->min_type( tag->getOption< std::string >( "type", "lbfgs_armijo_nonmonotone" ) );
	options_->minimize_tolerance( tag->getOption< core::Real >( "tolerance", 0.01 ) );

	// nma-min specific
	dampen_ = tag->getOption< core::Real >( "dampen", 0.05 );

	// high-level movemap factory
	core::select::movemap::MoveMapFactoryOP mmf( new core::select::movemap::MoveMapFactory );

	bool const chi( tag->getOption< bool >( "chi" ) ), bb( tag->getOption< bool >( "bb" ) );
	mmf->all_chi( chi );
	mmf->all_bb( bb );

	if ( tag->hasOption("jump") ) {
		utility::vector1<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );

		// string 'ALL' makes all jumps movable
		if ( jumps.size() == 1 && (jumps[1] == "ALL" || jumps[1] == "All" || jumps[1] == "all") ) {
			mmf->all_jumps( true );
		} else if ( tag->getOption< core::Size > ( "jump" ) == 0 ) {
			mmf->all_jumps( false );
		} else {
			for ( std::string const & jump : jumps ) {
				Size const value = std::atoi( jump.c_str() );
				core::select::jump_selector::JumpIndexSelectorOP jumpselect( new core::select::jump_selector::JumpIndexSelector( value ) );
				mmf->add_jump_action( core::select::movemap::mm_enable, jumpselect );
			}
		}
	}

	// fine-grained movemap control
	movemap_factory_ = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data, false, mmf );
}


NumericalDerivCheckResultOP
NormalModeMinimizer::deriv_check_result() const {
	return deriv_check_result_;
}

void
NormalModeMinimizer::deriv_check_local( pose::Pose const &pose,
	protocols::normalmode::NormalModeMultifunc f ) const {

	pose::Pose pose_work( pose );
	Real const dvar( 1.0e-2 );

	Multivec const vars0( f.nvar(), 0.0 );
	Multivec analytical_deriv( f.nvar(), 0.0 );

	// First calculate function and derivative for starting pose
	f.dfunc( vars0, analytical_deriv );

	// Then enumerate over degree of freedom
	Multivec vars;
	nma_tr << LJ(3,"DOF") << " "
		<< LJ(6,"dval") << " "
		<< LJ(10,"g_ana") << " "
		<< LJ(10,"g_num") << " "
		<< LJ(10,"dg") << " "
		<< LJ(10,"dg/g_num") << std::endl;

	for ( Size i_var = 1; i_var <= f.nvar(); ++i_var ) {

		Real f1, f2;
		// +dvar
		vars = vars0; vars[i_var] += dvar;
		f1 = f( vars );

		// -dvar
		vars = vars0; vars[i_var] -= dvar;
		f2 = f( vars );

		Real numeric_deriv = 0.5*(f1 - f2) /dvar;

		// Report
		nma_tr << I(3,i_var) << " "
			<< E(6,1,dvar) << " "
			<< F(10,5,analytical_deriv[i_var]) << " "
			<< F(10,5,numeric_deriv ) << " "
			<< F(10,5,analytical_deriv[i_var] - numeric_deriv  ) << " "
			<< F(10,5,(analytical_deriv[i_var] - numeric_deriv) / numeric_deriv   ) << std::endl;
	}

}

std::string NormalModeMinimizer::get_name() const {
	return mover_name();
}

std::string NormalModeMinimizer::mover_name() {
	return "NormalModeMinimizer";
}

void NormalModeMinimizer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( [] (std::string const& name) {
		return mover_name() + "_subelement_" + name + "Type";
		});
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy(xsd, subelements);

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"scorefxns", xs_string,
		"XRW TO DO",
		"score12");
	attlist + XMLSchemaAttribute::required_attribute(
		"chi", xsct_rosetta_bool,
		"XRW TO DO");
	attlist + XMLSchemaAttribute::required_attribute(
		"bb", xsct_rosetta_bool,
		"XRW TO DO");
	attlist + XMLSchemaAttribute::required_attribute(
		"jump", xs_string,
		"comma separated list or core::Size XRW TO DO");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_iter", xs_integer,
		"XRW TO DO",
		"200");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"min_type", xs_string,
		"XRW TO DO",
		"lbfgs_armijo_nonmonotone");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"tolerance", xsct_real,
		"XRW TO DO",
		"0.01");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"dampen", xsct_real,
		"XRW TO DO",
		"0.05");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"XRW TO DO",
		attlist, subelements );
}

std::string NormalModeMinimizerCreator::keyname() const {
	return NormalModeMinimizer::mover_name();
}

protocols::moves::MoverOP
NormalModeMinimizerCreator::create_mover() const {
	return protocols::moves::MoverOP( new NormalModeMinimizer );
}

void NormalModeMinimizerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NormalModeMinimizer::provide_xml_schema( xsd );
}


} // namespace normalmode
} // namespace protocols
