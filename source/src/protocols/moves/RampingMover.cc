// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/RampingMover.cc
/// @brief Mover class for ramping the repulsive term of a score function
/// over the course of several apply() evaluations.
/// @author Monica Berrondo

// Unit Headers
#include <protocols/moves/RampingMover.hh>
#include <protocols/moves/RampingMoverCreator.hh>

// Project headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/MonteCarlo.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// tracer
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <complex>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;
// C++ Headers

// ObjexxFCL Headers
//#include <ObjexxFCL/string.functions.hh>

namespace protocols {
namespace moves {

static basic::Tracer TR( "protocols.moves.RampingMover" );

RampingFunc::~RampingFunc() = default;

LinearFunc::~LinearFunc() = default;

FastLinearFunc::FastLinearFunc(
	Real xval_start_ramp,
	Real xval_end_ramp
) :
	xval_start_ramp_( xval_start_ramp ),
	xval_end_ramp_( xval_end_ramp )
{
	runtime_assert( 0 <= xval_start_ramp && xval_start_ramp <= 1 );
	runtime_assert( 0 <= xval_end_ramp   && xval_end_ramp   <= 1 );
	runtime_assert( xval_start_ramp <= xval_end_ramp );
}

FastLinearFunc::~FastLinearFunc() = default;

RampingFunc::Real
FastLinearFunc::func( Real x ) const
{
	/// If xval_start_ramp_ == xval_end_ramp_; step function.
	/// Make sure 1/(xval_end - xval_start ) is never evaluated!
	if ( x <= xval_start_ramp_ ) return 0.0;
	else if ( x >= xval_end_ramp_ ) return 1.0;
	else {
		return ( x - xval_start_ramp_ ) / (xval_end_ramp_ - xval_start_ramp_ );
	}
}

/// @brief Ramps rapidly from the starting value to the final value.
/// Not 1 at x=1.  Doesn't really finish at (1,1).
GeometricFunc::GeometricFunc( Real xval_at_0p5 ) :
	inv_xval_at_0p5_( 1 / xval_at_0p5 )
{
	runtime_assert( xval_at_0p5 != 0.0 );
}

/// @default ctor -- val of 0.5 produces 0.5
GeometricFunc::GeometricFunc() :
	inv_xval_at_0p5_( 3 )
{

}

GeometricFunc::~GeometricFunc() = default;

/// func(x) = 1 - exp( -1 * x * inv_xval_at_0p5 * 0.6931 );
RampingFunc::Real
GeometricFunc::func( Real x ) const
{
	static Real const ln_0p5 = std::log( 0.5 ); // ~0.6931
	return 1 - std::exp( -1 * x * inv_xval_at_0p5_ * ln_0p5 );

}

/// @brief Ramps slowly from the starting value to the final value
/// Non-zero for x = 0.  Doesn't really start at (0,0).
/// func(x) = exp( -1 * ( 1 - x ) / ( 1 - xval_at_0p5 ) * 0.6931 );
InvGeometricFunc::InvGeometricFunc( Real xval_at_0p5 )
:
	inv_one_minus_xval_at_0p5_( 1 / ( 1 - xval_at_0p5 ))
{}

InvGeometricFunc::InvGeometricFunc() :
	inv_one_minus_xval_at_0p5_( 1 / ( 1 - 0.66 ))
{}

InvGeometricFunc::~InvGeometricFunc() = default;

/// func(x) = exp( -1 * ( 1 - x ) / ( 1 - xval_at_0p5 ) * 0.6931 );
RampingFunc::Real
InvGeometricFunc::func( Real x ) const
{
	static Real const ln_0p5 = std::log( 0.5 );
	return std::exp( -1 * ( 1 - x ) * inv_one_minus_xval_at_0p5_ * ln_0p5 );
}

// XRW TEMP MoverOP RampingMoverCreator::create_mover() const { return MoverOP( new RampingMover ); }
// XRW TEMP std::string RampingMoverCreator::keyname() const { return RampingMover::mover_name(); }
// XRW TEMP std::string RampingMover::mover_name() { return "RampingMover"; }


/// RampingMover

RampingMover::RampingMover() :
	Mover( RampingMover::mover_name() ),
	mover_( /* 0 */ ),
	scorefxn_( /* 0 */ ),
	ramp_one_weight_( true ),
	score_type_( core::scoring::fa_rep ),
	ramping_funcs_for_weights_( core::scoring::n_score_types, nullptr ),
	outer_cycles_( 0 ),
	inner_cycles_( 0 ),
	mc_( /* 0 */ )
{
	ramping_funcs_for_weights_[ score_type_ ] = protocols::moves::RampingFuncOP( new LinearFunc );
}

RampingMover::RampingMover(
	MoverOP mover_in,
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::scoring::ScoreType score_type_in,
	int outer_cycles_in,
	int inner_cycles_in,
	MonteCarloOP  mc_in,
	bool geometric_in
) :
	Mover( RampingMover::mover_name() ),
	mover_(std::move( mover_in )),
	scorefxn_( scorefxn_in ), // replace this with clone() when symmetry comes online.
	ramp_one_weight_( true ),
	score_type_( score_type_in ),
	start_weights_( scorefxn_in->weights() ),
	ramping_funcs_for_weights_( core::scoring::n_score_types, nullptr ),
	outer_cycles_( outer_cycles_in ),
	inner_cycles_( inner_cycles_in ),
	mc_(std::move( mc_in ))
{
	start_weights_[ score_type_ ] = 0.2;
	end_weights_[   score_type_ ] = 1.0;
	ramping_funcs_for_weights_[ score_type_ ] = protocols::moves::RampingFuncOP( geometric_in ? (RampingFunc * ) new GeometricFunc : (RampingFunc * )  new LinearFunc );
}

RampingMover::RampingMover(
	MoverOP mover_in,
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::scoring::EnergyMap start_weights,
	core::scoring::EnergyMap end_weights,
	int outer_cycles_in,
	int inner_cycles_in,
	MonteCarloOP  mc_in
) :
	Mover( RampingMover::mover_name() ),
	mover_(std::move( mover_in )),
	scorefxn_(std::move( scorefxn_in )), // replace this with clone() when symmetry comes online.
	ramp_one_weight_( false ),
	score_type_( core::scoring::fa_rep ), // unused
	start_weights_( start_weights ),
	end_weights_( end_weights ),
	ramping_funcs_for_weights_( core::scoring::n_score_types, nullptr ),
	outer_cycles_( outer_cycles_in ),
	inner_cycles_( inner_cycles_in ),
	mc_(std::move( mc_in ))
{
	for ( Size ii = 1; ii <= core::scoring::n_score_types; ++ii ) {
		core::scoring::ScoreType iist = (core::scoring::ScoreType) ii;
		if ( start_weights_[ iist ] != end_weights_[ iist ] ) {
			ramping_funcs_for_weights_[ iist ] = protocols::moves::RampingFuncCOP( protocols::moves::RampingFuncOP( new LinearFunc ) );
		}
	}
}

RampingMover::~RampingMover() = default;

MoverOP RampingMover::clone() const
{
	return MoverOP( new RampingMover(*this) );
}

void
RampingMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & /*filters*/,
	Movers_map const & movers,
	Pose const & /*pose*/
)
{

	if ( !tag->hasOption("outer_cycles") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "You must specify the option outer_cycles in RampingMover");
	}
	if ( !tag->hasOption("inner_cycles") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "You must specify the option inner_cycles in RampingMover");
	}
	if ( !tag->hasOption("mover") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "You must specify the option mover in RampingMover");
	}

	//Get options values out of the tag
	outer_cycles_ = tag->getOption<int>("outer_cycles");
	inner_cycles_ = tag->getOption<int>("inner_cycles");
	std::string montecarlo_name(tag->getOption<std::string>("montecarlo","none"));
	std::string mover_name(tag->getOption<std::string>("mover"));

	core::scoring::ScoreFunctionOP sfxn( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( sfxn == nullptr ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "scorefxn required to create a RampingMover");
	}
	scorefxn_ = sfxn;


	//get the mover to ramp out of the movemap
	auto find_mover(movers.find(mover_name));
	if ( find_mover == movers.end() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "cannot find "+mover_name+" in mover map.");
	}
	mover_ = find_mover->second;

	//get the montecarlo object out of the datamap
	if ( montecarlo_name != "none" ) {
		mc_ = datamap.get_ptr<protocols::moves::MonteCarlo>("montecarlos", montecarlo_name);
	}

	// Two modes for the Ramping mover:
	// either we ramp one weight, in which case start-weight, end-weight,
	// ramp_func, and scoretype must be provided,
	// OR we ramp several weights, in which case those four should not be provided.
	// CASE 1: Ramp one weight

	if ( tag->hasOption("start_weight") ||
			tag->hasOption("end_weight") ||
			tag->hasOption("scoretype") ||
			tag->hasOption("ramp_func") ) {

		// 1st: Make sure all the required options have been provided.
		if ( !tag->hasOption("start_weight") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "One-weight ramping mode: you must specify the option start_weight in RampingMover");
		}
		if ( !tag->hasOption("end_weight") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "One-weight ramping mode: you must specify the option end_weight in RampingMover");
		}
		if ( !tag->hasOption("score_type") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "One-weight ramping mode: you must specify the option score_type in RampingMover");
		}
		if ( !tag->hasOption("ramp_func") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "One-weight ramping mode: you must specify the option ramp_func in RampingMover");
		}
		ramp_one_weight_ = true;

		core::Real start_weight(tag->getOption<core::Real>("start_weight"));
		core::Real end_weight(tag->getOption<core::Real>("end_weight"));
		std::string score_type(tag->getOption<std::string>("score_type"));
		std::string ramping_func(tag->getOption<std::string>("ramp_func"));

		//set up the start and end weights
		this->start_weight(start_weight);
		this->end_weight(end_weight);


		//set up the score type
		if ( ! core::scoring::ScoreTypeManager::is_score_type( score_type ) ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "The value for option score_type \"" + score_type + "\" is not a valid name for a score_type" );
		}
		score_type_ = core::scoring::ScoreTypeManager::score_type_from_name(score_type);
		ramping_funcs_for_weights_[score_type_] = instantiate_rampfunc( ramping_func, tag );

	} else {
		ramp_one_weight_ = false;
		if ( tag->hasOption("unramped_weights_from_sfxn") ) {
			std::string const scorefxn_key( tag->getOption<std::string>("unramped_weights_from_sfxn") );
			if ( ! datamap.has( "scorefxns", scorefxn_key ) ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
			}
			core::scoring::ScoreFunctionCOP sfxn = datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_key );
			start_weights_ = end_weights_ = sfxn->weights();

		}


		utility::vector0< TagCOP > const & rampterm_tags( tag->getTags() );
		for ( auto tag_ptr : rampterm_tags ) {
			if ( ! tag_ptr->hasOption("score_type") ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Ramping mover Add statement requires the score_type option");
			}
			if ( ! tag_ptr->hasOption("start_weight") ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Ramping mover Add statement requires the start_weight option");
			}
			if ( ! tag_ptr->hasOption("end_weight") ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Ramping mover Add statement requires the end_weight option");
			}
			if ( ! tag_ptr->hasOption("ramp_func") ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Ramping mover Add statement requires the ramp_func option");
			}
			core::Real start_weight(tag_ptr->getOption<core::Real>("start_weight"));
			core::Real end_weight(tag_ptr->getOption<core::Real>("end_weight"));
			std::string score_type(tag_ptr->getOption<std::string>("score_type"));
			std::string ramping_func(tag_ptr->getOption<std::string>("ramp_func"));

			//read the score type
			if ( ! core::scoring::ScoreTypeManager::is_score_type( score_type ) ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "The value for option score_type \"" + score_type + "\" is not a valid name for a score_type" );
			}
			core::scoring::ScoreType st = core::scoring::ScoreTypeManager::score_type_from_name(score_type);
			start_weights_[ st ] = start_weight;
			end_weights_[ st ] = end_weight;

			ramping_funcs_for_weights_[ st ] = instantiate_rampfunc( ramping_func, tag_ptr );
			TR << "Ramping weight " << score_type << std::endl;
		}
	}
}

void
RampingMover::apply( core::pose::Pose & pose )
{
	set_weights( start_weights_ );
	for ( int i = 0; i <= outer_cycles_; ++i ) {
		//T("protocols.moves.RampingMover") << "Move: " << i << "/" << outer_cycles_ << std::endl;
		update_weights( i );
		(*scorefxn_)(pose);

		// reset the score function for the monte carlo object, if any.
		if ( mc_ ) {
			mc_->reset_scorefxn( pose, *scorefxn_ );
		}
		//T("protocols.moves.RampingMover") << "Move: " << i << "/" << outer_cycles_ << " " << scorefxn_->weights()[ core::scoring::fa_rep ] << std::endl;

		// pose.dump_pdb("ramp" + ObjexxFCL::right_string_of(i,2,'0') + "_before.pdb" );
		for ( int j = 0; j < inner_cycles_; ++j ) {
			mover_->apply( pose );
		}
		// pose.dump_pdb("ramp" + ObjexxFCL::right_string_of(i,2,'0') + "_after.pdb" );
		(*scorefxn_)(pose);
		// scorefxn_->show(std::cout, pose);

	}
}

// @details The Ramping mover performs a shallow copy of the input
// score function pointer, so that multiple objects can share the
// same score function and be influenced by this mover's change to
// that score function.
void RampingMover::sfxn( core::scoring::ScoreFunctionOP scorefxn )
{
	scorefxn_ = scorefxn;
	start_weights_ = scorefxn->weights();
}


void
RampingMover::start_weight( core::Real start_weight_in )
{
	runtime_assert( ramp_one_weight_ );
	start_weights_[ score_type_ ] = start_weight_in;
}

void
RampingMover::end_weight( core::Real end_weight_in )
{
	runtime_assert( ramp_one_weight_ );
	end_weights_[ score_type_ ] = end_weight_in;
}

void
RampingMover::set_func_for_weight(
	core::scoring::ScoreType scoretype,
	RampingFuncOP func
)
{
	ramping_funcs_for_weights_[ scoretype ] = func;
}


void
RampingMover::update_weights( int round )
{
	core::Real progress = ((core::Real) round) / outer_cycles_;
	intermediate_weights_ = start_weights_;
	for ( Size ii = 1; ii <= core::scoring::n_score_types; ++ii ) {
		core::scoring::ScoreType iist = (core::scoring::ScoreType) ii;
		if ( ! ramping_funcs_for_weights_[ ii ] ) continue;
		core::Real const alpha = ramping_funcs_for_weights_[ ii ]->func( progress );
		core::Real const beta =  1 - alpha;
		/// if alpha  is outside [0..1], then extrapolate.
		intermediate_weights_[ iist ] = start_weights_[ iist ] * beta + end_weights_[ iist ] * alpha;
	}
	set_weights( intermediate_weights_ );
}

void
RampingMover::set_weights( core::scoring::EnergyMap const & emap )
{
	for ( Size ii = 1; ii <= core::scoring::n_score_types; ++ii ) {
		core::scoring::ScoreType iist = (core::scoring::ScoreType) ii;
		scorefxn_->set_weight( iist, emap[ iist ] );
		if ( emap[ iist ] == 0.0 ) continue;
		//T.Debug("protocols.moves.RampingMover") << iist << " " << emap[ iist ] << " ";
	}
	//T.Debug("protocols.moves.RampingMover") << std::endl;
}

RampingFuncOP
RampingMover::instantiate_rampfunc(
	std::string const & func_name,
	utility::tag::TagCOP tag_ptr
) const
{

	if ( func_name == "geometric" ) {
		core::Real xval_at_0p5 = tag_ptr->getOption< core::Real >( "xval_at_0p5", 0.75 );
		return RampingFuncOP( new GeometricFunc( xval_at_0p5 ) );
	} else if ( func_name == "linear" ) {
		return RampingFuncOP( new LinearFunc );
	} else if ( func_name == "fast_linear" ) {
		core::Real xval_start_ramp = tag_ptr->getOption<core::Real>("xval_start_ramp",0.0);
		core::Real xval_end_ramp = tag_ptr->getOption<core::Real>("xval_end_ramp",1.0);
		return RampingFuncOP( new FastLinearFunc( xval_start_ramp, xval_end_ramp ) );
	} else if ( func_name == "inverse_geometric" ) {
		core::Real xval_at_0p5 = tag_ptr->getOption<core::Real>("xval_at_0p5",0.75);
		return RampingFuncOP( new InvGeometricFunc( xval_at_0p5 ) );
	} else {
		utility_exit_with_message("option ramping_func in RampingMover must be: geometric, linear, fast_linear, or inverse_geometric");
	}

	return RampingFuncOP( new LinearFunc ); // appease compiler.
}

std::string RampingMover::get_name() const {
	return mover_name();
}

std::string RampingMover::mover_name() {
	return "RampingMover";
}

void RampingMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"start_weight", xsct_real,
		"starting weight for ramp");
	attlist + XMLSchemaAttribute::required_attribute(
		"end_weight", xsct_real,
		"ending weight for ramp");
	attlist + XMLSchemaAttribute::required_attribute(
		"scoretype", xs_string,
		"name of the score term to ramp");

	XMLSchemaRestriction ram_func_enumeration;
	ram_func_enumeration.name( "ramp_func_name" );
	ram_func_enumeration.base_type( xs_string );
	ram_func_enumeration.add_restriction( xsr_enumeration, "linear" );
	ram_func_enumeration.add_restriction( xsr_enumeration, "fast_linear" );
	ram_func_enumeration.add_restriction( xsr_enumeration, "geometric" );
	ram_func_enumeration.add_restriction( xsr_enumeration, "inverse_geometric" );
	xsd.add_top_level_element(ram_func_enumeration);
	attlist + XMLSchemaAttribute::required_attribute(
		"ramp_func", "ramp_func_name",
		"the ramp funct to use");

	attlist + XMLSchemaAttribute::required_attribute(
		"inner_cycles", xsct_non_negative_integer,
		"number of times to call inner mover in each score ramp increment");

	attlist + XMLSchemaAttribute::required_attribute(
		"outer_cycles", xsct_non_negative_integer,
		"number of increments to ramp score in");

	attlist + XMLSchemaAttribute::required_attribute(
		"mover", xs_string,
		"name of the inner mover to use");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"montecarlo", xs_string,
		"the name of the montecarlo object to use", "none");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Repeatedly applies a given mover while ramping the score from "
		"a low value to a high value",
		attlist );
}

std::string RampingMoverCreator::keyname() const {
	return RampingMover::mover_name();
}

protocols::moves::MoverOP
RampingMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RampingMover );
}

void RampingMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RampingMover::provide_xml_schema( xsd );
}



} // moves
} // protocols

