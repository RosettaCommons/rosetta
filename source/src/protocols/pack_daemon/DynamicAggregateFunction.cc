// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_daemon/DynamicAggregateFunction.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef USEMPI
// MPI Headers
#include <mpi.h>
#endif

// Unit headers
#include <protocols/pack_daemon/DynamicAggregateFunction.hh>

// Package headers
#include <protocols/pack_daemon/PackDaemon.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>
#include <numeric/expression_parser/Arithmetic.hh>
#include <protocols/multistate_design/MultiStatePacker.hh>

// Utility headers
#include <utility>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/mpi_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/FileContentsMap.hh>


namespace protocols {
namespace pack_daemon {

static basic::Tracer TR( "protocols.pack_daemon.DynamicAggregateFunction" );

using namespace numeric::expression_parser;

///// class VectorExpression : public Expression

VectorExpression::VectorExpression( std::string const & name ) : parent(), name_( name ) {}
VectorExpression::~VectorExpression() = default;

core::Real
VectorExpression::operator() () const
{
	utility_exit_with_message( "Illegal call to operator() on VectorExpression named " + name() );
	return 0.0;
}

ExpressionCOP
VectorExpression::differentiate( std::string const & ) const
{
	utility_exit_with_message( "Illegal call to differentiate() on VectorExpression named " + name() );
	return nullptr;
}

std::string const &
VectorExpression::name() const
{
	return name_;
}

/////// class VariableVectorExpression

VariableVectorExpression::VariableVectorExpression(
	std::string const & name,
	VariableExpressions const & vars
) :
	parent( name ),
	vars_( vars )
{}

VariableVectorExpression::~VariableVectorExpression() = default;

VariableVectorExpression::values
VariableVectorExpression::vector_values() const {
	values vals( vars_.size() );
	for ( Size ii = 1; ii <= vars_.size(); ++ii ) {
		vals[ ii ] = (*vars_[ ii ])();
	}
	return vals;
}

core::Size
VariableVectorExpression::size() const {
	return vars_.size();
}

/// @brief DO NOT call this funcion.
std::list< std::string >
VariableVectorExpression::active_variables() const
{
	utility_exit_with_message( "VariableVectorExpression::active_variables should not be called" );
	return std::list< std::string >();
}

utility::vector1< std::list< std::string > >
VariableVectorExpression::active_variables_vector() const
{
	utility::vector1< std::list< std::string > > active_vars_vector( size() );
	for ( Size ii = 1; ii <= vars_.size(); ++ii ) {
		active_vars_vector[ ii ] = vars_[ ii ]->active_variables();
		//for ( std::list< std::string >::const_iterator
		//  iter = active_vars_vector[ ii ].begin(), iter_end = active_vars_vector[ ii ].end();
		//  iter != iter_end; ++iter ) {
		// TR << "DEBUG: active_variables_vector() " << ii << " includes " << *iter << std::endl;
		//}
	}
	return active_vars_vector;
}

///////// class VectorFunction
VectorFunction::VectorFunction( VectorExpressionCOP ex ) : parent( ex ), vec_ex_( ex ) {}
VectorFunction::~VectorFunction() = default;

VectorExpressionCOP VectorFunction::vec_ex() const { return vec_ex_; }


VectorFunction2::VectorFunction2( VectorExpressionCOP ex1, VectorExpressionCOP ex2 ) :
	vec_ex1_(std::move( ex1 )),
	vec_ex2_(std::move( ex2 ))
{}

VectorFunction2::~VectorFunction2() = default;

VectorExpressionCOP VectorFunction2::vec_ex1() const { return vec_ex1_; }
VectorExpressionCOP VectorFunction2::vec_ex2() const { return vec_ex2_; }


IterativeVectorExpression::IterativeVectorExpression( std::string const & name ) : parent( name ) {}

IterativeVectorExpression::~IterativeVectorExpression() = default;

void IterativeVectorExpression::initialize(
	std::map< std::string, VectorExpressionCOP > const & vector_varnames,
	ArithmeticASTExpression const & expresion_ast,
	VectorExpressionCreator & expression_creator // holds a reference to my owning DynamicAggregateFunction
)
{
	Size s = vector_varnames.size();
	input_vector_expressions_.resize( s );
	local_variables_.resize( s );
	Size count = 0;
	for ( auto const & vector_varname : vector_varnames ) {
		++count;
		VariableExpressionOP varex( new VariableExpression( vector_varname.first, 0.0 ) );
		input_vector_expressions_[ count ] = vector_varname.second;
		local_variables_[ count ] = varex;
		local_variable_map_[ vector_varname.first ] = varex;
	}
	for ( Size ii = 2; ii <= s; ++ii ) {
		if ( input_vector_expressions_[ ii ]->size() != input_vector_expressions_[ ii - 1 ]->size() ) {
			utility_exit_with_message( "IterativeVectorExpression " + name() + " initialized with vector-expressions of uneven sizes: "
				+ input_vector_expressions_[ ii   ]->name() + " with size of " + utility::to_string( input_vector_expressions_[ ii   ]->size()) + ", and "
				+ input_vector_expressions_[ ii-1 ]->name() + " with size of " + utility::to_string( input_vector_expressions_[ ii-1 ]->size()) );
		}
	}

	/// This call will construct an expression tree from the ExpressionAST.  Control of flow
	/// will actually return to this object when the expression_creator encounters local variables.
	/// The expression_creator will hand control of flow to its owning DynamicAggregateFunction,
	/// and the DAF will hand control of flow to this object.  A pre-condition of this
	/// function call is that the DAF which is initializing this object must point to it with its
	/// focused_iterative_vector_expression_ member variable.

	expression_ = expression_creator.create_expression_tree( expresion_ast );

}

IterativeVectorExpression::values
IterativeVectorExpression::vector_values() const
{
	Size const n_vals = size();
	Size const n_inputs = input_vector_expressions_.size();

	if ( n_vals == 0 ) { values empty; return empty; }

	utility::vector1< values > vec_of_values( n_inputs );

	for ( Size ii = 1; ii <= n_inputs; ++ii ) {
		vec_of_values[ ii ] = input_vector_expressions_[ ii ]->vector_values();
		if ( ii != 1 ) {
			if ( vec_of_values[ ii ].size() != vec_of_values[ ii - 1 ].size() ) {
				utility_exit_with_message( "IterativeVectorExpression " + name() + " recieved vectors of uneven sizes from its children: "
					+ input_vector_expressions_[ ii   ]->name() + " with " + utility::to_string( vec_of_values[ ii   ].size()) + " values, and "
					+ input_vector_expressions_[ ii-1 ]->name() + " with " + utility::to_string( vec_of_values[ ii-1 ].size()) + " values." );
			}
		}
	}
	values return_vals( n_vals, 0.0 );
	for ( Size ii = 1; ii <= n_vals; ++ii ) {
		for ( Size jj = 1; jj <= n_inputs; ++jj ) {
			local_variables_[ jj ]->set_value( vec_of_values[ jj ][ ii ] );
		}
		return_vals[ ii ] = (*expression_)();
	}
	return return_vals;
}

core::Size
IterativeVectorExpression::size() const
{
	if ( input_vector_expressions_.size() == 0 ) {
		return 0;
	}
	return input_vector_expressions_[ 1 ]->size();
}

/// @details returns 0 if there is no local variable with name varname.
/// Do not throw an error if there is no such local variable -- the
/// initializing DAF will check its set of variable names after this object has
/// looked through its set of local variables.
VariableExpressionCOP
IterativeVectorExpression::local_variable( std::string const & varname ) const
{
	auto varit = local_variable_map_.find( varname );
	if ( varit != local_variable_map_.end() ) {
		return varit->second;
	}
	return nullptr;
}

std::list< std::string >
IterativeVectorExpression::active_variables() const
{
	utility_exit_with_message( "IterativeVectorExpression::active_variables should not be called" );
	return std::list< std::string >();
}

utility::vector1< std::list< std::string > >
IterativeVectorExpression::active_variables_vector() const
{
	core::Size const mysize = size();
	utility::vector1< std::list< std::string > > active_varibles_vector( mysize );
	for ( Size ii = 1; ii <= input_vector_expressions_.size(); ++ii ) {
		utility::vector1< std::list< std::string > > ii_active_varibles_vector = input_vector_expressions_[ ii ]->active_variables_vector();
		for ( Size jj = 1; jj <= mysize; ++jj ) {
			active_varibles_vector[ jj ].splice( active_varibles_vector[ jj ].end(), ii_active_varibles_vector[ jj ] );
		}
	}
	return active_varibles_vector;
}

/////// class Mean : public VectorFunction

Mean::Mean( VectorExpressionCOP ex ) : parent( ex ) {}
Mean::~Mean() = default;

core::Real
Mean::operator() () const
{
	VectorExpression::values vals = vec_ex()->vector_values();
	core::Real sum = 0;
	for ( core::Size ii = 1; ii <= vals.size(); ++ii ) {
		sum += vals[ ii ];
	}
	return sum / vals.size();
}

ExpressionCOP
Mean::differentiate( std::string const & ) const
{
	utility_exit_with_message( "Mean cannot be differentiated" );
	return nullptr;
}

std::list< std::string >
Mean::active_variables() const
{
	utility::vector1< std::list< std::string > > act_vars_vect = vec_ex()->active_variables_vector();
	std::list< std::string > all_vars_union;
	for ( core::Size ii = 1; ii <= act_vars_vect.size(); ++ii ) {
		all_vars_union.splice( all_vars_union.end(), act_vars_vect[ ii ] );
	}
	return all_vars_union;
}


/////// class VMax : public VectorFunction

VMax::VMax( VectorExpressionCOP ex ) : parent( ex ) {}
VMax::~VMax() = default;

core::Real
VMax::operator() () const
{
	VectorExpression::values vals = vec_ex()->vector_values();
	return utility::max( vals );
}

ExpressionCOP
VMax::differentiate( std::string const & ) const
{
	utility_exit_with_message( "VMax cannot be differentiated" );
	return nullptr;
}

std::list< std::string >
VMax::active_variables() const
{
	utility::vector1< std::list< std::string > > act_vars_vect = vec_ex()->active_variables_vector();
	VectorExpression::values vals = vec_ex()->vector_values();
	Size index = utility::arg_max( vals );
	return act_vars_vect[ index ];
}

/////// class VMin : public VectorFunction

VMin::VMin( VectorExpressionCOP ex ) : parent( ex ) {}

VMin::~VMin() = default;

core::Real
VMin::operator() () const
{
	VectorExpression::values vals = vec_ex()->vector_values();
	return utility::min( vals );
}

ExpressionCOP
VMin::differentiate( std::string const & ) const
{
	utility_exit_with_message( "VMin cannot be differentiated" );
	return nullptr;
}

std::list< std::string >
VMin::active_variables() const
{
	utility::vector1< std::list< std::string > > act_vars_vect = vec_ex()->active_variables_vector();
	VectorExpression::values vals = vec_ex()->vector_values();
	Size index = utility::arg_min( vals );
	//for ( std::list< std::string >::const_iterator iter = act_vars_vect[ index ].begin(),
	//    iter_end = act_vars_vect[ index ].end(); iter != iter_end; ++iter ) {
	// TR << "DEBUG Vmin active variables for index " << index << " " << *iter << std::endl;
	//}
	//TR << "DEBUG Vmin end active variables" << std::endl;
	return act_vars_vect[ index ];
}


VMaxBy::VMaxBy( VectorExpressionCOP ex1, VectorExpressionCOP ex2 ) : parent( ex1, ex2 ) {}
VMaxBy::~VMaxBy() = default;

core::Real
VMaxBy::operator() () const
{
	VectorExpression::values vals1 = vec_ex1()->vector_values();
	VectorExpression::values vals2 = vec_ex2()->vector_values();
	debug_assert( vals2.size() == vals1.size() );

	Size index = utility::arg_max( vals1 );
	return vals2[ index ];
}

ExpressionCOP
VMaxBy::differentiate( std::string const & ) const
{
	utility_exit_with_message( "VMaxBy cannot be differentiated" );
	return nullptr;
}

std::list< std::string >
VMaxBy::active_variables() const
{
	VectorExpression::values vals1 = vec_ex1()->vector_values();
	utility::vector1< std::list< std::string > > act_vars_vect2 = vec_ex2()->active_variables_vector();
	debug_assert( vals1.size() == act_vars_vect2.size() );

	Size index = utility::arg_max( vals1 );
	return act_vars_vect2[ index ];

}

VMinBy::VMinBy( VectorExpressionCOP ex1, VectorExpressionCOP ex2 ) : parent( ex1, ex2 ) {}

VMinBy::~VMinBy() = default;

core::Real
VMinBy::operator() () const
{
	VectorExpression::values vals1 = vec_ex1()->vector_values();
	VectorExpression::values vals2 = vec_ex2()->vector_values();
	debug_assert( vals1.size() == vals2.size() );

	Size index = utility::arg_min( vals1 );
	return vals2[ index ];
}

ExpressionCOP
VMinBy::differentiate( std::string const & ) const
{
	utility_exit_with_message( "VMinBy cannot be differentiated" );
	return nullptr;
}

std::list< std::string >
VMinBy::active_variables() const
{
	VectorExpression::values vals1 = vec_ex1()->vector_values();
	utility::vector1< std::list< std::string > > act_vars_vect2 = vec_ex2()->active_variables_vector();
	debug_assert( vals1.size() == act_vars_vect2.size() );

	Size index = utility::arg_min( vals1 );
	return act_vars_vect2[ index ];

}

/////// class PowExpression : public BinaryExpression

PowExpression::PowExpression( ExpressionCOP base, ExpressionCOP exponent ) : parent( base, exponent ) {}
PowExpression::~PowExpression() = default;

core::Real
PowExpression::operator() () const
{
	return std::pow( (*e1())(), (*e2())() );
}

ExpressionCOP
PowExpression::differentiate( std::string const & ) const
{
	utility_exit_with_message( "PowExpression::diffrentiate is unimplemented!" );
	return nullptr;
}

/////// class ExpExpression : public UnaryExpression

ExpExpression::ExpExpression( ExpressionCOP ex ) : parent( ex ) {}
ExpExpression::~ExpExpression() = default;

core::Real
ExpExpression::operator() () const
{
	return std::exp( (*ex())() );
}

ExpressionCOP
ExpExpression::differentiate( std::string const & ) const
{
	utility_exit_with_message( "ExpExpression::diffrentiate is unimplemented!" );
	return nullptr;
}


/////// class LnExpression : public UnaryExpression

LnExpression::LnExpression( ExpressionCOP ex ) : parent( ex ) {}
LnExpression::~LnExpression() = default;

core::Real
LnExpression::operator() () const
{
	return std::log( (*ex())() );
}

ExpressionCOP
LnExpression::differentiate( std::string const & ) const
{
	utility_exit_with_message( "LnExpression::diffrentiate is unimplemented!" );
	return nullptr;
}


InSetExpression::InSetExpression( ExpressionCOP ex ) : parent( ex ) {}
InSetExpression::~InSetExpression() = default;

void InSetExpression::value_set( utility::vector1< core::Real > const & values )
{
	value_set_ = values;
}

core::Real
InSetExpression::operator() () const
{
	core::Real val = (*ex())();
	for ( Size ii = 1; ii <= value_set_.size(); ++ii ) {
		if ( val == value_set_[ ii ] ) {
			return 1.0;
		}
	}
	return 0.0;
}

ExpressionCOP
InSetExpression::differentiate( std::string const & ) const
{
	utility_exit_with_message( "InSetExpression::diffrentiate is unimplemented!" );
	return nullptr;
}


/////// class VectorExpressionCreator : public ExpressionCreator

VectorExpressionCreator::VectorExpressionCreator( DynamicAggregateFunction const & owner ) : owner_( owner ) {}
VectorExpressionCreator::~VectorExpressionCreator() = default;

VectorExpressionCreator::ExpressionCOP
VectorExpressionCreator::handle_variable_expression( ArithmeticASTValue const & node )
{
	return owner_.variable_expression( node );
}

VectorExpressionCreator::ExpressionCOP
VectorExpressionCreator::handle_function_expression(
	FunctionTokenCOP function,
	utility::vector1< ExpressionCOP > const & args
)
{
	ExpressionCOP parent_expression = parent::handle_function_expression( function, args );
	if ( parent_expression ) {
		return parent_expression; // handles sqrt only
	}
	return owner_.function_expression( function, args );
}


///////

SurrogateVariableExpression::SurrogateVariableExpression(
	std::string const & name
) :
	parent( name ),
	root_expression_( /* 0 */ )
{}

SurrogateVariableExpression::SurrogateVariableExpression(
	std::string const & name,
	core::Real value
) :
	parent( name, value ),
	root_expression_( /* 0 */ )
{}

/// @details If you don't want to consider any of the variables
/// to be active in a SurrogateVariableExpression, then do not
/// set the root_expression_ for the SurrogateVariableExpression.
/// This allows non-expressions (e.g. the EntityFuncs) to still
/// contribute to the FITNESS without having to correspond to a
/// particular state.
std::list< std::string >
SurrogateVariableExpression::active_variables() const
{
	if ( root_expression_ ) {
		return root_expression_->active_variables();
	} else {
		std::list< std::string > empty_list;
		return empty_list;
	}
}

void
SurrogateVariableExpression::root_expression( ExpressionCOP setting )
{
	root_expression_ = setting;
}

ExpressionCOP
SurrogateVariableExpression::differentiate( std::string const & varname ) const
{
	return root_expression_->differentiate( varname );
}


///////

DynamicAggregateFunction::DynamicAggregateFunction() :
	num_entity_elements_( 0 ),
	file_contents_( utility::io::FileContentsMapOP( new utility::io::FileContentsMap ) )
{}

DynamicAggregateFunction::~DynamicAggregateFunction() = default;

void
DynamicAggregateFunction::set_num_entity_elements( Size setting ) {
	if ( setting == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "DynamicAggregateFunction::set_num_entity_elements may not be"
			"passed a number of entity elements == 0" );
	}
	if ( num_entity_elements_ != 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "DynamicAggregateFunction::set_num_entity_elements may only"
			"be called once" );
	}

	// set this once
	num_entity_elements_ = setting;
}

/// @details Required for processing the POSE_ENERGY and POSE_ENERGY_VECTOR commands
void
DynamicAggregateFunction::set_score_function( core::scoring::ScoreFunction const & sfxn )
{
	sfxn_ = sfxn.clone();
}


core::Size
DynamicAggregateFunction::num_states() const
{
	return variable_expressions_for_states_.size();
}

core::Size
DynamicAggregateFunction::num_npd_properties() const
{
	return variable_expressions_for_npd_properties_.size();
}

core::Size
DynamicAggregateFunction::num_npd_properties_for_state( core::Size state_id ) const
{
	return npd_variable_indices_for_states_[ state_id ].size();
}


core::Real
DynamicAggregateFunction::evaluate( StateEnergies const & state_energies, StateEnergies const & npd_properties, Entity const & entity )
{
	/// Assign all variables, evaluate scalar_expressions and vector expressions, then evaluate the final fitness expression
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "DAF::eval";
		for ( Size ii = 1; ii <= state_energies.size(); ++ii ) {
			TR.Debug << " " << state_energies[ ii ];
		}
	}
	assign_state_energies_to_variables_and_subexpressions( state_energies, npd_properties, entity );
	//TR << "Finished sub expression assignment " << std::endl;
	//TR << "Fitness expression " << fitness_exp_() << std::endl;

	core::Real score = (*fitness_exp_)();
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << " s: " << score << std::endl;
	}

	return score;
}

/// @details Relies on the Expression "active_variables" method
/// to determine which variables are selected and contribute to
/// the fitness function
DynamicAggregateFunction::StateIndices
DynamicAggregateFunction::select_relevant_states(
	StateEnergies const & state_energies,
	StateEnergies const & npd_properties,
	Entity const & entity
)
{
	TR << "Recovering relevant states" << std::endl;
	assign_state_energies_to_variables_and_subexpressions( state_energies, npd_properties, entity, true );
	std::list< std::string > active_varnames = fitness_exp_->active_variables();
	StateIndices active_indices; active_indices.reserve( active_varnames.size() );
	for ( std::list< std::string >::const_iterator
			iter = active_varnames.begin(), iter_end = active_varnames.end();
			iter != iter_end; ++iter ) {
		if ( state_variable_name_2_state_index_.find( *iter ) == state_variable_name_2_state_index_.end() ) {
			// Ignore active variables that are not in the state_variable_2_state_index_ map
			// These include POSE_ENERGY and POSE_ENERGY_VECTOR variables as well as
			// SCALAR_EXPRESSION variables that simply hold constants (e.g. "SCALAR_EXPRESSION five = 2 + 3" )
			TR << "The variable named '" << *iter << "' contributes to the fitness but does not correspond to a state." << std::endl;
			continue;
		}
		TR << "  fitness_exp_ identifies active variable: " << *iter << " " << state_variable_name_2_state_index_.find( *iter )->second << std::endl;
		active_indices.push_back( state_variable_name_2_state_index_.find( *iter )->second );
	}
	std::sort( active_indices.begin(), active_indices.end() );
	StateIndices::const_iterator end_iter = std::unique( active_indices.begin(), active_indices.end() );
	StateIndices return_indices; return_indices.reserve( active_indices.size() );
	StateIndices::const_iterator iter = active_indices.begin();
	while ( iter != end_iter ) {
		return_indices.push_back( *iter );
		++iter;
	}
	for ( core::Size ii = 1; ii <= return_indices.size(); ++ii ) {
		TR << "  FitnessFunction relevant state indices " << return_indices[ ii ] << std::endl;
	}
	return return_indices;
}

/// @details returns a VariableExpressionOP or VectorVariableExpressionOP given a particular
/// variable.  Throws an exception if a variable is requested that is not a member of
/// any of the following maps: named_state_expression_map_, state_vector_variables_, or
/// sub_expression_map_.
ExpressionCOP
DynamicAggregateFunction::variable_expression( ArithmeticASTValue const & var_node ) const
{
	if ( var_node.is_literal() ) {
		utility_exit_with_message( "Error in DynamicAggregateFunction::variable_expression; non-variable (literal) node recieved" +
			utility::to_string( var_node.literal_value() ));
	}

	if ( focused_iterative_vector_expression_ ) {
		ExpressionCOP iterative_vector_expression_local_variable =
			focused_iterative_vector_expression_->local_variable( var_node.variable_name() );
		if ( iterative_vector_expression_local_variable ) return iterative_vector_expression_local_variable;
	}

	auto named_state_iter =
		named_state_expression_map_.find( var_node.variable_name() );
	if ( named_state_iter != named_state_expression_map_.end() ) {
		return named_state_iter->second;
	}

	auto vector_variable_iter =
		state_vector_variables_.find( var_node.variable_name() );
	if ( vector_variable_iter != state_vector_variables_.end() ) {
		return vector_variable_iter->second;
	}

	auto scalar_expression_iter =
		scalar_expression_map_.find( var_node.variable_name() );
	if ( scalar_expression_iter != scalar_expression_map_.end() ) {
		return scalar_expression_iter->second;
	}

	auto vector_expression_iter =
		vector_expression_map_.find( var_node.variable_name() );
	if ( vector_expression_iter != vector_expression_map_.end() ) {
		return vector_expression_iter->second;
	}

	auto
		entity_funcs_iter = entity_funcs_.find( var_node.variable_name() );
	if ( entity_funcs_iter != entity_funcs_.end() ) {
		return entity_funcs_iter->second.second;
	}

	throw CREATE_EXCEPTION(utility::excn::Exception,  "Unexpected variable expression encountered while "
		"forming an expression from an ExpressionAST (after scanning/parsing completed!): " + var_node.variable_name() );

	return nullptr;
}

/// @details Handles the functions that the ExpressionCreator base class does not.
/// Ensures that the functions that are expecting vector arguments are actually given
/// vector arguments.  This is not guaranteed by the scanning or parsing of the input
/// file.
ExpressionCOP
DynamicAggregateFunction::function_expression(
	FunctionTokenCOP function,
	utility::vector1< ExpressionCOP > const & args
) const
{

	std::string const fname = function->name();
	if ( fname == "mean" ) {
		utility::vector1< VectorExpressionCOP > vector_expressions = verify_vector_arguments( fname, args, 1 );
		return ExpressionCOP( ExpressionOP( new Mean( vector_expressions[1] ) ) );
	} else if ( fname == "vmax" ) {
		utility::vector1< VectorExpressionCOP > vector_expressions = verify_vector_arguments( fname, args, 1 );
		return ExpressionCOP( ExpressionOP( new VMax( vector_expressions[1] ) ) );
	} else if ( fname == "vmin" ) {
		utility::vector1< VectorExpressionCOP > vector_expressions = verify_vector_arguments( fname, args, 1 );
		return ExpressionCOP( ExpressionOP( new VMin( vector_expressions[1] ) ) );
	} else if ( fname == "vmax_by" ) {
		utility::vector1< VectorExpressionCOP > vector_expressions = verify_vector_arguments( fname, args, 2 );
		return ExpressionCOP( ExpressionOP( new VMaxBy( vector_expressions[1], vector_expressions[2] ) ) );
	} else if ( fname == "vmin_by" ) {
		utility::vector1< VectorExpressionCOP > vector_expressions = verify_vector_arguments( fname, args, 2 );
		return ExpressionCOP( ExpressionOP( new VMinBy( vector_expressions[1], vector_expressions[2] ) ) );
	} else if ( fname == "exp" ) {
		if ( args.size() != 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "exp expression construction requested with more than one argument: " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new ExpExpression( args[ 1 ] ) ) );
	} else if ( fname == "ln" ) {
		if ( args.size() != 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "ln expression construction requested with more than one argument: " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new LnExpression( args[ 1 ] ) ) );
	} else if ( fname == "pow" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "pow expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new PowExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "ite" ) {
		if ( args.size() != 3 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "ite expression construction requested with nargs != 3. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new ITEExpression( args[ 1 ], args[ 2 ], args[ 3 ] ) ) );
	} else if ( fname == "abs" ) {
		if ( args.size() != 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "abs expression construction requested with nargs != 1. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new AbsoluteValueExpression( args[ 1 ] ) ) );
	} else if ( fname == "gt" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "gt expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new GT_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "lt" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "lt expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new LT_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "gte" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "gte expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new GTE_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "lte" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "lte expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new LTE_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "and" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "and expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new AndExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "or" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "or expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new OrExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "not" ) {
		if ( args.size() != 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "not expression construction requested with nargs != 1. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new NotExpression( args[ 1 ] ) ) );
	}
	throw CREATE_EXCEPTION(utility::excn::Exception,  "Unrecognized function requested of DynamicAggregateFunction: " + fname );
	return nullptr;
}

void DynamicAggregateFunction::add_file_contents(
	std::string const & fname,
	std::string const & contents
){
	file_contents_->set_file_contents( fname, contents );
}

std::string
DynamicAggregateFunction::state_name( Size state_index ) const
{
	return variable_expressions_for_states_[ state_index ]->name();
}

StructureFileNames const &
DynamicAggregateFunction::file_inputs_for_job( int job_index ) const
{
	return files_for_state_[ job_index ];
}

void DynamicAggregateFunction::read_all_variables_from_input_file( std::istream & input )
{

	initialize_scanner();

	function_names_.clear();
	function_names_.insert( "mean" );
	function_names_.insert( "vmax" );
	function_names_.insert( "vmin" );
	function_names_.insert( "vmax_by" );
	function_names_.insert( "vmin_by" );
	function_names_.insert( "ln" );
	function_names_.insert( "pow" );
	function_names_.insert( "exp" );
	function_names_.insert( "sqrt" );
	function_names_.insert( "ite" );
	function_names_.insert( "abs" );
	function_names_.insert( "gt" );
	function_names_.insert( "lt" );
	function_names_.insert( "gte" );
	function_names_.insert( "lte" );
	function_names_.insert( "and" );
	function_names_.insert( "or" );
	function_names_.insert( "not" );

	illegal_variable_names_.clear();
	illegal_variable_names_.insert( "min" );
	illegal_variable_names_.insert( "max" );

	utility::vector1< std::pair< std::string, std::string > > strucvec_filenames;
	std::map< std::string, ArithmeticASTExpressionOP > scalar_expression_asts;
	std::map< std::string, std::list< std::string > > vector_variables;
	std::map< std::string, std::pair< std::map< std::string, std::string >, ArithmeticASTExpressionOP > > vector_expression_asts;

	ArithmeticASTExpressionOP fitness_expression_ast( nullptr );

	named_state_data_file_names_.clear();
	state_vector_data_file_names_.clear();

	Size count_line( 0 );
	while ( input ) {
		++count_line;
		std::string line;
		std::getline( ( std::istream & ) input, line );
		if ( line.size() == 0 ) continue; // skip blank lines

		std::istringstream input_line( line );
		if ( input_line.peek() == '#' ) continue; // skip comment lines

		std::string command;
		input_line >> command;
		if ( command == "STATE" ) {
			process_STATE_line( line, count_line, input_line );
		} else if ( command == "STATE_VECTOR" ) {
			process_STATE_VECTOR_line( line, count_line, input_line, strucvec_filenames );
		} else if ( command == "POSE_ENERGY" ) {
			process_POSE_ENERGY_line( line, count_line, input_line );
		} else if ( command == "POSE_ENERGY_VECTOR" ) {
			process_POSE_ENERGY_VECTOR_line( line, count_line, input_line );
		} else if ( command == "NPD_PROPERTY" ) {
			process_NPD_PROPERTY_line( line, count_line, input_line );
		} else if ( command == "VECTOR_VARIABLE" ) {
			process_VECTOR_VARIABLE_line( line, count_line, input_line, vector_variables );
		} else if ( command == "SCALAR_EXPRESSION" ) {
			process_SCALAR_EXPRESSION_line( line, count_line, input_line, scalar_expression_asts );
		} else if ( command == "VECTOR_EXPRESSION" ) {
			process_VECTOR_EXPRESSION_line( line, count_line, input_line, vector_expression_asts );
		} else if ( command == "ENTITY_FUNCTION" ) {
			process_ENTITY_FUNCTION_line( line, count_line, input_line );
		} else if ( command == "FITNESS" ) {
			process_FITNESS_line( line, count_line, input_line, fitness_expression_ast );
		} else {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Unable to recognize command '" + command + "' while reading DynamicAggregateFunction file" );
		}
	}

	if ( ! fitness_expression_ast ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Unable to find FITNESS command in the file describing the DynamicAggretateFunction.\n"
			"The FITNESS command must exist  " );
	}

	/// NOW:
	/// 1. Read the STATE_VECTOR input files into the state_vector_data_file_names_ variable.
	///    Make sure that each input file is properly formatted.
	/// 2. Turn the ExpressionASTs into Expressions.
	///    a. Count the numer of states (structures)
	///    b. Create variables for the scalar- and vector-expressions.
	///    c. Create Expression trees for the scalar- and vector-expressions in the order they were defined.
	///    d. Create the FITNESS expression.

	Size n_vector_states( 0 );
	for ( Size ii = 1; ii <= strucvec_filenames.size(); ++ii ) {
		read_state_vector_file( strucvec_filenames[ ii ].first, strucvec_filenames[ ii ].second, n_vector_states );
	}

	Size const n_states = named_state_data_file_names_.size() + n_vector_states;
	Size const n_npd_properties = count_num_npd_properties();
	Size const n_variables = n_states + n_npd_properties + scalar_expression_asts.size();

	//state_names_.resize( n_states );
	variable_expressions_for_states_.resize( n_states );
	npd_variable_indices_for_states_.resize( n_states );
	variable_expressions_for_npd_properties_.resize( n_npd_properties );
	files_for_state_.resize( n_states );
	variable_expressions_.resize( n_variables );

	Size count_state( 0 ), count_npd_index( 0 ), count_variable_index( 0 );
	create_state_variable_expressions( count_state, count_npd_index, count_variable_index );
	create_variable_vector_expressions( count_state, count_npd_index, count_variable_index );
	create_scalar_and_vector_expression_variable_expressions( scalar_expression_asts, vector_variables, count_variable_index );

	turn_expression_ASTs_into_expressions( scalar_expression_asts, vector_expression_asts, fitness_expression_ast );

	/// Basically done.
	/// At this point, I should verify that the expressions can be evaulated -- that is,
	/// that the vector-variables are only handed to vector functions.
	/// I don't know the best way to do that.
}

std::list< std::pair< core::Size, std::string > >::const_iterator
DynamicAggregateFunction::npd_variable_indices_for_state_begin( core::Size state_id ) const
{
	return npd_variable_indices_for_states_[ state_id ].begin();
}

std::list< std::pair< core::Size, std::string > >::const_iterator
DynamicAggregateFunction::npd_variable_indices_for_state_end( core::Size state_id ) const {
	return npd_variable_indices_for_states_[ state_id ].end();
}


void
DynamicAggregateFunction::initialize_scanner()
{
	scanner_ = numeric::expression_parser::ArithmeticScannerOP( new ArithmeticScanner( false ) ); /// constructor without adding the "standard" functions
	scanner_->add_function( "sqrt", 1 ); // neiter min nor max are allowed functions.
	scanner_->add_function( "mean", 1 );
	scanner_->add_function( "vmax", 1 );
	scanner_->add_function( "vmin", 1 );
	scanner_->add_function( "vmax_by", 2 );
	scanner_->add_function( "vmin_by", 2 );
	scanner_->add_function( "exp", 1 );
	scanner_->add_function( "pow", 2 );
	scanner_->add_function( "ite", 3 );
	scanner_->add_function( "ln", 1 );
	scanner_->add_function( "abs", 1 );
	scanner_->add_function( "gt", 2 );
	scanner_->add_function( "lt", 2 );
	scanner_->add_function( "gte", 2 );
	scanner_->add_function( "lte", 2 );
	scanner_->add_function( "and", 2 );
	scanner_->add_function( "or", 2 );
	scanner_->add_function( "not", 1 );
}

/// @details  Reads the contents of the line which should include four things:
/// 1. an as-of-yet unused variable name to identify this state,
/// 2. a pdb name giving the file containing the state to be redesigned
/// 3. a correspondence file name which describes the connection between residues
///    in the state from 2. which should correspond to particular residues
///    that are being designed.
/// 4. a secondary resfile that describes how other residues in 2. should be repacked
///    or redesigned.
/// The variable name in 1. is checked against the variable_names_dec_line_ member variable,
/// and, if it is absent, be added to both variable_names_dec_line_ and scanner_.  This method
/// then adds the file names from 2, 3, and 4 to the named_state_data_file_names_
/// map.
void
DynamicAggregateFunction::process_STATE_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line
)
{
	std::string state_name, structure_file, correspondence_file, secondary_resfile;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read state name in the DynamicAggregateFunction"
			" input file after reading STATE on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> state_name;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read state variable name in the DynamicAggregateFunction"
			" input file after reading STATE name on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( state_name.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read state variable name in the DynamicAggregateFunction"
			" input file after reading STATE on line " + utility::to_string( line_number ) + "\n" + line );
	}

	verify_variable_name_or_throw( state_name, "STATE", line, line_number );

	input_line >> structure_file;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read pdb file name in the DynamicAggregateFunction"
			" input file after reading STATE pdb file on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( structure_file.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read pdb file name in the DynamicAggregateFunction"
			" input file after reading STATE name on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> correspondence_file;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read secondary resfile name in the DynamicAggregateFunction"
			" input file after reading STATE correspondence file on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( correspondence_file.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read correspondence file name in the DynamicAggregateFunction"
			" input file after reading STATE pdb file on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> secondary_resfile;
	if ( secondary_resfile.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read secondary resfile name in the DynamicAggregateFunction"
			" input file after reading STATE correspondence file on line " + utility::to_string( line_number ) + "\n" + line );
	}

	// Success.
	TR << "Read STATE line.  Name= " << state_name
		<< " PDB= " << structure_file << " Corr= " << correspondence_file
		<< " 2RF= " << secondary_resfile << std::endl;


	StructureFileNames sfn;
	sfn.pdb_name_ = structure_file;
	sfn.correspondence_file_name_ = correspondence_file;
	sfn.resfile_name_ = secondary_resfile;

	save_scalar_variable( state_name, line_number );

	state_variable_names_.insert( state_name );
	named_state_data_file_names_[ state_name ] = sfn;

}

/// @details This method reads a line beginning with the command
/// STATE_VECTOR looking for two things on this line:
/// 1. an as-of-yet unused variable name identifying the vector
/// 2. a file name which contains a list of states.
/// If the variable name is absent from the variables_ set, then
/// this method appends the name to the variables_ set and to the
/// scanner_, and appends the pair (1,2) to the input
/// strucvec_filenames parameter.  No other member variables are modified.
void
DynamicAggregateFunction::process_STATE_VECTOR_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line,
	utility::vector1< std::pair< std::string, std::string > > & strucvec_filenames
)
{
	std::string state_vector_variable_name, state_vector_filename;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read state-vector variable name in the DynamicAggregateFunction"
			" input file after reading STATE_VECTOR on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> state_vector_variable_name;
	if ( state_vector_variable_name.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read state-vector variable name in the DynamicAggregateFunction"
			" input file after reading STATE_VECTOR on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read state-vector file name in the DynamicAggregateFunction"
			" input file after reading STATE_VECTOR variable name on line " + utility::to_string( line_number ) + "\n" + line );
	}

	verify_variable_name_or_throw( state_vector_variable_name, "STATE_VECTOR", line, line_number );

	input_line >> state_vector_filename;
	if ( state_vector_filename.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read state-vector file name in the DynamicAggregateFunction"
			" input file after reading STATE_VECTOR variable name on line " + utility::to_string( line_number ) + "\n" + line );
	}

	// Success
	TR << "Read STATE_VECTOR line.  Name= " << state_vector_variable_name
		<< " File= " << state_vector_filename << std::endl;

	save_vector_variable( state_vector_variable_name, line_number );

	state_vector_variable_names_.insert( state_vector_variable_name );
	strucvec_filenames.push_back( std::make_pair( state_vector_variable_name, state_vector_filename ) );
}

void
DynamicAggregateFunction::process_POSE_ENERGY_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line
)
{
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read variable name in the DynamicAggregateFunction"
			" input file after reading POSE_ENERGY on line " + utility::to_string( line_number ) + "\n" + line );
	}
	// 1. read the new variable name that's being declared on this line
	std::string varname;
	input_line >> varname;
	if ( varname.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read variable name in the DynamicAggregateFunction"
			" input file after reading POSE_ENERGY on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a PDB name in the DynamicAggregateFunction"
			" input file after reading the variable name '" + varname + "' on line " + utility::to_string( line_number ) + "\n" + line );
	}
	verify_variable_name_or_throw( varname, "POSE_ENERGY", line, line_number );

	std::string pdb_name;
	input_line >> pdb_name;

	TR << "  Importing pose from pdb file '" << pdb_name << "'" << std::endl;
	//core::import_pose::pose_from_file( pose, pdb_name , core::import_pose::PDB_file);
	std::string pdb_string;
	try {
		pdb_string = file_contents_->get_file_contents( pdb_name );
	} catch ( utility::excn::Exception & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to open pdb file named '"
			+ pdb_name + "' given in the POSE_ENERGY command on line " + utility::to_string( line_number)
			+ "of the DynamicAggregateFunction fitness file" );
	}

	core::pose::Pose pose;
	core::import_pose::pose_from_pdbstring( pose, pdb_string, pdb_name );

	if ( pose.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Input pose given in file '"
			+ pdb_name + "' has zero residues.  Encountered while processing the '" + varname + "' variable in the DynamicAggregateFunction"
			" input file on line " + utility::to_string( line_number ) + "\n" + line );
	}

	// ok -- go ahead and score the pose
	TR << "  Scoring pose from pdb file '" << pdb_name << "'" << std::endl;
	core::Real score = (*sfxn_)( pose );
	scalar_expression_map_[ varname ] = numeric::expression_parser::VariableExpressionCOP( numeric::expression_parser::VariableExpressionOP( new VariableExpression( varname, score ) ) );
	save_scalar_variable( varname, line_number );

	TR << "Saving POSE_ENERGY of " << score << " in variable " << varname << std::endl;
}

void
DynamicAggregateFunction::process_POSE_ENERGY_VECTOR_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line
)
{
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read variable name in the DynamicAggregateFunction"
			" input file after reading POSE_ENERGY_VECTOR on line " + utility::to_string( line_number ) + "\n" + line );
	}
	// 1. read the new variable name that's being declared on this line
	std::string varname;
	input_line >> varname;
	if ( varname.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read variable name in the DynamicAggregateFunction"
			" input file after reading POSE_ENERGY_VECTOR on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read the name of a list file in the DynamicAggregateFunction"
			" input file after reading the variable name '" + varname + "' on line " + utility::to_string( line_number ) + "\n" + line );
	}
	verify_variable_name_or_throw( varname, "POSE_ENERGY_VECTOR", line, line_number );

	std::string fname;
	input_line >> fname;

	std::istringstream pdbvec_file;
	try {
		std::string fc = file_contents_->get_file_contents( fname );
		pdbvec_file.str( fc );
	} catch ( utility::excn::Exception & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to open pdb list file named '"
			+ fname + "' given in the POSE_ENERGY_VECTOR command on line " + utility::to_string( line_number)
			+ "of the DynamicAggregateFunction fitness file" );
	}

	std::list< std::string > pdb_names;
	Size count_pdbs = 0;
	Size pdbvec_linenum = 0;
	while ( pdbvec_file ) {
		++pdbvec_linenum;
		std::string line;
		std::getline( (std::istream &) (pdbvec_file), line );
		if ( line.size() == 0 ) continue;

		std::istringstream input_line( line );
		if ( input_line.peek() == '#' ) continue; // allow lines beginning with # to be commented out
		std::string pdbname;
		input_line >> pdbname;
		++count_pdbs;
		pdb_names.push_back( pdbname );
	}

	utility::vector1< VariableExpressionOP > pose_energy_variables( count_pdbs );
	count_pdbs = 0;
	for ( std::list< std::string >::const_iterator iter = pdb_names.begin(), iter_end = pdb_names.end(); iter != iter_end; ++iter ) {
		++count_pdbs;
		debug_assert( count_pdbs <= pose_energy_variables.size() );
		TR << "  Importing pose from pdb file " << *iter << std::endl;
		//core::import_pose::pose_from_file( pose, pdb_name , core::import_pose::PDB_file);
		std::string pdb_string;
		try {
			pdb_string = file_contents_->get_file_contents( *iter );
		} catch ( utility::excn::Exception & e ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to open pdb file named '"
				+ *iter + "' given in the POSE_ENERGY_VECTOR command on line " + utility::to_string( line_number)
				+ "of the DynamicAggregateFunction fitness file" );
		}
		core::pose::Pose pose;
		core::import_pose::pose_from_pdbstring( pose, pdb_string, *iter );

		if ( pose.size() == 0 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Input pose given in file '"
				+ *iter + "' has zero residues.  Encountered while processing the '" + varname + "' variable in the DynamicAggregateFunction"
				" input file on line " + utility::to_string( line_number ) + "\n" + line );
		}

		TR << "  Scoring pose from pdb file '" << *iter << std::endl;
		core::Real score = (*sfxn_)( pose );
		std::string newvar = varname + "_" + *iter;
		TR << "  Saving score of " << score << " in variable " << newvar << std::endl;
		pose_energy_variables[ count_pdbs ] = numeric::expression_parser::VariableExpressionOP( new VariableExpression( newvar, score ) );
	}

	save_vector_variable( varname, line_number );
	vector_expression_map_[ varname ] = protocols::pack_daemon::VectorExpressionCOP( protocols::pack_daemon::VectorExpressionOP( new VariableVectorExpression( varname, pose_energy_variables ) ) );


}


void
DynamicAggregateFunction::process_NPD_PROPERTY_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line
)
{
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read variable name in the DynamicAggregateFunction"
			" input file after reading NPD_PROPERTY on line " + utility::to_string( line_number ) + "\n" + line );
	}
	// 1. read the new variable name that's being declared on this line
	std::string varname;
	input_line >> varname;
	if ( varname.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read variable name in the DynamicAggregateFunction"
			" input file after reading NPD_PROPERTY on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an existing state- or state-vector "
			" variable name in the DynamicAggregateFunction"
			" input file after reading NPD_PROPERTY on line " + utility::to_string( line_number ) + "\n" + line );
	}
	verify_variable_name_or_throw( varname, "NPD_PROPERTY", line, line_number );

	// 2. Read the name of the existing state or state-vector variable
	std::string original_varname;
	input_line >> original_varname;
	if ( original_varname.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an existing state- or state-vector"
			" variable in the DynamicAggregateFunction"
			" input file after reading the new variable name on line " + utility::to_string( line_number ) + "\n" + line );
	}
	bool original_was_state_vector( false );
	if ( state_variable_names_.find( original_varname ) != state_variable_names_.end() ) {
		/// we have a state-variable, so we will add a new scalar variable to
	} else if ( state_vector_variable_names_.find( original_varname ) != state_vector_variable_names_.end() ) {
		/// OK we have a state-vector variable so we will add a new vector variable
		original_was_state_vector = true;
	} else {
		/// Error condition: the user has given a variable name that does not correspond to
		/// either a state or a state vector.
		if ( variable_names_dec_line_.find( original_varname ) != variable_names_dec_line_.end() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read the name for the existing state or state vector"
				" in the DynamicAggregateFunction input file, but instead was given the name of a non-state variable, "
				+ original_varname + " that had previously been declared on line " +
				utility::to_string( variable_names_dec_line_[ original_varname ] ) +
				".\n Error while processing line " + utility::to_string( line_number ) + "\n" + line );
		} else {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Read the name '" + original_varname + "' but "
				" expected to read the variable name for the existing state- or state-vector"
				" variable in the DynamicAggregateFunction"
				" input file after reading the new variable name on line " + utility::to_string( line_number ) + "\n" + line );
		}
	}

	// 3. Read the NPD property that will be calculated
	std::string npd_property;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a non-pairwise decomposable property"
			" name after reading the existing state- or state-vector"
			" variable name in the DynamicAggregateFunction"
			" input file on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> npd_property;
	if ( npd_property.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a non-pairwise decomposable property"
			" name after reading the existing state- or state-vector"
			" variable name in the DynamicAggregateFunction"
			" input file on line " + utility::to_string( line_number ) + "\n" + line );
	}
	// TO DO. Error check: is this a supported property?

	if ( original_was_state_vector ) {
		save_vector_variable( varname, line_number );
	} else {
		save_scalar_variable( varname, line_number );
	}

	npd_properties_for_state_variables_[ original_varname ].push_back( std::make_pair( npd_property, varname ) );
}

void
DynamicAggregateFunction::process_VECTOR_VARIABLE_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line,
	std::map< std::string, std::list< std::string > > & vector_variables
)
{
	std::string vector_variable_name, equals_sign, scalar_variable_name;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read vector variable name in the DynamicAggregateFunction"
			" input file after reading VECTOR_VARIABLE on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> vector_variable_name;
	if ( vector_variable_name.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read vector variable name in the DynamicAggregateFunction"
			" input file after reading VECTOR_VARIABLE on line " + utility::to_string( line_number ) + "\n" + line );
	}

	verify_variable_name_or_throw( vector_variable_name, "VECTOR_VARIABLE", line, line_number );

	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read equals sign in the DynamicAggregateFunction"
			" input file after reading " + vector_variable_name + " varname in the VECTOR_VARIABLE command on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> equals_sign;
	if ( equals_sign != "=" ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read equals sign in the DynamicAggregateFunction"
			" input file after reading " + vector_variable_name + " varname, but found '" + equals_sign + "' instead in the VECTOR_VARIABLE command on line " + utility::to_string( line_number ) + "\n" + line );
	}

	std::list< std::string > varname_list;
	while ( input_line ) {
		input_line >> scalar_variable_name;
		if ( !input_line ) continue;
		if ( scalar_variable_names_dec_line_.find( scalar_variable_name ) == scalar_variable_names_dec_line_.end() ) {
			if ( vector_variable_names_dec_line_.find( scalar_variable_name ) != vector_variable_names_dec_line_.end() ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Variable '" + scalar_variable_name + "' is a vector variable and"
					" cannot be listed in the scalar-variable list in a VECTOR_VARIABLE command.\n"
					"Error processing VECTOR_VARIABLE command on line " + utility::to_string( line_number ) + "\n" + line );

			} else {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Unknown variable '" + scalar_variable_name + "' requested"
					" in the VECTOR_VARIABLE command on line " + utility::to_string( line_number ) + "\n" + line );
			}
		} else {
			varname_list.push_back( scalar_variable_name );
		}
	}

	save_vector_variable( vector_variable_name, line_number );
	vector_variables[ vector_variable_name ] = varname_list;
}


/// @details This methods reads the line beginning with the command SCALAR_EXPRESSION
/// which should be in the following format ( items in <> are explained below)
/// SCALAR_EXPRESSION <varname> = <expression>
/// 1. varname must be an as-of-yet unused varible name.
/// 2. The form of expression must return a value (it should not be a vector expression).  It
///    may contain functions and variables so long as those variables have been declared
///    earlier in the file.
/// After the expression has been tokenized, it's abstract syntax tree (ast) is created and
/// appended to the input parameter scalar_expression_asts in tandem with the variable name.
/// The variable name is appended to the variables_ and scanner_ member variables.  No other
/// member variables are modified in this function.
void
DynamicAggregateFunction::process_SCALAR_EXPRESSION_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line,
	std::map< std::string, ArithmeticASTExpressionOP > & scalar_expression_asts
)
{
	std::string vname, equals_sign;//, restofline;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read variable name in the DynamicAggregateFunction"
			" input file after reading SCALAR_EXPRESSION command on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> vname;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read equals sign in the DynamicAggregateFunction"
			" input file after reading SCALAR_EXPRESSION variable name on line " + utility::to_string( line_number ) + "\n" + line );
	}

	verify_variable_name_or_throw( vname, "SCALAR_EXPRESSION", line, line_number );

	input_line >> equals_sign;
	if ( equals_sign != "=" ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read equals sign in the DynamicAggregateFunction"
			" input file after reading SCALAR_EXPRESSION variable name on line " + utility::to_string( line_number )
			+ "\n" + line + "\nbut encountered " + equals_sign + " instead." );
	}
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an expression in the DynamicAggregateFunction"
			" input file after reading SCALAR_EXPRESSION equals sign on line " + utility::to_string( line_number ) + "\n" + line );
	}

	std::string rest_of_line;
	std::getline( input_line, rest_of_line );

	TR << "On line " << line_number << ", attempting to tokenize scalar expression: " << rest_of_line << std::endl;
	TokenSetOP tokens = scanner_->scan( rest_of_line );
	TR << "On line " << line_number << ", attempting to parse scalar expression: " << rest_of_line << std::endl;
	ArithmeticASTExpressionOP expression_ast( new ArithmeticASTExpression );
	expression_ast->parse( *tokens );

	// Assuming we made it here, the AST has been correctly parsed.
	save_scalar_variable( vname, line_number );
	scalar_expression_asts[ vname ] = expression_ast;
	expression_evaluation_order_by_name_.emplace_back(std::make_pair( 1, vname ));
	TR << "SCALAR_EXPRESSION on line " << line_number << " successfully parsed" << std::endl;
}

void
DynamicAggregateFunction::process_VECTOR_EXPRESSION_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line,
	std::map< std::string, std::pair< std::map< std::string, std::string >, ArithmeticASTExpressionOP > > & vector_expression_asts
)
{
	std::string for_string, local_var, in_string, existing_vector_varname,
		comma_or_colon_string, vector_varname, equals_sign, rest_of_line;

	std::map< std::string, std::string > localvar_map;
	/// make a local copy and add variables to this local copy. Local variables declared in
	/// a single VECTOR_EXPRESSION line do not cary ovr into other lines.
	ArithmeticScanner local_scanner( *scanner_ );

	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read FOR in the DynamicAggregateFunction"
			" input file after reading VECTOR_EXPRESSION command on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> for_string;
	if ( for_string != "FOR" ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read FOR in the DynamicAggregateFunction"
			" input file after reading VECTOR_EXPRESSION command, but read '" + for_string + "' on line "
			+ utility::to_string( line_number ) + "\n" + line );
	}
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read local variable name in the DynamicAggregateFunction"
			" input file after reading FOR on a VECTOR_EXPRESSION command on line "
			+ utility::to_string( line_number ) + "\n" + line );
	}
	bool colon_found = false;
	while ( ! colon_found ) {
		if ( ! input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read local variable name in the DynamicAggregateFunction input file"
				" while processing multiple local variables on a VECTOR_EXPRESSION command on line "
				+ utility::to_string( line_number ) + "\n" + line );
		}
		input_line >> local_var;

		verify_variable_name_or_throw( local_var, "VECTOR_EXPRESSION", line, line_number );

		if ( localvar_map.find( local_var ) != localvar_map.end() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Local variable, '" + local_var + "' in VECTOR_EXPRESSION command appears multiple times on line "
				+ utility::to_string( line_number ) + "\n" + line );
		}

		if ( ! input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read 'IN' in the DynamicAggregateFunction input file"
				" following the declaration of local variable '" + local_var + "' on a VECTOR_EXPRESSION command on line "
				+ utility::to_string( line_number ) + "\n" + line );
		}
		input_line >> in_string;
		if ( in_string != "IN" ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read 'IN' in the DynamicAggregateFunction input file"
				" following the declaration of local variable '" + local_var + "', but read '" + in_string + "' in the VECTOR_EXPRESSION command on line "
				+ utility::to_string( line_number ) + "\n" + line );
		}
		if ( ! input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a vector-variable name in the DynamicAggregateFunction input file"
				" following '" + local_var + " IN' in the VECTOR_EXPRESSION command on line "
				+ utility::to_string( line_number ) + "\n" + line );
		}
		input_line >> existing_vector_varname;

		if ( vector_variable_names_dec_line_.find( existing_vector_varname ) == vector_variable_names_dec_line_.end() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "The vector-variable name '" + existing_vector_varname
				+ "'  given for local-variable '" + local_var + "' does not belong to an already-declared vector variable"
				+ ".  Error in the VECTOR_EXPRESSION command on line "
				+ utility::to_string( line_number ) + "\n" + line );
		}
		localvar_map[ local_var ] = existing_vector_varname;
		local_scanner.add_variable( local_var );
		if ( ! input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a colon or a comma following"
				" the declaration of local variable " + local_var + " in the VECTOR_EXPRESSION command on line "
				+ utility::to_string( line_number ) + "\n" + line );
		}
		input_line >> comma_or_colon_string;
		if ( comma_or_colon_string == "," ) {
			continue;
		}
		if ( comma_or_colon_string == ":" ) {
			break;
		}
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a colon or a comma following"
			" the declaration of local variable " + local_var + ", but found '" + comma_or_colon_string + "' instead."
			"\nError in the VECTOR_EXPRESSION command on line "
			+ utility::to_string( line_number ) + "\n" + line );
	}
	/// OK -- we've gotten this far.  Now we can read the vector-variable name, the equals sign, and the
	/// rest of the expression describing how to compute that vector variable.

	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a variable name"
			" following the colon in the VECTOR_EXPRESSION command on line "
			+ utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> vector_varname;

	/// Now check, do we have a valid variable name?
	verify_variable_name_or_throw( vector_varname, "VECTOR_EXPRESSION", line, line_number );
	if ( localvar_map.find( vector_varname ) != localvar_map.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Vector variable, '" + vector_varname + "' in VECTOR_EXPRESSION command appears first as a local variable\n.Line "
			+ utility::to_string( line_number ) + "\n" + line );
	}


	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an equals sign after reading variable-vector name '"
			+ vector_varname + "' in the VECTOR_EXPRESSION command on line "
			+ utility::to_string( line_number ) + "\n" + line );
	}

	input_line >> equals_sign;
	if ( equals_sign != "=" ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an equals sign after reading variable-vector name '"
			+ vector_varname + "' in the VECTOR_EXPRESSION command on line "
			+ utility::to_string( line_number ) + " but found '" + equals_sign + "' instead\n" + line );
	}
	std::getline( input_line, rest_of_line );
	TR << "On line " << line_number << ", attempting to tokenize VECTOR_EXPRESSION expression: " << rest_of_line << std::endl;
	TokenSetOP tokens = local_scanner.scan( rest_of_line );
	TR << "On line " << line_number << ", attempting to parse expression: " << rest_of_line << std::endl;
	ArithmeticASTExpressionOP vector_expression_ast( new ArithmeticASTExpression );
	vector_expression_ast->parse( *tokens );
	TR << "VECTOR_EXPRESSION on line " << line_number << " successfully parsed" << std::endl;

	save_vector_variable( vector_varname, line_number );
	expression_evaluation_order_by_name_.emplace_back(std::make_pair( 2, vector_varname ));
	vector_expression_asts[ vector_varname ] = std::make_pair( localvar_map, vector_expression_ast );
}

void
DynamicAggregateFunction::process_ENTITY_FUNCTION_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line
)
{
	std::string entityfunc_name, entityfunc_file;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read state name in the DynamicAggregateFunction"
			" input file after reading ENTITY_FUNCTION command on line " + utility::to_string( line_number ) + "\n" + line );
	}
	input_line >> entityfunc_name;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to EnitytFunc file name in the DynamicAggregateFunction"
			" input file after reading ENTITY_FUNCTION command on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( entityfunc_name.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read EntityFunc variable name in the DynamicAggregateFunction"
			" input file after reading ENTITY_FUNCTION command on line " + utility::to_string( line_number ) + "\n" + line );
	}

	verify_variable_name_or_throw( entityfunc_name, "ENTITY_FUNCTION", line, line_number );

	input_line >> entityfunc_file;
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read correspondence file name in the DynamicAggregateFunction"
			" input file after reading ENTITY_FUNCTION pdb file on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( entityfunc_file.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read pdb file name in the DynamicAggregateFunction"
			" input file after reading ENTITY_FUNCTION variable named " + entityfunc_name +
			" on line " + utility::to_string( line_number ) + "\n" + line );
	}

	/// Double check that the number of entity elements has been set
	if ( num_entity_elements_ == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "DynamicAggregateFunction::set_num_entity_elements must be called before EntityFuncs may be created" );
	}


	EntityFuncOP entfunc( new EntityFunc );
	entfunc->set_num_entity_elements( num_entity_elements_ );
	std::string entfunc_contents;

	try {
		entfunc_contents = file_contents_->get_file_contents( entityfunc_file );
	} catch ( utility::excn::Exception & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to open the entity function file file named '"
			+ entityfunc_file + "' given in the ENTITY_FUNCTION command on line " + utility::to_string( line_number)
			+ "of the DynamicAggregateFunction fitness file" );
	}

	std::istringstream iss( entfunc_contents );
	entfunc->initialize_from_input_file( iss );

	save_scalar_variable( entityfunc_name, line_number );
	entity_funcs_[ entityfunc_name ] = std::make_pair( entfunc, SurrogateVariableExpressionOP( new SurrogateVariableExpression( entityfunc_name ) ) );
	entity_funcs_dec_line_[ entityfunc_name ] = line_number;
}


/// @details This method reads the line of the input file beginning with the command FITNESS.
/// Each DynamicAggregateFunction input file should contain exactly one FITNESS command.
/// The fitness line should be of the form:
/// FITNESS <expression>
/// where the expression must return a value (it should not be a vector expression) and
/// may refer to functions and variables so long as those variables have been declared
/// earlier in the file.  The abstract syntax tree (ast) for the fitness expression
/// is returned in the input parameter, fitness_expression_ast.
void
DynamicAggregateFunction::process_FITNESS_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line,
	ArithmeticASTExpressionOP & fitness_expression_ast
)
{
	if ( fitness_expression_ast ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "FITNESS command appears multiple"
			" times in the DynamicAggregateFunction file.\nSecond occurrance found while reading\n"
			+ line + "\nLine # " + utility::to_string( line_number )  );
	}
	std::string rest_of_line;
	std::getline( input_line, rest_of_line );
	TR << "On line " << line_number << ", attempting to tokenize FITNESS expression: " << rest_of_line << std::endl;
	TokenSetOP tokens = scanner_->scan( rest_of_line );
	TR << "On line " << line_number << ", attempting to parse expression: " << rest_of_line << std::endl;
	fitness_expression_ast = ArithmeticASTExpressionOP( new ArithmeticASTExpression );
	fitness_expression_ast->parse( *tokens );
	TR << "FITNESS on line " << line_number << " successfully parsed" << std::endl;
}

void
DynamicAggregateFunction::save_scalar_variable( std::string const & varname, core::Size line_number )
{
	scanner_->add_variable( varname );
	variable_names_dec_line_[ varname ] = line_number;
	scalar_variable_names_dec_line_[ varname ] = line_number;
}

void
DynamicAggregateFunction::save_vector_variable( std::string const & varname, core::Size line_number )
{
	scanner_->add_variable( varname );
	variable_names_dec_line_[ varname ] = line_number;
	vector_variable_names_dec_line_[ varname ] = line_number;
}


/// @details This function reads the contents of a state-vector file.  The input
/// parameter vec_varname refers to the variable name associated with the file named
/// by the second input parameter -- this pairing occurred somewhere in the
/// DynamicAggregateFunction input file in a STATE_VECTOR command.  The format
/// of the state vector should be a series of lines each with three strings.
/// String 1: the name of the pdb file containing a state that should be repacked.
/// String 2: the name of the correspondence file mapping residues in the input pdb
///           to the shared residues.
/// String 3: the name of the secondary resfile describing how the packer should treat
///           residues other than the correspondence residues in the protein.
/// Lines beginning with "#" are ignored.
/// The triples read from this file are appended to the state_vector_data_file_names_
/// member variable.  The method increments the input variable n_vector_states for
/// each triple that it reads from the state-vector file.
void
DynamicAggregateFunction::read_state_vector_file(
	std::string const & vec_varname,
	std::string const & fname,
	Size & n_vector_states
)
{
	std::istringstream strucvec_file;
	try {
		std::string fc = file_contents_->get_file_contents( fname );
		strucvec_file.str( fc );
	} catch ( utility::excn::Exception & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to open state vector file named '" + fname + "' which was listed on line " + utility::to_string( variable_names_dec_line_[ vec_varname ] ) );
	}
	Size svline_num( 0 );
	utility::vector1< StructureFileNames > state_triples;
	while ( strucvec_file ) {
		++svline_num;
		std::string line;
		std::getline( (std::istream &) ( strucvec_file ), line );
		if ( line.size() == 0 ) continue;

		std::istringstream input_line( line );
		if ( input_line.peek() == '#' ) continue;

		std::string structure_file, correspondence_file, secondary_resfile;
		if ( ! input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read pdb file name as the first argument on line "
				+ utility::to_string( svline_num ) + " of the state-vector file" + fname + "\n" + line );
		}
		input_line >> structure_file;
		if ( ! input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read correspondence file name as the second argument on line "
				+ utility::to_string( svline_num ) + " of the state-vector file" + fname + "\n" + line );
		}
		input_line >> correspondence_file;
		if ( ! input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read secondary resfile name as the third argument on line "
				+ utility::to_string( svline_num ) + " of the state-vector file" + fname + "\n" + line );
		}
		input_line >> secondary_resfile;

		StructureFileNames sfn;
		sfn.pdb_name_ = structure_file;
		sfn.correspondence_file_name_ = correspondence_file;
		sfn.resfile_name_ = secondary_resfile;
		state_triples.push_back( sfn );
		TR << "Read state from " << fname << " PDB= " << structure_file
			<< " Corr= " << correspondence_file << " 2RF= " << secondary_resfile << std::endl;
		++n_vector_states;
	}
	if ( state_triples.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to find any states in state-vector file " + fname );
	}
	state_vector_data_file_names_[ vec_varname ] = state_triples;

}

/// @details For each element of the named_state_data_file_names_ data member,
/// creates a VariableExpression (named with the variable name given in the STATE
/// command from the DynamicAggregateFunction input file).  This variable expression
/// is put into three member variables:
/// the variable_expressions_for_states_ vector,
/// the variable_expressions_ vector, and
/// the named_state_expression_map_.
/// This method increments the count_states input parameter for each VariableExpression
/// that it processes, pairing the value of the count_states variable with the index
/// for this VariableExpression in both the variable_expressions_for_states_ vector and
/// in the variable_expressions_ vector.
void
DynamicAggregateFunction::create_state_variable_expressions(
	Size & count_state,
	Size & count_npd_index,
	Size & count_variable_index
)
{
	for ( std::map< std::string, StructureFileNames >::const_iterator
			iter = named_state_data_file_names_.begin(), iter_end = named_state_data_file_names_.end();
			iter != iter_end; ++iter ) {
		++count_state;
		++count_variable_index;
		variable_expressions_for_states_[ count_state ] = numeric::expression_parser::VariableExpressionOP( new VariableExpression( iter->first, 0.0 ) );
		variable_expressions_[ count_variable_index ] = variable_expressions_for_states_[ count_state ];
		files_for_state_[ count_state ] = iter->second;
		named_state_expression_map_[ iter->first ] = variable_expressions_for_states_[ count_state ];
		scalar_expression_map_[ iter->first ] = variable_expressions_for_states_[ count_state ];
		state_variable_name_2_state_index_[ iter->first ] = count_state;

		if ( npd_properties_for_state_variables_.find( iter->first ) != npd_properties_for_state_variables_.end() ) {
			std::list< std::pair< std::string, std::string > > const & npdlist( npd_properties_for_state_variables_[ iter->first ]);
			for ( auto const & npditer : npdlist ) {
				++count_npd_index;
				++count_variable_index;
				SurrogateVariableExpressionOP surrogate( new SurrogateVariableExpression( npditer.second, 0.0 ) );
				surrogate->root_expression( variable_expressions_for_states_[ count_state ] );

				npd_variable_indices_for_states_[ count_state ].push_back( std::make_pair( count_npd_index, npditer.first ) );
				variable_expressions_for_npd_properties_[ count_npd_index ] =
					variable_expressions_[ count_variable_index ] = surrogate;
				variable_name_2_variable_exp_index_[ npditer.second ] = count_variable_index;
				scalar_expression_map_[ npditer.second ] = variable_expressions_[ count_variable_index ];
			}
		}

	}
}

/// @details This method iterates across state_vector_data_file_names_ map,
/// and creates a VariableVectorExpression for each entry.  It adds this VariableVectorExpression
/// to the state_vector_variables_ map.  For each element in the state_vector_data_file_names_
/// map, and for each state-file-triple held in each particular element, it creates
/// a VariableExpression and pairs the index in the input count_state vector.  This
/// VariableExpression is added to three member variables:
/// the variable_expressions_for_states_ vector,
/// the variable_expressions_ vector.
/// Each VariableVectorExpression is given a vector of its set of VariableExpressions
/// in its constructor.  It will access these variables directly when it is queried for
/// its vector_values.  The indices of the state varaibles are inserted into the
/// state_indices_for_state_vector_ map.
void
DynamicAggregateFunction::create_variable_vector_expressions(
	Size & count_state,
	Size & count_npd_index,
	Size & count_variable_index
)
{
	for ( std::map< std::string, utility::vector1< StructureFileNames > >::const_iterator
			iter = state_vector_data_file_names_.begin(), iter_end = state_vector_data_file_names_.end();
			iter != iter_end; ++iter ) {
		utility::vector1< Size > indices( iter->second.size() );
		utility::vector1< VariableExpressionCOP > variables( iter->second.size() );
		std::map< std::string, utility::vector1< VariableExpressionCOP > > npd_property_variables;
		bool has_npd_properties = npd_properties_for_state_variables_.find( iter->first ) != npd_properties_for_state_variables_.end();

		for ( Size ii = 1; ii <= iter->second.size(); ++ii ) {
			++count_state;
			++count_variable_index;
			std::string ii_varname( iter->first + "_" + utility::to_string( ii ) );
			TR << "Adding state " << ii_varname << " with state index " << count_state << std::endl;
			variable_expressions_for_states_[ count_state ] = numeric::expression_parser::VariableExpressionOP( new VariableExpression( ii_varname, 0.0 ) );
			indices[ ii ] = count_state;
			variable_expressions_[ count_variable_index ] = variable_expressions_for_states_[ count_state ];
			files_for_state_[ count_state ] = iter->second[ ii ];
			variables[ ii ] = variable_expressions_for_states_[ count_state ];
			state_variable_name_2_state_index_[ ii_varname ] = count_state;
			variable_name_2_variable_exp_index_[ ii_varname ] = count_variable_index;

			/// Create vector expressions for this state and store them in the variable_expressions_for_npd_properties_ array
			if ( has_npd_properties ) {
				std::list< std::pair< std::string, std::string > > const & npdlist( npd_properties_for_state_variables_[ iter->first ]);
				for ( auto const & npditer : npdlist ) {
					++count_npd_index;
					++count_variable_index;
					npd_variable_indices_for_states_[ count_state ].push_back( std::make_pair( count_npd_index, npditer.first ) );
					std::string ii_npd_varname = npditer.second + "_" + utility::to_string( ii );
					SurrogateVariableExpressionOP surrogate( new SurrogateVariableExpression( ii_npd_varname, 0.0 ) );
					surrogate->root_expression( variable_expressions_for_states_[ count_state ] );
					variable_expressions_for_npd_properties_[ count_npd_index ] =
						variable_expressions_[ count_variable_index ] = surrogate;
					npd_property_variables[ npditer.first ].push_back( variable_expressions_[ count_variable_index ] );
					variable_name_2_variable_exp_index_[ ii_npd_varname ] = count_variable_index;
				}
			}
		}
		vector_expression_map_[ iter->first ] = state_vector_variables_[ iter->first ] = protocols::pack_daemon::VariableVectorExpressionOP( new VariableVectorExpression( iter->first, variables ) );
		if ( has_npd_properties ) {
			std::list< std::pair< std::string, std::string > > const & npdlist( npd_properties_for_state_variables_[ iter->first ]);
			for ( auto const & npditer : npdlist ) {
				//std::string ii_npd_varname = npditer->second + "_" + utility::to_string( ii );
				//npd_property_variables[ npditer->first ].push_back( new VariableExpressin( ii_npd_varname, 0.0 ) );
				vector_expression_map_[ npditer.second ] = protocols::pack_daemon::VectorExpressionCOP( protocols::pack_daemon::VectorExpressionOP( new VariableVectorExpression( npditer.second, npd_property_variables[ npditer.first ] ) ) );
			}
		}
		state_indices_for_state_vector_[ iter->first ] = indices;
	}

}

/// @details creates a variable expression for each sub-expression that will hold the
/// intermediate result of a sub-expression.  Increments the count_variable_index input parameter
/// and places the VariableExpression into the variable_expressions_ array into position
/// count_variable_index.
void
DynamicAggregateFunction::create_scalar_and_vector_expression_variable_expressions(
	std::map< std::string, ArithmeticASTExpressionOP > const & scalar_expression_asts,
	std::map< std::string, std::list< std::string > > const & vector_variables,
	Size & count_variable_index
)
{
	for ( auto const & scalar_expression_ast : scalar_expression_asts ) {
		++count_variable_index;
		SurrogateVariableExpressionOP surrogate( new SurrogateVariableExpression( scalar_expression_ast.first, 0.0 ) );
		variable_expressions_[ count_variable_index ] = surrogate;
		surrogate_expression_map_[ count_variable_index ] = surrogate;
		variable_name_2_variable_exp_index_[ scalar_expression_ast.first ] = count_variable_index;
		scalar_expression_map_[ scalar_expression_ast.first ] = variable_expressions_[ count_variable_index ];
	}

	for ( auto const & vector_variable : vector_variables ) {
		utility::vector1< VariableExpressionCOP > variables;
		variables.reserve( vector_variable.second.size() );
		for ( auto
				scalvar_iter = vector_variable.second.begin(), scalvar_iter_end = vector_variable.second.end();
				scalvar_iter != scalvar_iter_end; ++scalvar_iter ) {
			std::map< std::string, VariableExpressionCOP >::const_iterator scvar = scalar_expression_map_.find( *scalvar_iter );
			if ( scvar == scalar_expression_map_.end() ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Internal error: could not locate requested scalar variable '"
					+ *scalvar_iter + "' when constructing the VECTOR_VARIABLE '" + vector_variable.first );
			}
			variables.push_back( scvar->second );
		}
		vector_expression_map_[ vector_variable.first ] = protocols::pack_daemon::VectorExpressionCOP( protocols::pack_daemon::VectorExpressionOP( new VariableVectorExpression( vector_variable.first, variables ) ) );
	}
}

/// @details Turns expression abstract-syntax trees into Expression objects.  The
/// VectorExpressionCreator keeps a reference to me so that it can pass control of
/// flow to me when processing variable- and function-construction events.
/// The scalar_expressions_ array is updated to hold the ExpressionOPs for these sub-expressions
/// as pairs each Expression with the corresponding index of its VariableExpression
/// so that I can later update that VariableExpression (held in the variable_expressions_ array)
/// when evaluating that sub-expression.
/// The fitness_exp_ is assigned to the Expression coming from the fitness_expression_ast.
void
DynamicAggregateFunction::turn_expression_ASTs_into_expressions(
	std::map< std::string, ArithmeticASTExpressionOP > const & scalar_expression_asts,
	std::map< std::string, std::pair< std::map< std::string, std::string >, ArithmeticASTExpressionOP > >  const & vector_expression_asts,
	ArithmeticASTExpressionOP fitness_expression_ast
)
{
	VectorExpressionCreator expression_creator( *this );
	// start counting from the number of states (states) -- the sub-expressions start indexing
	// at num-states + 1.
	//Size count_variable_index = num_states();
	scalar_expressions_.clear();

	ASTPrinter printer;
	//printer.pretty( false );

	for ( std::list< std::pair< Size, std::string > >::const_iterator
			iter = expression_evaluation_order_by_name_.begin(), iter_end = expression_evaluation_order_by_name_.end();
			iter != iter_end; ++iter ) {
		if ( iter->first == 1 ) {
			TR << "Creating scalar expression for " << iter->second << std::endl;
			auto ast_iter = scalar_expression_asts.find( iter->second );
			if ( ast_iter == scalar_expression_asts.end() ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Internal error.  Unable to find scalar expression named '" + iter->second + "' in the scalar_expression_asts map" );
			}
			TR << "Creating expression from AST:\n" << printer.ast_string( *(ast_iter->second) ) << std::endl;

			ExpressionCOP var_expr = expression_creator.create_expression_tree( *(ast_iter->second) );
			Size const variable_index = variable_name_2_variable_exp_index_[ iter->second ];
			//++count_variable_index;
			scalar_expressions_.push_back( std::make_pair( variable_index, var_expr ));
			surrogate_expression_map_[ variable_index ]->root_expression( var_expr );
		} else if ( iter->first == 2 ) {
			TR << "Creating vector expression for " << iter->second << std::endl;
			auto
				ast_iter = vector_expression_asts.find( iter->second );
			if ( ast_iter == vector_expression_asts.end() ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Internal error.  Unable to find vector expression named '" + iter->second + "' in the vector_expression_asts map" );
			}

			std::pair< std::map< std::string, std::string >, ArithmeticASTExpressionOP > const & itvecexp_data =
				ast_iter->second;

			/// Create a mapping between the local-variable names and the vector-expression COPs
			/// that the local variables are intended to refer to.
			std::map< std::string, VectorExpressionCOP > vector_varnames;
			for ( auto const & vec_iter : itvecexp_data.first ) {
				//std::cout << "VectorExpression map: " << vec_iter->first << " " << vec_iter->second  << std::endl;
				if ( vector_expression_map_.find( vec_iter.second ) == vector_expression_map_.end() ) {
					throw CREATE_EXCEPTION(utility::excn::Exception,  "Variable " + vec_iter.second + " absent from the IterativeVectorExpression name map" );
				}
				vector_varnames[ vec_iter.first ] = vector_expression_map_[ vec_iter.second ];
			}
			IterativeVectorExpressionOP ivec_exp( new IterativeVectorExpression( iter->second ) );

			TR << "Creating expression from AST:\n" << printer.ast_string( *itvecexp_data.second ) << std::endl;

			/// Before initializing the iterative vector expression, I must hold a pointer to it.
			/// This allows me to pass control of flow to the iterative vector expression when I
			/// am asked to process variable expressions from the expression_creator.
			focused_iterative_vector_expression_ = ivec_exp;
			ivec_exp->initialize( vector_varnames, *itvecexp_data.second, expression_creator );
			focused_iterative_vector_expression_.reset();

			vector_expression_map_[ iter->second ] = ivec_exp;
		}
	}

	TR << "Creating fitness expression from AST:\n" << printer.ast_string( *fitness_expression_ast ) << std::endl;
	fitness_exp_ = expression_creator.create_expression_tree( *fitness_expression_ast );

}


/// @details Verifies that the args array for a vector function has size 1, and that the sole
/// element is a VectorExpression.  Returns an owning-pointer to the downcasted VectorExpression,
/// or throws an exception if the downcast fails.
utility::vector1< VectorExpressionCOP >
DynamicAggregateFunction::verify_vector_arguments(
	std::string const & fname,
	utility::vector1< ExpressionCOP > const & args,
	Size expected_nargs
) const
{
	if ( args.size() != expected_nargs ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "vector function expression " + fname + " construction requested with nargs != 1. Nargs= " + utility::to_string( args.size() )  );
	}
	utility::vector1< VectorExpressionCOP > vector_expressions( expected_nargs );
	for ( Size ii = 1; ii <= expected_nargs; ++ii ) {
		VectorExpressionCOP vec_ptr = utility::pointer::dynamic_pointer_cast< protocols::pack_daemon::VectorExpression const > ( args[ ii ] );
		if ( ! vec_ptr ) {
			VariableExpressionCOP var_ptr = utility::pointer::dynamic_pointer_cast< VariableExpression const > ( args[ ii ] );
			if ( var_ptr ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "vector function expression " + fname + " can only be constructed from a vector expression.\n"
					"Variable " + var_ptr->name() + " is not a vector variable." );
			} else {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "vector function expression " + fname + " can only be constructed from a vector expression." );
			}
		}
		vector_expressions[ ii ] = vec_ptr;
	}
	return vector_expressions;
}

void
DynamicAggregateFunction::verify_variable_name_or_throw(
	std::string const & vname,
	std::string const & command_name,
	std::string const & line,
	Size line_number
)
{
	if ( variable_names_dec_line_.find( vname ) != variable_names_dec_line_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Variable name " + vname + " appears multiple"
			" times in the DynamicAggregateFunction file.\nSecond occurrance found while reading a " +
			command_name + " command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}
	if ( function_names_.find( vname ) != function_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Declaration of variable '" + vname + "' in " + command_name + " command conflicts with a function name.  Line "
			+ utility::to_string( line_number ) + "\n" + line );
	}
	if ( entity_funcs_.find( vname ) != entity_funcs_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Declaration of variable '" + vname + "' in " + command_name + " command conflicts"
			" with a previously declared EntityFunc variable.\nPrevious declaration was on line "
			+ utility::to_string( entity_funcs_dec_line_[ vname ] ) + ".\nLine "
			+ utility::to_string( line_number ) + "\n" + line );

	}
	if ( illegal_variable_names_.find( vname ) != illegal_variable_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal name for variable, '" + vname + "' in the " +
			command_name + " command on line " + utility::to_string( line_number ) + "\n" + line );
	}

}

void
DynamicAggregateFunction::count_file_reads()
{
	file_contents_->delete_contents_at_nread_limit( true );
	file_contents_->refuse_unexpected_files( true );

	for ( std::map< std::string, StructureFileNames >::const_iterator
			iter = named_state_data_file_names_.begin(), iter_end = named_state_data_file_names_.end();
			iter != iter_end; ++iter ) {
		StructureFileNames const & sfn = iter->second;
		file_contents_->increment_nread_limit( sfn.pdb_name_ );
		file_contents_->increment_nread_limit( sfn.correspondence_file_name_ );
		file_contents_->increment_nread_limit( sfn.resfile_name_ );
	}
	for ( std::map< std::string, utility::vector1< StructureFileNames > >::const_iterator
			iter = state_vector_data_file_names_.begin(), iter_end = state_vector_data_file_names_.end();
			iter != iter_end; ++iter ) {
		for ( Size ii = 1; ii <= iter->second.size(); ++ii ) {
			StructureFileNames const & sfn = iter->second[ ii ];
			file_contents_->increment_nread_limit( sfn.pdb_name_ );
			file_contents_->increment_nread_limit( sfn.correspondence_file_name_ );
			file_contents_->increment_nread_limit( sfn.resfile_name_ );
		}
	}

}

std::string
DynamicAggregateFunction::get_file_contents(
	std::string const & filename
)
{
	return file_contents_->get_file_contents( filename );
}


void
DynamicAggregateFunction::assign_state_energies_to_variables_and_subexpressions(
	StateEnergies const & state_energies,
	StateEnergies const & npd_properties,
	Entity const & entity,
	bool verbose
)
{
	for ( std::map< std::string, std::pair< EntityFuncOP, VariableExpressionOP > >::const_iterator
			iter = entity_funcs_.begin(), iter_end = entity_funcs_.end();
			iter != iter_end; ++iter ) {
		core::Real entity_func_value = iter->second.first->evaluate( entity, verbose );
		iter->second.second->set_value( entity_func_value );
		if ( verbose ) {
			TR << "EntityFunc:: total for " << iter->second.second->name() << " " << entity_func_value << std::endl;
		}
	}

	for ( Size ii = 1; ii <= num_states(); ++ii ) {
		if ( verbose ) {
			TR << "Assigning value " << state_energies[ ii ] << " to VariableExpression " << ii << std::endl;
		}
		variable_expressions_for_states_[ ii ]->set_value( state_energies[ ii ] );
	}
	for ( Size ii = 1; ii <= npd_properties.size(); ++ii ) {
		variable_expressions_for_npd_properties_[ ii ]->set_value( npd_properties[ ii ] );
	}

	//TR << "Finished variable expression assignment " << std::endl;
	for ( Size ii = 1; ii <= scalar_expressions_.size(); ++ii ) {
		//TR << "Sub expression index " << ii << std::endl;
		//TR << "At " << scalar_expressions_[ ii ].second() << std::endl;
		//TR << "Evaluating sub_expression " << scalar_expressions_[ ii ].second->name() << std::endl;
		variable_expressions_[ scalar_expressions_[ ii ].first ]->set_value( (*(scalar_expressions_[ ii ].second))() );
		if ( verbose ) {
			TR << "sub expression " << variable_expressions_[ scalar_expressions_[ ii ].first ]->name()
				<< " evaluated to " << (*variable_expressions_[ scalar_expressions_[ ii ].first ]) () << std::endl;
		}
	}
}


core::Size
DynamicAggregateFunction::count_num_npd_properties() const
{
	Size count( 0 );
	/// 1. count the npd properties for states declared with the STATE command
	for ( auto const & named_state_data_file_name : named_state_data_file_names_ ) {
		if ( npd_properties_for_state_variables_.find( named_state_data_file_name.first ) != npd_properties_for_state_variables_.end() ) {
			count += npd_properties_for_state_variables_.find( named_state_data_file_name.first )->second.size();
		}
	}
	/// 2. count the npd properties for the states decalred with the STATE_VECTOR command
	for ( auto const & state_vector_data_file_name : state_vector_data_file_names_ ) {
		if ( npd_properties_for_state_variables_.find( state_vector_data_file_name.first ) != npd_properties_for_state_variables_.end() ) {
			count += npd_properties_for_state_variables_.find( state_vector_data_file_name.first )->second.size() * state_vector_data_file_name.second.size();
		}
	}
	return count;
}

////////////////////////////////////////////////
///////// DynamicAggregateFunctionDriver ///////
////////////////////////////////////////////////

void DynamicAggregateFunctionDriver::initialize_from_input_file(
	DaemonSetOP daemon_set,
	std::istream & input
)
{
	try {
		read_all_variables_from_input_file( input );
	} catch ( utility::excn::Exception & e ) {
		send_error_message_to_remote_daemon_sets();
		TR << "Initialization from input file failed with exception: " << e.msg() << std::endl;
		TR << "Remote daemon sets are spinning down" << std::endl;
		throw;// e;
	}
	initialize_pack_daemons( daemon_set );
}

void
DynamicAggregateFunctionDriver::initialize_pack_daemons( DaemonSetOP daemon_set )
{
	/// Look at all the files that need to be read.  If any file needs to be read
	/// multiple times, then keep the contents of that file in memory until its no longer needed.

	count_file_reads();

	/// Distribute the states across the nodes available for this job.  Some states will be assigned to this
	/// pack daemon, unless there are more nodes than states.

#ifdef USEMPI
	distribute_jobs_to_remote_daemons( daemon_set );
#else
	initialize_daemon_with_all_states( daemon_set );
#endif
	daemon_set->setup_daemons();
}

void
DynamicAggregateFunctionDriver::distribute_jobs_to_remote_daemons(
	DaemonSetOP daemon_set
)
{
	int MPI_nprocs = utility::mpi_nprocs();
	int nstructs = num_states();

	auto njobs_per_cpu = static_cast< int > ( std::ceil( (double) nstructs / MPI_nprocs ));
	/// Assign myself fewer states than to all other nodes if the number of states is not evenly divisible
	/// by the number of processors.

	int overhang = njobs_per_cpu * MPI_nprocs - nstructs;

	utility::vector0< std::list< int > > job_assignments( MPI_nprocs );
	if ( overhang == 0 ) {
		int state_count = 0;
		for ( int ii = 0; ii < MPI_nprocs; ++ii ) {
			for ( int jj = 1; jj <= njobs_per_cpu; ++jj ) {
				++state_count;
				job_assignments[ ii ].push_back( state_count );
			}
		}
	} else if ( overhang == 1 && njobs_per_cpu == 1 ) {
		int state_count = 0;
		for ( int ii = 1; ii < MPI_nprocs; ++ii ) {
			++state_count;
			job_assignments[ ii ].push_back( state_count );
		}
	} else {
		int state_count = 0;
		for ( int ii = 0; ii < MPI_nprocs; ++ii ) {
			int ii_n_jobs = ii < overhang ? njobs_per_cpu - 1 : njobs_per_cpu;
			for ( int jj = 1; jj <= ii_n_jobs; ++jj ) {
				++state_count;
				job_assignments[ ii ].push_back( state_count );
			}
		}
	}

	// Now assign jobs to nodes.
	for ( int ii = 1; ii < MPI_nprocs; ++ii ) {
		assign_jobs_to_remote_daemon_sets( ii, job_assignments[ ii ] );
	}
	bool jobs_successfully_initialized( true );

	std::string error_message;

	try {
		assign_jobs_to_local_daemon_set( job_assignments[ 0 ], daemon_set );
	} catch ( utility::excn::Exception & e ) {
		error_message += "Error from node 0\n";
		error_message += e.msg();
		TR << e.msg() << std::endl;
		jobs_successfully_initialized = false;
	}

	for ( int ii = 1; ii < MPI_nprocs; ++ii ) {
		if ( ! verify_remote_daemon_set_initialization_successful( ii ) ) {
			error_message += "Error from node " + utility::to_string( ii ) + "\n";
			std::string remote_message = utility::receive_string_from_node( ii );
			error_message += remote_message + "\n";
			TR << remote_message << std::endl;
			jobs_successfully_initialized = false;
		}
	}

	if ( jobs_successfully_initialized ) {
		for ( int ii = 1; ii < MPI_nprocs; ++ii ) {
			send_success_message_to_remote_daemon_set( ii );
		}
	} else {
		send_error_message_to_remote_daemon_sets();
#ifdef USEMPI
		MPI_Finalize();
#endif
		utility_exit_with_message( error_message );
	}
}

void
DynamicAggregateFunctionDriver::assign_jobs_to_local_daemon_set(
	std::list< int > const & job_indices,
	DaemonSetOP daemon_set
)
{
	for ( int state_id : job_indices ) {
		StructureFileNames const & sfn = file_inputs_for_job( state_id );
		core::pose::Pose pose;
		core::import_pose::pose_from_pdbstring( pose, get_file_contents( sfn.pdb_name_ ));
		std::istringstream corr_stream( get_file_contents( sfn.correspondence_file_name_ ) );
		std::istringstream resfile_stream( get_file_contents( sfn.resfile_name_ ) );
		daemon_set->add_pack_daemon(
			state_id,
			sfn.pdb_name_,
			pose,
			sfn.correspondence_file_name_,
			corr_stream,
			sfn.resfile_name_,
			resfile_stream );
		for ( auto
				npditer = npd_variable_indices_for_state_begin( state_id ),
				npditer_end = npd_variable_indices_for_state_end( state_id );
				npditer != npditer_end; ++npditer ) {
			daemon_set->add_npd_property_calculator_for_state( state_id, npditer->second, npditer->first );
		}
	}

}

void
DynamicAggregateFunctionDriver::assign_jobs_to_remote_daemon_sets(
	int proc_id,
	std::list< int > const & job_indices
)
{
	utility::send_integer_to_node( proc_id, add_daemon );

	int ndaemons = job_indices.size();
	utility::send_integer_to_node( proc_id, ndaemons );
	for ( int state_id : job_indices ) {
		StructureFileNames const & sfn = file_inputs_for_job( state_id );
		utility::send_integer_to_node( proc_id, state_id );
		utility::send_string_to_node( proc_id, sfn.pdb_name_ );
		utility::send_string_to_node( proc_id, get_file_contents( sfn.pdb_name_ ));
		utility::send_string_to_node( proc_id, sfn.correspondence_file_name_ );
		utility::send_string_to_node( proc_id, get_file_contents( sfn.correspondence_file_name_ ));
		utility::send_string_to_node( proc_id, sfn.resfile_name_ );
		utility::send_string_to_node( proc_id, get_file_contents( sfn.resfile_name_ ));
		int n_npd_properties_to_send = num_npd_properties_for_state( state_id );
		utility::send_integer_to_node( proc_id, n_npd_properties_to_send );
		//int count_npd = 0;
		for ( auto
				npditer = npd_variable_indices_for_state_begin( state_id ),
				npditer_end = npd_variable_indices_for_state_end( state_id );
				npditer != npditer_end; ++npditer ) {
			//++count_npd;
			utility::send_integer_to_node( proc_id, npditer->first );
			utility::send_string_to_node( proc_id, npditer->second );
		}
	}

}

bool DynamicAggregateFunctionDriver::verify_remote_daemon_set_initialization_successful( int proc_id ) const
{
	int message = utility::receive_integer_from_node( proc_id );
	if ( message == success_message ) return true;
	if ( message == error_message ) return false;

	throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected either a success_message or an error_message "
		"after initializing remote daemon sets from process " +
		utility::to_string( proc_id ) + " but received " + utility::to_string( message ) );
	return false;
}


void DynamicAggregateFunctionDriver::send_success_message_to_remote_daemon_set( int proc_id ) const
{
	utility::send_integer_to_node( proc_id, success_message );
}

void DynamicAggregateFunctionDriver::send_error_message_to_remote_daemon_sets() const
{
	int MPI_nprocs = utility::mpi_nprocs();
	for ( int ii = 1; ii < MPI_nprocs; ++ii ) {
		utility::send_integer_to_node( ii, error_message );
	}
}

void
DynamicAggregateFunctionDriver::initialize_daemon_with_all_states(
	DaemonSetOP daemon_set
)
{
	int nstructs = num_states();
	std::list< int > joblist;
	for ( int ii = 1; ii <= nstructs; ++ii ) joblist.push_back( ii );

	try {
		assign_jobs_to_local_daemon_set( joblist, daemon_set );
	} catch ( utility::excn::Exception & e ) {
		TR << "Error from daemon-set initialization \n";
		TR <<  e.msg();
		TR << std::endl;
		throw;// e;
	}
}


////////////////////////////////////////////////
///////// EntityFuncExpressionCreator //////////
////////////////////////////////////////////////

EntityFuncExpressionCreator::EntityFuncExpressionCreator(
	EntityFunc const & owner
) :
	owner_( owner )
{}

EntityFuncExpressionCreator::~EntityFuncExpressionCreator() = default;

ExpressionCOP
EntityFuncExpressionCreator::handle_variable_expression(
	ArithmeticASTValue const & node
)
{
	return owner_.variable_expression( node );
}

ExpressionCOP
EntityFuncExpressionCreator::handle_function_expression(
	FunctionTokenCOP function,
	utility::vector1< ExpressionCOP > const & args
)
{
	return owner_.function_expression( function, args );
}


////////////////////////////////////////////////
////////////////// Entity Funct ////////////////
////////////////////////////////////////////////

EntityFunc::EntityFunc() :
	num_entity_elements_( 0 )
{}

EntityFunc::~EntityFunc() = default;

void EntityFunc::set_num_entity_elements( Size num_ees )
{
	num_entity_elements_ = num_ees;
}

void EntityFunc::initialize_from_input_file( std::istream & input )
{
	if ( num_entity_elements_ == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "EntityFunc::initialize_from_input_file cannot be called without EntityFunc::set_num_entity_elements having first been called."  );
	}

	initialize_scanner_and_function_names();
	std::map< std::string, ArithmeticASTExpressionOP > expression_asts;
	ArithmeticASTExpressionOP score_expression_ast( nullptr );
	Size count_line( 0 );

	while ( input ) {
		++count_line;
		std::string line;
		std::getline( ( std::istream & ) input, line );
		if ( line.size() == 0 ) continue; // skip blank lines

		std::istringstream input_line( line );
		if ( input_line.peek() == '#' ) continue; // skip comment lines
		std::string command;
		input_line >> command;
		if ( command == "AA_SET" ) {
			process_AA_SET_line( line, count_line, input_line );
		} else if ( command == "SET_CONDITION" ) {
			process_SET_CONDITION_line( line, count_line, input_line );
		} else if ( command == "SUB_EXPRESSION" ) {
			process_SUB_EXPRESSION_line( line, count_line, input_line, expression_asts );
		} else if ( command == "SCORE" ) {
			process_SCORE_line( line, count_line, input_line, score_expression_ast );
		} else {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Unable to recognize command '" + command + "' while reading EntityFunc file" );
		}
	}

	if ( ! score_expression_ast ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "SCORE command was not found in EntityFunc file." );
	}
	turn_expression_ASTs_into_expressions( expression_asts, score_expression_ast );

}

core::Real
EntityFunc::evaluate( Entity const & entity, bool verbose )
{
	assign_entity_sequence_to_variables( entity );
	for ( Size ii = 1; ii <= expression_evaluation_order_.size(); ++ii ) {
		core::Real iival = (*expression_evaluation_order_[ ii ].first )();
		expression_evaluation_order_[ ii ].second->set_value( iival );
		if ( verbose ) {
			TR << "EntityFunc::evaluate " << expression_evaluation_order_[ ii ].second->name() << " " << iival << std::endl;
		}
	}
	return (*score_expression_)();
}

ExpressionCOP
EntityFunc::variable_expression( ArithmeticASTValue const & var_node ) const
{
	if ( var_node.is_literal() ) {
		utility_exit_with_message( "Error in EntityFunc::variable_expression; non-variable (literal) node recieved" +
			utility::to_string( var_node.literal_value() ));
	}

	auto var_iter =
		variable_expression_map_.find( var_node.variable_name() );
	if ( var_iter != variable_expression_map_.end() ) {
		return var_iter->second;
	}

	throw CREATE_EXCEPTION(utility::excn::Exception,  "Unable to find variable with name " + var_node.variable_name() + " in"
		" EntityFunc while trying to build expression trees for the expressions already tokenized and parsed.\n"
		"Cannot continue.  Contact Andrew Leaver-Fay (aleaverfay@gmail.com)" );
	return nullptr;
}

ExpressionCOP
EntityFunc::function_expression(
	FunctionTokenCOP function,
	utility::vector1< ExpressionCOP > const & args
) const
{
	std::string const fname = function->name();
	if ( fname == "max" ) {
		return ExpressionCOP( ExpressionOP( new MaxExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "min" ) {
		return ExpressionCOP( ExpressionOP( new MinExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "exp" ) {
		if ( args.size() != 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "exp expression construction requested with more than one argument: " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new ExpExpression( args[ 1 ] ) ) );
	} else if ( fname == "ln" ) {
		if ( args.size() != 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "ln expression construction requested with more than one argument: " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new LnExpression( args[ 1 ] ) ) );
	} else if ( fname == "pow" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "pow expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new PowExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "ite" ) {
		if ( args.size() != 3 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "ite expression construction requested with nargs != 3. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new ITEExpression( args[ 1 ], args[ 2 ], args[ 3 ] ) ) );
	} else if ( fname == "abs" ) {
		if ( args.size() != 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "abs expression construction requested with nargs != 1. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new AbsoluteValueExpression( args[ 1 ] ) ) );
	} else if ( fname == "gt" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "gt expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new GT_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "lt" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "lt expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new LT_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "gte" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "gte expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new GTE_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "lte" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "lte expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new LTE_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "eq" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "eq expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new EqualsExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "and" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "and expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new AndExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "or" ) {
		if ( args.size() != 2 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "or expression construction requested with nargs != 2. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new OrExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "not" ) {
		if ( args.size() != 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "not expression construction requested with nargs != 1. Nargs= " + utility::to_string( args.size() )  );
		}
		return ExpressionCOP( ExpressionOP( new NotExpression( args[ 1 ] ) ) );
	}

	throw CREATE_EXCEPTION(utility::excn::Exception,  "Unable to find function with name " + fname + " in"
		" EntityFunc while trying to build expression trees for the expressions already tokenized and parsed.\n"
		"Cannot continue.  Contact Andrew Leaver-Fay (aleaverfay@gmail.com)" );
	return nullptr;

}

void
EntityFunc::initialize_scanner_and_function_names()
{
	scanner_ = numeric::expression_parser::ArithmeticScannerOP( new ArithmeticScanner( false ) );
	scanner_->add_function( "sqrt", 1 );
	scanner_->add_function( "max", 2 );
	scanner_->add_function( "min", 2 );
	scanner_->add_function( "exp", 1 );
	scanner_->add_function( "pow", 2 );
	scanner_->add_function( "ite", 3 );
	scanner_->add_function( "ln", 1 );
	scanner_->add_function( "abs", 1 );
	scanner_->add_function( "gt", 2 );
	scanner_->add_function( "lt", 2 );
	scanner_->add_function( "gte", 2 );
	scanner_->add_function( "lte", 2 );
	scanner_->add_function( "eq", 2 );
	scanner_->add_function( "and", 2 );
	scanner_->add_function( "or", 2 );
	scanner_->add_function( "not", 1 );

	function_names_.clear();
	function_names_.insert( "max" );
	function_names_.insert( "min" );
	function_names_.insert( "ln" );
	function_names_.insert( "pow" );
	function_names_.insert( "exp" );
	function_names_.insert( "sqrt" );
	function_names_.insert( "ite" );
	function_names_.insert( "abs" );
	function_names_.insert( "gt" );
	function_names_.insert( "lt" );
	function_names_.insert( "gte" );
	function_names_.insert( "lte" );
	function_names_.insert( "eq" );
	function_names_.insert( "and" );
	function_names_.insert( "or" );
	function_names_.insert( "not" );

	illegal_variable_names_.clear();
	illegal_variable_names_.insert( "in" );
	//illegal_variable_names_.insert( "" ); // what else don't I want to exclude?

	// Add variable names for the entity elements
	entity_aas_.resize( num_entity_elements_ );
	for ( Size ii = 1; ii <= num_entity_elements_; ++ii ) {
		std::string ii_name = "ee_" + utility::to_string( ii );
		scanner_->add_variable( ii_name );
		VariableExpressionOP ee_var_expression( new VariableExpression( ii_name, 0.0 ) );
		entity_aas_[ ii ] = ee_var_expression;
		variable_expression_map_[ ii_name ] = ee_var_expression;
		subexpression_name_map_[ ii_name ] = ee_var_expression;
		subexpression_name_dec_line_[ ii_name ] = 0;
	}
}

void
EntityFunc::process_AA_SET_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line
)
{
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read amino acid set name in the EntityFunc"
			" input file after reading AA_SET on line " + utility::to_string( line_number ) + "\n" + line );
	}
	std::string aa_setname;
	input_line >> aa_setname;
	if ( aa_setname.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read amino acid set name in the EntityFunc"
			" input file after reading AA_SET on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( aa_sets_name_map_.find( aa_setname ) != aa_sets_name_map_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Amino-acid-set name '" + aa_setname + "' appears multiple"
			" times in the EntityFunc file.\nFirst occurrance was found on line " +
			utility::to_string( aa_sets_dec_line_[ aa_setname ] ) +
			".  Second occurrance found while reading"
			" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	if ( subexpression_name_map_.find( aa_setname ) != subexpression_name_map_.end() ) {
		if ( subexpression_name_dec_line_[ aa_setname ] != 0 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal amino-acid-set name '" + aa_setname + "' which was previously"
				" declared in the EntityFunc file as a sub-expression.\nFirst occurrance was found on line " +
				utility::to_string( subexpression_name_dec_line_[ aa_setname ] ) +
				".  Second occurrance found while reading"
				" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		} else {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal amino-acid-set name '" + aa_setname + "' which is reserved"
				" as a name for an entity-element variable.\nError found while reading"
				" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
	}

	if ( function_names_.find( aa_setname ) != function_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal amino-acid-set name '" + aa_setname + "' which is reserved"
			" as a name for an funcion.\nError encountered while reading"
			" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	if ( illegal_variable_names_.find( aa_setname ) != illegal_variable_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal amino-acid-set name '" + aa_setname + "' which is a reserved"
			" word for this input file format.\nError encountered while reading"
			" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	/// whitespace is fine to separate the aaset from the equals sign.
	char equals_sign( ' ' );
	while ( equals_sign == ' ' || equals_sign == '\t' ) {
		if ( !input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an equals sign after reading amino-acid-set name'"
				+ aa_setname + "' but found an end-of-line.\nError encountered while reading"
				" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
		input_line.get( equals_sign );
	}
	if ( equals_sign != '=' ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an equals sign after reading amino-acid-set name'"
			+ aa_setname + "' but found '" + equals_sign + "'\nError encountered while reading"
			" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	/// whitespace is fine to separate the equals sign from the left curly bracket.
	char lcurly = ' ';
	while ( lcurly == ' ' || lcurly == '\t' ) {
		if ( !input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a left curly bracket ('{') after reading"
				" an equals sign, but found an end-of-line.\nError encountered while reading"
				" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
		input_line.get( lcurly );
	}
	if ( lcurly != '{' ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a left curly bracket ('{') after reading"
			" an equals sign, but found '"  + utility::to_string( lcurly ) + "'\nError encountered while reading"
			" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	utility::vector1< core::Real > aas_in_set;

	char next_aa = ' ';
	while ( next_aa != '}' ) {
		while ( next_aa == ' ' || next_aa == '\t' ) {
			if ( !input_line ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a right curly bracket ('}') or a 1-letter"
					"  amino acid code, but found an end-of-line.\nError encountered while reading"
					" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
			}
			input_line.get( next_aa );
		}
		if ( next_aa == '}' ) break;
		if ( core::chemical::oneletter_code_specifies_aa( toupper( next_aa )) ) {
			aas_in_set.push_back( core::chemical::aa_from_oneletter_code( toupper( next_aa )) );
		} else {
			/// maybe later, if we want to include ncaa's then they're parsing will happen right here.
			/// but for now, no go.
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a 1-letter"
				"  amino acid code or a right curly brace, but found '" + utility::to_string( next_aa ) + "'.\n"
				"Error encountered while reading AA_SET command\n" + line + "\n"
				"Line # " + utility::to_string( line_number )  );
		}
		next_aa = ' ';
		while ( next_aa == ' ' || next_aa == '\t' ) {
			if ( !input_line ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a right curly bracket ('}') or a 1-letter"
					"  amino acid code, but found an end-of-line.\nError encountered while reading"
					" AA_SET command\n" + line + "\nLine # " + utility::to_string( line_number )  );
			}
			input_line.get( next_aa );
		}
		if ( next_aa == '}' ) break;
		if ( next_aa == ',' ) {
			next_aa = ' ';
			continue;
		}
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read comma"
			"  or a right curly brace, but found '" + utility::to_string( next_aa ) + "'.\n"
			"Error encountered while reading AA_SET command\n" + line + "\n"
			"Line # " + utility::to_string( line_number )  );
	}

	/// The rest of the line is ignored.

	/// OK we made it this far -- now keep track of the aas that define this set.
	aa_sets_name_map_[ aa_setname ] = aas_in_set;
	aa_sets_dec_line_[ aa_setname ] = line_number;
}

void
EntityFunc::process_SET_CONDITION_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line
)
{
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a set-condition name in the EntityFunc"
			" input file after reading SET_CONDITION on line " + utility::to_string( line_number ) + "\n" + line );
	}
	std::string condition_name;
	input_line >> condition_name;
	if ( condition_name.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a set-condition name in the EntityFunc"
			" input file after reading SET_CONDITION on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( aa_sets_name_map_.find( condition_name ) != aa_sets_name_map_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Set-condition name '" + condition_name + "' appears multiple"
			" times in the EntityFunc file.\nFirst occurrance was found on line " +
			utility::to_string( aa_sets_dec_line_[ condition_name ] ) +
			".  Second occurrance found while reading"
			" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	if ( subexpression_name_map_.find( condition_name ) != subexpression_name_map_.end() ) {
		if ( subexpression_name_dec_line_[ condition_name ] != 0 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal set-condition name '" + condition_name + "' which was previously"
				" declared in the EntityFunc file as a sub-expression.\nFirst occurrance was found on line " +
				utility::to_string( subexpression_name_dec_line_[ condition_name ] ) +
				".  Second occurrance found while reading"
				" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		} else {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal set-condition name '" + condition_name + "' which is reserved"
				" as a name for an entity-element variable.\nError found while reading"
				" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
	}

	if ( function_names_.find( condition_name ) != function_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal set-condition name '" + condition_name + "' which is reserved"
			" as a name for an funcion.\nError encountered while reading"
			" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	if ( illegal_variable_names_.find( condition_name ) != illegal_variable_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal set-condition name '" + condition_name + "' which is a reserved"
			" word for this input file format.\nError encountered while reading"
			" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}
	/// now read an equals sign.
	/// whitespace is fine to separate the set-condition name from the equals sign.
	char equals_sign( ' ' );
	while ( equals_sign == ' ' || equals_sign == '\t' ) {
		if ( !input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an equals sign after reading the set-condition name'"
				+ condition_name + "' but found an end-of-line.\nError encountered while reading"
				" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
		input_line >> equals_sign;
	}
	if ( equals_sign != '=' ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an equals sign after reading amino-acid-set name'"
			+ condition_name + "' but found '" + utility::to_string( equals_sign ) + "'\nError encountered while reading"
			" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	/// now read the name of the entity-element that we're going to be examining.
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read the entity-element name in the EntityFunc"
			" input file after reading the set condition name, \"" + condition_name + ",\" but found an end-of-line.\n"
			"Error encountered while reading SET_CONDITION command on line " + utility::to_string( line_number ) + "\n" + line );
	}
	std::string eename;
	input_line >> eename;
	if ( subexpression_name_dec_line_.find( eename ) == subexpression_name_dec_line_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an entity-element name in the EntityFunc"
			" input file after reading an equals sign, while processing the set-condition named \"" + condition_name + ",\" but"
			" found \"" + eename +"\" which is not a valid entity element name.\n"
			"There are " + utility::to_string( num_entity_elements_ ) + " defined for this EntityFunc.\n"
			"Entity-element variables should be named ee_1 through ee_" + utility::to_string( num_entity_elements_ ) + ".\n"
			"Error encountered while reading SET_CONDITION command on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( subexpression_name_dec_line_[ eename ] != 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an entity-element name in the EntityFunc"
			" input file after reading an equals sign, while processing the set-condition named \"" + condition_name + ",\" but"
			" found the named subexpression \"" + eename +"\" instead.\n"
			"Error encountered while reading SET_CONDITION command on line " + utility::to_string( line_number ) + "\n" + line );
	}

	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read the string \"in\" in the EntityFunc"
			" input file after reading the set condition name, \"" + condition_name + ",\" but found an end-of-line.\n"
			"Error encountered while reading SET_CONDITION command on line " + utility::to_string( line_number ) + "\n" + line );
	}

	/// now read the "in" string
	std::string in_string;
	input_line >> in_string;
	if ( in_string != "in" ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read the string \"in\" in the EntityFunc"
			" input file after reading the set condition name, \"" + condition_name + ",\" but found \"" +
			in_string + "\" instead.\n"
			"Error encountered while reading SET_CONDITION command on line " + utility::to_string( line_number ) + "\n" + line );
	}

	char next_char = ' ';
	while ( input_line.peek() == ' ' || input_line.peek()  == '\t' ) {
		if ( !input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a left curly bracket ('}') or an"
				"  amino acid set name, but found an end-of-line.\nError encountered while reading"
				" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
		input_line.get( next_char );
	}

	// Now we either read an amino-acid-set name, or we read a '{'
	utility::vector1< core::Real > aas_in_set;

	if ( input_line.peek() == '{' ) {
		char next_aa;
		input_line.get( next_aa ); // wipe away the left-curly bracket
		next_aa = ' ';
		while ( next_aa != '}' ) {
			while ( next_aa == ' ' || next_aa == '\t' ) {
				if ( !input_line ) {
					throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a right curly bracket ('}') or a 1-letter"
						"  amino acid code, but found an end-of-line.\nError encountered while reading"
						" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
				}
				input_line.get( next_aa );
			}
			if ( next_aa == '}' ) break;
			if ( core::chemical::oneletter_code_specifies_aa( toupper( next_aa ) ) ) {
				aas_in_set.push_back( core::chemical::aa_from_oneletter_code( toupper( next_aa )) );
			} else {
				/// maybe later, if we want to include ncaa's then they're parsing will happen right here.
				/// but for now, no go.
				throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a 1-letter"
					"  amino acid code or a right curly brace, but found '" + utility::to_string( next_aa ) + "'.\n"
					"Error encountered while reading SET_CONDITION command\n" + line + "\n"
					"Line # " + utility::to_string( line_number )  );
			}
			next_aa = ' ';
			while ( next_aa == ' ' || next_aa == '\t' ) {
				if ( !input_line ) {
					throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a right curly bracket ('}') or a 1-letter"
						"  amino acid code, but found an end-of-line.\nError encountered while reading"
						" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
				}
				input_line.get( next_aa );
			}
			if ( next_aa == '}' ) { break; }
			if ( next_aa == ',' ) {
				next_aa = ' ';
				continue;
			}

			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read comma"
				"  or a right curly brace, but found '" + utility::to_string( next_aa ) + "'.\n"
				"Error encountered while reading SET_CONDITION command\n" + line + "\n"
				"Line # " + utility::to_string( line_number )  );
		}

	} else {
		std::string aa_setname;
		input_line >> aa_setname;
		if ( aa_sets_name_map_.find( aa_setname ) == aa_sets_name_map_.end() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Amino-acid-set name '" + aa_setname +
				"' has not previously been declared.\nError encountered while reading"
				" SET_CONDITION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
		aas_in_set = aa_sets_name_map_[ aa_setname ];
	}

	/// Alright -- we've made it this far.  Therefore, we're now ready to create a set condition expression.
	InSetExpressionOP inset_expression( new InSetExpression( subexpression_name_map_[ eename ] ) );
	inset_expression->value_set( aas_in_set );

	SurrogateVariableExpressionOP surrogate_expression( new SurrogateVariableExpression( condition_name ) );
	surrogate_expression->root_expression( inset_expression );
	expression_evaluation_order_.push_back( std::make_pair( inset_expression, surrogate_expression ) );

	variable_expression_map_[ condition_name ] = surrogate_expression;
	subexpression_name_map_[ condition_name ] = inset_expression;
	subexpression_name_dec_line_[ condition_name ] = line_number;
	scanner_->add_variable( condition_name );
}

void
EntityFunc::process_SUB_EXPRESSION_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line,
	std::map< std::string, ArithmeticASTExpressionOP > & expression_asts
)
{
	if ( ! input_line ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a sub-expression name in the EntityFunc"
			" input file after reading SUB_EXPRESSION on line " + utility::to_string( line_number ) + "\n" + line );
	}
	std::string subexpression_name;
	input_line >> subexpression_name;
	if ( subexpression_name.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read a sub-expression name in the EntityFunc"
			" input file after reading SUB_EXPRESSION on line " + utility::to_string( line_number ) + "\n" + line );
	}
	if ( aa_sets_name_map_.find( subexpression_name ) != aa_sets_name_map_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Sub-expression name '" + subexpression_name + "' appears multiple"
			" times in the EntityFunc file.\nFirst occurrance was found on line " +
			utility::to_string( aa_sets_dec_line_[ subexpression_name ] ) +
			".  Second occurrance found while reading"
			" SUB_EXPRESSION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	if ( subexpression_name_map_.find( subexpression_name ) != subexpression_name_map_.end() ) {
		if ( subexpression_name_dec_line_[ subexpression_name ] != 0 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal sub-expression name '" + subexpression_name + "' which was previously"
				" declared in the EntityFunc file as a sub-expression.\nFirst occurrance was found on line " +
				utility::to_string( subexpression_name_dec_line_[ subexpression_name ] ) +
				".  Second occurrance found while reading"
				" SUB_EXPRESSION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		} else {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal sub-expression name '" + subexpression_name + "' which is reserved"
				" as a name for an entity-element variable.\nError found while reading"
				" SUB_EXPRESSION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
	}

	if ( function_names_.find( subexpression_name ) != function_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal sub-expression name '" + subexpression_name + "' which is reserved"
			" as a name for an funcion.\nError encountered while reading"
			" SUB_EXPRESSION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	if ( illegal_variable_names_.find( subexpression_name ) != illegal_variable_names_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Illegal sub-expression name '" + subexpression_name + "' which is a reserved"
			" word for this input file format.\nError encountered while reading"
			" SUB_EXPRESSION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}
	/// now read an equals sign.
	/// whitespace is fine to separate the set-condition name from the equals sign.
	char equals_sign( ' ' );
	while ( equals_sign == ' ' || equals_sign == '\t' ) {
		if ( !input_line ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an equals sign after reading the sub-expression name '"
				+ subexpression_name + "' but found an end-of-line.\nError encountered while reading"
				" SUB_EXPRESSION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
		}
		input_line >> equals_sign;
	}
	if ( equals_sign != '=' ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Expected to read an equals sign after reading the sub-expression name '"
			+ subexpression_name + "' but found '" + utility::to_string( equals_sign ) + "'\nError encountered while reading"
			" SUB_EXPRESSION command\n" + line + "\nLine # " + utility::to_string( line_number )  );
	}

	std::string rest_of_line;
	std::getline( input_line, rest_of_line );

	TR << "On line " << line_number << ", attempting to tokenize scalar expression: " << rest_of_line << std::endl;
	TokenSetOP tokens = scanner_->scan( rest_of_line );
	TR << "On line " << line_number << ", attempting to parse scalar expression: " << rest_of_line << std::endl;
	ArithmeticASTExpressionOP expression_ast( new ArithmeticASTExpression );
	expression_ast->parse( *tokens );

	expression_asts[ subexpression_name ] = expression_ast;
	SurrogateVariableExpressionOP surrogate_expression( new SurrogateVariableExpression( subexpression_name ) );
	expression_evaluation_order_.push_back( std::make_pair( ExpressionCOP( nullptr ), surrogate_expression ) );

	variable_expression_map_[ subexpression_name ] = surrogate_expression;
	subexpression_name_map_[ subexpression_name ] = nullptr;
	subexpression_name_dec_line_[ subexpression_name ] = line_number;
	scanner_->add_variable( subexpression_name );

}

void
EntityFunc::process_SCORE_line(
	std::string const & line,
	Size line_number,
	std::istream & input_line,
	ArithmeticASTExpressionOP & score_expression_ast
)
{
	if ( score_expression_ast ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Encountered a second SCORE line"
			"while processing the EntityFunc file\n"
			+ line + "\nLine # " + utility::to_string( line_number )  );
	}

	std::string rest_of_line;
	std::getline( input_line, rest_of_line );

	TR << "On line " << line_number << ", attempting to tokenize score expression: " << rest_of_line << std::endl;
	TokenSetOP tokens = scanner_->scan( rest_of_line );
	TR << "On line " << line_number << ", attempting to parse score expression: " << rest_of_line << std::endl;
	score_expression_ast = ArithmeticASTExpressionOP( new ArithmeticASTExpression );
	score_expression_ast->parse( *tokens );

}

void
EntityFunc::turn_expression_ASTs_into_expressions(
	std::map< std::string, ArithmeticASTExpressionOP > const & expression_asts,
	ArithmeticASTExpressionOP score_expression_ast
)
{
	EntityFuncExpressionCreator expression_creator( *this );
	for ( Size ii = 1; ii <= expression_evaluation_order_.size(); ++ii ) {
		if ( expression_evaluation_order_[ ii ].first == nullptr ) {
			std::string ii_name = expression_evaluation_order_[ ii ].second->name();
			auto iter = expression_asts.find( ii_name );
			ExpressionCOP var_expr = expression_creator.create_expression_tree( * iter->second );
			subexpression_name_map_[ ii_name ] = var_expr;
			expression_evaluation_order_[ ii ].first = var_expr;
			expression_evaluation_order_[ ii ].second->root_expression( var_expr );
		}
	}

	score_expression_ = expression_creator.create_expression_tree( * score_expression_ast );
}

void EntityFunc::assign_entity_sequence_to_variables( Entity const & entity )
{
	using namespace protocols::multistate_design;

	runtime_assert( entity.traits().size() == num_entity_elements_ );
	for ( Size ii = 1; ii <= entity.traits().size(); ++ii ) {
		auto const & pt_ptr( dynamic_cast< PosType const & > ( * entity.traits()[ ii ] ));
		core::chemical::AA entity_aa( pt_ptr.type() );
		entity_aas_[ ii ]->set_value( entity_aa );
	}
}


}
}
