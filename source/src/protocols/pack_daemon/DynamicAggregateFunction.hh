// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/DynamicAggregateFunction.hh
/// @brief  Declaration for a file-driven definition for the active states in a multistate design
///         and the set of mathematical expressions that together define the fitness function
///         for the sequence being designed.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_pack_daemon_DynamicAggregateFunction_hh
#define INCLUDED_protocols_pack_daemon_DynamicAggregateFunction_hh

// Unit headers
#include <protocols/pack_daemon/DynamicAggregateFunction.fwd.hh>

// Package headers
#include <protocols/pack_daemon/MultistateAggregateFunction.hh>

// Project headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <numeric/expression_parser/Arithmetic.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/FileContentsMap.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <set>

#include <protocols/pack_daemon/PackDaemon.fwd.hh>
#include <utility/vector0.hh>


namespace protocols {
namespace pack_daemon {

class VectorExpression : public numeric::expression_parser::Expression
{
public:
	typedef numeric::expression_parser::Expression    parent;
	typedef utility::vector1< core::Real >             values;
	typedef numeric::expression_parser::ExpressionCOP ExpressionCOP;

public:
	VectorExpression( std::string const & name );
	virtual ~VectorExpression();

	/// @brief DO NOT CALL THIS FUNCTION.  Vector expressions return
	/// vectors of values instead of a singular value.
	virtual
	core::Real
	operator() () const;

	/// @brief DO NOT CALL THIS FUNCTION.  Vector expressions cannot
	/// be differentiated.
	virtual
	numeric::expression_parser::ExpressionCOP
	differentiate( std::string const & varname ) const;

	virtual
	values
	vector_values() const = 0;

	virtual
	utility::vector1< std::list< std::string > >
	active_variables_vector() const = 0;


	/// @brief Returns the number of vector values that this Expression returns without
	/// computing those values.
	virtual
	core::Size
	size() const = 0;

	std::string const & name() const;

private:
	std::string name_;
};

class VariableVectorExpression : public VectorExpression
{
public:
	typedef VectorExpression parent;
	typedef utility::vector1< numeric::expression_parser::VariableExpressionCOP > VariableExpressions;

public:
	VariableVectorExpression( std::string const & name, VariableExpressions const & vars );
	~VariableVectorExpression();

	virtual
	values
	vector_values() const;

	/// @brief Returns the number of variable expressions this VectorExpression points at
	virtual
	core::Size
	size() const;

	virtual
	std::list< std::string >
	active_variables() const;

	virtual
	utility::vector1< std::list< std::string > >
	active_variables_vector() const;


private:
	VariableExpressions vars_;

};

class IterativeVectorExpression : public VectorExpression
{
public:
	typedef VectorExpression parent;
	typedef numeric::expression_parser::ArithmeticASTExpression ArithmeticASTExpression;
	typedef numeric::expression_parser::VariableExpressionOP     VariableExpressionOP;
	typedef numeric::expression_parser::VariableExpressionCOP    VariableExpressionCOP;
public:
	IterativeVectorExpression( std::string const & name );
	~IterativeVectorExpression();

	void initialize(
		std::map< std::string, VectorExpressionCOP > const & vector_varnames,
		ArithmeticASTExpression const & expresion_ast,
		VectorExpressionCreator & expression_creator // holds a reference to my owning DynamicAggregateFunction
	);

	virtual
	values
	vector_values() const;

	virtual
	core::Size
	size() const;

	numeric::expression_parser::VariableExpressionCOP
	local_variable( std::string const & varname ) const;

	virtual
	std::list< std::string >
	active_variables() const;

	virtual
	utility::vector1< std::list< std::string > >
	active_variables_vector() const;


private:
	utility::vector1< VectorExpressionCOP > input_vector_expressions_;
	utility::vector1< VariableExpressionOP > local_variables_; // op data members can have non-const operations performed on them
	std::map< std::string, VariableExpressionCOP > local_variable_map_;
	ExpressionCOP expression_;
};

class VectorFunction : public numeric::expression_parser::UnaryExpression
{
public:
	typedef numeric::expression_parser::UnaryExpression parent;
	typedef numeric::expression_parser::ExpressionCOP   ExpressionCOP;
	typedef utility::vector1< core::Size >              ArgIndices;

public:

	VectorFunction( VectorExpressionCOP ex );
	virtual ~VectorFunction();

protected:
	VectorExpressionCOP vec_ex() const;

private:
	VectorExpressionCOP vec_ex_;
};

class VectorFunction2 : public numeric::expression_parser::BinaryExpression
{
public:
	typedef numeric::expression_parser::BinaryExpression parent;
	typedef numeric::expression_parser::ExpressionCOP    ExpressionCOP;
	typedef utility::vector1< core::Size >               ArgIndices;

public:

	VectorFunction2( VectorExpressionCOP ex1, VectorExpressionCOP ex2 );
	virtual ~VectorFunction2();

protected:
	VectorExpressionCOP vec_ex1() const;
	VectorExpressionCOP vec_ex2() const;

private:
	VectorExpressionCOP vec_ex1_;
	VectorExpressionCOP vec_ex2_;
};

class Mean : public VectorFunction
{
public:
	typedef VectorFunction parent;

public:

	Mean( VectorExpressionCOP ex );
	virtual ~Mean();

	virtual
	core::Real
	operator() () const;

	virtual
	numeric::expression_parser::ExpressionCOP
	differentiate( std::string const & varname ) const;

	virtual
	std::list< std::string >
	active_variables() const;

};


class VMax : public VectorFunction
{
public:
	typedef VectorFunction parent;

public:

	VMax( VectorExpressionCOP ex );
	virtual ~VMax();

	virtual
	core::Real
	operator() () const;

	virtual
	numeric::expression_parser::ExpressionCOP
	differentiate( std::string const & varname ) const;

	virtual
	std::list< std::string >
	active_variables() const;

};

class VMin : public VectorFunction
{
public:
	typedef VectorFunction parent;

public:
	VMin( VectorExpressionCOP ex );

	virtual ~VMin();

	virtual
	core::Real
	operator() () const;

	virtual
	numeric::expression_parser::ExpressionCOP
	differentiate( std::string const & varname ) const;

	virtual
	std::list< std::string >
	active_variables() const;

};

/// @brief Take two vector expressions of equal length; returns the value from position i in
/// expression 2 where position i is the position with the largest value in expression 1
class VMaxBy : public VectorFunction2
{
public:
	typedef VectorFunction2 parent;

public:

	VMaxBy( VectorExpressionCOP ex1, VectorExpressionCOP ex2 );
	virtual ~VMaxBy();

	virtual
	core::Real
	operator() () const;

	virtual
	numeric::expression_parser::ExpressionCOP
	differentiate( std::string const & varname ) const;

	virtual
	std::list< std::string >
	active_variables() const;

};

/// @brief Take two vector expressions of equal length; returns the value from position i in
/// expression 2 where position i is the position with the smallest value in expression 1
class VMinBy : public VectorFunction2
{
public:
	typedef VectorFunction2 parent;

public:
	VMinBy( VectorExpressionCOP ex1, VectorExpressionCOP ex2 );

	virtual ~VMinBy();

	virtual
	core::Real
	operator() () const;

	virtual
	numeric::expression_parser::ExpressionCOP
	differentiate( std::string const & varname ) const;

	virtual
	std::list< std::string >
	active_variables() const;

};

class PowExpression : public numeric::expression_parser::BinaryExpression
{
public:
	typedef numeric::expression_parser::BinaryExpression parent;
	typedef numeric::expression_parser::ExpressionCOP    ExpressionCOP;

public:
	PowExpression( ExpressionCOP base, ExpressionCOP exponent );
	virtual ~PowExpression();

	virtual
	core::Real
	operator() () const;

	virtual
	ExpressionCOP
	differentiate( std::string const & varname ) const;

};

class ExpExpression : public numeric::expression_parser::UnaryExpression
{
public:
	typedef numeric::expression_parser::UnaryExpression parent;
	typedef numeric::expression_parser::ExpressionCOP   ExpressionCOP;

public:
	ExpExpression( ExpressionCOP ex );
	virtual ~ExpExpression();

	virtual
	core::Real
	operator() () const;

	virtual
	ExpressionCOP
	differentiate( std::string const & varname ) const;

};

class LnExpression : public numeric::expression_parser::UnaryExpression
{
public:
	typedef numeric::expression_parser::UnaryExpression parent;
	typedef numeric::expression_parser::ExpressionCOP   ExpressionCOP;

public:
	LnExpression( ExpressionCOP ex );
	virtual ~LnExpression();

	virtual
	core::Real
	operator() () const;

	virtual
	ExpressionCOP
	differentiate( std::string const & varname ) const;

};


/// @brief Returns "true" if the expression ex evaluates to one
/// of a set of indicated values.
class InSetExpression : public numeric::expression_parser::UnaryExpression
{
public:
	typedef numeric::expression_parser::UnaryExpression parent;
	typedef numeric::expression_parser::ExpressionCOP   ExpressionCOP;

public:
	InSetExpression( ExpressionCOP ex );
	virtual ~InSetExpression();

	void value_set( utility::vector1< core::Real > const & values );

	virtual
	core::Real
	operator() () const;

	virtual
	ExpressionCOP
	differentiate( std::string const & varname ) const;

private:
	utility::vector1< core::Real > value_set_;
};

/// @brief Stores the result of the surragate expression as if this
/// expression were a variable, but defers to the root expression for
/// questions of deriviatives and which variables are active
/// (this is not a real variable).
class SurrogateVariableExpression : public numeric::expression_parser::VariableExpression
{
public:
	typedef numeric::expression_parser::VariableExpression parent;
	typedef numeric::expression_parser::ExpressionCOP      ExpressionCOP;
public:
	SurrogateVariableExpression( std::string const & );
	SurrogateVariableExpression( std::string const & , core::Real value );

	virtual
	std::list< std::string >
	active_variables() const;

	void
	root_expression( ExpressionCOP setting );

	virtual
	ExpressionCOP
	differentiate( std::string const & varname ) const;

private:
	ExpressionCOP root_expression_;

};


class VectorExpressionCreator : public numeric::expression_parser::ExpressionCreator
{
public:
	typedef numeric::expression_parser::ExpressionCreator parent;
	typedef numeric::expression_parser::ExpressionCOP     ExpressionCOP;
public:
	VectorExpressionCreator( DynamicAggregateFunction const & owner );
	virtual ~VectorExpressionCreator();

	/// @brief Override the parent-class definition of this function to trap
	/// one specific kind of function call: max or min with a vector variable
	/// as one of their function argument
	//virtual
	//void
	//visit( numeric::expression_parser::ArithmeticASTFunction const & );

	virtual
	ExpressionCOP
	handle_variable_expression( numeric::expression_parser::ArithmeticASTValue const & );

	virtual
	ExpressionCOP
	handle_function_expression(
		numeric::expression_parser::FunctionTokenCOP function,
		utility::vector1< ExpressionCOP > const & args
	);

private:
	DynamicAggregateFunction const & owner_;

};

class StructureFileNames
{
public:
	std::string pdb_name_;
	std::string correspondence_file_name_;
	std::string resfile_name_;
};

class DynamicAggregateFunction : public MultistateAggregateFunction
{
public:
	typedef MultistateAggregateFunction                parent;

	typedef protocols::genetic_algorithm::Entity           Entity;

	typedef numeric::expression_parser::Expression    Expression;
	typedef numeric::expression_parser::ExpressionOP  ExpressionOP;
	typedef numeric::expression_parser::ExpressionCOP ExpressionCOP;

	typedef numeric::expression_parser::VariableExpression    VariableExpression;
	typedef numeric::expression_parser::VariableExpressionOP  VariableExpressionOP;
	typedef numeric::expression_parser::VariableExpressionCOP VariableExpressionCOP;

	typedef numeric::expression_parser::ArithmeticASTExpressionOP ArithmeticASTExpressionOP;

	typedef utility::vector1< core::Real >   ExpressionValues;
	typedef core::Size                       Size;

public:
	DynamicAggregateFunction();
	virtual ~DynamicAggregateFunction();

	void set_num_entity_elements( Size setting );
	void set_score_function( core::scoring::ScoreFunction const & sfxn );

	core::Size num_states() const;
	core::Size num_npd_properties() const;
	core::Size num_npd_properties_for_state( core::Size state_id ) const;

	virtual core::Real   evaluate( StateEnergies const &, StateEnergies const & npd_properties, Entity const &  );
	virtual StateIndices select_relevant_states( StateEnergies const & en, StateEnergies const & npd, Entity const & );


	ExpressionCOP
	variable_expression( numeric::expression_parser::ArithmeticASTValue const & ) const;

	ExpressionCOP
	function_expression(
		numeric::expression_parser::FunctionTokenCOP function,
		utility::vector1< ExpressionCOP > const & args
	) const;

	/// @brief Pair a file name with a string -- instead of opening a file
	/// when asked to read this file, the DAF will use the contents of the file
	/// as provided.  Useful for testing purposes.
	void add_file_contents( std::string const & fname, std::string const & contents );

	std::string state_name( Size state_index ) const;

	StructureFileNames const &
	file_inputs_for_job( int job_index ) const;

	/// @brief Initialize from an input stream without trying to communicate with
	/// remote nodes
	void read_all_variables_from_input_file( std::istream & input );

	std::list< std::pair< Size, std::string > >::const_iterator
	npd_variable_indices_for_state_begin( core::Size state_id ) const;

	std::list< std::pair< Size, std::string > >::const_iterator
	npd_variable_indices_for_state_end( core::Size state_id ) const;

private:

	void
	initialize_scanner();

	void
	process_STATE_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line
	);

	void
	process_STATE_VECTOR_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line,
		utility::vector1< std::pair< std::string, std::string > > & strucvec_filenames
	);

	void
	process_POSE_ENERGY_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line
	);

	void
	process_POSE_ENERGY_VECTOR_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line
	);


	void
	process_NPD_PROPERTY_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line
	);

	void
	process_VECTOR_VARIABLE_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line,
		std::map< std::string, std::list< std::string > > & vector_variables
	);


	void
	process_SCALAR_EXPRESSION_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line,
		std::map< std::string, ArithmeticASTExpressionOP > & scalar_expression_asts
	);

	void
	process_VECTOR_EXPRESSION_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line,
		std::map< std::string, std::pair< std::map< std::string, std::string >, ArithmeticASTExpressionOP > > & vector_expression_list
	);

	void
	process_ENTITY_FUNCTION_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line
	);

	void
	process_FITNESS_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line,
		ArithmeticASTExpressionOP & fitness_expression_ast
	);

	void
	save_scalar_variable( std::string const & varname, core::Size line_number );

	void
	save_vector_variable( std::string const & varname, core::Size line_number );

	void
	read_state_vector_file(
		std::string const & vec_varname,
		std::string const & fname,
		Size & n_vector_states
	);

	void
	create_state_variable_expressions(
		Size & count_state,
		Size & count_npd_index,
		Size & count_variable_index
	);

	void
	create_variable_vector_expressions(
		Size & count_state,
		Size & count_npd_index,
		Size & count_variable_index
	);

	void
	create_scalar_and_vector_expression_variable_expressions(
		std::map< std::string, ArithmeticASTExpressionOP > const & scalar_expression_asts,
		std::map< std::string, std::list< std::string > > const & vector_variables,
		Size & count_state
	);

	void
	turn_expression_ASTs_into_expressions(
		std::map< std::string, ArithmeticASTExpressionOP > const & scalar_expression_asts,
		std::map< std::string, std::pair< std::map< std::string, std::string >, ArithmeticASTExpressionOP > >  const & vector_expression_asts,
		ArithmeticASTExpressionOP fitness_expression_ast
	);

	utility::vector1< VectorExpressionCOP >
	verify_vector_arguments(
		std::string const & fname,
		utility::vector1< ExpressionCOP > const & args,
		Size expected_nargs
	) const;

	void
	verify_variable_name_or_throw(
		std::string const & vname,
		std::string const & command_name,
		std::string const & line,
		Size line_number
	);


	void
	assign_state_energies_to_variables_and_subexpressions(
		StateEnergies const & state_energies,
		StateEnergies const & npd_properties,
		Entity const & entity,
		bool verbose = false
	);

	/// @brief used to determine the number of requested NPD properties that will be
	/// calculated; used to size the variable_expressions_for_npd_properties_ array.
	Size count_num_npd_properties() const;

protected:

	void
	count_file_reads();

	std::string
	get_file_contents(
		std::string const & filename
	);

private:

	Size num_entity_elements_;
	core::scoring::ScoreFunctionOP sfxn_;

	ExpressionCOP fitness_exp_;
	utility::vector1< VariableExpressionOP > variable_expressions_for_states_;
	utility::vector1< std::list< std::pair< Size, std::string > > > npd_variable_indices_for_states_;
	utility::vector1< VariableExpressionOP > variable_expressions_for_npd_properties_;
	//utility::vector1< utility::vector1< VariableExpressionOP > > npd_properties_for_states_;
	utility::vector1< VariableExpressionOP > variable_expressions_;

	numeric::expression_parser::ArithmeticScannerOP scanner_;

	std::list< VectorFunctionOP > vfuncs_;
	std::list< std::pair< Size, std::string > > expression_evaluation_order_by_name_; // first == 1 for scalar, 2 for vector, second = name
	utility::vector1< std::pair< Size, ExpressionCOP > > scalar_expressions_;
	std::map< Size, SurrogateVariableExpressionOP > surrogate_expression_map_;
	std::map< std::string, VariableExpressionCOP > scalar_expression_map_;
	std::map< std::string, VectorExpressionCOP > vector_expression_map_;

	// Variables that correspond to protein states
	std::set< std::string > state_variable_names_;
	std::set< std::string > state_vector_variable_names_;
	std::map< std::string, ExpressionCOP > named_state_expression_map_;
	std::map< std::string, VariableVectorExpressionOP > state_vector_variables_;

	// first = property, second = name of variable for that property
	std::map< std::string, std::list< std::pair< std::string, std::string > > > npd_properties_for_state_variables_;

	std::map< std::string, ExpressionCOP > npd_property_expression_map_;
	std::map< std::string, VariableVectorExpressionOP > npd_property_vector_expression_map_;

	/// Function names may not be used as variable names.
	std::set< std::string > function_names_;

	/// This set represents strings that may not be used as variable names
	/// but which are not funciton names
	std::set< std::string > illegal_variable_names_;

	//// The following variables hold variable names and the lines on which they were declared.
	std::map< std::string, Size > variable_names_dec_line_;
	std::map< std::string, Size > scalar_variable_names_dec_line_;
	std::map< std::string, Size > vector_variable_names_dec_line_;

	std::map< std::string, std::pair< EntityFuncOP, VariableExpressionOP > > entity_funcs_;
	std::map< std::string, Size > entity_funcs_dec_line_;

	/// Keep track for each varible name (or scalar expression name) its index in the ?? vector
	std::map< std::string, Size > variable_name_2_variable_exp_index_;
	/// Keep track for each state-variable name its state's index
	std::map< std::string, Size > state_variable_name_2_state_index_;

	std::map< std::string, StructureFileNames > named_state_data_file_names_;
	std::map< std::string, utility::vector1< StructureFileNames > > state_vector_data_file_names_;

	utility::vector1< StructureFileNames > files_for_state_;
	std::map< std::string, utility::vector1< Size > > state_indices_for_state_vector_;

	IterativeVectorExpressionOP focused_iterative_vector_expression_;

	utility::io::FileContentsMapOP file_contents_;

};

class DynamicAggregateFunctionDriver : public DynamicAggregateFunction
{
public:
	void initialize_from_input_file( DaemonSetOP daemon_set, std::istream & input );

private:

	/// @brief After the input stream has been read, assign the states to particular nodes and to the
	/// the pack daemon on this node.
	void initialize_pack_daemons( DaemonSetOP daemon_set );

	void distribute_jobs_to_remote_daemons( DaemonSetOP daemon_set );

	void assign_jobs_to_local_daemon_set( std::list< int > const & job_indices, DaemonSetOP daemon_set );

	void
	assign_jobs_to_remote_daemon_sets(
		int proc_id,
		std::list< int > const & job_indices
	);

	bool verify_remote_daemon_set_initialization_successful( int proc_id ) const;
	void send_success_message_to_remote_daemon_set( int proc_id ) const;
	void send_error_message_to_remote_daemon_sets() const;

	void initialize_daemon_with_all_states(
		DaemonSetOP daemon_set
	);


};


class EntityFuncExpressionCreator : public numeric::expression_parser::ExpressionCreator
{
public:
	typedef numeric::expression_parser::ExpressionCreator parent;
	typedef numeric::expression_parser::ExpressionCOP     ExpressionCOP;
public:
	EntityFuncExpressionCreator( EntityFunc const & owner );
	virtual ~EntityFuncExpressionCreator();

	virtual
	ExpressionCOP
	handle_variable_expression( numeric::expression_parser::ArithmeticASTValue const & );

	virtual
	ExpressionCOP
	handle_function_expression(
		numeric::expression_parser::FunctionTokenCOP function,
		utility::vector1< ExpressionCOP > const & args
	);

private:
	EntityFunc const & owner_;

};


class EntityFunc : public utility::pointer::ReferenceCount
{
public:
	typedef utility::pointer::ReferenceCount parent;

	typedef protocols::genetic_algorithm::Entity           Entity;

	typedef numeric::expression_parser::Expression    Expression;
	typedef numeric::expression_parser::ExpressionOP  ExpressionOP;
	typedef numeric::expression_parser::ExpressionCOP ExpressionCOP;

	typedef numeric::expression_parser::VariableExpression    VariableExpression;
	typedef numeric::expression_parser::VariableExpressionOP  VariableExpressionOP;
	typedef numeric::expression_parser::VariableExpressionCOP VariableExpressionCOP;

	typedef numeric::expression_parser::ArithmeticASTExpressionOP ArithmeticASTExpressionOP;

	typedef core::Size Size;

public:
	EntityFunc();
	virtual ~EntityFunc();

	void set_num_entity_elements( Size num_ee );

	void initialize_from_input_file( std::istream & input );

	core::Real
	evaluate( Entity const & entity, bool verbose = false );

	ExpressionCOP
	variable_expression( numeric::expression_parser::ArithmeticASTValue const & ) const;

	ExpressionCOP
	function_expression(
		numeric::expression_parser::FunctionTokenCOP function,
		utility::vector1< ExpressionCOP > const & args
	) const;

private:
	void
	initialize_scanner_and_function_names();

	void
	process_AA_SET_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line
	);

	void
	process_SET_CONDITION_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line
	);

	void
	process_SUB_EXPRESSION_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line,
		std::map< std::string, ArithmeticASTExpressionOP > & expression_asts
	);

	void
	process_SCORE_line(
		std::string const & line,
		Size line_number,
		std::istream & input_line,
		ArithmeticASTExpressionOP & score_expression_ast
	);

	void
	turn_expression_ASTs_into_expressions(
		std::map< std::string, ArithmeticASTExpressionOP > const & expression_asts,
		ArithmeticASTExpressionOP score_expression_ast
	);

	void assign_entity_sequence_to_variables( Entity const & entity );

private:
	Size num_entity_elements_;

	utility::vector1< VariableExpressionOP > entity_aas_;

	std::map< std::string, utility::vector1< core::Real > > aa_sets_name_map_;
	std::map< std::string, Size > aa_sets_dec_line_;

	std::map< std::string, ExpressionCOP > subexpression_name_map_;
	std::map< std::string, Size > subexpression_name_dec_line_;

	std::map< std::string, VariableExpressionOP > variable_expression_map_;

	/// Evaluate the subexpressions and then store them in the surrogate variable expressions
	/// in order.
	utility::vector1< std::pair< ExpressionCOP, SurrogateVariableExpressionOP > > expression_evaluation_order_;

	ExpressionCOP score_expression_;

	numeric::expression_parser::ArithmeticScannerOP scanner_;
	/// Function names may not be used as variable names.
	std::set< std::string > function_names_;
	/// This set represents strings that may not be used as variable names
	/// but which are not funciton names
	std::set< std::string > illegal_variable_names_;


};

}
}

#endif
