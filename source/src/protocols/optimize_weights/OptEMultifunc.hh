// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/OptEMultifunc.hh
/// @brief  OptE mode multifunction class
/// @author Jim Havranek


#ifndef INCLUDED_protocols_optimize_weights_OptEMultifunc_hh
#define INCLUDED_protocols_optimize_weights_OptEMultifunc_hh

// Unit headers
#include <protocols/optimize_weights/OptEMultifunc.fwd.hh>

// Package headers
#include <numeric/expression_parser/Arithmetic.hh>

// Project headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/scoring/ScoreType.hh> // there is no .fwd.hh for this

/// C++ headers
#include <map>
#include <set>
#include <iosfwd>

#include <core/scoring/EnergyMap.hh>
#include <protocols/optimize_weights/OptEData.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace optimize_weights {

enum OptEMultifuncMPIMessages {
	EVAL_FUNC,
	EVAL_DFUNC,
	END_OF_MINIMIZATION
};

/// @brief OptE mode multifunction class
class OptEMultifunc : public core::optimization::Multifunc
{

public:
	typedef core::scoring::ScoreTypes ScoreTypes;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::optimization::Multivec Multivec;
	typedef core::Real Real;
	typedef core::Size Size;

public: // Creation

	// c-tor
	OptEMultifunc(
		OptEData & opte_data_in,
		EnergyMap const & fixed_terms_in,
		int num_free_in,
		ScoreTypes & score_list_in,
		ScoreTypes & fixed_score_list_in,
		Multivec const & component_weights
	);

	// c-tor
	OptEMultifunc(
		OptEData & opte_data_in,
		EnergyMap const & fixed_terms_in,
		int num_free_in,
		ScoreTypes const & score_list_in,
		ScoreTypes const & fixed_score_list_in,
		utility::vector1< Real > const & reference_energies_in,
		Multivec const & component_weights
	);

	/// @brief Destructor
	virtual
	~OptEMultifunc()
	{}


public: // Methods


	// @brief OptE func
	virtual
	Real
	operator ()(  Multivec const & vars ) const;

	/// @brief OptE dfunc
	virtual
	void
	dfunc(
		Multivec const & vars,
		Multivec & dE_dvars
	) const;

	/// @brief Does actual work for OptE minimization
	/* Real
	get_score_at_single_position(
		OptEPositionDataOP const this_pos,
		Multivec const & vars,
		Multivec & dvars
	) const; */

	/// @brief Extract variable weights from an Energy Map
	Multivec
	get_dofs_from_energy_map( EnergyMap const & start_vals ) const;

	/// @brief Expand free variables and combine with fixed to make an Energy Map
	EnergyMap
	get_energy_map_from_dofs(
		Multivec const & dofs
	) const;

	utility::vector1< Real >
	get_reference_energies_from_dofs(
		Multivec const & dofs
	) const;

	void set_starting_reference_energies( utility::vector1< Real > const & values )
	{
		starting_reference_energies_ = values;
	}

	/// @brief Non-driver node wait for MPI vars to evaluate either the func or the dfunc.
	void wait_for_remote_vars() const;

	/// @brief For driver node: inform the non-driver nodes that minimization is over.  Must
	/// be called before object is destructed (Should not be called in the destructor, as
	/// dstors should not throw exceptions, and MPI communication can absolutely result in exceptions).
	void declare_minimization_over() const;

	void fix_reference_energies( bool setting ) {
		fix_reference_energies_ = setting;
	}

	/// @brief Are the reference energies being optimized at all, or are they being held fixed?
	bool fix_reference_energies() const {
		return fix_reference_energies_;
	}

private:

	/// @brief send out messages over MPI for remote nodes to evaluate their func given the input vars.
	void mpi_broadcast_eval_func( Multivec const & vars ) const;
	/// @brief send out messages over MPI for remote nodes to evaluate their dfunc given the input vars.
	void mpi_broadcast_eval_dfunc( Multivec const & vars ) const;

	void mpi_broadcast_send_vars( Multivec const & vars ) const;

	void mpi_broadcast_receive_vars( Multivec & vars ) const;


	/// @brief collect func values from remote nodes and return their sum.
	Real mpi_receive_func() const;
	/// @brief collect dfunc valresultsues from remote nodes and increment the values in the dE_dvars input array.
	void mpi_receive_dfunc( Multivec & dE_vars ) const;

private: // data

	Size const num_energy_dofs_;
	int const num_ref_dofs_;
	int const num_total_dofs_;

	/// Rotamer energy components for all positions
	OptEData const & opte_data_;
	EnergyMap const & fixed_terms_;
	ScoreTypes const & score_list_;
	ScoreTypes const & fixed_score_list_;


	bool fix_reference_energies_;
	utility::vector1< Real > starting_reference_energies_;

	Multivec component_weights_;

	int mpi_rank_;
	int mpi_nprocs_;
	bool distribute_over_mpi_;

}; // OptEMultifunc


/// @brief OptE mode multifunction class
/*
class NoRefernceEnergyOptEMultifunc : public core::optimization::Multifunc
{
public:
	typedef core::scoring::ScoreTypes ScoreTypes;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::optimization::Multivec Multivec;
	typedef core::Real Real;
	typedef core::Size Size;

public: // Creation

	// c-tor
	OptEMultifunc(
		OptEData & opte_data_in,
		EnergyMap const & fixed_terms_in,
		int num_free_in,
		ScoreTypes & score_list_in,
		ScoreTypes & fixed_score_list_in,
		Multivec const & component_weights
	);

	// c-tor
	OptEMultifunc(
		OptEData & opte_data_in,
		EnergyMap const & fixed_terms_in,
		int num_free_in,
		ScoreTypes const & score_list_in,
		ScoreTypes const & fixed_score_list_in,
		utility::vector1< Real > const & reference_energies_in,
		Multivec const & component_weights
	);

	/// @brief Destructor
	virtual
	~OptEMultifunc()
	{}


public: // Methods


	// @brief OptE func
	virtual
	Real
	operator ()(  Multivec const & vars ) const;

	/// @brief OptE dfunc
	virtual
	void
	dfunc(
		Multivec const & vars,
		Multivec & dE_dvars
	) const;

	/// @brief Does actual work for OptE minimization
	Real
	get_score_at_single_position(
		OptEPositionDataOP const this_pos,
		Multivec const & vars,
		Multivec & dvars
	) const;

	/// @brief Extract variable weights from an Energy Map
	Multivec
	get_dofs_from_energy_map( EnergyMap const & start_vals ) const;

	/// @brief Expand free variables and combine with fixed to make an Energy Map
	EnergyMap
	get_energy_map_from_dofs(
		Multivec const & dofs
	) const;

	void set_starting_reference_energies( utility::vector1< Real > const & values )
	{
		starting_reference_energies_ = values;
	}

private: // data

	core::optimization::MultifuncOP multifunc_;
	utility::vector1< Real > reference_energies_;
	utility::vector1< Real > expanded_vars_;

};*/


/// DANGER DANGER DANGER
/// This class must never be allocated on the stack.  Instead, it should be allocated
/// (with "new") on the heap.  This class hands an owning-pointer to itself to another
/// class to create a call-back mechanism; this owning pointer will be invalid and
/// result in stack corruption if this class is allocated on the stack.
class WrapperOptEMultifunc : public core::optimization::Multifunc, public utility::pointer::enable_shared_from_this< WrapperOptEMultifunc >
{
public:
	typedef core::scoring::ScoreTypes ScoreTypes;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::optimization::Multivec Multivec;
	typedef core::Real Real;
	typedef core::Size Size;

public:
	WrapperOptEMultifunc();
	
	void init(
		ScoreTypes const & free_score_list,
		Size free_count, // the number of named dofs (score types) + the number of reference energies
		ScoreTypes const & fixed_score_list,
		EnergyMap  const & fixed_scores,
		OptEMultifuncOP optEfunc
	);

	~WrapperOptEMultifunc();

	/// self pointers
	inline WrapperOptEMultifuncCOP get_self_ptr() const { return shared_from_this(); }
	inline WrapperOptEMultifuncOP get_self_ptr() { return shared_from_this(); }
	inline WrapperOptEMultifuncCAP get_self_weak_ptr() const { return WrapperOptEMultifuncCAP( shared_from_this() ); }
	inline WrapperOptEMultifuncAP get_self_weak_ptr() { return WrapperOptEMultifuncAP( shared_from_this() ); }

	// @brief OptE func
	virtual
	Real
	operator ()( Multivec const & vars ) const;

	/// @brief OptE dfunc
	virtual
	void
	dfunc(
		Multivec const & vars,
		Multivec & dE_dvars
	) const;

	Size
	n_real_dofs() const;

	utility::vector1< Real >
	derived_dofs( Multivec const & vars ) const;

	void
	print_dofs( Multivec const & vars, std::ostream & ostr ) const;

	numeric::expression_parser::VariableExpressionOP
	register_variable_expression( std::string varname );

	void
	set_multifunc( OptEMultifuncOP multifunc );

private:
	ScoreTypes free_score_list_;
	ScoreTypes fixed_score_list_;
	EnergyMap  fixed_scores_;

	std::set< std::string > free_score_names_;
	std::set< std::string > fixed_score_names_;

	/// Dofs added by "NEW_DOF" lines in the wrapper input file
	std::set< std::string > new_dof_names_;

	/// Dofs marked as dependent by "DEPENDENT_DOF" lines in the wrapper input file.
	std::set< std::string > dependent_dof_names_;

	/// union of ( free_score_names_  - dependent_dof_names_ ) + fixed_score_names_ + new_dof_names
	std::set< std::string > valid_variable_names_;

	/// During expression creation; keep track of which variables are active.  This set
	/// Eventually ends up in the active_varibles_ array, and will be used to create
	/// derivative expressions.
	std::set< std::string > active_variables_this_dependent_dof_;

	/// Not every VariableExpression object registered by the ExpressionVisitor
	/// will end up in this map -- variables for fixed_score terms are not kept here, they
	/// are simply given a fixed value and never asked to update that value.
	/// Every vaiable in this list has its value updated at the beginning of scoring.
	std::map< std::string, OptEVariableExpressionOP > dof_variables_;

	// what order should dofs appear in for the OptEMultifunc?
	std::map< std::string, Size > optEmultifunc_dof_order_;

	/// non-depdendent dofs are handled by VariableExpressions that map
	/// from the real_dof id to the optE-percieved-dof id.
	utility::vector1< numeric::expression_parser::ExpressionCOP > optE_dof_expressions_;

	/// Which variables influence each optE-percieved-dofs?
	utility::vector1< std::set< std::string > > active_variables_;
	/// For each real dof, a list of each derivative expression and associated optE-percieved-dof id.
	utility::vector1< std::list< std::pair< Size, numeric::expression_parser::ExpressionCOP > > > real_dof_deriviative_expressions_;

	OptEMultifuncOP multifunc_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Size n_new_dofs_;
	Size n_real_dofs_;

};

class OptEVariableExpression : public numeric::expression_parser::VariableExpression
{
public:
	OptEVariableExpression( std::string const & name );
	OptEVariableExpression( std::string const & name, core::Real value );

	void set_id( core::Size id );
	void update_value_from_list( utility::vector1< core::Real > const & value_vector );
	core::Size id() const { return id_; }

private:
	core::Size id_;
};


class WrappedOptEExpressionCreator : public numeric::expression_parser::ExpressionCreator
{
public:
	WrappedOptEExpressionCreator();
	WrappedOptEExpressionCreator( WrapperOptEMultifuncAP multifunc );
	virtual ~WrappedOptEExpressionCreator();

	virtual
	numeric::expression_parser::ExpressionCOP
	handle_variable_expression( numeric::expression_parser::ArithmeticASTValue const & );

	virtual
	numeric::expression_parser::ExpressionCOP
	handle_function_expression(
		numeric::expression_parser::FunctionTokenCOP function,
		utility::vector1< numeric::expression_parser::ExpressionCOP > const & args
	);

	void
	set_wrapping_optE_multifunc( WrapperOptEMultifuncAP multifunc );

private:

	WrapperOptEMultifuncAP multifunc_;

};


} // namespace optimization
} // namespace core


#endif // INCLUDED_protocols_optimize_weights_OptEMultifunc_HH
