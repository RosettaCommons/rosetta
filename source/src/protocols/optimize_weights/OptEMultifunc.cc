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

// MPI headers
#ifdef USEMPI
#include <mpi.h>
#endif

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <protocols/optimize_weights/OptEMultifunc.hh>

// Project headers
#include <core/types.hh>

// AUTO-REMOVED #include <basic/options/util.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>

#include <basic/Tracer.hh>

#include <protocols/optimize_weights/OptEData.hh>
// AUTO-REMOVED #include <protocols/optimize_weights/NestedEnergyTermOptEData.hh>

/// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/format.hh> // F

// C++ headers
// AUTO-REMOVED #include <cmath>
#include <fstream>
#include <ostream>

// option key includes
#include <basic/options/keys/optE.OptionKeys.gen.hh>

//Auto Headers
#include <basic/options/option.hh>



using namespace core;
using namespace scoring;
using namespace optimization;

namespace protocols {
namespace optimize_weights {

using namespace ObjexxFCL::format;
using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.optimize_weights.OptEMultifunc" );

using namespace numeric::expression_parser;
//#undef NDEBUG


///
/// @begin OptEMultifunc::OptEMultifunc
///
OptEMultifunc::OptEMultifunc(
	OptEData & opte_data_in,
	EnergyMap const & fixed_terms_in,
	int num_free_in,
	ScoreTypes & score_list_in,
	ScoreTypes & fixed_score_list_in,
	Multivec const & component_weights
) :
	num_energy_dofs_( num_free_in ),
	num_ref_dofs_( chemical::num_canonical_aas ),
	num_total_dofs_( num_free_in + chemical::num_canonical_aas ),
	opte_data_( opte_data_in ),
	fixed_terms_( fixed_terms_in ),
	score_list_( score_list_in ),
	fixed_score_list_( fixed_score_list_in ),
	fix_reference_energies_( false ),
	starting_reference_energies_( chemical::num_canonical_aas, 0.0 ),
	component_weights_( component_weights ),
	mpi_rank_( 0 ),
	mpi_nprocs_( 0 ),
	distribute_over_mpi_( false )
{
#ifdef USEMPI
	MPI_Comm_rank( MPI_COMM_WORLD, & mpi_rank_);  /* get mpi rank */
	MPI_Comm_size( MPI_COMM_WORLD, & mpi_nprocs_);/* get number of processes */
	if ( basic::options::option[ basic::options::OptionKeys::optE::mpi_weight_minimization ] ) {
		distribute_over_mpi_ = true;
	}
	//TR << "OptEMultifunc created on node " << mpi_rank_ << std::endl;
#else
	(void) mpi_nprocs_;
#endif
}


///
/// @begin OptEMultifunc::OptEMultifunc
///
OptEMultifunc::OptEMultifunc(
	OptEData & opte_data_in,
	EnergyMap const & fixed_terms_in,
	int num_free_in,
	ScoreTypes const & score_list_in,
	ScoreTypes const & fixed_score_list_in,
	utility::vector1< Real > const & reference_energies_in,
	Multivec const & component_weights
) :
	num_energy_dofs_( num_free_in ),
	num_ref_dofs_( reference_energies_in.size() ),
	num_total_dofs_( num_free_in + reference_energies_in.size() ),
	opte_data_( opte_data_in ),
	fixed_terms_( fixed_terms_in ),
	score_list_( score_list_in ),
	fixed_score_list_( fixed_score_list_in ),
	fix_reference_energies_( false ),
	starting_reference_energies_( reference_energies_in ),
	component_weights_( component_weights ),
	mpi_rank_( 0 ),
	mpi_nprocs_( 0 ),
	distribute_over_mpi_( false )
{
#ifdef USEMPI
	MPI_Comm_rank( MPI_COMM_WORLD, & mpi_rank_);  /* get mpi rank */
	MPI_Comm_size( MPI_COMM_WORLD, & mpi_nprocs_);/* get number of processes */
	if ( basic::options::option[ basic::options::OptionKeys::optE::mpi_weight_minimization ] ) {
		distribute_over_mpi_ = true;
	}
	//TR << "OptEMultifunc created on node " << mpi_rank_ << std::endl;
#endif
}


///
/// @begin OptEMultifunc::operator()
///
/// @brief
/// The objective function for optE. Called in IterativeOptEDriver when optimizing the weights. Sums over all of the
/// PositionData objects in the OptEData object, and calls get_score() on each of them.  Each PositionData object implements
/// print_score() and get_score() methods that print/return how good a weight set is for optimizing the metric that PositionData
/// object represents.
///
Real
OptEMultifunc::operator()( Multivec const & vars ) const
{

	if ( distribute_over_mpi_ && mpi_rank_ == 0 ) {
		mpi_broadcast_eval_func( vars );
	}

	Multivec local_vars( num_total_dofs_ );
	if ( fix_reference_energies_ ) {
		for ( Size ii = 1; ii <= vars.size(); ++ii ) {
			local_vars[ ii ] = vars[ ii ];
		}
		for ( Size ii = 1, iipnfree = num_energy_dofs_ + 1; ii <= starting_reference_energies_.size(); ++ii, ++iipnfree ) {
			local_vars[ iipnfree ] = starting_reference_energies_[ ii ];
		}
	} else  {
		local_vars = vars;
	}

	// Apply this energy map to the rotamer energies to get the 'score'
	Real score( 0.0 );
	Multivec dummy( local_vars.size(), 0.0 );

	// Summing over positions
	for( OptEPositionDataOPs::const_iterator itr = opte_data_.position_data_begin(),
			e_itr = opte_data_.position_data_end() ; itr != e_itr ; ++itr  ) {
		Real const s = (*itr)->get_score( component_weights_, local_vars, dummy,
			num_energy_dofs_, num_ref_dofs_, num_total_dofs_,
			fixed_terms_, score_list_, fixed_score_list_ );
#ifndef WIN32
		bool bad = false;
		if(      std::isinf(s) ) { std::cerr << "Introduced INF score with " << OptEPositionDataFactory::optE_type_name( (*itr)->type() ) << " " << (*itr)->tag() << std::endl; bad=true; }
		else if( std::isnan(s) ) { std::cerr << "Introduced NAN score with " << OptEPositionDataFactory::optE_type_name( (*itr)->type() ) << " " << (*itr)->tag() << std::endl; bad=true; }
		else
#endif
			score += s;
#ifndef WIN32
		if ( bad ) {
			if ( basic::options::option[ basic::options::OptionKeys::optE::limit_bad_scores ].user() )					//NaN and inf errors can accumulate into the gigabytes if optE is left unattended
			{	
				static Size count = 0;																
				++count;
				if (count > 100000)
					utility::exit(__FILE__,__LINE__, "Counted over 100,000 inf/NaN scores. Admitting defeat now.");
			}
			std::cerr << "vars: " << std::endl;
			for ( Size ii = 1; ii <= local_vars.size(); ++ii ) {
				if ( ii != 1 ) std::cerr << ", ";
				std::cerr << ii << " " << local_vars[ ii ];
			}
			std::cerr << std::endl;
		}
#endif
	}

	if ( distribute_over_mpi_ && mpi_rank_ == 0 ) {
		score += mpi_receive_func();
	}

#ifndef NDEBUG
	if ( TR.visible() && mpi_rank_ == 0 ) {
		TR << "OptEMultifunc " << score << "\n";
		TR << "Vars: ";
		for ( Size ii = 1; ii <= local_vars.size(); ++ii ) {
			TR << " " << local_vars[ ii ];
		}
		TR << std::endl;
		TR << "dVars: ";
		for ( Size ii = 1; ii <= local_vars.size(); ++ii ) {
			TR << " " << dummy[ ii ];
		}
		TR << std::endl;
	}
#endif


	return score;
}

/// @brief OptE dfunc -- gets the partial derivatives of func for each dimension being minimized
void
OptEMultifunc::dfunc( Multivec const & vars, Multivec & dE_dvars ) const
{
	if ( distribute_over_mpi_ && mpi_rank_ == 0 ) {
		mpi_broadcast_eval_dfunc( vars );
	}

	Multivec local_vars( num_total_dofs_ );
	Multivec local_dE_dvars( num_total_dofs_, 0.0 );
	if ( fix_reference_energies_ ) {
		for ( Size ii = 1; ii <= vars.size(); ++ii ) {
			local_vars[ ii ] = vars[ ii ];
		}
		for ( Size ii = 1, iipnfree = num_energy_dofs_ + 1; ii <= starting_reference_energies_.size(); ++ii, ++iipnfree ) {
			local_vars[ iipnfree ] = starting_reference_energies_[ ii ];
		}
	} else  {
		local_vars = vars;
	}

	for ( Size ii(1); ii <= dE_dvars.size(); ++ii ) dE_dvars[ ii ] = 0.0;

	// over positions
	Real score( 0.0 );
	for( OptEPositionDataOPs::const_iterator itr = opte_data_.position_data_begin(),
			e_itr = opte_data_.position_data_end() ; itr != e_itr ; ++itr  ) {
		score += (*itr)->get_score( component_weights_, local_vars, local_dE_dvars,
			num_energy_dofs_, num_ref_dofs_, num_total_dofs_,
			fixed_terms_, score_list_, fixed_score_list_ );
		for( Size ii = 1 ; ii <= local_dE_dvars.size() ; ++ii ) {
#ifndef WIN32
			if(      std::isinf(local_dE_dvars[ ii ]) ) std::cerr << "Introduced INF deriv at " << ii << " with " << OptEPositionDataFactory::optE_type_name( (*itr)->type() ) << " " << (*itr)->tag() << std::endl;
			else if( std::isnan(local_dE_dvars[ ii ]) ) std::cerr << "Introduced NAN deriv at " << ii << " with " << OptEPositionDataFactory::optE_type_name( (*itr)->type() ) << " " << (*itr)->tag() << std::endl;
#endif
		}
	}

	/// Only store the vars that are actually being optimized; local vars will contain the reference energies
	for ( Size ii = 1; ii <= dE_dvars.size(); ++ii ) {
		dE_dvars[ ii ] = local_dE_dvars[ ii ];
	}

	if ( distribute_over_mpi_ && mpi_rank_ == 0 ) {
		mpi_receive_dfunc( dE_dvars );
	}


	if ( mpi_rank_ == 0 && TR.visible() ) { /// only the master node should output
		TR << "score: " << F(9,5,score) << std::endl; /// deceptive when scoring across multiple nodes
		TR << "dfuncs: ";
		for ( Size ii = 1 ; ii <= dE_dvars.size() ; ++ii ) TR << " " << F(9,3,dE_dvars[ ii ]);
		TR << std::endl;
		TR << "vars: ";
		for ( Size ii = 1; ii <= vars.size(); ++ii ) {
			TR << " " << F(9,3,vars[ ii ]);
		}
		TR << std::endl;
	}

}


///
/// @begin OptEMultifunc::get_dofs_from_energy_map()
///
/// @brief Extract variable weights from an Energy Map
///
Multivec
OptEMultifunc::get_dofs_from_energy_map( EnergyMap const & start_vals ) const
{
	Multivec dofs( fix_reference_energies_ ? num_energy_dofs_ : num_total_dofs_, 0.0 );

	Size dof_index( 1 );
	for( ScoreTypes::const_iterator itr = score_list_.begin(),
			end_itr = score_list_.end(); itr != end_itr; ++itr ) {
		dofs[ dof_index++ ] = start_vals[ *itr ];
	}
	if ( ! fix_reference_energies_ ) {
		for ( Size ii = 1; ii <= starting_reference_energies_.size(); ++ii ) {
			dofs[ dof_index++ ] = starting_reference_energies_[ ii ];
		}
	}
	return dofs;
}

///
/// @begin OptEMultifunc::get_energy_map_from_dofs()
///
/// @brief Expand free variables and combine with fixed to make an Energy Map.  Used by the IterativeOptEDriver at the
/// end of weight minimization to create an EnergyMap that uses the new weight set.  This EnergyMap then gets output
/// into a optE log file.
///
EnergyMap
OptEMultifunc::get_energy_map_from_dofs( Multivec const & dofs) const
{
	EnergyMap return_map( fixed_terms_ );

	// This covers the variable weights
	Size dof_index( 1 );
	for( ScoreTypes::const_iterator itr = score_list_.begin(),
			end_itr = score_list_.end() ;
			itr != end_itr ; ++itr ) {
		return_map[ *itr ] = dofs[ dof_index++ ];
	}

	// Now the fixed weights - but I think we should have them from the
	// copy constructor above, so this is probably a total waste....
	for( ScoreTypes::const_iterator itr = fixed_score_list_.begin(),
			end_itr = fixed_score_list_.end() ;
			itr != end_itr ; ++itr ) {
		return_map[ *itr ] = fixed_terms_[ *itr ];
	}

	return return_map;
}

utility::vector1< Real >
OptEMultifunc::get_reference_energies_from_dofs(
	Multivec const & dofs
) const
{
	if ( fix_reference_energies_ ) {
		return starting_reference_energies_;
	} else {
		utility::vector1< Real > refEs( num_ref_dofs_, 0 );
		for ( int ii = 1; ii <= num_ref_dofs_; ++ii ) {
			refEs[ ii ] = dofs[ ii + num_energy_dofs_ ];
		}
		return refEs;
	}
}

/// @brief Non-driver node wait for MPI vars to evaluate either the func or the dfunc.
void
OptEMultifunc::wait_for_remote_vars() const
{
#ifdef USEMPI
	Multivec vars, dE_dvars;
	while ( true ) {
		int message;
		/// Wait for a message from node 0.
		MPI_Bcast( &message, 1, MPI_INT, 0, MPI_COMM_WORLD );

		//TR << "OptEMultifunc::wait_for_remote_vars " << mpi_rank_ << " " << message << std::endl;


		if ( message == EVAL_FUNC ) {
			mpi_broadcast_receive_vars( vars );
			double my_func = operator() ( vars );
			MPI_Send( & my_func, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
		} else if ( message == EVAL_DFUNC ) {
			mpi_broadcast_receive_vars( vars );
			dE_dvars.resize( vars.size() );
			std::fill( dE_dvars.begin(), dE_dvars.end(), 0.0 );
			dfunc( vars, dE_dvars );
			double * dE_dvars_raw = new double[ vars.size() ];
			for ( Size ii = 1; ii <= vars.size(); ++ii ) { dE_dvars_raw[ ii - 1 ] = dE_dvars[ ii ]; }
			MPI_Send( dE_dvars_raw, vars.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
		} else if ( message == END_OF_MINIMIZATION ) {
			break;
		} else {
			std::cerr << "ERROR:  Unrecognized message from root node: " << message << std::endl;
			break;
		}
	}
#endif

}

/// @brief For driver node: inform the non-driver nodes that minimization is over.  Must
/// be called before object is destructed (Should not be called in the destructor, as
/// dstors should not throw exceptions, and MPI communication can absolutely result in exceptions).
void OptEMultifunc::declare_minimization_over() const
{
#ifdef USEMPI
	int message = END_OF_MINIMIZATION;
	MPI_Bcast( &message, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

}


/// @brief send out messages over MPI for remote nodes to evaluate their func given the input vars.
void OptEMultifunc::mpi_broadcast_eval_func(
#ifdef USEMPI
	Multivec const & vars
#else
	Multivec const &
#endif
) const
{
#ifdef USEMPI
	int message = EVAL_FUNC;
	MPI_Bcast( &message, 1, MPI_INT, 0, MPI_COMM_WORLD );

	mpi_broadcast_send_vars( vars );

	//TR << "OptEMultifunc::broadcast eval func" << mpi_rank_ << std::endl;

#endif
}


/// @brief send out messages over MPI for remote nodes to evaluate their dfunc given the input vars.
void OptEMultifunc::mpi_broadcast_eval_dfunc(
#ifdef USEMPI
	Multivec const & vars
#else
	Multivec const &
#endif
) const
{
#ifdef USEMPI
	int message = EVAL_DFUNC;
	MPI_Bcast( &message, 1, MPI_INT, 0, MPI_COMM_WORLD );

	mpi_broadcast_send_vars( vars );
	//TR << "OptEMultifunc::broadcast eval dfunc" << mpi_rank_ << std::endl;
#endif
}

void OptEMultifunc::mpi_broadcast_send_vars(
#ifdef USEMPI
	Multivec const & vars
#else
	Multivec const &
#endif
) const
{
#ifdef USEMPI
	int vars_size = vars.size();
	MPI_Bcast( & vars_size, 1, MPI_INT, 0, MPI_COMM_WORLD );

	double * raw_vars = new double[ vars_size ];
	for ( int ii = 1; ii <= vars_size; ++ii ) {
		raw_vars[ ii - 1 ] = vars[ ii ];
	}
	MPI_Bcast( raw_vars, vars_size, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	delete [] raw_vars;
#endif
}

void OptEMultifunc::mpi_broadcast_receive_vars(
#ifdef USEMPI
	Multivec & vars
#else
	Multivec  &
#endif
) const
{
#ifdef USEMPI
	int vars_size;
	MPI_Bcast( & vars_size, 1, MPI_INT, 0, MPI_COMM_WORLD );
	vars.resize( vars_size );

	double * raw_vars = new double[ vars_size ];
	MPI_Bcast( raw_vars, vars_size, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	for ( int ii = 1; ii <= vars_size; ++ii ) {
		vars[ ii ] = raw_vars[ ii - 1 ];
	}
	delete [] raw_vars;

#endif
}


/// @brief collect func values from remote nodes and return their sum.
Real
OptEMultifunc::mpi_receive_func() const
{
	Real total( 0.0 );
#ifdef USEMPI
	MPI_Status stat;

	for ( int ii = 1; ii < mpi_nprocs_; ++ii ) {
		double ii_func( 0 );
		MPI_Recv( & ii_func, 1, MPI_DOUBLE, ii, 1, MPI_COMM_WORLD, & stat );
		total += ii_func;
	}
#endif
	return total;
}

/// @brief collect dfunc valresultsues from remote nodes and increment the values in the dE_dvars input array.
void
OptEMultifunc::mpi_receive_dfunc(
#ifdef USEMPI
	Multivec & dE_dvars
#else
	Multivec  &
#endif
) const
{
#ifdef USEMPI
	MPI_Status stat;

	double * dE_dvars_raw = new double[ dE_dvars.size() ];
	for ( int ii = 1; ii < mpi_nprocs_; ++ii ) {
		MPI_Recv( dE_dvars_raw, dE_dvars.size(), MPI_DOUBLE, ii, 1, MPI_COMM_WORLD, & stat );
		for ( Size jj = 1; jj <= dE_dvars.size(); ++jj ) {
			dE_dvars[ jj ] += dE_dvars_raw[ jj - 1 ];
		}
	}
	delete [] dE_dvars_raw;
#endif
}

WrapperOptEMultifunc::WrapperOptEMultifunc() {}

void WrapperOptEMultifunc::init(
	ScoreTypes const & free_score_list,
	Size free_count,
	ScoreTypes const & fixed_score_list,
	EnergyMap  const & fixed_scores,
	OptEMultifuncOP optEfunc
) {
	//optE_dof_expressions_( optEfunc->fix_reference_energies() ? free_count : free_count + chemical::num_canonical_aas, 0 ),
	//active_variables_( optE_dof_expressions_.size() ),
	free_score_list_ = free_score_list;
	n_real_dofs_ = free_count;
	fixed_score_list_ = fixed_score_list;
	fixed_scores_ = fixed_scores;
	multifunc_ = optEfunc;
	//n_new_dofs_ = 0;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;

	// XCode doesn't like shorthand if statements in the initializer list
	optE_dof_expressions_.resize( optEfunc->fix_reference_energies() ? free_count : free_count + chemical::num_canonical_aas, 0 );
	active_variables_.resize( optE_dof_expressions_.size() );

	//std::cout << "WrapperOptEMultifunc ctor: free_count= " << free_count << " free_score_list.size()= " << free_score_list.size() << std::endl;

	if ( ! option[ optE::wrap_dof_optimization ].user() ) {
		utility_exit_with_message( "Error in WrapperOptEMultifunc constructor.  Cannot create WrapperOptEMultifunc if optE::wrap_dof_optimization is not on the command line");
	}

	ArithmeticScanner as; // comes built-in with min, max and sqrt.

	for ( Size ii = 1; ii <= free_score_list.size(); ++ii ) {
		std::string iiname = ScoreTypeManager::name_from_score_type( free_score_list[ ii ] );
		free_score_names_.insert( iiname );
		valid_variable_names_.insert( iiname );
		optEmultifunc_dof_order_[ iiname ] = ii;
		as.add_variable( iiname );
	}

	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		std::string iiname = ScoreTypeManager::name_from_score_type( fixed_score_list[ ii ] );
		valid_variable_names_.insert( iiname );
		as.add_variable( iiname );
	}

	/// This variable absolutely should not ever be held as an owning pointer
	/// through a member variable of the class, or you will have a cycle
	/// in the ownership graph and leak memory.
	WrappedOptEExpressionCreator expression_creator( get_self_weak_ptr() );

	std::ifstream wrapper_file( option[ optE::wrap_dof_optimization ]()().c_str() );
	bool finished_new_dof_header( false );


	while ( wrapper_file ) {
		std::string dof_dec;
		std::string dof_name;
		wrapper_file >> dof_dec;
		if ( dof_dec == "" ) {
			if ( !wrapper_file ) {
				break;
			} else {
				utility_exit_with_message( "Expected NEW_DOF or DEPENDENT_DOF from " + option[ optE::wrap_dof_optimization ]()() + " but got empty string" );
			}
		}
		wrapper_file  >> dof_name;
		std::cout << "READ: " << dof_dec << " " << dof_name << std::endl;
		if ( dof_dec == "NEW_DOF" ) {
			if ( finished_new_dof_header ) {
				utility_exit_with_message( "Encountered NEW_DOF declaration after a DEPENDENT_DOF delcaration.  All NEW_DOF declarations must be at the top of the file" );
			} else {
				valid_variable_names_.insert( dof_name );
				new_dof_names_.insert( dof_name );
				as.add_variable( dof_name );
			}
		} else if ( dof_dec == "DEPENDENT_DOF" ) {
			if ( ! finished_new_dof_header ) {
				finished_new_dof_header = true;
				/// prepare to parse dof dependencies.
			}

			if ( valid_variable_names_.find( dof_name ) == valid_variable_names_.end() ) {
				utility_exit_with_message( "Error in WrapperOptEMultifunc::WrapperOptEMultifunc()\nDid not find dof " + dof_name + " in valid_variable_names_ set; either it is not a free dof or is listed as a dependent dof twice" );
			}
			valid_variable_names_.erase( dof_name );
			std::string equals_sign;
			wrapper_file >> equals_sign;
			if ( equals_sign != "=" ) {
				utility_exit_with_message( "Expected an equals sign after reading 'DEPENDENT_DOF " + dof_name + "' but instead read " + equals_sign );
			}

			std::string expression;
			getline( wrapper_file, expression );
			std::cout << "READ EXPRESSION: " << expression << std::endl;
			TokenSetOP tokens = as.scan( expression );
			ArithmeticASTExpression ast_expression;
			ast_expression.parse( *tokens );

			active_variables_this_dependent_dof_.clear();
			numeric::expression_parser::ExpressionCOP derived_dof_expression = expression_creator.create_expression_tree( ast_expression );

			Size derived_dof_index = optEmultifunc_dof_order_[ dof_name ];
			optE_dof_expressions_[ derived_dof_index ] = derived_dof_expression;
			active_variables_[ derived_dof_index ] = active_variables_this_dependent_dof_;
			std::cout << "Created expression for " << dof_name << " index# " << derived_dof_index << std::endl;
		} else {
			utility_exit_with_message( "Expected either NEW_DOF or DEPENDENT_DOF, but got " + dof_dec + " from file " + option[ optE::wrap_dof_optimization ]()() );
		}

	}

	for ( Size ii = 1; ii <= optE_dof_expressions_.size(); ++ii ) {
		if ( optE_dof_expressions_[ ii ] == 0 ) {
			/// Need to create a variable expression for this dof so that it may be updated
			/// in each function evaluation
			/// Two cases: ii is a named dof, or ii is a reference energy.
			if ( ii <= free_score_list_.size() ) {
				std::string iiname = ScoreTypeManager::name_from_score_type( free_score_list_[ ii ] );
				OptEVariableExpressionOP varexp;
				if ( dof_variables_.find( iiname ) == dof_variables_.end() ) {
					varexp = OptEVariableExpressionOP( new OptEVariableExpression( iiname ) );
					dof_variables_[ iiname ] = varexp;
				} else {
					varexp = (dof_variables_.find( iiname ))->second;
				}
				active_variables_[ ii ].insert( iiname );
				optE_dof_expressions_[ ii ] = varexp;
			} else {
				std::string iiname = "ref" + utility::to_string( ii - free_score_list_.size() );
				//std::cout << "Adding reference energy dof: " << ii << " " << iiname << std::endl;
				// Assumption: the reference energy should not alrady be in the dof list
				assert( dof_variables_.find( iiname ) == dof_variables_.end() );
				OptEVariableExpressionOP varexp( new OptEVariableExpression( iiname ) );
				dof_variables_[ iiname ] = varexp;
				optE_dof_expressions_[ ii ] = varexp;
				active_variables_[ ii ].insert( iiname );
			}
		}
	}

	n_real_dofs_ = dof_variables_.size();
	Size count_real_dofs( 1 );
	for ( std::map< std::string, OptEVariableExpressionOP >::iterator
			iter = dof_variables_.begin(), iter_end = dof_variables_.end();
			iter != iter_end; ++iter ) {
		iter->second->set_id( count_real_dofs );
		++count_real_dofs;
	}
	real_dof_deriviative_expressions_.resize( n_real_dofs_ );

	for ( Size ii = 1; ii <= optE_dof_expressions_.size(); ++ii ) {
		numeric::expression_parser::ExpressionCOP iiexp = optE_dof_expressions_[ ii ];
		for ( std::set< std::string >::const_iterator
				variter = active_variables_[ ii ].begin(),
				variter_end = active_variables_[ ii ].end();
				variter != variter_end; ++variter ) {
			numeric::expression_parser::ExpressionCOP iiexp_dvar = iiexp->differentiate( *variter );
			if ( iiexp_dvar == 0 ) {
				utility_exit_with_message( "Error constructing parital derivative for '" +
					name_from_score_type( free_score_list[ ii ] ) +
					"' by variable '" + *variter + "'.  Null pointer returned." );
			}
			Size varindex = dof_variables_.find( *variter )->second->id();
			std::cout << "Adding dof derivative expression for " << *variter << " index#: " << varindex << " which appears in the expression for optEdof # " << ii << std::endl;
			real_dof_deriviative_expressions_[ varindex ].push_back( std::make_pair( ii, iiexp_dvar ) );
		}
	}
}

WrapperOptEMultifunc::~WrapperOptEMultifunc()
{
	//std::cout << "WrapperOptEMultifunc::~WrapperOptEMultifunc()" << std::endl;
}

// @brief OptE func
core::Real
WrapperOptEMultifunc::operator ()( Multivec const & vars ) const
{
	utility::vector1< Real > optEvars = derived_dofs( vars );
	Real score = (*multifunc_)( optEvars );

	TR << "WrapperOptEMultifunc func: " << F(7,2,score) << std::endl;
	TR << "Vars: ";
	for ( Size ii = 1; ii <= vars.size(); ++ii ) {
		TR << " " << vars[ ii ];
	}
	TR << std::endl;
	return score;
}

/// @brief OptE dfunc
void
WrapperOptEMultifunc::dfunc(
	Multivec const & vars,
	Multivec & dE_dvars
) const
{
	utility::vector1< Real > optEvars = derived_dofs( vars );
	Multivec dmultifunc_dvars( optEvars.size() );
	std::fill( dE_dvars.begin(), dE_dvars.end(), 0.0 );
	multifunc_->dfunc( optEvars, dmultifunc_dvars );

	for ( Size ii = 1; ii <= real_dof_deriviative_expressions_.size(); ++ii ) {
		for ( std::list< std::pair< Size, numeric::expression_parser::ExpressionCOP > >::const_iterator
				iter = real_dof_deriviative_expressions_[ ii ].begin(),
				iter_end = real_dof_deriviative_expressions_[ ii ].end();
				iter != iter_end; ++iter ) {
			dE_dvars[ ii ] += (dmultifunc_dvars[ iter->first ]) * ( (*(iter->second))());
			//std::cout << "Dof " << ii << " " << iter->first << ": " << (dmultifunc_dvars[ iter->first ])  << " * " << ( (*(iter->second))()) << std::endl;
		}
	}
	if ( TR.visible() ) {
		TR << "WrapperOptEMultifunc dfuncs:";
		for( Size ii = 1 ; ii <= dE_dvars.size() ; ++ii ) TR << " " << F(7,2,dE_dvars[ ii ]);
		TR << std::endl;
	}
}

Size
WrapperOptEMultifunc::n_real_dofs() const
{
	return n_real_dofs_;
}

utility::vector1< Real >
WrapperOptEMultifunc::derived_dofs( Multivec const & vars ) const
{
	//std::cout << "WrapperOptEMultifunc::derived_dofs()\n";
	utility::vector1< Real > optEvars( optE_dof_expressions_.size() );
	for ( std::map< std::string, OptEVariableExpressionOP >::const_iterator
			iter = dof_variables_.begin(), iter_end = dof_variables_.end();
			iter != iter_end; ++iter ) {
		iter->second->update_value_from_list( vars );
		//std::cout << "variable: " << iter->first << " " << (*iter->second)() << "\n";
	}
	for ( Size ii = 1; ii <= optE_dof_expressions_.size(); ++ii ) {
		optEvars[ ii ] = (*optE_dof_expressions_[ ii ])();
		//std::cout << "optEvar " << ii << " = " << optEvars[ ii ] << "\n";
	}
	//std::cout << std::endl;

	return optEvars;
}


void
WrapperOptEMultifunc::print_dofs(
	Multivec const & vars,
	std::ostream & ostr
) const
{
	ostr << "WrapperOptEMultifunc dofs:\n";
	for ( std::map< std::string, OptEVariableExpressionOP >::const_iterator
			iter = dof_variables_.begin(), iter_end = dof_variables_.end();
			iter != iter_end; ++iter ) {
		iter->second->update_value_from_list( vars );
		ostr << iter->first << " : " << (*(iter->second))() << "\n";
	}
}


numeric::expression_parser::VariableExpressionOP
WrapperOptEMultifunc::register_variable_expression( std::string varname )
{
	if ( valid_variable_names_.find( varname ) != valid_variable_names_.end() ) {
		/// This is a valid variable name.

		if ( dof_variables_.find( varname ) != dof_variables_.end() ) {
			return dof_variables_[ varname ];
		} else {
			OptEVariableExpressionOP varexp( new OptEVariableExpression( varname ) );
			if ( fixed_score_names_.find( varname ) == fixed_score_names_.end() ) {
				/// This is not a fixed dof, therefor it must be updated at each score evaluation.
				dof_variables_[ varname ] = varexp;
			} else {
				/// This is a fixed dof, so it does not need updating during scoring.  Set its value now,
				/// and refrain from adding this variable to the dof_variables_ map.
				ScoreType fixed_variable_scoretype = core::scoring::ScoreTypeManager::score_type_from_name( varname );
				varexp->set_value( fixed_scores_[ fixed_variable_scoretype ] );
			}
			active_variables_this_dependent_dof_.insert( varname );
			return varexp;
		}
	} else {
		/// We don't have a valid variable name.  Why not?
		if ( dependent_dof_names_.find( varname ) != dependent_dof_names_.end() ) {
			/// It's no longer a valid variable name because it's a former free-score name but has already
			/// been listed as a dependent dof.
			utility_exit_with_message( "Variable '" + varname +"' appearing on the right hand side of a DEPENDENT_DOF statement\nwas previously listed as a dependent variable" );
		} else {
			std::cerr << "Error: variable expression with name '" << varname << "' is not a valid variable name." << std::endl;
			std::cerr << "Free variables:" << std::endl;
			for ( std::set< std::string >::const_iterator
					iter = free_score_names_.begin(), iter_end = free_score_names_.end();
					iter != iter_end; ++iter ) {
				std::cerr << *iter << std::endl;
			}
			std::cerr << "Fixed variables:" << std::endl;
			for ( std::set< std::string >::const_iterator
					iter = fixed_score_names_.begin(), iter_end = fixed_score_names_.end();
					iter != iter_end; ++iter ) {
				std::cerr << *iter << std::endl;
			}
			std::cerr << "New variables:" << std::endl;
			for ( std::set< std::string >::const_iterator
					iter = new_dof_names_.begin(), iter_end = new_dof_names_.end();
					iter != iter_end; ++iter ) {
				std::cerr << *iter << std::endl;
			}
		}

		utility_exit_with_message("Could not register variable '" + varname + "' as it is neither a valid free, fixed nor new DOF" );
	}
	return 0;
}

void
WrapperOptEMultifunc::set_multifunc( OptEMultifuncOP multifunc )
{
	multifunc_ = multifunc;
}

OptEVariableExpression::OptEVariableExpression(std::string const & name ) :
	numeric::expression_parser::VariableExpression( name ),
	id_( 0 )
{}

OptEVariableExpression::OptEVariableExpression(std::string const & name, core::Real value ) :
	numeric::expression_parser::VariableExpression( name, value ),
	id_( 0 )
{}

void
OptEVariableExpression::set_id( core::Size id )
{
	id_ = id;
}

void OptEVariableExpression::update_value_from_list(
	utility::vector1< core::Real > const & value_vector
)
{
	set_value( value_vector[ id_ ] );
}

WrappedOptEExpressionCreator::WrappedOptEExpressionCreator()
{}

WrappedOptEExpressionCreator::WrappedOptEExpressionCreator(
	WrapperOptEMultifuncAP multifunc
)
	: multifunc_( multifunc )
{}

WrappedOptEExpressionCreator::~WrappedOptEExpressionCreator() {}

numeric::expression_parser::ExpressionCOP
WrappedOptEExpressionCreator::handle_variable_expression( ArithmeticASTValue const & node )
{
	WrapperOptEMultifuncOP multifunc( multifunc_ );
	return multifunc->register_variable_expression( node.variable_name() );
}

numeric::expression_parser::ExpressionCOP
WrappedOptEExpressionCreator::handle_function_expression(
	FunctionTokenCOP function,
	utility::vector1< numeric::expression_parser::ExpressionCOP > const & /*args*/
)
{
	utility_exit_with_message( "WrappedOptEExpressionCreator cannot process function " + function->name() );
	return 0;
}

void
WrappedOptEExpressionCreator::set_wrapping_optE_multifunc( WrapperOptEMultifuncAP multifunc )
{
	multifunc_ = multifunc;
}


} // namespace optimize_weights
} // namespace protocols



