// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/util/cart_deriv_funcs.hh
/// @brief  Classes for automating the testing of derivative evaluation
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_util_cart_deriv_funcs_HH
#define INCLUDED_util_cart_deriv_funcs_HH

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethod.hh>

#include <core/optimization/cartesian_minimize.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh> // for make_asymmetric_movemap


// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <list>


static basic::Tracer TR_cart_deriv_funcs("util.cart_deriv_funcs");

class CartAtomDeriv
{
public:
	inline
	CartAtomDeriv(
		core::id::AtomID const & at,
		core::Real dE_dx,
		core::Real dE_dy,
		core::Real dE_dz
	) :
		atid_( at ),
		dE_dxyz_( dE_dx, dE_dy, dE_dz )
	{}

	inline core::id::AtomID const & atid()    const { return atid_; }
	inline core::Vector     const & dE_dxyz() const { return dE_dxyz_;   }

private:
	core::id::AtomID atid_;
	core::Vector     dE_dxyz_;
};

class CartAtomDerivList
{
public:
	inline CartAtomDerivList() {}
	inline std::list< CartAtomDeriv >::const_iterator begin() const { return atom_derivs_.begin(); }
	inline std::list< CartAtomDeriv >::const_iterator end()   const { return atom_derivs_.end();   }
	inline void add( CartAtomDeriv const & atderiv ) { atom_derivs_.push_back( atderiv ); }

private:
	std::list< CartAtomDeriv > atom_derivs_;
};

class CartAtomDerivValidator
{
public:
	typedef core::id::AtomID_Map< core::Vector > AtomDerivMap;
	typedef core::optimization::symmetry::SymAtomTreeMinimizer SymAtomTreeMinimizer;
	typedef core::optimization::Multivec Multivec;

private:
	core::pose::PoseOP pose_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::kinematics::MoveMapOP move_map_;
	core::optimization::CartesianMinimizerMapOP min_map_;
	bool auto_update_;
	bool nonzero_deriv_only_; // in compute_pose_atom_derivs, only output atoms with nonzero derivative vectors
	std::list< core::Size > res_for_derivs_list_;

public:
	inline CartAtomDerivValidator() : auto_update_( false ), nonzero_deriv_only_( false ) {}
	inline CartAtomDerivValidator(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::kinematics::MoveMap const & move_map
	) :
		pose_( new core::pose::Pose( pose )),
		sfxn_( sfxn.clone() ),
		move_map_( new core::kinematics::MoveMap( move_map )),
		auto_update_( false ),
		nonzero_deriv_only_( false )
	{}
	inline void set_score_function( core::scoring::ScoreFunction const & sfxn ) { sfxn_ = sfxn.clone(); }
	inline void set_pose( core::pose::Pose const & p ) { pose_ = utility::pointer::make_shared< core::pose::Pose >( p ); }
	inline void set_movemap( core::kinematics::MoveMap const & mm ) {
		move_map_ = utility::pointer::make_shared< core::kinematics::MoveMap >( mm );
	}
	inline void set_nblist_auto_update( bool setting ) { auto_update_ = setting; }
	inline void set_nonzero_deriv_only( bool setting ) { nonzero_deriv_only_ = setting; }
	inline void add_res_for_deriv( core::Size resid ) { res_for_derivs_list_.push_back( resid ); }

	/// @brief Evaluate the derivatives for all the atoms in the list, and verify that they
	/// are within tolerance of the precomputed values stored in the AtomDerivList
	inline
	void validate_atom_deriv_list(
		CartAtomDerivList const & atderivs,
		core::Real tolerance = 1e-12
	)
	{
		using namespace core;

		TS_ASSERT( pose_ != 0 );
		TS_ASSERT( sfxn_ != 0 );
		TS_ASSERT( move_map_ != 0 );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			TR_cart_deriv_funcs << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return;
		}

		(*sfxn_)( *pose_ );

		kinematics::MoveMap asym_move_map;
		if ( pose::symmetry::is_symmetric( *pose_ ) ) {
			SymAtomTreeMinimizer::make_asymmetric_movemap( *pose_, *move_map_, asym_move_map );
		} else {
			asym_move_map = *move_map_;
		}

		core::optimization::CartesianMinimizerMap min_map;
		min_map.setup( *pose_, asym_move_map );

		pose_->energies().set_use_nblist( *pose_, min_map.domain_map(), auto_update_ );

		sfxn_->setup_for_minimizing( *pose_, min_map );
		activate_dof_deriv_terms_for_cart_min( *pose_, *sfxn_, min_map );

		sfxn_->setup_for_derivatives( *pose_ );

		AtomDerivMap curr_atom_derivs = eval_mingraph_derivs( min_map );

		for ( auto const & gold_atom_deriv : atderivs ) {
			core::id::AtomID id( gold_atom_deriv.atid() );
			core::Vector gold_dE_dxyz = gold_atom_deriv.dE_dxyz();
			core::Vector curr_dE_dxyz = curr_atom_derivs[ id ];
			bool bad = false;
			for ( int ii = 0; ii < 3; ++ii ) {
				TS_ASSERT_DELTA( gold_dE_dxyz[ ii ], curr_dE_dxyz[ ii ], tolerance );
				if ( std::abs( gold_dE_dxyz[ ii ] - curr_dE_dxyz[ ii ] ) > tolerance ) {
					bad = true;
				}
			}
			if ( bad ) {
				TR_cart_deriv_funcs << "Derivative evaluation for atom " << id.atomno() << " on residue " << id.rsd() << " failed: " << std::endl;
				TR_cart_deriv_funcs << "Gold dE_dxyz: (" << gold_dE_dxyz << " )" << std::endl;
				TR_cart_deriv_funcs << "Curr dE_dxyz: (" << curr_dE_dxyz << " )" << std::endl;
			}
		}

		// turn off nblist
		pose_->energies().reset_nblist();

	}

	/// @brief Copy and paste the result of this function into a .cxxtest.hh file so that
	/// you can invoke the above function (validate_atom_deriv_list) to test that the code
	/// continues to evaluate the proper derivatives given a score function and a move map.
	inline
	void compute_pose_atom_derivs()
	{
		using namespace core;

		TS_ASSERT( pose_ != 0 );
		TS_ASSERT( sfxn_ != 0 );
		TS_ASSERT( move_map_ != 0 );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			TR_cart_deriv_funcs << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return;
		}

		(*sfxn_)( *pose_ );

		core::kinematics::MoveMap asym_move_map;
		if ( pose::symmetry::is_symmetric( *pose_ ) ) {
			SymAtomTreeMinimizer::make_asymmetric_movemap( *pose_, *move_map_, asym_move_map );
		} else {
			asym_move_map = *move_map_;
		}

		core::optimization::CartesianMinimizerMap min_map;
		min_map.setup( *pose_, asym_move_map );

		pose_->energies().set_use_nblist( *pose_, min_map.domain_map(), auto_update_ );

		sfxn_->setup_for_minimizing( *pose_, min_map );
		activate_dof_deriv_terms_for_cart_min( *pose_, *sfxn_, min_map );

		sfxn_->setup_for_derivatives( *pose_ );

		int precision_original = TR_cart_deriv_funcs.precision();
		TR_cart_deriv_funcs.precision( 16 ); // write out at high precision.

		TR_cart_deriv_funcs << "using namespace core;" << std::endl;
		TR_cart_deriv_funcs << "using namespace core::id;" << std::endl;
		TR_cart_deriv_funcs << "CartAtomDerivList adl;" << std::endl;

		if ( res_for_derivs_list_.empty() ) {
			for ( core::Size ii = 1; ii <= pose_->size(); ++ii ) {
				res_for_derivs_list_.push_back( ii );
			}
		}
		AtomDerivMap curr_atom_derivs = eval_mingraph_derivs( min_map );
		for ( Size const ii : res_for_derivs_list_ ) {
			for ( core::Size jj = 1; jj <= pose_->residue( ii ).natoms(); ++jj ) {
				core::id::AtomID id( jj, ii );
				core::Vector dE_dxyz( curr_atom_derivs[ id ] );
				if ( nonzero_deriv_only_ ) {
					if ( dE_dxyz.length() == 0.0 ) continue;
				}
				TR_cart_deriv_funcs << "adl.add( AtomDeriv( AtomID( " << jj << ", " << ii << ")";
				for ( int ii = 0; ii < 3; ++ii ) {
					TR_cart_deriv_funcs << ", " << dE_dxyz[ ii ];
				}
				TR_cart_deriv_funcs << "));" << std::endl;
			}
		}
		// restore the precision before leaving this function
		TR_cart_deriv_funcs.precision( precision_original );

		// turn off nblist
		pose_->energies().reset_nblist();

	}

	void validate_start_func_matches_start_score() {

		std::pair< core::Real, core::Real > start_score_func = setup_for_minimizing();

		core::Real start_score = start_score_func.first;
		core::Real start_func = start_score_func.second;
		TS_ASSERT_DELTA( start_score, start_func, 1e-12 );
		if ( std::abs( start_score - start_func ) > 1e-12 ) {
			TR_cart_deriv_funcs << "Failed to match start_score and start_func in AtomDerivValidator::validate_start_func_matches_start_score()" << std::endl;
			TR_cart_deriv_funcs << "Start score: " << start_score << " Start func: " << start_func << std::endl;
		}

	}

	void validate_start_func_matches_start_score(
		core::Real start_score_gold,
		bool output_start_score = false,
		core::Real tolerance = 1e-12
	) {

		std::pair< core::Real, core::Real > start_score_func = setup_for_minimizing();

		core::Real start_score = start_score_func.first;
		TS_ASSERT_DELTA( start_score_gold, start_score, tolerance );
		if ( output_start_score ) {
			int precision_original = TR_cart_deriv_funcs.precision();
			TR_cart_deriv_funcs.precision( 16 ); // write out at high precision.
			TR_cart_deriv_funcs << "START SCORE: " << start_score << std::endl;
			TR_cart_deriv_funcs.precision( precision_original );
		}

		core::Real start_func = start_score_func.second;
		TS_ASSERT_DELTA( start_score, start_func, tolerance );
		if ( std::abs( start_score - start_func ) > tolerance ) {
			TR_cart_deriv_funcs << "Failed to match start_score and start_func in AtomDerivValidator::validate_start_func_matches_start_score()" << std::endl;
			TR_cart_deriv_funcs << "Start score: " << start_score << " Start func: " << start_func << std::endl;
		}

		// turn off nblist
		pose_->energies().reset_nblist();

	}

	inline
	core::id::AtomID_Map< core::Vector >
	eval_mingraph_derivs(
		core::optimization::CartesianMinimizerMap & min_map
	)
	{
		using namespace core;
		using namespace core::scoring;

		min_map.zero_stored_derivs();
		Multivec vars( min_map_->ndofs() );
		min_map_->copy_dofs_from_pose( *pose_, vars );

		Multivec dE_dvars( vars );
		optimization::cartesian_dfunc( *pose_, min_map, *sfxn_, vars, dE_dvars );

		id::AtomID_Map< Vector > atom_derivs( Vector( 0.0 ) );
		pose::initialize_atomid_map( atom_derivs, *pose_ );

		for ( Size ii = 1; ii <= pose_->total_residue(); ++ii ) {
			conformation::Residue const & ii_res( pose_->residue( ii ) );
			for ( Size jj = 1; jj <= ii_res.natoms(); ++jj ) {
				id::AtomID id( jj, ii );
				if ( min_map.atom_is_moving( id ) ) {
					atom_derivs[ id ] = dE_dvars[ min_map.get_atom_index(id) ];
				}
			}
		}
		return atom_derivs;
	}

	inline
	void
	simple_deriv_check( bool start_score_func_check, core::Real tolerance ) {
		using namespace core::optimization;

		if ( start_score_func_check ) {
			/// go ahead and make sure that the start score matches the start func
			validate_start_func_matches_start_score();
		} else {
			setup_for_minimizing();
		}
		sfxn_->setup_for_derivatives( *pose_ );

		// setup the function that we will pass to the simple deriv checker
		CartesianMultifunc f( *pose_, *min_map_, *sfxn_, false, false );

		Multivec vars( min_map_->ndofs() );
		min_map_->copy_dofs_from_pose( *pose_, vars );

		Multivec dE_dvars( vars );
		f.dfunc( vars, dE_dvars );

		NumericalDerivCheckResultOP result =
			utility::pointer::make_shared< NumericalDerivCheckResult >();
		result->send_to_stdout( false );

		cart_numerical_derivative_check(
			*min_map_, f, vars, dE_dvars, result, false );

		NumDerivCheckData const & result1( result->deriv_check_result( 1 ) );

		for ( int ii = 1; ii <= min_map_->ndofs(); ++ii ) {

			TS_ASSERT_DELTA( result1.step_data( ii, 1 ).num_deriv(), result1.step_data( ii, 1 ).ana_deriv(), tolerance );

			if ( std::abs( result1.step_data( ii, 1 ).num_deriv() - result1.step_data( ii, 1 ).ana_deriv() ) > tolerance ) {
				core::Size precision_old( TR_cart_deriv_funcs.precision() );
				TR_cart_deriv_funcs.precision( 16 );
				TR_cart_deriv_funcs << "Minmap dof " << ii << " from atom ";
				TR_cart_deriv_funcs << min_map_->get_atom( (ii-1) / 3 + 1 ) << " incorrectly computed:";
				TR_cart_deriv_funcs << " num_deriv: " << result1.step_data( ii, 1 ).num_deriv() << " ana_deriv: " << result1.step_data( ii, 1 ).ana_deriv();
				TR_cart_deriv_funcs << std::endl;
				TR_cart_deriv_funcs.precision( precision_old );
			}
		}

		// turn off nblist
		pose_->energies().reset_nblist();

	}

	/// @brief Setup the minimizer map and ready the pose for minimization.  Return the start score and the start func.
	inline
	std::pair< core::Real, core::Real >
	setup_for_minimizing() {
		using namespace core;
		using namespace core::optimization;

		TS_ASSERT( pose_ != 0 );
		TS_ASSERT( sfxn_ != 0 );
		TS_ASSERT( move_map_ != 0 );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			std::cerr << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return std::make_pair( -1234, -1234 );
		}

		core::Real start_score = (*sfxn_)(*pose_);

		/// BEGIN CartesianMinimizer setup block
		core::kinematics::MoveMap asym_move_map;
		if ( pose::symmetry::is_symmetric( *pose_ ) ) {
			SymAtomTreeMinimizer::make_asymmetric_movemap( *pose_, *move_map_, asym_move_map );
		} else {
			asym_move_map = *move_map_;
		}

		min_map_ = utility::pointer::make_shared< CartesianMinimizerMap >();
		min_map_->setup( *pose_, asym_move_map );

		pose_->energies().set_use_nblist( *pose_, min_map_->domain_map(), auto_update_ );
		sfxn_->setup_for_minimizing( *pose_, *min_map_ );
		activate_dof_deriv_terms_for_cart_min( *pose_, *sfxn_, *min_map_ );

		/// END CartesianMinimizer setup block

		core::Real start_func = (*sfxn_)(*pose_);

		return std::make_pair( start_score, start_func );

	}
};

#endif
