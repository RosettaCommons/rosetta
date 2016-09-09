// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/util/deriv_funcs.hh
/// @brief  Classes for automating the testing the derivative evaluation
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_util_symmetric_deriv_funcs_HH
#define INCLUDED_util_symmetric_deriv_funcs_HH

// Test headers
#include <test/util/deriv_funcs.hh>

// Project headers
#include <core/scoring/symmetry/SymmetricEnergies.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/optimization/symmetry/sym_atom_tree_minimize.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

class SymmetricAtomDerivValidator
{
private:
	core::pose::PoseOP pose_;
	core::scoring::symmetry::SymmetricScoreFunctionOP sfxn_;
	core::kinematics::MoveMapOP move_map_;
	core::optimization::symmetry::SymMinimizerMapOP sym_min_map_;
	//bool auto_update_; -- no autoupdate in symmetry
	bool nonzero_deriv_only_; // in compute_pose_atom_derivs, only output atoms with nonzero derivative vectors
	std::list< core::Size > res_for_derivs_list_;

public:
	inline SymmetricAtomDerivValidator() : nonzero_deriv_only_( false ) {}
	inline SymmetricAtomDerivValidator(
		core::pose::Pose const & pose,
		core::scoring::symmetry::SymmetricScoreFunction const & sfxn,
		core::kinematics::MoveMap const & move_map
	) :
		pose_( new core::pose::Pose( pose )),
		sfxn_( core::scoring::symmetry::symmetrize_scorefunction( sfxn ) ),
		move_map_( new core::kinematics::MoveMap( move_map )),
		nonzero_deriv_only_( false )
	{}
	inline void set_score_function( core::scoring::symmetry::SymmetricScoreFunction const & sfxn ) { sfxn_ = core::scoring::symmetry::symmetrize_scorefunction(sfxn); }
	inline void set_pose( core::pose::Pose const & p ) { pose_ = core::pose::PoseOP( new core::pose::Pose( p ) ); }
	inline void set_movemap( core::kinematics::MoveMap const & mm ) { move_map_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap( mm ) ); }
	inline void set_nonzero_deriv_only( bool setting ) { nonzero_deriv_only_ = setting; }
	inline void add_res_for_deriv( core::Size resid ) { res_for_derivs_list_.push_back( resid ); }

	/// Evaluate the derivatives for all the atoms in the list, and verify that they are within 10^-12
	/// of the precomputed values stored in the AtomDerivList
	inline
	void validate_atom_deriv_list( AtomDerivList const & atderivs )
	{
		using namespace core;
		using namespace kinematics;
		using namespace scoring::symmetry;
		using namespace conformation::symmetry;
		using namespace optimization::symmetry;

		TS_ASSERT( pose_ );
		TS_ASSERT( sfxn_ );
		TS_ASSERT( move_map_ );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			std::cout << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return;
		}

		(*sfxn_)( *pose_ );

		MoveMap semisym_move_map;
		SymAtomTreeMinimizer::make_semisymmetric_movemap( *pose_, *move_map_, semisym_move_map );

		SymmetricConformation const & symm_conf ( dynamic_cast< SymmetricConformation const & > ( pose_->conformation()) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		core::optimization::symmetry::SymMinimizerMap sym_min_map( *pose_, semisym_move_map, symm_info );

		pose_->energies().set_use_nblist( *pose_, sym_min_map.domain_map(), false );

		sfxn_->setup_for_minimizing( *pose_, sym_min_map );
		sfxn_->setup_for_derivatives( *pose_ );

		eval_mingraph_derivs( sym_min_map );

		for ( std::list< AtomDeriv >::const_iterator
				iter = atderivs.begin(), iter_end = atderivs.end();
				iter != iter_end; ++iter ) {
			core::id::AtomID id( iter->atid() );
			core::scoring::DerivVectorPair const & atderivs( sym_min_map.atom_derivatives( id.rsd() )[ id.atomno() ] );
			core::Vector F1( atderivs.f1() ), F2( atderivs.f2() );
			sfxn_->eval_npd_atom_derivative( id, *pose_, sym_min_map.domain_map(), F1, F2 );
			bool f1bad( false ), f2bad( false );
			for ( int ii = 0; ii < 3; ++ii ) {
				TS_ASSERT_DELTA( F1[ ii ], iter->f1()[ ii ], 1e-12 );
				if ( std::abs( F1[ ii ] - iter->f1()[ ii ] ) > 1e-12 ) {
					f1bad = true;
				}
			}
			if ( f1bad ) {
				std::cout << "Derivative evaluation for atom " << id.atomno() << " on residue " << id.rsd() << " failed: " << std::endl;
				std::cout << "Gold F1: (" << iter->f1()[ 0 ] << " " << iter->f1()[ 1 ] << " " << iter->f1()[ 2 ] << " )" << std::endl;
				std::cout << "Curr F1: (" << F1[ 0 ] << " " << F1[ 1 ] << " " << F1[ 2 ] << " )" << std::endl;
			}

			for ( int ii = 0; ii < 3; ++ii ) {
				TS_ASSERT_DELTA( F2[ ii ], iter->f2()[ ii ], 1e-12 );
				if ( std::abs( F2[ ii ] - iter->f2()[ ii ] ) > 1e-12 ) {
					f2bad = true;
				}
			}
			if ( f2bad ) {
				std::cout << "Derivative evaluation for atom " << id.atomno() << " on residue " << id.rsd() << " failed: " << std::endl;
				std::cout << "Gold F2: (" << iter->f2()[ 0 ] << " " << iter->f2()[ 1 ] << " " << iter->f2()[ 2 ] << " )" << std::endl;
				std::cout << "Curr F2: (" << F2[ 0 ] << " " << F2[ 1 ] << " " << F2[ 2 ] << " )" << std::endl;
			}
		}
	}

	/// Copy and paste the ouput generated by this function into a .cxxtest.hh file so that
	/// you can invoke the above function (validate_atom_deriv_list) to test that the code
	/// continues to evaluate the proper derivatives given a score function and a move map.
	inline
	void compute_pose_atom_derivs()
	{
		using namespace core;
		using namespace kinematics;
		using namespace scoring::symmetry;
		using namespace conformation::symmetry;
		using namespace optimization::symmetry;


		TS_ASSERT( pose_ );
		TS_ASSERT( sfxn_ );
		TS_ASSERT( move_map_ );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			std::cout << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return;
		}

		(*sfxn_)( *pose_ );

		MoveMap semisym_move_map;
		SymAtomTreeMinimizer::make_semisymmetric_movemap( *pose_, *move_map_, semisym_move_map );

		SymmetricConformation const & symm_conf ( dynamic_cast<SymmetricConformation const &> ( pose_->conformation()) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		core::optimization::symmetry::SymMinimizerMap sym_min_map( *pose_, semisym_move_map, symm_info );

		pose_->energies().set_use_nblist( *pose_, sym_min_map.domain_map(), false );

		sfxn_->setup_for_minimizing( *pose_, sym_min_map );
		sfxn_->setup_for_derivatives( *pose_ );

		int precision_original = std::cout.precision();
		std::cout.precision( 16 ); // write out at high precision.

		std::cout << "using namespace core;" << std::endl;
		std::cout << "using namespace core::id;" << std::endl;
		std::cout << "AtomDerivList adl;" << std::endl;

		if ( res_for_derivs_list_.empty() ) {
			for ( core::Size ii = 1; ii <= pose_->size(); ++ii ) {
				res_for_derivs_list_.push_back( ii );
			}
		}
		eval_mingraph_derivs( sym_min_map );
		for ( std::list< core::Size >::const_iterator iter = res_for_derivs_list_.begin(),
				iter_end = res_for_derivs_list_.end(); iter != iter_end; ++iter ) {
			core::Size ii = *iter;
			for ( core::Size jj = 1; jj <= pose_->residue( ii ).natoms(); ++jj ) {
				core::id::AtomID id( jj, ii );
				core::Vector F1( 0 ), F2( 0 );
				sfxn_->eval_npd_atom_derivative( id, *pose_, sym_min_map.domain_map(), F1, F2 );
				if ( nonzero_deriv_only_ ) {
					if ( F1.length() == 0.0 && F2.length() == 0.0 ) continue;
				}
				std::cout << "adl.add( AtomDeriv( AtomID( " << jj << ", " << ii << "), ";
				for ( int ii = 0; ii < 3; ++ii ) {
					std::cout << F1[ ii ] << ",";
				}
				for ( int ii = 0; ii < 2; ++ii ) {
					std::cout << F2[ ii ] << ",";
				}
				std::cout << F2[2] << "));" << std::endl;
			}
		}
		// restore the precision before leaving this function
		std::cout.precision( precision_original );
	}

	void validate_start_func_matches_start_score() {

		std::pair< core::Real, core::Real > start_score_func = setup_for_minimizing();

		core::Real start_score = start_score_func.first;
		core::Real start_func = start_score_func.second;
		TS_ASSERT_DELTA( start_score, start_func, 1e-12 );
		if ( std::abs( start_score - start_func ) > 1e-12 ) {
			std::cout << "Failed to match start_score and start_func in AtomDerivValidator::validate_start_func_matches_start_score()" << std::endl;
			std::cout << "Start score: " << start_score << " Start func: " << start_func << std::endl;
		}

	}

	void validate_start_func_matches_start_score( core::Real start_score_gold, bool output_start_score = false ) {

		std::pair< core::Real, core::Real > start_score_func = setup_for_minimizing();

		core::Real start_score = start_score_func.first;
		TS_ASSERT_DELTA( start_score_gold, start_score, 1e-12 );
		if ( output_start_score ) {
			int precision_original = std::cout.precision();
			std::cout.precision( 16 ); // write out at high precision.
			std::cout << "START SCORE: " << start_score << std::endl;
			std::cout.precision( precision_original );
		}

		core::Real start_func = start_score_func.second;
		TS_ASSERT_DELTA( start_score, start_func, 1e-12 );
		if ( std::abs( start_score - start_func ) > 1e-12 ) {
			std::cout << "Failed to match start_score and start_func in AtomDerivValidator::validate_start_func_matches_start_score()" << std::endl;
			std::cout << "Start score: " << start_score << " Start func: " << start_func << std::endl;
		}

	}

	inline
	void
	eval_mingraph_derivs(
		core::optimization::symmetry::SymMinimizerMap & sym_min_map
	)
	{
		using namespace core;
		using namespace core::scoring;
		using namespace core::scoring::symmetry;

		SymmetricEnergies const & symm_energies( dynamic_cast< SymmetricEnergies const & > ( pose_->energies()) );

		assert( symm_energies.minimization_graph() );
		assert( symm_energies.derivative_graph() );

		MinimizationGraphCOP mingraph  = symm_energies.minimization_graph();
		MinimizationGraphCOP dmingraph = symm_energies.derivative_graph();


		for ( Size ii = 1; ii <= pose_->size(); ++ii ) {
			MinimizationNode const & minnode =  * mingraph->get_minimization_node( ii );
			/// 1. eval intra-residue derivatives
			eval_atom_derivatives_for_minnode( minnode, pose_->residue( ii ), *pose_, sfxn_->weights(), sym_min_map.atom_derivatives( ii ) );
		}

		/// 2a. eval inter-residue derivatives
		for ( utility::graph::Node::EdgeListConstIter
				edgeit = mingraph->const_edge_list_begin(), edgeit_end = mingraph->const_edge_list_end();
				edgeit != edgeit_end; ++edgeit ) {
			MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
			Size const rsd1ind = minedge.get_first_node_ind();
			Size const rsd2ind = minedge.get_second_node_ind();
			conformation::Residue const & rsd1( pose_->residue( rsd1ind ));
			conformation::Residue const & rsd2( pose_->residue( rsd2ind ));
			ResSingleMinimizationData const & r1_min_data( mingraph->get_minimization_node( rsd1ind )->res_min_data() );
			ResSingleMinimizationData const & r2_min_data( mingraph->get_minimization_node( rsd2ind )->res_min_data() );

			eval_atom_derivatives_for_minedge( minedge, rsd1, rsd2,
				r1_min_data, r2_min_data, *pose_, sfxn_->weights(),
				sym_min_map.atom_derivatives( rsd1ind ), sym_min_map.atom_derivatives( rsd2ind ));
		}
		/// 2b. eval inter-residue derivatives
		for ( utility::graph::Node::EdgeListConstIter
				edgeit = dmingraph->const_edge_list_begin(), edgeit_end = dmingraph->const_edge_list_end();
				edgeit != edgeit_end; ++edgeit ) {
			MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
			Size const rsd1ind = minedge.get_first_node_ind();
			Size const rsd2ind = minedge.get_second_node_ind();
			conformation::Residue const & rsd1( pose_->residue( rsd1ind ));
			conformation::Residue const & rsd2( pose_->residue( rsd2ind ));
			ResSingleMinimizationData const & r1_min_data( dmingraph->get_minimization_node( rsd1ind )->res_min_data() );
			ResSingleMinimizationData const & r2_min_data( dmingraph->get_minimization_node( rsd2ind )->res_min_data() );

			eval_atom_derivatives_for_minedge( minedge, rsd1, rsd2,
				r1_min_data, r2_min_data, *pose_, sfxn_->weights(),
				sym_min_map.atom_derivatives( rsd1ind ), sym_min_map.atom_derivatives( rsd2ind ));
		}
	}

	inline
	bool // return true for success, false for failure
	simple_deriv_check( bool start_score_func_check, core::Real tolerance ) {
		using namespace core::optimization;
		using namespace core::optimization::symmetry;

		if ( start_score_func_check ) {
			/// go ahead and make sure that the start score matches the start func
			validate_start_func_matches_start_score();
		} else {
			setup_for_minimizing();
		}
		sfxn_->setup_for_derivatives( *pose_ );

		// setup the function that we will pass to the simple deriv checker
		SymAtomTreeMultifunc f( *pose_, *sym_min_map_, *sfxn_, false, false );

		// starting position -- "dofs" = Degrees Of Freedom
		Multivec vars( sym_min_map_->nangles() );
		sym_min_map_->copy_dofs_from_pose( *pose_, vars );

		Multivec dE_dvars( vars );
		f.dfunc( vars, dE_dvars );

		SimpleDerivCheckResult result =  simple_numeric_deriv_check( f, vars, dE_dvars, false, false, 1 );

		core::optimization::symmetry::SymMinimizerMap::const_iterator dof_iterator( sym_min_map_->begin() );
		bool full_success = true;

		for ( core::Size ii = 1; ii <= sym_min_map_->nangles(); ++ii, ++dof_iterator ) {

			TS_ASSERT_DELTA( result.step_data( ii, 1 ).num_deriv(), result.step_data( ii, 1 ).ana_deriv(), tolerance );
			if ( false ) { /// re-enable to look at all derivatives
				std::cout << "dof  " << ii << " " << (*dof_iterator)->dof_id() << std::endl;
				std::cout << "    F1: " << (*dof_iterator)->F1().x() << " " <<
					(*dof_iterator)->F1().y() << " " <<
					(*dof_iterator)->F1().z() << std::endl;
				std::cout << "    F2: " << (*dof_iterator)->F2().x() << " " <<
					(*dof_iterator)->F2().y() << " " <<
					(*dof_iterator)->F2().z() << std::endl;
				for ( core::Size jj = 1; jj <= (*dof_iterator)->atoms().size(); ++jj ) {
					core::id::AtomID const & id( (*dof_iterator)->atoms()[ jj ] );
					std::cout << "    Atom: " << id.rsd() << " " << id.atomno() << " " <<
						pose_->residue( id.rsd() ).name() << " " <<
						pose_->residue( id.rsd() ).atom_name( id.atomno() ) << std::endl;
				}
				std::cout << "    Numeric deriv: " << result.step_data( ii, 1 ).num_deriv() <<
					" analytic deriv: " << result.step_data( ii, 1 ).ana_deriv() << std::endl;
			}
			if ( std::abs( result.step_data( ii, 1 ).num_deriv() - result.step_data( ii, 1 ).ana_deriv() ) > tolerance ) {
				full_success = false;
				core::id::DOF_ID dofid( (*dof_iterator)->dof_id() );
				core::Size precision_old( std::cout.precision() );
				std::cout.precision( 16 );
				std::cout << "Minmap dof " << ii << " incorrectly computed for DOF: " << dofid.rsd() << " " <<
					dofid.atomno() << " " << dofid.type();
				if ( (*dof_iterator)->torsion_id().valid() ) {
					std::cout << "( " << (*dof_iterator)->torsion_id() << " )";
				}
				std::cout << " num_deriv: " << result.step_data( ii, 1 ).num_deriv() << " ana_deriv: " << result.step_data( ii, 1 ).ana_deriv();
				std::cout << std::endl;
				std::cout.precision( precision_old );
				for ( core::Size jj = 1; jj <= (*dof_iterator)->atoms().size(); ++jj ) {
					core::id::AtomID const & id( (*dof_iterator)->atoms()[ jj ] );
					std::cout << "  Atom: " << id.rsd() << " " << id.atomno() << " " <<
						pose_->residue( id.rsd() ).name() << " " <<
						pose_->residue( id.rsd() ).atom_name( id.atomno() ) << std::endl;
				}
			}
		}
		return full_success;
	}

	/// @brief Setup the minimizer map and ready the pose for minimization.  Return the start score and the start func.
	inline
	std::pair< core::Real, core::Real >
	setup_for_minimizing() {
		using namespace core;
		using namespace kinematics;
		using namespace scoring::symmetry;
		using namespace conformation::symmetry;
		using namespace optimization::symmetry;

		TS_ASSERT( pose_ );
		TS_ASSERT( sfxn_ );
		TS_ASSERT( move_map_ );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			std::cerr << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return std::make_pair( -1234, -1234 );
		}

		core::Real start_score = (*sfxn_)(*pose_);

		MoveMap semisym_move_map;
		SymAtomTreeMinimizer::make_asymmetric_movemap( *pose_, *move_map_, semisym_move_map ); //fd: new minimizer

		SymmetricConformation const & symm_conf ( dynamic_cast<SymmetricConformation const &> ( pose_->conformation()) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		sym_min_map_ = core::optimization::symmetry::SymMinimizerMapOP(
			new core::optimization::symmetry::SymMinimizerMap( *pose_, semisym_move_map, symm_info, true ) ); //fd: new minimizer

		pose_->energies().set_use_nblist( *pose_, sym_min_map_->domain_map(), false );
		sfxn_->setup_for_minimizing( *pose_, *sym_min_map_ );
		/// END AtomTreeMinimizer setup block

		core::Real start_func = (*sfxn_)(*pose_);

		return std::make_pair( start_score, start_func );
	}
};

#endif
