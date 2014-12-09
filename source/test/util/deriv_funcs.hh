// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/util/deriv_funcs.hh
/// @brief  Classes for automating the testing the derivative evaluation
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_util_deriv_funcs_HH
#define INCLUDED_util_deriv_funcs_HH

// Project headers
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/atom_tree_minimize.hh>
#include <core/optimization/AtomTreeMultifunc.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>
#include <core/types.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <list>

inline
core::kinematics::MoveMap
create_movemap_to_allow_all_torsions() {
	core::kinematics::MoveMap movemap;
	movemap.set_bb( true );
	movemap.set_chi( true );
	return movemap;
}

inline
core::kinematics::MoveMap
create_trpcage_movemap_to_allow_bb10_freedom() {
	core::kinematics::MoveMap movemap;
	movemap.set_bb( 10, true );
	return movemap;
}

class AtomDeriv
{
public:
	inline
	AtomDeriv(
		core::id::AtomID const & at,
		core::Real f1x,
		core::Real f1y,
		core::Real f1z,
		core::Real f2x,
		core::Real f2y,
		core::Real f2z
	) :
		atid_( at ),
		f1_( f1x, f1y, f1z ),
		f2_( f2x, f2y, f2z )
	{}

	inline core::id::AtomID const & atid() const { return atid_; }
	inline core::Vector     const & f1()   const { return f1_;   }
	inline core::Vector     const & f2()   const { return f2_;   }

private:
	core::id::AtomID atid_;
	core::Vector     f1_;
	core::Vector     f2_;
};

class AtomDerivList
{
public:
	inline AtomDerivList() {}
	inline std::list< AtomDeriv >::const_iterator begin() const { return atom_derivs_.begin(); }
	inline std::list< AtomDeriv >::const_iterator end()   const { return atom_derivs_.end();   }
	inline void add( AtomDeriv const & atderiv ) { atom_derivs_.push_back( atderiv ); }

private:
	std::list< AtomDeriv > atom_derivs_;
};

class AtomDerivValidator
{
private:
	core::pose::PoseOP pose_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::kinematics::MoveMapOP move_map_;
	core::optimization::MinimizerMapOP min_map_;
	bool auto_update_;
	bool nonzero_deriv_only_; // in compute_pose_atom_derivs, only output atoms with nonzero derivative vectors
	std::list< core::Size > res_for_derivs_list_;

public:
	inline AtomDerivValidator() : auto_update_( false ), nonzero_deriv_only_( false ) {}
	inline AtomDerivValidator(
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
	inline void set_pose( core::pose::Pose const & p ) { pose_ = core::pose::PoseOP( new core::pose::Pose( p ) ); }
	inline void set_movemap( core::kinematics::MoveMap const & mm ) { move_map_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap( mm ) ); }
	inline void set_nblist_auto_update( bool setting ) { auto_update_ = setting; }
	inline void set_nonzero_deriv_only( bool setting ) { nonzero_deriv_only_ = setting; }
	inline void add_res_for_deriv( core::Size resid ) { res_for_derivs_list_.push_back( resid ); }

	/// Evaluate the derivatives for all the atoms in the list, and verify that they are within 10^-12
	/// of the precomputed values stored in the AtomDerivList
	inline
	void validate_atom_deriv_list( AtomDerivList const & atderivs, core::Real tolerance = 1e-12 )
	{

		TS_ASSERT( pose_ != 0 );
		TS_ASSERT( sfxn_ != 0 );
		TS_ASSERT( move_map_ != 0 );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			std::cout << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return;
		}

		(*sfxn_)( *pose_ );

		core::optimization::MinimizerMap min_map;
		min_map.setup( *pose_, *move_map_ );

		pose_->energies().set_use_nblist( *pose_, min_map.domain_map(), auto_update_ );

		sfxn_->setup_for_minimizing( *pose_, min_map );
		sfxn_->setup_for_derivatives( *pose_ );

		eval_mingraph_derivs( min_map );

		for ( std::list< AtomDeriv >::const_iterator
				iter = atderivs.begin(), iter_end = atderivs.end();
				iter != iter_end; ++iter ) {
			core::id::AtomID id( iter->atid() );
			core::scoring::DerivVectorPair const & atderivs( min_map.atom_derivatives( id.rsd() )[ id.atomno() ] );
			core::Vector F1( atderivs.f1() ), F2( atderivs.f2() );
			sfxn_->eval_npd_atom_derivative( id, *pose_, min_map.domain_map(), F1, F2 );
			bool f1bad( false ), f2bad( false );
			for ( int ii = 0; ii < 3; ++ii ) {
				TS_ASSERT_DELTA( F1[ ii ], iter->f1()[ ii ], tolerance );
				if ( std::abs( F1[ ii ] - iter->f1()[ ii ] ) > tolerance ) {
					f1bad = true;
				}
			}
			if ( f1bad ) {
				std::cout << "Derivative evaluation for atom " << id.atomno() << " on residue " << id.rsd() << " failed: " << std::endl;
				std::cout << "Gold F1: (" << iter->f1()[ 0 ] << " " << iter->f1()[ 1 ] << " " << iter->f1()[ 2 ] << " )" << std::endl;
				std::cout << "Curr F1: (" << F1[ 0 ] << " " << F1[ 1 ] << " " << F1[ 2 ] << " )" << std::endl;
			}

			for ( int ii = 0; ii < 3; ++ii ) {
				TS_ASSERT_DELTA( F2[ ii ], iter->f2()[ ii ], tolerance );
				if ( std::abs( F2[ ii ] - iter->f2()[ ii ] ) > tolerance ) {
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

	/// Copy and paste the result of this function into a .cxxtest.hh file so that
	/// you can invoke the above function (validate_atom_deriv_list) to test that the code
	/// continues to evaluate the proper derivatives given a score function and a move map.
	inline
	void compute_pose_atom_derivs()
	{

		TS_ASSERT( pose_ != 0 );
		TS_ASSERT( sfxn_ != 0 );
		TS_ASSERT( move_map_ != 0 );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			std::cout << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return;
		}

		(*sfxn_)( *pose_ );

		core::optimization::MinimizerMap min_map;
		min_map.setup( *pose_, *move_map_ );

		pose_->energies().set_use_nblist( *pose_, min_map.domain_map(), auto_update_ );

		sfxn_->setup_for_minimizing( *pose_, min_map );
		sfxn_->setup_for_derivatives( *pose_ );

		int precision_original = std::cout.precision();
		std::cout.precision( 16 ); // write out at high precision.

		std::cout << "using namespace core;" << std::endl;
		std::cout << "using namespace core::id;" << std::endl;
		std::cout << "AtomDerivList adl;" << std::endl;

		if ( res_for_derivs_list_.empty() ) {
			for ( core::Size ii = 1; ii <= pose_->total_residue(); ++ii ) {
				res_for_derivs_list_.push_back( ii );
			}
		}
		eval_mingraph_derivs( min_map );
		for ( std::list< core::Size >::const_iterator iter = res_for_derivs_list_.begin(),
				iter_end = res_for_derivs_list_.end(); iter != iter_end; ++iter ) {
			core::Size ii = *iter;
			for ( core::Size jj = 1; jj <= pose_->residue( ii ).natoms(); ++jj ) {
				core::id::AtomID id( jj, ii );
				core::Vector F1( min_map.atom_derivatives( ii )[ jj ].f1() ), F2( min_map.atom_derivatives( ii )[ jj ].f2() );
				sfxn_->eval_npd_atom_derivative( id, *pose_, min_map.domain_map(), F1, F2 );
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

	void validate_start_func_matches_start_score( core::Real start_score_gold, bool output_start_score = false, core::Real tolerance = 1e-12 ) {

		std::pair< core::Real, core::Real > start_score_func = setup_for_minimizing();

		core::Real start_score = start_score_func.first;
		TS_ASSERT_DELTA( start_score_gold, start_score, tolerance );
		if ( output_start_score ) {
			int precision_original = std::cout.precision();
			std::cout.precision( 16 ); // write out at high precision.
			std::cout << "START SCORE: " << start_score << std::endl;
			std::cout.precision( precision_original );
		}

		core::Real start_func = start_score_func.second;
		TS_ASSERT_DELTA( start_score, start_func, tolerance );
		if ( std::abs( start_score - start_func ) > tolerance ) {
			std::cout << "Failed to match start_score and start_func in AtomDerivValidator::validate_start_func_matches_start_score()" << std::endl;
			std::cout << "Start score: " << start_score << " Start func: " << start_func << std::endl;
		}

	}

	inline
	void
	eval_mingraph_derivs(
		core::optimization::MinimizerMap & min_map
	)
	{
		using namespace core;
		using namespace core::scoring;

		assert( pose_->energies().minimization_graph() );
		MinimizationGraphCOP mingraph = pose_->energies().minimization_graph();


		for ( Size ii = 1; ii <= pose_->total_residue(); ++ii ) {
			MinimizationNode const & minnode =  * mingraph->get_minimization_node( ii );
			/// 1. eval intra-residue derivatives
			eval_atom_derivatives_for_minnode( minnode, pose_->residue( ii ), *pose_, sfxn_->weights(), min_map.atom_derivatives( ii ) );
		}

		/// 2. eval inter-residue derivatives
		for ( graph::Node::EdgeListConstIter
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
				min_map.atom_derivatives( rsd1ind ), min_map.atom_derivatives( rsd2ind ));
		}
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
		AtomTreeMultifunc f( *pose_, *min_map_, *sfxn_, false, false );

		// starting position -- "dofs" = Degrees Of Freedom
		Multivec vars( min_map_->nangles() );
		min_map_->copy_dofs_from_pose( *pose_, vars );

		Multivec dE_dvars( vars );
		f.dfunc( vars, dE_dvars );

		SimpleDerivCheckResult result =  simple_numeric_deriv_check( f, vars, dE_dvars, false, false, 1 );

		core::optimization::MinimizerMap::const_iterator dof_iterator( min_map_->begin() );
		for ( int ii = 1; ii <= min_map_->nangles(); ++ii, ++dof_iterator ) {

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
	}

	/// @brief Setup the minimizer map and ready the pose for minimization.  Return the start score and the start func.
	inline
	std::pair< core::Real, core::Real >
	setup_for_minimizing() {

		TS_ASSERT( pose_ != 0 );
		TS_ASSERT( sfxn_ != 0 );
		TS_ASSERT( move_map_ != 0 );
		if ( ! pose_ || ! sfxn_ || ! move_map_ ) {
			std::cerr << "ERROR: AtomDerivValidator incorrectly initialized" << std::endl;
			return std::make_pair( -1234, -1234 );
		}

		core::Real start_score = (*sfxn_)(*pose_);

		/// BEGIN AtomTreeMinimizer setup block
		min_map_ = core::optimization::MinimizerMapOP( new core::optimization::MinimizerMap );
		min_map_->setup( *pose_, *move_map_ );

		pose_->energies().set_use_nblist( *pose_, min_map_->domain_map(), auto_update_ );
		sfxn_->setup_for_minimizing( *pose_, *min_map_ );
		/// END AtomTreeMinimizer setup block

		core::Real start_func = (*sfxn_)(*pose_);

		return std::make_pair( start_score, start_func );
	}
};

#endif
