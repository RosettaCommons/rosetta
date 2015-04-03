// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/symmetry/SymAtomTreeMinimizer.cc
/// @brief  High-level atom tree minimizer class for symmetrical minimization
/// @author Ingemar Andre

// Unit headers
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

// Symmetry headers
#include <core/conformation/symmetry/util.hh>


// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>

#include <core/optimization/NumericalDerivCheckResult.fwd.hh>

// Project headers
#include <core/id/DOF_ID.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

//Auto Headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


using namespace ObjexxFCL::format;

namespace core {
namespace optimization {
namespace symmetry {

///////////////////////////////////////////////////////////////////////////////
Real
SymAtomTreeMinimizer::run(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map,
	scoring::ScoreFunction const & scorefxn,
	MinimizerOptions const & options
)
{
	check_setup( pose, move_map, scorefxn, options );

	using namespace core::conformation::symmetry;

	//typedef SymmetryInfo::DOF_IDs DOF_IDs;

	bool const use_nblist( options.use_nblist() );

	// it's important that the structure be scored prior to nblist setup
	Real const start_score( scorefxn( pose ) );

	// optionally try out a new approach to symmetric minimization
	bool const new_sym_min( basic::options::option[ basic::options::OptionKeys::optimization::new_sym_min ] );

	kinematics::MoveMap semisym_move_map;
	// I worry about this routine (make_assymetric_movemap) since it does not see dofs that were set individually (FIX!)
	if ( new_sym_min ) make_assymetric_movemap( pose, move_map, semisym_move_map ); // this is not quite right...
	else make_semisymmetric_movemap( pose, move_map, semisym_move_map );

	SymmetricConformation const & symm_conf ( dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
debug_assert( conformation::symmetry::is_symmetric( symm_conf ) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	// setup the minimizer map using the semi-symetric min map
	SymMinimizerMap sym_min_map( pose, semisym_move_map, symm_info, new_sym_min );
	//kinematics::DomainMap const & dm( sym_min_map.domain_map() );
	//for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
	//	std::cout << "(" << ii << ", " << dm(ii) << "), ";
	//	if ( ii % 10 == 0 ) std::cout << std::endl;
	//}
	//std::cout << std::endl;

	// if we are using the nblist, set it up
	if ( use_nblist ) {
		// setup a mask of the moving dofs
		pose.energies().set_use_nblist( pose, sym_min_map.domain_map(), options.nblist_auto_update() );
	}

	// etable functions must be initialized with asymm domain map
	scorefxn.setup_for_minimizing( pose, sym_min_map );

	// setup the function that we will pass to the low-level minimizer
	SymAtomTreeMultifunc f( pose, sym_min_map, scorefxn,
		options.deriv_check(), options.deriv_check_verbose() );

	// starting position -- "dofs" = Degrees Of Freedom
	Multivec dofs( sym_min_map.nangles() );
	sym_min_map.copy_dofs_from_pose( pose, dofs );

	Real const start_func( f( dofs ) );

	// std::cout << "Start score: " << start_score << " start_func " << start_func << std::endl;
	// APL Note: there might be a discrepancy between start_score and start_func if bb/bb hydrogen bonds
	// are being used.  Why?  Well, during minimization, bb/bb hbonds are calculated in the residue_pair_energy_ext
	// calls, but in regular scoring, bb/bb hbonds are calculated in setup_for_scoring.  This means that
	// bb/bb hbond energies are not stored in the EnergyGraph. Hbonding residue pairs that are not moving wrt
	// each other are not rescored during minimization, so their energies would be stored in the "fixed_energies_"
	// EnergyMap in the MinimizationGraph -- except for the fact that this EnergyMap is filled using energies
	// stored in the EnergyGraph, and the bb/bb hbond energies are not stored there.
	// The delta between start_score and start_func is not worrisome, however, because these energies are constant.
	// The delta has no impact on the minimizer's behavior.

	// now do the optimization with the low-level minimizer function
	Minimizer minimizer( f, options );
	minimizer.run( dofs );

	Real const end_func( f( dofs ) );

	//std::cout << "end_func:" << std::endl;
	//pose.energies().show( std::cout );

	// turn off nblist
	if ( use_nblist ) pose.energies().reset_nblist();

	// if we were doing rigid-body minimization, fold the rotation and
	// translation offsets into the jump transforms
	//
	// also sets rb dofs to 0.0, so in principle func value should be the same
	//
	sym_min_map.reset_jump_rb_deltas( pose, dofs );

	// rescore
	Real const end_score( scorefxn( pose ) );

	//std::cout << "end_score:" << std::endl;
	//pose.energies().show( std::cout );

	// we may not really need all these extra function evaluations
	// good for diagnostics though

	basic::Tracer core_optimize( "core.optimize", basic::t_debug );
	core_optimize << "SymAtomTreeMinimizer::run: nangles= " << sym_min_map.nangles() <<
		" start_score: " << F(12,3,start_score) <<
		" start_func: "  << F(12,3,start_func ) <<
		" end_score: "   << F(12,3,end_score  ) <<
		" end_func: "    << F(12,3,end_func   ) << std::endl;


	return end_score;
}

void
SymAtomTreeMinimizer::make_semisymmetric_movemap(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map_sym,
	kinematics::MoveMap & move_map_semisym
)
{
	using namespace core::conformation::symmetry;
	//typedef SymmetryInfo::DOF_IDs DOF_IDs;

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
debug_assert( conformation::symmetry::is_symmetric( SymmConf ) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	// copy bb/chi DOFs from symm movemap
	move_map_semisym = move_map_sym;

	for ( Size i=1; i<= pose.conformation().fold_tree().num_jump(); ++i ) {
		id::DOF_ID const & null_id
		( pose.conformation().dof_id_from_torsion_id(id::TorsionID(i,id::JUMP,1)));
		if ( symm_info->get_dof_derivative_weight( null_id , SymmConf ) > 0 ) {
			// if this is not the master get the master
			if ( symm_info->jump_is_independent( i ) ) continue;

			Size master_i = symm_info->jump_follows( i );

			for ( int j=1; j<= 6; ++j ) {
			id::DOF_ID const & id_master
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(master_i,id::JUMP,j)));
				id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(i,id::JUMP,j)));

				bool allow ( move_map_sym.get( id_master ) );
				move_map_semisym.set(id, allow );
			}
		}
	}

}


/// this routine will miss dofs that are set individually as opposed to the set_bb set_jump, etc,  instructions
void
SymAtomTreeMinimizer::make_assymetric_movemap(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map_sym,
	kinematics::MoveMap & move_map_asym
)
{
	using namespace core::conformation::symmetry;
	typedef SymmetryInfo::DOF_IDs DOF_IDs;

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
debug_assert( conformation::symmetry::is_symmetric( SymmConf ) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	for ( Size i=1; i<= pose.conformation().size(); ++i ) {
		if ( symm_info->bb_is_independent( i ) ) {
			bool bb ( move_map_sym.get_bb(i) );
			bool chi ( move_map_sym.get_chi(i) );
			move_map_asym.set_bb ( i, bb );
			move_map_asym.set_chi( i, chi );
			for ( std::vector< Size>::const_iterator
					clone     = symm_info->bb_clones( i ).begin(),
					clone_end = symm_info->bb_clones( i ).end();
					clone != clone_end; ++clone ){
				move_map_asym.set_bb ( *clone, bb );
				move_map_asym.set_chi( *clone, chi );
			}
		}
	}
	for ( Size i=1; i<= pose.conformation().fold_tree().num_jump(); ++i ) {
		if ( symm_info->jump_is_independent( i ) ) {
			for ( int j=1; j<= 6; ++j ) {
				id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(i,id::JUMP,j)));
				DOF_IDs const & dofs( symm_info->dependent_dofs( id, SymmConf ) );
				bool allow ( move_map_sym.get( id ) );
				move_map_asym.set(id, allow );
				for ( DOF_IDs::const_iterator dof =dofs.begin(), dofe= dofs.end(); dof != dofe; ++dof ) {
					move_map_asym.set( *dof, allow );
				}
			}
		}
	}
}

} // symmetry
} // namespace optimization
} // namespace core
