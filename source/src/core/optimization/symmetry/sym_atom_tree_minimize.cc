// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/optimization/atom_tree_minimize.hh
/// @brief  Atom tree minimization functions
/// @author Ingemar Andre

// Unit headers
#include <core/optimization/symmetry/sym_atom_tree_minimize.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/atom_tree_minimize.hh>

// Symmetry headers
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricEnergies.hh>
// AUTO-REMOVED #include <core/scoring/NeighborList.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/hbonds.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondSet.hh>

// // ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

// // Numeric headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>

// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/optimization/MinimizerMap.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>



using basic::T;
using basic::Error;
using basic::Warning;

using namespace ObjexxFCL::format;

namespace core {
namespace optimization {
namespace symmetry {

static basic::Tracer TR("core.optimization.symmetry");

typedef id::DOF_ID DOF_ID;


/////////////////////////////////////////////////////////////////////////////
/// @detailed
///car note that this calculates the deriv for all torsion angles even
///car those that are fixed. Because of the way that the derivative is
///car calculated, I don't believe this is a significant slow down (ie
///car have to run over all the atom pairs twice, regardless of the number
///car of torsion angles)
///
///car multiple neighborlists:
///car cendist       centroid distances in current structure, cutoff for vdw
///car dis2_tether   centroid distances in tether structure, cutoff for tether
///
///db computes the derivative of E  with respect to each
///db of the torsion angles.  Using the chain rule, we have
///db
///db dE/d phi = dE/dr *  dr/dphi
///db
///db dr/dphi  = Eab x (V-Vb) . (V' - V)/|V-V'|
///db
///db (the first cross product is the displacement of V upon a rotation dphi
///db around the unit vector Eab, Vb is the coordinates of the second atom in
///db the bond)
///db
///car dE/dR = 2r  (for vdw at least)
///db since | V-V'| = r,
///db
///db dE/ dphi = 2 Eab x (V-Vb) . (V' - V)
///db
///db  note that Eab and Vb are different for each torsion angle, but V'
///db and V are the same.  rearranging:
///db
///db = -  2 Eab X Vb . (V' - V) - 2 Eab . (V' x V).
///db
///db now we need the averages over all Vi of the difference and the
///db crossproduct of V and V'.
///
///car below, Eab x Vb is 'vec'
///car        Eab      is 'unit'
///car        (V'-V)   is 'f2'
///car        'F2tot' = f2*dE_dR (cumulative)
///car        (V' X V) is 'f1' ('F1_xxxE' is cumulative for potential xxx)
///car        eval_dE_dR actually returns dE_dR/r
///
///car if two atoms are fixed relatively in cartesian space, then dr/dphi = 0
///car and there is no contribution to the derivative
///

void
atom_tree_dfunc(
	pose::Pose & pose,
	SymMinimizerMap & symm_min_map,
	/*MinimizerMap & semisymm_min_map,
	MinimizerMap & asymm_min_map,*/
	scoring::ScoreFunction const & scorefxn,
	Multivec const & vars,
	Multivec & dE_dvars
)
{
	using namespace conformation::symmetry;
	// Initialize symmetry
	assert (pose::symmetry::is_symmetric( pose ) );
	SymmetricConformation & symm_conf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	//Real const deg2rad( numeric::conversions::radians(1.0) );

	dE_dvars.resize( symm_min_map.nangles() );

	// clear all the F1's and F2's
	symm_min_map.zero_torsion_vectors();


	// puts the degrees of freedom from vars into pose
	symm_min_map.copy_dofs_to_pose( pose, vars );


	/////////////////////////////////////////////////////////////////////////////
	// do some pre-computation prior to looping over the torsions

	// this will stash necessary information in the pose's energies object
	//
	scorefxn.setup_for_derivatives( pose );

	/////////////////////////////////////////////////////////////////////////////
	// get derivative of all atom pair potentials
	// this includes fa_pair and hbonds
	//
	// this call fills the F1's and F2's with contributions from their
	// immediately downstream atoms
	atom_tree_get_atompairE_deriv( pose, symm_min_map, scorefxn );


	/////////////////////////////////////////////////////////////////////////////
	// this should only be done once, after all torsion F1,F2's have
	// been filled in
	//
	// this sums all the F1,F2 contributions down the tree from leaves to root
	symm_min_map.link_torsion_vectors();


	/////////////////////////////////////////////////////////////////////////////
	// now loop over the torsions in the map
	int imap( 1 ); // for indexing into de_dvars( imap )
	for ( SymMinimizerMap::const_iterator it=symm_min_map.begin(), ite=symm_min_map.end();
			it != ite; ++it, ++imap ) {

		DOF_Node const & dof_node( **it );
		kinematics::tree::Atom const & atom( pose.atom_tree().atom( dof_node.atom_id() ) );
		/////////////////////////////////////////////////////////////////
		// derivatives of this particular degree of freedom
		//
		// eg rama,Paa,dunbrack,and torsional constraints
		Real sfxn_dof_deriv = scorefxn.eval_dof_derivative( dof_node.dof_id(), dof_node.torsion_id(), pose );
		Real scale = symm_min_map.torsion_scale_factor( dof_node ) / symm_info->score_multiply_factor();

		dE_dvars[ imap ] = torsional_derivative_from_cartesian_derivatives( atom, dof_node, sfxn_dof_deriv, scale );

	} // loop over map


	scorefxn.finalize_after_derivatives( pose );
}


///////////////////////////////////////////////////////////////////////////////

void
atom_tree_get_atompairE_deriv(
	pose::Pose & pose,
	SymMinimizerMap & symm_min_map,
	//MinimizerMap & semisymm_min_map,
	//MinimizerMap & asymm_min_map,
	scoring::ScoreFunction const & scorefxn
)
{
	using namespace conformation::symmetry;
	using namespace scoring;
	using namespace scoring::symmetry;

	SymmetricConformation const & symm_conf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	assert( conformation::symmetry::is_symmetric( symm_conf ) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );


	SymmetricEnergies const & symm_energies( dynamic_cast< SymmetricEnergies const & > (pose.energies()) );

	assert( symm_energies.minimization_graph() );
	assert( symm_energies.derivative_graph() );

	MinimizationGraphCOP mingraph  = symm_energies.minimization_graph();
	MinimizationGraphCOP dmingraph = symm_energies.derivative_graph();

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		MinimizationNode const & minnode =  * mingraph->get_minimization_node( ii );
		/// 1. eval intra-residue derivatives
		eval_atom_derivatives_for_minnode( minnode, pose.residue( ii ), pose, scorefxn.weights(), symm_min_map.atom_derivatives( ii ) );
	}

	/// 2a. eval inter-residue derivatives from the regular minimization graph
	for ( graph::Node::EdgeListConstIter
			edgeit = mingraph->const_edge_list_begin(), edgeit_end = mingraph->const_edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
		Size const rsd1ind = minedge.get_first_node_ind();
		Size const rsd2ind = minedge.get_second_node_ind();
		conformation::Residue const & rsd1( pose.residue( rsd1ind ));
		conformation::Residue const & rsd2( pose.residue( rsd2ind ));
		ResSingleMinimizationData const & r1_min_data( mingraph->get_minimization_node( rsd1ind )->res_min_data() );
		ResSingleMinimizationData const & r2_min_data( mingraph->get_minimization_node( rsd2ind )->res_min_data() );

		eval_weighted_atom_derivatives_for_minedge( minedge, rsd1, rsd2,
			r1_min_data, r2_min_data, pose, scorefxn.weights(),
			symm_min_map.atom_derivatives( rsd1ind ), symm_min_map.atom_derivatives( rsd2ind ));
	}

	/// 2b. eval inter-residue derivatives from derivative minimization graph
	for ( graph::Node::EdgeListConstIter
			edgeit = dmingraph->const_edge_list_begin(), edgeit_end = dmingraph->const_edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
		Size const rsd1ind = minedge.get_first_node_ind();
		Size const rsd2ind = minedge.get_second_node_ind();
		conformation::Residue const & rsd1( pose.residue( rsd1ind ));
		conformation::Residue const & rsd2( pose.residue( rsd2ind ));
		ResSingleMinimizationData const & r1_min_data( dmingraph->get_minimization_node( rsd1ind )->res_min_data() );
		ResSingleMinimizationData const & r2_min_data( dmingraph->get_minimization_node( rsd2ind )->res_min_data() );

		eval_weighted_atom_derivatives_for_minedge( minedge, rsd1, rsd2,
			r1_min_data, r2_min_data, pose, scorefxn.weights(),
			symm_min_map.atom_derivatives( rsd1ind ), symm_min_map.atom_derivatives( rsd2ind ));
	}


	//std::cerr << pose.fold_tree();

	// Loop over all dofs in the symmetric movemap
	// use atom lists from the semisymmetric movemap
	//
	for ( MinimizerMap::const_iterator iter = symm_min_map.begin(), iter_e = symm_min_map.end();
			iter != iter_e; ++iter ) {
		DOF_Node & dof_node( **iter );

		core::Real dof_wt_i =  symm_info->get_dof_derivative_weight( dof_node.dof_id(), symm_conf );
		// if( dof_wt_i == 0 ) continue;
		//std::cout << "  dof_node: " << dof_node.rsd() << " " << dof_node.atomno() << " " << dof_node.type() << " " << dof_wt_i << std::endl;

		Vector f1(0,0,0), f2(0,0,0);
		// loop through atoms first moved by this torsion
		for ( DOF_Node::AtomIDs::const_iterator it1=dof_node.atoms().begin(),
				it1e = dof_node.atoms().end();	it1 != it1e; ++it1 ) {
			id::AtomID const & atom_id( *it1 );

			/// Most of the derivative evaluation has already taken place by the time we get here.
			dof_node.F1() += symm_min_map.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f1();
			dof_node.F2() += symm_min_map.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f2();

			scorefxn.eval_npd_atom_derivative( atom_id, pose, symm_min_map.domain_map(), dof_node.F1(), dof_node.F2() );

		} // atom1

		//std::cout << "   ... summing " << dof_node.atom_id() << " with weight " << dof_wt_i << std::endl;
		//std::cout << "   " << f1.x() << " " << f1.y() << " " << f1.z() << " " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;

		dof_node.F1() += dof_wt_i * f1;
		dof_node.F2() += dof_wt_i * f2;
	}
	for ( MinimizerMap::const_iterator iter = symm_min_map.dependent_begin(), iter_e = symm_min_map.dependent_end();
			iter != iter_e; ++iter ) {
		DOF_Node & dof_node( **iter );

		core::Real dof_wt_i =  symm_info->get_dof_derivative_weight( dof_node.dof_id(), symm_conf );
		// if( dof_wt_i == 0 ) continue;
		//std::cout << "  dof_node: " << dof_node.rsd() << " " << dof_node.atomno() << " " << dof_node.type() << " " << dof_wt_i << std::endl;

		Vector f1(0,0,0), f2(0,0,0);
		// loop through atoms first moved by this torsion
		for ( DOF_Node::AtomIDs::const_iterator it1=dof_node.atoms().begin(),
				it1e = dof_node.atoms().end();	it1 != it1e; ++it1 ) {
			id::AtomID const & atom_id( *it1 );

			/// Most of the derivative evaluation has already taken place by the time we get here.
			f1 += symm_min_map.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f1();
			f2 += symm_min_map.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f2();

			scorefxn.eval_npd_atom_derivative( atom_id, pose, symm_min_map.domain_map(), f1, f2 );

		} // atom1

		Size this_jump = pose.fold_tree().get_jump_that_builds_residue( dof_node.dof_id().rsd() );
		Size master_jump = symm_info->jump_follows( this_jump );
		DOF_ID symm_dof_id(
			id::AtomID( dof_node.atomno(), pose.fold_tree().downstream_jump_residue( master_jump )  ), dof_node.type() );
		// get the equiv node in the symm min map
		DOF_NodeOP symm_dof_node = symm_min_map.dof_node_from_id( symm_dof_id );

		//std::cout << "   ... summing " << dof_node.atom_id() << " at dof node " << symm_dof_node->atom_id() << " with weight " << dof_wt_i << std::endl;
		//std::cout << "   ... summing " << dof_node.atom_id() << " at dof node " << symm_dof_node->atom_id() << " with weight " << dof_wt_i << std::endl;
		//std::cout << "   " << f1.x() << " " << f1.y() << " " << f1.z() << " " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;

		// add f1 f2 using dof deriv wt
		symm_dof_node->F1() += dof_wt_i * f1;
		symm_dof_node->F2() += dof_wt_i * f2;

	} // tor
}

///////////////////////////////////////////////////////////////////////////////
// temporary!
class MinDebug {
public:
	MinDebug( Size const nangles ):
		abs_deriv_dev( nangles ),
		rel_deriv_dev( nangles )
	{}

	utility::vector1< Real > abs_deriv_dev;
	utility::vector1< Real > rel_deriv_dev;

	Real best_cos_theta;
	Real best_abs_log_norm_ratio;
	Real best_norm_analytic;
	Real best_norm_numeric;
};

///////////////////////////////////////////////////////////////////////////////

void
numerical_derivative_check(
	SymMinimizerMap const & min_map,
	Multifunc const & func,
	Multivec const & start_vars,
	Multivec const & dE_dvars,
	bool const verbose // = true
)
{
	/////////////////////////////////////////////////////////////////////////////
	// NUMERICAL DERIVATIVE CHECK
	/////////////////////////////////////////////////////////////////////////////
	// how to analyze this:
	//
	// in gnuplot, look at numerical vs analytical derivs:
	// plot '< grep "^ratio" a3.log ' u 10:11,x
	//
	// also sort by deriv_dev lines:
	//
	// by magnitude of deviation
	// grep deriv_dev a3.log | sort -g +8
	//
	// or by ratio of deviation to actual
	//
	// grep deriv_dev a3.log | sort -g +9

	Size const nangles( min_map.nangles() );


	Real const increment = 0.0005; // PB -- 3/02
	Size const n_increment = 5;
	utility::vector1< Multivec > dE_dvars_numeric( n_increment );
	for ( Size i=1; i<= n_increment; ++i ) {
		dE_dvars_numeric[i].resize( nangles, 0.0 );
	}
	//FArray2D_Real dE_dvars_numeric( nangles, n_increment );

	// setup for saving diagnostics
	MinDebug min_debug( nangles );

// 	min_debug.nangles = nangles;
// 	if ( nangles > int( min_debug.abs_deriv_dev.size1() ) ) {
// 		min_debug.abs_deriv_dev.dimension( nangles );
// 		min_debug.rel_deriv_dev.dimension( nangles );
// 	}

	Multivec vars( start_vars );

	Real const f00 = func( vars );
	Size ii( 1 ); // for indexing into dE_dvars[ ii ]

	for ( MinimizerMap::const_iterator iter= min_map.begin(),
					iter_end= min_map.end(); iter != iter_end; ++iter, ++ii ) {
		DOF_Node const & dof_node( **iter );

		Real deriv_dev = 10000.0;
		for ( Size j = 1,factor=1; j <= n_increment; ++j ) {
			factor*=2;

			vars[ii] = start_vars[ii] + factor * increment;
			Real const f11 = func( vars );

			vars[ii] = start_vars[ii] - factor * increment;
			Real const f22 = func( vars );

			Real const deriv = ( f11 - f22 ) / ( factor * 2 * increment );

			dE_dvars_numeric[j][ii] = deriv;

			deriv_dev = std::min( deriv_dev, std::abs( deriv  - dE_dvars[ii] ) );

			vars[ii] = start_vars[ii];

			Real const ratio( std::abs( dE_dvars[ii] ) < 0.001 ? 0.0 :
												 deriv / dE_dvars[ii] );

			if ( verbose &&
					 ( std::abs(dE_dvars[ii]) > 0.001 || std::abs(deriv) > 0.001 ) ) {
				// if you change this output, please also change the comments
				// at the beginning of this section
				static bool ratio_header_output( false );
				if ( !ratio_header_output ) {
					ratio_header_output = true;
					TR << "ratio" <<
						A( 4, "inc" ) <<
						A( 4, "rsd" ) <<
						A( 4, "typ" ) <<
						A( 4, "atm" ) <<
						A( 5, "prsd" ) <<
						A( 5, "ptyp" ) <<
						A( 5, "patm" ) <<
						A( 5, "natm" ) <<
						A( 10, "numeric" ) <<
						A( 10, "analytic" ) <<
						A( 10, "ratio" ) <<
						A( 10, "f11" ) <<
						A( 10, "f00" ) <<
						A( 10, "f22" ) <<
						A( 10, "vars[ii]" ) << std::endl;
				}

				id::DOF_ID parent_id( id::BOGUS_DOF_ID );
				if ( dof_node.parent() ) {
					parent_id = dof_node.parent()->dof_id();
				}

				TR << "ratio" <<
					I( 4, j ) <<
					I( 4, dof_node.rsd() ) <<
					I( 4, dof_node.type() ) <<
					I( 4, dof_node.atomno() ) <<
					I( 5, parent_id.rsd() ) <<
					I( 5, parent_id.type() ) <<
					I( 5, parent_id.atomno() ) <<
					I( 5, dof_node.atoms().size()) <<
					F( 10, 4, deriv ) <<
					F( 10, 4, dE_dvars[ii] ) <<
					F( 10, 4, ratio ) <<
					F( 10, 4, f11 ) <<
					F( 10, 4, f00 ) <<
					F( 10, 4, f22 ) <<
					F( 10, 4, start_vars[ii] ) << std::endl;
			}
		}
		if ( true ) {

			Real const ratio( std::abs( dE_dvars[ii] ) < 0.001 ? 0.0 :
												 deriv_dev / std::abs( dE_dvars[ii] ) );

			min_debug.rel_deriv_dev[ ii ] = ratio;
			min_debug.abs_deriv_dev[ ii ] = deriv_dev;

			if ( verbose ) {
				// if you change this output, please also change the comments
				// at the beginning of this section
				TR << "deriv_dev:" << SS(ii) << SS(nangles) << SS(f00) <<
					SS( dof_node.type() ) << SS( dof_node.atomno() ) <<
					SS( dof_node.rsd() ) <<
					SS( dE_dvars[ii] ) << SS( deriv_dev ) << SS(ratio) << std::endl;
			}
		}
	}

	// calculate magnitudes, dot products of gradient vectors
	utility::vector1< Real > norm_numeric(n_increment,0.0), dot(n_increment,0.0);
	Real norm(0.0);
	for ( Size i=1; i<= nangles; ++i ) {
		norm += dE_dvars[i] * dE_dvars[i];
		for ( Size j=1; j<= n_increment; ++j ) {
			dot[j] += dE_dvars[i] * dE_dvars_numeric[j][i];
			norm_numeric[j] += dE_dvars_numeric[j][i] * dE_dvars_numeric[j][i];
		}
	}
	norm = std::sqrt( norm );

	min_debug.best_cos_theta = -10.0;
	min_debug.best_abs_log_norm_ratio = 200.0;
	min_debug.best_norm_analytic = 999.9;
	min_debug.best_norm_numeric  = 999.9;

	for ( Size j=1; j<= n_increment; ++j ) {
		norm_numeric[j] = std::sqrt( norm_numeric[j] );

		// handle strange cases
		Real log_norm_ratio;
		if ( norm < 0.001 && norm_numeric[j] < 0.001 ) {
			log_norm_ratio = 1.0;
		} else if ( norm < 0.001 ) {
			log_norm_ratio = 100.0;
		} else if ( norm_numeric[j] < 0.001 ) {
			log_norm_ratio = -100.0;
		} else {
			log_norm_ratio = std::log( norm_numeric[j] / norm );
		}

		Real const cos_theta( dot[j] / ( norm * norm_numeric[j]) );

		TR <<
			" norm: " << j << ' ' << F(12,4,norm) <<
			" norm_numeric: " << F(12,4,norm_numeric[j]) <<
			" cos_theta: " << F(7,4,cos_theta) <<
			" log_norm_ratio: " << F(9,4,log_norm_ratio) << std::endl;

		min_debug.best_cos_theta = std::max( min_debug.best_cos_theta,
																					cos_theta );
		if ( std::abs( log_norm_ratio ) < min_debug.best_abs_log_norm_ratio ) {
			min_debug.best_abs_log_norm_ratio = std::abs( log_norm_ratio );
			min_debug.best_norm_analytic = norm;
			min_debug.best_norm_numeric = norm_numeric[j];
		}
	}
}

} // namespace symmetry
} // namespace optimization
} // namespace core
// 	{ // dont want to modify the pose if we can avoid it
// 		Size const nangles( min_map.nangles() );
// 		Multivec tmp_vars( nangles );
// 		min_map.copy_dofs_from_pose( pose, tmp_vars );
// 		Real dev(0.0);
// 		for ( Size i=1; i<= nangles; ++i ) {
// 			dev += std::abs( tmp_vars[i] - vars[i] );
// 		}
// 		std::cout << "[ DEBUG ] vars dev in atom_tree_dfunc: " << dev << std::endl;
// 		if ( dev > 1e-2 ) {
// 			min_map.copy_dofs_to_pose( pose, vars );
// 			scorefxn( pose );
// 		}
// 	}
