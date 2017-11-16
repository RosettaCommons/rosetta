// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/atom_tree_minimize.cc
/// @brief  Atom tree minimization functions
/// @author Phil Bradley


// Unit headers
#include <core/optimization/atom_tree_minimize.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

// Project headers
#include <core/kinematics/AtomTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationGraph.hh>

// // ObjexxFCL headers
// #include <ObjexxFCL/FArray1A.hh>
// #include <ObjexxFCL/FArray3A.hh>
// #include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/format.hh>
// #include <ObjexxFCL/string.functions.hh>

// // Numeric headers
#include <numeric/constants.hh>
// #include <numeric/xyzVector.hh>
// #include <numeric/xyz.functions.hh>

// // Utility headers
// #include <utility/exit.hh>

// // C++ headers
// #include <algorithm>
// #include <cmath>
// #include <cstdlib>
// #include <iostream>

#include <basic/Tracer.hh>

#include <core/kinematics/tree/Atom.hh>
#include <core/optimization/MinimizerMap.hh>
#include <utility/vector1.hh>

using basic::Error;
using basic::Warning;
using namespace ObjexxFCL::format;

namespace core {
namespace optimization {

static basic::Tracer TR( "core.optimization" );

/////////////////////////////////////////////////////////////////////////////
/// @details
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


void
atom_tree_dfunc(
	pose::Pose & pose,
	MinimizerMap & min_map,
	scoring::ScoreFunction const & scorefxn,
	Multivec const & vars,
	Multivec & dE_dvars
)
{

	dE_dvars.resize( min_map.nangles() );

	// clear all the F1's and F2's
	min_map.zero_torsion_vectors();


	// puts the degrees of freedom from vars into pose
	min_map.copy_dofs_to_pose( pose, vars );


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
	atom_tree_get_atompairE_deriv( pose, min_map, scorefxn );


	/////////////////////////////////////////////////////////////////////////////
	// this should only be done once, after all torsion F1,F2's have
	// been filled in
	//
	// this sums all the F1,F2 contributions down the tree from leaves to root
	min_map.link_torsion_vectors();


	/////////////////////////////////////////////////////////////////////////////
	// now loop over the torsions in the map
	int imap( 1 ); // for indexing into de_dvars( imap )
	for ( auto it=min_map.begin(), ite=min_map.end();
			it != ite; ++it, ++imap ) {

		DOF_Node const & dof_node( **it );
		kinematics::tree::Atom const & atom( pose.atom_tree().atom( dof_node.atom_id() ) );
		/////////////////////////////////////////////////////////////////
		// derivatives of this particular degree of freedom
		//
		// eg rama,Paa,dunbrack,and torsional constraints
		Real sfxn_dof_deriv = scorefxn.eval_dof_derivative( dof_node.dof_id(), dof_node.torsion_id(), pose );
		Real scale = min_map.torsion_scale_factor( dof_node );

		dE_dvars[ imap ] = torsional_derivative_from_cartesian_derivatives( atom, dof_node, sfxn_dof_deriv, scale );

	} // loop over map


	scorefxn.finalize_after_derivatives( pose );

}

/// @details Refactored from atom_tree_dfunc above for use in alternate minimization contexts
Real
torsional_derivative_from_cartesian_derivatives(
	kinematics::tree::Atom const & atom,
	optimization::DOF_Node const & dof_node,
	Real dof_deriv, // derivatives applied directly to this DOF, if any
	Real torsion_scale_factor
)
{
	using namespace id;

	// NOTE: deriv is in the units of the degree of freedom as
	// represented internally, without any scale factors applied
	// ie the units returned by pose.get_atom_tree_torsion(...)
	//
	//
	// type  -- units
	// -------------------
	// PHI   -- radians
	// THETA -- radians
	// D     -- angstroms
	// RB1-3 -- angstroms
	// RB4-6 -- degrees    (!)

	Real deriv = 0.0;

	// this function determines the axis and center of rotation/translation
	Vector axis, end_pos;
	DOF_Type const type( dof_node.type() );
	atom.get_dof_axis_and_end_pos( axis, end_pos, type );


	// convert combined F1,F2 to an angular derivative
	//  using Abe Go trick
	//
	if ( type == PHI || type == THETA || type == RB4 || type == RB5 ||
			type == RB6 ) {
		// rotation about an axis
		// note: assumes we are dealing with RADIANS
		// scale factor below handles the fact that RB4-6 are degrees
		Real scale_factor( ( type == PHI || type == THETA ) ? 1 : numeric::constants::d::deg2rad );

		if ( type == THETA ) {
			// need to think about this more carefully
			using numeric::constants::f::pi; // silly -- why not numeric::constants::d?
			Real const theta( atom.dof( type ) );
			int const theta_mod
				( ( static_cast< int >( std::floor( theta/pi )))%2);
			if ( theta_mod == 1 || theta_mod == -1 ) {
				scale_factor *= -1.0f;
			}
		}

		deriv -= scale_factor * ( dot( axis, dof_node.F1() ) +
			dot( cross( axis, end_pos ), dof_node.F2() ) );
	} else {
		// translation along an axis
		deriv += dot( axis, dof_node.F2() );
	}

	deriv += dof_deriv;

	deriv /= torsion_scale_factor;

	return deriv;
}

///////////////////////////////////////////////////////////////////////////////

/// @details First evaluate all derivatives for the residue- and residue-pair decomposable
/// terms, using the minimization graph as a guide.  Then, evaluate any further derivatives
/// for those terms that are not residue- or residue-pair decomposable
void
atom_tree_get_atompairE_deriv(
	pose::Pose & pose,
	MinimizerMap & min_map,
	scoring::ScoreFunction const & scorefxn
)
{
	using namespace scoring;

	debug_assert( pose.energies().minimization_graph() );
	MinimizationGraphCOP mingraph = pose.energies().minimization_graph();

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		MinimizationNode const & minnode =  * mingraph->get_minimization_node( ii );
		/// 1. eval intra-residue derivatives
		eval_atom_derivatives_for_minnode( minnode, pose.residue( ii ), pose, scorefxn.weights(), min_map.atom_derivatives( ii ) );
	}

	/// 2. eval inter-residue derivatives
	for ( utility::graph::Node::EdgeListConstIter
			edgeit = mingraph->const_edge_list_begin(), edgeit_end = mingraph->const_edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
		Size const rsd1ind = minedge.get_first_node_ind();
		Size const rsd2ind = minedge.get_second_node_ind();
		conformation::Residue const & rsd1( pose.residue( rsd1ind ));
		conformation::Residue const & rsd2( pose.residue( rsd2ind ));
		ResSingleMinimizationData const & r1_min_data( mingraph->get_minimization_node( rsd1ind )->res_min_data() );
		ResSingleMinimizationData const & r2_min_data( mingraph->get_minimization_node( rsd2ind )->res_min_data() );

		eval_atom_derivatives_for_minedge( minedge, rsd1, rsd2,
			r1_min_data, r2_min_data, pose, scorefxn.weights(),
			min_map.atom_derivatives( rsd1ind ), min_map.atom_derivatives( rsd2ind ));
	}

	for ( auto iter = min_map.begin(), iter_e = min_map.end();
			iter != iter_e; ++iter ) {
		DOF_Node & dof_node( **iter );

		// loop through atoms first moved by this torsion
		for ( auto it1=dof_node.atoms().begin(),
				it1e = dof_node.atoms().end(); it1 != it1e; ++it1 ) {
			id::AtomID const & atom_id( *it1 );
			dof_node.F1() += min_map.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f1();
			dof_node.F2() += min_map.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f2();
			scorefxn.eval_npd_atom_derivative( atom_id, pose, min_map.domain_map(), dof_node.F1(), dof_node.F2() );
		} // atom1
	} // tor
}


SimpleDerivCheckResult
simple_numeric_deriv_check(
	Multifunc const & func,
	Multivec const & start_vars,
	Multivec const & dE_dvars,
	bool send_to_stdout,
	bool verbose,
	Size nsteps /* = 5 */
)
{

	/////////////////////////////////////////////////////////////////////////////
	// NUMERICAL DERIVATIVE CHECK
	/////////////////////////////////////////////////////////////////////////////
	// how to analyze this:
	//
	// in gnuplot, look at numerical vs analytical derivs:
	// plot '< grep "^core.optimization: ratio" a3.log ' u 4:5,x
	//
	// also sort by deriv_dev lines:
	//
	// by magnitude of deviation
	// grep deriv_dev a3.log | sort -g +6
	//
	// or by ratio of deviation to actual
	//
	// grep deriv_dev a3.log | sort -g +7

	Size const nangles( start_vars.size() );


	Real const increment = 0.0005; // PB -- 3/02
	Size const n_increment = nsteps;
	utility::vector1< Multivec > dE_dvars_numeric( n_increment );

	SimpleDerivCheckResult num_deriv_result(nangles, n_increment);

	for ( Size i=1; i<= n_increment; ++i ) {
		dE_dvars_numeric[i].resize( nangles, 0.0 );
	}

	Multivec vars( start_vars );

	Real const f00 = func( vars );

	for ( Size ii = 1; ii <= vars.size(); ++ii ) {
		Real deriv_dev = 10000.0;
		for ( Size j = 1,factor=1; j <= n_increment; ++j ) {
			factor*=2;

			vars[ii] = start_vars[ii] + factor * increment;
			Real const f11 = func( vars );

			//std::cout << "Vars, f11:";
			//for ( Size jj = 1; jj <= vars.size(); ++jj ) std::cout << " " << vars[ jj ];
			//std::cout << " " << f11 << std::endl;

			vars[ii] = start_vars[ii] - factor * increment;
			Real const f22 = func( vars );

			//std::cout << "Vars, f22:";
			//for ( Size jj = 1; jj <= vars.size(); ++jj ) std::cout << " " << vars[ jj ];
			//std::cout << " " << f22 << std::endl;


			Real const deriv = ( f11 - f22 ) / ( factor * 2 * increment );

			//std::cout << "Numeric deriv " << deriv << " anayltic: " << dE_dvars[ ii ] << std::endl;

			dE_dvars_numeric[j][ii] = deriv;

			deriv_dev = std::min( deriv_dev, std::abs( deriv  - dE_dvars[ii] ) );

			vars[ii] = start_vars[ii];

			Real const ratio( std::abs( dE_dvars[ii] ) < 0.001 ? 0.0 : deriv / dE_dvars[ii] );

			if ( std::abs(dE_dvars[ii]) > 0.001 || std::abs(deriv) > 0.001 ) {
				DerivCheckDataPoint dp( deriv, dE_dvars[ii], ratio, f11, f00, f22, start_vars[ ii ] );
				num_deriv_result.step_data( ii, j, dp );

				if ( verbose && send_to_stdout ) {
					// if you change this output, please also change the comments
					// at the beginning of this section
					static bool ratio_header_output( false );
					if ( !ratio_header_output ) {
						ratio_header_output = true;
						TR << "ratio" <<
							A( 4, "inc" ) <<
							A( 10, "numeric" ) <<
							A( 10, "analytic" ) <<
							A( 10, "ratio" ) <<
							A( 10, "f11" ) <<
							A( 10, "f00" ) <<
							A( 10, "f22" ) <<
							A( 10, "vars[ii]" ) << std::endl;
					}

					TR << "ratio" <<
						I( 4, j ) <<
						F( 10, 4, deriv ) <<                // column 4
						F( 10, 4, dE_dvars[ii] ) <<         // column 5
						F( 10, 4, ratio ) <<
						F( 10, 4, f11 ) <<
						F( 10, 4, f00 ) <<
						F( 10, 4, f22 ) <<
						F( 10, 4, start_vars[ii] ) << std::endl;
				}

			}
		}
		if ( true ) {

			Real const ratio( std::abs( dE_dvars[ii] ) < 0.001 ? 0.0 :
				deriv_dev / std::abs( dE_dvars[ii] ) );

			num_deriv_result.rel_deriv_dev( ii, ratio );
			num_deriv_result.abs_deriv_dev( ii, deriv_dev );

			if ( verbose && send_to_stdout ) {
				// if you change this output, please also change the comments
				// at the beginning of this section
				TR << "deriv_dev:" << SS(ii) << SS(nangles) << SS(f00) <<
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

	Real const too_small( 0.0001 );

	num_deriv_result.best_cos_theta( -10.0 );
	num_deriv_result.best_abs_log_norm_ratio( 200.0 );
	num_deriv_result.best_norm_analytic( 999.9 );
	num_deriv_result.best_norm_numeric( 999.9 );

	for ( Size j=1; j<= n_increment; ++j ) {
		norm_numeric[j] = std::sqrt( norm_numeric[j] );

		// handle strange cases
		Real log_norm_ratio;
		if ( norm < too_small && norm_numeric[j] < too_small ) {
			log_norm_ratio = 1.0;
		} else if ( norm < 0.125 * too_small ) {
			log_norm_ratio = 100.0;
		} else if ( norm_numeric[j] < 0.125 * too_small ) {
			log_norm_ratio = -100.0;
		} else {
			log_norm_ratio = std::log( norm_numeric[j] / norm );
		}

		Real const cos_theta( dot[j] / ( norm * norm_numeric[j]) );

		if ( send_to_stdout ) {
			TR <<
				" norm: " << j << ' ' << F(12,4,norm) <<
				" norm_numeric: " << F(12,4,norm_numeric[j]) <<
				" cos_theta: " << F(7,4,cos_theta) <<
				" log_norm_ratio: " << F(9,4,log_norm_ratio) << std::endl;
		}

		num_deriv_result.best_cos_theta( std::max( num_deriv_result.best_cos_theta(), cos_theta ) );
		if ( std::abs( log_norm_ratio ) < num_deriv_result.best_abs_log_norm_ratio() ) {
			num_deriv_result.best_abs_log_norm_ratio( std::abs( log_norm_ratio ) ) ;
			num_deriv_result.best_norm_analytic( norm );
			num_deriv_result.best_norm_numeric( norm_numeric[j] );
		}
	}
	if ( send_to_stdout ) {
		TR << "best : cos_theta " << num_deriv_result.best_cos_theta()
			<< " log_norm_ratio: " << num_deriv_result.best_abs_log_norm_ratio()
			<< " norm " << num_deriv_result.best_norm_analytic()
			<< " norm numeric: " << num_deriv_result.best_norm_numeric() << std::endl;
	}
	//deriv_check_result->add_deriv_data( min_debug );
	return num_deriv_result;
}

///////////////////////////////////////////////////////////////////////////////

void
numerical_derivative_check(
	MinimizerMap const & min_map,
	Multifunc const & func,
	Multivec const & start_vars,
	Multivec const & dE_dvars,
	NumericalDerivCheckResultOP deriv_check_result,
	bool const verbose // = true
)
{
	/////////////////////////////////////////////////////////////////////////////
	// NUMERICAL DERIVATIVE CHECK
	/////////////////////////////////////////////////////////////////////////////
	// how to analyze this:
	//
	// in gnuplot, look at numerical vs analytical derivs:
	// plot '< grep "^core.optimization: ratio" a3.log ' u 11:12,x
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

	bool const write_to_stdout( ! deriv_check_result || deriv_check_result->send_to_stdout() );

	Real const increment = 0.0005; // PB -- 3/02
	Size const n_increment = 5;
	utility::vector1< Multivec > dE_dvars_numeric( n_increment );
	for ( Size i=1; i<= n_increment; ++i ) {
		dE_dvars_numeric[i].resize( nangles, 0.0 );
	}
	//FArray2D_Real dE_dvars_numeric( nangles, n_increment );

	// setup for saving diagnostics
	NumDerivCheckDataOP min_debug;
	if ( deriv_check_result ) min_debug = NumDerivCheckDataOP( new NumDerivCheckData( nangles, n_increment ) );

	//  min_debug->nangles = nangles;
	//  if ( nangles > int( min_debug->abs_deriv_dev.size1() ) ) {
	//   min_debug->abs_deriv_dev.dimension( nangles );
	//   min_debug->rel_deriv_dev.dimension( nangles );
	//  }

	Multivec vars( start_vars );

	Real const f00 = func( vars );
	Size ii( 1 ); // for indexing into dE_dvars[ ii ]

	for ( auto iter= min_map.begin(),
			iter_end= min_map.end(); iter != iter_end; ++iter, ++ii ) {
		DOF_Node const & dof_node( **iter );

		if ( dof_node.dependent() ) { --ii; continue; } // special case for symmetry

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

			if ( std::abs(dE_dvars[ii]) > 0.001 || std::abs(deriv) > 0.001 ) {
				id::DOF_ID parent_id( id::DOF_ID::BOGUS_DOF_ID() );
				if ( dof_node.parent() ) {
					parent_id = dof_node.parent()->dof_id();
				}
				DOF_DataPoint dp( dof_node.dof_id(), parent_id, dof_node.atoms().size(),
					deriv, dE_dvars[ii], ratio, f11, f00, f22, start_vars[ ii ] );
				if ( min_debug ) min_debug->dof_step_data( ii, j, dp );
				if ( verbose && write_to_stdout ) {
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


					TR << "ratio" <<
						I( 4, j ) <<
						I( 4, dof_node.rsd() ) <<
						I( 4, dof_node.type() ) <<
						I( 4, dof_node.atomno() ) <<
						I( 5, parent_id.rsd() ) <<
						I( 5, parent_id.type() ) <<
						I( 5, parent_id.atomno() ) <<
						I( 5, dof_node.atoms().size()) <<
						F( 10, 4, deriv ) <<                // column 11
						F( 10, 4, dE_dvars[ii] ) <<         // column 12
						F( 10, 4, ratio ) <<
						F( 10, 4, f11 ) <<
						F( 10, 4, f00 ) <<
						F( 10, 4, f22 ) <<
						F( 10, 4, start_vars[ii] ) << std::endl;
				}

			}
		}
		if ( true ) {

			Real const ratio( std::abs( dE_dvars[ii] ) < 0.001 ? 0.0 :
				deriv_dev / std::abs( dE_dvars[ii] ) );

			if ( min_debug ) min_debug->rel_deriv_dev( ii, ratio );
			if ( min_debug ) min_debug->abs_deriv_dev( ii, deriv_dev );

			if ( verbose && write_to_stdout ) {
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

	Real lbest_cos_theta( -10.0 );
	Real lbest_abs_log_norm_ratio( 200.0 );
	Real lbest_norm_analytic( 999.9 );
	Real lbest_norm_numeric( 999.9 );

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

		if ( write_to_stdout ) {
			TR <<
				" norm: " << j << ' ' << F(12,4,norm) <<
				" norm_numeric: " << F(12,4,norm_numeric[j]) <<
				" cos_theta: " << F(7,4,cos_theta) <<
				" log_norm_ratio: " << F(9,4,log_norm_ratio) << std::endl;
		}

		lbest_cos_theta = std::max( lbest_cos_theta, cos_theta );
		if ( std::abs( log_norm_ratio ) < lbest_abs_log_norm_ratio ) {
			lbest_abs_log_norm_ratio = std::abs( log_norm_ratio );
			lbest_norm_analytic = norm;
			lbest_norm_numeric = norm_numeric[j];
		}
	}

	if ( min_debug ) min_debug->best_cos_theta( lbest_cos_theta );
	if ( min_debug ) min_debug->best_abs_log_norm_ratio( lbest_abs_log_norm_ratio );
	if ( min_debug ) min_debug->best_norm_analytic( lbest_norm_analytic );
	if ( min_debug ) min_debug->best_norm_numeric( lbest_norm_numeric );

	if ( deriv_check_result ) deriv_check_result->add_deriv_data( min_debug );
}


} // namespace optimization
} // namespace core
//  { // dont want to modify the pose if we can avoid it
//   Size const nangles( min_map.nangles() );
//   Multivec tmp_vars( nangles );
//   min_map.copy_dofs_from_pose( pose, tmp_vars );
//   Real dev(0.0);
//   for ( Size i=1; i<= nangles; ++i ) {
//    dev += std::abs( tmp_vars[i] - vars[i] );
//   }
//   std::cout << "[ DEBUG ] vars dev in atom_tree_dfunc: " << dev << std::endl;
//   if ( dev > 1e-2 ) {
//    min_map.copy_dofs_to_pose( pose, vars );
//    scorefxn( pose );
//   }
//  }
