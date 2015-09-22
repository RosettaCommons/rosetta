// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/cartesian_minimize.cc
/// @brief  Atom tree minimization functions
/// @author Frank DiMaio


// Unit headers
#include <core/optimization/cartesian_minimize.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/scoring/symmetry/SymmetricEnergies.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationGraph.hh>

#include <ObjexxFCL/format.hh>

// // Numeric headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>

#include <basic/Tracer.hh>

#include <core/optimization/CartesianMinimizerMap.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace ObjexxFCL::format;

namespace core {
namespace optimization {

static THREAD_LOCAL basic::Tracer TR( "core.optimization" );

/////////////////////////////////////////////////////////////////////////////
/// @details


void
cartesian_dfunc(
	pose::Pose & pose,
	CartesianMinimizerMap & min_map,
	scoring::ScoreFunction const & scorefxn,
	Multivec const & vars,
	Multivec & dE_dvars
) {
	dE_dvars.resize( min_map.ndofs() );

	// clear stored F1's and F2's
	min_map.zero_stored_derivs();

	// puts the degrees of freedom from vars into pose
	min_map.copy_dofs_to_pose( pose, vars );

	// do some pre-computation prior to looping over the torsions
	// this will stash necessary information in the pose's energies object
	scorefxn.setup_for_derivatives( pose );

	//fpd  scale derivatives for symmetry
	core::Real scale = 1;
	if ( pose::symmetry::is_symmetric( pose ) ) {
		using namespace conformation::symmetry;
		SymmetricConformation & symm_conf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		scale = symm_info->score_multiply_factor();
	}

	// get derivative of all atom pair potentials
	// this includes fa_pair and hbonds
	//
	// this call fills the F1's and F2's with contributions from their
	// immediately downstream atoms
	cartesian_collect_atompairE_deriv( pose, min_map, scorefxn, dE_dvars, scale );

	// now loop over the torsions in the map
	//Size ntorsions=min_map.ntorsions();
	cartesian_collect_torsional_deriv( pose, min_map, scorefxn, dE_dvars, scale );

	scorefxn.finalize_after_derivatives( pose );
}

///////////////////////////////////////////////////////////////////////////////
void
cartesian_collect_atompairE_deriv(
	pose::Pose & pose,
	CartesianMinimizerMap & min_map,
	scoring::ScoreFunction const & scorefxn,
	Multivec & dE_dvars,
	core::Real scale
) {
	using namespace scoring;
	using namespace scoring::symmetry;
	using namespace conformation::symmetry;

	debug_assert( pose.energies().minimization_graph() );
	MinimizationGraphCOP mingraph = pose.energies().minimization_graph();

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		MinimizationNode const & minnode =  * mingraph->get_minimization_node( ii );
		/// 1. eval intra-residue derivatives
		eval_atom_derivatives_for_minnode( minnode, pose.residue( ii ), pose, scorefxn.weights(), min_map.atom_derivatives( ii ) );
	}

	/// 2. eval inter-residue derivatives
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
			min_map.atom_derivatives( rsd1ind ), min_map.atom_derivatives( rsd2ind ));
	}

	// if we're symmetric loop over other edges
	if ( pose::symmetry::is_symmetric( pose ) ) {
		SymmetricEnergies const & symm_energies( dynamic_cast< SymmetricEnergies const & > (pose.energies()) );
		MinimizationGraphCOP dmingraph = symm_energies.derivative_graph();

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
				min_map.atom_derivatives( rsd1ind ), min_map.atom_derivatives( rsd2ind ));
		}
	}

	Size natoms=min_map.natoms();
	for ( Size i=1; i<=natoms; ++i ) {
		id::AtomID const & atom_id( min_map.get_atom(i) );
		core::Vector F1(0,0,0),F2(0,0,0);
		scorefxn.eval_npd_atom_derivative( atom_id, pose, min_map.domain_map(), F1, F2 );
		//F1 += min_map.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f1();
		F2 += min_map.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f2();
		dE_dvars[3*i-2] = F2[0]*scale;
		dE_dvars[3*i-1] = F2[1]*scale;
		dE_dvars[3*i  ] = F2[2]*scale;
	}
}

///////////////////////////////////////////////////////////////////////////////
void cartesian_collect_torsional_deriv(
	pose::Pose & pose,
	CartesianMinimizerMap & min_map,
	core::scoring::ScoreFunction const & scorefxn,
	Multivec & dE_dvars,
	core::Real scale
)
{
	using namespace core;
	using namespace core::optimization;
	using namespace id;

	// loop over the torsions
	Size ntorsions=min_map.ntorsions();
	for ( Size i=1; i<=ntorsions; ++i ) {
		id::DOF_ID const & dof_id( min_map.get_dof_id(i) );
		id::TorsionID const & TorsionID( min_map.get_TorsionID(i) );

		// eg rama,Paa,dunbrack,and torsional constraints
		Real deriv = scorefxn.eval_dof_derivative( dof_id, TorsionID, pose );

		if ( deriv == 0 ) continue;

		// work out which atoms are involved !
		AtomID id1, id2, id3, id4;
		pose.conformation().get_torsion_angle_atom_ids( TorsionID, id1,id2,id3,id4 );

		VectorQuad coords( pose.xyz(id1) , pose.xyz(id2) , pose.xyz(id3) , pose.xyz(id4) );
		VectorQuad grads;

		tors_deriv_to_cartesian( deriv, coords, grads );

		core::Size atmidx1 = min_map.get_atom_index(id1);
		core::Size atmidx2 = min_map.get_atom_index(id2);
		core::Size atmidx3 = min_map.get_atom_index(id3);
		core::Size atmidx4 = min_map.get_atom_index(id4);

		if ( atmidx1>0 ) {
			dE_dvars[3*atmidx1-2] += scale*(grads.get<0>().x());
			dE_dvars[3*atmidx1-1] += scale*(grads.get<0>().y());
			dE_dvars[3*atmidx1  ] += scale*(grads.get<0>().z());
		}
		if ( atmidx2>0 ) {
			dE_dvars[3*atmidx2-2] += scale*(grads.get<1>().x());
			dE_dvars[3*atmidx2-1] += scale*(grads.get<1>().y());
			dE_dvars[3*atmidx2  ] += scale*(grads.get<1>().z());
		}
		if ( atmidx3>0 ) {
			dE_dvars[3*atmidx3-2] += scale*(grads.get<2>().x());
			dE_dvars[3*atmidx3-1] += scale*(grads.get<2>().y());
			dE_dvars[3*atmidx3  ] += scale*(grads.get<2>().z());
		}
		if ( atmidx4>0 ) {
			dE_dvars[3*atmidx4-2] += scale*(grads.get<3>().x());
			dE_dvars[3*atmidx4-1] += scale*(grads.get<3>().y());
			dE_dvars[3*atmidx4  ] += scale*(grads.get<3>().z());
		}

	} // loop over torsions
}


///////////////////////////////////////////////////////////////////////////////

void
cart_numerical_derivative_check(
	CartesianMinimizerMap const & min_map,
	CartesianMultifunc const & func,
	Multivec const & start_vars,
	Multivec const & dE_dvars,
	NumericalDerivCheckResultOP deriv_check_result,
	bool const verbose // = true
)
{
	/////////////////////////////////////////////////////////////////////////////
	// NUMERICAL DERIVATIVE CHECK
	/////////////////////////////////////////////////////////////////////////////

	Size const ndofs( min_map.ndofs() );

	bool const write_to_stdout( ! deriv_check_result || deriv_check_result->send_to_stdout() );

	Real const increment = 0.0005; //fpd
	Size const n_increment = 1;
	utility::vector1< Multivec > dE_dvars_numeric( n_increment );
	for ( Size i=1; i<= n_increment; ++i ) {
		dE_dvars_numeric[i].resize( ndofs, 0.0 );
	}

	// setup for saving diagnostics
	NumDerivCheckDataOP min_debug;
	if ( deriv_check_result ) min_debug = NumDerivCheckDataOP( new NumDerivCheckData( ndofs, n_increment ) );

	Multivec vars( start_vars );

	//Real const f00 = func( vars );

	Size natoms=min_map.natoms();
	for ( Size i=1; i<=natoms; ++i ) {
		id::AtomID const & atm_id( min_map.get_atom(i) );
		for ( int jj=0; jj<3; ++jj ) {
			Size ii = (i-1)*3+jj+1;

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
					if ( verbose && write_to_stdout ) {
						// if you change this output, please also change the comments
						// at the beginning of this section
						static bool ratio_header_output( false );
						if ( !ratio_header_output ) {
							ratio_header_output = true;
							TR << "ratio" <<
								A( 4, "inc" ) <<
								A( 4, "rsd" ) <<
								A( 4, "atm" ) <<
								A( 4, "axs" ) <<
								A( 5, "natm" ) <<
								A( 10, "numeric" ) <<
								A( 10, "analytic" ) <<
								A( 10, "ratio" ) <<
								A( 10, "vars[ii]" ) << std::endl;
						}


						TR << "ratio" <<
							I( 4, j ) <<
							I( 4, atm_id.rsd() ) <<
							I( 4, atm_id.atomno() ) <<
							A( 4, (jj==0)?"x":((jj==1)?"y":"z") ) <<
							I( 5, natoms ) <<
							F( 10, 4, deriv ) <<                // column 11
							F( 10, 4, dE_dvars[ii] ) <<         // column 12
							F( 10, 4, ratio ) <<
							F( 10, 4, start_vars[ii] ) << std::endl;
					}

				}
			}
			if ( true ) {
				Real const ratio( std::abs( dE_dvars[ii] ) < 0.001 ? 0.0 :
					deriv_dev / std::abs( dE_dvars[ii] ) );

				if ( min_debug ) min_debug->rel_deriv_dev( ii, ratio );
				if ( min_debug ) min_debug->abs_deriv_dev( ii, deriv_dev );

				//if ( verbose && write_to_stdout ) {
				// if you change this output, please also change the comments
				// at the beginning of this section
				// TR << "deriv_dev:" << SS(ii) << SS(ndofs) << SS(f00) << SS( atm_id.atomno() ) <<
				//  SS( atm_id.rsd() ) << SS( dE_dvars[ii] ) << SS( deriv_dev ) << SS(ratio) << std::endl;
				//}
			}
		}
	}

	// calculate magnitudes, dot products of gradient vectors
	utility::vector1< Real > norm_numeric(n_increment,0.0), dot(n_increment,0.0);
	Real norm(0.0);
	for ( Size i=1; i<= ndofs; ++i ) {
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


// fpd convert torsional derivs to cartesian
void
tors_deriv_to_cartesian( Real dE_dtor, VectorQuad const & coords, VectorQuad & dE_dxs)
{
	// convert to cartesian derivatives on those atoms using the usual MD code - lift from PD
	//core::Real epot_dihedral_this=0;
	core::Real sin_phi, cos_phi;
	core::Vector vti_vta, vta_vtb, vtb_vtj;
	core::Vector nrml1, nrml2, nrml3;
	core::Real inv_nrml1_mag, inv_nrml2_mag, inv_nrml3_mag;

	core::Vector dcosdnrml1(0,0,0), dcosdnrml2(0,0,0), dsindnrml3(0,0,0), dsindnrml2(0,0,0);
	core::Vector f(0,0,0), fi(0,0,0), fab(0,0,0), fj(0,0,0);

	vti_vta = coords.get<0>() - coords.get<1>();
	vta_vtb = coords.get<1>() - coords.get<2>();
	vtb_vtj = coords.get<2>() - coords.get<3>();

	nrml1.x() = (vti_vta.y() * vta_vtb.z() - vti_vta.z() * vta_vtb.y());
	nrml1.y() = (vti_vta.z() * vta_vtb.x() - vti_vta.x() * vta_vtb.z());
	nrml1.z() = (vti_vta.x() * vta_vtb.y() - vti_vta.y() * vta_vtb.x());

	nrml2.x() = (vta_vtb.y() * vtb_vtj.z() - vta_vtb.z() * vtb_vtj.y());
	nrml2.y() = (vta_vtb.z() * vtb_vtj.x() - vta_vtb.x() * vtb_vtj.z());
	nrml2.z() = (vta_vtb.x() * vtb_vtj.y() - vta_vtb.y() * vtb_vtj.x());

	nrml3.x() = (vta_vtb.y() * nrml1.z() - vta_vtb.z() * nrml1.y());
	nrml3.y() = (vta_vtb.z() * nrml1.x() - vta_vtb.x() * nrml1.z());
	nrml3.z() = (vta_vtb.x() * nrml1.y() - vta_vtb.y() * nrml1.x());

	inv_nrml1_mag = 1.0 / nrml1.length();
	inv_nrml2_mag = 1.0 / nrml2.length();
	inv_nrml3_mag = 1.0 / nrml3.length();

	cos_phi = (nrml1.x() * nrml2.x() + nrml1.y() * nrml2.y() + nrml1.z() * nrml2.z()) * inv_nrml1_mag * inv_nrml2_mag;
	sin_phi = (nrml3.x() * nrml2.x() + nrml3.y() * nrml2.y() + nrml3.z() * nrml2.z()) * inv_nrml3_mag * inv_nrml2_mag;

	nrml2.x() *= inv_nrml2_mag;
	nrml2.y() *= inv_nrml2_mag;
	nrml2.z() *= inv_nrml2_mag;

	//phi = -atan2(sin_phi, cos_phi);

	if ( fabs(sin_phi) > 0.1 ) {
		nrml1.x() *= inv_nrml1_mag;
		nrml1.y() *= inv_nrml1_mag;
		nrml1.z() *= inv_nrml1_mag;

		dcosdnrml1.x() = inv_nrml1_mag * (nrml1.x() * cos_phi - nrml2.x());
		dcosdnrml1.y() = inv_nrml1_mag * (nrml1.y() * cos_phi - nrml2.y());
		dcosdnrml1.z() = inv_nrml1_mag * (nrml1.z() * cos_phi - nrml2.z());

		dcosdnrml2.x() = inv_nrml2_mag * (nrml2.x() * cos_phi - nrml1.x());
		dcosdnrml2.y() = inv_nrml2_mag * (nrml2.y() * cos_phi - nrml1.y());
		dcosdnrml2.z() = inv_nrml2_mag * (nrml2.z() * cos_phi - nrml1.z());

	} else {
		nrml3.x() *= inv_nrml3_mag;
		nrml3.y() *= inv_nrml3_mag;
		nrml3.z() *= inv_nrml3_mag;

		dsindnrml3.x() = inv_nrml3_mag * (nrml3.x() * sin_phi - nrml2.x());
		dsindnrml3.y() = inv_nrml3_mag * (nrml3.y() * sin_phi - nrml2.y());
		dsindnrml3.z() = inv_nrml3_mag * (nrml3.z() * sin_phi - nrml2.z());

		dsindnrml2.x() = inv_nrml2_mag * (nrml2.x() * sin_phi - nrml3.x());
		dsindnrml2.y() = inv_nrml2_mag * (nrml2.y() * sin_phi - nrml3.y());
		dsindnrml2.z() = inv_nrml2_mag * (nrml2.z() * sin_phi - nrml3.z());
	}
	dE_dtor *= -1;

	// forces
	if ( fabs(sin_phi) > 0.1 ) {
		dE_dtor /= sin_phi;
		fi.x() += dE_dtor * (vta_vtb.y() * dcosdnrml1.z() - vta_vtb.z() * dcosdnrml1.y());
		fi.y() += dE_dtor * (vta_vtb.z() * dcosdnrml1.x() - vta_vtb.x() * dcosdnrml1.z());
		fi.z() += dE_dtor * (vta_vtb.x() * dcosdnrml1.y() - vta_vtb.y() * dcosdnrml1.x());

		fj.x() += dE_dtor * (vta_vtb.z() * dcosdnrml2.y() - vta_vtb.y() * dcosdnrml2.z());
		fj.y() += dE_dtor * (vta_vtb.x() * dcosdnrml2.z() - vta_vtb.z() * dcosdnrml2.x());
		fj.z() += dE_dtor * (vta_vtb.y() * dcosdnrml2.x() - vta_vtb.x() * dcosdnrml2.y());

		fab.x() += dE_dtor * (vti_vta.z() * dcosdnrml1.y() - vti_vta.y() * dcosdnrml1.z() + vtb_vtj.y() * dcosdnrml2.z() - vtb_vtj.z() * dcosdnrml2.y());
		fab.y() += dE_dtor * (vti_vta.x() * dcosdnrml1.z() - vti_vta.z() * dcosdnrml1.x() + vtb_vtj.z() * dcosdnrml2.x() - vtb_vtj.x() * dcosdnrml2.z());
		fab.z() += dE_dtor * (vti_vta.y() * dcosdnrml1.x() - vti_vta.x() * dcosdnrml1.y() + vtb_vtj.x() * dcosdnrml2.y() - vtb_vtj.y() * dcosdnrml2.x());

	} else {
		dE_dtor /= -cos_phi;

		fi.x() += dE_dtor *
			((vta_vtb.y()*vta_vtb.y() + vta_vtb.z()*vta_vtb.z())*dsindnrml3.x()
			- vta_vtb.x()*vta_vtb.y()*dsindnrml3.y() - vta_vtb.x()* vta_vtb.z()*dsindnrml3.z());
		fi.y() += dE_dtor *
			((vta_vtb.z()*vta_vtb.z() + vta_vtb.x()*vta_vtb.x())*dsindnrml3.y()
			- vta_vtb.y()*vta_vtb.z()*dsindnrml3.z() - vta_vtb.y()* vta_vtb.x()*dsindnrml3.x());
		fi.z() += dE_dtor *
			((vta_vtb.x()*vta_vtb.x() + vta_vtb.y()*vta_vtb.y())*dsindnrml3.z()
			- vta_vtb.z()*vta_vtb.x()*dsindnrml3.x() - vta_vtb.z()*vta_vtb.y()*dsindnrml3.y());

		fj.x() += dE_dtor * (dsindnrml2.y() * vta_vtb.z() - dsindnrml2.z() * vta_vtb.y());
		fj.y() += dE_dtor * (dsindnrml2.z() * vta_vtb.x() - dsindnrml2.x() * vta_vtb.z());
		fj.z() += dE_dtor * (dsindnrml2.x() * vta_vtb.y() - dsindnrml2.y() * vta_vtb.x());

		fab.x() += dE_dtor * (-(vta_vtb.y() * vti_vta.y() + vta_vtb.z() * vti_vta.z()) * dsindnrml3.x()
			+ (2.0 * vta_vtb.x() * vti_vta.y() - vti_vta.x() * vta_vtb.y()) * dsindnrml3.y()
			+ (2.0 * vta_vtb.x() * vti_vta.z() - vti_vta.x() * vta_vtb.z()) * dsindnrml3.z() + dsindnrml2.z() * vtb_vtj.y() - dsindnrml2.y() * vtb_vtj.z());
		fab.y() += dE_dtor * (-(vta_vtb.z() * vti_vta.z() + vta_vtb.x() * vti_vta.x()) * dsindnrml3.y()
			+ (2.0 * vta_vtb.y() * vti_vta.z() - vti_vta.y() * vta_vtb.z()) * dsindnrml3.z()
			+ (2.0 * vta_vtb.y() * vti_vta.x() - vti_vta.y() * vta_vtb.x()) * dsindnrml3.x() + dsindnrml2.x() * vtb_vtj.z() - dsindnrml2.z() * vtb_vtj.x());
		fab.z() += dE_dtor * (-(vta_vtb.x() * vti_vta.x() + vta_vtb.y() * vti_vta.y()) * dsindnrml3.z()
			+ (2.0 * vta_vtb.z() * vti_vta.x() - vti_vta.z() * vta_vtb.x()) * dsindnrml3.x()
			+ (2.0 * vta_vtb.z() * vti_vta.y() - vti_vta.z() * vta_vtb.y()) * dsindnrml3.y() + dsindnrml2.y() * vtb_vtj.x() - dsindnrml2.x() * vtb_vtj.y());
	}

	boost::get<0>(dE_dxs) = fi;
	boost::get<1>(dE_dxs) = fab-fi;
	boost::get<2>(dE_dxs) = fj-fab;
	boost::get<3>(dE_dxs) = -fj;
}


} // namespace optimization
} // namespace core
