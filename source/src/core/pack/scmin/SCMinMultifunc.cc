// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/SCMinMultifunc.cc
/// @brief  Class for interfacing with the minimizer during sidechain minimization.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/scmin/SCMinMultifunc.hh>

// Package Headers
#include <core/pack/scmin/SCMinMinimizerMap.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/atom_tree_minimize.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/constants.hh>

#include <core/kinematics/Jump.hh>
#include <core/optimization/DOF_Node.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
using namespace ObjexxFCL::format;

#include <basic/Tracer.hh>

namespace core {
namespace pack {
namespace scmin {

static THREAD_LOCAL basic::Tracer TR( "core.pack.scmin.SCMinMultifunc" );


SCMinMultifunc::SCMinMultifunc(
	pose::Pose & p,
	utility::vector1< conformation::ResidueCOP > const & bg_residues,
	scoring::ScoreFunction const & sfxn,
	scoring::MinimizationGraph & mingraph,
	SCMinMinimizerMap & scminmap
) :
	pose_( p ),
	bg_residues_( bg_residues ),
	sfxn_( sfxn ),
	g_( mingraph ),
	scminmap_( scminmap ),
	scoretypes_( sfxn_.get_nonzero_weighted_scoretypes() )
{}

SCMinMultifunc::~SCMinMultifunc() {}

Real
SCMinMultifunc::operator ()( Multivec const & chi ) const
{
	using namespace scoring;

	scminmap_.assign_dofs_to_mobile_residues( chi );

	/// 1. setup for scoring -- nodes first
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		g_.get_minimization_node( iiresid )->setup_for_scoring( scminmap_.residue( iiresid ), pose_, sfxn_ );
	}
	/// edges second.
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		conformation::Residue const & iires( scminmap_.residue( iiresid ) );
		for ( utility::graph::Graph::EdgeListIter
				edgeit = g_.get_node( iiresid )->edge_list_begin(),
				edgeit_end = g_.get_node( iiresid )->edge_list_end();
				edgeit != edgeit_end; ++edgeit ) {
			Size const jjresid = (*edgeit)->get_other_ind( iiresid );
			bool const jjmoving = scminmap_.dm()( jjresid ) == 0;
			if ( jjmoving && iiresid > jjresid ) continue;
			conformation::Residue const & jjres( jjmoving ? scminmap_.residue( jjresid ) : *bg_residues_[ jjresid ] );

			MinimizationEdge & minedge = static_cast< MinimizationEdge & > ( (**edgeit) );
			if ( iiresid < jjresid ) {
				minedge.setup_for_scoring( iires, jjres, pose_, sfxn_ );
			} else {
				minedge.setup_for_scoring( jjres, iires, pose_, sfxn_ );
			}

		}
	}

	//emap.zero( scoretypes_ );
	EnergyMap emap;
	/// 2. score function evaluation -- 1body & 2body energies only
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		eval_res_onebody_energies_for_minnode( * g_.get_minimization_node( iiresid ), scminmap_.residue( iiresid ), pose_, sfxn_, emap );
	}
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		conformation::Residue const & iires( scminmap_.residue( iiresid ) );
		for ( utility::graph::Graph::EdgeListConstIter
				edgeit = g_.get_node( iiresid )->const_edge_list_begin(),
				edgeit_end = g_.get_node( iiresid )->const_edge_list_end();
				edgeit != edgeit_end; ++edgeit ) {
			Size const jjresid = (*edgeit)->get_other_ind( iiresid );
			bool const jjmoving = scminmap_.dm()( jjresid ) == 0;
			if ( jjmoving && iiresid > jjresid ) continue;
			conformation::Residue const & jjres( jjmoving ? scminmap_.residue( jjresid ) : *bg_residues_[ jjresid ] );
			MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );

			if ( iiresid < jjresid ) {
				eval_res_pair_energy_for_minedge( minedge, iires, jjres, pose_, sfxn_, emap );
			} else {
				eval_res_pair_energy_for_minedge( minedge, jjres, iires, pose_, sfxn_, emap );
			}
		}
	}
	return sfxn_.weights().dot( emap, scoretypes_ );
}


void
SCMinMultifunc::dfunc( Multivec const & chi, Multivec & dE_dchi ) const
{
	//using numeric::constants::d::pi;
	using namespace scoring;
	scminmap_.assign_dofs_to_mobile_residues( chi );
	scminmap_.zero_atom_derivative_vectors();

	/// 1. setup for derivatives -- nodes first
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		g_.get_minimization_node( iiresid )->setup_for_derivatives( scminmap_.residue( iiresid ), pose_, sfxn_ );
	}
	/// edges second.
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		conformation::Residue const & iires( scminmap_.residue( iiresid ) );
		for ( utility::graph::Graph::EdgeListIter
				edgeit = g_.get_node( iiresid )->edge_list_begin(),
				edgeit_end = g_.get_node( iiresid )->edge_list_end();
				edgeit != edgeit_end; ++edgeit ) {
			Size const jjresid = (*edgeit)->get_other_ind( iiresid );
			bool const jjmoving = scminmap_.dm()( jjresid ) == 0;
			if ( jjmoving && iiresid > jjresid ) continue;
			conformation::Residue const & jjres( jjmoving ? scminmap_.residue( jjresid ) : *bg_residues_[ jjresid ] );
			MinimizationEdge & minedge = static_cast< MinimizationEdge & > ( (**edgeit) );
			if ( iiresid < jjresid ) {
				minedge.setup_for_derivatives( iires, jjres, pose_, sfxn_ );
			} else {
				minedge.setup_for_derivatives( jjres, iires, pose_, sfxn_ );
			}
		}
	}

	/// 2a. Evaluate 1body atom derivatives for the active nodes
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		conformation::Residue const & iirsd = scminmap_.residue( iiresid );
		MinimizationNode const & minnode( * g_.get_minimization_node( iiresid ) );
		eval_atom_derivatives_for_minnode( minnode, iirsd, pose_,
			sfxn_.weights(), scminmap_.atom_derivatives( iiresid ) );
	}
	/// 2b. Evaluate 2body atom derivatives for the active nodes
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		MinimizationNode const & iiminnode( * g_.get_minimization_node( iiresid ) );
		conformation::Residue const & iires( scminmap_.residue( iiresid ) );
		for ( utility::graph::Graph::EdgeListConstIter
				edgeit = g_.get_node( iiresid )->edge_list_begin(),
				edgeit_end = g_.get_node( iiresid )->edge_list_end();
				edgeit != edgeit_end; ++edgeit ) {
			Size const jjresid = (*edgeit)->get_other_ind( iiresid );
			bool const jjmoving = scminmap_.dm()( jjresid ) == 0;
			if ( jjmoving && iiresid > jjresid ) continue;

			conformation::Residue const & jjres( jjmoving ? scminmap_.residue( jjresid ) : *bg_residues_[ jjresid ] );
			MinimizationNode const & jjminnode( * g_.get_minimization_node( jjresid ) );
			MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
			if ( iiresid < jjresid ) {
				eval_atom_derivatives_for_minedge( minedge, iires, jjres,
					iiminnode.res_min_data(), jjminnode.res_min_data(), pose_, sfxn_.weights(),
					scminmap_.atom_derivatives( iiresid ), scminmap_.atom_derivatives( jjresid ));
			} else {
				eval_atom_derivatives_for_minedge( minedge, jjres, iires,
					jjminnode.res_min_data(), iiminnode.res_min_data(), pose_, sfxn_.weights(),
					scminmap_.atom_derivatives( jjresid ), scminmap_.atom_derivatives( iiresid ));
			}
		}
	}


	/// 3. Accumulate the atom derivatives into the DOF_Node f1/f2s
	for ( Size ii = 1, iiend = scminmap_.n_dof_nodes(); ii <= iiend; ++ii ) {
		optimization::DOF_Node & iidofnode( scminmap_.dof_node( ii ) );
		iidofnode.F1() = 0; iidofnode.F2() = 0;
		optimization::DOF_Node::AtomIDs iiatoms( iidofnode.atoms() );
		for ( Size jj = 1, jjend = iiatoms.size(); jj <= jjend; ++jj ) {
			iidofnode.F1() += scminmap_.atom_derivatives( iiatoms[ jj ].rsd() )[ iiatoms[ jj ].atomno() ].f1();
			iidofnode.F2() += scminmap_.atom_derivatives( iiatoms[ jj ].rsd() )[ iiatoms[ jj ].atomno() ].f2();
		}
	}

	/// 3. Link torsion vectors
	scminmap_.link_torsion_vectors();

	id::DOF_ID junk;
	for ( Size ii = 1, iiend = scminmap_.n_dof_nodes(); ii <= iiend; ++ii ) {
		optimization::DOF_Node & iidofnode( scminmap_.dof_node( ii ) );
		//fpd  we assume that torsion-defined funcs have no derivatives w.r.t. d or theta DOFs
		Real dofderiv = 0.0;
		if ( iidofnode.type() == core::id::PHI ) {
			id::TorsionID torid = scminmap_.tor_for_dof( iidofnode.dof_id() );
			/// 4. Evaluate chi-dof derivatives
			dofderiv = eval_dof_deriv_for_minnode( * g_.get_minimization_node( iidofnode.rsd() ),
				scminmap_.residue( iidofnode.rsd() ), pose_, junk, torid, sfxn_, sfxn_.weights() );
		}

		// scalefactor
		Real factor( 1.0 );
		if ( iidofnode.type() == id::PHI ) {
			factor = numeric::constants::d::rad2deg;
		} else if ( iidofnode.type() == id::THETA ) {
			factor = numeric::constants::d::rad2deg * basic::options::option[ basic::options::OptionKeys::optimization::scale_theta ]();
			Real const theta( chi[ii] );
			int const theta_mod ( ( static_cast< int >( std::floor( theta/180 )))%2);
			if ( theta_mod == 1 || theta_mod == -1 ) factor *= -1.0;
		} else if ( iidofnode.type() == id::D ) {
			//fpd   The '0.01' here is odd
			//fpd     the reason is that scmin works in degree-space rather than radian-space (as atom-tree-min does)
			//fpd     we want to keep the scaling ~ the same as in atom-tree min
			factor = (1./0.01)*basic::options::option[ basic::options::OptionKeys::optimization::scale_d ]();
		}

		dE_dchi[ ii ] = optimization::torsional_derivative_from_cartesian_derivatives(
			scminmap_.atom( iidofnode.atom_id() ), iidofnode,
			dofderiv, factor );
	}

	//fpd uncomment to debug
	//scmin_numerical_derivative_check( chi, dE_dchi );
}

bool
SCMinMultifunc::abort_min( Multivec const & ) const {
	return false;
}


void
SCMinMultifunc::dump( Multivec const & /*vars*/, Multivec const & /*vars2*/ ) const
{}


void
SCMinMultifunc::scmin_numerical_derivative_check( Multivec const & start_vars, Multivec & dE_dvars ) const {
	/////////////////////////////////////////////////////////////////////////////
	// NUMERICAL DERIVATIVE CHECK
	/////////////////////////////////////////////////////////////////////////////

	Size const ndofs( scminmap_.n_dof_nodes() );

	Real const factor = 0.001; //fpd
	Multivec dE_dvars_numeric;
	dE_dvars_numeric.resize( ndofs, 0.0 );

	Multivec vars( start_vars );

	for ( Size ii=1; ii<=ndofs; ++ii ) {
		Real deriv_dev = 10000.0;

		vars[ii] = start_vars[ii] + factor;
		Real const f11 = (*this)( vars );
		vars[ii] = start_vars[ii] - factor;
		Real const f22 = (*this)( vars );
		Real const deriv = ( f11 - f22 ) / ( factor * 2 );

		dE_dvars_numeric[ii] = deriv;
		deriv_dev = std::min( deriv_dev, std::abs( deriv  - dE_dvars[ii] ) );
		vars[ii] = start_vars[ii];

		Real const ratio( std::abs( dE_dvars[ii] ) < 0.001 ? 0.0 : deriv / dE_dvars[ii] );
		if ( std::abs(dE_dvars[ii]) > 0.001 || std::abs(deriv) > 0.001 ) {
			static bool ratio_header_output( false );
			if ( !ratio_header_output ) {
				ratio_header_output = true;
				TR << "ratio" << A( 6, "dofid" ) << A( 10, "doftype" ) <<
					A( 10, "numeric" ) << A( 10, "analytic" ) << A( 10, "ratio" ) << A( 10, "vars[ii]" ) << std::endl;
			}
			TR << "ratio" <<  I( 6, ii ) << I( 10, scminmap_.dof_node(ii).type() ) <<
				F( 10, 4, deriv ) <<  F( 10, 4, dE_dvars[ii] ) <<  F( 10, 4, ratio ) <<
				F( 10, 4, start_vars[ii] ) << std::endl;
		}
	}

	// calculate magnitudes, dot products of gradient vectors
	Real norm_numeric(0.0), dot(0.0);
	Real norm(0.0);
	for ( Size i=1; i<= ndofs; ++i ) {
		norm += dE_dvars[i] * dE_dvars[i];
		dot += dE_dvars[i] * dE_dvars_numeric[i];
		norm_numeric += dE_dvars_numeric[i] * dE_dvars_numeric[i];
	}
	norm = std::sqrt( norm );
	norm_numeric = std::sqrt( norm_numeric );

	// handle strange cases
	Real log_norm_ratio;
	if ( norm < 0.001 && norm_numeric < 0.001 ) {
		log_norm_ratio = 1.0;
	} else if ( norm < 0.001 ) {
		log_norm_ratio = 100.0;
	} else if ( norm_numeric < 0.001 ) {
		log_norm_ratio = -100.0;
	} else {
		log_norm_ratio = std::log( norm_numeric / norm );
	}
	Real const cos_theta( dot / ( norm * norm_numeric) );

	TR << " norm: " << F(12,4,norm) <<
		" norm_numeric: " << F(12,4,norm_numeric) <<
		" cos_theta: " << F(7,4,cos_theta) <<
		" log_norm_ratio: " << F(9,4,log_norm_ratio) << std::endl;
}

} // namespace scmin
} // namespace pack
} // namespace core
