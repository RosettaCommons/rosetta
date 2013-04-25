// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/SCMinMultifunc.cc
/// @brief  Class for interfacing with the minimizer during sidechain minimization.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/scmin/SCMinMultifunc.hh>

// Package Headers
// AUTO-REMOVED #include <core/pack/scmin/AtomTreeCollection.hh>
#include <core/pack/scmin/SCMinMinimizerMap.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
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


namespace core {
namespace pack {
namespace scmin {


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
		for ( graph::Graph::EdgeListIter
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
		for ( graph::Graph::EdgeListConstIter
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
		for ( graph::Graph::EdgeListIter
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
		for ( graph::Graph::EdgeListConstIter
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
		if (iidofnode.type() == core::id::PHI) {
			id::TorsionID torid = scminmap_.tor_for_dof( iidofnode.dof_id() );
			/// 4. Evaluate chi-dof derivatives
			Real dofderiv = eval_dof_deriv_for_minnode( * g_.get_minimization_node( iidofnode.rsd() ),
				scminmap_.residue( iidofnode.rsd() ), pose_, junk, torid, sfxn_, sfxn_.weights() );

			dE_dchi[ ii ] = optimization::torsional_derivative_from_cartesian_derivatives(
				scminmap_.atom( iidofnode.atom_id() ), iidofnode,
				dofderiv, numeric::constants::d::rad2deg );
		}
	}
}

bool
SCMinMultifunc::abort_min( Multivec const & ) const {
	return false;
}


void
SCMinMultifunc::dump( Multivec const & /*vars*/, Multivec const & /*vars2*/ ) const
{}


/*void
SCMinMultifunc::eval_atom_deriv(
	id::AtomID const & atom,
	Vector & F1,
	Vector & F2
) const
{
	using namespace scoring;

	Size const rsdno = atom.rsd();
	Size const atomno = atom.atomno();

	conformation::Residue const & rsd = scminmap_.residue( rsdno );

	MinimizationNode const & minnode =  * g_.get_minimization_node( rsdno );
	/// 1. eval intra-residue derivatives
	eval_atom_derivative_for_minnode( minnode, atomno, rsd, pose_, scminmap_.dm(), sfxn_, sfxn_.weights(), F1, F2 );

	ResSingleMinimizationData const & ressingle_min_data( minnode.res_min_data() );
	/// 2. eval inter-residue derivatives
	for ( graph::Node::EdgeListConstIter
			edgeit = minnode.const_edge_list_begin(), edgeit_end = minnode.const_edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
		Size const other_rsdno = minedge.get_other_ind( rsdno );
		conformation::Residue const & other_rsd( scminmap_.dm()( other_rsdno ) == 0 ? scminmap_.residue( other_rsdno ) : *bg_residues_[ other_rsdno ]);
		ResSingleMinimizationData const & other_ressingle_min_data( g_.get_minimization_node( other_rsdno )->res_min_data() );

		eval_atom_deriv_for_minedge( minedge, atomno, rsd, other_rsd,
			ressingle_min_data, other_ressingle_min_data,
			pose_, scminmap_.dm(), sfxn_, sfxn_.weights(), F1, F2 );
	}
}*/


} // namespace scmin
} // namespace pack
} // namespace core
