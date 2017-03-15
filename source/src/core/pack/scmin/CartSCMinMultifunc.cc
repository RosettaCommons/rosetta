// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/CartSCMinMultifunc.cc
/// @brief  Class for interfacing with the minimizer during sidechain minimization.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/scmin/CartSCMinMultifunc.hh>

// Package Headers
#include <core/pack/scmin/CartSCMinMinimizerMap.hh>
#include <core/pack/scmin/SCMinMinimizerMap.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/cartesian_minimize.hh> // VectorQuad, tors->cart conversions

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/constants.hh>

#include <core/kinematics/Jump.hh>
#include <core/optimization/DOF_Node.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/chemical/ResidueType.hh>


namespace core {
namespace pack {
namespace scmin {


CartSCMinMultifunc::CartSCMinMultifunc(
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
	scminmap_(dynamic_cast<CartSCMinMinimizerMap&>( scminmap )),
	scoretypes_( sfxn_.get_nonzero_weighted_scoretypes() )
{
}

CartSCMinMultifunc::~CartSCMinMultifunc() {}


//fpd  same as AtomTree version
//fpd    - new logic is in assign_dofs_to_mobile_residues
Real
CartSCMinMultifunc::operator ()( Multivec const & xs ) const
{
	using namespace scoring;

	scminmap_.assign_dofs_to_mobile_residues( xs );

	/// 1. setup for scoring -- nodes first
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		g_.get_minimization_node( iiresid )->setup_for_scoring( scminmap_.residue( iiresid ), scminmap_.residue_data( iiresid ), pose_, sfxn_ );
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


//fpd only logic change in stages 3 & 4
void
CartSCMinMultifunc::dfunc( Multivec const & chi, Multivec & dE_dx ) const
{
	using namespace scoring;
	scminmap_.assign_dofs_to_mobile_residues( chi );
	scminmap_.zero_atom_derivative_vectors();

	/// 1. setup for derivatives -- nodes first
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		g_.get_minimization_node( iiresid )->setup_for_derivatives( scminmap_.residue( iiresid ), scminmap_.residue_data( iiresid ), pose_, sfxn_ );
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
		id::AtomID const & atom_id( scminmap_.get_atom(ii) );
		Vector F2 = scminmap_.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f2();
		dE_dx[3*ii-2] = F2[0];
		dE_dx[3*ii-1] = F2[1];
		dE_dx[3*ii  ] = F2[2];
	}

	/// 4. get torsional derivs & convert to cartesian
	id::DOF_ID junk;
	for ( Size ii = 1, iiend = scminmap_.nactive_residues(); ii <= iiend; ++ii ) {
		Size const iiresid = scminmap_.active_residue( ii );
		conformation::Residue const & iires( scminmap_.residue( iiresid ) );

		// loop through chis -- some of this logic might be better moved to scminminimizermap
		for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
			id::TorsionID torid(iiresid, id::CHI, jj);
			chemical::AtomIndices const &chi_atoms = iires.chi_atoms( jj );

			Real dofderiv = eval_dof_deriv_for_minnode(
				*g_.get_minimization_node( iiresid ),
				scminmap_.residue( iiresid ),
				pose_, junk, torid, sfxn_, sfxn_.weights() );

			optimization::VectorQuad coords(
				iires.xyz(chi_atoms[1]),
				iires.xyz(chi_atoms[2]),
				iires.xyz(chi_atoms[3]),
				iires.xyz(chi_atoms[4]) );
			optimization::VectorQuad grads;
			optimization::tors_deriv_to_cartesian( dofderiv, coords, grads );

			Size kkatmidx1 = scminmap_.get_atom_index(id::AtomID(chi_atoms[1], iiresid)); // returns 0 if this atom is fixed
			Size kkatmidx2 = scminmap_.get_atom_index(id::AtomID(chi_atoms[2], iiresid));
			Size kkatmidx3 = scminmap_.get_atom_index(id::AtomID(chi_atoms[3], iiresid));
			Size kkatmidx4 = scminmap_.get_atom_index(id::AtomID(chi_atoms[4], iiresid));
			if ( kkatmidx1>0 ) {
				dE_dx[3*kkatmidx1-2] += grads.get<0>().x();
				dE_dx[3*kkatmidx1-1] += grads.get<0>().y();
				dE_dx[3*kkatmidx1  ] += grads.get<0>().z();
			}
			if ( kkatmidx2>0 ) {
				dE_dx[3*kkatmidx2-2] += grads.get<1>().x();
				dE_dx[3*kkatmidx2-1] += grads.get<1>().y();
				dE_dx[3*kkatmidx2  ] += grads.get<1>().z();
			}
			if ( kkatmidx3>0 ) {
				dE_dx[3*kkatmidx3-2] += grads.get<2>().x();
				dE_dx[3*kkatmidx3-1] += grads.get<2>().y();
				dE_dx[3*kkatmidx3  ] += grads.get<2>().z();
			}
			if ( kkatmidx4>0 ) {
				dE_dx[3*kkatmidx4-2] += grads.get<3>().x();
				dE_dx[3*kkatmidx4-1] += grads.get<3>().y();
				dE_dx[3*kkatmidx4  ] += grads.get<3>().z();
			}
		}

		// TODO if bb is allowed to move, add in contributions from bb torsions
	}
}

bool
CartSCMinMultifunc::abort_min( Multivec const & ) const {
	return false;
}


void
CartSCMinMultifunc::dump( Multivec const & /*vars*/, Multivec const & /*vars2*/ ) const
{}


} // namespace scmin
} // namespace pack
} // namespace core
