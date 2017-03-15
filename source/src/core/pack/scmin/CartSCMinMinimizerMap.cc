// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/CartSCMinMinimizerMap.cc
/// @brief  Class for identifying the sidechain DOFs in the AtomTree which are free during
///         any particular call to the minimizer.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/pack/scmin/CartSCMinMinimizerMap.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/optimization/DOF_Node.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <numeric/constants.hh>

namespace core {
namespace pack {
namespace scmin {

optimization::DOF_NodeOP dummy_nodeop( new optimization::DOF_Node(id::DOF_ID(),NULL) );


CartSCMinMinimizerMap::CartSCMinMinimizerMap() :
	SCMinMinimizerMap(),
	nactive_moving_atoms_total_(0)
{
	// try to avoid reallocations by making this vector big enough for all AAs
	residue_coord_workspace_.reserve(40);
}

CartSCMinMinimizerMap::~CartSCMinMinimizerMap() {}

void CartSCMinMinimizerMap::set_total_residue( Size total_residue )
{
	reset_dof_nodes();
	nactive_residues_ = 0;

	atcs_for_residues_.resize( total_residue );
	active_residues_.resize( total_residue );
	active_residue_index_for_res_.resize( total_residue );
	moving_atoms_.resize( total_residue );
	nactive_moving_atoms_.resize( total_residue );
	atom_derivatives_.resize( total_residue );

	atoms_to_dofid_.resize( total_residue );
	dofid_to_atoms_.resize( total_residue ); // guess ... will grow as needed

	domain_map_.dimension( total_residue );

	std::fill( active_residues_.begin(), active_residues_.end(), 0 );
	std::fill( active_residue_index_for_res_.begin(), active_residue_index_for_res_.end(), 0 );

	for ( Size ii = 1; ii <= domain_map_.size(); ++ii ) domain_map_( ii ) = 1;
}

/// @brief Disable the minimization for all residues.  Ammortized O(1).
void CartSCMinMinimizerMap::clear_active_dofs()
{
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		domain_map_( active_residues_[ ii ] ) = 1;
		atcs_for_residues_[ active_residues_[ ii ]] = 0;
		active_residue_index_for_res_[ active_residues_[ ii ] ] = 0;
		active_residues_[ ii ] = 0;
	}
	nactive_residues_ = 0;
	reset_dof_nodes();
}

/// @details This should be called at most once per residue between calls to "clear_active_chi"
void CartSCMinMinimizerMap::activate_residue_dofs( Size resindex )
{
	debug_assert( domain_map_( resindex ) == 1 ); // activate_residue_chi should not have already been called.
	debug_assert( active_residue_index_for_res_[ resindex ] == 0 ); // activate_residue_chi should not have already been called.

	domain_map_( resindex ) = 0;
	active_residues_[ ++nactive_residues_ ] = resindex;
	active_residue_index_for_res_[ resindex ] = nactive_residues_;

}

/// @brief Convenience lookup -- turns over the request to the AtomTreeCollection
conformation::Residue const &
CartSCMinMinimizerMap::residue( Size seqpos ) const
{
	return atcs_for_residues_[ seqpos ]->active_residue();
}

basic::datacache::BasicDataCache &
CartSCMinMinimizerMap::residue_data( Size seqpos ) const
{
	return atcs_for_residues_[ seqpos ]->active_residue_data();
}


kinematics::tree::Atom const &
CartSCMinMinimizerMap::atom( AtomID const & atid ) const
{
	return atcs_for_residues_[ atid.rsd() ]->active_atom_tree().atom( id::AtomID( atid.atomno(), 1 ) );
}


void CartSCMinMinimizerMap::set_natoms_for_residue( Size resid, Size natoms )
{
	if ( atom_derivatives_[ resid ].size() < natoms ) {
		atom_derivatives_[ resid ].resize( natoms );
	}
}

/// @brief Invoked during the depth-first traversal through the AtomTree.  The AtomTree
/// is indicating that a particular torsion is dependent on another torsion.  Record
/// that fact.
void
CartSCMinMinimizerMap::add_torsion(
	DOF_ID const & /*new_torsion*/,
	DOF_ID const & /*parent*/
)
{
	; // no op
}

/// @brief Invoked during the depth-first traversal through the AtomTree; the atom
/// tree is indicating that a given atom is controlled by a particular DOF.  Record
/// that fact.
void
CartSCMinMinimizerMap::add_atom(
	AtomID const & /*atom_id*/,
	DOF_ID const & /*dof_id*/
)
{
	; // no op
}

/// @brief Traverse the atom trees in preparation for minimization to tie together all the
/// DOFs and the atoms they control.
void
CartSCMinMinimizerMap::setup( AtomTreeCollectionOP trees )
{
	reset_dof_nodes();
	atom_tree_collection_ = trees;

	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		atcs_for_residues_[ iiresid ] = trees->residue_atomtree_collection_op( iiresid );

		focused_residue_ = iiresid;
		conformation::Residue const & iires( atcs_for_residues_[ iiresid ]->active_residue() );


		/////////////////////
		// (1) cartesian dofs ...
		int stop1  = (int)iires.nheavyatoms();
		int start1 = std::min( (int)iires.first_sidechain_atom(), stop1+1 );
		int start2 = (int)iires.first_sidechain_hydrogen();
		int stop2  = (int)iires.natoms();
		int natoms_to_add = (stop1-start1+1)+(stop2-start2+1);

		if ( (int)moving_atoms_[iiresid].size() < natoms_to_add ) {
			moving_atoms_[iiresid].resize(natoms_to_add, id::AtomID( 0,0 ));
		}
		if ( (int)dofid_to_atoms_.size() < (int)(nactive_moving_atoms_total_+natoms_to_add) ) {
			dofid_to_atoms_.resize( nactive_moving_atoms_total_+natoms_to_add );
		}
		if ( atoms_to_dofid_[iiresid].size() < iires.natoms() ) {
			atoms_to_dofid_[iiresid].resize( iires.natoms() );
		}

		for ( int jj=start1; jj<=stop1; ++jj ) {
			++nactive_moving_atoms_[iiresid];
			++nactive_moving_atoms_total_;
			moving_atoms_[iiresid][nactive_moving_atoms_[iiresid]] = id::AtomID( jj,iiresid );
			atoms_to_dofid_[iiresid][jj] = nactive_moving_atoms_total_;
			dofid_to_atoms_[nactive_moving_atoms_total_] = id::AtomID( jj,iiresid );
		}
		for ( int jj=start2; jj<=stop2; ++jj ) {
			++nactive_moving_atoms_[iiresid];
			++nactive_moving_atoms_total_;
			moving_atoms_[iiresid][nactive_moving_atoms_[iiresid]] = id::AtomID( jj,iiresid );
			atoms_to_dofid_[iiresid][jj] = nactive_moving_atoms_total_;
			dofid_to_atoms_[nactive_moving_atoms_total_] = id::AtomID( jj,iiresid );
		}
	}
}

void
CartSCMinMinimizerMap::starting_dofs( optimization::Multivec & dof ) const
{
	dof.resize(3*nactive_moving_atoms_total_);

	Size ctr = 0;
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		conformation::Residue const & iires( atcs_for_residues_[ iiresid ]->active_residue() );
		core::Size nmoving_atoms_ii = moving_atoms_[iiresid].size();
		for ( Size jj = 1; jj <= nmoving_atoms_ii; ++jj ) {
			Vector const &jjxyz = iires.xyz(moving_atoms_[iiresid][jj].atomno());
			dof[++ctr] = jjxyz[0];
			dof[++ctr] = jjxyz[1];
			dof[++ctr] = jjxyz[2];
		}
	}

	debug_assert ( ctr == dof.size() );
}

void
CartSCMinMinimizerMap::assign_dofs_to_mobile_residues( optimization::Multivec const & dofs )
{
	debug_assert( dofs.size() == 3*nactive_moving_atoms_total_ );

	Size ctr = 0;
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		//conformation::Residue const & iires( atcs_for_residues_[ iiresid ]->active_residue() );
		core::Size nmoving_atoms_ii = moving_atoms_[iiresid].size();
		residue_coord_workspace_.resize(nmoving_atoms_ii);   // potentially a reallocation
		for ( Size jj = 1; jj <= nmoving_atoms_ii; ++jj ) {
			residue_coord_workspace_[jj][0] = dofs[++ctr];
			residue_coord_workspace_[jj][1] = dofs[++ctr];
			residue_coord_workspace_[jj][2] = dofs[++ctr];
		}
		atcs_for_residues_[ iiresid ]->set_rescoords( moving_atoms_[iiresid] , residue_coord_workspace_ );
	}

	//fpd -- set_rescoords() calls update_residue()
	//for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
	// atcs_for_residues_[ active_residues_[ ii ] ]->update_residue();
	//}
}

optimization::DOF_Node &
CartSCMinMinimizerMap::dof_node( Size /*index*/ ) {
	return (*dummy_nodeop);  // dummy
}

CartSCMinMinimizerMap::DOF_Node const &
CartSCMinMinimizerMap::dof_node_for_chi( Size /*resid*/, Size /*chiid*/ ) const
{
	return (*dummy_nodeop);  // dummy
}

id::TorsionID
CartSCMinMinimizerMap::tor_for_dof( DOF_ID const & dofid ) const
{
	if ( dofid.type() != core::id::PHI ) return core::id::BOGUS_TORSION_ID;

	Size const rsd( dofid.rsd() );
	Size const chi( residue( rsd ).type().last_controlling_chi( dofid.atomno() ) );
	id::TorsionID torid( rsd, id::CHI, chi );
	return torid;
}

void CartSCMinMinimizerMap::link_torsion_vectors()
{
	; // no op
}

void CartSCMinMinimizerMap::zero_atom_derivative_vectors()
{
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		for ( Size jj = 1, jjend = atom_derivatives_[ iiresid ].size(); jj <= jjend; ++jj ) {
			atom_derivatives_[ iiresid ][ jj ].f1() = 0.0;
			atom_derivatives_[ iiresid ][ jj ].f2() = 0.0;
		}
	}
}

void CartSCMinMinimizerMap::reset_dof_nodes()
{
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		nactive_moving_atoms_[iiresid] = 0;
	}
	nactive_moving_atoms_total_ = 0;
}


} // namespace scmin
} // namespace pack
} // namespace core

