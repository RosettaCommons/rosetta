// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/AtomTreeSCMinMinimizerMap.cc
/// @brief  Class for identifying the sidechain DOFs in the AtomTree which are free during
///         any particular call to the minimizer.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/pack/scmin/AtomTreeSCMinMinimizerMap.hh>

// Package Headers
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSet.fwd.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/scmin/AtomTreeCollection.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/optimization/DOF_Node.hh>

//#include <core/pose/Pose.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <numeric/constants.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

namespace core {
namespace pack {
namespace scmin {

AtomTreeSCMinMinimizerMap::AtomTreeSCMinMinimizerMap() :
	SCMinMinimizerMap(),
	dof_mask_( 1, false ),
  n_active_dof_nodes_( 0 ),
  dof_start_for_focused_residue_( 0 ),
  ndofs_added_for_focused_residue_( 0 )
{

}

AtomTreeSCMinMinimizerMap::~AtomTreeSCMinMinimizerMap() {}

void AtomTreeSCMinMinimizerMap::set_total_residue( Size total_residue )
{
	reset_dof_nodes();
	nactive_residues_ = 0;
	active_residues_.resize( total_residue );
	active_residue_index_for_res_.resize( total_residue );
	chi_start_for_active_residue_.resize( total_residue );
	active_residue_atom_to_dofnode_index_.resize( total_residue );
	atom_derivatives_.resize( total_residue );

	domain_map_.dimension( total_residue );
	atcs_for_residues_.resize( total_residue );
	if ( dof_nodes_.size() < total_residue ) dof_nodes_.reserve( total_residue ); // just a guess
	std::fill( active_residues_.begin(), active_residues_.end(), 0 );
	std::fill( active_residue_index_for_res_.begin(), active_residue_index_for_res_.end(), 0 );
	std::fill( chi_start_for_active_residue_.begin(), chi_start_for_active_residue_.end(), 0 );
	//std::fill( domain_map_.begin(), domain_map_.end(), 1 );
	for ( Size ii = 1; ii <= domain_map_.size(); ++ii ) domain_map_( ii ) = 1;
}

/// @brief Disable the minimization for all residues.  Ammortized O(1).
void AtomTreeSCMinMinimizerMap::clear_active_dofs()
{
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		domain_map_( active_residues_[ ii ] ) = 1;
		atcs_for_residues_[ active_residues_[ ii ]] = 0;
		active_residue_index_for_res_[ active_residues_[ ii ] ] = 0;

		active_residues_[ ii ] = 0;
		std::fill( active_residue_atom_to_dofnode_index_[ ii ].begin(), active_residue_atom_to_dofnode_index_[ ii ].end(), 0 );
		chi_start_for_active_residue_[ ii ] = 0;
	}
	nactive_residues_ = 0;
	reset_dof_nodes();
}

/// @details This should be called at most once per residue between calls to "clear_active_chi"
void AtomTreeSCMinMinimizerMap::activate_residue_dofs( Size resindex )
{
	assert( domain_map_( resindex ) == 1 ); // activate_residue_chi should not have already been called.
	assert( active_residue_index_for_res_[ resindex ] == 0 ); // activate_residue_chi should not have already been called.

	domain_map_( resindex ) = 0;
	active_residues_[ ++nactive_residues_ ] = resindex;
	active_residue_index_for_res_[ resindex ] = nactive_residues_;
}

/// @brief Convenience lookup -- turns over the request to the AtomTreeCollection
conformation::Residue const &
AtomTreeSCMinMinimizerMap::residue( Size seqpos ) const
{
	return atcs_for_residues_[ seqpos ]->active_residue();
}

kinematics::tree::Atom const &
AtomTreeSCMinMinimizerMap::atom( AtomID const & atid ) const
{
	return atcs_for_residues_[ atid.rsd() ]->active_atom_tree().atom( id::AtomID( atid.atomno(), 1 ) );
}


void AtomTreeSCMinMinimizerMap::set_natoms_for_residue( Size resid, Size natoms )
{
	if ( atom_derivatives_[ resid ].size() < natoms ) {
		atom_derivatives_[ resid ].resize( natoms );
	}
}

/// @brief Invoked during the depth-first traversal through the AtomTree.  The AtomTree
/// is indicating that a particular torsion is dependent on another torsion.  Record
/// that fact.
void
AtomTreeSCMinMinimizerMap::add_torsion(
	DOF_ID const & new_torsion,
	DOF_ID const & parent
)
{
	DOF_ID new_torsion_corrected( id::AtomID( new_torsion.atomno(), focused_residue_ ), new_torsion.type() );

	//std::cerr << "add_torsion   " << new_torsion << " " << parent << std::endl;
	DOF_NodeOP pnode( 0 );
	if ( parent.valid() ) {
		Size dofind=0;
		if (parent.type() == core::id::PHI) {
			dofind = atoms_representing_chis_[ parent.atomno() ];
		} else if (parent.type() == core::id::D) {
			dofind = atoms_representing_ds_[ parent.atomno() ];
		} else if (parent.type() == core::id::THETA) {
			dofind = atoms_representing_thetas_[ parent.atomno() ];
		}
		assert( dofind != 0 );
		Size const pind = dof_start_for_focused_residue_ + dofind - 1;
		assert( pind != 0 && pind <= n_active_dof_nodes_ );
		pnode = dof_nodes_[ pind ];
	}

	++ndofs_added_for_focused_residue_;
	if (new_torsion.type() == core::id::PHI) {
		assert( atoms_representing_chis_[ new_torsion.atomno() ] == 0 );
		atoms_representing_chis_[ new_torsion.atomno() ] = ndofs_added_for_focused_residue_;
	} else if (new_torsion.type() == core::id::D) {
		assert( atoms_representing_ds_[ new_torsion.atomno() ] == 0 );
		atoms_representing_ds_[ new_torsion.atomno() ] = ndofs_added_for_focused_residue_;
	} else if (new_torsion.type() == core::id::THETA) {
		assert( atoms_representing_thetas_[ new_torsion.atomno() ] == 0 );
		atoms_representing_thetas_[ new_torsion.atomno() ] = ndofs_added_for_focused_residue_;
	}

	if ( n_active_dof_nodes_ >= dof_nodes_.size() ) {
		dof_nodes_.push_back( optimization::DOF_NodeOP( new optimization::DOF_Node( new_torsion_corrected, pnode ) ) );
		++n_active_dof_nodes_;
	} else {
		++n_active_dof_nodes_;
		dof_nodes_[ n_active_dof_nodes_ ]->set_id( new_torsion_corrected );
		dof_nodes_[ n_active_dof_nodes_ ]->set_parent( pnode );
	}
}

/// @brief Invoked during the depth-first traversal through the AtomTree; the atom
/// tree is indicating that a given atom is controlled by a particular DOF.  Record
/// that fact.
void
AtomTreeSCMinMinimizerMap::add_atom(
	AtomID const & atom_id,
	DOF_ID const & dof_id
)
{
	//std::cerr << "add_atom   " << atom_id << " " << dof_id << std::endl;

	if ( dof_id.valid() ) {
		Size dofind;

		//		if (nonideal_) {
		if (dof_id.type() == core::id::D ) {
			Size const d = atoms_representing_ds_[ dof_id.atomno() ];
			assert( d != 0 );
			dofind = dof_start_for_focused_residue_ + d - 1;
			assert( dofind <= n_active_dof_nodes_ );
			dof_nodes_[ dofind ]->add_atom( id::AtomID( atom_id.atomno(), focused_residue_ ) );

		} else if (dof_id.type() == core::id::THETA) {
			Size const theta = atoms_representing_thetas_[ dof_id.atomno() ];
			assert( theta != 0 );
			dofind = dof_start_for_focused_residue_ + theta - 1;
			assert( dofind <= n_active_dof_nodes_ );
			dof_nodes_[ dofind ]->add_atom( id::AtomID( atom_id.atomno(), focused_residue_ ) );
		} else {
			Size const chi = atoms_representing_chis_[ dof_id.atomno() ];
			//assert( chi != 0 || nonideal_);  //fpd  with nonideality we may have a case where theta or d is movable but chi is not
			if (chi != 0) {
				dofind = dof_start_for_focused_residue_ + chi - 1;
				assert( dofind <= n_active_dof_nodes_ );
				/// add the atom, but give it the correct residue id, since the input atom_id will be 1.
				dof_nodes_[ dofind ]->add_atom( id::AtomID( atom_id.atomno(), focused_residue_ ) );
			}
		}
	}
}

/// @brief Traverse the atom trees in preparation for minimization to tie together all the
/// DOFs and the atoms they control.
void
AtomTreeSCMinMinimizerMap::setup( AtomTreeCollectionOP trees )
{
	reset_dof_nodes();
	atom_tree_collection_ = trees;
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		//assert( atcs_for_residues_[ iiresid ] == 0 );
		atcs_for_residues_[ iiresid ] = trees->residue_atomtree_collection_op( iiresid );
		active_residue_atom_to_dofnode_index_[ ii ].resize( atcs_for_residues_[ iiresid ]->active_residue().natoms(), 0 );

		focused_residue_ = iiresid;

		conformation::Residue const & iires( atcs_for_residues_[ iiresid ]->active_residue() );

		/// mark the DOFs in the DOF_ID_Mask for the chi in this residue as being free
		assert( dof_mask_.n_atom( 1 ) == ((nonideal_) ?
		        atoms_representing_chis_.size()+atoms_representing_ds_.size()+atoms_representing_thetas_.size() : atoms_representing_chis_.size()) );
		if ( dof_mask_.n_atom( 1 ) < ((nonideal_?3:1)*iires.natoms()) ) {
			atoms_representing_chis_.resize( iires.natoms(), 0 );
			if (nonideal_) {
				dof_mask_.resize( 1, 3*iires.natoms(), false );  //fpd at most 3 DOFS per atom
				atoms_representing_ds_.resize( iires.natoms(), 0 );
				atoms_representing_thetas_.resize( iires.natoms(), 0 );
			} else {
				dof_mask_.resize( 1, iires.natoms(), false );
			}
		}
		for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
			dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 4 ], 1 ), id::PHI ) ] = true;
		}
		dof_start_for_focused_residue_ = n_active_dof_nodes_ + 1;
		ndofs_added_for_focused_residue_ = 0;

		// descend through the atom tree for this residue.
		// In this call, add_atom() and add_torsion() wil be invoked on this MinimizerMap.
		DOF_ID invalid;
		atcs_for_residues_[ iiresid ]->active_atom_tree().root()->setup_min_map( invalid, dof_mask_, *this );

		// now mark the DOFs in the DOF_ID_Mask for the chi in this residue as being fixed
		for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
			dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 4 ], 1 ), id::PHI ) ] = false;
		}


		if (nonideal_) {
			for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
				if (jj == 1) {
					dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 3 ], 1 ), id::THETA ) ] = true;
				}
				dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 4 ], 1 ), id::THETA ) ] = true;
			}

			DOF_ID invalid2;
			atcs_for_residues_[ iiresid ]->active_atom_tree().root()->setup_min_map( invalid2, dof_mask_, *this );

			for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
				if (jj == 1) {
					dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 3 ], 1 ), id::THETA ) ] = false;
				}
				dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 4 ], 1 ), id::THETA ) ] = false;
			}

			for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
				if (jj == 1) {
					dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 3 ], 1 ), id::D ) ] = true;
				}
				dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 4 ], 1 ), id::D ) ] = true;
			}

			DOF_ID invalid3;
			atcs_for_residues_[ iiresid ]->active_atom_tree().root()->setup_min_map( invalid3, dof_mask_, *this );

			for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
				if (jj == 1) {
					dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 3 ], 1 ), id::D ) ] = false;
				}
				dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 4 ], 1 ), id::D ) ] = false;
			}
		}


		std::fill( atoms_representing_chis_.begin(), atoms_representing_chis_.end(), 0 ); // overkill
		if (nonideal_) {
			std::fill( atoms_representing_thetas_.begin(), atoms_representing_thetas_.end(), 0 ); //fpd  also overkill (?)
			std::fill( atoms_representing_ds_.begin(), atoms_representing_ds_.end(), 0 ); //fpd  also overkill (?)
		}
	}
}

void
AtomTreeSCMinMinimizerMap::starting_dofs( optimization::Multivec & dof ) const
{
	dof.resize( n_active_dof_nodes_ );

	for ( Size ii = 1; ii <= n_active_dof_nodes_; ++ii ) {
		DOF_Node const & iinode( * dof_nodes_[ ii ] );
		Size const iirsd = iinode.rsd();

		DOF_ID new_torsion_uncorrected( id::AtomID( iinode.atomno(), 1 ), iinode.type() );

		if (iinode.type() == core::id::PHI) {
			Size const iichi = atcs_for_residues_[ iirsd ]->active_restype().last_controlling_chi( iinode.atomno() );
			assert( atcs_for_residues_[ iirsd ]->active_restype().chi_atoms( iichi )[ 4 ] == (Size) iinode.atomno() );
			dof[ ii ] = atcs_for_residues_[ iirsd ]->active_residue().chi( iichi );
		} else if (iinode.type() == core::id::D) {
			dof[ ii ] = atcs_for_residues_[ iirsd ]->dof( new_torsion_uncorrected );
			//fpd   The '0.01' here is odd
			//fpd     the reason is that scmin works in degree-space rather than radian-space (as atom-tree-min does)
			//fpd     we want to keep the scaling the same as in atom-tree min
			dof[ ii ] *= (1./0.01)*basic::options::option[ basic::options::OptionKeys::optimization::scale_d ]();
		} else if (iinode.type() == core::id::THETA) {
			dof[ ii ] = numeric::constants::d::radians_to_degrees * atcs_for_residues_[ iirsd ]->dof( new_torsion_uncorrected );
			dof[ ii ] *= basic::options::option[ basic::options::OptionKeys::optimization::scale_theta ]();
		}
	}

	//fpd  previous version
	//for ( Size ii = 1; ii <= n_active_dof_nodes_; ++ii ) {
	//	DOF_Node const & iinode( * dof_nodes_[ ii ] );
	//	Size const iirsd = iinode.rsd();
	//	Size const iichi = atcs_for_residues_[ iirsd ]->active_restype().last_controlling_chi( iinode.atomno() );
	//	assert( atcs_for_residues_[ iirsd ]->active_restype().chi_atoms( iichi )[ 4 ] == (Size) iinode.atomno() );
	//	dof[ ii ] = atcs_for_residues_[ iirsd ]->active_residue().chi( iichi );
	//}
}

void
AtomTreeSCMinMinimizerMap::assign_dofs_to_mobile_residues( optimization::Multivec const & dofs )
{
	assert( dofs.size() == n_active_dof_nodes_ );

	for ( Size ii = 1; ii <= n_active_dof_nodes_; ++ii ) {
		DOF_Node const & iinode( * dof_nodes_[ ii ] );
		Size const iirsd = iinode.rsd();

		//fpd  for thetas and ds this is a bit of a hack
		//     the only d- and theta- minimizable atom not controlled by a chi is CB
		//     so we'll treat chi==0 in ResidueAtomTreeCollection::set_d/set_theta
		//         as the dist/angle completed by atoms 2 & 3 of chi 1
		Size const iichi = atcs_for_residues_[ iirsd ]->active_restype().last_controlling_chi( iinode.atomno() );
		assert( atcs_for_residues_[ iirsd ]->active_restype().chi_atoms( iichi==0?1:iichi )[ iichi==0?3:4 ] == (Size) iinode.atomno() );

		if (iinode.type() == core::id::PHI) {
			atcs_for_residues_[ iirsd ]->set_chi( iichi, dofs[ ii ] );
		} else if (iinode.type() == core::id::D) {
			//fpd   The '0.01' here is odd
			//fpd     the reason is that scmin works in degree-space rather than radian-space (as atom-tree-min does)
			//fpd     we want to keep the scaling the same as in atom-tree min
			atcs_for_residues_[ iirsd ]->set_d( iichi, 0.01 * dofs[ ii ] / basic::options::option[ basic::options::OptionKeys::optimization::scale_d ]() );
		} else if (iinode.type() == core::id::THETA) {
			atcs_for_residues_[ iirsd ]->set_theta( iichi, dofs[ ii ] / basic::options::option[ basic::options::OptionKeys::optimization::scale_theta ]() );
		}
	}

	//fpd  previous version
	//for ( Size ii = 1; ii <= n_active_dof_nodes_; ++ii ) {
	//	DOF_Node const & iinode( * dof_nodes_[ ii ] );
	//	Size const iirsd = iinode.rsd();
	//	Size const iichi = atcs_for_residues_[ iirsd ]->active_restype().last_controlling_chi( iinode.atomno() );
	//	assert( atcs_for_residues_[ iirsd ]->active_restype().chi_atoms( iichi )[ 4 ] == (Size) iinode.atomno() );
	//	atcs_for_residues_[ iirsd ]->set_chi( iichi, chi[ ii ] );
	//}

	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		atcs_for_residues_[ active_residues_[ ii ] ]->update_residue();
	}
}

AtomTreeSCMinMinimizerMap::DOF_Node const &
AtomTreeSCMinMinimizerMap::dof_node_for_chi( Size resid, Size chiid ) const
{
	Size atno = atcs_for_residues_[ resid ]->active_residue().chi_atoms( chiid )[ 4 ];
	Size actid = active_residue_index_for_res_[ resid ];
	Size dofnode_ind = active_residue_atom_to_dofnode_index_[ actid ][ atno ];
	return * dof_nodes_[ dofnode_ind ];
}

id::TorsionID
AtomTreeSCMinMinimizerMap::tor_for_dof( DOF_ID const & dofid ) const
{
	if (dofid.type() != core::id::PHI) return core::id::BOGUS_TORSION_ID;

	Size const rsd( dofid.rsd() );
	Size const chi( residue( rsd ).type().last_controlling_chi( dofid.atomno() ) );
	id::TorsionID torid( rsd, id::CHI, chi );
	return torid;
}

/// @details super simple -- no need for a sort (nor is there a need in the optimization::MinimizerMap for that matter).
void AtomTreeSCMinMinimizerMap::link_torsion_vectors()
{
	/// by construction, DOFs are added such that index( parent ) < index( child ).
	/// (i.e. we have a partial order of DOFs and DOF_Nodes are added in a monotonically increasing manner in that partial order)
	//fpd  this still holds for Ds/thetas
	for ( Size ii = n_active_dof_nodes_; ii >= 1; --ii ) {
		dof_nodes_[ ii ]->link_vectors();
	}
}

void AtomTreeSCMinMinimizerMap::zero_atom_derivative_vectors()
{
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		for ( Size jj = 1, jjend = atom_derivatives_[ iiresid ].size(); jj <= jjend; ++jj ) {
			atom_derivatives_[ iiresid ][ jj ].f1() = 0.0;
			atom_derivatives_[ iiresid ][ jj ].f2() = 0.0;
		}
	}
}

void AtomTreeSCMinMinimizerMap::reset_dof_nodes()
{
	for ( Size ii = 1; ii <= n_active_dof_nodes_; ++ii ) {
		dof_nodes_[ ii ]->set_parent( 0 );
		dof_nodes_[ ii ]->clear_atoms();
	}
	n_active_dof_nodes_ = 0;
}



} // namespace scmin
} // namespace pack
} // namespace core

