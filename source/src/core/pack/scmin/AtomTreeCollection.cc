// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/AtomTreeCollection.cc
/// @brief  Implementation for the classes holding sets of AtomTrees used during variuos packing+minimizing schemes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/scmin/AtomTreeCollection.hh>

// Package Headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/DOF_ID.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/util.hh>
#include <core/kinematics/tree/Atom.hh>
//#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

// Numeric headers
#include <numeric/constants.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace scmin {

ResidueAtomTreeCollectionMomento::ResidueAtomTreeCollectionMomento() :
	restype_index_( 0 ),
	natoms_( 0 )
{}

ResidueAtomTreeCollectionMomento::~ResidueAtomTreeCollectionMomento(){}

ResidueAtomTreeCollectionMomento::ResidueAtomTreeCollectionMomento(
	ResidueAtomTreeCollectionMomento const & src
) :
	utility::pointer::ReferenceCount(),
	restype_index_( src.restype_index_ ),
	natoms_( src.natoms_ ),
	coords_( src.coords_ )
{
}

ResidueAtomTreeCollectionMomento &
ResidueAtomTreeCollectionMomento::operator = ( ResidueAtomTreeCollectionMomento const & rhs )
{
	if ( this != & rhs ) {
		restype_index_ = rhs.restype_index_;
		natoms_ = rhs.natoms_;
		if ( coords_.size() < rhs.coords_.size() ) coords_.resize( rhs.coords_.size() );
		for ( Size ii = 1, iiend = natoms_; ii <= iiend; ++ii ) coords_[ ii ] = rhs.coords_[ ii ];
	}
	return *this;
}

void ResidueAtomTreeCollectionMomento::set_restype_index( Size setting )
{
	restype_index_ = setting;
}

void ResidueAtomTreeCollectionMomento::copy_coords( conformation::Residue const & res )
{
	if ( coords_.size() < res.natoms() ) { coords_.resize( res.natoms() ); }
	natoms_ = res.natoms();
	for ( Size ii = 1, iiend = res.natoms(); ii <= iiend; ++ii ) {
		coords_[ ii ] = res.xyz( ii );
	}
}


ResidueAtomTreeCollection::ResidueAtomTreeCollection(
	task::ResidueLevelTask const & rltask,
	conformation::Conformation const & conformation,
	conformation::Residue const & orig_res
) :
	active_restype_( 0 ),
	residue_uptodate_( true ),
	atom_tree_uptodate_( true ),
	atom_tree_representatives_( rltask.allowed_residue_types().size() ),
	residue_representatives_( atom_tree_representatives_.size() )
{
	//Size const nrestypes( atom_tree_representatives_.size() );
	Size const resid( orig_res.seqpos() );

	Size ii( 0 );
	for ( task::ResidueLevelTask::ResidueTypeCOPListConstIter
			allowed_iter = rltask.allowed_residue_types_begin(),
			allowed_end = rltask.allowed_residue_types_end();
			allowed_iter != allowed_end; ++allowed_iter ) {
		++ii;

		residue_representatives_[ ii ] = conformation::ResidueOP( new conformation::Residue(
			**allowed_iter, orig_res, conformation, rltask.preserve_c_beta() ) );
		residue_representatives_[ ii ]->seqpos( 1 ); // temporary -- while we construct the atom tree, pretend we're residue 1.
		kinematics::AtomPointer2D tree_atoms2d( 1 ); // 1 residue AtomTree.
		tree_atoms2d[ 1 ].resize( residue_representatives_[ ii ]->natoms() );
		kinematics::AtomPointer1D tree_atoms1d( residue_representatives_[ ii ]->natoms() );

		conformation::build_residue_tree( 1, *residue_representatives_[ ii ], tree_atoms1d, true );
		std::copy( tree_atoms1d.begin(), tree_atoms1d.end(), tree_atoms2d[1].begin() );
		atom_tree_representatives_[ ii ] = kinematics::AtomTreeOP( new kinematics::AtomTree( tree_atoms2d ) );

		residue_representatives_[ ii ]->seqpos( resid ); // restore the original sequence position
	}
}

ResidueAtomTreeCollection::ResidueAtomTreeCollection(
	rotamer_set::RotamerSet const & rset,
	core::Size resid
) :
	active_restype_( 0 ),
	residue_uptodate_( true ),
	atom_tree_uptodate_( true ),
	atom_tree_representatives_( rset.get_n_residue_types() ),
	residue_representatives_( atom_tree_representatives_.size() )
{

	for ( Size ii = 1; ii <= rset.get_n_residue_types(); ++ii ) {
		residue_representatives_[ ii ] = rset.rotamer( rset.get_residue_type_begin(ii) )->clone();
		residue_representatives_[ ii ]->seqpos( 1 ); // temporary -- while we construct the atom tree, pretend we're residue 1.
		kinematics::AtomPointer2D tree_atoms2d( 1 ); // 1 residue AtomTree.
		tree_atoms2d[ 1 ].resize( residue_representatives_[ ii ]->natoms() );
		kinematics::AtomPointer1D tree_atoms1d( residue_representatives_[ ii ]->natoms() );

		conformation::build_residue_tree( 1, *residue_representatives_[ ii ], tree_atoms1d, true );
		std::copy( tree_atoms1d.begin(), tree_atoms1d.end(), tree_atoms2d[1].begin() );
		atom_tree_representatives_[ ii ] = kinematics::AtomTreeOP( new kinematics::AtomTree( tree_atoms2d ) );

		residue_representatives_[ ii ]->seqpos( resid ); // restore the original sequence position
	}
}

ResidueAtomTreeCollection::~ResidueAtomTreeCollection()
{}


void ResidueAtomTreeCollection::set_active_restype_index( Size restype_index )
{
	debug_assert( restype_index > 0 && restype_index <= atom_tree_representatives_.size() );
	active_restype_ = restype_index;
	active_atom_tree_ = atom_tree_representatives_[ active_restype_ ];
	active_residue_ = residue_representatives_[ active_restype_ ];
}


chemical::ResidueType const &
ResidueAtomTreeCollection::active_restype() const
{
	return active_residue_->type();
}

conformation::ResidueCOP
ResidueAtomTreeCollection::active_residue_cop() const
{
	return active_residue_;
}

/// @brief The responsibility for making sure that the active residue and the active atomtree
/// are in synch is offloaded to an external class so that the calls to "active_residue()" and
/// "active_atom_tree()" are as fast as possible (and bitwise const for future multithreaded use).
/// After a round of set_chi() calls, the user for this class must update the residue coordinates.
void ResidueAtomTreeCollection::update_residue()
{
	using namespace id;
	using namespace conformation;
	using namespace kinematics;

	residue_uptodate_ = true;
	AtomTree const & tree( * atom_tree_representatives_[ active_restype_ ] );
	Residue & rsd( *residue_representatives_[ active_restype_ ] );
	if ( tree.residue_xyz_change_list_begin() != tree.residue_xyz_change_list_end() ) {
		debug_assert( *( tree.residue_xyz_change_list_begin() ) == 1 );
		for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
			rsd.set_xyz( ii, tree.xyz( id::AtomID( ii, 1 ) ));
		}

		tree.note_coordinate_change_registered();
	}

	// copy chi angles from atomtree->res
	for ( Size ii = 1, ii_end = rsd.nchi(); ii <= ii_end; ++ii ) {
		//std::cout << "Chi: " << ii << " " << tree.dof( DOF_ID( AtomID( rsd.chi_atoms(ii)[4], 1), PHI )) << " and ";
		rsd.chi()[ ii ] = numeric::constants::d::radians_to_degrees * tree.dof( DOF_ID( AtomID( rsd.chi_atoms(ii)[4], 1), PHI ));
		//std::cout << rsd.chi()[ ii ] << std::endl;
	}
	active_residue_->update_actcoord();
}

/// @brief See comments for update_residue().  After a call to "set_rescoords", the user must
/// call update_atomtree() to make sure the atomtree and the residue are in synch.
void ResidueAtomTreeCollection::update_atom_tree()
{
	atom_tree_uptodate_ = true;

	kinematics::AtomTree & tree( * atom_tree_representatives_[ active_restype_ ] );
	conformation::Residue & rsd( *residue_representatives_[ active_restype_ ] );

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		tree.set_xyz( id::AtomID( ii, 1 ), rsd.xyz( ii ) );
	}

	// force a dof update
	id::DOF_ID dofid( id::AtomID( 1, 1 ), id::RB1 );
	tree.dof( dofid );
}

core::Real ResidueAtomTreeCollection::dof( core::id::DOF_ID const &dofid )
{
	return atom_tree_representatives_[ active_restype_ ]->dof( dofid );
}

/// @brief Assigns the chi dihedral for the active restype.  Must be followed by a call to
/// update_residue() before the next call to active_residue()
void ResidueAtomTreeCollection::set_chi( Size chi_index, Real value )
{
	debug_assert( chi_index > 0 && chi_index <= residue_representatives_[ active_restype_ ]->nchi() );
	residue_uptodate_ = false;
	atom_tree_representatives_[ active_restype_ ]->set_dof(
		id::DOF_ID( id::AtomID( residue_representatives_[ active_restype_ ]->chi_atoms( chi_index )[ 4 ], 1 ), id::PHI ),
		numeric::constants::d::degrees_to_radians * value );
}

void ResidueAtomTreeCollection::set_d( Size chi_index, Real value )
{
	Size const effchi = (chi_index==0)? 1 : chi_index;
	Size const baseatom = (chi_index==0)? 3 : 4;

	debug_assert( effchi > 0 && effchi <= residue_representatives_[ active_restype_ ]->nchi() );
	residue_uptodate_ = false;
	atom_tree_representatives_[ active_restype_ ]->set_dof(
		id::DOF_ID( id::AtomID( residue_representatives_[ active_restype_ ]->chi_atoms( effchi )[ baseatom ], 1 ), id::D ),
		value  );
}

void ResidueAtomTreeCollection::set_theta( Size chi_index, Real value )
{
	Size const effchi = (chi_index==0)? 1 : chi_index;
	Size const baseatom = (chi_index==0)? 3 : 4;

	debug_assert( effchi > 0 && effchi <= residue_representatives_[ active_restype_ ]->nchi() );
	residue_uptodate_ = false;
	atom_tree_representatives_[ active_restype_ ]->set_dof(
		id::DOF_ID( id::AtomID( residue_representatives_[ active_restype_ ]->chi_atoms( effchi )[ baseatom ], 1 ), id::THETA ),
		numeric::constants::d::degrees_to_radians * value );
}

/// @brief Assigns the coordinates for a residue.  Must be followed by a call to
/// update_atom_tree() before the next cal to active_atom_tree().
void ResidueAtomTreeCollection::set_rescoords( conformation::Residue const & res )
{
	/// trust the the input residue's chi angles are correct
	debug_assert( & res.type() == & residue_representatives_[ active_restype_ ]->type() );
	for ( Size ii = 1; ii <= res.natoms(); ++ii ) {
		residue_representatives_[ active_restype_ ]->set_xyz( ii, res.xyz( ii ) );
	}
	for ( Size ii = 1, ii_end = res.nchi(); ii <= ii_end; ++ii ) {
		residue_representatives_[ active_restype_ ]->set_chi( ii, res.chi( ii ));
	}
	residue_representatives_[ active_restype_ ]->update_actcoord();
	atom_tree_uptodate_ = false;
}

/// @brief
void
ResidueAtomTreeCollection::set_rescoords( utility::vector1< Vector > const & coords )
{
	/// trust the the input residue's chi angles are correct
	debug_assert( coords.size() == residue_representatives_[ active_restype_ ]->natoms() );
	for ( Size ii = 1; ii <= residue_representatives_[ active_restype_ ]->natoms(); ++ii ) {
		residue_representatives_[ active_restype_ ]->set_xyz( ii, coords[ ii ] );
	}
	//residue_representatives_[ active_restype_ ]->update_actcoord();
	update_atom_tree();
	update_residue(); // now, copy the chi angles out of the atom tree
}

/// @brief
void
ResidueAtomTreeCollection::set_rescoords(
	utility::vector1< id::AtomID > const & atms, utility::vector1< Vector > const & coords )
{
	debug_assert( atms.size() == coords.size() );
	Size natoms = atms.size();

	/// trust the the input residue's chi angles are correct
	for ( Size ii = 1; ii <= natoms; ++ii ) {
		residue_representatives_[ active_restype_ ]->set_xyz( atms[ii].atomno(), coords[ ii ] );
	}
	//residue_representatives_[ active_restype_ ]->update_actcoord();
	update_atom_tree();
	update_residue(); // now, copy the chi angles out of the atom tree
}

void
ResidueAtomTreeCollection::save_momento( ResidueAtomTreeCollectionMomento & momento ) const
{
	debug_assert( residue_uptodate_ && atom_tree_uptodate_ );
	momento.set_restype_index( active_restype_ );
	momento.copy_coords( *residue_representatives_[ active_restype_ ] );
}

void ResidueAtomTreeCollection::update_from_momento( ResidueAtomTreeCollectionMomento const & momento )
{
	set_active_restype_index( momento.restype_index() );
	for ( Size ii = 1; ii <= active_residue_->natoms(); ++ii ) {
		active_residue_->set_xyz( ii, momento.coord( ii ) );
	}
	update_atom_tree();
	//update_residue();
	for ( Size ii = 1; ii <= active_residue_->natoms(); ++ii ) {
		debug_assert( active_residue_->xyz( ii ).distance_squared( momento.coord( ii ) ) < 1e-6 );
	}
	using namespace id;
	for ( Size ii = 1, ii_end = active_residue_->nchi(); ii <= ii_end; ++ii ) {
		//std::cout << "Chi: " << ii << " " << tree.dof( DOF_ID( AtomID( rsd.chi_atoms(ii)[4], 1), PHI )) << " and ";
		active_residue_->chi()[ ii ] = numeric::constants::d::radians_to_degrees * atom_tree_representatives_[ active_restype_ ]->dof( DOF_ID( AtomID( active_residue_->chi_atoms(ii)[4], 1), PHI ));
		/*std::cout << "   chi reported:" << active_residue_->chi()[ ii ] << " actual ";
		std::cout << numeric::dihedral_degrees( active_residue_->xyz( active_residue_->chi_atoms(ii)[ 1 ] ),
		active_residue_->xyz( active_residue_->chi_atoms(ii)[ 2 ] ),
		active_residue_->xyz( active_residue_->chi_atoms(ii)[ 3 ] ),
		active_residue_->xyz( active_residue_->chi_atoms(ii)[ 4 ] ) );
		std::cout << std::endl;
		debug_assert( std::abs( active_residue_->chi()[ ii ] -
		numeric::dihedral_degrees( active_residue_->xyz( active_residue_->chi_atoms(ii)[ 1 ] ),
		active_residue_->xyz( active_residue_->chi_atoms(ii)[ 2 ] ),
		active_residue_->xyz( active_residue_->chi_atoms(ii)[ 3 ] ),
		active_residue_->xyz( active_residue_->chi_atoms(ii)[ 4 ] )) < 1e-6 ) );*/
	}

}


AtomTreeCollection::AtomTreeCollection(
	pose::Pose const & pose,
	rotamer_set::RotamerSets const & rsets
) :
	//rotsets_( rotsets ),
	resid_2_moltenresid_( pose.size(), 0 ),
	moltenresid_2_resid_( rsets.nmoltenres() ),
	res_collections_( rsets.nmoltenres() )
{
	for ( Size ii = 1; ii <= rsets.nmoltenres(); ++ii ) {
		resid_2_moltenresid_[ rsets.moltenres_2_resid( ii ) ] = ii;
		moltenresid_2_resid_[ ii ] = rsets.moltenres_2_resid( ii );
		res_collections_[ ii ] = ResidueAtomTreeCollectionOP( new ResidueAtomTreeCollection(
			*rsets.rotamer_set_for_moltenresidue( ii ),
			moltenresid_2_resid_[ ii ]
			) );
	}
}


AtomTreeCollection::AtomTreeCollection(
	pose::Pose const & pose,
	task::PackerTask const & task
) :
	//rotsets_( rotsets ),
	resid_2_moltenresid_( pose.size(), 0 ),
	moltenresid_2_resid_( task.num_to_be_packed() ),
	res_collections_( task.num_to_be_packed() )
{
	Size count_moltenres( 0 );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( ! task.being_packed( ii ) ) continue;
		++count_moltenres;
		resid_2_moltenresid_[ ii ] = count_moltenres;
		moltenresid_2_resid_[ count_moltenres ] = ii;
		res_collections_[ count_moltenres ] = ResidueAtomTreeCollectionOP( new ResidueAtomTreeCollection( task.residue_task( ii ),  pose.conformation(), pose.residue( ii ) ) );
	}
}

AtomTreeCollection::AtomTreeCollection(
	pose::Pose const & /*pose*/,
	rotamer_set::RotamerSet const & rset,
	Size resid
)
{
	res_collections_.resize( 1 );
	res_collections_[ 1 ] = ResidueAtomTreeCollectionOP( new ResidueAtomTreeCollection( rset, resid  ) );
}

AtomTreeCollection::AtomTreeCollection(
	pose::Pose const & pose,
	task::ResidueLevelTask const & rltask,
	Size resid
)
{
	res_collections_.resize( 1 );
	res_collections_[ 1 ] = ResidueAtomTreeCollectionOP( new ResidueAtomTreeCollection( rltask, pose.conformation(), pose.residue( resid ) ) );
}

AtomTreeCollection::~AtomTreeCollection() {}

ResidueAtomTreeCollection &
AtomTreeCollection::moltenres_atomtree_collection( Size moltenresid )
{
	return * res_collections_[ moltenresid ];
}

ResidueAtomTreeCollection &
AtomTreeCollection::residue_atomtree_collection( Size resid )
{
	return * residue_atomtree_collection_op( resid );
}

ResidueAtomTreeCollectionOP
AtomTreeCollection::residue_atomtree_collection_op( Size resid )
{
	if ( resid_2_moltenresid_.size() != 0 ) {
		debug_assert( resid_2_moltenresid_[ resid ] != 0 );
		return res_collections_[ resid_2_moltenresid_[ resid ] ];
	} else {
		debug_assert( res_collections_.size() == 1 );
		return res_collections_[ 1 ];
	}
}

} // namespace scmin
} // namespace pack
} // namespace core
