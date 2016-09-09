// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/SymMinimizerMap.cc
/// @brief  MinimizerMap for symmetric minimization implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/optimization/symmetry/SymMinimizerMap.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/DOF_Node.hh>

// Project headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>


/// Numeric headers
#include <numeric/conversions.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <utility/options/BooleanVectorOption.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace core {
namespace optimization {
namespace symmetry {

/////////////////////////////////////////////////////////////////////////////
/// Stolen from MinimizerMap.cc
bool
DOF_Node_sorter( DOF_NodeCOP a, DOF_NodeCOP b )
{
	return *a < *b;
}

/////////////////////////////////////////////////////////////////////////////

SymMinimizerMap::SymMinimizerMap(
	pose::Pose const & pose,
	kinematics::MoveMap const & asymm_movemap,
	SymmetryInfoCOP symm_info,
	bool const new_sym_min // =false
) :
	pose_( pose ),
	symm_info_( symm_info ),
	res_interacts_with_asymmetric_unit_( pose.size(), false ),
	// n_dof_nodes_( 0 ),
	n_independent_dof_nodes_( 0 ),
	atom_derivatives_( pose.size() ),
	new_sym_min_( new_sym_min )
{

	//std::cout << "SymMinimizerMap ctor: " << std::endl;
	//std::cout << pose.conformation().fold_tree() << std::endl;

	DOF_NodeOP tmp(0);
	pose::initialize_dof_id_map( dof_node_pointer_, pose, tmp );

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		atom_derivatives_[ ii ].resize( pose.residue( ii ).natoms() );
		if ( symm_info_->bb_follows( ii ) == 0 ) {
			res_interacts_with_asymmetric_unit_[ ii ] = true; // residues in the asymmetric unit are always included
			continue;
		}
		/// ASSUMPTION the energy graph only includes edges where at least one node belongs to the asymmetric unit
		for ( utility::graph::Node::EdgeListConstIter
				iter = pose.energies().energy_graph().get_node( ii )->const_edge_list_begin(),
				iter_end = pose.energies().energy_graph().get_node( ii )->const_edge_list_end();
				iter != iter_end; ++iter ) {
			Size jj = (*iter)->get_other_ind( ii );
			if ( symm_info_->bb_follows( jj ) == 0 ) {
				res_interacts_with_asymmetric_unit_[ ii ] = true;
				break;
			}
		}
	}

	/// tell the atom tree about every DOF that will be changing, even if some of the DOFs are clones
	id::DOF_ID_Mask dof_mask( false );
	pose::setup_dof_mask_from_move_map( asymm_movemap, pose, dof_mask );

	/// figure out the mapping between dof_id's and torsion_id's, this is needed to correctly tell if a dof is
	/// independent or dependent
	pose::setup_dof_to_torsion_map( pose, dof_id2torsion_id_ );

	// this fills the torsion and atom lists
	DOF_ID temp( id::BOGUS_DOF_ID );
	pose.atom_tree().root()->setup_min_map( temp, dof_mask, *this );

	/// sort dof nodes by tree depth.
	dof_nodes_.sort( DOF_Node_sorter );


	// identify phi/psi/omega...
	//
	this->assign_rosetta_torsions(); // uses dof_id2torsion_id_ mapping that was set up above


	/// STOLEN directly from MinimizerMap.cc
	// setup the domain_map which indicates what rsd pairs are fixed/moving
	id::AtomID_Mask moving_dof, moving_xyz;
	core::pose::initialize_atomid_map( moving_xyz, pose, false );
	core::pose::initialize_atomid_map( moving_dof, pose, false );
	for ( const_iterator it = dof_nodes_.begin(), it_end = dof_nodes_.end();
			it != it_end; ++it ) {
		moving_dof[ (**it).atom_id() ] = true;
	}

	domain_map_.dimension( pose.size() );
	pose.conformation().atom_tree().update_domain_map
		( domain_map_, moving_dof, moving_xyz );

}

SymMinimizerMap::~SymMinimizerMap() {}

void
SymMinimizerMap::add_torsion(
	DOF_ID const & new_torsion,
	DOF_ID const & parent
)
{
	// std::cout << "add torsion: " << new_torsion.atomno() << " " << new_torsion.rsd() << " " << new_torsion.type();// << std::endl;
	// std::cout << "add torsion parent: " << parent.atomno() << " " << parent.rsd() << " " << parent.type() << std::endl;

	/// which kind of torsion is this guy
	id::TorsionID const & tor_id( dof_id2torsion_id_[ new_torsion ] );

	bool const independent( tor_id.valid() ? symm_info_->torsion_is_independent( tor_id ) :
		symm_info_->dof_is_independent( new_torsion, pose_.conformation() ) );

	// std::cout << "add torsion: " << new_torsion.atomno() << " " << new_torsion.rsd() << " " << new_torsion.type() <<
	//  ' ' << tor_id <<
	//  " ind: " << independent << " parent: " << parent.atomno() << " " << parent.rsd() << " " << parent.type() << std::endl;

	if ( new_sym_min_ ) {
		// add everything in the new approach
		add_new_dof_node( new_torsion, parent, !independent );
	} else {
		if ( independent ) {
			//std::cout << " ind"  << std::endl;
			add_new_dof_node( new_torsion, parent, false );
			//last_cloned_jump_ = DOF_ID( id::BOGUS_DOF_ID ); // set invalid
		} else {

			if ( new_torsion.type() >= id::RB1 ) { // We have a jump
				conformation::symmetry::SymmetricConformation const & symm_conf (
					dynamic_cast< conformation::symmetry::SymmetricConformation const & > ( pose_.conformation()) );
				debug_assert( conformation::symmetry::is_symmetric( symm_conf ) );
				if ( symm_info_->get_dof_derivative_weight( new_torsion, symm_conf ) != 0.0 ) {
					DOF_ID parent_dof( id::BOGUS_DOF_ID );
					add_new_dof_node( new_torsion, parent_dof, true );
				}
			}
		}
	}
}

void
SymMinimizerMap::add_atom(
	AtomID const & atom_id,
	DOF_ID const & dof_id
)
{
	if ( ! dof_id.valid() ) return;

	//std::cout << "add atom? atom " << atom_id << " dof " << dof_id.rsd() << " " << dof_id.atomno() << " " << dof_id.type() << std::endl;
	//std::cout << "    last_cloned_jump_?: "<< last_cloned_jump_.rsd() << " " << last_cloned_jump_.atomno() << " " << last_cloned_jump_.type() << std::endl;

	debug_assert( dof_node_pointer_[ dof_id ] );
	dof_node_pointer_[ dof_id ]->add_atom( atom_id );

	/*if ( symm_info_->dof_is_independent( dof_id, pose_.conformation() ) ) {
	dof_node_pointer_[ dof_id ]->add_atom( atom_id );
	std::cout << "    added to independ dof" << std::endl;
	} else if ( dof_node_pointer_[ dof_id ] ) {
	//debug_assert( dof_id == last_cloned_jump_ );
	debug_assert( ! symm_info_->dof_is_independent( dof_id, pose_.conformation() ) );
	//debug_assert( dof_node_pointer_[ last_cloned_jump_ ] );
	dof_node_pointer_[ dof_id ] ->add_atom( atom_id );
	std::cout << "    added to last_cloned_jump_: "<< last_cloned_jump_.rsd() << " " << last_cloned_jump_.atomno() << " " << last_cloned_jump_.type() << std::endl;
	}  */
}

kinematics::DomainMap const &
SymMinimizerMap::domain_map() const
{
	return domain_map_;
}

void
SymMinimizerMap::copy_dofs_from_pose(
	pose::Pose const & pose,
	Multivec & dofs
) const
{
	int imap = 1;
	for ( const_iterator it=dof_nodes_.begin(), it_end = dof_nodes_.end();
			it != it_end; ++it, ++imap ) {
		DOF_Node const & dof_node( **it );
		if ( new_sym_min_ && dof_node.dependent() ) { --imap; continue; }
		dofs[ imap ] = torsion_scale_factor( dof_node ) *
			pose.dof( dof_node.dof_id() );
	}
}

void
SymMinimizerMap::copy_dofs_to_pose(
	pose::Pose & pose,
	Multivec const & dofs
) const
{
	int imap = 1;
	for ( const_iterator it=dof_nodes_.begin(), it_end = dof_nodes_.end();
			it != it_end; ++it, ++imap ) {
		DOF_Node const & dof_node( **it );
		if ( new_sym_min_ && dof_node.dependent() ) { --imap; continue; }
		pose.set_dof( dof_node.dof_id(),
			dofs[ imap ] / torsion_scale_factor( dof_node ));
	}
}

DOF_NodeOP
SymMinimizerMap::dof_node_from_id( DOF_ID const &id ) const
{
	DOF_NodeOP node = 0;
	if ( id.valid() ) {
		node = dof_node_pointer_[ id ];
		if ( node == 0 ) {
			std::cerr << "DOF_ID does not exist in map! torsion= " << id << std::endl;
			utility_exit();
		}
	}
	return node;
}


void SymMinimizerMap::zero_torsion_vectors()
{
	for ( const_iterator iter = dof_nodes_.begin(),
			iter_end = dof_nodes_.end(); iter != iter_end; ++iter ) {
		(*iter)->F1() = 0;
		(*iter)->F2() = 0;
	}
	for ( Size ii = 1, iiend = atom_derivatives_.size(); ii <= iiend; ++ii ) {
		for ( Size jj = 1, jjend = atom_derivatives_[ ii ].size(); jj <= jjend; ++jj ) {
			atom_derivatives_[ ii ][ jj ].f1() = 0.0;
			atom_derivatives_[ ii ][ jj ].f2() = 0.0;
		}
	}
}

void
SymMinimizerMap::link_torsion_vectors()
{
#ifndef NDEBUG
	int last_depth = -1;
#endif
	for ( const_iterator it=dof_nodes_.begin(),
			it_end=dof_nodes_.end(); it != it_end; ++it ) {
		DOF_Node & dof_node( **it );
		debug_assert( last_depth == -1 || dof_node.depth() <= last_depth );
#ifndef NDEBUG
		last_depth = dof_node.depth();
#endif
		dof_node.link_vectors();
	}
}

/////////////////////////////////////////////////////////////////////////////
Real
SymMinimizerMap::torsion_scale_factor(
	DOF_Node const & dof_node
) const
{
	static Real const rad2deg( numeric::conversions::degrees(1.0) );
	DOF_Type const type( dof_node.type() );
	Real factor( 1.0 );
	if ( type == id::PHI ) {
		// bond torsion
		factor = rad2deg;
	} else if ( type == id::THETA ) {
		// bond angle
		factor = rad2deg * basic::options::option[ basic::options::OptionKeys::optimization::scale_theta ]();
	} else if ( type == id::D ) {
		// bond length
		factor = basic::options::option[ basic::options::OptionKeys::optimization::scale_d ]();
	} else if ( type == id::RB4 ||
			type == id::RB5 ||
			type == id::RB6 ) {
		// the jump_rb_delta's are stored in degrees!!!
		factor = basic::options::option[ basic::options::OptionKeys::optimization::scale_rbangle ]();
	} else if ( type == id::RB1 ||
			type == id::RB2 ||
			type == id::RB3 ) {
		// rigid body translation
		factor = basic::options::option[ basic::options::OptionKeys::optimization::scale_rb ]();
	}
	return factor;
}

/////////////////////////////////////////////////////////////////////////////
void
SymMinimizerMap::reset_jump_rb_deltas(
	pose::Pose & pose,
	Multivec & dofs
) const
{
	int imap = 1;
	for ( const_iterator it=dof_nodes_.begin(), it_end = dof_nodes_.end();
			it != it_end; ++it, ++imap ) {
		DOF_Node const & dof_node( **it );
		if ( new_sym_min_ && dof_node.dependent() ) { --imap; continue; }
		if ( DOF_type_is_rb( dof_node.type() ) ) {
			// will do this multiple times for each jump, but should be OK
			AtomID const & id( dof_node.atom_id() );
			kinematics::Jump jump( pose.jump( id ) );
			jump.fold_in_rb_deltas();
			pose.set_jump( id, jump );
			dofs[ imap ] = 0.0;
		}
	}
}


void
SymMinimizerMap::add_new_dof_node(
	DOF_ID const & new_torsion,
	DOF_ID const & parent,
	bool dependent
)
{
	DOF_NodeOP dof_node( new DOF_Node( new_torsion, DOF_NodeOP( 0 ) ) );
	dof_node->dependent( false ); // only used if new_sym_min_

	if ( parent.valid() ) {
		DOF_NodeOP parent_ptr = dof_node_pointer_[ parent ];
		if ( parent_ptr ) {
			dof_node->set_parent( parent_ptr );
		}
	}

	dof_node_pointer_[ new_torsion ] = dof_node;
	if ( new_sym_min_ ) { // new way
		dof_node->dependent( dependent );
		dof_nodes_.push_back( dof_node );
		if ( !dependent ) ++n_independent_dof_nodes_;
	} else { // old way
		if ( dependent ) {
			dependent_dof_nodes_.push_back( dof_node );// ++n_dof_nodes_;
		} else {
			dof_nodes_.push_back( dof_node ); ++n_independent_dof_nodes_;// ++n_dof_nodes_;
		}
	}
}

/// @details use the bb_follows mapping to find the residue in the asymmetric unit that
/// this cloned atom corresponds to
id::DOF_ID
SymMinimizerMap::asymmetric_dof( DOF_ID const & cloned_dof ) const
{
	if ( cloned_dof.type() >= id::RB1 && cloned_dof.type() <= id::RB6 ) {
		/// we have a jump
		Size cloned_jumpno = pose_.conformation().fold_tree().get_jump_that_builds_residue( cloned_dof.rsd() );
		debug_assert( symm_info_->jump_follows( cloned_jumpno ) != 0 );
		Size asymm_jumpno = symm_info_->jump_follows( cloned_jumpno );
		Size asymm_resno = pose_.conformation().fold_tree().jump_edge( asymm_jumpno ).stop();
		id::AtomID asymmatom( cloned_dof.atomno(), asymm_resno );
		return DOF_ID( asymmatom, cloned_dof.type() );
	} else {
		debug_assert( symm_info_->bb_follows( cloned_dof.rsd() ) != 0 );
		id::AtomID asymmatom( cloned_dof.atomno(), symm_info_->bb_follows( cloned_dof.rsd() ) );

		return DOF_ID( asymmatom, cloned_dof.type() );
	}
}

void
SymMinimizerMap::assign_rosetta_torsions()
{
	// mapping from AtomTree DOF ID's to bb/chi torsion angle ids
	//id::DOF_ID_Map< id::TorsionID > dof_map
	// ( id::BOGUS_TORSION_ID);
	//pose::setup_dof_to_torsion_map( pose, dof_map ); // we already did this
	id::DOF_ID_Map< id::TorsionID > const & dof_map( dof_id2torsion_id_ );

	for ( const_iterator it = dof_nodes_.begin(), it_end = dof_nodes_.end();
			it != it_end; ++it ) {
		DOF_Node & dof_node( **it );

		if ( dof_node.type() == id::PHI ) {
			id::TorsionID const & id( dof_map[ dof_node.dof_id() ] );
			if ( id.valid() ) {
				dof_node.torsion_id( id );
			}
		}
	}
}

} // symmetry
} // namespace optimization
} // namespace core
