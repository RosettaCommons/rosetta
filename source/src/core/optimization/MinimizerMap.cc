// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/min.cc
/// @brief  Kinematics
/// @author Phil Bradley

// Unit headers
#include <core/optimization/MinimizerMap.hh>


// Project headers
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

// Numeric headers
// #include <numeric/all.fwd.hh>
// #include <numeric/constants.hh>
#include <numeric/conversions.hh>
// #include <numeric/xyzMatrix.hh>
// #include <numeric/xyzVector.hh>

// ObjexxFCL headers
// #include <ObjexxFCL/FArray1D.hh>
// #include <ObjexxFCL/FArray2D.hh>
// #include <ObjexxFCL/FArray3D.hh>
//#include <ObjexxFCL/FArray4D.h>
//#include <ObjexxFCL/formatted.io.h>

// Utility headers
// #include <utility/exit.hh>
// #include <utility/io/orstream.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

// C++ headers
#include <cstdlib>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>


namespace core {
namespace optimization {


/////////////////////////////////////////////////////////////////////////////
MinimizerMap::~MinimizerMap()
{
	clear_dof_nodes();
}

/////////////////////////////////////////////////////////////////////////////
// this will only work properly if the DOFs are sorted in decreasing
// depth
//
void
MinimizerMap::link_torsion_vectors()
{
#ifndef NDEBUG
	int last_depth = -1; // int depth( 100000 ); -- bad code buries assumptions deep deep into code where no one knows to find it
#endif
	for ( iterator it=dof_nodes_.begin(),
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
void
MinimizerMap::zero_torsion_vectors()
{
	for ( iterator iter=dof_nodes_.begin(), iter_end=dof_nodes_.end();
				iter != iter_end; ++iter ) {
		(*iter)->F1() = 0.0;
		(*iter)->F2() = 0.0;
	}
	for ( Size ii = 1; ii <= atom_derivatives_.size(); ++ii ) {
		for ( Size jj = 1; jj <= atom_derivatives_[ ii ].size(); ++jj ) {
			atom_derivatives_[ ii ][ jj ].f1() = 0.0;
			atom_derivatives_[ ii ][ jj ].f2() = 0.0;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////

void
MinimizerMap::add_torsion(
	DOF_ID const & dof_id,
	DOF_ID const & parent_id
)
{

	// assign the parent
	DOF_NodeOP parent;
	if ( parent_id.valid() ) {
		// not the root torsion
		parent = dof_node_pointer_[ parent_id ];
		if ( parent == 0 ) {
			std::cerr << "parent torsion does not exist in map! torsion= " <<
				dof_id << " parent= " << parent_id << std::endl;
			utility_exit();
		}
	} else {
		// root torsion!
		parent = 0;
	}

	DOF_NodeOP dof( new DOF_Node( dof_id, parent ) );
	dof_node_pointer_[ dof_id ] = dof;
	dof_nodes_.push_back( dof );

}

/////////////////////////////////////////////////////////////////////////////
void
MinimizerMap::add_atom(
	AtomID const & atom_id,
	DOF_ID const & dof_id
)
{
	// add atom to list of atoms first moved by this torsion,
	// unless we are in the root rigid body

	if ( dof_id.valid() ) {
		DOF_NodeOP n( dof_node_pointer_[ dof_id ] );
		if ( n == 0 ) {
			std::cout << "torsion does not exist in map! torsion= " <<
				dof_id << std::endl;
			utility_exit();
		}

		// add the atom to the torsion
		n->add_atom( atom_id );
	}
}

/////////////////////////////////////////////////////////////////////////////
Real
MinimizerMap::torsion_scale_factor(
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
MinimizerMap::clear_dof_nodes()
{
	// delete old memory -- deprecated 7/1/2010 following DOF_Node's inherretance from ReferenceCount
	//for ( iterator it=dof_nodes_.begin(), it_end = dof_nodes_.end();
	//			it != it_end; ++it ) {
	//	delete (*it);
	//}

	dof_nodes_.clear();
	dof_node_pointer_.clear(); // have to clear this map, too
}

/////////////////////////////////////////////////////////////////////////////
void
MinimizerMap::reset( pose::Pose const & pose )
{
	clear_dof_nodes();

	// dimension the node_pointers to allow fast lookup
	DOF_NodeOP tmp(0);
	pose::initialize_dof_id_map( dof_node_pointer_, pose, tmp );
	atom_derivatives_.resize( pose.total_residue() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		atom_derivatives_[ ii ].resize( pose.residue( ii ).natoms() );
	}
}

/////////////////////////////////////////////////////////////////////////////
void
MinimizerMap::copy_dofs_from_pose(
	pose::Pose const & pose,
	Multivec & dofs
) const
{
	int imap = 1;
	for ( const_iterator it=dof_nodes_.begin(), it_end = dof_nodes_.end();
				it != it_end; ++it, ++imap ) {
		DOF_Node const & dof_node( **it );
		dofs[ imap ] = torsion_scale_factor( dof_node ) *
			pose.dof( dof_node.dof_id() );
	}
}


/////////////////////////////////////////////////////////////////////////////
void
MinimizerMap::copy_dofs_to_pose(
	pose::Pose & pose,
	Multivec const & dofs
) const
{
	int imap = 1;
	for ( const_iterator it=dof_nodes_.begin(), it_end = dof_nodes_.end();
			it != it_end; ++it, ++imap ) {
		DOF_Node const & dof_node( **it );
		pose.set_dof( dof_node.dof_id(),
			dofs[ imap ] / torsion_scale_factor( dof_node ));
	}
}


/////////////////////////////////////////////////////////////////////////////
void
MinimizerMap::reset_jump_rb_deltas(
	pose::Pose & pose,
	Multivec & dofs
) const
{
	int imap = 1;
	for ( const_iterator it=dof_nodes_.begin(), it_end = dof_nodes_.end();
			it != it_end; ++it, ++imap ) {
		DOF_Node const & dof_node( **it );
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

/////////////////////////////////////////////////////////////////////////////
bool
DOF_Node_sorter( DOF_NodeCOP a, DOF_NodeCOP b )
{
	return *a < *b;
}

/////////////////////////////////////////////////////////////////////////////

void
MinimizerMap::setup(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map
)
{
	// this clears dof_nodes_ and dimensions dof_node_pointer_
	reset( pose );

	// convert the allow_bb,allow_chi,allow_jump information
	// in the MoveMap into a simple boolean mask over AtomTree
	// degrees of freedom
	// this is necessary because the AtomTree doesn't know which degrees
	// of freedom are chi angles, which are phi angles, etc
	// the dof_mask is a low-level type of allow-move information
	// that's used by the atomtree atoms inside Atom::setup_min_map
	//
	id::DOF_ID_Mask dof_mask( false );
	pose::setup_dof_mask_from_move_map( move_map, pose, dof_mask );

	// this fills the torsion and atom lists
	DOF_ID tmp( id::BOGUS_DOF_ID );
	pose.atom_tree().root()->setup_min_map( tmp, dof_mask, *this );

	// sort DOFs for proper linking later on
	dof_nodes_.sort( DOF_Node_sorter );

	// identify phi/psi/omega...
	//
	this->assign_rosetta_torsions( pose );


	// setup the domain_map which indicates what rsd pairs are fixed/moving
	id::AtomID_Mask moving_dof, moving_xyz;
	core::pose::initialize_atomid_map( moving_xyz, pose, false );
	core::pose::initialize_atomid_map( moving_dof, pose, false );
	for ( const_iterator it = dof_nodes_.begin(), it_end = dof_nodes_.end();
				it != it_end; ++it ) {
		moving_dof[ (**it).atom_id() ] = true;
	}

	domain_map_.dimension( pose.total_residue() );
	pose.conformation().atom_tree().update_domain_map
		( domain_map_, moving_dof, moving_xyz );
}


/////////////////////////////////////////////////////////////////////////////
// private
void
MinimizerMap::assign_rosetta_torsions( pose::Pose const & pose )
{
	// mapping from AtomTree DOF ID's to bb/chi torsion angle ids
	id::DOF_ID_Map< id::TorsionID > dof_map
		( id::BOGUS_TORSION_ID);

	pose::setup_dof_to_torsion_map( pose, dof_map );

	for ( iterator it = dof_nodes_.begin(), it_end = dof_nodes_.end();
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


} // namespace kinematics
} // namespace core
