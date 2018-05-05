// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/GridScorer.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/GriddedAtomTreeMultifunc.hh>

#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/DOF_Node.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>
#include <core/optimization/atom_tree_minimize.hh>
#include <core/pose/Pose.hh>
#include <core/id/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

//////////////////////////////
namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

//using basic::T;
using basic::Tracer;
static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.GriddedAtomTreeMultifunc" );


GriddedAtomTreeMultifunc::GriddedAtomTreeMultifunc(
	LigandConformer & conf_in,
	core::pose::Pose & pose_in,
	GridScorer & scorefxn_in,
	core::optimization::MinimizerMap & min_map_in
) :
	conf_( conf_in ),
	pose_( pose_in ),
	sf_( scorefxn_in ),
	min_map_( min_map_in )
{
}

GriddedAtomTreeMultifunc::~GriddedAtomTreeMultifunc() { }

//
core::Real
GriddedAtomTreeMultifunc::operator ()( core::optimization::Multivec const & vars ) const {
	// score ligand conformer
	min_map_.copy_dofs_to_pose( pose_, vars );
	core::Real score = sf_.score( pose_, conf_ );
	return score;
}

//
void
GriddedAtomTreeMultifunc::dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const {
	using numeric::constants::d::pi;

	dE_dvars.resize( min_map_.nangles() );

	min_map_.zero_torsion_vectors();
	min_map_.copy_dofs_to_pose( pose_, vars );

	// get atom derivatives
	sf_.derivatives( pose_, conf_, min_map_ );

	// DEBUG DERIVS[1]
	//sf_.debug_deriv( pose_, conf_, min_map_ );

	// store f1/f2 on DOF nodes
	for ( auto iter = min_map_.begin(), iter_e = min_map_.end(); iter != iter_e; ++iter ) {
		core::optimization::DOF_Node & dof_node( **iter );
		for ( auto it1=dof_node.atoms().begin(),
				it1e = dof_node.atoms().end(); it1 != it1e; ++it1 ) {
			core::id::AtomID const & atom_id( *it1 );
			dof_node.F1() += min_map_.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f1();
			dof_node.F2() += min_map_.atom_derivatives( atom_id.rsd() )[ atom_id.atomno() ].f2();
		}
	}

	// pass f1/f2 up tree
	min_map_.link_torsion_vectors();

	// convert f1s/f2s into torsional derivatives
	int imap( 1 );
	for ( auto it=min_map_.begin(), ite=min_map_.end(); it != ite; ++it, ++imap ) {
		core::optimization::DOF_Node const & dof_node( **it );
		core::kinematics::tree::Atom const & atom( pose_.atom_tree().atom( dof_node.atom_id() ) );
		dE_dvars[ imap ] = 0.0;

		numeric::xyzVector< core::Real > axis, end_pos;
		core::id::DOF_Type const type( dof_node.type() );
		atom.get_dof_axis_and_end_pos( axis, end_pos, type );

		if ( type == core::id::PHI || type == core::id::THETA || type == core::id::RB4 || type == core::id::RB5 || type == core::id::RB6 ) {
			core::Real scale_factor( ( type == core::id::PHI || type == core::id::THETA ) ? 1 : numeric::constants::d::deg2rad );
			if ( type == core::id::THETA ) {
				core::Real const theta( atom.dof( type ) );
				int const theta_mod ( ( static_cast< int >( std::floor( theta/pi )))%2);
				if ( theta_mod == 1 || theta_mod == -1 ) scale_factor *= -1.0f;
			}

			dE_dvars[ imap ] -= scale_factor * ( dot( axis, dof_node.F1() ) + dot( cross( axis, end_pos ), dof_node.F2() ) );
		} else {
			dE_dvars[ imap ] += dot( axis, dof_node.F2() );
		}

		core::Real tors_scalefactor = min_map_.torsion_scale_factor( dof_node );
		dE_dvars[ imap ] += sf_.dof_derivative( pose_, min_map_, dof_node.dof_id(), dof_node.torsion_id() );
		dE_dvars[ imap ] /= tors_scalefactor;
	}

	// DEBUG DERIVS[2]
	//core::optimization::NumericalDerivCheckResultOP deriv_check_result( new core::optimization::NumericalDerivCheckResult );
	//deriv_check_result->send_to_stdout( true );
	//core::optimization::numerical_derivative_check( min_map_, *this, vars, dE_dvars, deriv_check_result, true );
}

void
GriddedAtomTreeMultifunc::dump( core::optimization::Multivec const & /*vars*/, core::optimization::Multivec const & /*vars2*/ ) const {
	//TR << "Calling GriddedAtomTreeMultifunc::dump()" << std::endl;
}




} // ga_dock
} // ligand_docking
} // protocols
