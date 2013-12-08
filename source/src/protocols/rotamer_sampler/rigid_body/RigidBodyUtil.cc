// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyUtil.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rotamer_sampler/rigid_body/RigidBodyUtil.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.rigid_body.RigidBodyUtil" );

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_atom_coordinates( utility::vector1< std::pair < id::AtomID, numeric::xyzVector< core::Real > > > & xyz_list, Size const & seq_num, core::conformation::Residue const & rsd_at_origin, core::kinematics::Stub const & moving_res_base_stub ){

		xyz_list.clear();

		numeric::xyzVector< core::Real > const & new_centroid = moving_res_base_stub.v;
		numeric::xyzMatrix< core::Real > const & new_coordinate_matrix = moving_res_base_stub.M;

		for ( Size at = 1; at <= rsd_at_origin.natoms(); at++ ){

			id::AtomID const id( at, seq_num );

			numeric::xyzVector< core::Real > atom_pos;

			atom_pos = new_coordinate_matrix * rsd_at_origin.xyz( at ); //I think the order here does matter.
			atom_pos = atom_pos + new_centroid; //I think the order here does matter.

			xyz_list.push_back( std::make_pair( id, atom_pos ) );
		}
	}



} //rigid_body
} //rotamer_sampler
} //protocols
