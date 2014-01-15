// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/phosphate/PhosphateUtil.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/enumerate/rna/phosphate/PhosphateUtil.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.enumerate.rna.phosphate.PhosphateUtil" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace phosphate {

	/////////////////////////////////////////////////////////////////
	void
	remove_terminal_phosphates( pose::Pose & pose ){

		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( pose.residue(n).has_variant_type( "FIVE_PRIME_PHOSPHATE" ) ){
				remove_variant_type_from_pose_residue( pose, "FIVE_PRIME_PHOSPHATE", n );
				add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", n );
			}
			remove_variant_type_from_pose_residue( pose, "THREE_PRIME_PHOSPHATE", n );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	correctly_position_five_prime_phosphate_SLOW( pose::Pose & pose, Size const res ) {
		using namespace core::chemical;
		ResidueTypeSet const & rsd_set = pose.residue( res ).residue_type_set();
		conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *( rsd_set.aa_map( aa_from_name( "RAD") )[1] ) ) ;
		pose.prepend_polymer_residue_before_seqpos( *new_rsd, res, true );
		pose.delete_polymer_residue( res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	correctly_position_five_prime_phosphate( pose::Pose & pose, Size const res ) {
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::id;

		// reposition XO3' "manually".
		Residue const & rsd = pose.residue( res );
		Vector const & P_xyz = rsd.xyz( " P  " );
		Vector const & O5prime_xyz = rsd.xyz( " O5'" );
		Vector const & C5prime_xyz = rsd.xyz( " C5'" );
		Vector const & OP2_xyz     = rsd.xyz( " OP2" );
		Vector const & XO3prime_xyz  = rsd.xyz( "XO3'" );
		Real const OP2_dihedral      = dihedral_degrees( OP2_xyz,      P_xyz, O5prime_xyz, C5prime_xyz );
		Real const XO3prime_dihedral = dihedral_degrees( XO3prime_xyz, P_xyz, O5prime_xyz, C5prime_xyz );
		static Real const desired_XO3prime_OP2_offset = 114.6; // from rna_phenix files.
		Real const rotation_amount = ( desired_XO3prime_OP2_offset  - ( XO3prime_dihedral - OP2_dihedral ) );
		//		TR << "ROTATION AMOUNT " << rotation_amount << std::endl;
		AtomID XO3prime_ID = AtomID( rsd.atom_index( "XO3'"), res );
		DOF_ID XO3prime_DOF_ID(XO3prime_ID, PHI );
		//		pose.dump_pdb( "START.pdb" );
		//		TR << "DOF ORIGINAL: " << pose.dof( XO3prime_DOF_ID ) << std::endl;
		pose.set_dof( XO3prime_DOF_ID, pose.dof( XO3prime_DOF_ID ) + numeric::conversions::radians( rotation_amount ) );
		//		TR << "DOF NEW     : " << pose.dof( XO3prime_DOF_ID ) << std::endl;
		//		Vector const XO3prime_xyz_NEW  = pose.residue( res ).xyz( "XO3'" );
		//		Real const XO3prime_dihedral_NEW = dihedral_degrees( XO3prime_xyz_NEW, P_xyz, O5prime_xyz, C5prime_xyz );
		//		TR << "DESIRED OFFSET " << desired_XO3prime_OP2_offset << "  and actual: " << XO3prime_dihedral_NEW - OP2_dihedral << " from original: " << XO3prime_dihedral - OP2_dihedral << std::endl;
		//		pose.dump_pdb( "END.pdb" );
		//		exit( 0 );

	}

} //phosphate
} //rna
} //enumerate
} //stepwise
} //protocols
