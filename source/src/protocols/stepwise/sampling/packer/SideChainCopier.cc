// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/packer/SideChainCopier.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/packer/SideChainCopier.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/copydofs/util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.packer.SideChainCopier" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace packer {

	//Constructor
	SideChainCopier::SideChainCopier( core::pose::Pose const & reference_pose,
																		bool const copy_o2prime_hydrogens /* = false */ ):
		reference_pose_( reference_pose ),
		copy_o2prime_hydrogens_( copy_o2prime_hydrogens )
	{
		for ( Size n = 1; n <= reference_pose_.total_residue(); n++ ) copy_res_.push_back( n );
	}

	//Constructor
	SideChainCopier::SideChainCopier( core::pose::Pose const & reference_pose,
																		utility::vector1< Size > const & copy_res,
																		bool const copy_o2prime_hydrogens /* = false */ ):
		reference_pose_( reference_pose ),
		copy_res_( copy_res ),
		copy_o2prime_hydrogens_( copy_o2prime_hydrogens )
	{
	}

	//Destructor
	SideChainCopier::~SideChainCopier()
	{}

	//////////////////////////////////////////////////////////////////////////////////
	void
	SideChainCopier::apply( core::pose::Pose & viewer_pose ){

		Pose pose = viewer_pose; // to prevent crashes in graphics.

		runtime_assert( pose.total_residue() == reference_pose_.total_residue() );
		Real const ho2prime_tolerance( 1.0e-3 );
		Real const orient_backbone_tolerance( 1.0e-3 );

		for ( Size i = 1; i <= copy_res_.size(); i++ ){

			Size const n = copy_res_[i];
			if ( pose.residue_type( n ).is_RNA() ) {
				if ( !copy_o2prime_hydrogens_ ) continue;
				make_variants_match( pose, reference_pose_, n, "VIRTUAL_O2PRIME_HYDROGEN" );
				id::TorsionID torsion_id( i, id::CHI, 4 ); // 2'-OH
				if ( std::abs( pose.torsion( torsion_id ) - reference_pose_.torsion( torsion_id ) ) < ho2prime_tolerance ) continue;
				pose.set_torsion( torsion_id, reference_pose_.torsion( torsion_id ) );

			}	else if ( pose.residue_type( n ).is_protein() ) {

				ResidueOP rsd = reference_pose_.residue( n ).clone();

				bool const orient_backbone = ( ( reference_pose_.residue( n ).xyz( 1 ) - pose.residue( n ).xyz( 1 ) ).length() > orient_backbone_tolerance );
				if ( orient_backbone ) rsd->place( pose.residue( n ), pose.conformation() );

				// before putting into new pose, need to move backbone atoms back, so that only side chain atoms are changed.
				for ( Size k = 1; k <= rsd->type().last_backbone_atom(); k++ ){
					rsd->set_xyz( k, pose.residue( n ).xyz( rsd->atom_name( k ) ) );
				}
				for ( Size k = rsd->type().nheavyatoms() + 1; k < rsd->type().first_sidechain_hydrogen(); k++ ){
					rsd->set_xyz( k, pose.residue( n ).xyz( rsd->atom_name( k ) ) );
				}

				pose.conformation().replace_residue( n, *rsd, false /*orient_backbone*/ );
			}

		}

		viewer_pose = pose;
	}

} //packer
} //sampling
} //stepwise
} //protocols
