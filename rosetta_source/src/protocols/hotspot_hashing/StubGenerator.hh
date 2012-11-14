// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Alex Ford (fordas@u.washington.edu)

#ifndef INCLUDED_protocols_hotspot_hashing_StubGenerator_hh
#define INCLUDED_protocols_hotspot_hashing_StubGenerator_hh

// Unit headers
#include <protocols/hotspot_hashing/StubGenerator.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace hotspot_hashing {

class StubGenerator
{
  public:
    StubGenerator();

    static core::conformation::ResidueOP getStubByName( std::string name )
    {
      using namespace core::conformation;
      using namespace core::chemical;
      
      ResidueOP residue = ResidueFactory::create_residue(
          rsd_set_from_cmd_line()->name_map( name ) );

			applyRT(residue, residueStubCentroidTransform(residue));
    }

    static void applyRT( core::conformation::ResidueOP residue, RT transform )
		{
			for (core::Size i = 1; i <= residue->natoms(); i++) 
			{
				residue->set_xyz(
						i,
						transform.get_rotation() * (residue->xyz(i) + transform.get_translation()));
			}
		}

		static RT residueStubCentroidTransform(core::conformation::ResidueCOP const residue)
		{
			// Canonical transform aligns CA atom to <0, 0, 0>
			// CA->SC heavyatom centroid vector along <1,0,0>
			// CA->C vector on the XY plane (<CA->C> * <0,0,1> == 0)
			
			Vector position = -residue->xyz(residue->atom_index("CA"));
			
			Vector xunit = Vector(1, 0, 0);
			Vector yunit = Vector(0, 1, 0);

			Vector cacentroid_vector = residueStubCentroid(residue) + position;
			Matrix cacentroid_rotation = rotation_matrix( cacentroid_vector.cross(xunit), angle_of(cacentroid_vector, xunit));

			Vector cac_vector_prime = cacentroid_rotation * (residue->xyz(residue->atom_index("C")) + position);
			Vector cac_vector_zyprojection = Vector(0, cac_vector_prime.y(), cac_vector_prime.z());
			Matrix cac_rotation = rotation_matrix( cac_vector_zyprojection.cross(yunit), angle_of(cac_vector_zyprojection, yunit));

			return RT(cac_rotation * cacentroid_rotation, position);
		}

		static Vector residueStubCentroid(core::conformation::ResidueCOP const residue)
		{
			Vector centroid;
			centroid = 0;

			if (residue->first_sidechain_atom() > residue->nheavyatoms())
			{
				return residue->xyz(residue->nbr_atom());
			}

			for (core::Size i = residue->first_sidechain_atom(); i <= residue->nheavyatoms(); ++i)
			{
				centroid += residue->xyz(i);
			}

			centroid /= (1 + residue->nheavyatoms() - residue->first_sidechain_atom());

			return centroid;
		}
};

} // namespace hotspot_hashing
} // namespace protocols

#endif
