// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Alex Ford (fordas@u.washington.edu)

// Unit headers
#include <protocols/hotspot_hashing/StubGenerator.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>

#include <protocols/hotspot_hashing/SearchPattern.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.io.hh>


namespace protocols {
namespace hotspot_hashing {

static THREAD_LOCAL basic::Tracer tr( "protocols.hotspot_hashing.StubGenerator" );

core::conformation::ResidueOP StubGenerator::getStubByName( std::string name )
{
	using namespace core::conformation;
	using namespace core::chemical;

	ResidueOP residue = ResidueFactory::create_residue(
		rsd_set_from_cmd_line().lock()->name_map( name ) ); // FIXME: potential null pointer from lock()

	moveFromStubFrame(residue, residueStubCentroidFrame(residue));

	return residue;
}

void StubGenerator::placeResidueAtTransform( core::pose::Pose & pose, core::conformation::ResidueCOP sourceResidue, core::kinematics::Stub transform, core::Size & residuejumpindex, core::Size & residueindex)
{
	tr.Debug << "Placing at transform: " << transform << std::endl;

	core::conformation::ResidueOP new_residue( new core::conformation::Residue(*sourceResidue) );
	moveIntoStubFrame(new_residue, transform);

	// Places residue at last jump & residue number
	placeResidueOnPose(pose, new_residue);
	residueindex = pose.size();
	residuejumpindex = pose.num_jump();

	tr.Debug << "Placed residue at anchor location: " << pose.residue(residueindex).xyz("CA") << std::endl;
}

void StubGenerator::placeResidueOnPose(core::pose::Pose & pose, core::conformation::ResidueCOP residue)
{
	pose.append_residue_by_jump(
		*residue,
		pose.size(),
		"",
		residue->atom_name(residue->nbr_atom()),
		true );

	if ( pose.num_jump() != 0 ) {
		// Adjust jump to realign residue
		core::Size atom_center;
		core::Size atom_a;
		core::Size atom_b;
		residue->select_orient_atoms(atom_center, atom_a, atom_b);

		core::kinematics::Stub source_stub(
			residue->xyz(atom_center),
			residue->xyz(atom_a),
			residue->xyz(atom_b));

		core::kinematics::Stub current_stub(
			pose.residue(pose.size()).xyz(atom_center),
			pose.residue(pose.size()).xyz(atom_a),
			pose.residue(pose.size()).xyz(atom_b));

		core::kinematics::RT transform(source_stub, current_stub);

		core::kinematics::Stub upstreamstub = pose.conformation().upstream_jump_stub(pose.num_jump());
		core::kinematics::Jump residuejump = pose.jump(pose.num_jump());
		core::kinematics::Jump newjump(residuejump);

		newjump.rotation_by_matrix(upstreamstub, source_stub.v, transform.get_rotation());
		if ( transform.get_translation().length() != 0 ) {
			newjump.translation_along_axis(upstreamstub, transform.get_translation(), transform.get_translation().length());
		}

		pose.set_jump(pose.num_jump(), newjump);
	}


	tr.Debug << "Appended residue on pose. Residue: " << pose.size() << " Jump: " << pose.num_jump() << " Anchor atom: " << residue->atom_name(residue->nbr_atom()) << std::endl;
}

void StubGenerator::moveIntoStubFrame( core::conformation::ResidueOP residue, core::kinematics::Stub transform )
{
	for ( core::Size i = 1; i <= residue->natoms(); i++ ) {
		residue->set_xyz(i, transform.local2global(residue->xyz(i)));
	}
}

void StubGenerator::moveFromStubFrame( core::conformation::ResidueOP residue, core::kinematics::Stub transform )
{
	for ( core::Size i = 1; i <= residue->natoms(); i++ ) {
		residue->set_xyz(i, transform.global2local(residue->xyz(i)));
	}
}

core::kinematics::Stub StubGenerator::residueStubOrientFrame(core::conformation::ResidueCOP const residue)
{
	core::Size a, b, c;
	residue->select_orient_atoms( a, b, c );
	core::kinematics::Stub result(residue->xyz(a), residue->xyz(b), residue->xyz(c));

	return result;
}

core::kinematics::Stub StubGenerator::residueStubCentroidFrame(core::conformation::ResidueCOP const residue)
{
	// Obtain the stub needed to generate
	// Canonical transform aligns CA atom to <0, 0, 0>
	// CA->SC heavyatom centroid vector along <1,0,0>
	// CA->C vector on the XY plane (<CA->C> * <0,0,1> == 0)

	core::kinematics::Stub result;

	result.from_four_points(
		residue->xyz(residue->atom_index("CA")),
		residueStubCentroid(residue),
		residue->xyz(residue->atom_index("CA")),
		residue->xyz(residue->atom_index("C")));

	return result;
}

StubGenerator::Vector StubGenerator::residueStubCentroid(core::conformation::ResidueCOP const residue)
{
	Vector centroid;
	centroid = 0;

	if ( residue->first_sidechain_atom() > residue->nheavyatoms() ) {
		return residue->xyz(residue->nbr_atom());
	}

	for ( core::Size i = residue->first_sidechain_atom(); i <= residue->nheavyatoms(); ++i ) {
		centroid += residue->xyz(i);
	}

	centroid /= (1 + residue->nheavyatoms() - residue->first_sidechain_atom());

	return centroid;
}

} // namespace hotspot_hashing
} // namespace protocols
