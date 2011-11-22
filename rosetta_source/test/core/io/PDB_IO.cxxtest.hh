// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/PDB_IO.cxxtest.hh
/// @brief  test suite for PDB reader/writer
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include "platform/types.hh"

#include <test/core/init_util.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>


// Package Headers

#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.io.PDB_IO.cxxtest");

using namespace core;

class PDB_IO : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	PDB_IO() {}


	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-no_optH" );
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// ------------------------------------------ //
	/// @brief test how PDB IO
	void test_pdb_io() {
		core::pose::Pose pose;
		const std::string original_file_name("core/io/test_in.pdb");
		core::import_pose::pose_from_pdb(pose, original_file_name);

		// see if number of residues are correct
		TS_ASSERT_EQUALS( pose.total_residue(), 116u );

		// write pose to a new file...
		const std::string tmp_file_name("PDB_IO_cxxtest.pdb._tmp_");
		core::io::pdb::dump_pdb(pose, tmp_file_name);
		//pose.dump_pdb(tmp_file_name);

		TS_ASSERT_FILE_EQ(original_file_name.c_str(), tmp_file_name.c_str());

		// read writen file as new pose object
		core::pose::Pose P2;
		core::import_pose::pose_from_pdb(P2, tmp_file_name);

		// see if number of residues are correct
		TS_ASSERT_EQUALS( P2.total_residue(), 116u );


		// now test if residue information are the same
		bool isExit = false;
		for(Size i=1; i<=pose.total_residue(); i++) {
			core::conformation::Residue const & R1( pose.residue(i) );
			core::conformation::Residue const & R2( P2.residue(i) );

			// check if number of atoms are the same
			TS_ASSERT_EQUALS( R1.natoms(), R2.natoms());

			TS_ASSERT_EQUALS( R1.name3(), R2.name3());

			for(Size a=1; a<=R1.natoms(); a++) {
				conformation::Atom const & A1( R1.atom(a) );
				conformation::Atom const & A2( R2.atom(a) );

				TS_ASSERT_EQUALS( R1.atom_name(a), R2.atom_name(a));
				if( R1.atom_name(a) != R2.atom_name(a) ) isExit=true;

				TS_ASSERT_EQUALS( A1.xyz().x(), A2.xyz().x());
				TS_ASSERT_EQUALS( A1.xyz().y(), A2.xyz().y());
				TS_ASSERT_EQUALS( A1.xyz().z(), A2.xyz().z());

				if( A1.xyz().x() != A2.xyz().x() ) isExit=true;
				if( A1.xyz().y() != A2.xyz().y() ) isExit=true;
				if( A1.xyz().z() != A2.xyz().z() ) isExit=true;

				if( isExit ) break;
			}
			if( isExit ) break;
		}
		//TR << "Number of residues:" << pose.total_residue() << "\n";
		//TS_TRACE("some message");
	}

	void test_pdb_read_partial_residues() {
		core::pose::Pose pose;
		// This file has a leading fragment of an Arg residue that needs to be ignored.
		core::import_pose::pose_from_pdb(pose, "core/io/1ten.pdb");
		TS_ASSERT_EQUALS( pose.total_residue(), 89 );
		TR << pose.annotated_sequence() << std::endl;
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "L[LEU_p:NtermProteinFull]DAPSQIEVKDVTDTTALITWFKPLAEIDGIELTYGIKDVPGDRTTIDLTEDENQYSIGNLKPDTEYEVSLISRRGDMSSNPAKETFTT[THR_p:CtermProteinFull]" );
	}
};

