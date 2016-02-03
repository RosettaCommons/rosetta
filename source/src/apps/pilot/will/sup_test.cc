// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/sup_test.cc
/// @brief test rms_util superimpose stuff

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>


int main (int argc, char *argv[])
{

	try {


	devel::init(argc,argv);

	core::pose::Pose mod_pose,ref_pose;

	core::import_pose::pose_from_file(mod_pose,basic::options::option[basic::options::OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);
	core::import_pose::pose_from_file(ref_pose,basic::options::option[basic::options::OptionKeys::in::file::s]()[2], core::import_pose::PDB_file);

	using namespace core::id;
	AtomID_Map<AtomID> atom_map;
	core::pose::initialize_atomid_map(atom_map,mod_pose,BOGUS_ATOM_ID);
	for(core::Size ir = 1; ir <= 6; ++ir) {
		for(core::Size ia = 1; ia <= 4; ia++) {
			core::Size ref_rsd = ref_pose.n_residue() - 6 + ir;
			atom_map[ AtomID(ia,ir) ] = AtomID(ia,ref_rsd);
		}
	}

	core::Real rms = core::scoring::superimpose_pose(mod_pose,ref_pose,atom_map);

	std::cout << "aligned region rms: " << rms << std::endl;

	mod_pose.dump_pdb("moved_pose.pdb");
	ref_pose.dump_pdb("fixed_pose.pdb");

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


