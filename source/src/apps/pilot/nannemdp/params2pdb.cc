// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file param2pdb.cc
/// @brief  read in a param file print out a PDB file.
/// @author Gordon Lemmon (glemmon@gmail.com)

//////////////////////////////////////////////////////////////////////////////

#include <utility/vector0.hh>

#include <devel/init.hh>
#include <core/types.hh>
//#include <core/options/option.hh>
//#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <protocols/rigid/RigidBodyMover.hh>

#include <numeric/xyz.functions.hh>
#include <utility/io/ozstream.hh>


int main(int argc, char* argv[])
{
	try{
	devel::init(argc, argv);

	core::chemical::ResidueTypeSelector rs;
	rs.set_property("PHOSPHONATE");
	core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
	core::chemical::ResidueTypeSetCAP rsd_set= cm->residue_type_set( core::chemical::FA_STANDARD );
	core::chemical::ResidueTypeCOPs fragment_types= rs.select( *rsd_set );

	core::chemical::ResidueTypeCOPs::const_iterator begin= fragment_types.begin();
	for(; begin!= fragment_types.end(); ++begin){
		core::conformation::Residue temp= core::conformation::Residue(**begin, true);
		std::string name= temp.name();
		name.append(".pdb");
		utility::io::ozstream out(name.c_str(), std::ios::out | std::ios::binary);
		//if(!out) {
		//	basic::Tracer::Error() << "FileData::dump_pdb: Unable to open file:" << name << " for writing!!!" << std::endl;
		//	return false;
		//}

		core::Size start=1;
		core::io::pdb::dump_pdb_residue(temp, start, out);
		out.close();
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}


}
