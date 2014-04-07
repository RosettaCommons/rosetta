// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <protocols/rigid/RigidBodyMover.hh>

#include <numeric/xyz.functions.hh>
#include <utility/io/ozstream.hh>


int main(int argc, char* argv[])
{
	try{
	devel::init(argc, argv);

	core::chemical::ResidueSelector rs;
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
