
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

/// @file apps/pilot/smlewis/Extender.cc
/// @brief script to add residues to a PDB.  This script is NOT intended to make a good conformation - just the sequence-correct PDB for other modes (loop modeling, etc)
/// @author Steven Lewis

// Unit Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/chemical/ChemicalManager.hh> //CENTROID, FA_STANDARD
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

// Utility Headers
#include <devel/init.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>

// C++ headers
#include <string>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("apps.pilot.smlewis.Extender");

namespace basic{ namespace options{ namespace OptionKeys{
utility::options::FileOptionKey const extension("extension");
}}}//basic::options::OptionKeys


///@brief helper parses extension file
void read_in_extension(
											 core::pose::Pose const & pose,
											 core::Size & loop_start,
											 std::string & extension)
{
	std::string filename(basic::options::option[ basic::options::OptionKeys::extension ].value());

	//find the file, open it, handle error
	utility::io::izstream extension_stream( filename );
	if ( !extension_stream ) {
		Error() << "Can't open extension file, looked for: " << filename
						<< " use -extension <filename> to specify" << std::endl;
		utility_exit();
	}

	core::Size PDBloopstart;
	char chain;

	extension_stream >> chain >> PDBloopstart >> extension;

	//if the stream fails, bad formatting
	if( extension_stream.fail() ){
		Error() << "Can't parse extension file.  Using the PDB numbering, format is:\n    chain loopstart extension\n where extension is presented as CAPITAL SINGLE LETTER CODES out of ACDEFGHIKLMNPQRSTVWY" << std::endl;
		utility_exit();
	}

 	if (chain == '_') chain = ' ';

	core::pose::PDBPoseMap const & pose_map(pose.pdb_info()->pdb2pose());
	loop_start = pose_map.find(chain, PDBloopstart);

}//read_in_extension

int main( int argc, char* argv[] )
{

	try {

	using basic::options::option;
	using namespace basic::options::OptionKeys;
	option.add( extension, "extension file").def("extension");
	devel::init(argc, argv);

	//set up our pose
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, basic::options::start_file() );
	core::Size const poselength = pose.total_residue();

	//setup ultimate numberings in the combined pose
	core::Size loop_start;
	std::string extension;
	read_in_extension( pose, loop_start, extension);

	if( (loop_start > poselength) || (loop_start < 1) ){
		utility_exit_with_message("loop_start is not inside pose");
	}	else if (loop_start == 1) {
		utility_exit_with_message("N-terminal extension not supported yet");
	} else if (loop_start < poselength && loop_start != 1) {
		utility_exit_with_message("internal extension not supported yet");
	} else if (loop_start == poselength){
		TR << "C-terminal extension" << std::endl;
	} else {
		utility_exit_with_message("Whaaaat?");
	}

	core::chemical::ResidueTypeSetCAP typeset(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));

	core::Size const ext_length(extension.size());
	for(core::Size i(0); i<ext_length; ++i){
		core::Size insert_pos(loop_start+i);

		//name3 of new residue
		std::string name3(core::chemical::name_from_aa(core::chemical::aa_from_oneletter_code(extension[i])));

		//concrete new residue
		core::conformation::Residue newres(typeset->name_map(name3), true);

		pose.conformation().safely_append_polymer_residue_after_seqpos( newres, insert_pos, true );
		pose.set_phi(insert_pos, -150.0);
		pose.set_psi(insert_pos, 150.0);
		pose.set_omega(insert_pos, 180.0);

	} //finish inserting


	pose.dump_pdb("result.pdb");

	//internal loops ruminations
	//make loops fold tree; jump from  insert point-1 to insert point plus two
	//insert between insert point and insert point +1

	//feed loop into movemap, remember to leave omegas fixed at 180

	// 	protocols::loops::CcdLoopClosureMover close( interface_->loop(1), movemap );
	// 	close.apply(combined); //close the other gap
	// 	combined.dump_pdb("combined_nopretty.pdb");

	TR << "************************d**o**n**e**************************************" << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
