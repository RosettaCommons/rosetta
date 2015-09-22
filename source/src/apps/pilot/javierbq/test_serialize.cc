// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/dssp/Dssp.hh>


// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/PDBInfo.hh>

#include <string>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "serialization_test" );

class ThisApplication {
    public:
        ThisApplication(){};
        static void register_options();
		};

		OPT_KEY( String , pdb)

		void ThisApplication::register_options() {
				NEW_OPT( pdb, "pdb file", "");
		}


std::string
pose_to_string(const core::pose::Pose& pose){
	std::ostringstream ss;
	core::io::silent::SilentFileData sfd;
	core::io::silent::BinarySilentStruct sstruct;
	sstruct.fill_struct(pose, pose.pdb_info()->name());
	sstruct.print_header(ss);
	sfd.write_silent_struct(sstruct, ss);
	return ss.str();
}

core::pose::Pose
string_to_pose(const std::string& silent_string){
	std::istringstream ss;
	ss.str(silent_string);
	core::io::silent::SilentFileData sfd;
	core::pose::Pose pose;

	utility::vector1<std::string> readall;
	sfd.read_stream(ss, readall,false);
	sfd.get_structure(sfd.begin()->decoy_tag()).fill_pose(pose);
	return pose;
}

class SerializablePose: public core::pose::Pose {
	public:
		SerializablePose(const core::pose::Pose& p)
					: Pose(p)
			{ }

		SerializablePose();

		virtual core::pose::PoseOP clone() const {
				return new SerializablePose( *this );
		}
};

typedef utility::pointer::owning_ptr< SerializablePose > SerializablePoseOP;
typedef utility::pointer::owning_ptr< SerializablePose const > SerializablePoseCOP;

int
main ( int argc, char* argv[] ){
	try {

		ThisApplication::register_options();
		devel::init(argc, argv);
		TR << "Reaading File " << option[pdb]() << std::endl;
		//core::pose::PoseOP p = import_pose::pose_from_pdb( option[pdb]());
		core::pose::PoseOP pop = import_pose::pose_from_pdb( option[pdb]());
		SerializablePose p(*pop);


		TR << "size of pose to flatten " << p.total_residue() << std::endl;

		TR << "ss1:" << p.secstruct() << std::endl;
		p.set_secstruct(2,'H');
		TR << "ss2:" << p.secstruct() << std::endl;
		std::string buffer = pose_to_string(p);

		core::pose::Pose p2 = string_to_pose(buffer);
		TR << "size of pose to unflatten " << p2.total_residue() << std::endl;
		TR << "ss:" << p2.secstruct() << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

    return 0;
}
