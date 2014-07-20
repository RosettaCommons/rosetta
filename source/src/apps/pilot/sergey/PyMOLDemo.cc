// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
///
/// @brief
/// @author Sergey Lyskov

#include <basic/Tracer.hh>

#include <protocols/moves/PyMolMover.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/PDB_Info.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <devel/init.hh>


//#include <unistd.h>

#include <utility/io/zipstream.hpp>

#include <utility/tools/make_map.hh>
#include <utility/excn/Exceptions.hh>

#include <stdio.h>
#include <cstring>

basic::Tracer TR("PyMOLDemo");


#include <unistd.h>

//Auto Headers
#include <core/import_pose/import_pose.hh>


int main( int argc, char * argv [] )
{

	try {

	using namespace core;

	devel::init(argc, argv);

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb(pose, "src/python/bindings/test/data/test_in.pdb");

	protocols::moves::AddPyMolObserver( pose );
	for(int j=0; j<32; j++) {
		pose.set_phi(70, pose.phi(70) + 1. );
		usleep(1000000);
	}

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function_legacy( scoring::PRE_TALARIS_2013_STANDARD_WTS );
	scorefxn->score(pose);
	//pose.energies().residue_total_energies(1);
	//T("Scoring done!") << "---------------------" << std::endl;

	protocols::moves::PyMolMover pymol;

	//pymol.keep_history(true);
	//pymol.update_interval(.1);

	core::Real a = 0.;

	pymol.print("Hi PyMOL!\n");

	for(int j=0; j<16; j++) {
		pose.set_phi(50, a);
		a += 1.;

		std::ostringstream msg;
		msg << "Tottal energy=" << scorefxn->score(pose);
		pymol.print(msg.str());

		//TR << "Sending pose..." << std::endl;
		pymol.apply(pose);
		//TR << "Sending energies..." << std::endl;
		pymol.send_energy(pose);

		usleep(1000000);
		//TR << a << std::endl ;
	}

	for(int j=0; j<1000; j++) {
		for(unsigned int r=1; r<=pose.total_residue(); r++) {
			std::map<int, int> C = utility::tools::make_map(int(r), int(protocols::moves::XC_white),
															int(1+pose.total_residue() - r), int(protocols::moves::XC_red) );
			pymol.send_colors(pose, C );
			usleep(1000000);
		}
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}



int _main( int argc, char * argv [] )
{
	using namespace core;

	devel::init(argc, argv);

	protocols::moves::UDPSocketClient s;

	std::string msg;
	for(int i=0; i<65536; i++) msg.push_back('_');

	msg = msg + "Qwe...";

	std::ostringstream ostringstream_;
	zlib_stream::zip_ostream zipper(ostringstream_, true);
	zipper << msg;
	//zipper.zflush();
	zipper.zflush_finalize();

	s.sendMessage(ostringstream_.str());
	// ___Qwe...') 65542


	/*
	sockaddr_in addr;
	memset(&addr, '\0', sizeof(sockaddr_in));


	addr.sin_family = AF_INET;       // host byte order
	addr.sin_port = htons(65000);     // short, network byte order
	//addr.sin_addr.s_addr = INADDR_ANY;
	addr.sin_addr.s_addr = inet_addr("127.0.0.1");

	int s_id = socket(AF_INET, SOCK_DGRAM, 0);

	char * buffer = "Some message here...\0";
	int error;// = sendto(s_id, buffer, strlen(buffer),0 , (struct sockaddr *)&addr, sizeof(struct sockaddr_in));

	if( error == -1 ) {
		printf("Cannot connect... Exiting.\n");
		return(1);
	}*/
	return 0;
}
