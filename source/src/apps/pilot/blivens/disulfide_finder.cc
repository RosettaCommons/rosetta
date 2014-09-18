
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file disulfide_finder.cc
/// @brief Gives some stats on disulfide bonds
/// @author Spencer Bliven
/// @date Created October 2008
/// @details
/// @section cli Command Line
/// @code disulfide_finder -l pdbfiles -o outfile -database db @endcode


//Utility
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <fstream>

//Core Chemistry
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

//Command line Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//disulfides
#include <apps/pilot/blivens/disulfides.hh>

using namespace core;
using namespace std;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;


//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>
static thread_local basic::Tracer TR( "pilot_apps.blivens.disulfide_finder" );


int
usage(char* msg)
{
	TR	<< "usage: disulfide_finder -l pdbfiles... -o outfile -database db" << endl
	<< msg << endl;
	exit(1);
}



int
main( int argc, char * argv [] )
{
  try {

	//Load a known disulfide bond as a reference

	utility::vector1<string> infiles;
	ofstream out;

	//init options system
	devel::init(argc, argv);

	if( option[ in::file::s ].user() || option[ in::file::l ].user() ) {
		infiles = basic::options::start_files();
	}
	else {
		return usage("No in file given: Use -s or -l to designate pdb files to search for disulfides");
	}

	if( option[ out::file::o ].user() ) {

		string outfile = option[ out::file::o ]();
		out.open(outfile.c_str());

	}
	else {
		return usage("No out file given: Use -o to designate an output file");
	}

	if( ! option[in::path::database].user() ) {
		return usage("No Database specified: Use -database to designate the location");
	}

	out << "file\ttrueDS\tpossibleDS\tmissedDS\ttrueNeg"<<endl;

	//for each input pdb
	for(vector1<string>::const_iterator infile_it = infiles.begin();	infile_it != infiles.end(); ++infile_it) {

		TR << "Analyzing "<<*infile_it<<endl;
		pose::Pose pose;
		core::import_pose::pose_from_pdb(pose, *infile_it);

		//Look for disulfide bonds
		int trueDS=0; // actual && possible
		int possibleDS=0; // possible but not actual
		int missedDS=0; // actual but not possible - nature must be wrong!
		int trueNeg=0; //for completeness, neither possible nor actual
		//for each pair of cysteines i,j
		for(Size i(1); i<= pose.total_residue()-1;i++) {
			//only consider cysteines
			if( pose.residue_type(i).name1() != 'C' )
				continue;

			for(Size j(i+1); j<= pose.total_residue();j++) {
				if( pose.residue_type(j).name1() != 'C' )
					continue;

				bool possible = possible_disulfide(pose, i, j);
				bool actual = actual_disulfide(pose, i, j);

				if( possible ) {
					if (actual) {
						TR << "Correctly identified a DS between "<<i<<" and "<<j<<endl;
						trueDS++;
					}
					else {
						TR << "Found a possible DS between "<<i<<" and "<<j<<endl;
						possibleDS++;
					}
				}
				else { //not possible
					if (actual) {
						TR << "Missed an actual DS between "<<i<<" and "<<j<<endl;
						missedDS++;
					}
					else {
						trueNeg++;
					}
				}
			} }//end cysteine pair

			TR << "Searched "<< trueDS+possibleDS+missedDS+trueNeg <<" residue pairs for possible DS bonds."<<endl;
			TR << "Correctly identified "<<trueDS<<" of "<<trueDS+missedDS<<" real DS."<<endl;
			TR << "Identified "<<possibleDS<<" additional possible DS."<<endl;

			out << *infile_it << '\t'
				<< trueDS << '\t' << possibleDS << '\t'
				<< missedDS << '\t' << trueNeg << endl;
	}//finished with this pdb

	out.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;
}//end main
