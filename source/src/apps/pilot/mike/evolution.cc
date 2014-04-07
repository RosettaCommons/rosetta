
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


// libRosetta headers

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

// Auto-header: duplicate removed #include <numeric/random/random.hh>

static numeric::random::RandomGenerator test_RG(1254); // <- Magic number, do not change it!!!

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/evolution.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace protocols;

using utility::vector1;
using io::pdb::dump_pdb;

struct PoseAndData {
public:
	void set_file_name( std::string filename_ ) {
		filename = filename_;
		std::cout << "reading: " << filename << std::endl;
		utility::vector1< std::string > token;

		token = utility::string_split( filename, '.' );
		if ( token.size() < 2 ) {
			utility_exit_with_message( "Invalid filename: " + std::string(filename) + " \n" );
		}
		clusternumber = utility::string2int( token[2] );
		if ( clusternumber < 0 ) {
			utility_exit_with_message( "Invalid filename: " + std::string(filename) + " \n" );
		}
	}

public:
	std::string filename;
	int silent_index;
	pose::Pose pose;
	Real energy;
	int clusternumber;
	Real rms;
};

bool compareEnergies( const PoseAndData &p1, const PoseAndData &p2 ) {
	return p1.energy < p2.energy;
}

void readPoseAndData_PDB(
	PoseAndData &pad,
	const std::string &filename,
	core::scoring::ScoreFunctionOP scorefxn,
	pose::PoseOP native_pose,
	bool havenative
) {

	pad.set_file_name( filename );
	core::import_pose::pose_from_pdb( pad.pose, pad.filename );
	pad.energy = (*scorefxn)(pad.pose);
	pad.rms = -1;
	pad.silent_index = -1;
	if ( havenative ) {
		pad.rms = protocols::simple_filters::native_CA_rmsd( *native_pose, pad.pose );
	}
}

void readPoseAndData_SILENT(
	PoseAndData &pad,
	const std::string &filename,
	core::io::silent::SilentFileData::iterator &data,
	core::scoring::ScoreFunctionOP scorefxn,
	pose::PoseOP native_pose,
	bool havenative
) {

	pad.set_file_name( filename );
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	data->fill_pose( pad.pose, *rsd_set );
	pad.energy = (*scorefxn)(pad.pose);
	pad.rms = -1;
	pad.silent_index = -1;
	if ( havenative ) {
		pad.rms = protocols::simple_filters::native_CA_rmsd( *native_pose, pad.pose );
	}
}


void fillPoseAndDataList(
	std::vector< utility::file::FileName > &list,
	std::vector< PoseAndData > &poses,
	core::scoring::ScoreFunctionOP scorefxn
) {

	using basic::options::option;
	using namespace basic::options::OptionKeys;

	bool havenative = false;
	pose::PoseOP native_pose = new pose::Pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() ); // default is standard fullatom residue_set
		havenative = true;
	}

	using namespace utility::file;
	for (std::vector< FileName >::iterator i = list.begin(), i_end = list.end(); i != i_end; ++i) {
		std::string filename( i->name() );
		std::ifstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while( getline(data, line) ) {
			PoseAndData pad;
			readPoseAndData_PDB( pad, FileName(line), scorefxn, native_pose, havenative);
			poses.push_back(pad);
		}
		data.close();
	}

}



void processChild(
	PoseAndData &child,
	std::vector< PoseAndData > &parent_list,
	int size_limit
) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	Real rms_limit = option[ evolution::rms_threshold ]();
	Real rms_topmargin = option[ evolution::rms_topmargin ]();

	std::cout << "Processing Child:  " << child.filename << "    " << child.clusternumber << "  "
		<< child.energy << "  " << child.rms << std::endl;

	// now process the child.


	//	int bestcluster = -1;
	int nclusters=0;
	std::vector< PoseAndData >::iterator bestparent = parent_list.begin();
	Real bestrms = 1000000;
	for (std::vector< PoseAndData >::iterator parent = parent_list.begin(), i_end = parent_list.end(); parent != i_end; ++parent) {
		Real rms = core::scoring::CA_rmsd( child.pose, parent->pose );
		if ( parent->clusternumber > nclusters) nclusters = parent->clusternumber;
		if ( rms < bestrms ) {
			bestrms = rms;
			bestparent = parent;
		}


	}

	std::cout << "BestRMS: " << parent_list.size() << "   " << bestrms << "   " << rms_limit << "  " << rms_topmargin << "  "
		<< bestparent->clusternumber << "  " << child.clusternumber << std::endl;

	if ( bestrms > rms_topmargin ) return;

	if ( bestrms > rms_limit ) {
		std::cout << "Adding new group" << std::endl;

		// find highest energy in *original* cluster assigned.
		bool status = true;
		int clustersize=0;
		std::vector< PoseAndData >::iterator highenergy  = parent_list.begin();
		for (std::vector< PoseAndData >::iterator parent = parent_list.begin(), i_end = parent_list.end(); parent != i_end; ++parent) {
			if ( parent->clusternumber != child.clusternumber ) continue;
			clustersize++;
			if ( status || ( parent->energy > highenergy->energy ) ) {
				highenergy = parent;
				status = false;
			};
		}

		if (clustersize > 1 )
		if ( highenergy->energy > child.energy ) {
			// delete that structure !

			if ( (int)parent_list.size() > size_limit )
			 parent_list.erase( highenergy );

			// add structure as new cluster
			child.clusternumber = nclusters+1;
			parent_list.push_back( child );
			nclusters++;
		}

	} else {
		// add to original cluster

		// find highest energy member of the best cluster

		bool status = true;
		std::vector< PoseAndData >::iterator highenergy  = parent_list.begin();
		for (std::vector< PoseAndData >::iterator parent = parent_list.begin(), i_end = parent_list.end(); parent != i_end; ++parent) {
			//if ( parent->clusternumber != bestparent->clusternumber ) continue;
			if ( status || ( parent->energy > highenergy->energy ) ) {
				highenergy = parent;
				status = false;
			};
		}

		int saveclusternumber = highenergy->clusternumber;

		// Only if new guy actually has better energy then parent to be removed ! (otherwise groups with only 1 parent get screwed over!!)
		if ( child.energy < highenergy->energy ) {
			(*highenergy)=child;
			highenergy->clusternumber = saveclusternumber;
		}
	}
} // processChild

void processChildren(
	std::vector< utility::file::FileName > &list,
	core::scoring::ScoreFunctionOP scorefxn,
	std::vector< PoseAndData > &parent_list,
	int size_limit
) {

	using basic::options::option;
	using namespace basic::options::OptionKeys;

	bool havenative = false;
	pose::PoseOP native_pose = new pose::Pose;
	if ( option[ in::file::native ].active() ) {
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() ); // default is standard fullatom residue_set
		havenative = true;
	}

	using namespace utility::file;
	for (std::vector< FileName >::iterator i = list.begin(), i_end = list.end(); i != i_end; ++i) {
		std::string filename( i->name() );
		std::ifstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while( getline(data, line) ) {
			PoseAndData child;

			//Decide if it's a PDB or a Silent File

			std::string filename  = FileName(line);
			std::string extention = filename.substr( filename.length() - 3, 3);

			if ( extention == "pdb") {
				std::cout << filename << "  " << extention << "  " << "PDB !" << std::endl;
				readPoseAndData_PDB( child, FileName(line), scorefxn, native_pose, havenative);
				processChild(child, parent_list, size_limit );
			} else if ( extention == "out") {
				std::cout << filename << "  " << extention << "  " << "OUTFILE !" << std::endl;
				core::io::silent::SilentFileData sfd;
				sfd.read_file( filename );
				int count = 0;
				for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
					readPoseAndData_SILENT( child, FileName(line), iter, scorefxn, native_pose, havenative);
					processChild(child, parent_list, size_limit );
					child.silent_index = count;
					count++;
				}
			} else {
				 utility_exit_with_message( "Unable to open file - unknown format (neither .pdb nor .out): " + filename + '\n' );
			}
		} // getline
		data.close();
	} // for filenames
} // processChildren

struct EnergyAndFilename {
	Real energy;
	std::string filename;
	int clusternumber;
};

/// @brief Takes a list of filenames that yield input poses, scores the poses using
/// the provided ScoreFunction,
void processChildrenIntensification(
	std::vector< utility::file::FileName > & list,
	core::scoring::ScoreFunctionOP scorefxn,
	std::vector< PoseAndData > & parent_list,
	int size_limit = 200
) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;

	// make a copy pof the parent list
	std::vector< PoseAndData > full_list( parent_list );

	bool havenative = false;
	pose::PoseOP native_pose = new pose::Pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() ); // default is standard fullatom residue_set
		havenative = true;
	}

	// blitz all the poses
	for (std::vector< PoseAndData >::iterator structure = full_list.begin(),
			i_end = full_list.end(); structure != i_end; ++structure) {
		structure->pose = pose::Pose();
	}

	for (std::vector< FileName >::iterator i = list.begin(), i_end = list.end(); i != i_end; ++i) {
		std::string filename( i->name() );
		std::ifstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while( getline(data, line) ) {
			PoseAndData child;
			std::string filename =  FileName(line);
			std::string extention = filename.substr( filename.length() - 3, 3);
			if ( extention == "pdb") {
				std::cout << filename << "  " << extention << "  " << "PDB !" << std::endl;
				readPoseAndData_PDB( child, FileName(line), scorefxn, native_pose, havenative);
				std::cout << child.energy << "  " << child.rms << std::endl;
				child.pose = pose::Pose(); // zap pose
				full_list.push_back( child );
			} else
			if ( extention == "out") {
				std::cout << filename << "  " << extention << "  " << "OUTFILE !" << std::endl;
				core::io::silent::SilentFileData sfd;
				sfd.read_file( filename );
				int count = 0;
				for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
					readPoseAndData_SILENT( child, FileName(line), iter, scorefxn, native_pose, havenative);
					child.silent_index = count;
					std::cout << child.energy << "  " << child.rms << std::endl;
					child.pose = pose::Pose(); // zap pose
					full_list.push_back( child );
					count++;
				}
			} else {
				 utility_exit_with_message( "Unable to open file - unknown format (neither .pdb nor .out): " + filename + '\n' );
			}
		}
		data.close();
	}

	// sorting by energy
	std::sort( full_list.begin(), full_list.end(), compareEnergies );

	parent_list.clear();

	int counter = 0;
	for (std::vector< PoseAndData >::iterator structure = full_list.begin(),
			i_end = full_list.end(); structure != i_end; ++structure) {
		if ( counter >= size_limit ) continue;
		parent_list.push_back( *structure );
		counter ++;
	}

	// re-read in surviving parents
	for (std::vector< PoseAndData >::iterator parent = parent_list.begin(), i_end = parent_list.end(); parent != i_end; ++parent) {
		std::string filename =  parent->filename;
		std::string extention = filename.substr( filename.length() - 3, 3);
		if ( extention == "pdb" ) {
			std::cout << filename << "  " << extention << "  " << "PDB !" << std::endl;
			readPoseAndData_PDB( *parent, parent->filename, scorefxn, native_pose, havenative  );
		} else if ( extention == "out") {
			std::cout << filename << "  " << extention << "  " << "OUTFILE !" << std::endl;
			core::io::silent::SilentFileData sfd;
			sfd.read_file( filename );
			int count = 0;
			for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
				if ( parent->silent_index == count) {
					readPoseAndData_SILENT(  *parent,  parent->filename, iter, scorefxn, native_pose, havenative);
					parent->silent_index = count;
				}
				count++;
			}
		} else {
			 utility_exit_with_message( "Unable to open file - unknown format (neither .pdb nor .out): " + filename + '\n' );
		}
	} // for parent_list
} // void processChildrenIntensification

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
    try {
	devel::init(argc, argv);

	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace scoring;
	using namespace utility::file;
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = getScoreFunction();

	// get parent list
	std::vector< FileName > parent_list;
	if ( option[ evolution::parentlist ].user() ) {
		parent_list = option[ evolution::parentlist ]().vector(); // make a copy (-l)
	}	else {
		utility_exit_with_message( " -evolution::parentlist must be given ! \n");
	}

	// get child list
	std::vector< FileName > child_list;
	if ( option[ evolution::childlist ].user() ) {
		child_list = option[ evolution::childlist]().vector(); // make a copy (-l)
	}	else {
		utility_exit_with_message( " -evolution::childlist must be given ! \n");
	}

	std::vector< PoseAndData > parent_poses;
	fillPoseAndDataList( parent_list, parent_poses, scorefxn );

	int c=0;

	// decide if this round we diversify ir intensify
	if ( option[ evolution::action ]() == "diversify" ) {
//		processChildrenIntensification( child_list, scorefxn, parent_poses, int( 0.5* float( child_list.size() + parent_poses.size() ) ));
		processChildren( child_list, scorefxn, parent_poses, 200);
	} else if ( option[ evolution::action ]() == "intensify" ) {
		processChildrenIntensification( child_list, scorefxn, parent_poses, 200 );
	} else {
		utility_exit_with_message( " -evolution::action must one of ( diversify | intensify ) \n");
	}


	// print final population!

	c=0;
	for (std::vector< PoseAndData >::iterator i = parent_poses.begin(), i_end = parent_poses.end(); i != i_end; ++i) {
		c++;
		std::cout << "Parent:  " << c << "  " << i->filename << "    " << i->silent_index << "  "
			<< i->clusternumber << "  "
			<< i->energy << "  "
			<< i->rms <<  std::endl;
	}


	// write final population i.e. print out PDBs

	std::string targetdir =  option[ evolution::targetdir ]();
	int last_clusternumber=0;
	for (std::vector< PoseAndData >::iterator i = parent_poses.begin(), i_end = parent_poses.end(); i != i_end; ++i) {
		if ( i->clusternumber > last_clusternumber ) last_clusternumber = i->clusternumber;
	}

	std::vector < int > count_per_cluster ( last_clusternumber+1 ,0   );
	std::cout << "SIZE: " << count_per_cluster.size();

	std::string scorefile ( targetdir + "/initial_scorerms");
	std::ofstream outfile ( scorefile.c_str() );

	Real worst_energy = -100000000.0;
	for (std::vector< PoseAndData >::iterator i = parent_poses.begin(), i_end = parent_poses.end(); i != i_end; ++i) {
		i->filename = targetdir + "/"
		                       "tag." + right_string_of( i->clusternumber, 3, '0' ) + "." +
		                       right_string_of( count_per_cluster[ i->clusternumber ], 3, '0'  ) + ".pdb";
		i->pose.dump_pdb( i->filename );
		count_per_cluster[ i->clusternumber ]++;
		outfile << i->rms << "  " << i->energy << std::endl;
		if ( i->energy > worst_energy ) worst_energy = i->energy;
	}

	std::string filter_settings ( targetdir + "/filter_settings");
	std::ofstream ffile ( filter_settings.c_str() );

	Real padding_score_filter  =  option[ evolution::padding_score_filter ]();
	Real padding_stage2_filter =  option[ evolution::padding_stage2_filter ]();

	ffile << " -looprelax::final_score_filter "
				<< worst_energy + padding_score_filter << "\n"
	      << " -relax::filter_stage2_quarter	"
				<< worst_energy + padding_stage2_filter << "\n"
				<< std::endl;
	ffile.close();

	outfile.close();
	std::cout << "Normal termination." << std::endl;
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}






