// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file sheet_finder.cc
/// @brief protocol to search list of PDBs output beta sheets with low rmsd to starting structure
/// @author Ben Stranges

// Unit Headers
// Project Headers
#include <basic/options/util.hh>
#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
// Auto-header: duplicate removed #include <basic/options/util.hh>

#include <core/id/AtomID_Map.hh>
#include <utility/file/FileName.hh>

#include <numeric/xyzVector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <utility/vector1.functions.hh>
#include <utility/vector1.hh>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <utility/io/izstream.hh>
#include <string>


// option key includes

//#include <basic/options/keys/align.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>


//namespaces
using namespace core;
using namespace core::scoring;
using namespace core::pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility;

// protocol specific options
namespace align {
FileOptionKey const master_pose_range( "align::master_pose_range");
FileOptionKey const list_pose_range( "align::list_pose_range" );
}

/// Reads in file that specifies which residues in master pose will be aligned
/// Each line of the file master_pose_range should have start and end position to align
bool
read_master_pose_range(
	std::string const & filename,
	utility::vector1< std::pair < Size, Size > > & master_ranges
)
{
	utility::io::izstream data( filename );
	if ( !data ) {
		std::cerr << " can not open master pose range file: " << filename << std::endl;
		return false;
	}

	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
		Size begin, end;
		line_stream >> begin >> end;
		if ( line_stream.fail() ) {
			std::cout << " can not parse line in master pose range file: " << line << std::endl;
			return false;
		}
		std::pair < Size, Size > range ( begin, end );
		master_ranges.push_back( range );
	}
	data.close();
	data.clear();
	return true;
}

/// Reads in file that specifies which residues in list pose will be aligned
/// Each line of the file list_pose_range should have start and end position to align
bool
read_list_pose_range(
	std::string const & filename,
	utility::vector1< std::pair < Size, Size > > & list_ranges
)
{
	utility::io::izstream data( filename );
	if ( !data ) {
		std::cerr << " can not open list pose range file: " << filename << std::endl;
		return false;
	}

	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
		Size begin, end;
		line_stream >> begin >> end;
		if ( line_stream.fail() ) {
			std::cout << " can not parse line in list pose range file: " << line << std::endl;
			return false;
		}
		std::pair < Size, Size > range ( begin, end );
		list_ranges.push_back( range );
	}
	data.close();
	data.clear();
	return true;
}

///Create task factory for rotamer trials mover
core::pack::task::TaskFactoryOP
generate_factory(){
	using namespace core::pack::task;
	TaskFactoryOP task_factory = new TaskFactory();

	//this assumes that the two protein partners are chains 1 and 2 - this is dangerous!!!!
	task_factory->push_back( new protocols::toolbox::task_operations::RestrictToInterfaceOperation(1, 2));
	task_factory->push_back( new operation::RestrictToRepacking );
	task_factory->push_back( new operation::InitializeFromCommandline() );
	//std::cout << "using default TaskFactory ( prevent repacking at metal site, detect interface" << std::endl;
	return task_factory;
}

///Main function////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char* argv[] )
{
	try{
		std::cout << "blah" << std::endl;
		using basic::options::option;
		option.add( align::master_pose_range, "file with master pose range").def("master_range");
		option.add( align::list_pose_range, "file with list pose range").def("list_range");
		devel::init(argc, argv);

		std::cout << " Begin " << std::endl;

		//-s read in "master" PDB
		core::pose::Pose master_pose;
		std::string pdbname(basic::options::option[ basic::options::OptionKeys::in::file::s ].value()[1]);
		core::import_pose::pose_from_file( master_pose, pdbname , core::import_pose::PDB_file);

		//-l read in list of pdbs
		core::pose::Pose list_pose;
		utility::vector1< std::string > pdbs( basic::options::start_files() );

		//-native for rmsd calculations
		core::pose::Pose nativepose;
		if ( option[ in::file::native ].active() ) {
			core::import_pose::pose_from_file( nativepose, basic::options::option[ basic::options::OptionKeys::in::file::native ]().name() , core::import_pose::PDB_file);
		}

		//PDBs read in through -l
		for ( core::Size i=1; i<= pdbs.size(); ++i ) {
			std::cout  << "PDBs from list: " << pdbs[i] << std::endl;
		}

		// read in master range
		std::string filename_master_range = option[ align::master_pose_range ]();
		utility::vector1< std::pair < Size, Size > >  master_ranges;
		read_master_pose_range( filename_master_range, master_ranges );

		//read in list range
		std::string filename_list_range = option[ align::list_pose_range ]();
		utility::vector1< std::pair < Size, Size > >  list_ranges;
		read_list_pose_range( filename_list_range, list_ranges );

		//find chain to cut out of master_pose later
		char pep_chain ('A');
		utility::vector1< core::Size >  res_to_lose;
		for ( core::Size i=1; i<= master_pose.size(); ++i ) {
			if ( master_pose.pdb_info()->chain(i) == pep_chain ) {
				res_to_lose.push_back(i);
			}
		}

		//make a scorefunction
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

		//initialize PackRotamersMover with defined task factory and score function
		protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
		pack_mover->task_factory( generate_factory() );
		pack_mover->score_function( scorefxn );

		core::Size numpdbs(pdbs.size());
		std::cout << "number of pdbs set "<< numpdbs << std::endl;

		//Begin loop over -l pdbs
		for ( core::Size j = 1; j < numpdbs; ++j ) {
			std::cout << "Begin loop " << j << std::endl;
			core::import_pose::pose_from_file( list_pose, pdbs[j] , core::import_pose::PDB_file);

			std::cout << "master PDB read in " << pdbname << std::endl;
			std::cout << "list PDB read in " << pdbs[j] << std::endl;

			//find RMSD between list_pose and native_pose
			float rmsd;
			if ( option[ in::file::native ].active() ) {
				rmsd = core::scoring::CA_rmsd( list_pose, nativepose );
				std::cout << "RMSD " << rmsd << " " <<pdbs[j] << std::endl;
			} else {
				rmsd = core::scoring::CA_rmsd( list_pose, list_pose );
			}

			//set atom map for superimpose_pose
			core::id::AtomID_Map< id::AtomID > atom_map;
			core::pose::initialize_atomid_map( atom_map, master_pose, core::id::AtomID::BOGUS_ATOM_ID() ); // maps every atomid to bogus atom

			utility::vector1< core::id::AtomID > ids1;
			utility::vector1< core::id::AtomID > ids2;
			std::pair< Size, Size > p1;
			std::pair< Size, Size > p2;
			p1 = master_ranges[ 1 ];
			p2 = list_ranges[ 1 ];

			std::cout << p1.first << " " << p1.second << " master_range" << std::endl;
			std::cout << p2.first << " " << p2.second << " list_range" << std::endl;

			for ( core::Size i=p1.first; i<= p1.second; ++i ) {
				core::id::AtomID dummy_atomid (master_pose.residue(i).atom_index("CA"), i);
				ids1.push_back(dummy_atomid);
			}

			for ( core::Size i=p2.first; i<= p2.second; ++i ) {
				core::id::AtomID dummy_atomid (list_pose.residue(i).atom_index("CA"), i);
				ids2.push_back(dummy_atomid);
			}
			assert( ids1.size()== ids2.size());
			std::cout << "ready to set atom map" << std::endl;

			for ( core::Size i=1; i<=ids1.size(); ++i ) {
				atom_map[ ids1[i] ] = ids2[i];
			}
			superimpose_pose( master_pose, list_pose, atom_map );
			std::cout << pdbname <<" superimposed onto "<< pdbs[j] << std::endl;

			//remove peptide from edit_master_pose
			core::pose::Pose edit_master_pose(master_pose);
			edit_master_pose.conformation().delete_residue_range_slow( min(res_to_lose), max(res_to_lose));
			core::pose::Pose full_pose (list_pose);
			std::cout << "deleting residues"<< " of " << pdbname << std::endl;

			//append edit_master_pose residues into full_pose
			full_pose.append_residue_by_jump(edit_master_pose.residue( 1 ), full_pose.size(), "" , "",  true /*start new chain*/);
			for ( core::Size n = 2; n<= edit_master_pose.size(); ++n ) {
				full_pose.append_residue_by_bond(edit_master_pose.residue ( n ));
			}

			std::cout << "rotamer repacking" << std::endl;
			pack_mover->apply(full_pose);


			// final output
			std::ostringstream outputfilename;
			outputfilename << "crib_" << pdbs[j];
			std::string fname = outputfilename.str();
			std::ofstream out( fname.c_str() );
			out << "RMSD: " << rmsd << '\n';
			full_pose.dump_scored_pdb( fname, *scorefxn );
			out.close();
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
