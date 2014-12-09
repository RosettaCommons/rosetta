// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/io/pdb/file_data.cc
///
/// @brief
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)

#include <core/pose/annotated_sequence.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
//#include <core/chemical/residue_io.hh>

#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/Ramachandran.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>


#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>



#include <basic/database/open.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pose_io.hh>


#include <utility/stream_util.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <protocols/init/init.hh>
#include <core/init/init.hh>
#include <devel/init.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>


//silly using/typedef

using namespace core;
using namespace basic;

using utility::vector1;
using basic::T;
using basic::Error;
using basic::Warning;


typedef std::map< std::string, std::map< std::string, numeric::xyzVector< Real > > > Coords;

typedef vector1< std::string > Strings;

///////////////////////////////////////////////////////////////////////////////
// some silly helper routines:
//
///////////////////////////////////////////////////////////////////////////////

std::string readFile(std::string fname)
{
	Size fsize;
	std::string res;

	std::ifstream file(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
	if( file.is_open() )  {
		fsize = file.tellg();
		res.resize( fsize );
		file.seekg(0, std::ios::beg);
		file.read(&res[0], fsize);
		file.close();
	}
	else std::cout << "file not found!";
	return res;
}


// Pose Observer example
#include <core/pose/Pose.hh>
#include <core/pose/signals/GeneralEvent.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


class PosePyObserver : public utility::pointer::ReferenceCount
{
public:
	PosePyObserver() {};
	virtual ~PosePyObserver() {
		std::cout << "~PosePyObserver..." << std::endl;
	};


	void link(core::pose::Pose &p) {
		p.attach_general_obs(&PosePyObserver::generalEvent, this);
		//pose_ = p;
	}

	void unlink(core::pose::PoseOP &p) {
		p->detach_general_obs(&PosePyObserver::generalEvent, this);
		//pose_ = 0;
	}

	virtual void generalEvent( core::pose::signals::GeneralEvent const & ) {
		std::cout << "PosePyObserver::generalEvent... C++ version" << std::endl;
	};

private:
	core::pose::PoseOP pose_; // we want to keep OP to linked pose, so it cannot got deleted by Python
};


///////////////////////////////////////////////////////////////////////////////

int main( int argc, char * argv [] )
{

	try {


	using namespace core;

	//devel::init(argc, argv);
	//protocols::init::init(argc, argv);
	core::init::init(argc, argv);

	{
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb(pose, "src/python/bindings/test/data/test_in.pdb");
		//core::import_pose::pose_from_pdb(pose, "test_in.pdb");

		//core::scoring::ScoreFunction scorefxn;
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function_legacy( scoring::PRE_TALARIS_2013_STANDARD_WTS );
		T("Score:") << scorefxn->score(pose)  << std::endl;

		//scorefxn(pose);
		pose.energies().residue_total_energies(1);
		T("Scoring done!") << "---------------------" << std::endl;
    }
    {
		T("Testing pose_from_sequence...") << std::endl;
		std::string sequence(1000, 'V');
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, sequence, "fa_standard");
		T("Testing pose_from_sequence... Done!") << std::endl;

		return 0;
	}


	chemical::ResidueTypeSetCAP residue_set( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

	T("pdb_demo") << "---------------------" << std::endl;

	T("pdb_demo") << "\n\n\n\n\n>>>>> Reading file..." << std::endl;
	//std::string data = readFile("test_in.pdb");
	std::string data = readFile("src/python/bindings/test/data/test_in.pdb");

	T("pdb_demo") << ">>>>> Creating file data..." << std::endl;
	core::io::pdb::FileData fd = core::io::pdb::PDB_DReader::createFileData(data);

	T("pdb_demo") << "Remarks:" << std::endl;
	for(Size i=0; i<fd.remarks->size(); i++) {
		T("pdb_demo") << fd.remarks->at(i).num << " " << fd.remarks->at(i).value << std::endl;
	}


	T("pdb_demo") << ">>>>> Building Pose..." << std::endl;

	//T("pdb_demo") << fd;

	core::pose::Pose pose;
	core::import_pose::build_pose(fd, pose, *residue_set);
	T("pdb_demo") << " pose.total_residue()=" << pose.total_residue() << "\n";
	//T("pdb_demo") << fd;
	T("pdb_demo") << "Back to PDB now... -----------------" << std::endl;

	core::io::pdb::FileData fd1;
	fd1.init_from_pose(pose);

	//fd1 = fd;
	std::cout << io::pdb::PDB_DReader::createPDBData(fd);

	pose.dump_pdb("output.pdb");


	using namespace basic::options::OptionKeys;

	std::vector<int> A;  A.push_back(1);  A.push_back(2);  A.push_back(3);  A.push_back(5);
	utility::vector1<int> B;  B.push_back(10);  B.push_back(20);  B.push_back(30);  B.push_back(45);
	std::map<int, std::string> M;  M[1]="one";  M[2]="two";  M[3]="1+2";

	T("Demo") << "vector:" << A << " vector1:" << B << " map:" << M << "\n";

	T("Error", basic::t_error) << "Some error here!!!\n";
	T("core.pose") << "Some core pose message\n";
	T("core") << "Some core message\n";

	Error() << "Some error test...\n";
	Warning() << "Some warning test...\n";


	std::cout << "Scoring...\n";

	core::scoring::ScoreFunction scorefxn;
	scorefxn(pose);
	pose.energies().residue_total_energies(1);

	std::cout << "Creating PosePyObserver...\n";
	PosePyObserver PO;
	PO.link(pose);

	std::cout << "Scoring pose agagin...\n";
	scorefxn(pose);
	pose.energies().residue_total_energies(1);


	std::cout << "Done! -------------------------------\n";

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

