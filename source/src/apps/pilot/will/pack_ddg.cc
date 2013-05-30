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




#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoringManager.hh>
#include <basic/Tracer.hh>

#include <protocols/simple_moves/MinMover.hh>
//#include <protocols/moves/ResidueMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>
#include <core/scoring/packstat/compute_sasa.hh>

// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <map>
#include <ctime>

#include <core/import_pose/import_pose.hh>
using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace io::pdb;
using core::pose::Pose;
using utility::vector1;
using core::Real;

std::map<std::string,std::string>
get_resname1to3()
{
	std::map<std::string,std::string> m;
	m["G"] = "GLY";
	m["A"] = "ALA";
	m["L"] = "LEU";
	m["M"] = "MET";
	m["F"] = "PHE";
	m["W"] = "TRP";
	m["K"] = "LYS";
	m["Q"] = "GLN";
	m["E"] = "GLU";
	m["S"] = "SER";
	m["P"] = "PRO";
	m["V"] = "VAL";
	m["I"] = "ILE";
	m["C"] = "CYS";
	m["Y"] = "TYR";
	m["H"] = "HIS";
	m["R"] = "ARG";
	m["N"] = "ASN";
	m["D"] = "ASP";
	m["T"] = "THR";
	return m;
}



core::conformation::ResidueOP
get_residue( std::string name ) {
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	if( 1 == name.size() ) name = ((std::string)get_resname1to3()[name]);
	return ResidueFactory::create_residue( residue_set->name_map(name) );
}


void
refine_pose( Pose & pose, int seqpos = 0 ){
	using namespace core::scoring;

  ScoreFunctionOP sfstd( getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS ) );

	kinematics::MoveMapOP mm = new kinematics::MoveMap;
	mm->set_bb ( false ); mm->set_chi( true );	mm->set_jump( false );
	if ( 0 != seqpos ) {
		mm->set_bb( seqpos, true );
	}
	protocols::simple_moves::MinMover minstd( mm, sfstd , "dfpmin", 0.001, true, false, false );
	minstd.min_options()->nblist_auto_update(true);

	pack::task::PackerTaskOP taskstd = pack::task::TaskFactory::create_packer_task( pose );
	taskstd->restrict_to_repacking();
	taskstd->or_include_current(true);
	for( int i = 1; i <= (int)pose.total_residue(); ++i ) {
		taskstd->nonconst_residue_task(i).or_ex1(true);
		taskstd->nonconst_residue_task(i).or_ex2(true);
	}
	if( 0 != seqpos ) {
		taskstd->nonconst_residue_task(seqpos).or_ex3(true);
		taskstd->nonconst_residue_task(seqpos).or_ex4(true);
	}
	protocols::simple_moves::PackRotamersMover packstd( sfstd, taskstd );

	packstd.apply( pose );
	//minstd.apply( pose );

}

void
mutate_residue( Pose & pose, int const seqpos, std::string res_to ) {
	core::conformation::ResidueOP r = get_residue( res_to );
	pose.replace_residue( seqpos, *r, true );

}

///////////////////////////////////////////////////////////////////////////////
Real
packing_score( Pose const & pose, core::Size const rsd_id = 0 )
{
	if( 0 == rsd_id ) {
		Real ps = core::scoring::packstat::compute_packing_score( pose, 1 );
		return ps;
	} else {
		return core::scoring::packstat::compute_residue_packing_score( pose, rsd_id, 1 );
	}
}



void
test_ddg( Pose const orig, std::string pdb, int seqpos,
					std::string from_res, std::string to_res, Real ddg1, Real ddg2 ) {

	Pose pose = orig;

	if( get_resname1to3()[from_res] != pose.residue(seqpos).name3() ) {
		std::cerr << "bad from_res " << pdb << " " << seqpos << " " << from_res << " " << to_res << std::endl;
		// std::exit(-1);
	}

 	Real ps_pre = core::scoring::packstat::compute_residue_packing_scores( pose, 1 )[seqpos];
	// dump_pdb(pose,"before.pdb");

	refine_pose(pose,seqpos);

	std::cout << "DDG_PACKSTAT " << pdb    << " " << from_res << " " << seqpos << " "
						<< to_res << " " << ddg1     << " " << ddg2   << " before ";
	Real ps_before = packing_score(pose,seqpos);
	// dump_pdb(pose,"mid.pdb");

	pose = orig;
	mutate_residue( pose, seqpos, to_res );
	refine_pose(pose,seqpos);

	std::cout << " after ";
	Real ps_after = packing_score(pose,seqpos);
	std::cout << std::endl;
	// dump_pdb(pose,"after.pdb");
	//
	std::cerr << "DDG_PACKSCORE " << pdb
			<< " " << from_res
	 		<< " " << seqpos
			<< " " << to_res
	    << " ddg " << ddg1 << " " << ddg2 << " packing: " << ps_pre << " " << ps_before << " " << ps_after << std::endl;

}

void
read_ddg_file( std::string pdb, std::string fname ) {
	Pose orig;
	core::import_pose::pose_from_pdb( orig, pdb+".pdb" );
	using namespace std;
	ifstream in( fname.c_str() );
	char buf[999];
	while(true) {
		in.getline(buf,999);
		istringstream iss(buf);
		int num,seqpos;
		string from,to,tmpstr;
		Real ddg1,ddg2;
		if( !(iss >> num >> ddg1 >> ddg2 >> from >> seqpos >> to) )	break;
		if( from.size() > 1 ) {
			std::cerr << "error parsing line: " << pdb << " " << buf << std::endl;
			continue;
		}
		if( (iss >> tmpstr) ) {
			std::cerr << "can only handle one mutation: " << pdb << " " << buf << std::endl;
			continue; // can only handle one mutation
		}
		if( seqpos < 1 || seqpos > (int)orig.total_residue() ) {
			std::cerr << "seqpos out of range: " << pdb << " " << buf << std::endl;
			continue;
		}
		test_ddg( orig, pdb, seqpos, from, to, ddg1, ddg2 );
		// std::cerr << pdb << " " << seqpos << " " << from << " " << to << std::endl;

	}
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init(argc, argv);

	std::string d = "/Users/sheffler/project/packstat/ddg/training_set/";

	// test_ddg( "1bni", 87, "L", "T" );

	// read_ddg_file( "1bni", d+"1bni-ddG.txt" );
	// read_ddg_file( "1bvc", d+"1bvc-ddG.txt" );
	// read_ddg_file( "1fkj", d+"1fkj-ddG.txt" );
	// read_ddg_file( "1ftg", d+"1ftg-ddG.txt" );
	// read_ddg_file( "1shf", d+"1shf-ddG.txt" );
	// read_ddg_file( "1stn", d+"1stn-ddG.txt" );
	read_ddg_file( "1u5p", d+"1u5p-ddG.txt" );
	read_ddg_file( "2ci2", d+"2ci2-ddG.txt" );
	// read_ddg_file( "2lzm", d+"2lzm-ddG.txt" );
	// read_ddg_file( "5azu", d+"5azu-ddG.txt" );

	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

