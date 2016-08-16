// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Rescore membrane protein test
/// @author Bjorn Wallner

// libRosetta headers


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedStubID.hh>

#include <core/fragment/JumpSRFD.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/MembraneTopology.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <basic/database/open.hh>
//#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <ObjexxFCL/format.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <numeric/xyzVector.io.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>


using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

//Auto Headers
#include <core/import_pose/import_pose.hh>


int
main( int argc, char* argv [] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );
	using namespace core;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// setup residue types
//	core::chemical::ResidueTypeSetCAP rsd_set;
//	if ( option[ in::file::fullatom ]() ) {
//		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
//	} else {
//		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
//	}

	// configure score function
//	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("score_membrane");

	// configure silent-file data object
//	core::io::silent::SilentFileData sfd;
	Size p1=16;
	Size p2=209;
	Size cut=(p1+p2)/2;
	core::pose::Pose start_pose;
	std::string infile  = *(option[ in::file::silent ]().begin());
	core::import_pose::pose_from_file( start_pose, infile, core::import_pose::PDB_file);
	core::pose::Pose pose(start_pose);
	Size nres=pose.total_residue();
	core::kinematics::FoldTree f(nres);
	f.new_jump(p1,p2,cut);
//	f.set_jump_atoms(1,"N","N",true);
	f.put_jump_stubs_intra_residue();
	std::cout <<  f << "\n";
	pose.fold_tree(f);
	std::cout <<  pose.fold_tree() << "\n";

	/*
	Vector const & tmp_n1( pose.residue( p1 ).atom( "N" ).xyz());
	Vector const & tmp_ca1( pose.residue( p1 ).atom( "CA" ).xyz());
	Vector const & tmp_c1( pose.residue( p1 ).atom( "C" ).xyz());
	Vector const & tmp_o1( pose.residue( p1 ).atom( "O" ).xyz());

	Vector const  & tmp_n2( pose.residue( p2 ).atom( "N" ).xyz());
	Vector const & tmp_ca2( pose.residue( p2 ).atom( "CA" ).xyz());
	Vector const & tmp_c2( pose.residue( p2 ).atom( "C" ).xyz());
	Vector const & tmp_o2( pose.residue( p2 ).atom( "O" ).xyz());
	core::kinematics::RT rt( core::kinematics::Stub( tmp_ca1, tmp_n1, tmp_ca1, tmp_c1 ), core::kinematics::Stub( tmp_ca2, tmp_n2, tmp_ca2, tmp_c2) );
	std::cout << "RT " << rt << std::endl;
    core::kinematics::Jump jump(rt);
	std::cout << "JUMP " << jump << std::endl;

	*/


	core::kinematics::Stub stub1(
		pose.residue( p1 ).atom( "CA" ).xyz(),
		pose.residue( p1 ).atom( "N" ).xyz(),
		pose.residue( p1 ).atom( "CA" ).xyz(),
		pose.residue( p1 ).atom( "C" ).xyz());
	core::kinematics::Stub stub2(
		pose.residue( p2 ).atom( "CA" ).xyz(),
		pose.residue( p2 ).atom( "N" ).xyz(),
		pose.residue( p2 ).atom( "CA" ).xyz(),
		pose.residue( p2 ).atom( "C" ).xyz());

	std::cout << stub1 << std::endl;
	std::cout << stub2 << std::endl;

	core::kinematics::RT rt(
	core::kinematics::Stub(
		pose.residue( p1 ).atom( "CA" ).xyz(),
		pose.residue( p1 ).atom( "N" ).xyz(),
		pose.residue( p1 ).atom( "CA" ).xyz(),
		pose.residue( p1 ).atom( "C" ).xyz()),
	core::kinematics::Stub(
		pose.residue( p2 ).atom( "CA" ).xyz(),
		pose.residue( p2 ).atom( "N" ).xyz(),
		pose.residue( p2 ).atom( "CA" ).xyz(),
		pose.residue( p2).atom( "C" ).xyz()));

	core::kinematics::RT rt2(
	core::kinematics::Stub(
		pose.residue( p1 ).atom( "N" ).xyz(),
		pose.residue( p1 ).atom( "CA" ).xyz(),
		pose.residue( p1 ).atom( "C" ).xyz()),
	core::kinematics::Stub(
		pose.residue( p2 ).atom( "N" ).xyz(),
		pose.residue( p2 ).atom( "CA" ).xyz(),
		pose.residue( p2).atom( "C" ).xyz()));


	std::cout << "JUMP " << rt2 << std::endl;
	std::cout << "RT " << rt << std::endl;

	id::StubID up_stub,down_stub;
	//core::id::NamedStubID name_stub1("CA","N","CA","C", p1 );

    up_stub = core::id::StubID( core::id::NamedStubID( "CA","N","CA","C", p1 ), pose );
    down_stub = core::id::StubID( core::id::NamedStubID( "CA","N","CA","C", p2 ), pose );
    pose.conformation().set_stub_transform( up_stub, down_stub,rt );


	//pose.set_jump( 1, jump );
	core::Real rmsd(core::scoring::CA_rmsd( start_pose, pose ));
	std::cout << "START RMS " << rmsd <<"\n";
//	pose.fold_tree(f);
// pose.set_jump( 1, jump2.reverse() );
	//core::Real rmsd2(core::scoring::CA_rmsd( start_pose, pose ));
	//std::cout << "START RMS " << rmsd2 <<"\n";

	pose.dump_pdb("tmp.pdb");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
