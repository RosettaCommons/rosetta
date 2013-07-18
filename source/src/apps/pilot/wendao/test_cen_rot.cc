// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file test_bbmc.cc
/// @brief test centroid rot model
/// @author Yuan Liu


// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>

#include <core/chemical/AA.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

// 
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
//#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

//
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>

#include <sstream>
#include <iostream>
#include <string>

using namespace core;
using namespace core::import_pose;
using namespace core::io::pdb;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core::pose;
using namespace core::chemical;
using namespace core::scoring;
using namespace core::kinematics;
using namespace protocols::moves;

// create a TaskFactory with the resfile
using namespace core::pack::task;
using namespace core::pack::dunbrack;
using namespace core::pack::dunbrack::cenrot;

static numeric::random::RandomGenerator RG(62331911);
basic::Tracer TR("pilot.wendao.cenrot");
std::map<std::string, Real> masslst;

void switch_to_residue_type_set_cenrot( core::pose::Pose & pose ) {
	//clear old energy cache
	pose.energies().clear();
	//get restype
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( "centroid_rot" );

	//loop for each residue
	for ( core::Size i=1; i<= pose.total_residue(); ++i ) {
		//get the residue
		core::conformation::Residue const & rsd( pose.residue(i) );
		//check current restype
		std::string const & current_type_set_name ( rsd.type().residue_type_set().name() );
		if ( current_type_set_name == "centroid_rot" ) {
			TR.Warning << "core::util::switch_to_residue_type_set: residue " << i 
			<< " already in centroid_rot residue_type_set" << std::endl;
			continue;
		}

		//gen new residue
		core::conformation::ResidueOP new_rsd( 0 );
		if( ( rsd.aa() == aa_unk ) ){
			//skip
		}
		else if ( rsd.name().substr(0,3)=="CYD" ) {
			core::chemical::ResidueType const & new_rsd_type( rsd_set->name_map("CYS") );
			new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
		}
		else  {
			core::chemical::ResidueType const & new_rsd_type( rsd_set->name_map(rsd.name().substr(0,3)) );
			new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );

			if ( ! new_rsd ) {
				TR.Warning << "Did not find perfect match for residue: "  << rsd.name()
				<< "at position " << i << ". Trying to find acceptable match. " << std::endl;
				core::chemical::ResidueTypeCOPs const & rsd_types( rsd_set->name3_map( rsd.name().substr(0,3) ) );
				for ( core::Size j=1; j<= rsd_types.size(); ++j ) {
					core::chemical::ResidueType const & new_rsd_type( *rsd_types[j] );
					if ( rsd.type().name3()  == new_rsd_type.name3()  ) {
						new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
						break;
					}
				}
				if (  new_rsd ) {
					TR.Warning << "Found an acceptable match: " << rsd.type().name() << " --> " << new_rsd->name() << std::endl;
				}
			}
		}

		//find the centroid postion
		PointPosition cenrotxyz(0,0,0);
		Real mass = 0.0;
		if (rsd.name3()=="GLY") {
			cenrotxyz = rsd.atom("CA").xyz();
		}
		else if (rsd.name3()=="ALA") {
			cenrotxyz = rsd.atom("CB").xyz();
		}
		else {
			//std::cout<<rsd.name()<<std::endl;
			for (	Size na=rsd.type().first_sidechain_atom(); 
					na<=rsd.type().nheavyatoms();
					na++) {
				if (rsd.atom_name(na)==" CB ") continue;
				std::string elem = rsd.atom_name(na).substr(1,1);
				cenrotxyz += (rsd.atoms()[na].xyz()*masslst[elem]);
				mass += masslst[elem];
				TR.Debug << "|" << rsd.atom_name(na) << "|"
				<< " (" << masslst[elem] << ") "
				<< rsd.atoms()[na].xyz().x() << ","
				<< rsd.atoms()[na].xyz().y() << ","
				<< rsd.atoms()[na].xyz().z() << std::endl;
			}
			cenrotxyz = cenrotxyz/mass;
		}

		//replace
		if ( ! new_rsd ) {
			std::cerr << pose.sequence() << std::endl;
			std::cerr  << "can not find a residue type that matches the residue " << rsd.name()
			<< "at position " << i << std::endl;
			//utility_exit_with_message( "switch_to_cenrot_residue_type_set fails\n" );
			continue;
		}

		//replace it
		pose.replace_residue( i, *new_rsd, false );
		//set centroid_rot xyz
		pose.set_xyz(id::AtomID(pose.residue(i).atom_index("CEN"), i), cenrotxyz);
	}
}

int main( int argc, char * argv [] ) {
	devel::init(argc, argv);
	masslst["C"]=12.0107;
	masslst["O"]=15.9994;
	masslst["S"]=32.066;
	masslst["N"]=14.00674;
	masslst["H"]=1.00794;
	
	ResidueTypeSetCAP rsd_set;
	rsd_set=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom("fa_standard");
	protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot("centroid_rot");

	Pose p;
	core::import_pose::pose_from_pdb( p, *rsd_set, option[ in::file::native ]() );
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();

	p.dump_pdb("out_fa.pdb");

	//switch
	//std::cout << "switch" << std::endl;
	if (0) {
		to_centroid.apply(p);
	} else {
		to_cenrot.apply(p);
	}

	p.dump_pdb("out_cenrot.pdb");

	to_fullatom.apply(p);
	p.dump_pdb("out_fa2.pdb");

	switch_to_residue_type_set_cenrot(p);
	p.dump_pdb("out_cenrot2.pdb");

	to_centroid.apply(p);
	p.dump_pdb("out_cen.pdb");

	to_cenrot.apply(p);
	p.dump_pdb("out_cenrot3.pdb");	

	//std::cout << p.residue(1).type().name() << std::endl;
	//std::cout << "E0= " << (*score_fxn)(p) << std::endl;
	//score_fxn->show(TR,p);
	//TR.flush();



	return 0;

	if (0) {
		std::cout << "min_init" << std::endl;
		protocols::simple_moves::MinMover minmover;
		minmover.score_function(*score_fxn);
		minmover.min_type("dfpmin");
		minmover.tolerance(1e-4);
		core::kinematics::MoveMapOP final_mm = new core::kinematics::MoveMap();
		//final_mm->set_chi( true );
		final_mm->set_bb( true );
		minmover.movemap(final_mm);

		std::cout << "min_apply" << std::endl;
		minmover.apply(p);
		
	}
	else {
		core::kinematics::MoveMap mm;
		mm.set_bb( true );
		//only for sidechain
		//mm.set_chi ( true );
		//mm.set( core::id::THETA, true );
		//mm.set( core::id::D, true );

		core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.00001, true, true, true );
		core::optimization::CartesianMinimizer minimizer;
		std::cout << "CART MINTEST: " << "\n";
		std::cout << "start score: " << (*score_fxn)(p) << "\n";
		long t1=clock();
		minimizer.run( p, mm, *score_fxn, options );
		long t2=clock();
		double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
		std::cout << "end score: " << (*score_fxn)(p) << "\n";
		std::cout << "MIN TIME: " << time << " sec \n";
	}

	(*score_fxn)(p);
	score_fxn->show(TR,p);
	TR.flush();
	p.dump_pdb("minimized.pdb");
	return 0;
}

