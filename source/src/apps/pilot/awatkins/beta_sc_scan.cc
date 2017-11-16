// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <devel/init.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

#include <unistd.h>

/*#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path_traits.hpp>
#include <boost/filesystem/config.hpp>*/

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace id;
using namespace chemical;
using namespace conformation;
using namespace scoring;
using namespace func;
using namespace constraints;
using namespace pose;
using namespace protocols::simple_moves;

// tracer - used to replace cout
static basic::Tracer TR("BetaScScan");

void set_up_movemap( std::string const & n3, core::pose::Pose const & pose,  kinematics::MoveMapOP mm ) {
	// Everything can move except chi and bb torsions
	for ( Size ii = 1; ii <= pose.residue( 1 ).natoms(); ++ii ) {
		using namespace core::id;

		mm->set( DOF_ID( AtomID( ii, 1 ), D ), true );
		mm->set( DOF_ID( AtomID( ii, 1 ), THETA ), true );
		//mm->set( DOF_ID( AtomID( ii, 1 ), PHI ), true );
	}
	//mm->set_chi( false );
	//mm->set_bb( false );
	//mm->show();

	// OK, set chi 2 true bc it's a proton chi!
	if ( n3 == "B3S" ) {
		mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "HG" ), 1 ), PHI ), true );
		mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "OG" ), 1 ), PHI ), false );
	} else if ( n3 == "B3I" ) {
		mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "CG1" ), 1 ), PHI ), false );
		mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "CD1" ), 1 ), PHI ), false );
	} else if ( n3 == "B3L" ) {
		mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "CG" ), 1 ), PHI ), false );
		mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "CD1" ), 1 ), PHI ), false );
	} else if ( n3 == "B3W" ) {
		mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "CG" ), 1 ), PHI ), false );
		mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "CD1" ), 1 ), PHI ), false );
	}

	mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "CO" ), 1 ), PHI ), false );
	mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "C" ), 1 ), PHI ), false );
	mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "NM" ), 1 ), PHI ), false );
}

void scan_chi( core::pose::Pose & pose, core::scoring::ScoreFunctionOP score_fxn, kinematics::MoveMapOP mm, std::string const & ext ) {

	protocols::simple_moves::MinMoverOP min( new protocols::simple_moves::MinMover( mm, score_fxn, "linmin_iterated", 0.001, true ) );//lbfgs_armijo_nonmonotone", 0.001, true ) );


	utility::vector1< utility::vector1< Real > > results;
	if ( pose.residue_type( 1 ).name3() == "B3S" ) {
		for ( Real chi1 = 0; chi1 <= 350; chi1 += 10 ) {
			utility::vector1< Real > row;
			pose.conformation().set_torsion( TorsionID( 1, id::CHI, 1 ), chi1 );

			// score the pose
			//Real orig_ener( ( *score_fxn )( pose ) );

			min->apply( pose );

			std::stringstream ss;
			ss << "pose_" << chi1 << "_" /*<< chi2 << "_" */<< ext << ".pdb";
			pose.dump_pdb( ss.str() );
			// score the pose
			Real min_ener( pose.energies().total_energy() );
			row.push_back(min_ener);
			results.push_back( row );
		}

		for ( Size ii = 1; ii <= results.size(); ++ii ) {
			TR << (Real(ii * 10 - 10) ) << "\t";
			for ( Size jj = 1; jj <= results[ii].size(); ++jj ) {
				TR << results[ii][jj] << "\t";
			}
			TR << "\n";
		}
		TR << std::endl;
	} else {
		for ( Real chi1 = 0; chi1 <= 350; chi1 += 10 ) {
			utility::vector1< Real > row;
			for ( Real chi2 = 0; chi2 <= 350; chi2 += 10 ) {
				pose.conformation().set_torsion( TorsionID( 1, id::CHI, 1 ), chi1 );
				pose.conformation().set_torsion( TorsionID( 1, id::CHI, 2 ), chi2 );

				// score the pose
				//Real orig_ener( ( *score_fxn )( pose ) );

				min->apply( pose );

				std::stringstream ss;
				ss << "pose_" << chi1 << "_" << chi2 << "_" << ext << ".pdb";
				pose.dump_pdb( ss.str() );
				// score the pose
				Real min_ener( pose.energies().total_energy() );
				row.push_back(min_ener);
			}
			results.push_back( row );
		}

		TR << "\t";
		for ( Real chi2 = 0; chi2 <= 350; chi2 += 10 ) {
			TR << chi2 << "\t";
		}
		TR << "\n";
		for ( Size ii = 1; ii <= results.size(); ++ii ) {
			TR << (Real(ii * 10 - 10) ) << "\t";
			for ( Size jj = 1; jj <= results[ii].size(); ++jj ) {
				TR << results[ii][jj] << "\t";
			}
			TR << "\n";
		}
		TR << std::endl;
	}
}

void score_dir( std::string const & dir, core::scoring::ScoreFunctionOP score_fxn ) {

	TR << "Scoring all files in the dir " << dir << "..." << std::endl;
	PoseOP pose( new Pose );
	utility::vector1< utility::vector1< Real > > results;

	for ( Real chi1 = 0; chi1 <= 350; chi1 += 10 ) {
		utility::vector1< Real > row;
		TR << "Chi1 " << chi1 << "... " << std::endl;
		for ( Real chi2 = 0; chi2 <= 350; chi2 += 10 ) {

			std::stringstream fns;
			// This is just to transpose: we need to have chi1 columns and chi2 rows in final.
			fns << dir + "/B3W_" << chi2 << "_" << chi1 << ".pdb";

			try {
				int res = access(fns.str().c_str(), R_OK);
				if ( res < 0 ) {
					row.push_back( 100 );
					continue;
				}
				pose = import_pose::pose_from_file( fns.str(), import_pose::PDB_file );
				if ( chi1 == 0 && chi2 == 0 ) {
					pose->dump_pdb( "test.pdb" );
				}
				row.push_back( ( *score_fxn )( *pose ) );
			} catch ( ... ) {
				row.push_back( 100000 );
			}
		}
		results.push_back( row );
	}

	TR << "\t";
	for ( Real chi2 = 0; chi2 <= 350; chi2 += 10 ) {
		TR << chi2 << "\t";
	}
	TR << "\n";
	for ( Size ii = 1; ii <= results.size(); ++ii ) {
		TR << (Real(ii * 10 - 10) ) << "\t";
		for ( Size jj = 1; jj <= results[ii].size(); ++jj ) {
			TR << results[ii][jj] << "\t";
		}
		TR << "\n";
	}
	TR << std::endl;

	// Can't use this part of boost, aww
	/*using namespace boost::filesystem;

	directory_iterator end_itr;
	for ( directory_iterator itr( dir ); itr != end_itr; ++itr ) {
	pose = import_pose::pose_from_file( itr->path().string(), import_pose::PDB_file );
	TR << itr->path().filename() << ( *score_fxn )( *pose ) << std::endl;
	}*/
}

void score_dir_single( std::string const & dir, core::scoring::ScoreFunctionOP score_fxn ) {

	TR << "Scoring all files in the dir " << dir << "..." << std::endl;
	PoseOP pose( new Pose );
	utility::vector1< Real > results;

	for ( Real chi1 = 0; chi1 <= 350; chi1 += 10 ) {

		std::stringstream fns;
		// This is just to transpose: we need to have chi1 columns and chi2 rows in final.
		fns << dir + "/B3S_" << chi1 << ".pdb";

		try {
			int res = access(fns.str().c_str(), R_OK);
			if ( res < 0 ) {
				results.push_back( 100 );
				continue;
			}
			pose = import_pose::pose_from_file( fns.str(), import_pose::PDB_file );
			if ( chi1 == 0 ) {
				pose->dump_pdb( "test.pdb" );
			}
			results.push_back( ( *score_fxn )( *pose ) );
		} catch ( ... ) {
			results.push_back( 100000 );
		}
	}

	TR << "\n";
	for ( Size ii = 1; ii <= results.size(); ++ii ) {
		TR << (Real(ii * 10 - 10) ) << "\t" << results[ii] << "\n";
	}
	TR << std::endl;
}

int
main( int argc, char* argv[] )
{
	try {
		devel::init(argc, argv);


		// create score function
		scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );

		/*
		std::string const n3 = "B3W";
		chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		Pose pose;
		pose.append_residue_by_jump( Residue( restype_set->name_map( n3+":AcetylatedNtermProteinFull:MethylatedCtermProteinFull" ), true ), 1 );

		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		set_up_movemap( n3, pose, mm );

		// Extended

		pose.conformation().set_torsion( TorsionID( 1, id::BB, 1 ), -140 );
		pose.conformation().set_torsion( TorsionID( 1, id::BB, 2 ),   75 );
		pose.conformation().set_torsion( TorsionID( 1, id::BB, 3 ),   65 );

		scan_chi( pose, score_fxn, mm, "ext" );

		// Helical

		pose.conformation().set_torsion( TorsionID( 1, id::BB, 1 ), -140 );
		pose.conformation().set_torsion( TorsionID( 1, id::BB, 2 ),   60 );
		pose.conformation().set_torsion( TorsionID( 1, id::BB, 3 ), -120 );

		scan_chi( pose, score_fxn, mm, "hel" );
		*/
		score_dir_single( "renamed", score_fxn );
		//score_dir( "renamed", score_fxn );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
