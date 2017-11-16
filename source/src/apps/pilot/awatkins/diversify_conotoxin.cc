// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <numeric/xyz.functions.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/ncbb/util.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/CyclizationMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/string_util.hh>

// C++ headers
#include <string>
#include <sstream>

using namespace core;
using namespace conformation;
using namespace core::chemical;
using namespace scoring;
using namespace constraints;
using namespace func;
using namespace pose;
using namespace protocols;
using namespace protocols::ncbb;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("DiversifyConotoxin");

utility::vector1< utility::vector1< Size > >
power_set( utility::vector1< Size > thiol_positions ) {

	utility::vector1< utility::vector1< Size > > ret;

	Size size = thiol_positions.size();
	Size pow_set_size = pow(2, size);
	Size i, j;

	for ( i = 1; i <= pow_set_size; ++i ) {
		utility::vector1< Size > x;
		for ( j = 1; j <= size; ++j ) {
			if ( (i-1) & ( 1<<(j-1)) ) x.push_back( thiol_positions[j] );
		}
		ret.push_back( x );
	}
	return ret;
}

int
main( int argc, char* argv[] )
{
	try {
		devel::init(argc, argv);

		scoring::ScoreFunctionOP score_fxn = scoring::get_score_function();

		score_fxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		score_fxn->set_weight( core::scoring::angle_constraint, 1.0 );
		score_fxn->set_weight( core::scoring::dihedral_constraint, 1.0 );//10.0 );

		Pose pose;
		std::string filename = option[in::file::s].value()[1];
		import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
		TR << "Importing pose from " << filename << std::endl;
		core::chemical::ResidueTypeSetCOP residue_set_cap = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		// hard coded for now, yuck
		utility::vector1< Size > thiol_positions;
		thiol_positions.push_back( 5 );
		thiol_positions.push_back( 6 );
		thiol_positions.push_back( 11 );
		thiol_positions.push_back( 19 );
		utility::vector1< Size > thiol_partners;
		thiol_partners.push_back( 11 );
		thiol_partners.push_back( 19 );
		thiol_partners.push_back( 5 );
		thiol_partners.push_back( 6 );

		utility::vector1<std::string> names;
		//names.push_back( "DCYS" );
		names.push_back( "B3C" );
		//names.push_back( "C41" );
		//names.push_back( "F41" );


		utility::vector1< utility::vector1< Size > > thiol_options = power_set( thiol_positions );
		for ( Size ii = 1; ii <= thiol_options.size(); ++ii ) {
			TR << "{ ";
			for ( Size jj = 1; jj <= thiol_options[ ii ].size(); ++jj ) {
				TR <<thiol_options[ ii ][jj] << " ";
			}
			TR << "} " << std::endl;
		}
		TR << std::endl;
		filename.erase( filename.size()-4 );

		for ( Size ti = 1; ti <= names.size(); ++ti ) {
			std::string thiol_nm = names[ti];

			// skip empty set
			for ( Size ii = 2; ii <= thiol_options.size(); ++ii ) {
				TR << "Doing " << thiol_options[ ii ] << std::endl;

				std::string outfile = filename;
				outfile += "_" + thiol_nm + "_";
				for ( Size jj = 1; jj <= thiol_options[ ii ].size(); ++jj ) {
					outfile += utility::to_string( thiol_options[ ii ][ jj ] );
				}
				outfile += ".pdb";


				Pose copy_pose;
				if ( thiol_nm == "B3C" ) {
					// Add every residue from pose, one by one.
					ResidueOP new_b3c = ResidueFactory::create_residue( residue_set_cap->name_map( "B3C" ) );//,
					ResidueOP new_b3cd = ResidueFactory::create_residue( residue_set_cap->name_map( "B3C:disulfide" ) );//,

					ResidueOP new_cys = ResidueFactory::create_residue( residue_set_cap->name_map( "CYS" ) );//,

					copy_pose.append_residue_by_jump( pose.residue(1), 1 );
					Size counter = 1;
					for ( Size ri = 2; ri <= pose.size(); ++ri ) {
						if ( counter <= thiol_options[ ii ].size() &&  ri == thiol_options[ ii ][counter] ) {
							copy_pose.append_residue_by_bond( *new_b3c, false );
							counter++;
						} else {
							copy_pose.append_residue_by_bond( pose.residue( ri ), false );
						}
					}
					copy_pose.conformation().detect_disulfides();

					// Now add special constraints to get bond working...
					for ( Size jj = 1; jj <= thiol_options[ ii ].size(); ++jj ) {
						Size resi = thiol_options[ ii ][ jj ];
						if ( resi > 1 ) {
							copy_pose.add_constraint(
								AtomPairConstraintOP( new AtomPairConstraint(
								*new AtomID( copy_pose.residue( resi-1 ).atom_index( "C" ), resi-1 ),
								*new AtomID( copy_pose.residue( resi   ).atom_index( "N" ), resi   ),
								TopOutFuncOP( new TopOutFunc( 100, 1.33, 5 ) ) ) ) );
						}
						if ( resi < copy_pose.size() ) {
							copy_pose.add_constraint(
								AtomPairConstraintOP( new AtomPairConstraint(
								*new AtomID( copy_pose.residue( resi+1 ).atom_index( "N" ), resi+1 ),
								*new AtomID( copy_pose.residue( resi   ).atom_index( "C" ), resi   ),
								TopOutFuncOP( new TopOutFunc( 100, 1.33, 5 ) ) ) ) );
						}
					}

					kinematics::MoveMapOP mm( new kinematics::MoveMap() );

					for ( Size kk = 1; kk <= copy_pose.size(); ++kk ) {
						if ( thiol_nm == "B3C" || copy_pose.residue( kk ).chain() == copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_options[ ii ][1] ) ).chain() ) {
							mm->set_bb( kk, true );
							mm->set_chi( kk, true );
						}
					}


					// create minimization mover
					simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, "lbfgs_armijo_nonmonotone"/*option[ OptionKeys::run::min_type ].value()*/, 0.01, true ) );
					for ( Real mi = 0.01; mi <= 2; mi += 0.2 ) {
						TR << "Minimization with apc weight " << mi << std::endl;
						score_fxn->set_weight( atom_pair_constraint, mi );
						minM->apply( copy_pose );
					}

					for ( Real mi = 2; mi >= 0; mi -= 0.4 ) {
						TR << "Minimization with apc weight " << mi << std::endl;
						score_fxn->set_weight( atom_pair_constraint, mi );
						minM->apply( copy_pose );
					}

					for ( Size jj = 1; jj <= thiol_options[ ii ].size(); ++jj ) {
						copy_pose.replace_residue( thiol_options[ ii ][ jj ], *new_b3cd, true );
					}

					copy_pose.conformation().detect_disulfides();

				} else {


					copy_pose = pose;

					TR << "Copied pose" << std::endl;


					// replace this thiol position
					for ( Size jj = 1; jj <= thiol_options[ ii ].size(); ++jj ) {
						std::string rn = copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_options[ ii ][ jj ] ) ).name();
						rn.erase( 0, rn.find(':') );
						rn = thiol_nm + rn;
						ResidueOP new_hcys = ResidueFactory::create_residue(
							residue_set_cap->name_map( rn/*"C26:disulfide"*/ ),
							copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_options[ ii ][ jj ] ) ),
							copy_pose.conformation() );
						copy_pose.replace_residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_options[ ii ][ jj ] ), *new_hcys, true );
					}
				}
				/*
				core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms(
				copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_positions[ ii ] ) ),
				*new_hcys,
				copy_pose.conformation() );
				*/
				//copy_pose.conformation().replace_residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_positions[ ii ] ), *new_hcys, true );
				TR << "Residue replaced" << std::endl;

				// add constraint
				if ( thiol_nm == "B3C" ) {
					for ( Size jj = 1; jj <= thiol_positions.size()/2; ++jj ) {
						copy_pose.add_constraint(
							AtomPairConstraintOP( new AtomPairConstraint(
							*new AtomID( copy_pose.residue( thiol_positions[ jj ] ).atom_index(copy_pose.residue( thiol_positions[ jj ] ).type().get_disulfide_atom_name()), thiol_positions[ jj ] ),
							*new AtomID( copy_pose.residue( thiol_partners[ jj ]  ).atom_index(copy_pose.residue( thiol_partners[ jj ]  ).type().get_disulfide_atom_name()), thiol_partners[ jj ]  ),
							//HarmonicFuncOP( new HarmonicFunc( 2.02, 0.02 ) ) ) ) );
							TopOutFuncOP( new TopOutFunc( 100, 2.02, 5 ) ) ) ) );
					}
				} else {

					for ( Size jj = 1; jj <= thiol_positions.size()/2; ++jj ) {
						copy_pose.add_constraint(
							AtomPairConstraintOP( new AtomPairConstraint(
							*new AtomID( copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_positions[ jj ] ) ).atom_index(copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_positions[ jj ] ) ).type().get_disulfide_atom_name()), thiol_positions[ jj ] ),
							*new AtomID( copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_partners[ jj ]  ) ).atom_index(copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_partners[ jj ]  ) ).type().get_disulfide_atom_name()), thiol_partners[ jj ]  ),
							//HarmonicFuncOP( new HarmonicFunc( 2.02, 0.02 ) ) ) ) );
							TopOutFuncOP( new TopOutFunc( 100, 2.02, 5 ) ) ) ) );
					}
				}
				TR << "constraints added" << std::endl;
				copy_pose.dump_pdb( outfile /*"out.pdb"*/ );

				// create move map for minimization
				kinematics::MoveMapOP mm( new kinematics::MoveMap() );
				//mm->set_bb( true );
				//mm->set_chi( true  );
				//mm->set_bb( false );
				for ( Size kk = 1; kk <= copy_pose.size(); ++kk ) {
					if ( thiol_nm == "B3C" || copy_pose.residue( kk ).chain() == copy_pose.residue( copy_pose.pdb_info()->pdb2pose( 'F', thiol_options[ ii ][1] ) ).chain() ) {
						mm->set_bb( kk, true );
						mm->set_chi( kk, true );
					}
				}
				//mm->set_chi( copy_pose.pdb_info()->pdb2pose( 'F', thiol_positions[ ii ] ), true );
				//mm->set_chi( copy_pose.pdb_info()->pdb2pose( 'F', thiol_partners[ ii ]  ), true );

				// create minimization mover
				simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, "lbfgs_armijo_nonmonotone"/*option[ OptionKeys::run::min_type ].value()*/, 0.01, true ) );
				for ( Real mi = 0.01; mi <= 2; mi += 0.2 ) {
					score_fxn->set_weight( atom_pair_constraint, mi );
					minM->apply( copy_pose );
				}

				for ( Real mi = 2; mi >= 0; mi -= 0.4 ) {
					score_fxn->set_weight( atom_pair_constraint, mi );
					minM->apply( copy_pose );
				}

				copy_pose.dump_pdb( outfile /*"out.pdb"*/ );
			}
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main

