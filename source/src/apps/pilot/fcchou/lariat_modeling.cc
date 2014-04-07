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
#include <core/init/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/farna/RNA_SuiteAssign.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <protocols/farna/RNA_IdealCoord.hh>
#include <protocols/farna/RNA_ProtocolUtil.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>


#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <utility/excn/Exceptions.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

using namespace core;
using namespace core::pose;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;
using numeric::conversions::radians;
using numeric::conversions::degrees;

static numeric::random::RandomGenerator RG(5075);  // <- Magic number, do not change it!

utility::vector1 < utility::vector1 < Real > >
get_torsion_set( Pose const & pose )
{
	using namespace id;
	Size const total_res = pose.total_residue();
	utility::vector1 < utility::vector1 < Real > > torsion_set;
	utility::vector1 < Real > torsion;
	Size moving_suite;

	for (Size i = 1; i <= 3; ++i) {
		moving_suite = total_res - 4 + i;
		torsion.clear();
		torsion.push_back( pose.torsion( TorsionID( moving_suite, id::BB, 5 ) ) );
		torsion.push_back( pose.torsion( TorsionID( moving_suite, id::BB, 6 ) ) );
		torsion.push_back( pose.torsion( TorsionID( moving_suite+1, id::BB, 1 ) ) );
		torsion.push_back( pose.torsion( TorsionID( moving_suite+1, id::BB, 2 ) ) );
		torsion.push_back( pose.torsion( TorsionID( moving_suite+1, id::BB, 3 ) ) );
		torsion.push_back( pose.torsion( TorsionID( moving_suite+1, id::BB, 4 ) ) );
		torsion.push_back( pose.torsion( TorsionID( moving_suite+1, id::CHI, 1 ) ) );
		torsion.push_back( pose.torsion( TorsionID( moving_suite+1, id::CHI, 4 ) ) );
		torsion_set.push_back( torsion );
	}
	return torsion_set;
}

void
apply_torsion_set( Pose & pose, utility::vector1 < utility::vector1 < Real > > const & torsion_set )
{
	using namespace id;
	static protocols::farna::RNA_IdealCoord const ideal_coord_rna;
	Size const total_res = pose.total_residue();
	Size moving_suite;

	for (Size i = 1; i <= 3; ++i) {
		moving_suite = total_res - 4 + i;
		utility::vector1 < Real > const & torsion = torsion_set[i];
		pose.set_torsion( TorsionID( moving_suite, id::BB, 5 ), torsion[1] );
		pose.set_torsion( TorsionID( moving_suite, id::BB, 6 ), torsion[2]  );
		pose.set_torsion( TorsionID( moving_suite+1, id::BB, 1 ), torsion[3]  );
		pose.set_torsion( TorsionID( moving_suite+1, id::BB, 2 ), torsion[4]  );
		pose.set_torsion( TorsionID( moving_suite+1, id::BB, 3 ), torsion[5]  );
		pose.set_torsion( TorsionID( moving_suite+1, id::CHI, 1 ), torsion[7]  );
		pose.set_torsion( TorsionID( moving_suite+1, id::CHI, 4 ), torsion[8]  );
		if (torsion[6] < 115) {
			if ( pose.torsion( TorsionID( moving_suite+1, id::BB, 4 ) ) > 115) {
				ideal_coord_rna.apply(pose, moving_suite+1, true);
			}
		} else {
			if ( pose.torsion( TorsionID( moving_suite+1, id::BB, 4 ) ) < 115) {
				ideal_coord_rna.apply(pose, moving_suite+1, false);
			}
		}
	}
}

void
update_torsion_set( utility::vector1 < utility::vector1 < Real > > & torsion_set, Real const stdev )
{
	for (Size i = 1; i <= torsion_set.size(); ++i) {
		for (Size j = 1; j <= torsion_set[i].size(); ++j) {
			if (i == 1 && (j == 1 || j == 2) )  {
				continue;
			} else if ( j == 6 && RG.uniform() < 0.2) {
				torsion_set[i][j] = (RG.uniform() < 0.5) ? 82 : 130;
			} else {
				torsion_set[i][j] += RG.gaussian() * stdev;
			}
			if ( torsion_set[i][j] > 180 ) torsion_set[i][j] -= 360;
			if ( torsion_set[i][j] < -180 ) torsion_set[i][j] += 360;
		}
	}

}

bool
pose_list_compare( std::pair <Real, Pose> const & i, std::pair <Real, Pose> const & j)
{
	return (i.first < j.first);
}

void
lariat_modeling ()
{
	using namespace chemical;
	using namespace scoring;
	using namespace core::id;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace utility::io;
	using namespace scoring::constraints;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->
	          residue_type_set ( RNA );

	Pose pose;
	std::string pdb_name;
	if ( option[ in::file::native ].user() ) {
		pdb_name = option[in::file::native]();
		import_pose::pose_from_pdb ( pose, *rsd_set, pdb_name );
		protocols::farna::make_phosphate_nomenclature_matches_mini(pose);
	} else {
		utility_exit_with_message("User must specify -native option!");
	}

	//Setup score function.
	std::string score_weight_file = "rna_hires";
	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_weight_file= option[ basic::options::OptionKeys::score::weights ]();
		std::cout << "User passed in score:weight option: " << score_weight_file << std::endl;
	}
	core::scoring::ScoreFunctionOP scorefxn =
			ScoreFunctionFactory::create_score_function ( score_weight_file );

	protocols::farna::RNA_IdealCoord const ideal_coord_rna;
	Size const total_res = pose.total_residue();
	for (Size i = total_res - 2; i <= total_res; ++i) {
		ideal_coord_rna.apply( pose, i );
	}

	ConstraintSetOP cst_set = pose.constraint_set()->clone();
	Size const atm_indexO2prime = pose.residue(total_res).atom_index ( "O2'" );
	Vector const coord_P = pose.residue(total_res - 2).xyz("P");
	Vector const coord_O3prime = pose.residue(total_res - 3).xyz("O3'");
	Vector const PO3_bond_norm = (coord_O3prime - coord_P).normalize();
	Vector const coord_O2prime = coord_P - PO3_bond_norm * 3.7;
	cst_set->add_constraint( new CoordinateConstraint( AtomID(atm_indexO2prime, total_res),
																						AtomID(1, 1), coord_O2prime,
																						new HarmonicFunc( 0.0, 0.5 ) ) );
	pose.constraint_set ( cst_set );

	Real kT = 1.0;
	Size const n_steps_per_cycle = 10000;
	Real stdev = 3;
	utility::vector1 < utility::vector1 < Real > > torsion_set = get_torsion_set(pose);
	utility::vector1 < utility::vector1 < Real > > torsion_set_new = torsion_set;
	Real score, score_new;

	for (Size i = 0; i != 11; ++ i) {
		Real const cst_weight = i * 0.1;
		scorefxn -> set_weight( coordinate_constraint, cst_weight );
		score = (*scorefxn) (pose);
		Size n_accept = 0;
		for (Size j = 0; j != n_steps_per_cycle; ++j) {
			update_torsion_set( torsion_set_new, stdev );
			apply_torsion_set(pose, torsion_set_new);
			score_new = (*scorefxn) (pose);
			if (score_new < score || RG.uniform() < exp( (score - score_new) / kT) ){
				score = score_new;
				torsion_set = torsion_set_new;
				++n_accept;
			} else {
				torsion_set_new = torsion_set;
			}
		}
		std::cout << "Cycle " << i+1 << "/11 completed! accept rate = " << double(n_accept) / n_steps_per_cycle <<std::endl;
	}

	Size const n_step = 500000;
	Size const n_step_per_save = 10000;
	kT = 0.5;
	stdev = 1;
	Size n_step_after_save = 0;
	Size n_accept = 0;
	utility::vector1 < std::pair <Real, Pose> > pose_list;
	Pose lowest_pose;
	Real lowest_score = 99999;

	for (Size j = 0; j != n_step; ++j) {
		update_torsion_set( torsion_set_new, stdev );
		apply_torsion_set(pose, torsion_set_new);
		score_new = (*scorefxn) (pose);
		if (score_new < score || RG.uniform() < exp( (score - score_new) / kT) ){
			++n_accept;
			score = score_new;
			torsion_set = torsion_set_new;
			if (score < lowest_score) {
				lowest_score = score;
				lowest_pose = pose;
			}
		} else {
			torsion_set_new = torsion_set;
		}
		++n_step_after_save;
		if (n_step_after_save == n_step_per_save) {
			pose_list.push_back( std::pair <Real, Pose> (0.0, pose) );
			n_step_after_save = 0;
		}
	}
	std::cout << "accept rate = " << double(n_accept) / n_step << std::endl;

	pose_list.push_back( std::pair <Real, Pose> (0.0, lowest_pose) );

	AtomTreeMinimizer minimizer;
	float const dummy_tol ( 0.00000001 );
	MinimizerOptions min_options1 ( "dfpmin_armijo", dummy_tol, true, false, false );

	kinematics::MoveMap mm;
	mm.set_bb ( false );
	mm.set_chi ( false );
	mm.set_jump ( false );
	for (Size i = 1; i <= 3; ++i) {
		Size const moving_res = total_res - 3 + i;
		mm.set_bb ( moving_res, true );
		mm.set_chi ( moving_res, true );
	}

	for (Size i = 1; i <= pose_list.size(); ++i) {
		Pose & pose = pose_list[i].second;
		minimizer.run ( pose, mm, *scorefxn, min_options1 );
		Real const score = (*scorefxn) (pose);
		pose_list[i].first = score;
	}

	std::sort( pose_list.begin(), pose_list.end(), pose_list_compare );

	core::io::silent::SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	for (Size j = 1; j <= pose_list.size(); ++j) {
		std::ostringstream oss;
		oss << "decoy_" << j;
		core::io::silent::BinaryRNASilentStruct s( pose_list[j].second, oss.str() );
		silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );
	}

}

int
main( int argc, char * argv [] )
{
    try {
	using namespace core;

	/////////////////////////////
	// setup
	//////////////////////////////
	core::init::init ( argc, argv );
	//////////////////////////////
	// end of setup
	//////////////////////////////

	lariat_modeling();
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;

}
