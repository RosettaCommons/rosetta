// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilot/chrisk/rotamer_analysis.cc
/// @brief

//core library
#include <math.h>
#include <stdlib.h>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/chemical/util.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/types.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/sasa.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/tools/make_map.hh>
#include <utility/string_util.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rtmin.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>

//protocols library (Movers)
#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

//calculator stuff
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>
#include <basic/MetricValue.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/NumberWaterHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/SemiExplicitWaterUnsatisfiedPolarsCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/ExplicitWaterUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/SemiExplicitWaterUnsatisfiedHBondsCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/ExplicitWaterUnsatisfiedHBondsCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedHBondsCalculator.hh>

//utilities

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rot_anl.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

//local options
namespace basic{ namespace options{ namespace OptionKeys{
}}}//basic::options::OptionKeys

////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace basic;
using namespace id;
using namespace pose;
using namespace pack;
using namespace conformation;
using namespace chemical;
using namespace kinematics;
using namespace scoring;
using namespace options;
using namespace basic::options::OptionKeys;
using namespace optimization;
namespace OK = OptionKeys;
using utility::vector1;
using std::string;
using import_pose::pose_from_pdb;
using io::pdb::dump_pdb; // deprecated though
using namespace ObjexxFCL;
using basic::T;
using basic::Warning;
using basic::Error;

static thread_local basic::Tracer TR( "chrisk.byres_data" );
static thread_local basic::Tracer TR_unsat( "chrisk.unsat_calc" );


//local options
namespace byres_data
{
  basic::options::StringOptionKey tag( "byres_data:tag" );
//  basic::options::BooleanOptionKey calc_sasa( "byres_data:calc_sasa" );
  basic::options::BooleanOptionKey solvate( "byres_data:solvate" );
  basic::options::BooleanOptionKey calcs_only( "byres_data:calcs_only" );
  basic::options::BooleanOptionKey solv_unsat_calc( "byres_data:solv_unsat_calc" );
  basic::options::RealOptionKey repeat( "byres_data:repeat" );
  basic::options::IntegerOptionKey nloop_solvadd( "byres_data:nloop_solvadd" );
  basic::options::IntegerOptionKey nloop_solvdock( "byres_data:nloop_solvdock" );
  basic::options::IntegerOptionKey nloop_hbscan( "byres_data:nloop_hbscan" );
  basic::options::BooleanOptionKey align_native_seq( "byres_data:align_native_seq" );
}

/*
//need special function for water, has nor base2
//assuming O=1, H=4,5 for TP5 water
Stub
build_water_frame(
	Residue const & rsd,
	Size const & atm )
{
	//if atm is O
	Size base( 4 ), base2( 5 );
	//else
	if( atm == 4 ) base = 1;
	else if( atm == 5 ) base2 = 1;

	Vector const & Hxyz( rsd.atom( atm ).xyz() );
	Vector const & Dxyz( rsd.atom( base ).xyz() );
	Vector const & DBxyz( rsd.atom( base2 ).xyz() );
	return Stub(Hxyz Hxyz, Hxyz - Dxyz, Dxyz - DBxyz);
}

//return acc_base --> hydrogen angle close to ideal
// *actual* ideal angles would involve a db lookup, but we just wanna get in the neighborhood
Real
get_ideal_BAH_angle(
	Residue const & rsd,
	Size const & aatm)
{
	Hybridization const & acc_hybrid( rsd.atom_type(aatm).hybridization());

	switch(acc_hybrid){
		case SP2_HYBRID:
			return 122.6;
			break;
		case SP3_HYBRID:
			return 122.2;
			break;
		case RING_HYBRID:
			return 0.0;
			break;
		default:
			TR << "No ideal acc angle. Unrecognized acceptor hybridization" << std::endl;
	}
	return 120.0;
}

//build hbond frames
//graciously copied from momeara
Stub
build_donor_frame(
	Residue const & rsd,
	Size const & hatm )
{
	Vector const & Hxyz( rsd.atom( hatm ).xyz() );
	Vector const & Dxyz( rsd.atom(rsd.atom_base( hatm )).xyz());
	Vector const & DBxyz( rsd.atom(rsd.atom_base(rsd.atom_base(hatm))).xyz());
	return Stub(Hxyz Hxyz, Dxyz, DBxyz);
}

Stub
build_acceptor_frame(
	Residue const & rsd,
	Size const & aatm)
{
	Vector const & Axyz( rsd.atom( aatm ).xyz() );
	Vector const & ABxyz( rsd.atom( rsd.atom_base( aatm ) ).xyz() );
	Vector const & AB2xyz( rsd.atom(rsd.abase2(aatm)).xyz());
	Hybridization const & acc_hybrid( rsd.atom_type(aatm).hybridization());

	switch(acc_hybrid){
		case SP2_HYBRID:
			return Stub(Axyz, Axyz, ABxyz, AB2xyz );
			break;
		case SP3_HYBRID:
			//using hybrid type stub for sp3 because want to scan both sides
			//return Stub(Axyz, Axyz, AB2xyz, ABxyz );
			return Stub(Axyz, Axyz, Real(0.5)*(ABxyz+AB2xyz), ABxyz);
			break;
		case RING_HYBRID:
			return Stub(Axyz, Axyz, Real(0.5)*(ABxyz+AB2xyz), ABxyz);
			break;
		default:
			TR << "Cannot build acceptor frame: Unrecognized acceptor hybridization" << std::endl;
	}
	return Stub();
}

void
append_rsd_by_kinematic_hbond_jump_near_atom(
	pose::Pose & pose,
	Size seqpos,
	Size atomno,
	conformation::Residue new_rsd,
	Size new_atomno,
	Real dist_min,
	Real dist_max
)
{
	typedef  numeric::xyzMatrix< Real > Matrix;
	using namespace scoring::hbonds;

	Residue rsd( pose.residue( seqpos ) );
	//append by jump from seqpos atomno to new_rsd atom 1, maybe make random downstream atom?
	pose.append_residue_by_jump( new_rsd, seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), true );
	Size new_seqpos( pose.total_residue() );
	Size jump_number( pose.fold_tree().num_jump() );
	Jump jump( pose.jump( jump_number ) );

	//which is acceptor, donor?
	Size aatm( 0 ), acc_pos( 0 ), hatm( 0 ), don_pos( 0 );
	if( rsd.atom_type( atomno ).is_polar_hydrogen() &&
		new_rsd.atom_type( new_atomno ).is_acceptor() ){
		don_pos = seqpos;
		hatm = atomno;
		acc_pos = new_seqpos;
		aatm = new_atomno;
	}
	else if( new_rsd.atom_type( new_atomno ).is_polar_hydrogen() &&
		rsd.atom_type( atomno ).is_acceptor() ){
		don_pos = new_seqpos;
		hatm = new_atomno;
		acc_pos = seqpos;
		aatm = atomno;
	}
	else{ utility_exit_with_message( "ERROR: res " + string_of( seqpos ) + " atom " + string_of( atomno ) +
		" res " + string_of( new_seqpos ) + " atom " + string_of( new_atomno ) + " is not HB don/acc pair!!\n" );
	}
	//now get their base atoms to get datm and batm
	Size datm( pose.residue( don_pos ).atom_base( hatm ) );
	Size batm( pose.residue( acc_pos ).atom_base( aatm ) );

	//random distance within window
	//must import from hbonds/constants.hh
	Size steps( option[ byres_data::nloop_hbscan ] );
	Real AHdist_min(MIN_R), AHdist_max(MAX_R); Size AHdist_steps(steps);
	Real cosBAH_min(MIN_xH), cosBAH_max(MAX_xH); Size cosBAH_steps(steps);
	Real cosAHD_min(MIN_xD), cosAHD_max(MAX_xD); Size cosAHD_steps(steps);
	core::Real const MIN_xC = 0.;
	core::Real const MAX_xC = { numeric::constants::d::pi_2 }; // chi cutoff
	Real chi_min(MIN_xC), chi_max(MAX_xC); Size chi_steps(steps);

//debug
//	for( Real AHdist = AHdist_min; AHdist <= AHdist_max;
//					AHdist += (AHdist_max - AHdist_min)/static_cast<Real>(AHdist_steps-1)){
	for( Real cosBAH = cosBAH_min; cosBAH <= cosBAH_max;
					cosBAH += (cosBAH_max - cosBAH_min)/static_cast<Real>(cosBAH_steps-1)){
			for( Real cosAHD = cosAHD_min - .0001; cosAHD <= cosAHD_max;
							cosAHD += (cosAHD_max - cosAHD_min)/static_cast<Real>(cosAHD_steps-1)){
				for( Real chi = chi_min - .0001; chi <= chi_max;
								chi += (chi_max - chi_min)/static_cast<Real>(chi_steps-1)){

						//BAH angle depends on hybridization
	//					Real cosBAH( cos( get_ideal_BAH_angle( pose.residue( acc_pos ), aatm ) ) ); //angle from acc_base -> hydrogen
						Real AHdist( 2.2 );

						//trigonomixxx str8up stolen from momeara
						// sqrt( 1 - cos^2 ) = sin
						Vector Hxyz( cosBAH * AHdist,
								sin( chi ) * sqrt( 1 - cosBAH * cosBAH ) * AHdist,
								cos( chi ) * sqrt( 1 - cosBAH * cosBAH ) * AHdist );

						Real ADdist, cosBAD;
						Vector Dxyz;
						Stub acc_frame;

						// triangle AHD the AHD angle is
						ADdist = sqrt(AHdist*AHdist + 1 +2*AHdist*cosAHD);
						if ( cosBAH < .5 ) {
								cosBAD = (cosBAH*AHdist + cosBAH*cosAHD + sqrt(1-cosBAH*cosBAH)*sqrt(1-cosAHD*cosAHD))/ADdist;

						} else {
								cosBAD = (cosBAH*AHdist + cosBAH*cosAHD - sqrt(1-cosBAH*cosBAH)*sqrt(1-cosAHD*cosAHD))/ADdist;
						}
						Dxyz = Vector(cosBAD*ADdist,
								sin(chi)*sqrt(1-cosBAD*cosBAD)*ADdist,
								cos(chi)*sqrt(1-cosBAD*cosBAD)*ADdist);
	//debug

						TR << "hbond geometry:" << "\t";
						TR << "\tAHdist: " << AHdist << "\t";
						TR << "\tcosBAH: " << cosBAH << "\t";
						TR << "\tchi:    " << chi << "\t";
						TR << "\tcosAHD: " << cosAHD << "\t";
						TR << "\tHxyz:   " << "(" << Hxyz.x() << "," << Hxyz.y() << "," << Hxyz.z() << ")" << "\t";
						TR << "\tcosBAD: " << cosBAD << "\t";
						TR << "\tDxyz: " << Dxyz.x() << "," << Dxyz.y() << "," << Dxyz.z() << ")" << std::endl;



						//get the ideal hbond RT
						Stub don_frame( Hxyz, Hxyz, Dxyz, Real( 2 ) * Hxyz);

						if( don_pos == seqpos ){
							RT hbond_geometry( acc_frame, don_frame );

							//get atom ids
							AtomID aatm_id( aatm, acc_pos );
							AtomID hatm_id( hatm, don_pos );

							//get stubs for don and h
							AtomTree const & at = pose.conformation().atom_tree();
							Stub don_AT_frame( at.atom( hatm_id ).get_stub());
							Stub acc_AT_frame( at.atom( aatm_id ).get_stub());


							Stub don_hbond_frame;
							Stub acc_hbond_frame;

							don_hbond_frame = build_donor_frame( pose.residue( don_pos ), hatm );
							acc_hbond_frame = build_acceptor_frame( pose.residue( acc_pos ), aatm );

							RT don_AT_to_hbond( don_AT_frame, don_hbond_frame );
							RT acc_hbond_to_AT( acc_hbond_frame, acc_AT_frame );

							// Make hbond jump
							Stub a,b,c;
							//set stub a relative to don_AT_frame
							don_AT_to_hbond.make_jump( don_AT_frame, a );
							//set stub b relative to a
							hbond_geometry.make_jump( a, b );
							//define jump from actual don frame
							Jump hbond_jump( don_AT_frame, b);
							pose.set_jump( jump_number, hbond_jump );
						}
						else{
							RT hbond_geometry( acc_frame, don_frame );

							//get atom ids
							AtomID aatm_id( aatm, acc_pos );
							AtomID hatm_id( hatm, don_pos );

							//get stubs for acc and h
							AtomTree const & at = pose.conformation().atom_tree();
							Stub acc_AT_frame( at.atom( aatm_id ).get_stub());
							Stub don_AT_frame( at.atom( hatm_id ).get_stub());


							Stub acc_hbond_frame;
							Stub don_hbond_frame;

							acc_hbond_frame = build_acceptor_frame( pose.residue( acc_pos ), aatm );
							don_hbond_frame = build_donor_frame( pose.residue( don_pos ), hatm );

							RT acc_AT_to_hbond( acc_AT_frame, acc_hbond_frame );
							RT don_hbond_to_AT( don_hbond_frame, don_AT_frame );

							// Make hbond jump
							Stub a,b,c;
							//set stub a relative to acc_AT_frame
							acc_AT_to_hbond.make_jump( acc_AT_frame, a );
							//set stub b relative to a
							hbond_geometry.make_jump( a, b );
							//define jump from actual acc frame
							Jump hbond_jump(acc_AT_frame, b);
							pose.set_jump( jump_number, hbond_jump );
						}
					//debug
					//pose.dump_pdb( "watest." + string_of( atomno ) + "." + string_of( cosAHD ) + "." + string_of( cosBAH ) + "." + string_of( chi ) + ".pdb" );
				}
			}
	}
}

void
append_rsd_by_hbond_jump_near_atom(
	pose::Pose & pose,
	Size seqpos,
	Size atomno,
	conformation::Residue new_rsd,
	Size new_atomno,
	Real dist_min,
	Real dist_max
)
{
	typedef  numeric::xyzMatrix< Real > Matrix;

	Residue rsd( pose.residue( seqpos ) );
	//append by jump from seqpos atomno to new_rsd atom 1, maybe make random downstream atom?
	pose.append_residue_by_jump( new_rsd, seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), true );
	Size new_seqpos( pose.total_residue() );
	Size jump_number( pose.fold_tree().num_jump() );
	Jump jump( pose.jump( jump_number ) );

	//which is acceptor, donor?
	Size aatm( 0 ), acc_pos( 0 ), hatm( 0 ), don_pos( 0 );
	if( rsd.atom_type( atomno ).is_polar_hydrogen() &&
		new_rsd.atom_type( new_atomno ).is_acceptor() ){
		don_pos = seqpos;
		hatm = atomno;
		acc_pos = new_seqpos;
		aatm = new_atomno;
	}
	else if( new_rsd.atom_type( new_atomno ).is_polar_hydrogen() &&
		rsd.atom_type( atomno ).is_acceptor() ){
		don_pos = new_seqpos;
		hatm = new_atomno;
		acc_pos = seqpos;
		aatm = atomno;
	}
	else{ utility_exit_with_message( "ERROR: res " + string_of( seqpos ) + " atom " + string_of( atomno ) +
		" res " + string_of( new_seqpos ) + " atom " + string_of( new_atomno ) + " is not HB don/acc pair!!\n" );
	}
	//now get their base atoms to get datm and batm
	Size datm( pose.residue( don_pos ).atom_base( hatm ) );
	Size batm( pose.residue( acc_pos ).atom_base( aatm ) );

	//random distance within window
	Real AH_dist( dist_min + numeric::random::rg().uniform() * ( dist_max - dist_min ) );
	Real AHD_ang( 180.0 ); //angle from acc -> donor
	//BAH angle depends on hybridization
	Real BAH_ang( get_ideal_BAH_angle( pose.residue( acc_pos ), aatm ) ); //angle from acc_base -> hydrogen
	//pick a random chi angle
	Real BAHD_chi( numeric::random::rg().uniform() * 360.0 );

	//set translation
	jump.random_trans( AH_dist );
	pose.set_jump( jump_number, jump );
pose.dump_pdb( "ah.pdb" );
//	FoldTree f_jump( pose.fold_tree() );
	FoldTree f_jump( pose.total_residue() );
	FoldTree f_rot( pose.total_residue() );

	//mrrow?
	f_jump.new_jump( seqpos, new_seqpos, pose.total_residue() - 1 );
	f_jump.reorder( seqpos );
	f_jump.set_jump_atoms( jump_number, seqpos, rsd.atom_name( atomno ), new_seqpos, new_rsd.atom_name( new_atomno ) );
	pose.fold_tree( f_jump );

	//new cutpoint at end of chain
//	Size new_cutpoint( pose.conformation().chain_end( pose.chain( seqpos ) ) );
	//if we're already at end of chain, swap jump to first res
//	if( new_cutpoint == seqpos ){
//		f_rot.new_jump( 1, new_seqpos, seqpos );
//		new_cutpoint = 2;
//	}

	f_rot.new_chemical_bond( seqpos, new_seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), pose.total_residue() - 1 );
//	f_rot.reorder( seqpos );
	pose.fold_tree( f_rot );
	//to set rotation
	pose.conformation().set_bond_angle( AtomID( aatm, acc_pos ),
		AtomID( hatm, don_pos  ),
		AtomID( datm, don_pos  ),
		AHD_ang );
pose.dump_pdb( "ahd.pdb" );
	pose.conformation().set_bond_angle( AtomID( batm, acc_pos ),
		AtomID( aatm, acc_pos  ),
		AtomID( hatm, don_pos  ),
		BAH_ang );
pose.dump_pdb( "bah.pdb" );
	pose.conformation().set_torsion_angle( AtomID( batm, acc_pos ),
		AtomID( aatm, acc_pos  ),
		AtomID( hatm, don_pos  ),
		AtomID( datm, don_pos  ),
		BAH_ang );
pose.dump_pdb( "tor.pdb" );
	//now flip back to jump fold tree
	pose.fold_tree( f_jump );


//	pose.dump_pdb( "qqq." + string_of( seqpos ) + "." + string_of( atomno ) + ".pdb" );

}
*/

//check for clashes based on atom distances
//only check atoms in AtomID vector
//should be way faster than calculating entire score
bool
fast_clash_check(
	Pose const & pose,
	vector1< id::AtomID > const check_atids,
	Real const clash_dist_cut
)
{
	Real const clash_dist2_cut( clash_dist_cut * clash_dist_cut );
	for( Size iatid = 1; iatid <= check_atids.size(); ++iatid ){
		Vector const at1_xyz( pose.xyz( check_atids[ iatid ] ) );
		for( Size res2 = 1; res2 <= pose.total_residue(); ++res2 ){
			for( Size at2 = 1; at2 <= pose.residue( res2 ).natoms(); ++at2 ){
				//skip virtual atoms!
				if( pose.residue( res2 ).atom_type( at2 ).lj_wdepth() == 0.0 ) continue;
				id::AtomID atid2( at2, res2 );
				//skip if atid2 is in check_atids
				bool skip_at2( false );
				for( Size jatid = 1; jatid <= check_atids.size(); ++jatid ){
					if( atid2 == check_atids[ jatid ] ){ skip_at2 = true; break; }
				}
				if( skip_at2 ) continue;
				Real const dist2( at1_xyz.distance_squared( pose.xyz( atid2 ) ) );
				if( dist2 < clash_dist2_cut ){
					//TR_unsat << "CLASH!: " << check_atids[ iatid ] << " - " << atid2 <<
					//	 " = " << dist2 << std::endl;
					return true;
				}
			}
		}
	}
	return false;
}

/*
void
append_water_by_hbond_jump_near_atom(
	pose::Pose pose,
	ScoreFunctionOP scorefxn,
	Size seqpos,
	Size atomno,
	conformation::Residue new_rsd,
	Size new_atomno,
	Real dist_min,
	Real dist_max
)
{
	using namespace scoring::hbonds;
	Residue rsd( pose.residue( seqpos ) );
	Pose ref_pose( pose );

	//append by jump from seqpos atomno to new_rsd atom 1, maybe make random downstream atom?
	pose.append_residue_by_jump( new_rsd, seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), true );

	Size new_seqpos( pose.total_residue() );
	Size jump_number( pose.fold_tree().num_jump() );
	Jump jump( pose.jump( jump_number ) );
	//store water atom ids for clash check
	vector1< id::AtomID > clash_check_atids;
	for( Size iat = 1; iat <= new_rsd.natoms(); ++iat ){
		clash_check_atids.push_back( id::AtomID( iat, new_seqpos ) );
	}

	//which is acceptor, donor?
	Size aatm( 0 ), acc_pos( 0 ), hatm( 0 ), don_pos( 0 );
	bool wat_is_acc( false );
	//water is acceptor
	if( rsd.atom_type( atomno ).is_polar_hydrogen() &&
		new_rsd.atom_type( new_atomno ).is_acceptor() ){
		don_pos = seqpos;
		hatm = atomno;
		acc_pos = new_seqpos;
		aatm = new_atomno;
		wat_is_acc = true;
}
	//or water is donor
	else if( new_rsd.atom_type( new_atomno ).is_polar_hydrogen() &&
		rsd.atom_type( atomno ).is_acceptor() ){
		don_pos = new_seqpos;
		hatm = new_atomno;
		acc_pos = seqpos;
		aatm = atomno;
	}
	else{ utility_exit_with_message( "ERROR: res " + string_of( seqpos ) + " atom " + string_of( atomno ) +
		" res " + string_of( new_seqpos ) + " atom " + string_of( new_atomno ) + " is not HB don/acc pair!!\n" );
	}
	//need to reset actual

	//now get their base atoms to get datm and batm
	Size datm( pose.residue( don_pos ).atom_base( hatm ) );
	Size batm( pose.residue( acc_pos ).atom_base( aatm ) );
	Size b2atm( pose.residue( acc_pos ).abase2( aatm ) );
	Size dbatm( pose.residue( don_pos ).atom_base( datm ) ); //hpol base2

	//try adding vrts on the end?
	ResidueTypeSet const & rsd_set( rsd.residue_type_set() );
	conformation::ResidueOP vrt_rsd( conformation::ResidueFactory::create_residue( rsd_set.name_map( "VRT" ) ) );
	pose.append_residue_by_jump( *vrt_rsd, pose.total_residue() );
	FoldTree f_jump( pose.fold_tree() );
	//just min the new jump
	MoveMapOP mm = new MoveMap;
	mm->set_jump( jump_number, true );
	protocols::moves::MinMoverOP min_mover = new protocols::moves::MinMover( mm, scorefxn, "dfpmin", 0.01, true );

	//now go to chemical bond
	FoldTree f_rot( pose.total_residue() );

	//run through oxygen regardless
	f_rot.new_chemical_bond( seqpos, new_seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), pose.total_residue() - 2 );
	pose.fold_tree( f_rot );

	//must import from hbonds/constants.hh
	Size steps( option[ byres_data::nloop_hbscan ] ); //default 10 --> 10^5 total
	Real AHdist_min(MIN_R), AHdist_max(MAX_R); Size AHdist_steps(steps);
	Real cosBAH_min(MIN_xH), cosBAH_max(MAX_xH); Size cosBAH_steps(steps);
	Real cosAHD_min(MIN_xD), cosAHD_max(MAX_xD); Size cosAHD_steps(steps);
	Real B2BAHchi_min( 0 ), B2BAHchi_max( numeric::constants::f::pi_2 ); Size B2BAHchi_steps(steps);
	//last chi only needs [0,pi] to avoid degeneracy, right ?!?!
	Real BAHDchi_min( 0 ), BAHDchi_max( numeric::constants::f::pi ); Size BAHDchi_steps(steps);

	Pose start_pose( pose );
//	Real AHdist( 3.0 );
//	Real BAHDchi( 0 );
	for( Real AHdist = AHdist_min; AHdist <= AHdist_max;
			AHdist += (AHdist_max - AHdist_min)/static_cast<Real>(AHdist_steps-1)){
		for( Real cosBAH = cosBAH_min - 0.0001; cosBAH <= cosBAH_max;
				cosBAH += (cosBAH_max - cosBAH_min)/static_cast<Real>(cosBAH_steps-1)){
			for( Real cosAHD = cosAHD_min - .0001; cosAHD <= cosAHD_max;
					cosAHD += (cosAHD_max - cosAHD_min)/static_cast<Real>(cosAHD_steps-1)){
				for( Real B2BAHchi = B2BAHchi_min - .0001; B2BAHchi <= B2BAHchi_max;
						B2BAHchi += (B2BAHchi_max - B2BAHchi_min)/static_cast<Real>(B2BAHchi_steps-1)){
					for( Real BAHDchi = BAHDchi_min - .0001; BAHDchi <= BAHDchi_max;
							BAHDchi += (BAHDchi_max - BAHDchi_min)/static_cast<Real>(BAHDchi_steps-1)){

						//TODO: call hbonds::hbond_compute_energy( ... ) with these angles
						// and skip if is zero hb energy

						Real AHDang( numeric::constants::f::pi - std::acos( cosAHD ) );
						Real BAHang( numeric::constants::f::pi - std::acos( cosBAH ) );

						//can we replace this pose copy with just resetting the f_rot tree?
						//should be way faster, as long as not reinit'ing the angles is cool
						pose.fold_tree( f_rot );
						//pose = start_pose;

						pose.conformation().set_bond_angle( AtomID( batm, acc_pos ),
								AtomID( aatm, acc_pos  ),
								AtomID( hatm, don_pos  ),
								BAHang );

						pose.conformation().set_bond_angle( AtomID( aatm, acc_pos ),
								AtomID( hatm, don_pos  ),
								AtomID( datm, don_pos  ),
								AHDang );

						pose.conformation().set_torsion_angle( AtomID( batm, acc_pos ),
								AtomID( aatm, acc_pos  ),
								AtomID( hatm, don_pos  ),
								AtomID( datm, don_pos  ),
								BAHDchi );
						if( wat_is_acc ){
							//need to redefine hbond chi torsion if water acceptor (hbond chi undefined)
							pose.conformation().set_torsion_angle( AtomID( dbatm, don_pos ),
									AtomID( datm, don_pos  ),
									AtomID( hatm, don_pos  ),
									AtomID( aatm, acc_pos  ),
									B2BAHchi );
						}
						else{
							pose.conformation().set_torsion_angle( AtomID( b2atm, acc_pos ),
									AtomID( batm, acc_pos  ),
									AtomID( aatm, acc_pos  ),
									AtomID( hatm, don_pos  ),
									B2BAHchi );
						}

						pose.conformation().set_bond_length( AtomID( aatm, acc_pos ),
								AtomID( hatm, don_pos  ),
								AHdist );

						//do fast clash check, OH hbonds are only 0.8A!
						if( fast_clash_check( pose, clash_check_atids, 0.8 ) ) continue;

						pose.fold_tree( f_jump );
						//minimize the new jump, too slow!
						//min_mover->apply( pose );
						scorefxn->score( pose );
						//check if hbonded

						Real wat_score( pose.energies().residue_total_energies( new_seqpos ).dot( scorefxn->weights() ) );
						TR_unsat << "AHdist: " << AHdist << " ";
						TR_unsat << "\tBAHang: " << BAHang << " ";
						TR_unsat << "\tAHDang: " << AHDang << " ";
						TR_unsat << "\tB2BAHchi: " << B2BAHchi << " ";
						TR_unsat << "\tBAHDchi: " << BAHDchi << " ";
						TR_unsat << "\twat_score: " << wat_score << std::endl;

						//pose.dump_pdb( "watest." + string_of( atomno ) + "." + string_of( Size( numeric::conversions::degrees( AHDang ) ) ) + "." + string_of( Size( numeric::conversions::degrees( BAHang ) ) ) + "." + string_of( Size( numeric::conversions::degrees( B2BAHchi ) ) ) + ".pdb" );

					}
				}
			}
		}
	}

}

void
append_rsd_by_jump_near_atom(
	pose::Pose & pose,
	Size seqpos,
	Size atomno,
	conformation::Residue new_rsd,
	Size new_atomno,
	Real dist_min,
	Real dist_max
)
{
	typedef  numeric::xyzMatrix< Real > Matrix;

	Residue rsd( pose.residue( seqpos ) );
	//append by jump from seqpos atomno to new_rsd atom 1, maybe make random downstream atom?
	pose.append_residue_by_jump( new_rsd, seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), true );
	Size new_seqpos( pose.total_residue() );
	Size jump_number( pose.fold_tree().num_jump() );
	Jump jump( pose.jump( jump_number ) );

	//set jump distance as random val from dist_min to dist_max
	Real jump_dist( dist_min + numeric::random::rg().uniform() * ( dist_max - dist_min ) );
	jump.random_trans( jump_dist );
	//set jump rotation as random matrix
	jump.set_rotation( protocols::geometry::random_reorientation_matrix( 360, 360 ) );

	//and set jump in pose at our new jump
	pose.set_jump( jump_number, jump );

}

void
dock_water_to_atom(
	pose::Pose & pose,
	ScoreFunctionOP scorefxn,
	Size seqpos,
	Size atomno,
	conformation::Residue wat_rsd,
	Size new_atomno,
	Real dist_min,
	Real dist_max
)
{
	typedef  numeric::xyzMatrix< Real > Matrix;

	//attempt appending in 10 random orientations
	pose::Pose start_pose( pose );
	protocols::moves::MonteCarloOP mc_create( new protocols::moves::MonteCarlo( pose, *scorefxn, 0.8 ) );
	for( Size i = 1; i <= 10; ++i ){
		pose = start_pose;
		append_rsd_by_jump_near_atom( pose, seqpos, atomno, wat_rsd, new_atomno, dist_min, dist_max );

		Size jump_number( pose.fold_tree().num_jump() );
		//gaussian perturbations to RB dofs
		protocols::moves::MonteCarloOP mc_dock( new protocols::moves::MonteCarlo( pose, *scorefxn, 0.8 ) );
		for( Size i = 1; i <= 10; ++i ){
			Jump jump( pose.jump( jump_number ) );
			jump.gaussian_move( 1, 0.05, 90.0 );
			pose.set_jump( jump_number, jump );
			mc_dock->boltzmann( pose );
		}
		mc_dock->recover_low( pose );
		mc_create->boltzmann( pose );
	}
	mc_create->recover_low( pose );
}

void
solvate_residue_test(
	pose::Pose & pose,
	ScoreFunctionOP scorefxn,
	Size seqpos,
	Real shell_cutoff
)
{
	pose::Pose in_pose( pose );
	Real min_dist( 1.5 );
	Residue rsd( pose.residue( seqpos ) );
	ResidueTypeSet const & rsd_set( rsd.residue_type_set() );
	ResidueOP wat_rsd( ResidueFactory::create_residue( rsd_set.name_map( "TP3" ) ) );
	//turn off solvation score
	Real fa_sol_wt( scorefxn->get_weight( fa_sol ) );
	scorefxn->set_weight( fa_sol, 0.0 );

	for ( chemical::AtomIndices::const_iterator anum  = rsd.accpt_pos().begin(),
		anume = rsd.accpt_pos().end(); anum != anume; ++anum ){
		Size const iatom( *anum );
		append_rsd_by_kinematic_hbond_jump_near_atom( pose, seqpos, iatom, *wat_rsd, 2, min_dist, shell_cutoff );
		pose = in_pose;
	}
	for ( chemical::AtomIndices::const_iterator hnum  = rsd.Hpos_polar().begin(),
		hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const iatom( *hnum );
		append_rsd_by_kinematic_hbond_jump_near_atom( pose, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff );
		pose = in_pose;
	}
}

void
solvate_residue(
	pose::Pose & pose,
	ScoreFunctionOP scorefxn,
	Size seqpos,
	Real shell_cutoff
)
{
	pose::Pose in_pose( pose );
	Size const max_attempt_per_atom( 20 );
	Size const max_wat_per_atom( 5 );
	Real min_dist( 1.5 );
	Residue rsd( pose.residue( seqpos ) );
	ResidueTypeSet const & rsd_set( rsd.residue_type_set() );
	ResidueOP wat_rsd( ResidueFactory::create_residue( rsd_set.name_map( "TP3" ) ) );
	//turn off solvation score
	Real fa_sol_wt( scorefxn->get_weight( fa_sol ) );
	scorefxn->set_weight( fa_sol, 0.0 );

	//optimization stuff
	protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *scorefxn, 0.8 ) );

	//do polar hydrogens
	for ( chemical::AtomIndices::const_iterator hnum  = rsd.Hpos_polar().begin(),
		hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const iatom( *hnum );
		Size n_wat( 0 );
		//try a few times to append water molecules
		for( Size i_mc = 1; i_mc <= max_attempt_per_atom; ++i_mc ){
			if( n_wat >= max_wat_per_atom ) break;
			//append water molecule from iatom to water oxygen (atom 1)
//			append_rsd_by_hbond_jump_near_atom( pose, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff );
//			append_water_by_hbond_jump_near_atom( pose, scorefxn, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff );
//			append_rsd_by_jump_near_atom( pose, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff );
//			dock_water_to_atom( pose, scorefxn, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff );
			append_rsd_by_kinematic_hbond_jump_near_atom( pose, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff );
			pose = in_pose;
//debug
			break;
			scorefxn->score( pose );
			//just min the new jump
			MoveMapOP mm = new MoveMap;
			mm->set_jump( pose.fold_tree().num_jump(), true );
			protocols::moves::MinMoverOP min_mover = new protocols::moves::MinMover( mm, scorefxn, "dfpmin", 0.001, false );
			min_mover->apply( pose );
			if( mc->boltzmann( pose ) ) ++n_wat;
		}
	}
	//then do acceptors
	for ( chemical::AtomIndices::const_iterator anum  = rsd.accpt_pos().begin(),
		anume = rsd.accpt_pos().end(); anum != anume; ++anum ){
		Size const iatom( *anum );
		Size n_wat( 0 );
		//try a few times to append water molecules
		for( Size i_mc = 1; i_mc <= max_attempt_per_atom; ++i_mc ){
			if( n_wat >= max_wat_per_atom ) break;
			//append water molecule from iatom to water H (atom 4)
//			append_rsd_by_hbond_jump_near_atom( pose, seqpos, iatom, *wat_rsd, 4, min_dist, shell_cutoff );
//			append_water_by_hbond_jump_near_atom( pose, scorefxn, seqpos, iatom, *wat_rsd, 2, min_dist, shell_cutoff );
//			append_rsd_by_jump_near_atom( pose, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff );
//			dock_water_to_atom( pose, scorefxn, seqpos, iatom, *wat_rsd, 1, min_dist, shell_cutoff );
			append_rsd_by_kinematic_hbond_jump_near_atom( pose, seqpos, iatom, *wat_rsd, 4, min_dist, shell_cutoff );
			pose = in_pose;
//debug
			break;
			scorefxn->score( pose );
			//just min the new jump
			MoveMapOP mm = new MoveMap;
			mm->set_jump( pose.fold_tree().num_jump(), true );
			protocols::moves::MinMoverOP min_mover = new protocols::moves::MinMover( mm, scorefxn, "dfpmin", 0.001, false );
			min_mover->apply( pose );
			if( mc->boltzmann( pose ) ) ++n_wat;
		}
	}

	//reset fa_sol back to normal
	scorefxn->set_weight( fa_sol, fa_sol_wt );
}
*/


//set occupancy and bfactor data from native pose into pose
void
set_pose_occ_and_bfac(
	Pose & pose,
	Pose const native_pose,
	vector1< Size > const native_seqpos_map
)
{
	for( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ){
		if( seqpos > native_seqpos_map.size() ) continue;
		Size native_seqpos( native_seqpos_map[ seqpos ] );
//		if( pose.residue( seqpos ).name3().compare( native_pose.residue( native_seqpos_map[ seqpos ] ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
		//skip if mismatch
		if( native_seqpos == 0 ||
				pose.residue( seqpos ).name3().compare( native_pose.residue( native_seqpos ).name3() ) != 0 ||
				pose.residue( seqpos ).natoms() != native_pose.residue( native_seqpos ).natoms()
				) continue;
		Residue rsd( pose.residue( seqpos ) );
		for( Size ii = 1; ii <= rsd.natoms(); ++ii ){
			if( rsd.atom_is_hydrogen( ii ) ) continue;
			pose.pdb_info()->occupancy( seqpos, ii, native_pose.pdb_info()->occupancy( native_seqpos, ii ) );
			pose.pdb_info()->temperature( seqpos, ii, native_pose.pdb_info()->temperature( native_seqpos, ii ) );
		}
	}
}

//remove all water molecules
Pose
remove_pose_water(
	pose::Pose pose
)
{
	Size nres( pose.total_residue() );
	//delete from end to being so dont mess up seqpos numbering
	for( Size seqpos = nres; seqpos >= 1; --seqpos ){
		if( pose.residue( seqpos ).name1() == 'w' ){
			pose.delete_polymer_residue( seqpos );
		}
	}
	return pose;
}

//normalize residue sasa by exposed value
Real
normalize_residue_sasa(
	pose::Pose pose,
	Size const seqpos,
	Real const res_sasa
)
{
	Size aaidx( pose.residue( seqpos ).aa() );
	switch( aaidx ){
    case  1:
			return res_sasa / 170; // 1 A
    	break;
		case  2:
			return res_sasa / 170; // 2 C
    	break;
		case  3:
			return res_sasa / 210; // 3 D
    	break;
		case  4:
			return res_sasa / 250; // 4 E
    	break;
		case  5:
			return res_sasa / 290; // 5 F
    	break;
		case  6:
			return res_sasa / 170; // 6 G
    	break;
		case  7:
			return res_sasa / 220; // 7 H
    	break;
		case  8:
			return res_sasa / 230; // 8 I
    	break;
		case  9:
			return res_sasa / 260; // 9 K
    	break;
		case 10:
			return res_sasa / 230; // 10 L
    	break;
		case 11:
			return res_sasa / 240; // 11 M
    	break;
		case 12:
			return res_sasa / 190; // 12 N
    	break;
		case 13:
			return res_sasa / 220; // 13 P
    	break;
		case 14:
			return res_sasa / 220; // 14 Q
    	break;
		case 15:
			return res_sasa / 260; // 15 R
    	break;
		case 16:
			return res_sasa / 180; // 16 S
    	break;
		case 17:
			return res_sasa / 200; // 17 T
    	break;
		case 18:
			return res_sasa / 200; // 18 V
    	break;
		case 19:
			return res_sasa / 300; // 19 W
    	break;
		case 20:
			return res_sasa / 290; // 20 Y
			break;
		default:
			return res_sasa;
			break;
	}
}

//get lk energy of one atom by a given residue
//graciously stolen from pbradley
Real
get_atom_lk_energy_by_residue_no_count_pair(
		scoring::ScoreFunctionOP scorefxn,
		conformation::Residue const & rsd1,
		Size const atom1,
		conformation::Residue const & rsd2
)
{

	ObjexxFCL::FArray3D< Real > const & solv1( ScoringManager::get_instance()->etable(
		scorefxn->energy_method_options().etable_type() )->solv1() );
	Real const safe_max_dis2( ScoringManager::get_instance()->etable(
		scorefxn->energy_method_options().etable_type() )->get_safe_max_dis2() );
	Real const etable_bins_per_A2( ScoringManager::get_instance()->etable(
		scorefxn->energy_method_options().etable_type() )->get_bins_per_A2() );

	// setup residue information
	Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
	Size const atom1_type_index( rsd1.atom( atom1 ).type() );

	Real total_lk_energy( 0.0 );
	for ( Size atom2=1; atom2<= rsd2.nheavyatoms(); ++atom2 ) {
		Vector const & atom2_xyz( rsd2.xyz( atom2 ) );

		Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

		if ( ( d2 >= safe_max_dis2) || ( d2 < 1e-3 ) ) continue; // exclude self...

		// setup for solvation Etable lookups
		Size const atom2_type_index( rsd2.atom( atom2 ).type() );
		Real const d2_bin = d2 * etable_bins_per_A2;
		int disbin = static_cast< int >( d2_bin ) + 1;
		Real  frac = d2_bin - ( disbin - 1 );
		int const l1 = solv1.index( disbin, atom2_type_index, atom1_type_index );

		Real const lk_energy_of_atom1_by_atom2
			( ( ( 1. - frac ) * solv1[ l1 ] + frac * solv1[ l1+1 ] ) );

		total_lk_energy += lk_energy_of_atom1_by_atom2;
//		std::cout << rsd1.name3() + " " + rsd1.atom_name( atom1 ) + " <- " + rsd2.name3() + " " + rsd2.atom_name( atom2 ) + "\t" + string_of( lk_energy_of_atom1_by_atom2 ) + "\n";

	}
	return total_lk_energy;
}

/*
//total LK burial for single atom, sums over self residue atoms and energy graph nbr atoms
Real
calc_lk_burial_for_single_atom(
		Size const atom1,
		conformation::Residue const & rsd1,
		pose::Pose const & pose
		)
{


	PROF_START( util::CALC_LK_BURIAL_FOR_SINGLE_ATOM );

	/// this could be bad if atom type sets don't match up
	static chemical::ResidueTypeSet const * fa_standard_rsd_set
		( & ( *chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) ) );

	static methods::EnergyMethodOptions energy_method_options;

	static methods::LK_BallEnergy lk_ball_energy( energy_method_options );


	if ( & ( rsd1.residue_type_set() ) != fa_standard_rsd_set ) {
		TR.Trace << "skipping lk desolvation calculation for non-fa-standard residue" << std::endl;
		return 0.0;
	}

	/// what is this atom's desolvation parameter?

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	Size const pos1( rsd1.seqpos() );

	Real weighted_desolvation_no_count_pair( 0.0 );
	Real const lk_dgfree( -1.0 * rsd1.atom_type( atom1 ).lk_dgfree() );
	for ( graph::Graph::EdgeListConstIter
			ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			ire = energy_graph.get_node( pos1 )->const_edge_list_end();
			ir != ire; ++ir ) {
		EnergyEdge const * edge( static_cast< EnergyEdge const *> (*ir) );
		Size const pos2( edge->get_other_ind( pos1 ) );
		conformation::Residue const & rsd2( pose.residue( pos2 ) );
		if ( rsd2.aa() == chemical::aa_h2o || rsd2.aa() == chemical::aa_vrt ) continue;
		assert( pos2 != pos1 );
		weighted_desolvation_no_count_pair +=
			lk_ball_energy.calculate_lk_desolvation_of_single_atom_by_residue_no_count_pair( atom1, rsd1, rsd2 );
	}

	/// add something here to check for water rotamers built where there was a virtual residue -- will have no neighbors!
	if ( energy_graph.get_node( pos1 )->const_edge_list_begin() ==
			energy_graph.get_node( pos1 )->const_edge_list_end() ) {
		TR.Trace << "calc_lk_desolvation_for_single_atom: no nbrs!" << std::endl;
	}

	weighted_desolvation_no_count_pair +=
		lk_ball_energy.calculate_lk_desolvation_of_single_atom_by_residue_no_count_pair( atom1, rsd1, rsd1 );

	PROF_STOP( util::CALC_LK_BURIAL_FOR_SINGLE_ATOM );

	return weighted_desolvation_no_count_pair / lk_dgfree;

}
*/

Real
get_atom_lk_energy(
	pose::Pose const pose,
	scoring::ScoreFunctionOP scorefxn,
	Size const seqpos,
	Size const iatom
)
{
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	Real fa_sol_atom = 0.;
	core::conformation::Residue const & resl( pose.residue( seqpos ) );
//	core::conformation::Atom const & atoml( resl.atom(iatom) );
	scoring::etable::TableLookupEtableEnergy const etable_energy(
		*ScoringManager::get_instance()->etable( scorefxn->energy_method_options().etable_type() ),
		scorefxn->energy_method_options()
	);
	//over all other residue nbrs in energy graph
	for ( graph::Graph::EdgeListConstIter
			ir  = energy_graph.get_node( seqpos )->const_edge_list_begin(),
			ire = energy_graph.get_node( seqpos )->const_edge_list_end();
			ir != ire; ++ir ) {
		EnergyEdge const * edge( static_cast< EnergyEdge const *> (*ir) );
		Size const nbr_seqpos( edge->get_other_ind( seqpos ) );
		core::conformation::Residue const & resu( pose.residue( nbr_seqpos ) );
		//skip solvent atoms!
		if( resu.name1() == 'w' ) continue;

		/*
		//this just gets symmetrized energy!!
		//over all atoms in nbr residue
		for ( Size jatom=1, jatom_end = resu.nheavyatoms(); jatom <= jatom_end; ++jatom ) {
			core::conformation::Atom const & atomu( resu.atom(jatom) );
			core::DistanceSquared dsq;
			core::scoring::Weight weight = 1.0;
			core::Size path_dist(0);
			core::Energy atr, rep, solv, bb;
			dsq = atoml.xyz().distance_squared( atomu.xyz() );

			core::scoring::etable::count_pair::CountPairFunctionOP cpfxn = core::scoring::etable::count_pair::CountPairFactory::create_count_pair_function( resl, resu, scoring::etable::count_pair::CP_CROSSOVER_4 );
			if ( cpfxn->count( iatom, jatom, weight, path_dist ) ) {
				etable_energy.atom_pair_energy(atoml,atomu,weight,atr,rep,solv,bb,dsq);
				fa_sol_atom += solv;
			}
		}
		*/
		fa_sol_atom += get_atom_lk_energy_by_residue_no_count_pair( scorefxn, resl, iatom, resu );
	}

	//and get burial from other atoms in self residue
	fa_sol_atom += get_atom_lk_energy_by_residue_no_count_pair( scorefxn, resl, iatom, resl );

	return fa_sol_atom;
}

Real
get_atom_lk_burial(
	pose::Pose const pose,
	scoring::ScoreFunctionOP scorefxn,
	Size const seqpos,
	Size const iatom
)
{
	Real const atom_lk_energy( get_atom_lk_energy( pose, scorefxn, seqpos, iatom ) );
	Real const lk_dgfree( pose.residue( seqpos ).atom_type( iatom ).lk_dgfree() );
	//HACKATTACK! carbonyl carbons have 0 dgfree! returns bogus -1 if dgfree == 0
	if( lk_dgfree == 0 ) return -1;
	Real atom_lk_burial( atom_lk_energy / ( -1 * lk_dgfree ) );
//	std::cout << pose.residue( seqpos ).atom_name( iatom ) + "\t" << atom_lk_burial << "\t" + string_of( lk_dgfree ) + "\n";
	return atom_lk_burial;
}

//add up atomic lk burials, normalize by nheavyatoms
Real
get_res_avg_lk_burial(
	pose::Pose const pose,
	scoring::ScoreFunctionOP scorefxn,
	Size const seqpos,
	bool incl_bb,
	bool incl_sc
)
{
	Real res_lk_burial( 0. );
	Size n_atoms( 0 );
	//HACKATTACK! carbonyl carbons have 0 dgfree! skip val if is NAN
	for( Size iatom = 1; iatom <= pose.residue( seqpos ).nheavyatoms(); ++iatom ){
		//exclude bb or sc atoms?
		if( !incl_bb && iatom < pose.residue( seqpos ).first_sidechain_atom() ) continue;
		if( !incl_sc && iatom >= pose.residue( seqpos ).first_sidechain_atom() ) continue;
		Real atom_lk_burial( get_atom_lk_burial( pose, scorefxn, seqpos, iatom ) );
		//check for bogus neg value
		if( atom_lk_burial < 0 ) continue;
		res_lk_burial += get_atom_lk_burial( pose, scorefxn, seqpos, iatom );
		++n_atoms;
	}
	//return 0 if skipped all atoms
	if( n_atoms == 0 ) return 0;
	return res_lk_burial / n_atoms;
}

/*
//get sum of sr and lr hbond_bb by res
Real
get_res_hbond_bb_score_raw(
	pose::Pose const pose,
	Size const seqpos
)
{
	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	HBondSet hbond_bb_set;
	//exclude bbsc, scbb, and scsc
	fill_hbond_set( pose, false, hbond_bb_set, false, true, true, true );
	Real rsd_hbond_bb_energy( 0.0 );
	for( Size i_hb = 1; i_hb <= Size ( hbond_bb_set.nhbonds() ); ++i_hb ){
			HBond hb( hbond_bb_set.hbond( i_hb ) );
			if ( hb.don_res() == seqpos || hb.acc_res() == seqpos ) {
					//divide by 2 for half of interacting pair
					rsd_hbond_bb_energy += ( hb.energy() / 2.0 );
			}
	}
	return rsd_hbond_bb_energy;
}
*/
/*
//get # sidechain hbonds
Size
get_n_bb_hbonds(
		Pose const pose,
		Size const seqpos
		)
{
	using namespace core::scoring::hbonds;
	Size n_hbonds( 0 );
	HBondSet hbset;
	hbset.setup_for_residue_pair_energies( pose, false, false );
	for( Size i = 1; i <= hbset.nhbonds(); ++i ){
		HBond hb( hbset.hbond( i ) );
		if( hb.don_res() == seqpos && hb.don_hatm_is_protein_backbone() ) ++n_hbonds;
		else if( hb.acc_res() == seqpos && hb.acc_atm_is_protein_backbone() ) ++n_hbonds;
	}
	return n_hbonds;
}
*/
//get # hbonds by type
void
get_n_hbonds(
		Pose const pose,
		Size const seqpos,
		Size & n_sc_hbonds,
		Size & n_bb_hbonds,
		Size & n_sc_wat_hbonds,
		Size & n_bb_wat_hbonds,
		Size & n_hbonds,
		Size & n_unsat
		)
{
	using namespace core::scoring::hbonds;
	HBondSet hbset;
	hbset.setup_for_residue_pair_energies( pose, false, false );
	conformation::Residue rsd( pose.residue( seqpos ) );
	//count n hbonds for donors and acceptors
	vector1< Size > n_atom_hbonds( rsd.natoms(), 0 );
	for( Size i = 1; i <= hbset.nhbonds(); ++i ){
		HBond hb( hbset.hbond( i ) );
		if( hb.don_res() == seqpos ){
			Size don_atm( rsd.atom_base( hb.don_hatm() ) );
			n_atom_hbonds[ don_atm ] += 1;
			if( hb.acc_res_is_protein() || hb.acc_res_is_dna() ){
				if( hb.don_hatm_is_protein_backbone() ) ++n_bb_hbonds;
				else ++n_sc_hbonds;
			}
			else if( pose.residue( hb.acc_res() ).name1() == 'w' ){
				if( hb.don_hatm_is_protein_backbone() ) ++n_bb_wat_hbonds;
				else ++n_sc_wat_hbonds;
			}
		}
		else if( hb.acc_res() == seqpos ){
			n_atom_hbonds[ hb.acc_atm() ] += 1;
			if( hb.don_res_is_protein() || hb.don_res_is_dna() ){
				if( hb.acc_atm_is_protein_backbone() ) ++n_bb_hbonds;
				else ++n_sc_hbonds;
			}
			else if( pose.residue( hb.don_res() ).name1() == 'w' ){
					if( hb.acc_atm_is_protein_backbone() ) ++n_bb_wat_hbonds;
					else ++n_sc_wat_hbonds;
			}
		}
	}
	n_hbonds = n_sc_hbonds + n_bb_hbonds + n_sc_wat_hbonds + n_bb_wat_hbonds;
	//calc unsat hbonds a la BuriedUnsatisfied calculator
	for( Size atm = 1; atm <= rsd.nheavyatoms(); ++atm ){
		if( !( rsd.atom_type( atm ).is_acceptor() || rsd.atom_type( atm ).is_donor() ) ) continue;
		//we need at least one hbond
		if( n_atom_hbonds[ atm ] == 0 ){
			n_unsat += 1; continue;
		}
		//may need > 1 hbond, see BuriedUnsatCalc
		//behavior copied from Florian!
		std::string atom_type( rsd.type().atom_type( atm ).name() );
		Size satisfac_cut = 3;
		if( atom_type == "OH" || atom_type == "OCbb" || atom_type == "S" ){
			satisfac_cut = 2;
		}
		Size bonded_heavyatoms = rsd.n_bonded_neighbor_all_res( atm )
			- rsd.type().number_bonded_hydrogens( atm );
		if( bonded_heavyatoms + n_atom_hbonds[ atm ] < satisfac_cut ){
			n_unsat += 1;
		}
	}
}

//get total res-water energy, ignores solvation!
Real
get_res_water_energy(
	pose::Pose pose,
	ScoreFunction scorefxn,
	Size seqpos
)
{
	//turn off solvation!
	scorefxn.set_weight( fa_sol, 0.0 );
	Real total_water_energy( 0.0 );
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( graph::Graph::EdgeListConstIter
			ir  = energy_graph.get_node( seqpos )->const_edge_list_begin(),
			ire = energy_graph.get_node( seqpos )->const_edge_list_end();
			ir != ire; ++ir ){
		Size nbr( (*ir)->get_other_ind( seqpos ) );
		if( pose.residue( nbr ).name1() != 'w' ) continue;
		EnergyEdge const * edge( static_cast< const EnergyEdge *> (*ir) );
		EnergyMap const & emap( edge->fill_energy_map());
		total_water_energy += emap.dot( scorefxn.weights() );
	}
	return total_water_energy;
}


//get sidechain avg bfactor of seqpos, ignore hydrogens
Real
get_sc_bfactor(
		Pose const native_pose,
		Size const seqpos
		)
{
	if( seqpos == 0 ) return 0.0;
	Residue native_rsd( native_pose.residue( seqpos ) );
	Real sc_bfactor( 0.0 );
	Size count( 0 );
	//skip if no sidechain
	if( native_rsd.first_sidechain_atom() > native_rsd.natoms() ) return 0;
	//get bfacs
	for( Size ii = native_rsd.first_sidechain_atom(); ii <= native_rsd.natoms(); ++ii ){
		if( native_rsd.atom_is_hydrogen( ii ) ) continue;
		if( native_pose.pdb_info()->occupancy( seqpos, ii ) >= 0.5 ){
			sc_bfactor += native_pose.pdb_info()->temperature( seqpos, ii );
			++count;
		}
	}
	sc_bfactor /= count;
	return sc_bfactor;
}

/*
//get sidechain rmsd of seqpos, ignore hydrogens and low occ
Real
get_sc_rmsd(
		Pose const pose,
		Pose const native_pose,
		Size const seqpos
		)
{
	if( pose.residue( seqpos ).name3().compare( native_pose.residue( seqpos ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
	Residue rsd( pose.residue( seqpos ) );
	Residue native_rsd( native_pose.residue( seqpos ) );
	Real sc_rmsd( 0.0 );
	Size count( 0 );
	for( Size ii = rsd.first_sidechain_atom(); ii <= rsd.natoms(); ++ii ){
		if( rsd.atom_is_hydrogen( ii ) ) continue;
		if( native_pose.pdb_info()->occupancy( seqpos, ii ) >= 0.5 ){
			sc_rmsd += rsd.xyz( ii ).distance_squared( native_rsd.xyz( ii ) );
			++count;
		}
	}
	sc_rmsd = std::sqrt( sc_rmsd / count );
	return sc_rmsd;
}
*/

//get automorphic rmsd
Real
get_ca_distance(
		Pose pose,
		Pose ref_pose,
		Size const seqpos,
		Size const ref_seqpos
		)
{
	return pose.residue( seqpos ).atom( "CA" ).xyz().distance( ref_pose.residue( ref_seqpos ).atom( "CA" ).xyz() );
}

//get automorphic rmsd
Real
get_sc_automorphic_rmsd(
		Pose pose,
		Pose ref_pose,
		Size const seqpos,
		Size const ref_seqpos,
		bool const super
		)
{
	//return bogus value if mismatch
	if( ref_seqpos == 0 ||
			pose.residue( seqpos ).name3().compare( ref_pose.residue( ref_seqpos ).name3() ) != 0 ||
			pose.residue( seqpos ).natoms() != ref_pose.residue( ref_seqpos ).natoms()
			) return 9999.;
//	if( pose.residue( seqpos ).name3().compare( ref_pose.residue( ref_seqpos ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
	//add variant type to virt bb atoms
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_BB", seqpos );
	pose::add_variant_type_to_pose_residue( ref_pose, "VIRTUAL_BB", ref_seqpos );

	Residue rsd( pose.residue( seqpos ) );
	Residue ref_rsd( ref_pose.residue( ref_seqpos ) );

	Real sc_rmsd( scoring::automorphic_rmsd( pose.residue( seqpos ), ref_pose.residue( ref_seqpos ), super /*superpose*/ ) );
	return sc_rmsd;
}

/*
//get sidechain rmsd of seqpos, ignore hydrogens and low occ
Real
get_three_atom_sc_rmsd(
		Pose const pose,
		Pose const native_pose,
		Size const seqpos
		)
{
	if( pose.residue( seqpos ).name3().compare( native_pose.residue( seqpos ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
	Residue rsd( pose.residue( seqpos ) );
	Residue native_rsd( native_pose.residue( seqpos ) );
	vector1< std::string > atom_ids( sc_rmsd_AtomIDs[ rsd.name3() ] );
	Real sc_rmsd( 0.0 );
	Size count( 0 );
	for( Size ii = 1; ii <= atom_ids.size(); ++ii ){
		sc_rmsd += rsd.xyz( rsd.atom_index( atom_ids[ ii ] ) ).distance_squared( native_rsd.xyz( native_rsd.atom_index( atom_ids[ ii ] ) ) );
		++count;
	}
	sc_rmsd = std::sqrt( sc_rmsd / count );
	return sc_rmsd;
}
*/

void
split_fa_dun(
		Pose const pose,
		Size const seqpos,
		Real & fa_dun_rot,
		Real & fa_dun_dev,
		pack::dunbrack::RotVector & rotvec
		)
{

	using namespace scoring;
	using namespace pack::dunbrack;

	conformation::Residue const & rsd( pose.residue( seqpos ) );
	bool is_dun02( true );
	if( option[ corrections::score::dun10 ].user() ) is_dun02 = false;
	RotamerLibrary const & rlib( * RotamerLibrary::get_instance() );
	RotamerLibraryScratchSpace scratch;

	Real fa_dun_tot = rlib.rotamer_energy( rsd, scratch );
	fa_dun_rot = scratch.fa_dun_rot();
	fa_dun_dev = scratch.fa_dun_dev();
	//semirotameric returns 0 for fa_dun_rot!
	if( abs( fa_dun_tot - fa_dun_rot - fa_dun_dev ) > 0.01 ) fa_dun_rot = ( fa_dun_tot - fa_dun_dev );

	rotamer_from_chi( rsd, rotvec );
}

std::string
get_aa_torsion_string(
		Pose const pose,
		Size const seqpos
		)
{
	std::string tor_anl;
	tor_anl += ( "phi: " + string_of( pose.phi( seqpos ) ) + " " );
	tor_anl += ( "psi: " + string_of( pose.psi( seqpos ) ) + " " );
	if( pose.residue( seqpos ).nchi() < 1 ){
		tor_anl += ( "chi1: 0 chi2: 0 chi3: 0 chi4: 0 " );
		return tor_anl;
	}
	tor_anl += ( "chi1: " + string_of( pose.chi( 1, seqpos ) ) + " " );
	if( pose.residue( seqpos ).nchi() < 2 ){
		tor_anl += ( "chi2: 0 chi3: 0 chi4: 0 " );
		return tor_anl;
	}
	tor_anl += ( "chi2: " + string_of( pose.chi( 2, seqpos ) ) + " " );
	if( pose.residue( seqpos ).nchi() < 3 ){
		tor_anl += ( "chi3: 0 chi4: 0 " );
		return tor_anl;
	}
	tor_anl += ( "chi3: " + string_of( pose.chi( 3, seqpos ) ) + " " );
	if( pose.residue( seqpos ).nchi() < 4 ){
		tor_anl += ( "chi4: 0 " );
		return tor_anl;
	}
	tor_anl += ( "chi4: " + string_of( pose.chi( 4, seqpos ) ) + " " );
	return tor_anl;
}

//////OUTPUT///////
void
get_res_data_ss(
		io::silent::SilentStructOP & ss,
		Pose pose,
		Pose const native_pose,
		Size const native_seqpos,
		Size const seqpos,
		ScoreFunctionOP scorefxn,
		ScoreFunctionOP scorefxn_edens,
		bool const do_sc_rmsd
//		utility::vector1< core::Real > residue_sasa
		)
{
	using namespace core;
	using namespace core::scoring;


	std::map< std::string, Real > res_data_map;

	scorefxn->score( pose );
	EnergyMap weights( pose.energies().weights() );
	EnergyMap rsd_energies( pose.energies().residue_total_energies( seqpos )  );
	EnergyMap rsd_energies_nat( native_pose.energies().residue_total_energies( native_seqpos )  );

	Real total_score( pose.energies().residue_total_energies( seqpos ).dot( scorefxn->weights() ) );
	Real total_score_nat( native_pose.energies().residue_total_energies( native_seqpos ).dot( scorefxn->weights() ) );
	ss->add_energy( "score", total_score );
	ss->add_energy( "score_d", total_score - total_score_nat );

	for ( int ii = 1; ii <= scoring::n_score_types; ++ii ) {
		if ( weights[ ScoreType(ii) ] != 0.0 ) {
			Real const value( rsd_energies[ ScoreType(ii) ] );
			Real const value_nat( rsd_energies_nat[ ScoreType(ii) ] );
			std::string const scorename( name_from_score_type( ScoreType(ii) ) );
			ss->add_energy( scorename + "_d", value - value_nat );
		}
	}

	//dunbrack split
	Real fa_dun_rot, fa_dun_dev;
	pack::dunbrack::RotVector rotvec;
	Real fa_dun_rot_nat, fa_dun_dev_nat;
	pack::dunbrack::RotVector rotvec_nat;
	split_fa_dun( pose, seqpos, fa_dun_rot, fa_dun_dev, rotvec );
	split_fa_dun( native_pose, native_seqpos, fa_dun_rot_nat, fa_dun_dev_nat, rotvec_nat );
	Real fa_dun_wt( scorefxn->get_weight( fa_dun ) );
	fa_dun_rot *= fa_dun_wt;
	fa_dun_dev *= fa_dun_wt;
	fa_dun_rot_nat *= fa_dun_wt;
	fa_dun_dev_nat *= fa_dun_wt;
	ss->add_energy( "fa_dun_rot", fa_dun_rot );
	ss->add_energy( "fa_dun_dev", fa_dun_dev );
	ss->add_energy( "fa_dun_rot_d", fa_dun_rot - fa_dun_rot_nat );
	ss->add_energy( "fa_dun_dev_d", fa_dun_dev - fa_dun_dev_nat );

	//res hbond_bb energy
//	Real hbond_bb_score_raw( get_res_hbond_bb_score_raw( pose, seqpos ) );
//	ss->add_energy( "hbond_bb_raw", hbond_bb_score_raw );

	//avg atomic lk burial in residue
	Real res_lk_burial( get_res_avg_lk_burial( pose, scorefxn, seqpos, true, true ) );
	Real bb_lk_burial( get_res_avg_lk_burial( pose, scorefxn, seqpos, true, false ) );
	Real sc_lk_burial( get_res_avg_lk_burial( pose, scorefxn, seqpos, false, true ) );
	ss->add_energy( "res_lk_burial", res_lk_burial );
	ss->add_energy( "bb_lk_burial", bb_lk_burial );
	ss->add_energy( "sc_lk_burial", sc_lk_burial );

/*
	//residue sasa
	if( option[ byres_data::calc_sasa ] ){
			Real residue_sasa_norm( normalize_residue_sasa( pose, seqpos, residue_sasa[ seqpos ] ) );
			ss->add_energy( "res_sasa", residue_sasa_norm );
	}
*/

	//electron density res score
	Real edens_score( 0.0 );
	if( option[ edensity::mapfile ].user() ){
		scorefxn_edens->score( pose );
		edens_score = pose.energies().residue_total_energies( seqpos )[ elec_dens_window ];
	}
	ss->add_energy( "elec_dens_window", edens_score );

	//n hbonds
/*
	Size n_sc_hbonds( 0 ), n_bb_hbonds( 0 ), n_sc_wat_hbonds( 0 ), n_bb_wat_hbonds( 0 ), n_hbonds( 0 ), n_unsat_hbonds( 0 );
	get_n_hbonds( pose, seqpos, n_sc_hbonds, n_bb_hbonds, n_sc_wat_hbonds, n_bb_wat_hbonds, n_hbonds, n_unsat_hbonds );
	ss->add_energy( "n_sc_hbonds", static_cast< Real >( n_sc_hbonds ) );
	ss->add_energy( "n_bb_hbonds", static_cast< Real >( n_bb_hbonds ) );
	ss->add_energy( "n_sc_wat_hbonds", static_cast< Real >( n_sc_wat_hbonds ) );
	ss->add_energy( "n_bb_wat_hbonds", static_cast< Real >( n_bb_wat_hbonds ) );
	ss->add_energy( "n_hbonds", static_cast< Real >( n_hbonds ) );
	ss->add_energy( "n_unsat_hbonds", static_cast< Real >( n_unsat_hbonds ) );
*/

	//water energy
	Real res_water_energy( get_res_water_energy( pose, *scorefxn, seqpos ) );
	ss->add_energy( "res_water_energy", res_water_energy );

	//things after this point should not be delta'ed

	//this will xform 4 rotbin indices into a 4 digit number
	Real rotbin_val( 0 );
	Size n_rotbins( 4 );
	for( Size i_rotvec = 1; i_rotvec <= n_rotbins; ++i_rotvec ){
		if( i_rotvec <= rotvec.size() ) rotbin_val += ( ( rotvec[ i_rotvec ] + 1 ) * std::pow( 10.0, static_cast< Real >( 4 - i_rotvec ) ) );
		else rotbin_val += std::pow( 10.0, static_cast< Real >( 4 - i_rotvec ) );
	}
	ss->add_energy( "rotbin", rotbin_val );

	//torsion angles
	ss->add_energy( "phi", pose.phi( seqpos ) );
	ss->add_energy( "psi", pose.psi( seqpos ) );
  vector1< Real > chi_data( n_rotbins, 0. );
  vector1< Real > chis( pose.residue( seqpos ).chi() );
  for( Size i = 1; i <= chis.size(); ++i ) chi_data[ i ] = chis[ i ];
  ss->add_energy( "chi1", chi_data[ 1 ] );
  ss->add_energy( "chi2", chi_data[ 2 ] );
  ss->add_energy( "chi3", chi_data[ 3 ] );
  ss->add_energy( "chi4", chi_data[ 4 ] );

  vector1< Real > nat_chi_data( n_rotbins, 0. );
  vector1< Real > nat_chis( native_pose.residue( native_seqpos ).chi() );
  for( Size i = 1; i <= nat_chis.size(); ++i ) nat_chi_data[ i ] = nat_chis[ i ];
  vector1< Real > d_chi_data( n_rotbins, 0. );
  for( Size i = 1; i <= chi_data.size(); ++i ){
		d_chi_data[ i ] = std::abs( chi_data[ i ] - nat_chi_data[ i ] );
		if ( d_chi_data[ i ] > 180. ) d_chi_data[ i ] -= 180.;
	}
  ss->add_energy( "d_chi1", d_chi_data[ 1 ] );
  ss->add_energy( "d_chi2", d_chi_data[ 2 ] );
  ss->add_energy( "d_chi3", d_chi_data[ 3 ] );
  ss->add_energy( "d_chi4", d_chi_data[ 4 ] );

	//sidechain bfactor from native
	Real sc_bfactor( get_sc_bfactor( native_pose, native_seqpos ) );
	ss->add_energy( "sc_bfactor", sc_bfactor );

	//sidechain rmsd
	Real sc_auto_rmsd_nat( 0 );
	if( do_sc_rmsd ){
		sc_auto_rmsd_nat = get_sc_automorphic_rmsd( pose, native_pose, seqpos, native_seqpos, false );
		ss->add_energy( "sc_rmsd", sc_auto_rmsd_nat );
		sc_auto_rmsd_nat = get_sc_automorphic_rmsd( pose, native_pose, seqpos, native_seqpos, true );
		ss->add_energy( "sc_rmsd_super", sc_auto_rmsd_nat );
		Real bb_ca_dist( get_ca_distance( pose, native_pose, seqpos, native_seqpos ) );
		ss->add_energy( "bb_ca_dist", bb_ca_dist );
	}

	Size aaidx( pose.residue( seqpos ).aa() );
	ss->add_energy( "aaidx", static_cast< Real >( aaidx ) );

	ss->add_energy( "seqpos", static_cast< Real >( seqpos ) );
	ss->add_string_value( "aa", pose.residue( seqpos ).name3() );
	ss->add_string_value( "pdb_chain", utility::string_split( pose.pdb_info()->pose2pdb( seqpos ) )[ 1 ] );
	ss->add_string_value( "pdb_resnum", utility::string_split( pose.pdb_info()->pose2pdb( seqpos ) )[ 2 ] );

}

///////////////////////////////////////////////////////////////////////////////
	void
byres_analysis(
	pose::Pose pose
)
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//tag is "." or "tag."
	std::string tag( option[ byres_data::tag ] );
	tag = tag + ".";
	std::string const scname_byres( tag + "byres.sc" );
	std::string const scname_byatom( tag + "byatom.sc" );
	std::string const scname( tag + "sc" );

	//do rmsd calc?
	bool const calc_rmsd( option[ in::file::native ].user() );

	std::string pdbname( core::pose::tag_from_pose( pose ) );
	/*
	std::string pdbname( start_file() );
	Size pdbnamestart( 0 );
	if( pdbname.find_last_of( "/" ) < ( pdbname.size() - 1 ) ) pdbnamestart = pdbname.find_last_of( "/" ) + 1;
	std::string pdbname( pdbname, pdbnamestart, pdbname.size() - pdbnamestart - 4 );
	Pose pose;
	pose_from_pdb( pose, pdbname );
	*/

	//create a ScoreFunction from commandline options
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	//and add decompose bb_hbond energies option
	core::scoring::methods::EnergyMethodOptionsOP emopts(
		new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() )
	);
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
//	emopts->hbond_options().ignore_water_hb_env_dep( true );
	scorefxn->set_energy_method_options( *emopts );
	//clone into edens scorefxn
	core::scoring::ScoreFunctionOP scorefxn_edens( ( *scorefxn ).clone() );

	//add edens score?
	if( option[ edensity::mapfile ].user() ){
		protocols::electron_density::SetupForDensityScoringMoverOP edens_mover(
			new protocols::electron_density::SetupForDensityScoringMover );
		edens_mover->apply( pose );
		std::string const pdbname_out( pdbname + ".edock." + "pdb" );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn_edens );
		pose.dump_scored_pdb( pdbname_out, *scorefxn_edens );
	}

	//set non-edens scorefxn edens wt to 0.0 so can rtmin and stuff
	scorefxn->set_weight( elec_dens_window, 0.0 );

	std::string native_pdbname( pdbname );
	Pose native_pose( pose );
	//init native_seqpos_map to point to self
	vector1< Size > native_seqpos_map;
	for( Size iseq = 1; iseq <= pose.total_residue(); ++iseq ){ native_seqpos_map.push_back( iseq ); }
	//actually have a diff native pose?
	if( option[ in::file::native ].user() ){
		native_pdbname = option[ in::file::native ]();
		pose_from_pdb( native_pose, native_pdbname );
		//align structures for edensity cals?
		if( option[ edensity::mapfile ].user() ) core::scoring::calpha_superimpose_pose( native_pose, pose );
		//get sequence mapping?
		if( option[ byres_data::align_native_seq ] ){
			using namespace sequence;
			utility::file::FileName blosum62( "/work/chrisk/rosetta/rosetta_database/sequence/substitution_matrix/IDENTITY" ); //TODO: can we get the rosetta db env variable?
			ScoringSchemeOP ss( new MatrixScoringScheme( -30, -5, blosum62 ) );
			NWAligner nw_aligner;
			SequenceOP seq1( new Sequence( pose.sequence(), "target", 1 ) );
			SequenceOP seq2( new Sequence( native_pose.sequence(), "ref", 1 ) );
			SequenceAlignment global_align = nw_aligner.align( seq1, seq2, ss );
			SequenceMapping seq_map( global_align.sequence_mapping( 1, 2 ) );
			seq_map.show( TR );
			native_seqpos_map = seq_map.mapping();
		}
		set_pose_occ_and_bfac( pose, native_pose, native_seqpos_map );
	}

	scorefxn->score( pose );
	scorefxn->score( native_pose );

/*
	utility::vector1< core::Real > residue_sasa;
	//calc sasa, must remove water molecules first!
	core::id::AtomID_Map< core::Real > atom_sasa;
	if( option[ byres_data::calc_sasa ] ){
		pose::Pose dry_pose( remove_pose_water( pose ) );
		core::Real const probe_radius(1.4);
		core::scoring::calc_per_atom_sasa( dry_pose, atom_sasa, residue_sasa, probe_radius);
	}
*/

	//silent-type output
	io::silent::SilentFileData pose_silent_data;
	io::silent::SilentFileData res_silent_data;
	io::silent::SilentFileData atom_silent_data;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//create a task factory to get resfile data
	core::pack::task::TaskFactoryOP task_factory = new core::pack::task::TaskFactory;
	//use resfile to note residues to analyze...
	Pose start_pose( pose );
	if( option[ packing::resfile ].user() ) task_factory->push_back( new core::pack::task::operation::ReadResfile );
	else task_factory->push_back( new core::pack::task::operation::RestrictToRepacking() );
	core::pack::task::PackerTaskOP packer_task( task_factory->create_task_and_apply_taskoperations( pose ) );

	Size nres( pose.total_residue() );

	//solvate all polar groups?
/*
	if( option[ byres_data::solvate ] ){
//debug
		for( Size seqpos = 1; seqpos <= nres; ++seqpos ){
			solvate_residue_test( pose, scorefxn, seqpos, 4.0 );
		}
		pose.dump_scored_pdb( pdbname + ".solv.pdb", *scorefxn );
	}
*/

	// do pose_metric_calculators?
	basic::MetricValue< Size > all_bur_unsat_polars;
	basic::MetricValue< Size > all_solv_unsat_polars;
//	basic::MetricValue< Size > all_exsolv_unsat_polars;
	basic::MetricValue< vector1< Size > > res_bur_unsat_polars;
//	basic::MetricValue< vector1< Size > > res_exsolv_unsat_polars;
	basic::MetricValue< vector1< Size > > res_solv_unsat_polars;
	basic::MetricValue< vector1< Real > > res_semiex_score;
	basic::MetricValue<id::AtomID_Map< bool > > atom_bur_unsat_polars;
	basic::MetricValue<id::AtomID_Map< bool > > atom_semiex_unsat_polars;
//	basic::MetricValue<id::AtomID_Map< bool > > atom_exsolv_unsat_polars;
	basic::MetricValue<id::AtomID_Map< Real > > atom_semiex_score;
	basic::MetricValue<id::AtomID_Map< Real > > atom_sasa;
	basic::MetricValue<id::AtomID_Map< Size > > atom_hbonds;
//	basic::MetricValue<id::AtomID_Map< Size > > atom_water_hbonds;

	pose.metric( "num_hbonds", "atom_Hbonds", atom_hbonds );
//	pose.metric( "num_water_hbonds", "atom_Hbonds", atom_water_hbonds );
	pose.metric( "sasa", "atom_sasa", atom_sasa );

	pose.metric( "bur_unsat_polars", "all_bur_unsat_polars", all_bur_unsat_polars );
	pose.metric( "bur_unsat_polars", "residue_bur_unsat_polars", res_bur_unsat_polars );
	pose.metric( "bur_unsat_polars", "atom_bur_unsat", atom_bur_unsat_polars );

//	pose.metric( "exsolv_unsat_polars", "all_unsat_polars", all_exsolv_unsat_polars );
//	pose.metric( "exsolv_unsat_polars", "residue_unsat_polars", res_exsolv_unsat_polars );
//	pose.metric( "exsolv_unsat_polars", "atom_unsat", atom_exsolv_unsat_polars );

	if( option[ byres_data::solv_unsat_calc ] ){
		pose.metric( "solv_unsat_polars", "all_unsat_polars", all_solv_unsat_polars );
		pose.metric( "solv_unsat_polars", "residue_unsat_polars", res_solv_unsat_polars );
		pose.metric( "solv_unsat_polars", "atom_unsat", atom_semiex_unsat_polars );
		pose.metric( "solv_unsat_polars", "residue_semiex_score", res_semiex_score );
		pose.metric( "solv_unsat_polars", "atom_semiex_score", atom_semiex_score );
	}

	io::silent::SilentStructOP pose_ss( new io::silent::ScoreFileSilentStruct );
	core::scoring::hbonds::HBondSet hbset;
	hbset.setup_for_residue_pair_energies( pose, false, false );
	//each residue
	for( Size seqpos = 1; seqpos <= nres; ++seqpos ){
		if( packer_task->being_packed( seqpos ) ){
//			pose = start_pose;
			if( !option[ packing::resfile ].user() && !pose.residue( seqpos ).is_protein() && !pose.residue( seqpos ).is_DNA() ) continue;

			// WARNING TMP HACK BAIL IF NO NATIVE RES //
			if( seqpos > native_seqpos_map.size() ) continue;
			// WARNING TMP HACK BAIL IF NO NATIVE RES //

			//is the native seqpos the same, or coming from a seq alignment?
			Size native_seqpos( seqpos );
			if( option[ byres_data::align_native_seq ] ) native_seqpos = native_seqpos_map[ seqpos ];

			// WARNING TMP HACK BAIL IF NO NATIVE RES //
			if( native_seqpos < 1 ) continue;
			// WARNING TMP HACK BAIL IF NO NATIVE RES //

			io::silent::SilentStructOP res_ss( new io::silent::ScoreFileSilentStruct );
			res_ss->decoy_tag( pdbname );

			//print byatom data
			//only over non-water residues
			io::silent::SilentStructOP atom_ss( new io::silent::ScoreFileSilentStruct );
			Residue const rsd( pose.residue( seqpos ) );
			for( Size iatom = 1; iatom <= rsd.natoms(); ++iatom ){
				id::AtomID const id( iatom, seqpos );
				if( !( rsd.atom_type( iatom ).is_acceptor() || rsd.atom_type( iatom).is_donor() ) ) continue;
				std::string atom_name( rsd.atom_name( iatom ) );
				atom_name.erase( 0, atom_name.find_first_not_of( " " ) );
				atom_name.erase( atom_name.find_last_not_of( " " ) + 1 );
				atom_ss->decoy_tag( rsd.name3() + "_" + string_of( seqpos ) +
						"_" + atom_name );
				atom_ss->add_string_value( "atom_type", rsd.atom_type( iatom ).name() );
				atom_ss->add_energy( "seqpos", static_cast< Real >( seqpos ) );
				atom_ss->add_string_value( "pdbname", pdbname );
				atom_ss->add_string_value( "pdb_chain", utility::string_split( pose.pdb_info()->pose2pdb( seqpos ) )[ 1 ] );
				atom_ss->add_string_value( "pdb_resnum", utility::string_split( pose.pdb_info()->pose2pdb( seqpos ) )[ 2 ] );
				//			TR_unsat << "seqpos\t" << seqpos << "\tatom\t" << rsd.atom_name( iatom ) << "\t";
				//			TR_unsat << "wat_hbonds\t" << atom_water_hbonds.value()[ id ] << "\t";
				atom_ss->add_energy( "hbonds", static_cast< Real >( atom_hbonds.value()[ id ] ) );
//				atom_ss->add_energy( "wat_hbonds", static_cast< Real >( atom_water_hbonds.value()[ id ] ) );
				//need to sum attached H's for donors
				Real this_sasa( atom_sasa.value()[ id ] );
				for( Size hcount = rsd.type().attached_H_begin( iatom );
						hcount<= rsd.type().attached_H_end( iatom ); hcount++){
					this_sasa += atom_sasa.value()[ core::id::AtomID ( hcount, seqpos ) ];
				}
				atom_ss->add_energy( "sasa", this_sasa );
				Real atom_lk_burial( get_atom_lk_burial( pose, scorefxn, seqpos, iatom ) );
				atom_ss->add_energy( "atom_lk_burial", atom_lk_burial );
				atom_ss->add_energy( "bur_unsat", static_cast< Real >( atom_bur_unsat_polars.value()[ id ] ) );
//				atom_ss->add_energy( "exsolv_unsat", static_cast< Real >( atom_exsolv_unsat_polars.value()[ id ] ) );

				//iter thru atom's hbonds
				Real hbond_prot_prot( 0.0 );
				Real hbond_prot_water( 0.0 );
				Real hbond_prot_water_ref( 0.0 );
				Real hbond_water_ref( -0.353 );
				if( rsd.atom_type( iatom ).is_acceptor() ){
					vector1< hbonds::HBondCOP > const atom_hbset( hbset.atom_hbonds( id ) );
					for( Size ihb = 1; ihb <= atom_hbset.size(); ++ihb ){
						scoring::hbonds::HBond hb( *atom_hbset[ ihb ] );
						Size dres( hb.don_res() );
						Size ares( hb.acc_res() );
						Real hbe( hb.energy() * hb.weight() );
						if( pose.residue( dres ).is_protein() && pose.residue( ares ).is_protein() ){
							hbond_prot_prot += hbe;
						} else if( pose.residue( dres ).is_protein() && pose.residue( ares ).name1() == 'w' ||
								pose.residue( ares ).is_protein() && pose.residue( dres ).name1() == 'w' ){
							hbond_prot_water += hbe;
							hbond_prot_water_ref += ( hbe - hbond_water_ref );
						}
					}
				}
				//must sum up attached H's if donor
				if( rsd.atom_type( iatom ).is_donor() ){
					for( Size ihatom = rsd.type().attached_H_begin( iatom );
							ihatom <= rsd.type().attached_H_end( iatom ); ihatom++){
						id::AtomID const h_id( ihatom, seqpos );
						vector1< hbonds::HBondCOP > const atom_hbset( hbset.atom_hbonds( h_id ) );
						for( Size ihb = 1; ihb <= atom_hbset.size(); ++ihb ){
							scoring::hbonds::HBond hb( *atom_hbset[ ihb ] );
							Size dres( hb.don_res() );
							Size ares( hb.acc_res() );
							Real hbe( hb.energy() * hb.weight() );
							if( pose.residue( dres ).is_protein() && pose.residue( ares ).is_protein() ){
								hbond_prot_prot += hbe;
							} else if( pose.residue( dres ).is_protein() && pose.residue( ares ).name1() == 'w' ||
									pose.residue( ares ).is_protein() && pose.residue( dres ).name1() == 'w' ){
								hbond_prot_water += hbe;
								hbond_prot_water_ref += ( hbe - hbond_water_ref );
							}
						}
					}
				}

				Real hbond_total_ref( hbond_prot_prot + hbond_prot_water_ref );
				atom_ss->add_energy( "hbond_prot-prot", hbond_prot_prot );
//				atom_ss->add_energy( "hbond_prot-water", hbond_prot_water );
//				atom_ss->add_energy( "hbond_prot-water-ref", hbond_prot_water_ref );
				atom_ss->add_energy( "hbond_total-ref", hbond_total_ref );
				if( option[ byres_data::solv_unsat_calc ] ){
					atom_ss->add_energy( "semiex_unsat", static_cast< Real >( atom_semiex_unsat_polars.value()[ id ] ) );
					atom_ss->add_energy( "semiex_score", atom_semiex_score.value()[ id ] );
					atom_ss->add_energy( "hbond_total-semiex", ( atom_semiex_score.value()[ id ] + hbond_prot_prot ) );
				}


				//write scorefile line
				atom_silent_data.write_silent_struct( *atom_ss, scname_byatom );

				//accumulate atom data into res data
				vector1< io::silent::SilentEnergy > atom_ss_data( atom_ss->energies() );
				for( Size idata = 1; idata <= atom_ss_data.size(); ++idata ){
					std::string dataname( atom_ss_data[ idata ].name() );
					Real atom_dataval( atom_ss_data[ idata ].value() );
					Real res_dataval( atom_dataval + res_ss->get_energy( dataname ) );
					res_ss->add_energy( dataname, res_dataval );
				}
			}

			//get a bunch of data about residue
			if( !option[ byres_data::calcs_only ] ){
//				get_res_data_ss( res_ss, pose, native_pose, native_seqpos, seqpos, scorefxn, scorefxn_edens, calc_rmsd, residue_sasa );
				get_res_data_ss( res_ss, pose, native_pose, native_seqpos, seqpos, scorefxn, scorefxn_edens, calc_rmsd );
			}

			//now calculator vals
			res_ss->add_energy( "bur_unsat_polars", res_bur_unsat_polars.value()[ seqpos ] );
//			res_ss->add_energy( "exsolv_unsat_polars", res_exsolv_unsat_polars.value()[ seqpos ] );
			if( option[ byres_data::solv_unsat_calc ] ){
					res_ss->add_energy( "solv_unsat_polars", res_solv_unsat_polars.value()[ seqpos ] );
			}
/*
			res_ss->add_energy( "bur_unsat_hbonds", res_bur_unsat_hbonds.value()[ seqpos ] );
			res_ss->add_energy( "exsolv_unsat_hbonds", res_exsolv_unsat_hbonds.value()[ seqpos ] );
			if( option[ byres_data::solv_unsat_calc ] ){
					res_ss->add_energy( "solv_unsat_hbonds", res_solv_unsat_hbonds.value()[ seqpos ] );
			}
*/
			//write scorefile line
			res_silent_data.write_silent_struct( *res_ss, scname_byres );

			//accumulate res data into pose data
			vector1< io::silent::SilentEnergy > res_ss_data( res_ss->energies() );
			for( Size idata = 1; idata <= res_ss_data.size(); ++idata ){
				std::string dataname( res_ss_data[ idata ].name() );
				Real res_dataval( res_ss_data[ idata ].value() );
				Real pose_dataval( res_dataval + pose_ss->get_energy( dataname ) );
				pose_ss->add_energy( dataname, pose_dataval );
			}
		}
	}

	//do correction for hbond count
	if( !option[ byres_data::calcs_only ] ){
		Real n_bb_hbonds( pose_ss->get_energy( "n_bb_hbonds" ) / 2.0 );
		Real n_sc_hbonds( pose_ss->get_energy( "n_sc_hbonds" ) / 2.0 );
		Real n_bb_wat_hbonds( pose_ss->get_energy( "n_bb_wat_hbonds" ) );
		Real n_sc_wat_hbonds( pose_ss->get_energy( "n_sc_wat_hbonds" ) );
		Real n_hbonds( n_bb_hbonds + n_sc_hbonds + n_bb_wat_hbonds + n_sc_wat_hbonds );
		pose_ss->add_energy( "n_bb_hbonds", n_bb_hbonds );
		pose_ss->add_energy( "n_sc_hbonds", n_sc_hbonds );
		pose_ss->add_energy( "n_hbonds", n_hbonds );
	}

	pose_ss->decoy_tag( pdbname );
	pose_silent_data.write_silent_struct( *pose_ss, scname );
}

void
load_coords()
{
	using namespace core;

	import_pose::pose_stream::MetaPoseInputStream input =  import_pose::pose_stream::streams_from_cmd_line();
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	core::pose::Pose pose;
	while ( input.has_another_pose() ) {
		input.fill_pose( pose, *rsd_set );
		for( Size i = 1; i<= option[ byres_data::repeat ]; ++i ){
			byres_analysis( pose );
		}
	}
}

void
register_calcs()
{
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );
//	core::pose::metrics::PoseMetricCalculatorOP num_water_hbonds_calculator = new protocols::toolbox::pose_metric_calculators::NumberWaterHBondsCalculator();
//	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "num_water_hbonds", num_water_hbonds_calculator );
	core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculatorLegacy();
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );

	core::pose::metrics::PoseMetricCalculatorOP bur_unsat_calculator = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator( "sasa", "num_hbonds" );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "bur_unsat_polars", bur_unsat_calculator );
//	core::pose::metrics::PoseMetricCalculatorOP exsolv_unsat_calculator = new protocols::toolbox::pose_metric_calculators::ExplicitWaterUnsatisfiedPolarsCalculator( scorefxn );
//	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "exsolv_unsat_polars", exsolv_unsat_calculator );
	if( option[ byres_data::solv_unsat_calc ] ){
		core::pose::metrics::PoseMetricCalculatorOP solv_unsat_calculator = new protocols::toolbox::pose_metric_calculators::SemiExplicitWaterUnsatisfiedPolarsCalculator("num_hbonds", scorefxn );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "solv_unsat_polars", solv_unsat_calculator );
	}
/*
	core::pose::metrics::PoseMetricCalculatorOP bur_unsat_hb_calculator = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedHBondsCalculator( "sasa", "num_hbonds" );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "bur_unsat_hbonds", bur_unsat_hb_calculator );
	core::pose::metrics::PoseMetricCalculatorOP exsolv_unsat_hb_calculator = new protocols::toolbox::pose_metric_calculators::ExplicitWaterUnsatisfiedHBondsCalculator("num_hbonds" );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "exsolv_unsat_hbonds", exsolv_unsat_hb_calculator );
	if( option[ byres_data::solv_unsat_calc ] ){
		core::pose::metrics::PoseMetricCalculatorOP solv_unsat_hb_calculator = new protocols::toolbox::pose_metric_calculators::SemiExplicitWaterUnsatisfiedHBondsCalculator("num_hbonds", scorefxn );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "solv_unsat_hbonds", solv_unsat_hb_calculator );
	}
*/
}

void*
my_main( void*)
{
	register_calcs();
	load_coords();
	exit(0);

}

int
main( int argc, char * argv [] )
{
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility;

	option.add( byres_data::solvate, "solvate the whole protein" ).def( false );
	option.add( byres_data::calcs_only, "pose calcs only" ).def( false );
	option.add( byres_data::solv_unsat_calc, "add explicitwater pose metric calc" ).def( false );
//	option.add( byres_data::calc_sasa, "calc byatom and byres sasa" ).def( false );
	option.add( byres_data::tag, "nametag" ).def( "score" );
	option.add( byres_data::repeat, "rescore structure" ).def( 1 );
	option.add( byres_data::nloop_solvadd, "number of iter of append water molecule per atom" ).def( 10 );
	option.add( byres_data::nloop_solvdock, "number of iter of RBMover water molecule per atom" ).def( 10 );
	option.add( byres_data::nloop_hbscan, "number of iter of hbscan per polar h or acc" ).def( 10 );
	option.add( byres_data::align_native_seq, "do a sequence alignment to map native seqpos to pose seqpos" ).def( false );

	devel::init(argc, argv);

	protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


