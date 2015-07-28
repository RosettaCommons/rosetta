// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   cov_hbs.cc
/// @brief  Sidechain conjugation to acryl amides
/// @author Andy Watkins (amw579@nyu.edu)

// includes
#include <iostream>
#include <fstream>
#include <string>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/ncbb/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/func/SumFunc.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Numeric Headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::a3b_hbs;

using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR("Vancomycin");

void add_constraints(
	core::pose::Pose & vancomycin
) {

	core::id::AtomID const CZ2_2( vancomycin.residue(2).type().atom_index( "CZ2" ), 2 );
	core::id::AtomID const OT_2( vancomycin.residue(2).type().atom_index( "OT" ), 2 );
	core::id::AtomID const CZ2_6( vancomycin.residue(6).type().atom_index( "CZ2" ), 6 );
	core::id::AtomID const OT_6( vancomycin.residue(6).type().atom_index( "OT" ), 6 );
	core::id::AtomID const CD1( vancomycin.residue(4).type().atom_index( "CD1" ), 4 );
	core::id::AtomID const CE ( vancomycin.residue(4).type().atom_index( "CE"  ), 4 );
	core::id::AtomID const CG1( vancomycin.residue(4).type().atom_index( "CG1" ), 4 );
	core::id::AtomID const CG2( vancomycin.residue(4).type().atom_index( "CG2" ), 4 );
	core::id::AtomID const CD2( vancomycin.residue(4).type().atom_index( "CD2" ), 4 );
	core::id::AtomID const CD12( vancomycin.residue(5).type().atom_index( "CD1" ), 5 );
	core::id::AtomID const CG12( vancomycin.residue(7).type().atom_index( "CG1" ), 7 );
	core::id::AtomID const CB ( vancomycin.residue(7).type().atom_index( "CB" ), 7 );
	core::id::AtomID const CD13( vancomycin.residue(7).type().atom_index( "CD1" ), 7 );
	
	core::scoring::func::HarmonicFuncOP harm( new core::scoring::func::HarmonicFunc( 1.31, 0.02 ) );
	core::scoring::constraints::AtomPairConstraintOP bond1( new core::scoring::constraints::AtomPairConstraint( OT_2, CD1, harm ) );
	core::scoring::constraints::AtomPairConstraintOP bond2( new core::scoring::constraints::AtomPairConstraint( OT_6, CD2, harm ) );
	core::scoring::constraints::AtomPairConstraintOP bond3( new core::scoring::constraints::AtomPairConstraint( CD12, CG12, harm ) );
	
	vancomycin.add_constraint( bond1 );
	vancomycin.add_constraint( bond2 );
	vancomycin.add_constraint( bond3 );
	
	core::scoring::func::CircularHarmonicFuncOP aharm(  new core::scoring::func::CircularHarmonicFunc( 104.5*3.14159/180, 0.02 ) );
	core::scoring::func::CircularHarmonicFuncOP aharm2( new core::scoring::func::CircularHarmonicFunc( 120  *3.14159/180, 0.02 ) );
	
	core::scoring::constraints::AngleConstraintOP ang1( new core::scoring::constraints::AngleConstraint( CZ2_2, OT_2, CD1, aharm2 ) );
	core::scoring::constraints::AngleConstraintOP ang2( new core::scoring::constraints::AngleConstraint( CZ2_6, OT_6, CD2, aharm2 ) );
	core::scoring::constraints::AngleConstraintOP ang3( new core::scoring::constraints::AngleConstraint( OT_2, CD1, CE, aharm2 ) );
	core::scoring::constraints::AngleConstraintOP ang4( new core::scoring::constraints::AngleConstraint( OT_6, CD2, CE, aharm2 ) );
	core::scoring::constraints::AngleConstraintOP ang5( new core::scoring::constraints::AngleConstraint( CD12, CG12, CD13, aharm2 ) );
	
	vancomycin.add_constraint( ang1 );
	vancomycin.add_constraint( ang2 );
	vancomycin.add_constraint( ang3 );
	vancomycin.add_constraint( ang4 );
	vancomycin.add_constraint( ang5 );
	
	core::scoring::func::CircularHarmonicFuncOP charm( new core::scoring::func::CircularHarmonicFunc( 3.14, 0.02 ) );
	core::scoring::constraints::DihedralConstraintOP dih1( new core::scoring::constraints::DihedralConstraint( OT_2, CD1, CG1, CE, charm ) );
	core::scoring::constraints::DihedralConstraintOP dih2( new core::scoring::constraints::DihedralConstraint( OT_6, CD2, CG2, CE, charm ) );
	core::scoring::constraints::DihedralConstraintOP dih3( new core::scoring::constraints::DihedralConstraint( CD12, CG12, CB, CD13, charm ) );
	
	vancomycin.add_constraint( dih1 );
	vancomycin.add_constraint( dih2 );
	vancomycin.add_constraint( dih3 );
}

int main ( int argc, char* argv[] )
{
try {
	using namespace core;
	using namespace utility;
	using namespace scoring;
	using namespace pose;
	using namespace core::chemical;
	using namespace conformation;
	using namespace func;
	using namespace constraints;
	
	using namespace core::id;
	using namespace core::pack;
	using namespace core::pack::task;
	//option[ chemical::patch_selectors ].push_back( "CTERM_AMIDATION" );
		
	devel::init(argc, argv);
	
	ScoreFunctionOP scorefxn = get_score_function();
	if ( scorefxn->get_weight( atom_pair_constraint ) == 0.0 ) {
		scorefxn->set_weight( atom_pair_constraint, 0.1 );
	}
	if ( scorefxn->get_weight( angle_constraint ) == 0.0 ) {
		scorefxn->set_weight( angle_constraint, 1.0 );
	}
	if ( scorefxn->get_weight( dihedral_constraint ) == 0.0 ) {
		scorefxn->set_weight( dihedral_constraint, 1.0 );
	}
	
	// Vanc sequence: N-methyl D-leu
	
	core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	//ResidueType const & dleu_type = residue_set_cap->name_map( "DLEU:NtermProteinMethylated" );
	ResidueType const & dleu_type = residue_set_cap->name_map( "DTRP:NtermProteinMethylated" );
	ResidueType const & V01_type = residue_set_cap->name_map( "V01:aryl-O-conjugated" );
	//ResidueType const & asn_type = residue_set_cap->name_map( "ASN" );
	ResidueType const & asn_type = residue_set_cap->name_map( "LYS" );
	ResidueType const & V02_type = residue_set_cap->name_map( "V02:aryl-O-conjugated:phg_cd1_conjugation:phg_cd2_conjugation" );
	ResidueType const & V02_type2 = residue_set_cap->name_map( "V02:phg_cd1_conjugation" );
	ResidueType const & V04_type = residue_set_cap->name_map( "V04:aryl-O-conjugated" );
	//ResidueType const & V03_type = residue_set_cap->name_map( "V03:CtermProteinFull:aryl-C-conjugated" );
	ResidueType const & V03_type = residue_set_cap->name_map( "V03:aryl-C-conjugated" );
	
	Residue dleu( dleu_type, true );
	Residue V01( V01_type, true );
	Residue asn( asn_type, true );
	Residue V02( V02_type, true );
	Residue V022( V02_type2, true );
	Residue V04( V04_type, true );
	Residue V03( V03_type, true );
	
	Pose vancomycin;
	vancomycin.append_residue_by_jump( dleu, 1 );
	vancomycin.append_residue_by_bond( V01, true );
	vancomycin.append_residue_by_bond( asn, true );
	vancomycin.append_residue_by_bond( V02, true );
	vancomycin.append_residue_by_bond( V022, true );
	vancomycin.append_residue_by_bond( V04, true );
	vancomycin.append_residue_by_bond( V03, true );
	// extra!
	vancomycin.append_residue_by_bond( *new Residue( residue_set_cap->name_map( "GLN:CtermProteinFull" ), true ), true );
	
	//Size v02_connid1 = V02_type.residue_connection_id_for_atom( pose.residue( 4 ).atom_index( "CD1" ) );
	//Size v02_connid2 = V02_type.residue_connection_id_for_atom( pose.residue( 4 ).atom_index( "CD2" ) );
	//new_cys->residue_connection_partner( cys_connid, resi_vdp, vdp_connid );
	
	vancomycin.conformation().declare_chemical_bond( 2, "OT", 4, "CD1" );
	vancomycin.conformation().declare_chemical_bond( 6, "OT", 4, "CD2" );
	vancomycin.conformation().declare_chemical_bond( 5, "CD1", 7, "CG1" );
	
	add_constraints( vancomycin );
	
	for ( Size ii = 1; ii <= vancomycin.total_residue(); ++ii ) {
		vancomycin.set_phi( ii, -150 );
		vancomycin.set_psi( ii, 150 );
		vancomycin.set_omega( ii, 180 );
	}
	
	vancomycin.dump_pdb( "vancomycin.pdb" );
	
	kinematics::MoveMapOP pert_mm( new kinematics::MoveMap() );
	for ( Size i = 1; i <= vancomycin.total_residue(); ++i ) {
		pert_mm->set_bb( i, true );
		pert_mm->set_chi( i, true );
	}
	protocols::simple_moves::MinMoverOP min_mover( new simple_moves::MinMover( pert_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true ) );
	
	for ( Real apc = 0; apc <= 2; apc += 0.1 ) {
		scorefxn->set_weight( atom_pair_constraint, apc );
		scorefxn->set_weight( dihedral_constraint, apc );
		
		min_mover->apply( vancomycin );
	}
	vancomycin.dump_pdb( "vancomycin_min1.pdb" );
	for ( Real apc = 2; apc >= 0; apc -= 0.1 ) {
		scorefxn->set_weight( atom_pair_constraint, apc );
		scorefxn->set_weight( dihedral_constraint, apc );
		
		min_mover->apply( vancomycin );
	}
	vancomycin.dump_pdb( "vancomycin_min2.pdb" );


	
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
