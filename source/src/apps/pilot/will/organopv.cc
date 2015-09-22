// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/coiled_coil.cc
/// @brief samples coiled coils

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <devel/init.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/optimizeH.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/ResDepAngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ResDepAtomPairConstraint.hh>
#include <core/scoring/constraints/ResDepDihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <basic/Tracer.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/flxbb/DesignLayerOperation.fwd.hh>
#include <protocols/flxbb/DesignLayerOperation.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/relax/SimpleMultiRelax.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/mynamespaces.ihh>


static THREAD_LOCAL basic::Tracer TR( "organopv" );


void add_apc(core::pose::Pose & pose, core::id::AtomID aid1, core::id::AtomID aid2, core::Real mean, core::Real sd,
			 core::scoring::ScoreType st = core::scoring::atom_pair_constraint)
{
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::AtomPairConstraint(
		aid1, aid2, new core::scoring::constraints::HarmonicFunc(mean,sd), st );
	pose.add_constraint(cc);
}


void
add_agc(
	core::pose::Pose & pose,
	core::id::AtomID aid1,
	core::id::AtomID aid2,
	core::id::AtomID aid3,
	core::Real mean,
	core::Real sd
) {
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::AngleConstraint(
		aid1, aid2, aid3, new core::scoring::constraints::HarmonicFunc(mean,sd) );
	pose.add_constraint(cc);
}


void minimize(Pose & pose, ScoreFunctionOP sf) {
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_chi(true);
	movemap->set_bb(false);
	movemap->set_jump(true);
	core::conformation::symmetry::make_symmetric_movemap( pose, *movemap );
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
	m.apply(pose);
}

// void moverb(Pose & pose, ScoreFunctionOP sf)

void design(Pose & pose, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	vector1< bool > aas(20,false);
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	aas[core::chemical::aa_ala] = true;
	aas[core::chemical::aa_asp] = true;
	aas[core::chemical::aa_asn] = true;
	aas[core::chemical::aa_gln] = true;
	aas[core::chemical::aa_glu] = true;
	aas[core::chemical::aa_ile] = true;
	aas[core::chemical::aa_leu] = true;
	aas[core::chemical::aa_lys] = true;
	aas[core::chemical::aa_ser] = true;
	aas[core::chemical::aa_thr] = true;
	aas[core::chemical::aa_tyr] = true;
	aas[core::chemical::aa_val] = true;
	task->nonconst_residue_task(1).restrict_absent_canonical_aas(aas);
	task->nonconst_residue_task(2).restrict_absent_canonical_aas(aas);
	task->nonconst_residue_task(3).restrict_to_repacking();
	task->nonconst_residue_task(4).restrict_to_repacking();
	// task->nonconst_residue_task(2).prevent_repacking();
	// task->nonconst_residue_task(3).prevent_repacking();
	// task->nonconst_residue_task(4).prevent_repacking();
	task->or_include_current(false);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);


}


void* doit(void* /*x = NULL*/) {
	using namespace core;
	using namespace pose;
	using namespace protocols;
	using namespace moves;
	using namespace ObjexxFCL::format;
	using numeric::random::uniform;


	Pose pose;
	// chemical::ResidueTypeSetCAP rs( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
	// chemical::make_pose_from_sequence(pose,std::string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),*rs,false);
	// for(Size i = 1; i <= pose.n_residue(); ++i) {
	// 	pose.set_phi  (i,-60); pose.set_psi  (i,-45); pose.set_omega(i,180);
	// }
	// pose.dump_pdb("helix.pdb");
	// std::exit(-1);

	// io::pdb::pose_from_pdb(pose,"params/pyr_0001.pdb");
	// chemical::remove_lower_terminus_type_from_pose_residue(pose,1);
	// chemical::remove_upper_terminus_type_from_pose_residue(pose,1);
	// for(Size i = 1; i <= 36; ++i) {
	// 	pose.set_chi(1,1,10.0*i);
	// 	pose.dump_pdb("pyr_chi_"+string_of(i)+".pdb");
	// }
	// std::exit(-1);


	io::pdb::pose_from_pdb(pose,"input/helix_pyr.pdb");
	chemical::remove_lower_terminus_type_from_pose_residue(pose,1);
	chemical::remove_upper_terminus_type_from_pose_residue(pose,2);
	// chemical::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, 2 );
	// chemical::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, 1 );
	// chemical::add_variant_type_to_pose_residue( pose, "VIRTUAL_N", 1 );
	// chemical::add_variant_type_to_pose_residue( pose, "VIRTUAL_C", 2 );
	pose.set_phi  (1,-58.4193);
	pose.set_psi  (1,-46.555);
	pose.set_omega(1,-179.774);
	pose.set_phi  (2,-58.4193);
	pose.set_psi  (2,-46.555);
	pose.set_omega(2,-179.774);
	pose.set_xyz(AtomID(8,1),numeric::xyzVector<Real>(1.404,0.540,-3.183));

	// for(Size i=1; i<= pose.residue(1).natoms(); ++i) {
	// 	std::cerr << "ATOM " << i << " " << pose.residue(1).atom_name(i) << std::endl;
	// }
	// for(Size i=1; i<= pose.residue(2).natoms(); ++i) {
	// 	std::cerr << "ATOM " << i << " " << pose.residue(2).atom_name(i) << std::endl;
	// }
	// for(Size i=1; i<= pose.residue(3).natoms(); ++i) {
	// 	std::cerr << "ATOM " << i << " " << pose.residue(3).atom_name(i) << std::endl;
	// }
	// std::exit(-1);


	core::conformation::symmetry::make_symmetric_pose( pose );

	// pose.dump_pdb("sym_test.pdb");
	// std::exit(-1);

	// std::cerr << "TMP " << pose.residue(3).atom_index("C1") << std::endl;
	using namespace core::scoring::constraints;
	std::map<std::string,std::string> atomRO,atomLO,atomL,atomR,atomLO2;
	atomRO["ALA"]= "CB";	atomR["ALA"]= "CA";
	atomRO["ASN"]="OD1";	atomR["ASN"]= "CG";
	atomRO["ASP"]="OD1";	atomR["ASP"]= "CG";
	atomRO["GLN"]="OE1";	atomR["GLN"]= "CD";
	atomRO["GLU"]="OE1";	atomR["GLU"]= "CD";
	atomRO["ILE"]="CD1";	atomR["ILE"]="CG1";
	atomRO["LEU"]="CD2";	atomR["LEU"]= "CG";
	atomRO["LYS"]= "NZ";	atomR["LYS"]= "CE";
	atomRO["SER"]= "OG";	atomR["SER"]= "CB";
	atomRO["THR"]="OG1";	atomR["THR"]= "CB";
	atomRO["TYR"]= "OH";	atomR["TYR"]= "CZ";
	atomRO["VAL"]="CG2";	atomR["VAL"]= "CB";
	atomLO ["PYR"] = "B";
	atomL  ["PYR"] = "C";
	atomLO2["PYR"] = "P";

	pose.add_constraint(new ResDepAtomPairConstraint( atomRO,atomLO,1,3,new HarmonicFunc(0.0,0.03) ));
	pose.add_constraint(new ResDepAtomPairConstraint( atomRO,atomLO,2,4,new HarmonicFunc(0.0,0.03) ));
	pose.add_constraint(new ResDepAngleConstraint(atomR,atomLO,atomL,1,3,3, new HarmonicFunc(1.884956,0.003) ));
	pose.add_constraint(new ResDepAngleConstraint(atomR,atomLO,atomL,2,4,4, new HarmonicFunc(1.884956,0.003) ));
	// pose.add_constraint(new ResDepDihedralConstraint(atomR,atomLO,atomL,atomLO2,1,3,3,3, new HarmonicFunc(0.0,0.003) ));
	// pose.add_constraint(new ResDepDihedralConstraint(atomR,atomLO,atomL,atomLO2,2,4,4,4, new HarmonicFunc(0.0,0.003) ));

	pose.add_constraint(new AtomPairConstraint( AtomID(2,3), AtomID(3,4 ), new HarmonicFunc(2.8,0.08) ));
	pose.add_constraint(new AtomPairConstraint( AtomID(2,4), AtomID(3,11), new HarmonicFunc(2.8,0.08) ));
	pose.add_constraint(new AngleConstraint( AtomID(1,3), AtomID(2,3), AtomID(3,4 ), new HarmonicFunc(2.094395,0.03) ));
	pose.add_constraint(new AngleConstraint( AtomID(1,4), AtomID(2,4), AtomID(3,11), new HarmonicFunc(2.094395,0.03) ));

	// hbonding
	pose.add_constraint(new ResDepAngleConstraint(atomR,atomLO,atomLO2,1,3, 4, new HarmonicFunc(2.094395,0.03) ));
	pose.add_constraint(new ResDepAngleConstraint(atomR,atomLO,atomLO2,2,4,11, new HarmonicFunc(2.094395,0.03) ));

	pose.add_constraint(new AngleConstraint( AtomID(2,3), AtomID(3,3), AtomID(1,4 ), new HarmonicFunc(2.094395,0.03) ));
	pose.add_constraint(new AngleConstraint( AtomID(2,4), AtomID(3,4), AtomID(1,11), new HarmonicFunc(2.094395,0.03) ));

	for(Size i = 1; i <= 16; i++) {
		// std::cerr << "ATOM " << i << " " << pose.residue(3).atom_name(i) << std::endl;
		add_apc(pose,AtomID(i,3),AtomID(i, 4),3.0,1.0);
		add_apc(pose,AtomID(i,4),AtomID(i,11),3.0,1.0);
	}

	ScoreFunctionOP sf;
	sf = core::scoring::get_score_function();
	// sf = new core::scoring::ScoreFunction();
	// sf->set_weight(core::scoring::chainbreak,100.00);
	sf->set_weight(core::scoring::atom_pair_constraint,2.00);
	sf->set_weight(core::scoring::angle_constraint,2.00);
	sf->set_weight(core::scoring::dihedral_constraint,2.00);
	sf->set_weight(core::scoring::fa_sol,0.00);
	sf = new core::scoring::symmetry::SymmetricScoreFunction(*sf);

	sf->show(pose);

	std::string tag = string_of(uniform());

	if (option[ basic::options::OptionKeys::parser::view ]()) {
		protocols::viewer::add_conformation_viewer(pose.conformation(),"smhybrid",1200,1200);
	}

	Pose min = pose;
	for(Size i = 1; i < 9999; ++i) {
		for(Size k = 1; k < 9; ++k) {
			for(Size j = 1; j <= pose.fold_tree().num_jump(); ++j) {
				if( !( pose.fold_tree().downstream_jump_residue(j)==3 || pose.fold_tree().downstream_jump_residue(j)==4 ) ) continue;
				core::kinematics::Jump jmp(pose.jump(j));
				jmp.set_rotation(numeric::rotation_matrix_degrees(numeric::xyzVector<Real>(1,0,0),gaussian())*jmp.get_rotation());
				jmp.set_rotation(numeric::rotation_matrix_degrees(numeric::xyzVector<Real>(0,1,0),gaussian())*jmp.get_rotation());
				jmp.set_rotation(numeric::rotation_matrix_degrees(numeric::xyzVector<Real>(0,0,1),gaussian())*jmp.get_rotation());
				jmp.set_translation(jmp.get_translation()+(gaussian()/3.0)*numeric::xyzVector<Real>(1,0,0));
				jmp.set_translation(jmp.get_translation()+(gaussian()/3.0)*numeric::xyzVector<Real>(0,1,0));
				jmp.set_translation(jmp.get_translation()+(gaussian()/3.0)*numeric::xyzVector<Real>(0,0,1));
				pose.set_jump(j,jmp);
			}
			design(pose,sf);
			minimize(pose,sf);
			if( (*sf)(min) > (*sf)(pose) ) min = pose;
		}
		// pose = min;
	}
	std::cerr << "SCORE " << (*sf)(pose) << " " << tag << std::endl;
	sf->show(pose);

	pose.dump_pdb("organopv"+tag+".pdb");

	// std::exit(-1);
	return NULL;
}


int
main( int argc, char * argv [] )
{

	try {

	devel::init(argc,argv);

	if (option[ basic::options::OptionKeys::parser::view ]()) {
		protocols::viewer::viewer_main( &doit );
	} else {
		while(1) doit(NULL);
	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
