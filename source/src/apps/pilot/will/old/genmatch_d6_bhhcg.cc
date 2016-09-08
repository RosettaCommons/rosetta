// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <boost/tuple/tuple.hpp>
#include <basic/database/open.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
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
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static THREAD_LOCAL basic::Tracer TR( "genmatch_d6_bpy" );


void myoptH(Pose & pose, ScoreFunctionOP sf) {
	add_lower_terminus_type_to_pose_residue(pose,1);
	add_upper_terminus_type_to_pose_residue(pose,pose.size());
	core::pack::optimizeH(pose,*sf);
	remove_lower_terminus_type_from_pose_residue(pose,1);
	remove_upper_terminus_type_from_pose_residue(pose,pose.size());
}


// find 2 HIS that can chelate a tetrahedral metal
void run() {
	using namespace basic::options::OptionKeys;
	using namespace core::id;

	core::chemical::ResidueTypeSetCAP cen_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
	core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Real tmpdis = option[willmatch::max_dis_metal]();
	Real const MXDSMTL = tmpdis*tmpdis;
	Real const MXAGMTL = option[willmatch::max_ang_metal]();
	Real const MATCH_OVERLAP_DOT = cos(numeric::conversions::radians(option[willmatch::match_overlap_ang]()));

	vector1<string> infiles;
	if( option[in::file::l].user() ) {
		utility::io::izstream in(option[in::file::l]()[1]);
		string tmp;
		while(in >> tmp) infiles.push_back(tmp);
	} else if(option[in::file::s].user()) {
		infiles = option[in::file::s]();
	} else {
		utility_exit_with_message("no input!");
	}

	for(Size ifile = 1; ifile <= infiles.size(); ifile++) {
		string infile = infiles[ifile];
		Pose in_cen,in_fa;
		pose_from_file(in_fa, *fa_residue_set,infile, core::import_pose::PDB_file);
		Pose native = in_fa;
		pose_from_file(in_cen,*cen_residue_set,infile, core::import_pose::PDB_file);
		Size nres = in_cen.size();
		core::chemical::ResidueType const & ala( in_cen.residue(1).residue_type_set().name_map("ALA") );
		core::chemical::ResidueType const & alafa( in_fa.residue(1).residue_type_set().name_map("ALA") );
		// core::chemical::ResidueType const & hise( in_fa.residue(1).residue_type_set().name_map("HIS") );
		// core::chemical::ResidueType const & hisd( in_fa.residue(1).residue_type_set().name_map("HIS_D") );
		for(Size i = 1; i <= nres; ++i) {
			core::pose::replace_pose_residue_copying_existing_coordinates(in_cen,i,ala);
			core::pose::replace_pose_residue_copying_existing_coordinates(in_fa,i,alafa);
		}
		Pose const init_pose = in_cen;
		Pose const fa_pose = in_fa;
		ImplicitFastClashCheck clashcheck(init_pose,basic::options::option[basic::options::OptionKeys::willmatch::clash_dis]());

		ScoreFunctionOP sf = core::scoring::get_score_function();

		Real chi1incr = option[willmatch::chi1_increment]();
		Real chi2incr = option[willmatch::chi2_increment]();
		vector1<Real> CHI1,CHI2;
		for(Real i = 0; i < 360; i+= chi1incr) CHI1.push_back(i);
		for(Real i = 0; i < 360; i+= chi2incr) CHI2.push_back(i);

		// setup HIS residues for checking
		Pose hse,hsd,bpy,glu,cys;
		core::pose::make_pose_from_sequence(hse,"H[HIS]"  ,*fa_residue_set,false);
		core::pose::make_pose_from_sequence(hsd,"H[HIS_D]",*fa_residue_set,false);
		core::pose::make_pose_from_sequence(bpy,"X[BPY]"  ,*fa_residue_set,false);
		core::pose::make_pose_from_sequence(glu,"E"       ,*fa_residue_set,false);
		core::pose::make_pose_from_sequence(cys,"C"       ,*fa_residue_set,false);
		hsd.set_dof(DOF_ID(AtomID(hsd.residue(1).atom_index("HD1"),1),D),2.1);
		hse.set_dof(DOF_ID(AtomID(hse.residue(1).atom_index("HE2"),1),D),2.1);
		// hse.dump_pdb("hse.pdb");
		// hsd.dump_pdb("hsd.pdb");
		ObjexxFCL::FArray3D<Vec> chi2cen(2,CHI1.size(),CHI2.size()),chi2ori(2,CHI1.size(),CHI2.size());
		ObjexxFCL::FArray2D<Vec> bci2cen(  CHI1.size(),CHI2.size()),bci2ori(  CHI1.size(),CHI2.size()),bci2end(  CHI1.size(),CHI2.size());
		ObjexxFCL::FArray2D<Vec> gci2cen(  CHI1.size(),CHI2.size()),gci2ori(  CHI1.size(),CHI2.size()),gci2cd(  CHI1.size(),CHI2.size());
		ObjexxFCL::FArray2D<Vec> mcentroid(2,CHI1.size());
		Stub stub(hse.xyz(AtomID(5,1)),hse.xyz(AtomID(2,1)),hse.xyz(AtomID(1,1)));
		// precompute cen and ori for hisd and hise
		for(Size i = 1; i <= CHI1.size(); ++i) {
			hsd.set_chi(1,1,CHI1[i]);
			hse.set_chi(1,1,CHI1[i]);
			bpy.set_chi(1,1,CHI1[i]);
			glu.set_chi(1,1,CHI1[i]);
			mcentroid(1,i) = stub.global2local(  hse.residue(1).xyz("CG") );
			mcentroid(2,i) = stub.global2local( (hse.residue(1).xyz("CG")-hse.residue(1).xyz("CB")).normalized()*5.4 + hse.residue(1).xyz("CB") );
			for(Size j = 1; j <= CHI2.size(); ++j) {
				hse.set_chi(2,1,CHI2[j]);
				hsd.set_chi(2,1,CHI2[j]);
				bpy.set_chi(2,1,CHI2[j]);
				glu.set_chi(2,1,CHI2[j]);
				chi2cen(1,i,j) = stub.global2local( hsd.residue(1).xyz("HD1"));
				chi2ori(1,i,j) = stub.global2local( hsd.residue(1).xyz("ND1"));
				chi2cen(2,i,j) = stub.global2local( hse.residue(1).xyz("HE2"));
				chi2ori(2,i,j) = stub.global2local( hse.residue(1).xyz("NE2"));
				bci2cen(  i,j) = stub.global2local( bpy.residue(1).xyz("ZN"));
				bci2end(  i,j) = stub.global2local( bpy.residue(1).xyz("CO2"));
				Vec cd = glu.residue(1).xyz("CD");
				Vec cg = glu.residue(1).xyz("CG");
				Vec ori = (cg-cd).normalized();
				Vec cen = cd - 2.362957*ori;
				ori = cen + ori;
				gci2cen(  i,j) = stub.global2local( cen );
				gci2ori(  i,j) = stub.global2local( ori );
				gci2cd (  i,j) = stub.global2local( cd );
			}
		}

		// string fname = utility::file_basename(infile)+"_d6_bpy_matches.pdb";
		// utility::io::ozstream out(fname);
		// out << "MODEL BASE"	<< endl;
		// fa_pose.dump_pdb(out);
		// out << "ENDMDL" << endl;

		Pose pose = fa_pose;

		Size count = 0;

		vector1<Size> batm;
		TR << "BATM:";
		for(Size i = 6; i <= bpy.residue(1).nheavyatoms(); ++i) {
			// TR << "'" << bpy.residue(1).atom_name(i) << "'" << std::endl;
			if(bpy.residue(1).atom_name(i)==" NE1") { continue; }
			if(bpy.residue(1).atom_name(i)==" NN1") { continue; }
			if(bpy.residue(1).atom_name(i)==" ZN ") { continue; }
			batm.push_back(i);
			TR << " " << i;
		}
		TR << std::endl;
		vector1<Size> hatm;
		TR << "HATM:";
		for(Size i = 6; i <= hsd.residue(1).nheavyatoms(); ++i) {
			if(hsd.residue(1).atom_name(i)==" ND1") { continue; }
			if(hsd.residue(1).atom_name(i)==" NE2") { continue; }
			hatm.push_back(i);
			TR << " " << i;
		}
		TR << std::endl;
		vector1<Size> eatm;
		TR << "EATM:";
		for(Size i = 6; i <= glu.residue(1).nheavyatoms(); ++i) {
			if(glu.residue(1).atom_name(i)==" ND1") { continue; }
			if(glu.residue(1).atom_name(i)==" NE2") { continue; }
			eatm.push_back(i);
			TR << " " << i;
		}
		TR << std::endl;
		vector1<Size> scanres;
		if(option[willmatch::residues].user()) {
			TR << "input scanres!!!!!!" << std::endl;
			scanres = option[willmatch::residues]();
		} else {
			for(Size i = 2; i <= pose.size()-1; ++i) {
				if(!pose.residue(i).has("N" )) { continue; }
				if(!pose.residue(i).has("CA")) { continue; }
				if(!pose.residue(i).has("C" )) { continue; }
				if(!pose.residue(i).has("O" )) { continue; }
				if(!pose.residue(i).has("CB")) { continue; }
				if(pose.residue(i).name3()=="PRO") { continue; }
				scanres.push_back(i);
			}
		}

		vector1<Vec> foundcen,foundori1,foundori2,foundori3,foundori4,foundori5,foundori6;
		Size lastb=0,lasti=0,lastj=0,laste=0;
		for(vector1<Size>::const_iterator biter = scanres.begin(); biter != scanres.end(); ++biter) {
			Size brsd = *biter;

			TR << "scanning bpy rsd " << infile << " " << brsd << std::endl;
			if(brsd != lastb) {	lastb = brsd; foundcen.clear(); foundori1.clear(); foundori2.clear(); foundori3.clear(); foundori4.clear(); foundori5.clear(); foundori6.clear();	}
			Stub s3(init_pose.xyz(AtomID(5,brsd)),init_pose.xyz(AtomID(2,brsd)),init_pose.xyz(AtomID(1,brsd)));
			for(Size kch1 = 1; kch1 <= CHI1.size(); ++kch1) {
				for(Size kch2 = 1; kch2 <= CHI2.size(); ++kch2) {
					if(!clashcheck.clash_check(s3.local2global(bci2end(kch1,kch2)))) { continue; }
					if(!clashcheck.clash_check(s3.local2global(bci2cen(kch1,kch2)))) { continue; }
					pose.replace_residue(brsd,bpy.residue(1),true);
					pose.set_chi(1,brsd,CHI1[kch1]);
					pose.set_chi(2,brsd,CHI2[kch2]);

					bool clash = false;
					for(Size k = 6; k <= pose.residue(brsd).nheavyatoms(); ++k) {
						if(!clashcheck.clash_check(pose.residue(brsd).xyz(k),brsd)) clash=true;
						if(clash) break;
					}
					if(clash) { continue; }

					Vec const cenb = pose.residue(brsd).xyz("ZN");
					Vec orik1 = (pose.residue(brsd).xyz("NE1")-cenb).normalized();
					Vec orik2 = (pose.residue(brsd).xyz("NN1")-cenb).normalized();
					Vec const ligk1 = pose.residue(brsd).xyz("NE1");
					Vec const ligk2 = pose.residue(brsd).xyz("NN1");
					{
						Real ang = (90.0-angle_degrees(orik1,Vec(0,0,0),orik2))/2.0;
						orik1 = rotation_matrix_degrees(orik2.cross(orik1),ang)*orik1;
						orik2 = rotation_matrix_degrees(orik1.cross(orik2),ang)*orik2;
					}
					for(vector1<Size>::const_iterator iiter = scanres.begin(); iiter != scanres.end(); ++iiter) {
						Size irsd = *iiter;
						if(irsd==brsd) { continue; }
						// if(irsd != lasti) {	lasti = irsd; foundcen.clear(); foundori1.clear(); foundori2.clear(); foundori3.clear(); foundori4.clear(); foundori5.clear(); foundori6.clear();	}
						if(pose.residue(irsd).xyz("CB").distance_squared(cenb) > 36.0) { continue; }
						Stub si(init_pose.xyz(AtomID(5,irsd)),init_pose.xyz(AtomID(2,irsd)),init_pose.xyz(AtomID(1,irsd)));
						for(Size ich1 = 1; ich1 <= CHI1.size(); ++ich1) {
							for(Size ich2 = 1; ich2 <= CHI2.size(); ++ich2) {
								for(Size ide = 1; ide <= 2; ide++) {
									Vec const ceni = si.local2global(chi2cen(ide,ich1,ich2));
									if(ceni.distance_squared(cenb) > MXDSMTL) { continue; }
									Vec const orii = (si.local2global(chi2ori(ide,ich1,ich2))-ceni).normalized();
									bool his180 = false;
									{
										Real ang1 = numeric::conversions::degrees(acos(orii.dot(orik1)));
										Real ang2 = numeric::conversions::degrees(acos(orii.dot(orik2)));
										if( ang1 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL ) { continue; }
										if( ang2 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang2 && ang2 < 180-MXAGMTL ) { continue; }
										his180 = (ang1 > 135.0 || ang2 > 135.0);
									}
									if(!clashcheck.clash_check(si.local2global(chi2cen(ide,ich1,ich2)))) { continue; }

									pose.replace_residue(irsd,( ide==1 ? hsd : hse ).residue(1),true);
									pose.set_chi(1,irsd,CHI1[ich1]);
									pose.set_chi(2,irsd,CHI2[ich2]);
									Vec const ligi = pose.residue(irsd).xyz( (ide==1) ? "ND1" : "NE2" );

									clash = false;
									for(vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b)
										for(vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h)
											if( pose.residue(brsd).xyz(*b).distance_squared(pose.residue(irsd).xyz(*h)) < 7.0 ) clash=true;
									if(clash) { continue; }

									for(Size a = 6; a <= pose.residue(irsd).nheavyatoms(); ++a) {
										if(!clashcheck.clash_check(pose.residue(irsd).xyz(a),irsd)) { clash=true; if(clash) break; }
									}
									if(clash) { continue; }

									for(vector1<Size>::const_iterator jiter = scanres.begin(); jiter != scanres.end(); ++jiter) {
										Size jrsd = *jiter;
										if(jrsd <= irsd) { continue; }
										if(jrsd==brsd || jrsd==irsd) { continue; }
										// if(jrsd != lastj) {	lastj = jrsd; foundcen.clear(); foundori1.clear(); foundori2.clear(); foundori3.clear(); foundori4.clear(); foundori5.clear(); foundori6.clear();	}
										if(pose.residue(jrsd).xyz("CB").distance_squared(cenb) > 36.0) { continue; }
										Stub sj(init_pose.xyz(AtomID(5,jrsd)),init_pose.xyz(AtomID(2,jrsd)),init_pose.xyz(AtomID(1,jrsd)));
										for(Size jch1 = 1; jch1 <= CHI1.size(); ++jch1) {
											for(Size jch2 = 1; jch2 <= CHI2.size(); ++jch2) {
												for(Size jde = 1; jde <= 2; jde++) {
													Vec const cenj = sj.local2global(chi2cen(jde,jch1,jch2));
													if(cenj.distance_squared(cenb) > MXDSMTL) { continue; }
													clash = false;
													Vec const orij = (sj.local2global(chi2ori(jde,jch1,jch2))-cenj).normalized();
													bool localhis180 = his180;
													bool hishis180 = false;
													{
														Real ang1 = numeric::conversions::degrees(acos(orij.dot(orii)));
														if( ang1 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL ) { continue; }
														if( ang1 > 135.0 ) hishis180 = true;
													}
													if( !localhis180 && !hishis180 ) {
														// TR << "NOT his180, hishis180" << std::endl;
														Real ang1 = numeric::conversions::degrees(acos(orij.dot(orik1)));
														Real ang2 = numeric::conversions::degrees(acos(orij.dot(orik2)));
														if( ang1 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL ) { continue; }
														if( ang2 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang2 && ang2 < 180-MXAGMTL ) { continue; }
														localhis180 = (ang1 > 135.0 || ang2 > 135.0);
													} else {
														// TR << "his180 OR hishis180" << std::endl;
														Real ang1 = numeric::conversions::degrees(acos(orij.dot(orik1)));
														Real ang2 = numeric::conversions::degrees(acos(orij.dot(orik2)));
														if( ang1 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang1) { continue; }
														if( ang2 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang2) { continue; }
													}
													if(!clashcheck.clash_check(sj.local2global(chi2cen(jde,jch1,jch2)))) { continue; }

													pose.replace_residue(jrsd,( jde==1 ? hsd : hse ).residue(1),true);
													pose.set_chi(1,jrsd,CHI1[jch1]);
													pose.set_chi(2,jrsd,CHI2[jch2]);
													Vec const ligj = pose.residue(jrsd).xyz( (jde==1) ? "ND1" : "NE2" );

													for(vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h) {
														for(vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b) {
															if( pose.residue(brsd).xyz(*b).distance_squared(pose.residue(jrsd).xyz(*h)) < 7.0 ) {
																clash=true;
															}
														}
														for(vector1<Size>::const_iterator h2 = hatm.begin(); h2 != hatm.end(); ++h2) {
															if( pose.residue(irsd).xyz(*h2).distance_squared(pose.residue(jrsd).xyz(*h)) < 7.0 ) {
																clash=true;
															}
														}
													}
													if(clash) { continue; }

													for(Size a = 6; a <= pose.residue(jrsd).nheavyatoms(); ++a) {
														if(!clashcheck.clash_check(pose.residue(jrsd).xyz(a),jrsd)) { clash=true; if(clash) break; }
													}
													if(clash) { continue; }

													// {
													// 	bool output = true;
													// 	if(ceni.distance_squared(cenb) > MXDSMTL/2.0) output = false;
													// 	if(cenj.distance_squared(cenb) > MXDSMTL/2.0) output = false;
													// 	{
													// 		Real ang1 = numeric::conversions::degrees(acos(orii.dot(orik1)));
													// 		Real ang2 = numeric::conversions::degrees(acos(orii.dot(orik2)));
													// 		if( ang1 < 90.0-MXAGMTL/2.0 || 90.0+MXAGMTL/2.0 < ang1 && ang1 < 180-MXAGMTL/2.0 ) output = false;
													// 		if( ang2 < 90.0-MXAGMTL/2.0 || 90.0+MXAGMTL/2.0 < ang2 && ang2 < 180-MXAGMTL/2.0 ) output = false;
													// 	}
													// 	{
													// 		Real ang1 = numeric::conversions::degrees(acos(orij.dot(orii)));
													// 		if( ang1 < 90.0-MXAGMTL/2.0 || 90.0+MXAGMTL/2.0 < ang1 && ang1 < 180-MXAGMTL/2.0 ) output = false;
													// 	}
													// 	{
													// 		Real ang1 = numeric::conversions::degrees(acos(orij.dot(orik1)));
													// 		Real ang2 = numeric::conversions::degrees(acos(orij.dot(orik2)));
													// 		if( ang1 < 90.0-MXAGMTL/2.0 || 90.0+MXAGMTL/2.0 < ang1 && ang1 < 180-MXAGMTL/2.0 ) output = false;
													// 		if( ang2 < 90.0-MXAGMTL/2.0 || 90.0+MXAGMTL/2.0 < ang2 && ang2 < 180-MXAGMTL/2.0 ) output = false;
													// 	}
													// 	if(output) {
													// 		pose = fa_pose;
													// 		pose.replace_residue(brsd,bpy.residue(1),true);
													// 		pose.set_chi(1,brsd,CHI1[kch1]);
													// 		pose.set_chi(2,brsd,CHI2[kch2]);
													// 		pose.replace_residue(irsd,( ide==1 ? hsd : hse ).residue(1),true);
													// 		pose.set_chi(1,irsd,CHI1[ich1]);
													// 		pose.set_chi(2,irsd,CHI2[ich2]);
													// 		pose.replace_residue(jrsd,( jde==1 ? hsd : hse ).residue(1),true);
													// 		pose.set_chi(1,jrsd,CHI1[jch1]);
													// 		pose.set_chi(2,jrsd,CHI2[jch2]);
													// 		string outfname = utility::file_basename(infile)+"_"+string_of(brsd)+"_"+string_of(irsd)+"_"+string_of(jrsd)+"_"+string_of(numeric::random::uniform())+".pdb.gz";
													// 		TR << "HIT! " << outfname << std::endl;
													// 		pose.dump_pdb(option[out::file::o]()+"/"+outfname);
													// 	}
													// }


													for(vector1<Size>::const_iterator citer = scanres.begin(); citer != scanres.end(); ++citer) {
														Size crsd = *citer;
														if(crsd==brsd||crsd==irsd||crsd==jrsd) { continue; }
														if(pose.residue(crsd).xyz("CB").distance_squared(cenb) > 16.0) { continue; }
														clash = false;
														Vec cbca = (pose.residue(crsd).xyz(5)-pose.residue(crsd).xyz(2)).normalized();
														Vec CB = pose.residue(crsd).xyz(5);
														Vec SG = CB - cbca*1.8;
														SG = rotation_matrix_degrees(cbca.cross(pose.residue(crsd).xyz(2)),115.8)*(SG-CB)+CB;
														Mat Rc = rotation_matrix_degrees(cbca,chi1incr);
														for(Size cch1 = 1; cch1 <= CHI1.size(); ++cch1) {
															SG = Rc*(SG-CB)+CB;
															if(fabs(SG.distance(cenb)-2.2) > MXDSMTL) { continue; }
															if(fabs(angle_degrees(CB,SG,cenb)-100) > MXAGMTL ) { continue; }

															for(vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b)
																if( pose.residue(brsd).xyz(*b).distance_squared(SG) < 8.0 ) clash=true;
															for(vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h) {
																if( pose.residue(irsd).xyz(*h).distance_squared(SG) < 8.0 ) clash=true;
																if( pose.residue(jrsd).xyz(*h).distance_squared(SG) < 8.0 ) clash=true;
															}
															if(clash) { continue; }

															if(!clashcheck.clash_check(SG,crsd)) { clash=true; if(clash) break; }
															if(clash) { continue; }

															Vec oric = (cenb-SG).normalized();
															{
																Real ang = angle_degrees( SG, cenb, ligi );
																if( ang < 90.0-MXAGMTL || 90.0+MXAGMTL < ang && ang < 180-MXAGMTL ) { continue; }
																ang = angle_degrees( SG, cenb, ligj );
																if( ang < 90.0-MXAGMTL || 90.0+MXAGMTL < ang && ang < 180-MXAGMTL ) { continue; }
																ang = angle_degrees( SG, cenb, ligk1 );
																if( ang < 90.0-MXAGMTL || 90.0+MXAGMTL < ang && ang < 180-MXAGMTL ) { continue; }
																ang = angle_degrees( SG, cenb, ligk2 );
																if( ang < 90.0-MXAGMTL || 90.0+MXAGMTL < ang && ang < 180-MXAGMTL ) { continue; }
															}

															TR << "BHHC HIT!!! looking for GLU" << std::endl;

															for(vector1<Size>::const_iterator eiter = scanres.begin(); eiter != scanres.end(); ++eiter) {
																Size ersd = *eiter;
																if(ersd==brsd||ersd==irsd||ersd==jrsd||ersd==crsd) { continue; }
																if(pose.residue(ersd).xyz("CB").distance_squared(cenb) > 49.0) { continue; }
																// if(ersd != laste) {	laste = ersd; foundcen.clear(); foundori1.clear(); foundori2.clear(); foundori3.clear(); foundori4.clear(); foundori5.clear(); foundori6.clear();	}
																clash = false;
																Stub se(init_pose.xyz(AtomID(5,ersd)),init_pose.xyz(AtomID(2,ersd)),init_pose.xyz(AtomID(1,ersd)));
																for(Size ech1 = 1; ech1 <= CHI1.size(); ++ech1) {
																	for(Size ech2 = 1; ech2 <= CHI2.size(); ++ech2) {
																		//Vec const cene = se.local2global(gci2cen(ech1,ech2));
																		//if(cene.distance_squared(cenb) > MXDSMTL) { continue; }
																		Vec const CD = se.local2global(gci2cd(ech1,ech2));
																		Vec const orie = (cenb-CD).normalized();
																		if( CD.distance_squared(cenb) > 16.0 ) continue;

																		pose.replace_residue(ersd,glu.residue(1),true);
																		pose.set_chi(1,ersd,CHI1[ech1]);
																		pose.set_chi(2,ersd,CHI2[ech2]);

																		for(vector1<Size>::const_iterator e = eatm.begin(); e != eatm.end(); ++e) {
																			for(vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b)
																				if( pose.residue(brsd).xyz(*b).distance_squared(pose.residue(ersd).xyz(*e)) < 7.0 ) clash=true;
																			for(vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h)
																				if( pose.residue(irsd).xyz(*h).distance_squared(pose.residue(ersd).xyz(*e)) < 7.0 ) clash=true;
																			for(vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h)
																				if( pose.residue(jrsd).xyz(*h).distance_squared(pose.residue(ersd).xyz(*e)) < 7.0 ) clash=true;
																			if( SG.distance_squared(pose.residue(ersd).xyz(*e)) < 7.0 ) clash=true;
																		}
																		if(clash) { continue; }

																		for(Size a = 6; a <= pose.residue(ersd).nheavyatoms()-2; ++a) {
																			if(!clashcheck.clash_check(pose.residue(ersd).xyz(a),ersd)) { clash=true; if(clash) break; }
																		}
																		if(clash) { continue; }


																		Vec cen = (2*cenb + ceni + cenj)/4.0;
																		bool overlap = false;
																		Vec ori1 = orik1;
																		Vec ori2 = orik2;
																		Vec ori3 = orii;
																		Vec ori4 = orij;
																		Vec ori5 = oric;
																		Vec ori6 = orie;
																		for(Size i = 1; i <= foundcen.size(); ++i) {
																			if( cen.distance_squared(foundcen[i]) < option[willmatch::match_overlap_dis]() &&
																				ori1.dot(foundori1[i])            > MATCH_OVERLAP_DOT                      &&
																				ori2.dot(foundori2[i])            > MATCH_OVERLAP_DOT                      &&
																				ori3.dot(foundori3[i])            > MATCH_OVERLAP_DOT                      &&
																				ori4.dot(foundori4[i])            > MATCH_OVERLAP_DOT                      &&
																				ori5.dot(foundori5[i])            > MATCH_OVERLAP_DOT                      &&
																				ori6.dot(foundori6[i])            > MATCH_OVERLAP_DOT                      ){
																				 // TR << "overlap!" << std::endl;
																				overlap = true;
																				break;
																			}
																		}
																		if(overlap) continue;
																		count++;
																		foundcen.push_back(cen);
																		foundori1.push_back(ori1);
																		foundori2.push_back(ori2);
																		foundori3.push_back(ori3);
																		foundori4.push_back(ori4);
																		foundori5.push_back(ori5);
																		foundori6.push_back(ori6);
																		{
																			Pose opose = native;
																			opose.replace_residue(brsd,bpy.residue(1),true);
																			opose.set_chi(1,brsd,CHI1[kch1]);
																			opose.set_chi(2,brsd,CHI2[kch2]);
																			opose.replace_residue(irsd,( ide==1 ? hsd : hse ).residue(1),true);
																			opose.set_chi(1,irsd,CHI1[ich1]);
																			opose.set_chi(2,irsd,CHI2[ich2]);
																			opose.replace_residue(jrsd,( jde==1 ? hsd : hse ).residue(1),true);
																			opose.set_chi(1,jrsd,CHI1[jch1]);
																			opose.set_chi(2,jrsd,CHI2[jch2]);
																			opose.replace_residue(crsd,cys.residue(1),true);
																			opose.set_xyz(AtomID(6,crsd),SG);
																			opose.replace_residue(ersd,glu.residue(1),true);
																			opose.set_chi(1,ersd,CHI1[ech1]);
																			opose.set_chi(2,ersd,CHI2[ech2]);
																			string outfname = utility::file_basename(infile)+"_FULL_"+string_of(brsd)+"_"+string_of(irsd)+"_"+string_of(jrsd)+"_"+string_of(crsd)+"_"+string_of(kch1)+"_"+string_of(kch2)+"_"+string_of(ich1)+"_"+string_of(ich2)+"_"+string_of(ide)+"_"+string_of(jch1)+"_"+string_of(jch2)+"_"+string_of(jde)+"_"+string_of(cch1)+".pdb.gz";
																			TR << "HIT! " << outfname << std::endl;
																			opose.dump_pdb(option[out::file::o]()+"/"+outfname);
																			utility_exit_with_message("testing");
																		}
																	}
																}
															}

															//}
														}
													}

												}
											}
										}
									}

								}
							}
						}
					}

				}
			}
		}

		// out.close();
	}
}


int main (int argc, char *argv[]) {

	try {

	devel::init(argc,argv);
	run();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


