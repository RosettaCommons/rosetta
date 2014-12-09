// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <apps/pilot/will/will_util.ihh>

using core::Real;
using core::Size;
using core::pose::Pose;
using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using std::string;
using utility::vector1;
using ObjexxFCL::string_of;
using ObjexxFCL::lead_zero_string_of;
using numeric::min;
using core::import_pose::pose_from_pdb;
using basic::options::option;
using numeric::min;
using numeric::max;

typedef utility::vector1<core::Real> Reals;
typedef utility::vector1<core::Size> Sizes;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

static thread_local basic::Tracer TR( "genmatch_d6_bpy" );



void myoptH(Pose & pose, ScoreFunctionOP sf) {
	add_lower_terminus_type_to_pose_residue(pose,1);
	add_upper_terminus_type_to_pose_residue(pose,pose.n_residue());
	core::pack::optimizeH(pose,*sf);
	remove_lower_terminus_type_from_pose_residue(pose,1);
	remove_upper_terminus_type_from_pose_residue(pose,pose.n_residue());
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

	core::io::silent::SilentFileData sfd;

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
		Pose in_fa;
		pose_from_pdb(in_fa, *fa_residue_set,infile);
		for(Size ir = 1; ir <= in_fa.size(); ++ir) {
			if(in_fa.residue(ir).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(in_fa,ir);
			if(in_fa.residue(ir).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(in_fa,ir);
		}
		Pose native = in_fa;
		Size nres = in_fa.n_residue();
		core::chemical::ResidueType const & rtala( in_fa.residue(1).residue_type_set().name_map("ALA") );
		core::chemical::ResidueType const & rthis( in_fa.residue(1).residue_type_set().name_map("HIS") );
		core::chemical::ResidueType const & rtglu( in_fa.residue(1).residue_type_set().name_map("GLU") );
		core::chemical::ResidueType const & rtbpy( in_fa.residue(1).residue_type_set().name_map("BPY") );
		core::chemical::ResidueType const & rtphe( in_fa.residue(1).residue_type_set().name_map("PHE") );
		for(Size i = 1; i <= nres; ++i) core::pose::replace_pose_residue_copying_existing_coordinates(in_fa,i,rtala);
		Pose const fa_pose = in_fa;
		ImplicitFastClashCheck clashcheck(fa_pose,basic::options::option[basic::options::OptionKeys::willmatch::clash_dis]());

		ScoreFunctionOP sf = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
		sf->set_weight(core::scoring::fa_dun,1.0);

		Real chi1incr = option[willmatch::chi1_increment]();
		Real chi2incr = option[willmatch::chi2_increment]();
		vector1<Real> CHI1,CHI2;
		for(Real i = 0; i < 360; i+= chi1incr) CHI1.push_back( i>180.0 ? i-360.0 : i );
		for(Real i = 0; i < 360; i+= chi2incr) CHI2.push_back( i>180.0 ? i-360.0 : i );

		utility::io::ozstream bhhout(option[out::file::o]()+"/"+utility::file_basename(infile)+".bhh_match");
		// setup HIS residues for checking
		Pose hse,hsd,bpy,glu;
		core::pose::make_pose_from_sequence(hse,"H[HIS]"  ,*fa_residue_set,false);
		core::pose::make_pose_from_sequence(hsd,"H[HIS_D]",*fa_residue_set,false);
		core::pose::make_pose_from_sequence(bpy,"X[BPY]"  ,*fa_residue_set,false);
		core::pose::make_pose_from_sequence(glu,"E"       ,*fa_residue_set,false);
		core::pose::remove_lower_terminus_type_from_pose_residue(hse,1);
		core::pose::remove_lower_terminus_type_from_pose_residue(hsd,1);
		core::pose::remove_lower_terminus_type_from_pose_residue(bpy,1);
		core::pose::remove_lower_terminus_type_from_pose_residue(glu,1);
		core::pose::remove_upper_terminus_type_from_pose_residue(hse,1);
		core::pose::remove_upper_terminus_type_from_pose_residue(hsd,1);
		core::pose::remove_upper_terminus_type_from_pose_residue(bpy,1);
		core::pose::remove_upper_terminus_type_from_pose_residue(glu,1);
		hsd.set_dof(DOF_ID(AtomID(hsd.residue(1).atom_index("HD1"),1),D),2.1);
		hse.set_dof(DOF_ID(AtomID(hse.residue(1).atom_index("HE2"),1),D),2.1);
		// hse.dump_pdb("hse.pdb");
		// hsd.dump_pdb("hsd.pdb");
		ObjexxFCL::FArray3D<Vec> chi2cen(2,CHI1.size(),CHI2.size()),chi2ori(2,CHI1.size(),CHI2.size());
		ObjexxFCL::FArray2D<Vec> bci2cen(  CHI1.size(),CHI2.size()),bci2ori(  CHI1.size(),CHI2.size()),bci2end(  CHI1.size(),CHI2.size());
		ObjexxFCL::FArray2D<Vec> gci2cen(  CHI1.size(),CHI2.size()),gci2ori(  CHI1.size(),CHI2.size());
		ObjexxFCL::FArray2D<Vec> hci1cen(2,CHI1.size());
		ObjexxFCL::FArray1D<Vec> bci1cen(  CHI1.size());
		ObjexxFCL::FArray1D<Vec> eci1cen(  CHI1.size());
		Stub stub(hse.residue(1).xyz("CB"),hse.residue(1).xyz("CA"),hse.residue(1).xyz("CA")+(hse.residue(1).xyz("N")-hse.residue(1).xyz("C")));
		// precompute cen and ori for hisd and hise
	//Real bdscbmx2 = (MXDSMTL+5.64+0.10)*(MXDSMTL+5.64+0.10);
		Real hddcbmx2 = (MXDSMTL+3.54+0.10)*(MXDSMTL+3.54+0.10);
		Real hedcbmx2 = (MXDSMTL+5.74+0.10)*(MXDSMTL+5.74+0.10);
		Real edscbmx2 = (MXDSMTL+5.24+0.10)*(MXDSMTL+5.24+0.10);
	//Real bdsc1mx2 = (MXDSMTL+3.24+0.10)*(MXDSMTL+3.24+0.10);
		Real hddc1mx2 = (MXDSMTL+3.30+0.10)*(MXDSMTL+3.30+0.10);
		Real hedc1mx2 = (MXDSMTL+1.64+0.10)*(MXDSMTL+1.64+0.10);
		Real edsc1mx2 = (MXDSMTL+3.87+0.10)*(MXDSMTL+3.87+0.10);
		for(Size i = 1; i <= CHI1.size(); ++i) {
			hsd.set_chi(1,1,CHI1[i]);
			hse.set_chi(1,1,CHI1[i]);
			bpy.set_chi(1,1,CHI1[i]);
			glu.set_chi(1,1,CHI1[i]);
			hci1cen(1,i) = stub.global2local( hsd.residue(1).xyz("CG")-0.007*((hsd.residue(1).xyz("CG")-hsd.residue(1).xyz("CB")).normalized()) ); // CG - 0.006699*(CG-CB) 3.11925602854 from actual cen
			hci1cen(2,i) = stub.global2local( hse.residue(1).xyz("CG")+3.947*((hse.residue(1).xyz("CG")-hse.residue(1).xyz("CB")).normalized()) ); //	CG + 3.946941*(CG-CB) 1.53906602847 from actual cen
			bci1cen(  i) = stub.global2local( bpy.residue(1).xyz("CG")+3.034*((bpy.residue(1).xyz("CG")-bpy.residue(1).xyz("CB")).normalized()) ); // GC + 3.033686*(CG-CB) 3.23823647147 from actual cen
			eci1cen(  i) = stub.global2local( glu.residue(1).xyz("CG")+1.991*((glu.residue(1).xyz("CG")-glu.residue(1).xyz("CB")).normalized()) ); // CG + 1.991210*(CG-CB) 3.86805118484 from actual cen
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
			}
		}

		vector1<Size> scanres;
		if(option[willmatch::residues].user()) {
			TR << "input scanres!!!!!!" << std::endl;
			scanres = option[willmatch::residues]();
		} else {
			for(Size i = 1; i <= in_fa.n_residue(); ++i) {
				if(!in_fa.residue(i).has("N" )) { continue; }
				if(!in_fa.residue(i).has("CA")) { continue; }
				if(!in_fa.residue(i).has("C" )) { continue; }
				if(!in_fa.residue(i).has("O" )) { continue; }
				if(!in_fa.residue(i).has("CB")) { continue; }
				if(in_fa.residue(i).name3()=="PRO") { continue; }
				scanres.push_back(i);
			}
		}

		// precompute acceptable chis
		ObjexxFCL::FArray1D<Real> ballow0(                        in_fa.size(),9e9),eallow0(                        in_fa.size(),9e9),hallow0(                        in_fa.size(),9e9);
		ObjexxFCL::FArray2D<Real> ballow1(            CHI1.size(),in_fa.size(),9e9),eallow1(            CHI1.size(),in_fa.size(),9e9),hallow1(            CHI1.size(),in_fa.size(),9e9);
		ObjexxFCL::FArray3D<Real> ballow2(CHI2.size(),CHI1.size(),in_fa.size(),9e9),eallow2(CHI2.size(),CHI1.size(),in_fa.size(),9e9),hallow2(CHI2.size(),CHI1.size(),in_fa.size(),9e9);
		{
			ScoreFunctionOP sftmp = new core::scoring::ScoreFunction;
			sftmp->set_weight(core::scoring::fa_dun,1.0);
			TR << "precomputing dun scores and clashes" << std::endl;
			// BPY
			Pose pose = in_fa,pose3=bpy;
			pose3.append_residue_by_bond(pose.residue(2));
			pose3.append_residue_by_bond(pose.residue(3));
			for(Sizes::const_iterator iri = scanres.begin(); iri != scanres.end(); ++iri) {
				Size const ir(*iri);
				core::pose::replace_pose_residue_copying_existing_coordinates(pose,ir,rtphe); // clash w/BPY

				pose3.replace_residue(1,pose.residue(((ir-1+pose.size()-1)%(pose.size()))+1),false);
				pose3.replace_residue(2,pose.residue(((ir+0+pose.size()-1)%(pose.size()))+1),false);
				pose3.replace_residue(3,pose.residue(((ir+1+pose.size()-1)%(pose.size()))+1),false);
				for(Size i1 = 1; i1 <= CHI1.size(); ++i1) {
					pose3.set_chi(1,2,CHI1[i1]);
					for(Size i2 = 1; i2 <= CHI2.size(); ++i2) {
						pose3.set_chi(2,2,CHI2[i2]);
						sftmp->score(pose3);
						ballow2(i2,i1,ir) =	pose3.energies().residue_total_energies(2)[core::scoring::fa_dun]; // dun score w/PHE
					}
				}
				pose.replace_residue(ir,bpy.residue(1),true);
				pose3.replace_residue(2,pose.residue(ir),false);
				for(Size i1 = 1; i1 <= CHI1.size(); ++i1) {
					pose3.set_chi(1,2,CHI1[i1]);
					for(Size i2 = 1; i2 <= CHI2.size(); ++i2) {
						pose3.set_chi(2,2,CHI2[i2]);
						for(Size a = 6; a <= pose3.residue(2).nheavyatoms(); ++a) {
							if(!clashcheck.clash_check(pose3.residue(2).xyz(a),ir)) {
								ballow2(i2,i1,ir) += 9e9;
								break;
							}
						}
						if( ballow1(i1,ir) > ballow2(i2,i1,ir) ) ballow1(i1,ir) = ballow2(i2,i1,ir);

						// pose.set_chi(1,ir,CHI1[i1]);
						// pose.set_chi(2,ir,CHI2[i2]);
						// pose.dump_pdb("test1.pdb");
						// pose3.dump_pdb("test3.pdb");
						// utility_exit_with_message("");

					}
					if( ballow0(ir) > ballow1(i1,ir) ) ballow0(ir) = ballow1(i1,ir);
				}
				// GLU
				pose.replace_residue(ir,glu.residue(1),true);
				pose3.replace_residue(2,pose.residue(ir),false);
				pose3.set_phi(2,pose.phi(ir)); pose3.set_psi(2,pose.psi(ir)); pose3.set_omega(2,pose.omega(ir));
				pose3.set_chi(3,2,0.0);
				for(Size i1 = 1; i1 <= CHI1.size(); ++i1) {
					pose3.set_chi(1,2,CHI1[i1]);
					for(Size i2 = 1; i2 <= CHI2.size(); ++i2) {
						pose3.set_chi(2,2,CHI2[i2]);
						sftmp->score(pose3);
						eallow2(i2,i1,ir) =	pose3.energies().residue_total_energies(2)[core::scoring::fa_dun];
						for(Size a = 6; a <= pose3.residue(2).nheavyatoms()-2; ++a) {
							if(!clashcheck.clash_check(pose3.residue(2).xyz(a),ir)) {
								eallow2(i2,i1,ir) += 9e9;
								break;
							}
						}
						if( eallow1(i1,ir) > eallow2(i2,i1,ir) ) eallow1(i1,ir) = eallow2(i2,i1,ir);
					}
					if( eallow0(ir) > eallow1(i1,ir) ) eallow0(ir) = eallow1(i1,ir);
				}
				// HIS
				pose.replace_residue(ir,hsd.residue(1),true);
				pose3.replace_residue(2,pose.residue(ir),false);
				pose.replace_residue(ir,hse.residue(1),true);
				pose3.set_phi  (2,pose.phi  (ir)); pose3.set_psi  (2,pose.psi  (ir));	pose3.set_omega(2,pose.omega(ir));
				for(Size i1 = 1; i1 <= CHI1.size(); ++i1) {
					pose3.set_chi(1,2,CHI1[i1]);
					for(Size i2 = 1; i2 <= CHI2.size(); ++i2) {
						pose3.set_chi(2,2,CHI2[i2]);
						pose3.energies().clear_energies();
						sftmp->score(pose3);
						hallow2(i2,i1,ir) =	pose3.energies().residue_total_energies(2)[core::scoring::fa_dun];

						if(21==ir&&38==i1&&25==i2) {
							TR << "PRE " << pose3.phi(2) << " " << pose3.psi(2) << " " << pose3.omega(2) << " " << pose3.chi(1,2) << " " << pose3.chi(2,2) << std::endl;
							pose3.dump_pdb("test3.pdb");
							TR << pose3.chi(1,2) << " " << pose3.chi(2,2) << " " << pose3.residue(2).name() << " " << pose3.energies().residue_total_energies(2)[core::scoring::fa_dun] << std::endl;
							pose3.set_chi(1,2,pose3.chi(1,2)+360.0);
							pose3.set_chi(2,2,pose3.chi(2,2)+360.0);
							pose3.energies().clear_energies();
							sftmp->score(pose3);
							TR << pose3.chi(1,2) << " " << pose3.chi(2,2) << " " << pose3.residue(2).name() << " " << pose3.energies().residue_total_energies(2)[core::scoring::fa_dun] << std::endl;
							pose3.set_chi(1,2,pose3.chi(1,2)-360.0);
							pose3.set_chi(2,2,pose3.chi(2,2)-360.0);
							pose3.energies().clear_energies();
							sftmp->score(pose3);
							TR << pose3.chi(1,2) << " " << pose3.chi(2,2) << " " << pose3.residue(2).name() << " " << pose3.energies().residue_total_energies(2)[core::scoring::fa_dun] << std::endl;
							pose3.set_chi(1,2,pose3.chi(1,2)-360.0);
							pose3.set_chi(2,2,pose3.chi(2,2)-360.0);
							pose3.energies().clear_energies();
							sftmp->score(pose3);
							TR << pose3.chi(1,2) << " " << pose3.chi(2,2) << " " << pose3.residue(2).name() << " " << pose3.energies().residue_total_energies(2)[core::scoring::fa_dun] << std::endl;
						}
						// pose.set_chi(1,ir,CHI1[i1]);
						// pose.set_chi(2,ir,CHI2[i2]);
						// //pose.set_chi(3,ir,0.0);
						// sftmp->score(pose);
						// Real db2 = pose.energies().residue_total_energies(ir)[core::scoring::fa_dun];
						// if( fabs(hallow2(i2,i1,ir)-db2) > 0.001 ) {
						// 	TR << "DUN CMP " << ir << " " << hallow2(i2,i1,ir) << " " << db2 << std::endl;
						// 	TR << "PHI " << pose.phi(  ir) << " " << pose3.phi(  2) << " " << std::endl;
						// 	TR << "PSI " << pose.psi(  ir) << " " << pose3.psi(  2) << " " << std::endl;
						// 	TR << "OMG " << pose.omega(ir) << " " << pose3.omega(2) << " " << std::endl;
						// 	TR << "CHI " << pose.chi(1,ir) << " " << pose3.chi(1,2) << " " << std::endl;
						// 	TR << "CHI " << pose.chi(2,ir) << " " << pose3.chi(2,2) << " " << std::endl;
						// 	pose3.dump_pdb("pose3.pdb");
						// 	pose.dump_pdb("pose.pdb");
						// 	utility_exit_with_message("lkjasldfj");
						// }

						for(Size a = 6; a <= pose3.residue(2).nheavyatoms(); ++a) {
							if(!clashcheck.clash_check(pose3.residue(2).xyz(a),ir)) {
								hallow2(i2,i1,ir) += 9e9;
								break;
							}
						}
						if( hallow1(i1,ir) > hallow2(i2,i1,ir) ) hallow1(i1,ir) = hallow2(i2,i1,ir);
					}
					if( hallow0(ir) > hallow1(i1,ir) ) hallow0(ir) = hallow1(i1,ir);
				}
				core::pose::replace_pose_residue_copying_existing_coordinates(pose,ir,rtala);
			}
			TR << "DONE precomputing dun scores and clashes" << std::endl;
		}
		Real const DUN_THRESH = option[willmatch::fa_dun_thresh]();


		// string fname = utility::file_basename(infile)+"_d6_bpy_matches.pdb";
		// utility::io::ozstream out(fname);
		// out << "MODEL BASE"	<< endl;
		// fa_pose.dump_pdb(out);
		// out << "ENDMDL" << endl;

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

		vector1<Vec> Ns(in_fa.n_residue()),CAs(in_fa.n_residue()),Cs(in_fa.n_residue()),CBs(in_fa.n_residue());
		for(Size i = 1; i <= in_fa.n_residue(); ++i) Ns [i] = in_fa.residue(i).xyz("N" );
		for(Size i = 1; i <= in_fa.n_residue(); ++i) CAs[i] = in_fa.residue(i).xyz("CA");
		for(Size i = 1; i <= in_fa.n_residue(); ++i) Cs [i] = in_fa.residue(i).xyz("C" );
		for(Size i = 1; i <= in_fa.n_residue(); ++i) CBs[i] = in_fa.residue(i).xyz("CB");

		Pose pose = fa_pose;
		//		Size lastb=0,lasti=0,lastj=0,laste=0;
		for(vector1<Size>::const_iterator biter = scanres.begin(); biter != scanres.end(); ++biter) {
			vector1<Vec> foundcenbg,foundori1bg,foundori2bg,foundori3bg,foundori4bg,foundori5bg,foundori6bg;
			vector1<Vec> foundcensg,foundori1sg,foundori2sg,foundori3sg,foundori4sg,foundori5sg,foundori6sg;
			Size brsd = *biter;
			//			if(brsd != 110) continue;
			TR << "scanning bpy rsd " << infile << " " << brsd << std::endl;
			Stub s3(CBs[brsd],CAs[brsd],CAs[brsd]+(Ns[brsd]-Cs[brsd]));
			if(ballow0(brsd) > DUN_THRESH) continue; // if all chi1 / chi2 are clash or bad rot
			for(Size kch1 = 1; kch1 <= CHI1.size(); ++kch1) {
				if(ballow1(kch1,brsd) > DUN_THRESH) continue; // if all chi2 are clash or bad rot
				for(Size kch2 = 1; kch2 <= CHI2.size(); ++kch2) {
					if(ballow2(kch2,kch1,brsd) > DUN_THRESH) continue;  // chi2 is clash or bad rot
					if(!clashcheck.clash_check(s3.local2global(bci2end(kch1,kch2)))) { continue; }
					if(!clashcheck.clash_check(s3.local2global(bci2cen(kch1,kch2)))) { continue; }
					pose.replace_residue(brsd,bpy.residue(1),true);
					pose.set_chi(1,brsd,CHI1[kch1]);
					pose.set_chi(2,brsd,CHI2[kch2]);

					bool clash = false;
					// already done above
					// for(Size k = 6; k <= pose.residue(brsd).nheavyatoms(); ++k) {
					// 	if(!clashcheck.clash_check(pose.residue(brsd).xyz(k),brsd)) clash=true;
					// 	if(clash) break;
					// }
					// if(clash) { continue; }

					Vec const cenb = pose.residue(brsd).xyz("ZN");
					Vec orik1 = (pose.residue(brsd).xyz("NE1")-cenb).normalized();
					Vec orik2 = (pose.residue(brsd).xyz("NN1")-cenb).normalized();
					Vec const ligb1 = pose.residue(brsd).xyz("NE1");
					Vec const ligb2 = pose.residue(brsd).xyz("NN1");
					{
						Real ang = (90.0-angle_degrees(orik1,Vec(0,0,0),orik2))/2.0;
						orik1 = rotation_matrix_degrees(orik2.cross(orik1),ang)*orik1;
						orik2 = rotation_matrix_degrees(orik1.cross(orik2),ang)*orik2;
					}
					for(vector1<Size>::const_iterator iiter = scanres.begin(); iiter != scanres.end(); ++iiter) {
						Size irsd = *iiter;
						if(irsd==brsd) { continue; }
						Stub si(CBs[irsd],CAs[irsd],CAs[irsd]+(Ns[irsd]-Cs[irsd]));
						for(Size ide = 1; ide <= 2; ide++) {

							if(hallow0(irsd) > DUN_THRESH) continue; // if all chi1 / chi2 are clash or bad rot
							if(CBs[irsd].distance_squared(cenb) > ((ide==1)?hddcbmx2:hedcbmx2) ) { continue; }							// if CB too far, skip
							for(Size ich1 = 1; ich1 <= CHI1.size(); ++ich1) {
								if(hallow1(ich1,irsd) > DUN_THRESH) continue; // if all chi2 are clash or bad rot
								if(si.local2global(hci1cen(ide,ich1)).distance_squared(cenb) > ((ide==1)?hddc1mx2:hedc1mx2)) { continue; }							// if chi1 cen too far, skip
								for(Size ich2 = 1; ich2 <= CHI2.size(); ++ich2) {
									if(hallow2(ich2,ich1,irsd) > DUN_THRESH) continue;  // chi2 is clash or bad rot

									Vec const ceni = si.local2global(chi2cen(ide,ich1,ich2));
									if(ceni.distance_squared(cenb) > MXDSMTL) { continue; }
									Vec const orii = (si.local2global(chi2ori(ide,ich1,ich2))-ceni).normalized();
									bool his180 = false;
									{
										Real ang1 = numeric::conversions::degrees(acos(orii.dot(orik1)));
										Real ang2 = numeric::conversions::degrees(acos(orii.dot(orik2)));
										if( ang1 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL) ) { continue; }
										if( ang2 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang2 && ang2 < 180-MXAGMTL) ) { continue; }
										his180 = (ang1 > 135.0 || ang2 > 135.0);
									}
									if(!clashcheck.clash_check(si.local2global(chi2cen(ide,ich1,ich2)))) { continue; }

									pose.replace_residue(irsd,( ide==1 ? hsd : hse ).residue(1),true);
									pose.set_chi(1,irsd,CHI1[ich1]);
									pose.set_chi(2,irsd,CHI2[ich2]);
									Vec const ligi = pose.residue(irsd).xyz( (ide==1)?"ND1":"NE2" );

									clash = false;
									for(vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b)
										for(vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h)
											if( pose.residue(brsd).xyz(*b).distance_squared(pose.residue(irsd).xyz(*h)) < 7.0 ) clash=true;
									if(clash) { continue; }

									// for(Size a = 6; a <= pose.residue(irsd).nheavyatoms(); ++a) {
									// 	if(!clashcheck.clash_check(pose.residue(irsd).xyz(a),irsd)) { clash=true; if(clash) break; }
									// }
									// if(clash) { continue; }

									for(vector1<Size>::const_iterator jiter = scanres.begin(); jiter != scanres.end(); ++jiter) {
										Size jrsd = *jiter;
										if(jrsd <= irsd) { continue; }
										if(jrsd==brsd || jrsd==irsd) { continue; }
										Stub sj(CBs[jrsd],CAs[jrsd],CAs[jrsd]+(Ns[jrsd]-Cs[jrsd]));
										for(Size jde = 1; jde <= 2; jde++) {

											if(hallow0(jrsd) > DUN_THRESH) continue; // if all chi1 / chi2 are clash or bad rot
											if(CBs[jrsd].distance_squared(cenb) > ((jde==1)?hddcbmx2:hedcbmx2) ) { continue; }							// if CB too far, skip
											for(Size jch1 = 1; jch1 <= CHI1.size(); ++jch1) {
												if(hallow1(jch1,jrsd) > DUN_THRESH) continue; // if all chi2 are clash or bad rot
												if(sj.local2global(hci1cen(jde,jch1)).distance_squared(cenb) > ((jde==1)?hddc1mx2:hedc1mx2)) { continue; }							// if chi1 cen too far, skip
												for(Size jch2 = 1; jch2 <= CHI2.size(); ++jch2) {
													if(hallow2(jch2,jch1,jrsd) > DUN_THRESH) continue;  // chi2 is clash or bad rot

													Vec const cenj = sj.local2global(chi2cen(jde,jch1,jch2));
													if(cenj.distance_squared(cenb) > MXDSMTL) { continue; }
													Vec const orij = (sj.local2global(chi2ori(jde,jch1,jch2))-cenj).normalized();
													bool localhis180 = his180;
													bool hishis180 = false;
													{
														Real ang1 = numeric::conversions::degrees(acos(orij.dot(orii)));
														if( ang1 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL) ) { continue; }
														if( ang1 > 135.0 ) hishis180 = true;
													}
													if( !localhis180 && !hishis180 ) {
														// TR << "NOT his180, hishis180" << std::endl;
														Real ang1 = numeric::conversions::degrees(acos(orij.dot(orik1)));
														Real ang2 = numeric::conversions::degrees(acos(orij.dot(orik2)));
														if( ang1 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL) ) { continue; }
														if( ang2 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang2 && ang2 < 180-MXAGMTL) ) { continue; }
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
													Vec const ligj = pose.residue(jrsd).xyz( (jde==1)?"ND1":"NE2" );

													clash = false;
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

													// for(Size a = 6; a <= pose.residue(jrsd).nheavyatoms(); ++a) {
													// 	if(!clashcheck.clash_check(pose.residue(jrsd).xyz(a),jrsd)) { clash=true; if(clash) break; }
													// }
													// if(clash) { continue; }

													bhhout<< brsd << " " << kch1 << " " << kch2 <<
														" " << irsd << " " << ich1 << " " << ich2 << " " << ide <<
														" " << jrsd << " " << jch1 << " " << jch2 << " " << jde << std::endl;

													for(vector1<Size>::const_iterator eiter = scanres.begin(); eiter != scanres.end(); ++eiter) {
														Size ersd = *eiter;
														if(ersd==brsd||ersd==irsd||ersd==jrsd) { continue; }

														clash = false;
														Stub se(CBs[ersd],CAs[ersd],CAs[ersd]+(Ns[ersd]-Cs[ersd]));

														//if(eallow0(ersd) > DUN_THRESH) continue; // if all chi1 / chi2 are clash or bad rot
														if(CBs[ersd].distance_squared(cenb) > edscbmx2 ) { continue; }							// if CB too far, skip
														for(Size ech1 = 1; ech1 <= CHI1.size(); ++ech1) {
															//if(eallow1(ech1,ersd) > DUN_THRESH) continue; // if all chi2 are clash or bad rot
															if( se.local2global(eci1cen(ech1)).distance_squared(cenb) > edsc1mx2 ) { continue; }							// if chi1 cen too far, skip
															for(Size ech2 = 1; ech2 <= CHI2.size(); ++ech2) {
																//if(eallow2(ech2,ech1,ersd) > DUN_THRESH) continue;  // chi2 is clash or bad rot

																Vec const cene = se.local2global(gci2cen(ech1,ech2));
																if(cene.distance_squared(cenb) > MXDSMTL) { continue; }
																pose.replace_residue(ersd,glu.residue(1),true);
																pose.set_chi(1,ersd,CHI1[ech1]);
																pose.set_chi(2,ersd,CHI2[ech2]);
																Vec const CD = pose.residue(ersd).xyz("CD");
																Real ange1 = angle_degrees(CD,cenb,ligb1);
																Real ange2 = angle_degrees(CD,cenb,ligb2);
																Real ange3 = angle_degrees(CD,cenb,ligi);
																Real ange4 = angle_degrees(CD,cenb,ligj);
																Real MXAGtmp = min(MXAGMTL,22.5);
																if( ange1 < 90.0-MXAGMTL || (90.0+MXAGtmp < ange1 && ange1 < 135.0-MXAGtmp) || 135.0+MXAGMTL < ange1 ) continue;
																if( ange2 < 90.0-MXAGMTL || (90.0+MXAGtmp < ange2 && ange2 < 135.0-MXAGtmp) || 135.0+MXAGMTL < ange2 ) continue;
																if( ange3 < 90.0-MXAGMTL || (90.0+MXAGtmp < ange3 && ange3 < 135.0-MXAGtmp) || 135.0+MXAGMTL < ange3 ) continue;
																if( ange4 < 90.0-MXAGMTL || (90.0+MXAGtmp < ange4 && ange4 < 135.0-MXAGtmp) || 135.0+MXAGMTL < ange4 ) continue;
																int n90 = 0, n135 = 0;
																if( ange1 < 112.5 ) n90++; else n135++;
																if( ange2 < 112.5 ) n90++; else n135++;
																if( ange3 < 112.5 ) n90++; else n135++;
																if( ange4 < 112.5 ) n90++; else n135++;
																bool biglu = true;
																if( n90 != 2 || n135 != 2 ) biglu = false;
																Vec ea1(0,0,0),ea2(0,0,0);
																if(biglu) {
																	if( ange1 > 112.5 && ange2 > 112.5 ) { ea1 = ligb1; ea2 = ligb2; }
																	if( ange1 > 112.5 && ange3 > 112.5 ) { ea1 = ligb1; ea2 = ligi ; }
																	if( ange1 > 112.5 && ange4 > 112.5 ) { ea1 = ligb1; ea2 = ligj ; }
																	if( ange2 > 112.5 && ange3 > 112.5 ) { ea1 = ligb2; ea2 = ligi ; }
																	if( ange2 > 112.5 && ange4 > 112.5 ) { ea1 = ligb2; ea2 = ligj ; }
																	if( ange3 > 112.5 && ange4 > 112.5 ) { ea1 = ligi ; ea2 = ligj ; }
																	Real aer1 = dihedral_degrees(pose.residue(ersd).xyz("OE1"),pose.residue(ersd).xyz("CG"),CD,ea1);
																	Real aer2 = dihedral_degrees(pose.residue(ersd).xyz("OE2"),pose.residue(ersd).xyz("CG"),CD,ea2);
																	// TR << ange1 << " " << ange2 << " " << ange3 << " " << ange4 << " " << aer1 << " " << aer2 << std::endl;
																	pose.set_chi(3,ersd,pose.chi(3,ersd)+(aer1+aer2)/2.0);
																} else {
																	if( ange1 > 112.5 ) { ea1 = ligb1; }
																	if( ange2 > 112.5 ) { ea1 = ligb1; }
																	if( ange3 > 112.5 ) { ea1 = ligb1; }
																	if( ange4 > 112.5 ) { ea1 = ligb2; }
																	Real aer1 = dihedral_degrees(pose.residue(ersd).xyz("OE1"),pose.residue(ersd).xyz("CG"),CD,ea1);
																	pose.set_chi(3,ersd,pose.chi(3,ersd)+aer1);
																}

																clash = false;
																for(Size a = 6; a <= pose.residue(ersd).nheavyatoms(); ++a) {
																	Real TH = 8.0;
																	if( a > pose.residue(ersd).nheavyatoms()-2 ) TH = 4.0;
																	else if(!clashcheck.clash_check(pose.residue(ersd).xyz(a),ersd)) { clash=true; }
																	for(vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b)
																		if( pose.residue(brsd).xyz(*b).distance_squared(pose.residue(ersd).xyz(a)) < TH ) clash=true;
																	 for(vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h)
																		if( pose.residue(irsd).xyz(*h).distance_squared(pose.residue(ersd).xyz(a)) < TH ) clash=true;
																	for(vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h)
																		if( pose.residue(jrsd).xyz(*h).distance_squared(pose.residue(ersd).xyz(a)) < TH ) clash=true;

																	if(clash) break;
																}
																if(clash) { continue; }

																// if( ersd==76 && ech1==51 && ech2==29 ) {
																// 	pose.dump_pdb("test0.pdb");
																// 	utility_exit_with_message("tesljtl;sfj");
																// }

																bool overlap = false;
																Vec cen = (2*cenb + ceni + cenj + cene)/5.0;
																Vec ori1 = orik1;
																Vec ori2 = orik2;
																Vec ori3 = orii;
																Vec ori4 = orij;
																Vec ori5 = (CD-cene).normalized();
																Vec ori6 = (CD-cene).normalized();
																vector1<Vec> & foundcen ((biglu)? foundcenbg: foundcensg);
																vector1<Vec> & foundori1((biglu)?foundori1bg:foundori1sg);
																vector1<Vec> & foundori2((biglu)?foundori2bg:foundori2sg);
																vector1<Vec> & foundori3((biglu)?foundori3bg:foundori3sg);
																vector1<Vec> & foundori4((biglu)?foundori4bg:foundori4sg);
																vector1<Vec> & foundori5((biglu)?foundori5bg:foundori5sg);
															  vector1<Vec> & foundori6((biglu)?foundori6bg:foundori6sg);

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
																	opose.replace_residue(ersd,glu.residue(1),true);
																	opose.set_chi(1,ersd,CHI1[ech1]);
																	opose.set_chi(2,ersd,CHI2[ech2]);
																	opose.set_chi(3,ersd,pose.chi(3,ersd));

																	string outfname = utility::file_basename(infile)+"_FULL_"+string_of(brsd)+"_"+string_of(irsd)+"_"+string_of(jrsd)+"_"+string_of(ersd);
																	outfname += "_"+string_of(kch1)+"_"+string_of(kch2)+"_"+string_of(ich1)+"_"+string_of(ich2)+"_"+string_of(ide);
																	outfname += "_"+string_of(jch1)+"_"+string_of(jch2)+"_"+string_of(jde)+"_"+string_of(ech1)+"_"+string_of(ech2);
																	if(biglu) outfname += "_BIGLU";
																	outfname += ".pdb.gz";
																	TR << "HIT! " << outfname << std::endl;

																	sf->score(opose);
																	//Real dbb = opose.energies().residue_total_energies(brsd)[core::scoring::fa_dun];
																	Real dbi = opose.energies().residue_total_energies(irsd)[core::scoring::fa_dun];
																	Real dbj = opose.energies().residue_total_energies(jrsd)[core::scoring::fa_dun];
																	Real dbe = opose.energies().residue_total_energies(ersd)[core::scoring::fa_dun];
																	core::id::AtomID_Map<Real> atom_sasa;
																	vector1<Real> rsd_sasa;
																	core::scoring::calc_per_atom_sasa( opose, atom_sasa, rsd_sasa, 2.0, false );
																	TR << "DUN " << irsd << " " << ich1 << " " << ich2 << " "
																		 << dbi << " " << hallow2(ich2,ich1,irsd) << " "
																		 << dbj << " " << hallow2(jch2,jch1,jrsd) << " "
																		 << dbe << " " << eallow2(ech2,ech1,ersd) << " " << std::endl;
																	if( fabs(dbi-hallow2(ich2,ich1,irsd)) > 0.001 || fabs(dbj-hallow2(jch2,jch1,jrsd)) > 0.001 ) {
																		TR << "RUN " << opose.phi(irsd) << " " << opose.psi(irsd) << " " << opose.omega(irsd) << " " << opose.chi(1,irsd) << " " << opose.chi(2,irsd) << std::endl;
																		opose.set_chi(1,irsd,opose.chi(1,irsd)+360.0);
																		opose.set_chi(2,irsd,opose.chi(2,irsd)+360.0);
																		opose.energies().clear_energies();
																		sf->score(opose);
																		TR << opose.chi(1,irsd) << " " << opose.chi(2,irsd) << " " << opose.residue(irsd).name() << " " << opose.energies().residue_total_energies(irsd)[core::scoring::fa_dun] << std::endl;
																		opose.set_chi(1,irsd,opose.chi(1,irsd)-360.0);
																		opose.set_chi(2,irsd,opose.chi(2,irsd)-360.0);
																		opose.energies().clear_energies();
																		sf->score(opose);
																		TR << opose.chi(1,irsd) << " " << opose.chi(2,irsd) << " " << opose.residue(irsd).name() << " " << opose.energies().residue_total_energies(irsd)[core::scoring::fa_dun] << std::endl;
																		opose.set_chi(1,irsd,opose.chi(1,irsd)-360.0);
																		opose.set_chi(2,irsd,opose.chi(2,irsd)-360.0);
																		opose.energies().clear_energies();
																		sf->score(opose);
																		TR << opose.chi(1,irsd) << " " << opose.chi(2,irsd) << " " << opose.residue(irsd).name() << " " << opose.energies().residue_total_energies(irsd)[core::scoring::fa_dun] << std::endl;
																		opose.set_chi(1,irsd,opose.chi(1,irsd)+360.0);
																		opose.set_chi(2,irsd,opose.chi(2,irsd)+360.0);
																		opose.energies().clear_energies();
																		sf->score(opose);
																		TR << opose.chi(1,irsd) << " " << opose.chi(2,irsd) << " " << opose.residue(irsd).name() << " " << opose.energies().residue_total_energies(irsd)[core::scoring::fa_dun] << std::endl;




																		opose.dump_pdb("test.pdb");
																		utility_exit_with_message("DUN!");
																	}

																	core::io::silent::SilentStructOP ss_out_all( new core::io::silent::ScoreFileSilentStruct );
																	ss_out_all->fill_struct(opose,outfname);
																	ss_out_all->add_energy( "bpy_chi1" , CHI1[kch1] );
																	ss_out_all->add_energy( "bpy_chi2" , CHI2[kch2] );
																	ss_out_all->add_energy( "dun_his1", dbi );
																	ss_out_all->add_energy( "dun_his2", dbj );
																	ss_out_all->add_energy( "dun_glu" , dbe );
																	ss_out_all->add_energy( "sasa_bpy" , rsd_sasa[brsd] );
																	ss_out_all->add_energy( "sasa_his1", rsd_sasa[irsd] );
																	ss_out_all->add_energy( "sasa_his2", rsd_sasa[jrsd] );
																	ss_out_all->add_energy( "sasa_glu" , rsd_sasa[ersd] );
																	sfd.write_silent_struct( *ss_out_all, option[ out::file::silent ]() );

																	// utility::io::ozstream otmp(option[out::file::o]()+"/"+outfname);
																	// otmp << "REMARK 666 MATCH TEMPLATE A BPY "+ObjexxFCL::format::I(4,brsd)+" MATCH MOTIF A HIS "+ObjexxFCL::format::I(4,irsd)+"  1  1" << std::endl;
																	// otmp << "REMARK 666 MATCH TEMPLATE A BPY "+ObjexxFCL::format::I(4,brsd)+" MATCH MOTIF A HIS "+ObjexxFCL::format::I(4,jrsd)+"  2  1" << std::endl;
																	// otmp << "REMARK 666 MATCH TEMPLATE A BPY "+ObjexxFCL::format::I(4,brsd)+" MATCH MOTIF A GLU "+ObjexxFCL::format::I(4,ersd)+"  3  1" << std::endl;
																	// opose.dump_pdb(otmp);
																	// otmp.close();

																	// static int outcount = 0;
																	// outcount++;
																	// if(outcount == 2)	utility_exit_with_message("test");
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
			}
		}

		bhhout.close();
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


