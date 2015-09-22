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
#include <core/io/pdb/pose_io.hh>
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

static THREAD_LOCAL basic::Tracer TR( "sf4_hash" );

struct MyRT {
	Vec trans;
	Real w,x,y,z;
	MyRT() : trans(0,0,0),w(0),x(0),y(0),z(0) {}
	MyRT(Vec trans_in, Vec axis, Real ang) {
		from_axis(trans_in,axis,ang);
	}
	MyRT(Vec trans_in, Mat rot_in) {
		from_rot(trans_in,rot_in);
	}
	void from_rot(Vec trans_in, Mat rot) {
		Real angle;
		Vec axis = rotation_axis( rot, angle );
		from_axis(trans_in,axis,angle);
	}
	void from_axis(Vec trans_in, Vec axis, Real ang) {
		trans = trans_in;
		// w = cos(ang/2.0);
		// x = axis.x();
		// y = axis.y();
		// z = axis.z();
		// Real l = sqrt(w*w+x*x+y*y+z*z);
		// w /= l;
		// x /= l;
		// y /= l;
		// z /= l;
		w = ang;
		x = axis.x();
		y = axis.y();
		z = axis.z();
	}
	Real distance_squared(MyRT o) {
		return (w-o.w)*(w-o.w)+(x-o.x)*(x-o.x)+(y-o.y)*(y-o.y)+(z-o.z)*(z-o.z);
	}
};

std::string SAN1("CB"),SAN2("CA"),SAN3("N");

void alncys(core::pose::Pose & pose, Real c1, Real c2, Real c3) {
 	pose.set_chi(1,1,c1);
	pose.set_chi(2,1,c2);
	trans_pose(pose,Vec(2.225,2.225,2.225)-pose.residue(1).xyz("SG") );
	alignaxis (pose,Vec(2.225,2.225,2.225),pose.residue(1).xyz("SG")-pose.residue(1).xyz("HG"),pose.residue(1).xyz("SG"));
	Real ang = dihedral_degrees(Vec(1,-1,-1),Vec(0,0,0),pose.residue(1).xyz("SG"),pose.residue(1).xyz("CB"));
	rot_pose(pose, Vec(1,1,1), c3-ang );
}

void get_cys_rts(utility::vector1<utility::vector1<core::Real> > & CHI, utility::vector1<Stub> & STs, utility::vector1<utility::vector1<MyRT> > & RTs) {

	core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	Vec F1( 0.961, 0.961, 0.961), F2 (0.961,-0.961,-0.961), F3(-0.961,-0.961, 0.961), F4(-0.961, 0.961,-0.961);
	Vec S1(-1.279,-1.279,-1.279), S2(-1.279, 1.279, 1.279), S3( 1.279, 1.279,-1.279), S4( 1.279,-1.279, 1.279);
	cout<<"HETATM"<<I(5,1)<<' '<<" S1 "<<' ' <<	"SF4"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,S1.x())<<F(8,3,S1.y())<<F(8,3,S1.z())<<F(6,2,1.0)<<F(6,2,1.0)<<"           S\n";
	cout<<"HETATM"<<I(5,2)<<' '<<" S2 "<<' ' <<	"SF4"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,S2.x())<<F(8,3,S2.y())<<F(8,3,S2.z())<<F(6,2,1.0)<<F(6,2,1.0)<<"           S\n";
	cout<<"HETATM"<<I(5,3)<<' '<<" S3 "<<' ' <<	"SF4"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,S3.x())<<F(8,3,S3.y())<<F(8,3,S3.z())<<F(6,2,1.0)<<F(6,2,1.0)<<"           S\n";
	cout<<"HETATM"<<I(5,4)<<' '<<" S4 "<<' ' <<	"SF4"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,S4.x())<<F(8,3,S4.y())<<F(8,3,S4.z())<<F(6,2,1.0)<<F(6,2,1.0)<<"           S\n";
	cout<<"HETATM"<<I(5,5)<<' '<<" F1 "<<' ' <<	"SF4"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,F1.x())<<F(8,3,F1.y())<<F(8,3,F1.z())<<F(6,2,1.0)<<F(6,2,1.0)<<"          FE\n";
	cout<<"HETATM"<<I(5,6)<<' '<<" F2 "<<' ' <<	"SF4"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,F2.x())<<F(8,3,F2.y())<<F(8,3,F2.z())<<F(6,2,1.0)<<F(6,2,1.0)<<"          FE\n";
	cout<<"HETATM"<<I(5,7)<<' '<<" F3 "<<' ' <<	"SF4"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,F3.x())<<F(8,3,F3.y())<<F(8,3,F3.z())<<F(6,2,1.0)<<F(6,2,1.0)<<"          FE\n";
	cout<<"HETATM"<<I(5,8)<<' '<<" F4 "<<' ' <<	"SF4"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,F4.x())<<F(8,3,F4.y())<<F(8,3,F4.z())<<F(6,2,1.0)<<F(6,2,1.0)<<"          FE\n";

	vector1<Real> CHI1(3); CHI1[1] =  -60.0; CHI1[2] = 60.0; CHI1[3] = 180.0;
	vector1<Real> CHI2(3); CHI2[1] = -110.0; CHI2[2] =  3.0; CHI2[3] =  65.0;
	vector1<Real> CHI3(3); CHI3[1] = -120.0; CHI3[2] =  0.0; CHI3[3] = 120.0;

	Pose cys1,cys2;
	core::pose::make_pose_from_sequence(cys1,"C",*fa_residue_set);
	core::pose::make_pose_from_sequence(cys2,"C",*fa_residue_set);
	core::pose::remove_lower_terminus_type_from_pose_residue(cys1,1);
	core::pose::remove_lower_terminus_type_from_pose_residue(cys2,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(cys1,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(cys2,1);
	// cys1.dump_pdb("test.pdb");

	CHI.clear();
	for(Size ich = 1; ich <= 27; ich++) {
		Real ichi1 = CHI1[((ich-1)/1)%3+1];
		Real ichi2 = CHI2[((ich-1)/3)%3+1];
		Real ichi3 = CHI3[((ich-1)/9)%3+1];
		alncys(cys1,ichi1,ichi2,ichi3);
		bool clash = false;
		for(Size i = 1; i <= 3; ++i) {
			if(cys1.residue(1).xyz(i).distance_squared(S1) < 9.0 ) clash = true;
			if(cys1.residue(1).xyz(i).distance_squared(S2) < 9.0 ) clash = true;
			if(cys1.residue(1).xyz(i).distance_squared(S3) < 9.0 ) clash = true;
			if(cys1.residue(1).xyz(i).distance_squared(S4) < 9.0 ) clash = true;
			if(cys1.residue(1).xyz(i).distance_squared(F1) < 9.0 ) clash = true;
			if(cys1.residue(1).xyz(i).distance_squared(F2) < 9.0 ) clash = true;
			if(cys1.residue(1).xyz(i).distance_squared(F3) < 9.0 ) clash = true;
			if(cys1.residue(1).xyz(i).distance_squared(F4) < 9.0 ) clash = true;
		}
		if(clash) continue;
		vector1<Real> chi; chi.push_back(ichi1); chi.push_back(ichi2); chi.push_back(ichi3);
		CHI.push_back(chi);
		//		cys1.dump_pdb("startcys_"+string_of(CHI.size())+".pdb");
		TR << "rot " << CHI.size() << " " << ichi1 << " " << ichi2 << " " << ichi3 << std::endl;
	}


	RTs.resize(CHI.size());
	STs.resize(CHI.size());
	for(Size i = 1; i <= CHI.size(); i++) RTs[i].resize(CHI.size());

	for(Size ich = 1; ich <= CHI.size(); ich++) {
		Real ichi1 = CHI[ich][1];
		Real ichi2 = CHI[ich][2];
		Real ichi3 = CHI[ich][3];
		TR << "ich " << ich << " " << ichi1 << " " << ichi2 << " " << ichi3 << std::endl;
		{
			cys1.set_chi(1,1,ichi1);
			cys1.set_chi(2,1,ichi2);
			trans_pose(cys1,Vec(2.225,2.225,2.225)-cys1.residue(1).xyz("SG") );
			alignaxis (cys1,Vec(2.225,2.225,2.225),cys1.residue(1).xyz("SG")-cys1.residue(1).xyz("HG"),cys1.residue(1).xyz("SG"));
			Real ang = dihedral_degrees(Vec(1,-1,-1),Vec(0,0,0),cys1.residue(1).xyz("SG"),cys1.residue(1).xyz("CB"));
			rot_pose(cys1, Vec(1,1,1), ichi3-ang );
			//cys1.dump_pdb("cys1_"+string_of(ich)+".pdb");
		}
		Stub stub1(cys1.residue(1).xyz(SAN1),cys1.residue(1).xyz(SAN2),cys1.residue(1).xyz(SAN3));
		STs[ich] = stub1;
		for(Size jch = 1; jch <= CHI.size(); jch++) {
			Real jchi1 = CHI[jch][1];
			Real jchi2 = CHI[jch][2];
			Real jchi3 = CHI[jch][3];
			// TR << "ich " << ich << " " << chi1 << " " << chi2 << " " << chi3 << std::endl;
			{
				cys2.set_chi(1,1,jchi1);
				cys2.set_chi(2,1,jchi2);
				trans_pose(cys2,Vec(2.225,-2.225,-2.225)-cys2.residue(1).xyz("SG") );
				alignaxis (cys2,Vec(2.225,-2.225,-2.225),cys2.residue(1).xyz("SG")-cys2.residue(1).xyz("HG"),cys2.residue(1).xyz("SG"));
				Real ang = dihedral_degrees(Vec(1,1,1),Vec(0,0,0),cys2.residue(1).xyz("SG"),cys2.residue(1).xyz("CB"));
				rot_pose(cys2, Vec(1,-1,-1), jchi3-ang );
				//			cys2.dump_pdb("cys2_"+string_of(jch)+".pdb");
			}
			Stub stub2(cys2.residue(1).xyz(SAN1),cys2.residue(1).xyz(SAN2),cys2.residue(1).xyz(SAN3));
			core::kinematics::RT rt;
			rt.from_stubs(stub1,stub2);
			RTs[ich][jch].from_rot(rt.get_translation(),rt.get_rotation());
		// rt.get_rotation(), RTs[ich][jch].angle );
			// TR << "axisang " << RTs[ich][jch].angle << " " << RTs[ich][jch].axis << std::endl;
		}
	}
	//utility_exit_with_message("debug");
}

void run_sf4h() {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::id;

	vector1<Real> CHI1(3); CHI1[1] =  -60.0; CHI1[2] = 60.0; CHI1[3] = 180.0;
	vector1<Real> CHI2(3); CHI2[1] = -110.0; CHI2[2] =  3.0; CHI2[3] =  65.0;
	vector1<Real> CHI3(3); CHI3[1] = -120.0; CHI3[2] =  0.0; CHI3[3] = 120.0;

	core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	core::chemical::ResidueTypeSetCAP cen_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );

	vector1<vector1<Real> > CHI;
	vector1<Stub> STs;
	vector1<vector1<MyRT> > RTs;
	get_cys_rts(CHI,STs,RTs);
	Vec lb(9e9,9e9,9e9),ub(-9e9,-9e9,-9e9);
	for(Size i = 1; i <= RTs.size(); i++) {
		for(Size j = 1; j <= RTs.size(); j++) {
			Vec t = RTs[i][j].trans; lb.min(t);	ub.max(t);
		}
	}
	TR << "min " << lb << std::endl;
	TR << "max " << ub << std::endl;

	core::conformation::ResidueOP cys = core::conformation::ResidueFactory::create_residue(fa_residue_set->name_map("CYS"));

	vector1<string> infiles = option[in::file::s]();
	for(Size ifile = 1; ifile <= infiles.size(); ++ifile) {
		string infile = utility::file_basename(infiles[ifile]);
		TR << "scaning pdb " << infiles[ifile] << std::endl;
		Pose pose;
		core::import_pose::pose_from_pdb(pose,*fa_residue_set,infiles[ifile]);
		Size gap = 0;
		for(Size ir = 1; ir <= pose.n_residue()-2; ++ir) {
			if(pose.residue(ir).name3()=="GLY" || pose.residue(ir).name3()=="PRO") continue;
			Stub stub1(pose.residue(ir).xyz(SAN1),pose.residue(ir).xyz(SAN2),pose.residue(ir).xyz(SAN3));
			for(Size jr = ir+2; jr <= ir+15; ++jr) {
				if(jr > pose.n_residue()) break;
				if(pose.residue(jr).name3()=="GLY" || pose.residue(jr).name3()=="PRO") continue;
				if( pose.residue(jr-1).xyz("C").distance_squared(pose.residue(jr).xyz("N")) > 4.0 ) {
					gap = jr-1;
				}
				if( ir <= gap && jr > gap ) continue;
				Stub stub2(pose.residue(jr).xyz(SAN1),pose.residue(jr).xyz(SAN2),pose.residue(jr).xyz(SAN3));
				core::kinematics::RT rt;
				rt.from_stubs(stub1,stub2);
				Vec t = rt.get_translation();
				//				if( t.x()+0.5 < lb.x() || ub.x() < t.x()-0.5 || t.y()+0.5 < lb.y() || ub.y() < t.y()-0.5 || t.z()+0.5 < lb.z() || ub.z() < t.z()-0.5 ) continue;
				MyRT mrt(rt.get_translation(),rt.get_rotation());
				for(Size ich = 1; ich <= RTs.size(); ich++) {
					for(Size jch = 1; jch <= RTs.size(); jch++) {
						if(RTs[ich][jch].trans.distance_squared(mrt.trans) > 0.5) continue;
						// if(RTs[ich][jch].distance_squared(mrt) > 0.0301537) continue;

						Vec axis1 = Vec(RTs[ich][jch].x,RTs[ich][jch].y,RTs[ich][jch].z).normalized();
						Vec axis2 = Vec(          mrt.x,          mrt.y,          mrt.z).normalized();

						if(axis1.dot(axis2) < 0.95) continue;
						if(fabs(RTs[ich][jch].w-mrt.w) > 0.2) continue;
						TR << "HIT " << infile << " " << ir << " " << jr << " " << ich << " " << jch << std::endl;
						TR << "trans 1" << mrt.trans << std::endl;
						TR << "trans 2" << RTs[ich][jch].trans << std::endl;
						TR << "rot1 " << axis1 << std::endl;
						TR << "rot2 " << axis2 << std::endl;
						TR << "ANG " << RTs[ich][jch].w << " " << mrt.w << std::endl;
						Pose sub;
						sub.append_residue_by_jump(pose.residue(ir),1);
						for(Size kr = ir+1; kr <= jr; ++kr) sub.append_residue_by_bond(pose.residue(kr));
						sub.replace_residue(1,*cys,true);
						cys->set_chi(1,CHI[jch][1]);
						cys->set_chi(2,CHI[jch][2]);
						sub.replace_residue(sub.n_residue(),*cys,true);

						cys->set_chi(1,CHI[ich][1]);
						cys->set_chi(2,CHI[ich][2]);
						trans_pose(sub,Vec(2.225,2.225,2.225)-sub.residue(1).xyz("SG") );
						alignaxis (sub,Vec(2.225,2.225,2.225),sub.residue(1).xyz("SG")-sub.residue(1).xyz("HG"),sub.residue(1).xyz("SG"));
						Real ang = dihedral_degrees(Vec(1,-1,-1),Vec(0,0,0),sub.residue(1).xyz("SG"),sub.residue(1).xyz("CB"));
						rot_pose(sub, Vec(1,1,1), CHI[ich][3]-ang );

						std::string fn = string_of(ich)+"_"+infile+"_"+string_of(ir)+"_"+string_of(jr)+".pdb";
						TR << "dumping pdb " << fn << std::endl;
						sub.dump_pdb(fn);
						utility_exit_with_message("debug");
					}
				}
			}
		}
	}
}


int main (int argc, char *argv[]) {

	try {

	devel::init(argc,argv);
	run_sf4h();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


