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
#include <basic/options/option_macros.hh>
//#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <protocols/sic_dock/xyzStripeHashPoseWithMeta.hh>
#include <sstream>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>

using core::Size;
using core::Real;
using ObjexxFCL::string_of;


OPT_1GRP_KEY( Real , gendsf, rms_cutoff )

void register_options() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		NEW_OPT( gendsf::rms_cutoff, "rms rudundency cutoff", 2.0 );
}


Real norm_degrees(Real x) {
	if(x < 0.0) return x+360.0;
	return x;
}

core::pack::rotamer_set::RotamerSetOP
get_rotamers(
	core::pose::Pose & pose,
  	core::Size irsd,
  	bool ex1,
  	bool ex2
) {
	core::pack::rotamer_set::RotamerSetOP rotset;
	core::scoring::ScoreFunction dummy_sfxn;
	dummy_sfxn( pose );
	core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
	dummy_task->initialize_from_command_line();
	dummy_task->nonconst_residue_task( irsd ).restrict_to_repacking();
	dummy_task->nonconst_residue_task( irsd ).and_extrachi_cutoff(0);
	dummy_task->nonconst_residue_task( irsd ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
	dummy_task->nonconst_residue_task( irsd ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
	dummy_task->nonconst_residue_task( irsd ).or_ex1( ex1 ); //since we only want rotamers for the specified restype
	dummy_task->nonconst_residue_task( irsd ).or_ex2( ex2 ); //since we only want rotamers for the specified restype
	core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
	core::pack::rotamer_set::RotamerSetFactory rsf;
	rotset = rsf.create_rotamer_set( pose.residue( irsd ) );
	rotset->set_resid( irsd );
	rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
	return rotset;
}

// get chi lists with more canonical rotamers first
vector1<vector1<Real> >
get_chis_complicated(
  core::pose::Pose & pose,
  core::Size irsd
){
	vector1<vector1<Real> > chis;
	core::pack::rotamer_set::RotamerSetOP rotset1[16];
	core::pack::rotamer_set::RotamerSetOP rotset2[16];
	core::pack::rotamer_set::RotamerSetOP rotset00 = get_rotamers(pose,irsd,false,false);
	core::pack::rotamer_set::RotamerSetOP rotset10 = get_rotamers(pose,irsd, true,false);
	core::pack::rotamer_set::RotamerSetOP rotset01 = get_rotamers(pose,irsd,false, true);
	core::pack::rotamer_set::RotamerSetOP rotset11 = get_rotamers(pose,irsd, true, true);
	rotset1[ 0] = rotset00; rotset2[ 0] = rotset00;
	rotset1[ 1] = rotset00; rotset2[ 1] = rotset10;
	rotset1[ 2] = rotset10; rotset2[ 2] = rotset00;
	rotset1[ 3] = rotset10; rotset2[ 3] = rotset10;
	rotset1[ 4] = rotset00; rotset2[ 4] = rotset01;
	rotset1[ 5] = rotset01; rotset2[ 5] = rotset00;
	rotset1[ 6] = rotset10; rotset2[ 6] = rotset01;
	rotset1[ 7] = rotset01; rotset2[ 7] = rotset10;
	rotset1[ 8] = rotset01; rotset2[ 8] = rotset01;
	rotset1[ 9] = rotset00; rotset2[ 9] = rotset11;
	rotset1[10] = rotset11; rotset2[10] = rotset00;
	rotset1[11] = rotset10; rotset2[11] = rotset11;
	rotset1[12] = rotset11; rotset2[12] = rotset10;
	rotset1[13] = rotset01; rotset2[13] = rotset11;
	rotset1[14] = rotset11; rotset2[14] = rotset01;
	rotset1[15] = rotset11; rotset2[15] = rotset11;
	for(Size irs = 0; irs < 16; ++irs) {
		for(Size irot = 1; irot <= rotset1[irs]->num_rotamers(); ++irot) {
			core::Real chi1i = rotset1[irs]->rotamer(irot)->chi(1);
			core::Real chi2i = rotset1[irs]->rotamer(irot)->chi(2);
			if( fabs(chi2i) < 30.0 ) continue;
			for(Size jrot = 1; jrot <= rotset2[irs]->num_rotamers(); ++jrot) {
				core::Real chi1j = rotset2[irs]->rotamer(jrot)->chi(1);
				core::Real chi2j = rotset2[irs]->rotamer(jrot)->chi(2);
				if( fabs(chi2j) < 30.0 ) continue;
				bool seenit = false;
				for(Size i = 1; i <= chis.size(); ++i) {
					if( fabs(norm_degrees(chis[i][1])-norm_degrees(chi1i)) < 5.0 && fabs(norm_degrees(chis[i][2])-norm_degrees(chi2i)) < 5.0 &&
						fabs(norm_degrees(chis[i][3])-norm_degrees(chi1j)) < 5.0 && fabs(norm_degrees(chis[i][4])-norm_degrees(chi2j)) < 5.0 ){
						seenit = true;
					}
				}
				if(!seenit) {
					vector1<Real> tmp(4);
					tmp[1] = chi1i;
					tmp[2] = chi2i;
					tmp[3] = chi1j;
					tmp[4] = chi2j;
					chis.push_back(tmp);
				}
			}
		}
	}
	return chis;
}

vector1<vector1<Real> >
get_chis(
  core::pose::Pose & pose,
  core::Size irsd
){
	vector1<vector1<Real> > chis;
	core::pack::rotamer_set::RotamerSetOP rotset00 = get_rotamers(pose,irsd,false,false);
	for(Size irot = 1; irot <= rotset00->num_rotamers(); ++irot) {
		core::Real chi1i = rotset00->rotamer(irot)->chi(1);
		core::Real chi2i = rotset00->rotamer(irot)->chi(2);
		if( fabs(chi2i) < 30.0 ) continue;
		for(Size jrot = 1; jrot <= rotset00->num_rotamers(); ++jrot) {
			core::Real chi1j = rotset00->rotamer(jrot)->chi(1);
			core::Real chi2j = rotset00->rotamer(jrot)->chi(2);
			if( fabs(chi2j) < 30.0 ) continue;
			vector1<Real> tmp(4);
			tmp[1] = chi1i;
			tmp[2] = chi2i;
			tmp[3] = chi1j;
			tmp[4] = chi2j;
			chis.push_back(tmp);
		}
	}
	return chis;
}


int get_N(core::pose::Pose const & inp) {
	Size c = inp.chain(1);
	Size i = 2;
	for(i = 2; i <= inp.n_residue(); ++i) {
		if(inp.chain(i) != (int)c) {
			break;
		}
	}
	Size N = i-1;
	if( inp.n_residue() % N != 0 ) {
		utility_exit_with_message("N does not divide n_residue! "
			+ObjexxFCL::string_of(N)+" "+ObjexxFCL::string_of(inp.n_residue()) );
	} else if( inp.residue(N).name() != inp.residue(2*N).name() ) {
		utility_exit_with_message("res N & 2*N not same type!!!");
	} else if( inp.residue(N).nchi() > 0 && fabs(norm_degrees(inp.residue(N).chi(1)) - norm_degrees(inp.residue(2*N).chi(1))) > 1.0 ) {
		utility_exit_with_message("res N & 2*N not same chi1: "+
			ObjexxFCL::string_of(inp.residue(N).chi(1))+" "+ObjexxFCL::string_of(inp.residue(2*N).chi(1)));
	} else if( inp.residue(N).nchi() > 1 && fabs(norm_degrees(inp.residue(N).chi(2)) - norm_degrees(inp.residue(2*N).chi(2))) > 1.0 ) {
		utility_exit_with_message("res N & 2*N not same chi2: "+
			ObjexxFCL::string_of(inp.residue(N).chi(2))+" "+ObjexxFCL::string_of(inp.residue(2*N).chi(2)));
	} else if( inp.residue(N).nchi() > 2 && fabs(norm_degrees(inp.residue(N).chi(3)) - norm_degrees(inp.residue(2*N).chi(3))) > 1.0 ) {
		utility_exit_with_message("res N & 2*N not same chi3: "+
			ObjexxFCL::string_of(inp.residue(N).chi(3))+" "+ObjexxFCL::string_of(inp.residue(2*N).chi(3)));
	} else if( inp.residue(N).nchi() > 3 && fabs(norm_degrees(inp.residue(N).chi(4)) - norm_degrees(inp.residue(2*N).chi(4))) > 1.0 ) {
		utility_exit_with_message("res N & 2*N not same chi4: "+
			ObjexxFCL::string_of(inp.residue(N).chi(4))+" "+ObjexxFCL::string_of(inp.residue(2*N).chi(4)));
	}
	return N;
}

core::pose::Pose make_two_trimers(core::pose::Pose const & inp, Size icys, Size N) {
	core::pose::Pose pose;
	for(Size i = 1; i <= 3*N; ++i) {
		if(inp.residue(i).is_lower_terminus()||inp.residue(i).is_ligand()) {
			pose.append_residue_by_jump(inp.residue(i),1);
		} else {
			pose.append_residue_by_bond(inp.residue(i));
		}
	}

	// get sym xform from inp sub1 / sub4
	Vec cen = (center_of_geom(inp,1,N)+center_of_geom(inp,3*N+1,4*N))/2.0;
	Vec a(0,0,0);
	for(Size i = 1; i <= N; ++i) {
		Vec tmp = Vec(0,0,0);
		tmp += inp.xyz(AtomID(2,i))+inp.xyz(AtomID(2,3*N+i));
		tmp = tmp/(Real)2.0 - cen;
		if( a.length() > 0.0001 && a.dot(tmp) < 0 ) tmp = -tmp;
		a += tmp;
	}
	a.normalize();
	Mat rot = rotation_matrix_degrees(a,180.0);

	core::pose::Pose partner(pose);
	rot_pose(partner,rot,cen);
	for(Size i = 1; i <= partner.n_residue(); ++i) {
		if(partner.residue(i).is_lower_terminus()||partner.residue(i).is_ligand()) {
			pose.append_residue_by_jump(partner.residue(i),1);
		} else {
			pose.append_residue_by_bond(partner.residue(i));
		}
	}
	assert(pose.n_residue()==6*N);

	std::cout << pose.fold_tree() << std::endl;

	using namespace core::kinematics;
	FoldTree ft = pose.fold_tree();
	ObjexxFCL::FArray2D<int> jumps(2,5);
	ObjexxFCL::FArray1D<int> cuts(5);
	for(int i = 1; i <= 5; ++i) cuts(i) = N*i;
	jumps(1,1) = 1*N; jumps(2,1) = 2*N;
	jumps(1,2) = 1*N; jumps(2,2) = 3*N;
	jumps(1,3) = icys; jumps(2,3) = 3*N+icys;
	jumps(1,4) = 4*N; jumps(2,4) = 5*N;
	jumps(1,5) = 4*N; jumps(2,5) = 6*N;
	ft.tree_from_jumps_and_cuts(6*N,5,jumps,cuts);
	ft.set_jump_atoms(3,"SG","SG");

	pose.fold_tree(ft);

	std::cout << pose.fold_tree() << std::endl;

	return pose;
}

void set_disulf(core::pose::Pose & pose, Size icys, Size N, Size Nsub, Real chi1, Real chi2, Real chiss, Real chi2b, Real chi1b) {
	Real delta_chi1  = chi1  - dihedral_degrees(pose.xyz(AtomID(1,0*N+icys)),pose.xyz(AtomID(2,0*N+icys)),pose.xyz(AtomID(5,0*N+icys)),pose.xyz(AtomID(6,0*N+icys)));
	Real delta_chi2  = chi2  - dihedral_degrees(pose.xyz(AtomID(2,0*N+icys)),pose.xyz(AtomID(5,0*N+icys)),pose.xyz(AtomID(6,0*N+icys)),pose.xyz(AtomID(6,Nsub*N+icys)));
	Real delta_chiss = chiss - dihedral_degrees(pose.xyz(AtomID(5,0*N+icys)),pose.xyz(AtomID(6,0*N+icys)),pose.xyz(AtomID(6,Nsub*N+icys)),pose.xyz(AtomID(5,Nsub*N+icys)));
	Real delta_chi2b = chi2b - dihedral_degrees(pose.xyz(AtomID(2,Nsub*N+icys)),pose.xyz(AtomID(5,Nsub*N+icys)),pose.xyz(AtomID(6,Nsub*N+icys)),pose.xyz(AtomID(6,0*N+icys)));
	Real delta_chi1b = chi1b - dihedral_degrees(pose.xyz(AtomID(1,Nsub*N+icys)),pose.xyz(AtomID(2,Nsub*N+icys)),pose.xyz(AtomID(5,Nsub*N+icys)),pose.xyz(AtomID(6,Nsub*N+icys)));

	pose.set_chi(1,icys,chi1);
	Vec CA = pose.xyz(AtomID(2,icys));
	Vec CB = pose.xyz(AtomID(5,icys));
	rot_pose( pose, rotation_matrix_degrees(CB -CA ,delta_chi1 ), CA , Nsub*N+1, 2*Nsub*N );
	pose.set_chi(2,icys,chi2);
	Vec SG = pose.xyz(AtomID(6,icys));
	rot_pose( pose, rotation_matrix_degrees(SG -CB ,delta_chi2 ), CB , Nsub*N+1, 2*Nsub*N );
	Vec SGb = pose.xyz(AtomID(6,Nsub*N+icys));
	rot_pose( pose, rotation_matrix_degrees(SGb-SG ,delta_chiss), SG , Nsub*N+1, 2*Nsub*N );
	Vec CBb = pose.xyz(AtomID(5,Nsub*N+icys));
	pose.set_chi(2,Nsub*N+icys,chi2b);
	rot_pose( pose, rotation_matrix_degrees(CBb-SGb,delta_chi2b), SGb, Nsub*N+1, 2*Nsub*N );
	Vec CAb = pose.xyz(AtomID(2,Nsub*N+icys));
	pose.set_chi(1,Nsub*N+icys,chi1b);
	rot_pose( pose, rotation_matrix_degrees(CAb-CBb,delta_chi1b), CBb, Nsub*N+1, 2*Nsub*N );
}

bool clash_check_bb(core::pose::Pose & pose, Size icys, Size N, Size Nsub) {
	for(Size ir = 1; ir <= Nsub*N; ++ir) {
		Size inatom = pose.residue(ir).aa() == core::chemical::aa_gly ? 4 : 5;
		for(Size ia = 1; ia <= inatom; ++ia) {
			Vec iX = pose.xyz(AtomID(ia,ir));
			if( ir!=icys && pose.xyz(AtomID(6,Nsub*N+icys)).distance_squared(iX) < 7.0 ) return false;
			for(Size jr = Nsub*N+1; jr <= 2*Nsub*N; ++jr) {
				Size jnatom = pose.residue(jr).aa() == core::chemical::aa_gly ? 4 : 5;
				for(Size ja = 1; ja <= jnatom; ++ja) {
					if( pose.xyz(AtomID(ja,jr)).distance_squared(iX) < 8.0 ) return false;
				}
			}
		}
	}
	return true;
}

Real repackmin(core::pose::Pose & pose, core::scoring::ScoreFunction const & sfxn, Size icys, Size N) {
	core::conformation::ResidueOP old1 = pose.residue(0*N+icys).clone();
	core::conformation::ResidueOP old2 = pose.residue(3*N+icys).clone();
	// core::pose::replace_pose_residue_copying_existing_coordinates(pose,0*N+icys,pose.residue(1).residue_type_set().name_map("CYD"));
	// core::pose::replace_pose_residue_copying_existing_coordinates(pose,3*N+icys,pose.residue(1).residue_type_set().name_map("CYD"));
	pose.conformation().detect_disulfides();
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(pose);
	task->restrict_to_repacking();
	task->or_include_current(true);
	core::pack::pack_rotamers(pose,sfxn,task);

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
	movemap->set_bb(false); movemap->set_chi(true); movemap->set_jump(false);
	// movemap->set_dof();
	protocols::simple_moves::MinMover min(movemap, &sfxn, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
	min.apply(pose);

	Real s = sfxn(pose);
	pose.replace_residue(0*N+icys,*old1,false);
	pose.replace_residue(3*N+icys,*old2,false);
	return s;
}

struct HIT {
	Pose pose;
	Real ddg;
	string fn;
	HIT(Pose const & _pose, Real _ddg, string _fn) : pose(_pose), ddg(_ddg), fn(_fn) {}
};

void generate_disulfide_conformations(core::pose::Pose const & in_pose, Size icys, Size N) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::pose::Pose pose(in_pose);
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	scorefxn->set_weight(core::scoring::dslf_ss_dst,10.0); // S-S
	scorefxn->set_weight(core::scoring::dslf_cs_ang,10.0); // CB-S-S
	scorefxn->set_weight(core::scoring::dslf_ss_dih,10.0); //
	scorefxn->set_weight(core::scoring::dslf_ca_dih, 1.0);
	core::scoring::ScoreFunction sf;

	core::pose::Pose tmp(pose);
	tmp.dump_pdb("test.pdb");
	Real score_ref = repackmin(tmp,*scorefxn,icys,N);
	tmp.dump_pdb("testmin.pdb");


	if( pose.residue(icys).xyz("SG").distance(pose.residue(3*N+icys).xyz("SG")) > 2.2 ) {
		utility_exit_with_message("BAD DISULFIDE!!!");
	}
	Real start_chi1  = dihedral_degrees(pose.xyz(AtomID(1,0*N+icys)),pose.xyz(AtomID(2,0*N+icys)),pose.xyz(AtomID(5,0*N+icys)),pose.xyz(AtomID(6,0*N+icys)));
	Real start_chi2  = dihedral_degrees(pose.xyz(AtomID(2,0*N+icys)),pose.xyz(AtomID(5,0*N+icys)),pose.xyz(AtomID(6,0*N+icys)),pose.xyz(AtomID(6,3*N+icys)));
	// unused ~Labonte
	//Real start_chiss = dihedral_degrees(pose.xyz(AtomID(5,0*N+icys)),pose.xyz(AtomID(6,0*N+icys)),pose.xyz(AtomID(6,3*N+icys)),pose.xyz(AtomID(5,3*N+icys)));
	//Real start_chi2b = dihedral_degrees(pose.xyz(AtomID(2,3*N+icys)),pose.xyz(AtomID(5,3*N+icys)),pose.xyz(AtomID(6,3*N+icys)),pose.xyz(AtomID(5,0*N+icys)));
	//Real start_chi1b = dihedral_degrees(pose.xyz(AtomID(1,3*N+icys)),pose.xyz(AtomID(2,3*N+icys)),pose.xyz(AtomID(5,3*N+icys)),pose.xyz(AtomID(6,3*N+icys)));
	// swap in a CYS with HG placed where paired S goes
	Pose cys;
 	make_pose_from_sequence(cys,"C",core::chemical::FA_STANDARD,false);
	remove_lower_terminus_type_from_pose_residue(cys,1);
	remove_upper_terminus_type_from_pose_residue(cys,1);
	cys.set_dof(core::id::DOF_ID(core::id::AtomID(cys.residue(1).atom_index("HG"),1),core::id::D    ),2.02);
	cys.set_dof(core::id::DOF_ID(core::id::AtomID(cys.residue(1).atom_index("HG"),1),core::id::THETA),1.322473);
	cys.set_chi(1,1,start_chi1);
	cys.set_chi(2,1,start_chi2);
	pose.replace_residue(0*N+icys,cys.residue(1),true);
	pose.replace_residue(1*N+icys,cys.residue(1),true);
	pose.replace_residue(2*N+icys,cys.residue(1),true);
	pose.replace_residue(3*N+icys,cys.residue(1),true);
	pose.replace_residue(4*N+icys,cys.residue(1),true);
	pose.replace_residue(5*N+icys,cys.residue(1),true);

	core::pose::Pose small;
	Size smallsize = 1;
	Size nsmall = 0, smallicys;
	Size smst = (Size)max(1,(int)icys-(int)smallsize);
	for(Size ir = smst; ir <= min(N,icys+smallsize); ++ir) {
		nsmall++;
		if(ir==icys) smallicys = nsmall;
		if(ir==smst || pose.residue(ir).is_lower_terminus()||pose.residue(ir).is_ligand()) {
			   small.append_residue_by_jump(pose.residue(ir),1);
		} else small.append_residue_by_bond(pose.residue(ir));
	}
	for(Size ir = smst; ir <= min(N,icys+smallsize); ++ir) {
		if(ir==smst || pose.residue(ir).is_lower_terminus()||pose.residue(ir).is_ligand()) {
			   small.append_residue_by_jump(pose.residue(3*N+ir),1);
		} else small.append_residue_by_bond(pose.residue(3*N+ir));
	}
	assert( nsmall > 1);
	core::pose::Pose med;
	Size medsize = 20;
	Size nmed = 0, medicys;
	Size medst = (Size)max(1,(int)icys-(int)medsize);
	for(Size ir = medst; ir <= min(N,icys+medsize); ++ir) {
		nmed++;
		if(ir==icys) medicys = nmed;
		if(ir==medst || pose.residue(ir).is_lower_terminus()||pose.residue(ir).is_ligand()) {
			   med.append_residue_by_jump(pose.residue(ir),1);
		} else med.append_residue_by_bond(pose.residue(ir));
	}
	for(Size ir = medst; ir <= min(N,icys+medsize); ++ir) {
		if(ir==medst || pose.residue(ir).is_lower_terminus()||pose.residue(ir).is_ligand()) {
			   med.append_residue_by_jump(pose.residue(3*N+ir),1);
		} else med.append_residue_by_bond(pose.residue(3*N+ir));
	}
	assert( nmed > 21);

	vector1<HIT> alts;
	vector1<vector1<Real> > chis = get_chis(pose,icys);
	for(Size irot = 1; irot <= chis.size(); ++irot){
		Real chi1i = chis[irot][1];
		Real chi2i = chis[irot][2];
		Real chi1j = chis[irot][3];
		Real chi2j = chis[irot][4];
		std::cout << "progress " << 2*irot << " of " << 2*chis.size() << std::endl;
		for(Size idsf = 0; idsf < 2; ++idsf) {
			Real chiss = 90.0;
			if( idsf == 1 ) chiss =  -90.0;
			// if( idsf == 2 ) chiss =   80.0;
			// if( idsf == 3 ) chiss =  100.0;
			// if( idsf == 4 ) chiss =  -80.0;
			// if( idsf == 5 ) chiss = -100.0;

			string fn = "alt_disulf_"+string_of(chi1i)+"_"+string_of(chi2i)+"_"+string_of(chiss)+"_"+string_of(chi2j)+"_"+string_of(chi1j)+".pdb";
			// std::cout << "check " << string_of(chi1i)+"_"+string_of(chi2i)+"_"+string_of(chiss)+"_"+string_of(chi2j)+"_"+string_of(chi1j) << std::endl;
			// check small
			set_disulf(small,smallicys,nsmall,1,chi1i,chi2i,chiss,chi2j,chi1j);
			if( !clash_check_bb(small,smallicys,nsmall,1) ) {
				std::cout << "clash small  " << chi1i << " " << " " << chi2i << " " << chiss << " " << chi2j << " " << chi1j << std::endl;
				continue;
			}
			set_disulf(med,medicys,nmed,1,chi1i,chi2i,chiss,chi2j,chi1j);
			if( !clash_check_bb(med,medicys,nmed,1) ) {
				std::cout << "clash medium " << chi1i << " " << " " << chi2i << " " << chiss << " " << chi2j << " " << chi1j << std::endl;
				continue;
			}
			// now full
			set_disulf(pose,icys,N,3,chi1i,chi2i,chiss,chi2j,chi1j);
			if( !clash_check_bb(pose,icys,N,3) ) {
				std::cout << "clash full   " << chi1i << " " << " " << chi2i << " " << chiss << " " << chi2j << " " << chi1j << std::endl;
				continue;
			}
			if( core::scoring::CA_rmsd(in_pose,pose) < option[gendsf::rms_cutoff]() ) {
				std::cout << "same " << " " << string_of(chi1i)+"_"+string_of(chi2i)+"_"+string_of(chiss)+"_"+string_of(chi2j)+"_"+string_of(chi1j)+".pdb";
				// pose.dump_pdb("same_"+string_of(chi1i)+"_"+string_of(chi2i)+"_"+string_of(chiss)+"_"+string_of(chi2j)+"_"+string_of(chi1j)+".pdb");
				continue;
			}
			Real ddg = repackmin(pose,*scorefxn,icys,N) - score_ref;
			// if(ddg > 100.0) continue;
			bool seenit = false;
			for(Size ip = 1; ip <= alts.size(); ++ip) {
				if( core::scoring::CA_rmsd(pose,alts[ip].pose) < option[gendsf::rms_cutoff]() ) {
					if( alts[ip].ddg > ddg ) {
						// std::cout << "alt_disulfide ddG " << ddg << " " << fn << std::endl;
						alts[ip] = HIT(pose,ddg,fn); // if lower ddg, replace it
					}
					seenit = true;
				}
			}
			if(!seenit) {
				alts.push_back(HIT(pose,ddg,fn));
			}
			// utility_exit_with_message("arst");
		}
	}
	for(Size ihit = 1; ihit <= alts.size(); ++ihit) {
		std::cout << "alt_disulfide ddG " << alts[ihit].ddg << " " << alts[ihit].fn << std::endl;
		alts[ihit].pose.dump_pdb(option[out::file::o]+"/"+alts[ihit].fn);
	}
}


int main (int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;
	for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
		string fn = option[in::file::s]()[ifn];
		Pose pnat;
		core::import_pose::pose_from_file(pnat,fn, core::import_pose::PDB_file);
		Size N = get_N(pnat);
		Size icys = 1; for(icys = 1; icys < N; ++icys) if(pnat.residue(icys).aa() == core::chemical::aa_cys) break;
		Pose tri2 = make_two_trimers(pnat,icys,N);
		generate_disulfide_conformations(tri2,icys,N);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


