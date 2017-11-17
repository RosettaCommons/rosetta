// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/willmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <basic/database/open.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/will_util.ihh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


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
using core::import_pose::pose_from_file;
using basic::options::option;
using numeric::min;
using numeric::max;

typedef utility::vector1<core::Real> Reals;
typedef utility::vector1<core::Size> Sizes;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

static basic::Tracer TR( "willmatch_d6_bpy" );


void myoptH(Pose & pose, ScoreFunctionOP sf) {
	add_lower_terminus_type_to_pose_residue(pose,1);
	add_upper_terminus_type_to_pose_residue(pose,pose.size());
	core::pack::optimizeH(pose,*sf);
	remove_lower_terminus_type_from_pose_residue(pose,1);
	remove_upper_terminus_type_from_pose_residue(pose,pose.size());
}


void refine(Pose & pose, ScoreFunctionOP sf, Size r1, Size r2, Size r3, Size r4 ) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	Size N = pose.size();
	vector1< bool > aac(20,false);
	aac[core::chemical::aa_ala] = true;
	//aac[core::chemical::aa_cys] = true;
	//aac[core::chemical::aa_asp] = true;
	//aac[core::chemical::aa_glu] = true;
	aac[core::chemical::aa_phe] = true;
	//aac[core::chemical::aa_gly] = true;
	//aac[core::chemical::aa_his] = true;
	aac[core::chemical::aa_ile] = true;
	//aac[core::chemical::aa_lys] = true;
	aac[core::chemical::aa_leu] = true;
	aac[core::chemical::aa_met] = true;
	aac[core::chemical::aa_asn] = true;
	//aac[core::chemical::aa_pro] = true;
	aac[core::chemical::aa_gln] = true;
	//aac[core::chemical::aa_arg] = true;
	aac[core::chemical::aa_ser] = true;
	aac[core::chemical::aa_thr] = true;
	aac[core::chemical::aa_val] = true;
	aac[core::chemical::aa_trp] = true;
	aac[core::chemical::aa_tyr] = true;

	vector1<Size> iface(N,1);
	for ( Size ir = 1; ir <= N; ++ir ) {
		if ( ir==r1 || ir==r2 || ir==r3 || ir==r4 ) iface[ir] = 0;
		else if ( pose.residue(ir).name3()=="CYS" || pose.residue(ir).name3()=="GLY" || pose.residue(ir).name3()=="PRO" ) iface[ir] = 0;
		else {
			Real closestatom=9e9,closestcb=9e9;
			if ( pose.residue(r1).xyz("NE1").distance_squared( pose.xyz(AtomID(2,ir)) ) > 225.0 ) continue;
			//Real closestatom=9e9,closestcb=9e9;
			for ( Size jri = 1; jri <= 4; jri++ ) {
				Size jr = r1; if ( 2==jri ) jr = r2; else  if ( 3==jri ) jr = r3; else  if ( 4==jri ) jr = r4;
				for ( Size ja = 5; ja <= pose.residue(jr).nheavyatoms(); ++ja ) {
					Vec aj = pose.xyz(AtomID(ja,jr));
					for ( Size ia = 5; ia <= pose.residue(ir).nheavyatoms(); ++ia ) {
						Real d = aj.distance_squared( pose.xyz(AtomID(ia,ir)) );
						closestatom = min(closestatom,d);
						if ( ia==5 ) closestcb = min(closestcb,d);
					}
				}
			}
			closestcb = sqrt(closestcb);
			closestatom = sqrt(closestatom);
			if ( closestatom < 6.0 ) {
				iface[ir] = 3;
			}//  else if( closestatom < 8.0) {
			//  iface[ir] = 2;
			// }
		}
	}

	for ( Size i = 1; i <= N; ++i ) {
		if (        iface[i] == 3 ) {
			bool tmp = aac[pose.residue(i).aa()];
			aac[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aac);
			task->nonconst_residue_task(i).or_include_current(true);
			aac[pose.residue(i).aa()] = tmp;
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else if ( iface[i] == 2 ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_include_current(true);
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else {
			task->nonconst_residue_task(i).prevent_repacking();
		}

	}

	Real worig = sf->get_weight(core::scoring::res_type_constraint);
	if ( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);
	pose.add_constraints( res_cst );

	protocols::simple_moves::PackRotamersMover repack( sf, task );
	repack.apply(pose);

	// cleanup 2
	pose.remove_constraints( res_cst );
	sf->set_weight(core::scoring::res_type_constraint,worig);


	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	// movemap->set_jump(false);
	// movemap->set_bb(false);
	// movemap->set_chi(true);
	// protocols::simple_moves::MinMover( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);
	// movemap->set_jump(false);
	// movemap->set_bb(true);
	// movemap->set_chi(true);
	// protocols::simple_moves::MinMover( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);
	movemap->set_jump(true);
	movemap->set_bb(true);
	movemap->set_chi(true);
	protocols::simple_moves::MinMover( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);


}


// find 2 HIS that can chelate a tetrahedral metal
void run() {
	using namespace basic::options::OptionKeys;
	using namespace core::id;
	using namespace core;
	using core::scoring::func::FuncOP;

	core::chemical::ResidueTypeSetCOP cen_residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ) );
	core::chemical::ResidueTypeSetCOP  fa_residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	Real tmpdis = option[willmatch::max_dis_metal]();
	Real const MXDSMTL = tmpdis*tmpdis;
	Real const MXAGMTL = option[willmatch::max_ang_metal]();
	Real const MATCH_OVERLAP_DOT = cos(numeric::conversions::radians(option[willmatch::match_overlap_ang]()));

	core::io::silent::SilentFileData sfd;

	vector1<string> infiles;
	if ( option[in::file::l].user() ) {
		utility::io::izstream in(option[in::file::l]()[1]);
		string tmp;
		while ( in >> tmp ) infiles.push_back(tmp);
	} else if ( option[in::file::s].user() ) {
		infiles = option[in::file::s]();
	} else {
		utility_exit_with_message("no input!");
	}

	std::string startfile = "";
	Size startbpy = 0;
	{
		utility::io::izstream inz(option[out::file::o]()+"/willmatch_d6_bpy.progress");
		while ( inz >> startfile >> startbpy ) ;
		inz.close();
		if ( startfile != "" ) TR << "continuing from " << startfile << " " << startbpy << std::endl;
	}
	utility::io::ozstream oprogress(option[out::file::o]()+"/willmatch_d6_bpy.progress");

	for ( Size ifile = 1; ifile <= infiles.size(); ifile++ ) {
		string infile = infiles[ifile];
		if ( startfile != "" && startfile != infile ) continue; // CHECKPOINT
		Pose in_fa;
		pose_from_file(in_fa, *fa_residue_set,infile, core::import_pose::PDB_file);
		if ( in_fa.size() > 300 ) continue;
		for ( Size ir = 1; ir <= in_fa.size(); ++ir ) {
			if ( in_fa.residue(ir).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(in_fa,ir);
			if ( in_fa.residue(ir).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(in_fa,ir);
		}
		Pose native = in_fa;
		Size nres = in_fa.size();
		core::scoring::packing::HolesParams hp_dec15;
		hp_dec15.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
		Real natrholes = core::scoring::packing::compute_dec15_score(native);

		core::chemical::ResidueType const & rtala( in_fa.residue(1).residue_type_set().name_map("ALA") );
		//core::chemical::ResidueType const & rthis( in_fa.residue(1).residue_type_set().name_map("HIS") );
		//core::chemical::ResidueType const & rtglu( in_fa.residue(1).residue_type_set().name_map("GLU") );
		//core::chemical::ResidueType const & rtbpy( in_fa.residue(1).residue_type_set().name_map("BPY") );
		core::chemical::ResidueType const & rtphe( in_fa.residue(1).residue_type_set().name_map("PHE") );
		for ( Size i = 1; i <= nres; ++i ) core::pose::replace_pose_residue_copying_existing_coordinates(in_fa,i,rtala);
		Pose const fa_pose = in_fa;
		ImplicitFastClashCheck clashcheck(fa_pose,basic::options::option[basic::options::OptionKeys::willmatch::clash_dis]());

		ScoreFunctionOP sf     = core::scoring::get_score_function();
		ScoreFunctionOP sfhard = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
		//sf->set_weight(core::scoring::fa_dun,1.0);
		sf->set_weight(core::scoring::fa_sol,0.0);
		sf->set_weight(core::scoring::hbond_sr_bb,2.0);
		sf->set_weight(core::scoring::hbond_lr_bb,2.0);
		sf->set_weight(core::scoring::hbond_sc,3.0);
		sf->set_weight(core::scoring::hbond_bb_sc,3.0);
		sf->set_weight(core::scoring::fa_intra_rep,0.4);
		sf->set_weight(core::scoring::fa_dun,0.1);
		sf->set_weight(core::scoring::atom_pair_constraint,1.0);
		sf->set_weight(core::scoring::angle_constraint    ,1.0);
		sf->set_weight(core::scoring::dihedral_constraint ,1.0);
		//sfhard->set_weight(core::scoring::fa_atr,1.2);
		//sfhard->set_weight(core::scoring::fa_rep,0.2);
		//sfhard->set_weight(core::scoring::fa_sol,0.0);
		sfhard->set_weight(core::scoring::hbond_sr_bb,3.0);
		sfhard->set_weight(core::scoring::hbond_lr_bb,3.0);
		//sfhard->set_weight(core::scoring::hbond_sc   ,3.0);
		//sfhard->set_weight(core::scoring::hbond_bb_sc,3.0);
		sfhard->set_weight(core::scoring::fa_intra_rep,0.4);
		sfhard->set_weight(core::scoring::fa_dun,0.1);
		sfhard->set_weight(core::scoring::atom_pair_constraint,1.0);
		sfhard->set_weight(core::scoring::angle_constraint    ,1.0);
		sfhard->set_weight(core::scoring::dihedral_constraint ,1.0);

		Real chi1incr = option[willmatch::chi1_increment]();
		Real chi2incr = option[willmatch::chi2_increment]();
		vector1<Real> CHI1,CHI2;
		Real ofst1 = floor(chi1incr/2.0)+0.5, ofst2 = floor(chi2incr/2.0)+0.5;

		for ( Real i = ofst1; i < 360.0; i+= chi1incr ) CHI1.push_back( (i>180.0 ? i-360.0 : i));
		for ( Real i = ofst2; i < 360.0; i+= chi2incr ) CHI2.push_back( (i>180.0 ? i-360.0 : i));

		//utility::io::ozstream bhhout(option[out::file::o]()+"/"+utility::file_basename(infile)+".bhh_match");
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
		for ( Size i = 1; i <= CHI1.size(); ++i ) {
			hsd.set_chi(1,1,CHI1[i]);
			hse.set_chi(1,1,CHI1[i]);
			bpy.set_chi(1,1,CHI1[i]);
			glu.set_chi(1,1,CHI1[i]);
			hci1cen(1,i) = stub.global2local( hsd.residue(1).xyz("CG")-0.007*((hsd.residue(1).xyz("CG")-hsd.residue(1).xyz("CB")).normalized()) ); // CG - 0.006699*(CG-CB) 3.11925602854 from actual cen
			hci1cen(2,i) = stub.global2local( hse.residue(1).xyz("CG")+3.947*((hse.residue(1).xyz("CG")-hse.residue(1).xyz("CB")).normalized()) ); // CG + 3.946941*(CG-CB) 1.53906602847 from actual cen
			bci1cen(  i) = stub.global2local( bpy.residue(1).xyz("CG")+3.034*((bpy.residue(1).xyz("CG")-bpy.residue(1).xyz("CB")).normalized()) ); // GC + 3.033686*(CG-CB) 3.23823647147 from actual cen
			eci1cen(  i) = stub.global2local( glu.residue(1).xyz("CG")+1.991*((glu.residue(1).xyz("CG")-glu.residue(1).xyz("CB")).normalized()) ); // CG + 1.991210*(CG-CB) 3.86805118484 from actual cen
			for ( Size j = 1; j <= CHI2.size(); ++j ) {
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

		ObjexxFCL::FArray1D<Real> bsasa0(            in_fa.size(),9e9);
		ObjexxFCL::FArray2D<Real> bsasa1(CHI1.size(),in_fa.size(),9e9);

		vector1<Real> natsasa; {
			core::id::AtomID_Map<Real> atom_sasa;
			core::scoring::calc_per_atom_sasa( native, atom_sasa, natsasa, 2.0, false );   }
		//for(Size i = 1; i <= natsasa.size(); ++i) std::cout << i << " " << natsasa[i] << std::endl;

		vector1<Size> scanres;
		// scanres.push_back(4);
		// scanres.push_back(18);
		// scanres.push_back(21);
		// scanres.push_back(88);
		//std::cerr << "aoirstenarsoit" << std::endl;
		if ( option[willmatch::residues].user() ) {
			TR << "input scanres!!!!!!" << std::endl;
			scanres = option[willmatch::residues]();
		} else {
			for ( Size i = 1; i <= in_fa.size(); ++i ) {
				if ( !in_fa.residue(i).has("N" ) ) { continue; }
				if ( !in_fa.residue(i).has("CA") ) { continue; }
				if ( !in_fa.residue(i).has("C" ) ) { continue; }
				if ( !in_fa.residue(i).has("O" ) ) { continue; }
				if ( !in_fa.residue(i).has("CB") ) { continue; }
				if ( in_fa.residue(i).name3()=="PRO" ) { continue; }
				//if(natsasa[i] > 0) continue;
				scanres.push_back(i);
			}
		}


		// precompute acceptable chis
		ObjexxFCL::FArray1D<Real> ballow0(                        in_fa.size(),9e9),eallow0(                        in_fa.size(),9e9),hallow0(                        in_fa.size(),9e9);
		ObjexxFCL::FArray2D<Real> ballow1(            CHI1.size(),in_fa.size(),9e9),eallow1(            CHI1.size(),in_fa.size(),9e9),hallow1(            CHI1.size(),in_fa.size(),9e9);
		ObjexxFCL::FArray3D<Real> ballow2(CHI2.size(),CHI1.size(),in_fa.size(),9e9),eallow2(CHI2.size(),CHI1.size(),in_fa.size(),9e9),hallow2(CHI2.size(),CHI1.size(),in_fa.size(),9e9);
		//core::pack::rotamers::SingleResidueRotamerLibraryCAP plib = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( rtphe );
		//core::pack::rotamers::SingleResidueRotamerLibraryCAP hlib = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( hsd.residue(1).type() );
		//core::pack::rotamers::SingleResidueRotamerLibraryCAP elib = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( glu.residue(1).type() );
		//core::pack::dunbrack::RotamerLibraryScratchSpace scratch;
		{
			Real mbpyir = 9e9;
			TR << "precomputing dun scores and clashes, size " << scanres.size() << std::endl;
			// BPY
			Pose pose = in_fa,pose3=bpy;
			pose3.append_residue_by_bond(pose.residue(2));
			pose3.append_residue_by_bond(pose.residue(3));
			Size count = 0;
			TR << scanres.size() << " ";
			for ( Sizes::const_iterator iri = scanres.begin(); iri != scanres.end(); ++iri ) {
				TR << ++count << " ";
				Size const ir(*iri);
				if ( natsasa[ir] > 0.0 ) { ballow0(ir)=9e9; } else {
					core::pose::replace_pose_residue_copying_existing_coordinates(pose,ir,rtphe); // clash w/BPY
					pose3.replace_residue(1,pose.residue(((ir-1+pose.size()-1)%(pose.size()))+1),false);
					pose3.replace_residue(2,pose.residue(((ir+0+pose.size()-1)%(pose.size()))+1),false);
					pose3.replace_residue(3,pose.residue(((ir+1+pose.size()-1)%(pose.size()))+1),false);
					for ( Size i1 = 1; i1 <= CHI1.size(); ++i1 ) {
						pose3.set_chi(1,2,CHI1[i1]);
						for ( Size i2 = 1; i2 <= CHI2.size(); ++i2 ) {
							pose3.set_chi(2,2,CHI2[i2]);
							sf->score(pose3);
							ballow2(i2,i1,ir) = pose3.energies().residue_total_energies(2)[core::scoring::fa_intra_rep] - 7.0;
							mbpyir = min(mbpyir,ballow2(i2,i1,ir));
						}
					}
					pose.replace_residue(ir,bpy.residue(1),true);
					pose3.replace_residue(2,pose.residue(ir),false);
					for ( Size i1 = 1; i1 <= CHI1.size(); ++i1 ) {
						pose3.set_chi(1,2,CHI1[i1]);
						for ( Size i2 = 1; i2 <= CHI2.size(); ++i2 ) {
							pose3.set_chi(2,2,CHI2[i2]);
							for ( Size a = 6; a <= pose3.residue(2).nheavyatoms(); ++a ) {
								if ( !clashcheck.clash_check(pose3.residue(2).xyz(a),ir) ) {
									ballow2(i2,i1,ir) += 9e9;
									break;
								}
							}
							//if( bsasa1(i1,ir) > 0.1 ) ballow2(i2,i1,ir) += 9e9; // this doesn't seem to work...
							if ( ballow1(i1,ir) > ballow2(i2,i1,ir) ) ballow1(i1,ir) = ballow2(i2,i1,ir);
						}
						if ( ballow0(ir) > ballow1(i1,ir) ) ballow0(ir) = ballow1(i1,ir);
					}
				}
				// // GLU
				// pose.replace_residue(ir,glu.residue(1),true);
				// pose3.replace_residue(2,pose.residue(ir),false);
				// pose3.set_phi(2,pose.phi(ir)); pose3.set_psi(2,pose.psi(ir)); pose3.set_omega(2,pose.omega(ir));
				// for(Size i1 = 1; i1 <= CHI1.size(); ++i1) {
				//  pose3.set_chi(1,2,CHI1[i1]);
				//  for(Size i2 = 1; i2 <= CHI2.size(); ++i2) {
				//    pose3.set_chi(2,2,CHI2[i2]);
				//    Real dbmin = 9e9, argmn = 0.0;
				//    for(Real i3 = -177.5; i3 < 180.0; i3+=1) {
				//      pose3.set_chi(3,2,i3);
				//      sf->score(pose3);
				//      Real tmp = pose3.energies().residue_total_energies(2)[core::scoring::fa_intra_rep];// elib->rotamer_energy( pose3.residue(2), scratch );
				//      if( tmp < dbmin ) { dbmin = tmp; argmn = i3; }
				//    }
				//    //std::cerr << "TEST " << dbmin << " " << ir << " " << i1 << " " << i2 << " " << argmn << std::endl;
				//    eallow2(i2,i1,ir) = dbmin;
				//    for(Size a = 6; a <= pose3.residue(2).nheavyatoms()-2; ++a) {
				//      if(!clashcheck.clash_check(pose3.residue(2).xyz(a),ir)) {
				//        eallow2(i2,i1,ir) += 9e9;
				//        break;
				//      }
				//    }
				//    if( eallow1(i1,ir) > eallow2(i2,i1,ir) ) eallow1(i1,ir) = eallow2(i2,i1,ir);
				//  }
				//  if( eallow0(ir) > eallow1(i1,ir) ) eallow0(ir) = eallow1(i1,ir);
				// }
				// HIS
				pose.replace_residue(ir,hsd.residue(1),true);
				pose3.replace_residue(2,pose.residue(ir),false);
				pose.replace_residue(ir,hse.residue(1),true);
				pose3.set_phi  (2,pose.phi  (ir)); pose3.set_psi  (2,pose.psi  (ir)); pose3.set_omega(2,pose.omega(ir));
				for ( Size i1 = 1; i1 <= CHI1.size(); ++i1 ) {
					pose3.set_chi(1,2,CHI1[i1]);
					for ( Size i2 = 1; i2 <= CHI2.size(); ++i2 ) {
						pose3.set_chi(2,2,CHI2[i2]);
						sf->score(pose3);
						hallow2(i2,i1,ir) = pose3.energies().residue_total_energies(2)[core::scoring::fa_intra_rep];//hlib->rotamer_energy( pose3.residue(2), scratch );
						for ( Size a = 6; a <= pose3.residue(2).nheavyatoms(); ++a ) {
							if ( !clashcheck.clash_check(pose3.residue(2).xyz(a),ir) ) {
								hallow2(i2,i1,ir) += 9e9;
								break;
							}
						}
						if ( hallow1(i1,ir) > hallow2(i2,i1,ir) ) hallow1(i1,ir) = hallow2(i2,i1,ir);
					}
					if ( hallow0(ir) > hallow1(i1,ir) ) hallow0(ir) = hallow1(i1,ir);
				}
				core::pose::replace_pose_residue_copying_existing_coordinates(pose,ir,rtala);
			}
			TR << "DONE precomputing dun scores and clashes " << mbpyir << std::endl;
		}
		Real const DUN_THRESH = option[willmatch::fa_dun_thresh]();
		//utility_exit_with_message("aorsnt");
		// string fname = utility::file_basename(infile)+"_d6_bpy_matches.pdb";
		// utility::io::ozstream out(fname);
		// out << "MODEL BASE"  << endl;
		// fa_pose.dump_pdb(out);
		// out << "ENDMDL" << endl;

		Size count = 0;

		vector1<Size> batm;
		TR << "BATM:";
		for ( Size i = 6; i <= bpy.residue(1).nheavyatoms(); ++i ) {
			// TR << "'" << bpy.residue(1).atom_name(i) << "'" << std::endl;
			if ( bpy.residue(1).atom_name(i)==" NE1" ) { continue; }
			if ( bpy.residue(1).atom_name(i)==" NN1" ) { continue; }
			if ( bpy.residue(1).atom_name(i)==" ZN " ) { continue; }
			batm.push_back(i);
			TR << " " << i;
		}
		TR << std::endl;
		vector1<Size> hatm;
		TR << "HATM:";
		for ( Size i = 6; i <= hsd.residue(1).nheavyatoms(); ++i ) {
			if ( hsd.residue(1).atom_name(i)==" ND1" ) { continue; }
			if ( hsd.residue(1).atom_name(i)==" NE2" ) { continue; }
			hatm.push_back(i);
			TR << " " << i;
		}
		TR << std::endl;
		vector1<Size> eatm;
		TR << "EATM:";
		for ( Size i = 6; i <= glu.residue(1).nheavyatoms(); ++i ) {
			if ( glu.residue(1).atom_name(i)==" ND1" ) { continue; }
			if ( glu.residue(1).atom_name(i)==" NE2" ) { continue; }
			eatm.push_back(i);
			TR << " " << i;
		}
		TR << std::endl;

		vector1<Vec> Ns(in_fa.size()),CAs(in_fa.size()),Cs(in_fa.size()),CBs(in_fa.size());
		for ( Size i = 1; i <= in_fa.size(); ++i ) Ns [i] = in_fa.residue(i).xyz("N" );
		for ( Size i = 1; i <= in_fa.size(); ++i ) CAs[i] = in_fa.residue(i).xyz("CA");
		for ( Size i = 1; i <= in_fa.size(); ++i ) Cs [i] = in_fa.residue(i).xyz("C" );
		for ( Size i = 1; i <= in_fa.size(); ++i ) CBs[i] = in_fa.residue(i).xyz("CB");

		Pose pose = fa_pose;
		//    Size lastb=0,lasti=0,lastj=0,laste=0;
		//for(vector1<Size>::const_iterator biter = scanres.begin(); biter != scanres.end(); ++biter) {
		TR << "matching res, size " << scanres.size() << std::endl;
		//#pragma omp parallel for
		for ( Size ii = 1; ii <= scanres.size(); ++ii ) {
			if ( count%10==0 ) TR << ObjexxFCL::format::I(4,ii) << "/" << scanres.size() << std::endl;
			vector1<Vec> foundcenbg,foundori1bg,foundori2bg,foundori3bg,foundori4bg,foundori5bg,foundori6bg;
			vector1<Vec> foundcensg,foundori1sg,foundori2sg,foundori3sg,foundori4sg,foundori5sg,foundori6sg;
			//Size brsd = *biter;
			Size brsd = scanres[ ii ];
			//      if(brsd != 110) continue;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if ( natsasa[brsd] > 0 ) continue;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if ( startbpy != 0 && startbpy != brsd ) continue;        // CHECKPOINT
			if ( startbpy == brsd ) { startfile = ""; startbpy = 0; } // CHECKPOINT
			oprogress << infile << " " << brsd << std::endl;
			//TR << "scanning bpy rsd " << infile << " " << brsd << std::endl;
			if ( natsasa[brsd] > 0.1 ) continue;
			Stub s3(CBs[brsd],CAs[brsd],CAs[brsd]+(Ns[brsd]-Cs[brsd]));
			if ( ballow0(brsd) > DUN_THRESH ) continue; // if all chi1 / chi2 are clash or bad rot
			for ( Size kch1 = 1; kch1 <= CHI1.size(); ++kch1 ) {
				if ( ballow1(kch1,brsd) > DUN_THRESH ) continue; // if all chi2 are clash or bad rot
				for ( Size kch2 = 1; kch2 <= CHI2.size(); ++kch2 ) {
					if ( ballow2(kch2,kch1,brsd) > DUN_THRESH ) continue;  // chi2 is clash or bad rot
					if ( !clashcheck.clash_check(s3.local2global(bci2end(kch1,kch2))) ) { continue; }
					if ( !clashcheck.clash_check(s3.local2global(bci2cen(kch1,kch2))) ) { continue; }
					pose.replace_residue(brsd,bpy.residue(1),true);
					pose.set_chi(1,brsd,CHI1[kch1]);
					pose.set_chi(2,brsd,CHI2[kch2]);

					bool clash = false;
					// already done above
					for ( Size k = 6; k <= pose.residue(brsd).nheavyatoms(); ++k ) {
						if ( !clashcheck.clash_check(pose.residue(brsd).xyz(k),brsd) ) clash=true;
						if ( clash ) break;
					}
					if ( clash ) { continue; }

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
					for ( vector1<Size>::const_iterator iiter = scanres.begin(); iiter != scanres.end(); ++iiter ) {
						Size irsd = *iiter;
						if ( irsd==brsd ) { continue; }
						Stub si(CBs[irsd],CAs[irsd],CAs[irsd]+(Ns[irsd]-Cs[irsd]));
						/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						for ( Size ide = 1; ide <= 2; ide++ ) {
							/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							if ( hallow0(irsd) > DUN_THRESH ) continue; // if all chi1 / chi2 are clash or bad rot
							if ( CBs[irsd].distance_squared(cenb) > ((ide==1)?hddcbmx2:hedcbmx2) ) { continue; }              // if CB too far, skip
							for ( Size ich1 = 1; ich1 <= CHI1.size(); ++ich1 ) {
								if ( hallow1(ich1,irsd) > DUN_THRESH ) continue; // if all chi2 are clash or bad rot
								if ( si.local2global(hci1cen(ide,ich1)).distance_squared(cenb) > ((ide==1)?hddc1mx2:hedc1mx2) ) { continue; }              // if chi1 cen too far, skip
								for ( Size ich2 = 1; ich2 <= CHI2.size(); ++ich2 ) {
									if ( hallow2(ich2,ich1,irsd) > DUN_THRESH ) continue;  // chi2 is clash or bad rot

									Vec const ceni = si.local2global(chi2cen(ide,ich1,ich2));
									if ( ceni.distance_squared(cenb) > MXDSMTL ) { continue; }
									Vec const orii = (si.local2global(chi2ori(ide,ich1,ich2))-ceni).normalized();
									bool his180 = false;
									{
										Real ang1 = numeric::conversions::degrees(acos(orii.dot(orik1)));
										Real ang2 = numeric::conversions::degrees(acos(orii.dot(orik2)));
										if ( ang1 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL) ) { continue; }
										if ( ang2 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang2 && ang2 < 180-MXAGMTL) ) { continue; }
										his180 = (ang1 > 135.0 || ang2 > 135.0);
									}
									if ( !clashcheck.clash_check(si.local2global(chi2cen(ide,ich1,ich2))) ) { continue; }

									pose.replace_residue(irsd,( ide==1 ? hsd : hse ).residue(1),true);
									pose.set_chi(1,irsd,CHI1[ich1]);
									pose.set_chi(2,irsd,CHI2[ich2]);
									Vec const ligi = pose.residue(irsd).xyz( (ide==1)?"ND1":"NE2" );

									clash = false;
									for ( vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b ) {
										for ( vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h ) {
											if ( pose.residue(brsd).xyz(*b).distance_squared(pose.residue(irsd).xyz(*h)) < 7.0 ) clash=true;
										}
									}
									if ( clash ) { continue; }

									for ( Size a = 6; a <= pose.residue(irsd).nheavyatoms(); ++a ) {
										if ( !clashcheck.clash_check(pose.residue(irsd).xyz(a),irsd) ) { clash=true; if ( clash ) break; }
									}
									if ( clash ) { continue; }

									for ( vector1<Size>::const_iterator jiter = scanres.begin(); jiter != scanres.end(); ++jiter ) {
										Size jrsd = *jiter;
										if ( jrsd <= irsd ) { continue; }
										if ( jrsd==brsd || jrsd==irsd ) { continue; }
										Stub sj(CBs[jrsd],CAs[jrsd],CAs[jrsd]+(Ns[jrsd]-Cs[jrsd]));
										/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										for ( Size jde = 1; jde <= 2; jde++ ) {
											/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
											if ( hallow0(jrsd) > DUN_THRESH ) continue; // if all chi1 / chi2 are clash or bad rot
											if ( CBs[jrsd].distance_squared(cenb) > ((jde==1)?hddcbmx2:hedcbmx2) ) { continue; }  // if CB too far, skip
											for ( Size jch1 = 1; jch1 <= CHI1.size(); ++jch1 ) {
												if ( hallow1(jch1,jrsd) > DUN_THRESH ) continue; // if all chi2 are clash or bad rot
												if ( sj.local2global(hci1cen(jde,jch1)).distance_squared(cenb) > ((jde==1)?hddc1mx2:hedc1mx2) ) { continue; }              // if chi1 cen too far, skip
												for ( Size jch2 = 1; jch2 <= CHI2.size(); ++jch2 ) {
													if ( hallow2(jch2,jch1,jrsd) > DUN_THRESH ) continue;  // chi2 is clash or bad rot

													Vec const cenj = sj.local2global(chi2cen(jde,jch1,jch2));
													if ( cenj.distance_squared(cenb) > MXDSMTL ) { continue; }
													Vec const orij = (sj.local2global(chi2ori(jde,jch1,jch2))-cenj).normalized();
													bool localhis180 = his180;
													bool hishis180 = false;
													{
														Real ang1 = numeric::conversions::degrees(acos(orij.dot(orii)));
														if ( ang1 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL) ) { continue; }
														if ( ang1 > 135.0 ) hishis180 = true;
													}
													if ( !localhis180 && !hishis180 ) {
														// TR << "NOT his180, hishis180" << std::endl;
														Real ang1 = numeric::conversions::degrees(acos(orij.dot(orik1)));
														Real ang2 = numeric::conversions::degrees(acos(orij.dot(orik2)));
														if ( ang1 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang1 && ang1 < 180-MXAGMTL) ) { continue; }
														if ( ang2 < 90.0-MXAGMTL || (90.0+MXAGMTL < ang2 && ang2 < 180-MXAGMTL) ) { continue; }
														localhis180 = (ang1 > 135.0 || ang2 > 135.0);
													} else {
														// TR << "his180 OR hishis180" << std::endl;
														Real ang1 = numeric::conversions::degrees(acos(orij.dot(orik1)));
														Real ang2 = numeric::conversions::degrees(acos(orij.dot(orik2)));
														if ( ang1 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang1 ) { continue; }
														if ( ang2 < 90.0-MXAGMTL || 90.0+MXAGMTL < ang2 ) { continue; }
													}
													if ( !clashcheck.clash_check(sj.local2global(chi2cen(jde,jch1,jch2))) ) { continue; }

													pose.replace_residue(jrsd,( jde==1 ? hsd : hse ).residue(1),true);
													pose.set_chi(1,jrsd,CHI1[jch1]);
													pose.set_chi(2,jrsd,CHI2[jch2]);
													Vec const ligj = pose.residue(jrsd).xyz( (jde==1)?"ND1":"NE2" );

													clash = false;
													for ( vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h ) {
														for ( vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b ) {
															if ( pose.residue(brsd).xyz(*b).distance_squared(pose.residue(jrsd).xyz(*h)) < 7.0 ) {
																clash=true;
															}
														}
														for ( vector1<Size>::const_iterator h2 = hatm.begin(); h2 != hatm.end(); ++h2 ) {
															if ( pose.residue(irsd).xyz(*h2).distance_squared(pose.residue(jrsd).xyz(*h)) < 7.0 ) {
																clash=true;
															}
														}
													}
													if ( clash ) { continue; }

													for ( Size a = 6; a <= pose.residue(jrsd).nheavyatoms(); ++a ) {
														if ( !clashcheck.clash_check(pose.residue(jrsd).xyz(a),jrsd) ) { clash=true; if ( clash ) break; }
													}
													if ( clash ) { continue; }

													//bhhout<< brsd << " " << kch1 << " " << kch2 <<
													//  " " << irsd << " " << ich1 << " " << ich2 << " " << ide <<
													//  " " << jrsd << " " << jch1 << " " << jch2 << " " << jde << std::endl;

													for ( vector1<Size>::const_iterator eiter = scanres.begin(); eiter != scanres.end(); ++eiter ) {
														Size ersd = *eiter;
														if ( ersd==brsd||ersd==irsd||ersd==jrsd ) { continue; }

														clash = false;
														Stub se(CBs[ersd],CAs[ersd],CAs[ersd]+(Ns[ersd]-Cs[ersd]));

														//if(eallow0(ersd) > DUN_THRESH) continue; // if all chi1 / chi2 are clash or bad rot
														if ( CBs[ersd].distance_squared(cenb) > edscbmx2 ) { continue; }              // if CB too far, skip
														for ( Size ech1 = 1; ech1 <= CHI1.size(); ++ech1 ) {
															//if(eallow1(ech1,ersd) > DUN_THRESH) continue; // if all chi2 are clash or bad rot
															if ( se.local2global(eci1cen(ech1)).distance_squared(cenb) > edsc1mx2 ) { continue; }              // if chi1 cen too far, skip
															for ( Size ech2 = 1; ech2 <= CHI2.size(); ++ech2 ) {
																//if(eallow2(ech2,ech1,ersd) > DUN_THRESH) continue;  // chi2 is clash or bad rot

																Vec const cene = se.local2global(gci2cen(ech1,ech2));
																if ( cene.distance_squared(cenb) > MXDSMTL ) { continue; }
																pose.replace_residue(ersd,glu.residue(1),true);
																pose.set_chi(1,ersd,CHI1[ech1]);
																pose.set_chi(2,ersd,CHI2[ech2]);
																Vec const CD = pose.residue(ersd).xyz("CD");
																Real ange1 = angle_degrees(CD,cenb,ligb1);
																Real ange2 = angle_degrees(CD,cenb,ligb2);
																Real ange3 = angle_degrees(CD,cenb,ligi);
																Real ange4 = angle_degrees(CD,cenb,ligj);
																Real MXAGtmp = min(MXAGMTL,22.5);
																if ( ange1 < 90.0-MXAGMTL || (90.0+MXAGtmp < ange1 && ange1 < 135.0-MXAGtmp) || 135.0+MXAGMTL < ange1 ) continue;
																if ( ange2 < 90.0-MXAGMTL || (90.0+MXAGtmp < ange2 && ange2 < 135.0-MXAGtmp) || 135.0+MXAGMTL < ange2 ) continue;
																if ( ange3 < 90.0-MXAGMTL || (90.0+MXAGtmp < ange3 && ange3 < 135.0-MXAGtmp) || 135.0+MXAGMTL < ange3 ) continue;
																if ( ange4 < 90.0-MXAGMTL || (90.0+MXAGtmp < ange4 && ange4 < 135.0-MXAGtmp) || 135.0+MXAGMTL < ange4 ) continue;
																int n90 = 0, n135 = 0;
																if ( ange1 < 112.5 ) n90++; else n135++;
																if ( ange2 < 112.5 ) n90++; else n135++;
																if ( ange3 < 112.5 ) n90++; else n135++;
																if ( ange4 < 112.5 ) n90++; else n135++;
																bool biglu = true;
																if ( n90 != 2 || n135 != 2 ) biglu = false;
																Vec ea1(0,0,0),ea2(0,0,0);
																if ( biglu ) {
																	if ( ange1 > 112.5 && ange2 > 112.5 ) { ea1 = ligb1; ea2 = ligb2; }
																	if ( ange1 > 112.5 && ange3 > 112.5 ) { ea1 = ligb1; ea2 = ligi ; }
																	if ( ange1 > 112.5 && ange4 > 112.5 ) { ea1 = ligb1; ea2 = ligj ; }
																	if ( ange2 > 112.5 && ange3 > 112.5 ) { ea1 = ligb2; ea2 = ligi ; }
																	if ( ange2 > 112.5 && ange4 > 112.5 ) { ea1 = ligb2; ea2 = ligj ; }
																	if ( ange3 > 112.5 && ange4 > 112.5 ) { ea1 = ligi ; ea2 = ligj ; }
																	Real aer1 = dihedral_degrees(pose.residue(ersd).xyz("OE1"),pose.residue(ersd).xyz("CG"),CD,ea1);
																	Real aer2 = dihedral_degrees(pose.residue(ersd).xyz("OE2"),pose.residue(ersd).xyz("CG"),CD,ea2);
																	// TR << ange1 << " " << ange2 << " " << ange3 << " " << ange4 << " " << aer1 << " " << aer2 << std::endl;
																	pose.set_chi(3,ersd,pose.chi(3,ersd)+(aer1+aer2)/2.0);
																} else {
																	if ( ange1 > 112.5 ) { ea1 = ligb1; }
																	if ( ange2 > 112.5 ) { ea1 = ligb1; }
																	if ( ange3 > 112.5 ) { ea1 = ligb1; }
																	if ( ange4 > 112.5 ) { ea1 = ligb2; }
																	Real aer1 = dihedral_degrees(pose.residue(ersd).xyz("OE1"),pose.residue(ersd).xyz("CG"),CD,ea1);
																	pose.set_chi(3,ersd,pose.chi(3,ersd)+aer1);
																}

																clash = false;
																for ( Size a = 6; a <= pose.residue(ersd).nheavyatoms(); ++a ) {
																	Real TH = 8.0;
																	if ( a > pose.residue(ersd).nheavyatoms()-2 ) TH = 4.0;
																	else if ( !clashcheck.clash_check(pose.residue(ersd).xyz(a),ersd) ) { clash=true; }
																	for ( vector1<Size>::const_iterator b = batm.begin(); b != batm.end(); ++b ) {
																		if ( pose.residue(brsd).xyz(*b).distance_squared(pose.residue(ersd).xyz(a)) < TH ) clash=true;
																	}
																	for ( vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h ) {
																		if ( pose.residue(irsd).xyz(*h).distance_squared(pose.residue(ersd).xyz(a)) < TH ) clash=true;
																	}
																	for ( vector1<Size>::const_iterator h = hatm.begin(); h != hatm.end(); ++h ) {
																		if ( pose.residue(jrsd).xyz(*h).distance_squared(pose.residue(ersd).xyz(a)) < TH ) clash=true;
																	}

																	if ( clash ) break;
																}
																if ( clash ) { continue; }

																// if( ersd==76 && ech1==51 && ech2==29 ) {
																//  pose.dump_pdb("test0.pdb");
																//  utility_exit_with_message("tesljtl;sfj");
																// }

																//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

																for ( Size i = 1; i <= foundcen.size(); ++i ) {
																	if ( cen.distance_squared(foundcen[i]) < option[willmatch::match_overlap_dis]() &&
																			ori1.dot(foundori1[i])            > MATCH_OVERLAP_DOT                      &&
																			ori2.dot(foundori2[i])            > MATCH_OVERLAP_DOT                      &&
																			ori3.dot(foundori3[i])            > MATCH_OVERLAP_DOT                      &&
																			ori4.dot(foundori4[i])            > MATCH_OVERLAP_DOT                      &&
																			ori5.dot(foundori5[i])            > MATCH_OVERLAP_DOT                      &&
																			ori6.dot(foundori6[i])            > MATCH_OVERLAP_DOT                      ) {
																		// TR << "overlap!" << std::endl;
																		overlap = true;
																		break;
																	}
																}
																if ( overlap ) continue;
																count++;
																foundcen.push_back(cen);
																foundori1.push_back(ori1);
																foundori2.push_back(ori2);
																foundori3.push_back(ori3);
																foundori4.push_back(ori4);
																foundori5.push_back(ori5);
																foundori6.push_back(ori6);
																//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

																{
																	string outfname = utility::file_basename(infile)+"_FULL_"+string_of(brsd)+"_"+string_of(irsd)+"_"+string_of(jrsd)+"_"+string_of(ersd);
																	outfname += "_"+string_of(kch1)+"_"+string_of(kch2)+"_"+string_of(ich1)+"_"+string_of(ich2)+"_"+string_of(ide);
																	outfname += "_"+string_of(jch1)+"_"+string_of(jch2)+"_"+string_of(jde)+"_"+string_of(ech1)+"_"+string_of(ech2);
																	if ( biglu ) outfname += "_BIGLU";
																	outfname += ".pdb.gz";
																	TR << "HIT! " << outfname << std::endl;

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

																	//opose.dump_pdb( option[out::file::o]()+"/"+outfname+"_PRE.pdb.gz");
																	replace_pose_residue_copying_existing_coordinates(opose, irsd, opose.residue(1).residue_type_set().name_map(ide==1?"HIS_EEE":"HIS_DDD") );
																	replace_pose_residue_copying_existing_coordinates(opose, jrsd, opose.residue(1).residue_type_set().name_map(jde==1?"HIS_EEE":"HIS_DDD") );

																	AtomID bzn = AtomID(opose.residue(brsd).atom_index("ZN" ),brsd);
																	AtomID bne = AtomID(opose.residue(brsd).atom_index("NE1"),brsd);
																	AtomID bnn = AtomID(opose.residue(brsd).atom_index("NN1"),brsd);
																	AtomID inh = AtomID(opose.residue(irsd).atom_index(ide==1?"HD1":"HE2"),irsd);
																	AtomID inn = AtomID(opose.residue(irsd).atom_index(ide==1?"ND1":"NE2"),irsd);
																	AtomID jnh = AtomID(opose.residue(jrsd).atom_index(jde==1?"HD1":"HE2"),jrsd);
																	AtomID jnn = AtomID(opose.residue(jrsd).atom_index(jde==1?"ND1":"NE2"),jrsd);
																	AtomID ecg = AtomID(opose.residue(ersd).atom_index("CG" ),ersd);
																	AtomID ecd = AtomID(opose.residue(ersd).atom_index("CD" ),ersd);
																	AtomID oe1 = AtomID(opose.residue(ersd).atom_index("OE1"),ersd);
																	AtomID oe2 = AtomID(opose.residue(ersd).atom_index("OE2"),ersd);

																	Real d1 = pose.residue(brsd).xyz("ZN").distance(pose.residue(ersd).xyz("OE1"));
																	Real d2 = pose.residue(brsd).xyz("ZN").distance(pose.residue(ersd).xyz("OE2"));
																	core::scoring::func::FuncOP disfunc0( new core::scoring::func::HarmonicFunc( 0.0, 0.2 ) );
																	core::scoring::func::FuncOP disfunc( new core::scoring::func::HarmonicFunc( 2.2, 0.2 ) );
																	core::scoring::func::FuncOP angfunc0( new core::scoring::func::CircularHarmonicFunc(   0.0   , 0.2 ) );
																	core::scoring::func::FuncOP angfunc90( new core::scoring::func::CircularHarmonicFunc(   1.5707, 0.2 ) );
																	core::scoring::func::FuncOP angfunc180( new core::scoring::func::CircularHarmonicFunc( 2*1.5707, 0.2 ) );
																	core::scoring::func::FuncOP angfunc2( new core::scoring::func::CircularHarmonicFunc(   2.0943, 0.2 ) );

																	bool angbei = numeric::angle_degrees(opose.xyz(bne),opose.xyz(bzn),opose.xyz(inn)) < 135.0;
																	bool angbni = numeric::angle_degrees(opose.xyz(bnn),opose.xyz(bzn),opose.xyz(inn)) < 135.0;
																	bool angbej = numeric::angle_degrees(opose.xyz(bne),opose.xyz(bzn),opose.xyz(jnn)) < 135.0;
																	bool angbnj = numeric::angle_degrees(opose.xyz(bnn),opose.xyz(bzn),opose.xyz(jnn)) < 135.0;
																	bool anghij = numeric::angle_degrees(opose.xyz(inn),opose.xyz(bzn),opose.xyz(jnn)) < 135.0;

																	opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint(bzn,inh,disfunc0) ) ));
																	opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint(bzn,jnh,disfunc0) ) ));
																	opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(bne,bzn,inn,(angbei?angfunc90:angfunc180)) ) ));
																	opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(bnn,bzn,inn,(angbni?angfunc90:angfunc180)) ) ));
																	opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(bne,bzn,jnn,(angbej?angfunc90:angfunc180)) ) ));
																	opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(bnn,bzn,jnn,(angbnj?angfunc90:angfunc180)) ) ));
																	opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(inn,bzn,jnn,(anghij?angfunc90:angfunc180)) ) ));

																	// TR << "!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n" << std::endl;
																	// if(!biglu) continue;

																	Real omx1 = max( max( max( numeric::angle_degrees(opose.xyz(bne),opose.xyz(bzn),opose.xyz(oe1)),
																		numeric::angle_degrees(opose.xyz(bnn),opose.xyz(bzn),opose.xyz(oe1))),
																		numeric::angle_degrees(opose.xyz(inn),opose.xyz(bzn),opose.xyz(oe1))),
																		numeric::angle_degrees(opose.xyz(jnn),opose.xyz(bzn),opose.xyz(oe1)));
																	Real omx2 = max( max( max( numeric::angle_degrees(opose.xyz(bne),opose.xyz(bzn),opose.xyz(oe2)),
																		numeric::angle_degrees(opose.xyz(bnn),opose.xyz(bzn),opose.xyz(oe2))),
																		numeric::angle_degrees(opose.xyz(inn),opose.xyz(bzn),opose.xyz(oe2))),
																		numeric::angle_degrees(opose.xyz(jnn),opose.xyz(bzn),opose.xyz(oe2)));


																	bool angbe1 = numeric::angle_degrees(opose.xyz(bne),opose.xyz(bzn),opose.xyz(oe1)) != omx1;
																	bool angbn1 = numeric::angle_degrees(opose.xyz(bnn),opose.xyz(bzn),opose.xyz(oe1)) != omx1;
																	bool anghi1 = numeric::angle_degrees(opose.xyz(inn),opose.xyz(bzn),opose.xyz(oe1)) != omx1;
																	bool anghj1 = numeric::angle_degrees(opose.xyz(jnn),opose.xyz(bzn),opose.xyz(oe1)) != omx1;

																	bool angbe2 = numeric::angle_degrees(opose.xyz(bne),opose.xyz(bzn),opose.xyz(oe2)) != omx2;
																	bool angbn2 = numeric::angle_degrees(opose.xyz(bnn),opose.xyz(bzn),opose.xyz(oe2)) != omx2;
																	bool anghi2 = numeric::angle_degrees(opose.xyz(inn),opose.xyz(bzn),opose.xyz(oe2)) != omx2;
																	bool anghj2 = numeric::angle_degrees(opose.xyz(jnn),opose.xyz(bzn),opose.xyz(oe2)) != omx2;

																	if ( ((int)angbe1 + (int)angbn1 + (int)anghi1 + (int)anghj1) != 3 ) utility_exit_with_message("GLU1 ANG WRONG!!!!!!");
																	if ( ((int)angbe2 + (int)angbn2 + (int)anghi2 + (int)anghj2) != 3 ) utility_exit_with_message("GLU2 ANG WRONG!!!!!!");

																	if ( biglu || d1 < d2 ) {
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint(bzn,oe1,disfunc) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(bne,bzn,oe1,(angbe1?angfunc90:angfunc180)) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(bnn,bzn,oe1,(angbn1?angfunc90:angfunc180)) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(inn,bzn,oe1,(anghi1?angfunc90:angfunc180)) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(jnn,bzn,oe1,(anghj1?angfunc90:angfunc180)) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(    bzn,oe1,ecd,angfunc2) ) ));
																		if ( fabs(numeric::dihedral_degrees(opose.xyz(bzn),opose.xyz(oe1),opose.xyz(ecd),opose.xyz(ecg))) > 90.0 ) {
																			opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint(bzn,oe1,ecd,ecg,angfunc180) ) ));
																		} else opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint(bzn,oe1,ecd,ecg,angfunc0  ) ) ));
																	}
																	if ( biglu || d1 > d2 ) {
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint(bzn,oe2,disfunc) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(bne,bzn,oe2,(angbe2?angfunc90:angfunc180)) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(bnn,bzn,oe2,(angbn2?angfunc90:angfunc180)) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(inn,bzn,oe2,(anghi2?angfunc90:angfunc180)) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(jnn,bzn,oe2,(anghj2?angfunc90:angfunc180)) ) ));
																		opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint(    bzn,oe2,ecd,angfunc2) ) ));
																		if ( fabs(numeric::dihedral_degrees(opose.xyz(bzn),opose.xyz(oe2),opose.xyz(ecd),opose.xyz(ecg))) > 90.0 ) {
																			opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint(bzn,oe2,ecd,ecg,angfunc180) ) ));
																		} else opose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint(bzn,oe2,ecd,ecg,angfunc0  ) ) ));
																	}

																	//opose.dump_pdb("test.pdb");
																	//utility_exit_with_message("oarsetno");

																	//refine(opose,sf    ,brsd,irsd,jrsd,ersd);
																	refine(opose,sfhard,brsd,irsd,jrsd,ersd);
																	core::scoring::calpha_superimpose_pose(opose,native);
																	Real rms = core::scoring::CA_rmsd(opose,native);

																	sf->score(opose);
																	Real dbb = opose.energies().residue_total_energies(brsd)[core::scoring::fa_intra_rep];
																	Real dbi = opose.energies().residue_total_energies(irsd)[core::scoring::fa_intra_rep];//hlib->rotamer_energy( opose.residue(irsd), scratch );
																	Real dbj = opose.energies().residue_total_energies(jrsd)[core::scoring::fa_intra_rep];//hlib->rotamer_energy( opose.residue(jrsd), scratch );
																	Real dbe = opose.energies().residue_total_energies(ersd)[core::scoring::fa_intra_rep];//elib->rotamer_energy( opose.residue(ersd), scratch );
																	// if( dbe > DUN_THRESH ) continue;

																	Real hbi = opose.energies().residue_total_energies(irsd)[core::scoring::hbond_sc]+opose.energies().residue_total_energies(irsd)[core::scoring::hbond_bb_sc];
																	Real hbj = opose.energies().residue_total_energies(jrsd)[core::scoring::hbond_sc]+opose.energies().residue_total_energies(jrsd)[core::scoring::hbond_bb_sc];
																	Real hbe = opose.energies().residue_total_energies(ersd)[core::scoring::hbond_sc]+opose.energies().residue_total_energies(ersd)[core::scoring::hbond_bb_sc];

																	core::id::AtomID_Map<Real> atom_sasa;
																	vector1<Real> rsd_sasa;
																	core::scoring::calc_per_atom_sasa( opose, atom_sasa, rsd_sasa, 2.0, false );
																	// if( rsd_sasa[brsd] > 0.1 ) continue;

																	core::io::silent::SilentStructOP ss_out_all( new core::io::silent::ScoreFileSilentStruct );

																	Real decrholes = core::scoring::packing::compute_dec15_score(opose);
																	core::scoring::packing::HolesResult hr = core::scoring::packing::compute_holes_score(opose,hp_dec15);
																	Real locrholes = 0.0, count = 0.0;
																	for ( Size i = 1; i <= hr.atom_scores.n_atom(brsd); ++i ) { count += 1; locrholes += hr.atom_scores[AtomID(i,brsd)]; }
																	for ( Size i = 1; i <= hr.atom_scores.n_atom(irsd); ++i ) { count += 1; locrholes += hr.atom_scores[AtomID(i,irsd)]; }
																	for ( Size i = 1; i <= hr.atom_scores.n_atom(jrsd); ++i ) { count += 1; locrholes += hr.atom_scores[AtomID(i,jrsd)]; }
																	for ( Size i = 1; i <= hr.atom_scores.n_atom(ersd); ++i ) { count += 1; locrholes += hr.atom_scores[AtomID(i,ersd)]; }
																	locrholes /= count;

																	ss_out_all->fill_struct(opose,outfname);
																	ss_out_all->add_energy( "rms" , rms );
																	ss_out_all->add_energy( "rholes", decrholes );
																	ss_out_all->add_energy( "drholes", decrholes - natrholes);
																	ss_out_all->add_energy( "lrholes", locrholes );
																	ss_out_all->add_energy( "hb_his1", hbi );
																	ss_out_all->add_energy( "hb_his2", hbj );
																	ss_out_all->add_energy( "hb_glu" , hbe );
																	ss_out_all->add_energy( "bpy_chi1" , CHI1[kch1] );
																	ss_out_all->add_energy( "bpy_chi2" , CHI2[kch2] );
																	ss_out_all->add_energy( "dun_bpy" , dbb );
																	ss_out_all->add_energy( "dun_his1", dbi );
																	ss_out_all->add_energy( "dun_his2", dbj );
																	ss_out_all->add_energy( "dun_glu" , dbe );
																	ss_out_all->add_energy( "sasa_bpy" , rsd_sasa[brsd] );
																	ss_out_all->add_energy( "sasa_his1", rsd_sasa[irsd] );
																	ss_out_all->add_energy( "sasa_his2", rsd_sasa[jrsd] );
																	ss_out_all->add_energy( "sasa_glu" , rsd_sasa[ersd] );
																	sfd.write_silent_struct( *ss_out_all, option[out::file::o]()+"/"+option[ out::file::silent ]() );

																	utility::io::ozstream otmp(option[out::file::o]()+"/"+outfname);
																	otmp << "REMARK 666 MATCH TEMPLATE A BPY "+ObjexxFCL::format::I(4,brsd)+" MATCH MOTIF A HIS "+ObjexxFCL::format::I(4,irsd)+"  1  1" << std::endl;
																	otmp << "REMARK 666 MATCH TEMPLATE A BPY "+ObjexxFCL::format::I(4,brsd)+" MATCH MOTIF A HIS "+ObjexxFCL::format::I(4,jrsd)+"  2  1" << std::endl;
																	otmp << "REMARK 666 MATCH TEMPLATE A BPY "+ObjexxFCL::format::I(4,brsd)+" MATCH MOTIF A GLU "+ObjexxFCL::format::I(4,ersd)+"  3  1" << std::endl;
																	opose.dump_pdb(otmp);
																	otmp.close();

																	// static int outcount = 0;
																	// outcount++;
																	// if(outcount == 2)  utility_exit_with_message("test");
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

		//bhhout.close();
	}
	oprogress.close();
}


int main (int argc, char *argv[]) {

	try {

		devel::init(argc,argv);
		run();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
