// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/crossmatch.cc
/// @brief crosses matches fast

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/crossmatch.OptionKeys.gen.hh>
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
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/packstat/sasa_dot_locations.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>
#include <sstream>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>

#include <apps/pilot/will/will_util.ihh>

// #include <devel/init.hh>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

using core::id::AtomID;
using basic::options::option;
using core::kinematics::Stub;
using core::kinematics::FoldTree;
using core::pose::Pose;
using core::pose::PoseAP;
using core::pose::PoseCAP;
using core::pose::PoseCOP;
using core::pose::PoseOP;
using core::Real;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::conversions::radians;
using numeric::constants::r::pi;
using numeric::rotation_matrix_degrees;
using numeric::xyzVector;
using numeric::xyzMatrix;
using numeric::x_rotation_matrix_degrees;
using numeric::y_rotation_matrix_degrees;
using numeric::z_rotation_matrix_degrees;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::format::LJ;
using ObjexxFCL::format::SS;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using std::pair;
using utility::io::izstream;
using utility::io::ozstream;
using utility::pointer::owning_ptr;
using utility::pointer::access_ptr;
using utility::pointer::ReferenceCount;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;
typedef numeric::xyzMatrix<Real> Mat;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;


using core::pose::Pose;
using core::scoring::ScoreFunctionOP;

static thread_local basic::Tracer TR( "crossmatch_3e" );

static core::io::silent::SilentFileData sfd;

void design_homodimer(Pose & pose, ScoreFunctionOP sf, vector1<Size> const & matchres, bool designall = false ){
	using namespace core::pack::task;

	Real worig = sf->get_weight(core::scoring::res_type_constraint);
	if( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);

	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	// task->initialize_extra_rotamer_flags_from_command_line();
	vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys] = false;
	aas[core::chemical::aa_his] = false;
	// aas[core::chemical::aa_met] = false;
	aas[core::chemical::aa_pro] = false;
	aas[core::chemical::aa_gly] = false;
	if(option[basic::options::OptionKeys::willmatch::exclude_ala]()) aas[core::chemical::aa_ala] = false;
	if(!designall) {
		if(option[basic::options::OptionKeys::smhybrid::design_hydrophobic]()) {
			aas[core::chemical::aa_ser] = false;
			aas[core::chemical::aa_thr] = false;
			aas[core::chemical::aa_asp] = false;
			aas[core::chemical::aa_glu] = false;
			aas[core::chemical::aa_lys] = false;
			aas[core::chemical::aa_arg] = false;
			aas[core::chemical::aa_asn] = false;
			aas[core::chemical::aa_gln] = false;
		}
	}

	core::conformation::symmetry::SymmetryInfoCOP symminfo = core::pose::symmetry::symmetry_info(pose);
	Size end_of_prot_1 = symminfo->get_nres_subunit();
	Size end_of_prot_2 = 2*end_of_prot_1;
	vector1<Size> interface;
	for(Size i = 1; i <= end_of_prot_1; ++i) {
		// if( sasa.size() >= i && sasa[i] > option[basic::options::OptionKeys::willmatch::cb_sasa_thresh]() ) continue;
		AtomID aid(5,i);
		if(pose.residue(i).nheavyatoms() < 5) continue; // don't design GLY or VIRT
		for(Size j = end_of_prot_1+1; j <= end_of_prot_2; ++j) {
			// if( nres > end_of_prot_2 && sasa[i] > option[basic::options::OptionKeys::willmatch::cb_sasa_thresh]()/2.0 ) continue;
			AtomID aid2(5,j);
			if(pose.residue(j).nheavyatoms() < 5) continue; // don't design GLY or VIRT
			if(pose.xyz(aid).distance_squared(pose.xyz(aid2)) < 100.0) {
				interface.push_back(i);
			}
		}
	}

	// TR << "INTERFACE ";
	for(Size i = 1; i <= end_of_prot_1; ++i) {
		if(std::find(matchres.begin(),matchres.end(),i)!=matchres.end()) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if(std::find(interface.begin(),interface.end(),i)!=interface.end()){
			// TR << i << " ";
			bool tmp = aas[pose.residue(i).aa()];
			aas[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
			aas[pose.residue(i).aa()] = tmp;
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
			task->nonconst_residue_task(i).or_include_current(true);
		} else {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_include_current(true);
		}
	}
	// TR << std::endl;

	core::pack::make_symmetric_PackerTask_by_truncation(pose,task);
	// TR << *task << std::endl;

	// pose.dump_pdb("test.pdb");
	// if(uniform() > 0.2) std::exit(-1);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);

	// cleanup 2
	pose.remove_constraints( res_cst );
	sf->set_weight(core::scoring::res_type_constraint,worig);

}




int main (int argc, char *argv[])
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init(argc,argv);

	ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction( core::scoring::get_score_function() );

	Pose nat;
	core::import_pose::pose_from_pdb(nat,option[willmatch::native1]());
	vector1<Pose> m;
	core::import_pose::pose_from_pdb(m,*core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD),option[in::file::s]()[1]);
	for(Size ilig = 2; ilig <= m.size(); ++ilig) {
		Pose base(nat);
		Pose lig(m[ilig]);
		// TR << lig.pdb_info()->number(1) << " " << lig.pdb_info()->number(2) << " " << lig.pdb_info()->number(3) << std::endl;
		core::pose::remove_lower_terminus_type_from_pose_residue(lig,1);
		core::pose::remove_lower_terminus_type_from_pose_residue(lig,2);
		core::pose::remove_lower_terminus_type_from_pose_residue(lig,3);
		core::pose::remove_upper_terminus_type_from_pose_residue(lig,1);
		core::pose::remove_upper_terminus_type_from_pose_residue(lig,2);
		core::pose::remove_upper_terminus_type_from_pose_residue(lig,3);
		base.replace_residue( lig.pdb_info()->number(1), lig.residue(1), true );
		base.replace_residue( lig.pdb_info()->number(2), lig.residue(2), true );
		base.replace_residue( lig.pdb_info()->number(3), lig.residue(3), true );
		Vec cen = lig.residue(4).xyz(1);
		Vec ori = (lig.residue(5).xyz(1) - cen).normalized();
		Vec axs = (lig.residue(3).xyz("CD") - cen).normalized();
		axs -= projection_matrix(ori)*axs;
		axs = rotation_matrix_degrees(ori,90.0) * axs;
		Mat R = rotation_matrix_degrees(axs,180.0);
		core::conformation::symmetry::SymmDataOP symmdata = core::pose::symmetry::symm_data_from_axis(base.n_residue(),2,axs,cen);
		core::pose::symmetry::make_symmetric_pose(base,*symmdata);
		bool clash = false;
		for(Size ir = 1; ir <= nat.n_residue(); ++ir) {
			for(Size jr = nat.n_residue()+1; jr <= 2*nat.n_residue(); ++jr) {
				for(Size ia = 1; ia <= 4; ia++) {
					for(Size ja = 1; ja <= 4; ja++) {
						if( base.xyz(AtomID(ia,ir)).distance_squared(base.xyz(AtomID(ja,jr))) < 9.0 ) clash = true;
					}
				}
				if(clash) break;
			}
			if(clash) break;
		}
		if(clash) continue;

		vector1<Size> matchres; matchres.push_back(lig.pdb_info()->number(1)); matchres.push_back(lig.pdb_info()->number(2)); matchres.push_back(lig.pdb_info()->number(3));
		// base.dump_pdb("init.pdb");
		design_homodimer(base,sf,matchres);
		sf->score(base);
		// base.dump_pdb("desn.pdb");


		Size nmut = 0;
		for(Size ir = 1; ir <= nat.n_residue(); ++ir) {
			if( nat.residue(ir).name3() != base.residue(ir).name3() ) nmut++;
		}
		core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
		ss_out->fill_struct(base, "3E_homo_"+string_of(ilig-1));
	 	ss_out->add_energy( "iface", nmut );
		sfd.write_silent_struct( *ss_out, option[ basic::options::OptionKeys::out::file::silent ]() );

		// base.energies().show(cout);

		// utility_exit_with_message("testing");

		TR << "dumping pdb " << ilig << " " << "3E_homo_"+string_of(ilig-1)+".pdb" << std::endl;
		utility::io::ozstream out("3E_homo_"+string_of(ilig-1)+".pdb");
		base.dump_pdb(out);
		Vec viz;
		viz = cen+2*ori; out <<"HETATM"<<I(5,9998)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"B"<<I(4,base.n_residue()+1)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		viz = cen-2*ori; out <<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"B"<<I(4,base.n_residue()+2)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		out.close();
	}

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
















