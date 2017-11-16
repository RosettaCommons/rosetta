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

static basic::Tracer TR( "genmatch" );


inline Real const sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
	if ( sqdist > 9.*9. ) {
		return 0.0;
	} else if ( sqdist < 6.*6. ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0  - sqr( (dist - 6.) / (9. - 6.) ) );
	}
}


Real iface_check_c3(Pose & pose, Size nres, vector1<Size> const & iface_candidates) {
	Real num = 0;
	for ( vector1<Size>::const_iterator i=iface_candidates.begin(),ie=iface_candidates.end(); i != ie; ++i ) {
		for ( vector1<Size>::const_iterator j=iface_candidates.begin(),je=iface_candidates.end(); j != je; ++j ) {
			num += sigmoidish_neighbor(pose.residue(*i).xyz(5).distance_squared(pose.residue(*j+1*nres).xyz(5)));
			num += sigmoidish_neighbor(pose.residue(*i).xyz(5).distance_squared(pose.residue(*j+2*nres).xyz(5)));
		}
	}
	return num;
}

vector1<Size> read_res_list(string fn) {
	vector1<Size> l;
	if ( fn=="" ) return l;
	if ( fn=="_" ) return l;
	if ( fn.size()==1 && fn[0]==(char)0 ) return l;
	izstream in(fn);
	if ( !in.good() ) {
		utility_exit_with_message("can't open res list file '"+fn+"'");
	}
	Size r;
	while ( in >> r ) l.push_back(r);
	return l;
}

void repack(Pose & pose, Size nres, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->initialize_extra_rotamer_flags_from_command_line();
	for ( Size i = 1; i <= nres; ++i ) {
		if ( pose.residue(i).name3()=="BPY" ) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else {
			task->nonconst_residue_task(i).restrict_to_repacking();
		}
	}
	// TR << *task << std::endl;
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}

void design(Pose & pose, Size nres, ScoreFunctionOP sf) {
	core::id::AtomID_Map< bool > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, false );
	for ( Size ir = 1; ir <= pose.size(); ++ir ) {
		atom_map.set(AtomID(2,ir) , true );
		atom_map.set(AtomID(3,ir) , true );
		atom_map.set(AtomID(5,ir) , true );
	}
	core::id::AtomID_Map<Real> atom_sasa; utility::vector1<Real> sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, sasa, 2.3, false, atom_map );
	for ( Size i = 1; i <= sasa.size(); ++i ) if ( atom_sasa.n_atom(i) > 4 ) sasa[i] = atom_sasa[AtomID(5,i)];

	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys] = false;
	aas[core::chemical::aa_his] = false;
	// aas[core::chemical::aa_met] = false;
	aas[core::chemical::aa_pro] = false;
	aas[core::chemical::aa_gly] = false;
	aas[core::chemical::aa_gly] = false;
	if ( option[basic::options::OptionKeys::willmatch::exclude_ala]() ) aas[core::chemical::aa_ala] = false;
	if ( option[basic::options::OptionKeys::smhybrid::design_hydrophobic]() ) {
		aas[core::chemical::aa_ser] = false;
		aas[core::chemical::aa_thr] = false;
		aas[core::chemical::aa_asp] = false;
		aas[core::chemical::aa_glu] = false;
		aas[core::chemical::aa_lys] = false;
		aas[core::chemical::aa_arg] = false;
		aas[core::chemical::aa_asn] = false;
		aas[core::chemical::aa_gln] = false;
	}

	vector1<Size> fixed;
	if ( option[basic::options::OptionKeys::willmatch::fixed_res].user() ) {
		utility::io::izstream in(option[basic::options::OptionKeys::willmatch::fixed_res]());
		Size tmp;
		while ( in>>tmp ) fixed.push_back(tmp);
		in.close();
	}

	vector1<Size> interface;
	if ( option[basic::options::OptionKeys::willmatch::design_interface]() ) {
		for ( Size i = 1; i <= nres; ++i ) {
			if ( sasa.size() >= i && sasa[i] > 15.0 ) continue;
			AtomID aid(5,i);
			if ( pose.residue(i).nheavyatoms() < 5 ) aid.atomno() = 2;
			for ( Size j = nres+1; j <= 3*nres; ++j ) {
				AtomID aid2(5,j);
				if ( pose.residue(j).nheavyatoms() < 5 ) aid.atomno() = 2;
				if ( pose.xyz(aid).distance_squared(pose.xyz(aid2)) < 49 ) {
					interface.push_back(i);
				}
			}
		}
		for ( Size i = 1; i <= nres; ++i ) {
			Vec xyz = pose.xyz(AtomID(5,i));
			xyz.z() = 0;
			if ( xyz.length() < 8.0 ) interface.push_back(i);
		}
	}
	for ( Size i = 1; i <= nres; ++i ) {
		if ( pose.residue(i).name3()=="BPY" ) {
			task->nonconst_residue_task(i).prevent_repacking();
			// std::exit(-1);
		} else if ( std::find(fixed.begin(),fixed.end(),i)!=fixed.end() ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_ONE_STDDEV );
			task->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_ONE_STDDEV );
		} else if ( std::find(interface.begin(),interface.end(),i)!=interface.end() ) {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_ONE_STDDEV );
			task->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_ONE_STDDEV );
		}
	}
	for ( Size i = nres+1; i <= pose.size(); ++i ) {
		task->nonconst_residue_task(i).prevent_repacking();
	}
	TR << *task << std::endl;

	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}


void design_dyad(Pose & pose, Size const r1, Size const r2, ScoreFunctionOP sf, bool ex = false) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys] = false;
	aas[core::chemical::aa_his] = false;
	aas[core::chemical::aa_tyr] = false;

	aas[core::chemical::aa_ser] = false;
	aas[core::chemical::aa_thr] = false;
	aas[core::chemical::aa_asp] = false;
	aas[core::chemical::aa_glu] = false;
	aas[core::chemical::aa_lys] = false;
	aas[core::chemical::aa_arg] = false;
	aas[core::chemical::aa_asn] = false;
	aas[core::chemical::aa_gln] = false;
	// aas[core::chemical::aa_ala] = false;
	aas[core::chemical::aa_gly] = false;

	sf->score(pose);
	Real worig = sf->get_weight(core::scoring::res_type_constraint);
	if ( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);

	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);

	// aas[core::chemical::aa_trp] = false;
	if ( ex ) task->initialize_extra_rotamer_flags_from_command_line();

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue(i).name3()=="GLY" || pose.residue(i).name3()=="CYS" || pose.residue(i).name3()=="PRO" || pose.residue(i).name3()=="CYD" ) {
			task->nonconst_residue_task(i).prevent_repacking();
			continue;
		}
		Real d = 9e9; // set d by min dist to any heavy on r1,r2
		for ( Size ia = 1; ia <= pose.residue(i).nheavyatoms(); ++ia ) {
			numeric::xyzVector<Real> xyzi = pose.residue(i).xyz(ia);
			for ( Size ja = 1; ja <= pose.residue(r1).nheavyatoms(); ++ja ) {
				if ( xyzi.distance_squared(pose.residue(r1).xyz(ja)) < d ) d = xyzi.distance_squared( pose.residue(r1).xyz(ja) );
			}
			for ( Size ja = 1; ja <= pose.residue(r2).nheavyatoms(); ++ja ) {
				if ( xyzi.distance_squared(pose.residue(r2).xyz(ja)) < d ) d = xyzi.distance_squared( pose.residue(r2).xyz(ja) );
			}
		}
		if ( r1 == i || r2 == i ) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if ( d <  7.0* 7.0 ) {
			bool tmp = aas[pose.residue(i).aa()];
			aas[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
			aas[pose.residue(i).aa()] = tmp;
			task->nonconst_residue_task(i).or_include_current(true);
		} else if ( d < 13.0*13.0 ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_include_current(true);
		} else {
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	// TR << *task << std::endl;
	protocols::simple_moves::PackRotamersMover repack( sf, task );
	repack.apply(pose);

	sf->show(pose);

	pose.remove_constraints( res_cst );
	sf->set_weight(core::scoring::res_type_constraint,worig);
}


void minimize(Pose & pose, Size nres, Size , ScoreFunctionOP sf, int bb=0) {
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	// core::pose::symmetry::make_symmetric_movemap(pose,*movemap);
	movemap->set_chi(true);
	movemap->set_bb(false);
	movemap->set_jump(false);
	if ( bb==1 ) for ( Size i = 1; i <= nres; ++i ) if ( pose.secstruct(i)=='L' ) movemap->set_bb(i,true);
	if ( bb>=2 ) for ( Size i = 1; i <= nres; ++i ) movemap->set_bb(i,true);
	// movemap->set_chi(bpyres,false);

	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );

	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );

	m.apply(pose);

}


std::pair<Size,Size> makesplitwork_3bpy(Size total) {
	using namespace basic::options::OptionKeys;
	Size size1 = 1;
	Size size2 = total;
	if ( option[willmatch::splitwork].user() ) {
		Size part   = option[willmatch::splitwork]()[1];
		Size nparts = option[willmatch::splitwork]()[2];
		size1 = (part-1)*std::ceil(((Real)total)/(Real)nparts)+1;
		size2 = (part  )*std::ceil(((Real)total)/(Real)nparts);
		if ( option[in::file::s]().size() == 1 ) {
			Real frac1 = ((Real)part-1)/(Real)nparts;
			Real frac2 = ((Real)part  )/(Real)nparts;
			frac1 = 1.0 - sqrt(1.0-frac1);
			frac2 = 1.0 - sqrt(1.0-frac2);
			size1 = std::ceil(frac1*(Real)total)+1;
			size2 = std::ceil(frac2*(Real)total);
		}
	}
	TR << "SIZE " << size1 << " " << size2 << std::endl;
	return std::pair<Size,Size>(size1,size2);
}

std::pair<vector1<Size>,vector1<Size> > makesplitwork(Size total, Size total2 = 0) {
	using namespace basic::options::OptionKeys;
	vector1<Size> idx1;
	vector1<Size> idx2;
	Size part   = 1;
	Size nparts = 1;
	if ( option[willmatch::splitwork].user() ) {
		part   = option[willmatch::splitwork]()[1];
		nparts = option[willmatch::splitwork]()[2];
	}
	TR << "makesplitwork " << total << " " << total2 << " " << part << " " << nparts << std::endl;
	if ( total2 == 0 ) {
		for ( Size i = part; i <= total; i += nparts ) {
			idx1.push_back(i);
		}
	} else {
		if ( option[in::file::s]().size() == 1 ) {
			for ( Size i = part; i <= total; i += nparts ) {
				for ( Size j = i; j <= total; j += 1 ) {
					idx1.push_back(i);
					idx2.push_back(j);
				}
			}
		} else {
			for ( Size i = part; i <= total; i += nparts ) {
				for ( Size j = 1; j <= total2; j += 1 ) {
					idx1.push_back(i);
					idx2.push_back(j);
				}
			}
		}
	}
	return std::pair<vector1<Size>,vector1<Size> >(idx1,idx2);
}

std::pair<vector1<Size>,vector1<Size> > makesplitwork( utility::vector1<Size> res1, utility::vector1<Size> res2 ) {
	using namespace basic::options::OptionKeys;
	vector1<Size> idx1;
	vector1<Size> idx2;
	Size part   = 1;
	Size nparts = 1;
	if ( option[willmatch::splitwork].user() ) {
		part   = option[willmatch::splitwork]()[1];
		nparts = option[willmatch::splitwork]()[2];
	}
	TR << "makesplitwork " << res1.size() << " " << res2.size() << " " << part << " " << nparts << std::endl;
	for ( Size i = part; i <= res1.size(); i += nparts ) {
		for ( Size j = 1; j <= res2.size(); j += 1 ) {
			idx1.push_back(res1[i]);
			idx2.push_back(res2[j]);
		}
	}
	return std::pair<vector1<Size>,vector1<Size> >(idx1,idx2);
}


void myoptH(Pose & pose, ScoreFunctionOP sf) {
	add_lower_terminus_type_to_pose_residue(pose,1);
	add_upper_terminus_type_to_pose_residue(pose,pose.size());
	core::pack::optimizeH(pose,*sf);
	remove_lower_terminus_type_from_pose_residue(pose,1);
	remove_upper_terminus_type_from_pose_residue(pose,pose.size());
}


void run_3bpy() {
	using namespace basic::options::OptionKeys;

	// core::chemical::ResidueTypeSetCAP cen_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
	core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Pose init_fa;
	vector1<Pose> bpys(2);
	core::import_pose::pose_from_file(bpys[1] ,*fa_residue_set,option[in::file::s]()[2], core::import_pose::PDB_file);
	core::import_pose::pose_from_file(bpys[2] ,*fa_residue_set,option[in::file::s]()[3], core::import_pose::PDB_file);
	core::pose::remove_lower_terminus_type_from_pose_residue(bpys[1],1);
	core::pose::remove_lower_terminus_type_from_pose_residue(bpys[2],1);
	core::pose::remove_upper_terminus_type_from_pose_residue(bpys[1],1);
	core::pose::remove_upper_terminus_type_from_pose_residue(bpys[2],1);
	pose_from_file(init_fa,   *fa_residue_set,option[in::file::s]()[1], core::import_pose::PDB_file);
	Pose orig = init_fa;
	// pose_from_file(init_pose,*cen_residue_set,option[in::file::s]()[1], core::import_pose::PDB_file);
	Size nres = init_fa.size();
	core::chemical::ResidueType const & ala( init_fa.residue(1).residue_type_set().name_map("ALA") );
	for ( Size i = 1; i <= nres; ++i ) {
		core::pose::replace_pose_residue_copying_existing_coordinates(init_fa,i,ala);
	}
	vector1<Size> iface_candidates;
	vector1<Size> tmp = read_res_list(basic::options::option[basic::options::OptionKeys::willmatch::exclude_res1]());
	for ( Size i = 1; i <= init_fa.size(); ++i ) if ( std::find(tmp.begin(),tmp.end(),i)==tmp.end() ) iface_candidates.push_back(i);

	ImplicitFastClashCheck clashcheck(init_fa,basic::options::option[basic::options::OptionKeys::willmatch::clash_dis]());

	ScoreFunctionOP sf = core::scoring::get_score_function();

	string pref = "";
	vector1<numeric::xyzTriple<Size> > chis;
	if ( option[willmatch::chilist].user() ) {
		pref = "chilist.";
		utility::io::izstream in(option[willmatch::chilist]());
		Size tmp0,tmp1,tmp2;
		while ( in >> tmp0 >> tmp1 >> tmp2 ) {
			chis.push_back(numeric::xyzTriple<Size>(tmp0,tmp1,tmp2));
		}
		in.close();
	} else {
		for ( Size i = 2; i <= nres-1; ++i ) {
			for ( Real chi1 = 0.0; chi1 < 360.0; chi1+=option[willmatch::chi1_increment]() ) {
				for ( Real chi2 = 0.0; chi2 < 360.0; chi2+=option[willmatch::chi2_increment]() ) {
					chis.push_back(numeric::xyzTriple<Size>(i,chi1,chi2));
				}
			}
		}
	}
	std::pair<Size,Size> split = makesplitwork_3bpy(chis.size());
	Size last_rsd = 0;
	for ( Size ibpy = 1; ibpy <= 2; ibpy++ ) {
		Pose pose;
		Pose bpy = bpys[ibpy];
		for ( Size ichi = max(1lu,split.first); ichi <= min(chis.size(),split.second); ichi++ ) {
			Size irsd = chis[ichi].x();
			Size chi1 = chis[ichi].y();
			Size chi2 = chis[ichi].z();
			if ( irsd + nres*ibpy != last_rsd ) {
				TR << "PROGRESS: " << I(3,irsd) << " " << I(3,chi1) << " " << I(3,chi2) << std::endl;
				pose = init_fa;
				core::pose::replace_pose_residue_copying_existing_coordinates(pose,irsd,pose.residue(1).residue_type_set().name_map("BPY"));
				core::id::AtomID_Map< core::id::AtomID > atom_map;
				core::pose::initialize_atomid_map( atom_map, pose, core::id::AtomID::BOGUS_ATOM_ID() ); // maps every atomid to bogus atom
				for ( Size j = 6; j <= bpy.residue(1).natoms(); ++j ) atom_map[AtomID(j,irsd)] = AtomID(j,1);
				core::scoring::superimpose_pose( pose, bpy, atom_map );
				FoldTree ft;
				utility_exit_with_message("set_root_atomno not implemented!");
				// ft.set_root_atomno( bpy.residue(1).atom_index("NE1") );
				ft.add_edge(irsd,1,-1);
				ft.add_edge(irsd,pose.size(),-1);
				ft.reorder(irsd);
				pose.fold_tree(ft);
				core::pose::symmetry::make_symmetric_pose( pose );
				last_rsd = irsd + nres*ibpy;
			}
			pose.set_chi(1,irsd,chi1);
			pose.set_chi(2,irsd,chi2);
			if ( ! clashcheck.clash_check_trimer(pose,irsd) ) { /*TR << "clash fail" << std::endl;*/ continue; }
			Real iface = iface_check_c3(pose,nres,iface_candidates);
			if ( iface < (Size)option[willmatch::interface_size]() ) { /*TR << "iface fail " << iface << std::endl;*/ continue; }
			TR << "MATCH " << irsd << " " << chi1 << " " << chi2 << std::endl;
			Pose dpose = pose;
			for ( Size j = 1; j <= nres; ++j ) {
				if ( dpose.residue(j).name3()!="BPY" ) core::pose::replace_pose_residue_copying_existing_coordinates(dpose,j,orig.residue(j).type());
				else core::pose::replace_pose_residue_copying_existing_coordinates(dpose,j,fa_residue_set->name_map("BPY"));
			}

			design(dpose,nres,sf);
			minimize(dpose,nres,irsd,sf,0);
			design(dpose,nres,sf);

			core::id::AtomID_Map< bool > atom_map;
			core::pose::initialize_atomid_map( atom_map, dpose, false );
			for ( Size ir = 1; ir <= dpose.size(); ++ir ) {
				for ( Size ia = 1; ia <= dpose.residue(ir).nheavyatoms(); ++ia ) {
					Size t = dpose.residue(ir).atom_type_index(ia);
					if ( t==3 || t==4 || t==5 || t==6 || t==19 || t==20 ) atom_map.set(AtomID(ia,ir) , true );
				}
			}
			vector1<Size> reset_res;
			core::id::AtomID_Map<Real> atom_sasa; utility::vector1<Real> sasa;
			core::scoring::calc_per_atom_sasa( dpose, atom_sasa, sasa, 1.5, false, atom_map );
			for ( Size i = 1; i <= nres; ++i ) {
				if ( dpose.residue(i).name3()=="PHE" /*&& sasa[i] > 0.0*/ ) { TR << "SASA PHE " << sasa[i] << std::endl; reset_res.push_back(i); }
				if ( dpose.residue(i).name3()=="TRP" /*&& sasa[i] > 0.0*/ ) { TR << "SASA TRP " << sasa[i] << std::endl; reset_res.push_back(i); }
				if ( dpose.residue(i).name3()=="TYR" /*&& sasa[i] > 0.0*/ ) { TR << "SASA TYR " << sasa[i] << std::endl; reset_res.push_back(i); }
				if ( dpose.residue(i).name3()=="ILE" /*&& sasa[i] > 0.0*/ ) { TR << "SASA ILE " << sasa[i] << std::endl; reset_res.push_back(i); }
				if ( dpose.residue(i).name3()=="LEU" /*&& sasa[i] > 0.0*/ ) { TR << "SASA LEU " << sasa[i] << std::endl; reset_res.push_back(i); }
				if ( dpose.residue(i).name3()=="VAL" /*&& sasa[i] > 0.0*/ ) { TR << "SASA VAL " << sasa[i] << std::endl; reset_res.push_back(i); }
				if ( dpose.residue(i).name3()=="MET" /*&& sasa[i] > 0.0*/ ) { TR << "SASA MET " << sasa[i] << std::endl; reset_res.push_back(i); }
				if ( dpose.residue(i).name3()=="ALA" /*&& sasa[i] > 0.0*/ ) { TR << "SASA ALA " << sasa[i] << std::endl; reset_res.push_back(i); }
			}
			for ( Size i = 1; i <= reset_res.size(); ++i ) {
				if ( dpose.residue(reset_res[i]).name3()=="BPY" ) continue;
				core::pose::replace_pose_residue_copying_existing_coordinates(dpose,reset_res[i],orig.residue(reset_res[i]).type());
			}
			repack(dpose,nres,sf);

			Real rholes = core::scoring::packstat::compute_packing_score( dpose , 1 );
			string fname = pref+utility::file_basename(string(option[in::file::s]()[1])+"_"+"bpy"+string_of(ibpy)+"_"+string_of(irsd)+"_chi_"+string_of(chi1)+"_"+string_of(chi2)+".pdb");
			TR << "GENMATCH_SCORE " << I(3,irsd) << " " << I(3,chi1) << " " << I(3,chi2) << " " << I(3,iface) << " " << F(7,3,(*sf)(dpose)) << " " << F(7,3,dpose.energies().total_energies()[core::scoring::fa_atr]) << " " << F(5,4,rholes) << " " << fname << std::endl;
			utility::io::ozstream out(fname);
			dpose.dump_pdb(out);
			if ( basic::options::option[basic::options::OptionKeys::smhybrid::add_cavities]() ) {
				core::scoring::packstat::output_packstat_pdb(dpose,out);
			}
			out.close();

			// utility_exit_with_message("DEBUG");

			// for(Size i = 1; i <= nres; ++i) {
			//  for(Size j = 1; j <= dpose.residue(i).nchi(); ++j) {
			//    TR << "CHI " << i << " " << j << " " << dpose.chi(j,i) << std::endl;
			//  }
			// }
		}
	}

}

// find 2 HIS that can chelate a tetrahedral metal
void run_zn2his() {
	using namespace basic::options::OptionKeys;
	using namespace core::id;
	TR << "run_zn2his" << std::endl;

	core::chemical::ResidueTypeSetCAP cen_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
	core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Real tmpdis = option[willmatch::max_dis_metal]();
	Real const MXDSMTL = tmpdis*tmpdis;
	Real const MXAGMTL = option[willmatch::max_ang_metal]();
	Real const MATCH_OVERLAP_DOT = cos(numeric::conversions::radians(option[willmatch::match_overlap_ang]()));
	TR << "MATCH_OVERLAP_DOT: " << MATCH_OVERLAP_DOT << std::endl;
	string infile = option[in::file::s]()[1];
	Pose in_cen,in_fa;
	pose_from_file(in_fa, *fa_residue_set,infile, core::import_pose::PDB_file);
	pose_from_file(in_cen,*cen_residue_set,infile, core::import_pose::PDB_file);
	Size nres = in_cen.size();
	core::chemical::ResidueType const & ala( in_cen.residue(1).residue_type_set().name_map("ALA") );
	core::chemical::ResidueType const & alafa( in_fa.residue(1).residue_type_set().name_map("ALA") );
	// core::chemical::ResidueType const & hise( in_fa.residue(1).residue_type_set().name_map("HIS") );
	// core::chemical::ResidueType const & hisd( in_fa.residue(1).residue_type_set().name_map("HIS_D") );
	for ( Size i = 1; i <= nres; ++i ) {
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
	for ( Real i = 0; i < 360; i+= chi1incr ) CHI1.push_back(i);
	for ( Real i = 0; i < 360; i+= chi2incr ) CHI2.push_back(i);

	// setup HIS residues for checking
	Pose he,hd;
	core::pose::make_pose_from_sequence(he,"H[HIS]"  ,*fa_residue_set,false);
	core::pose::make_pose_from_sequence(hd,"H[HIS_D]",*fa_residue_set,false);
	hd.set_dof(DOF_ID(AtomID(hd.residue(1).atom_index("HD1"),1),D),2.1);
	he.set_dof(DOF_ID(AtomID(he.residue(1).atom_index("HE2"),1),D),2.1);
	// he.dump_pdb("he.pdb");
	// hd.dump_pdb("hd.pdb");
	ObjexxFCL::FArray3D<Vec> chi2cen(2,CHI1.size(),CHI2.size()),chi2ori(2,CHI1.size(),CHI2.size());
	ObjexxFCL::FArray2D<Vec> mcentroid(2,CHI1.size());
	Stub stub(he.xyz(AtomID(5,1)),he.xyz(AtomID(2,1)),he.xyz(AtomID(1,1)));
	// precompute cen and ori for hisd and hise
	for ( Size i = 1; i <= CHI1.size(); ++i ) {
		hd.set_chi(1,1,CHI1[i]);
		he.set_chi(1,1,CHI1[i]);
		mcentroid(1,i) = stub.global2local(  he.residue(1).xyz("CG") );
		mcentroid(2,i) = stub.global2local( (he.residue(1).xyz("CG")-he.residue(1).xyz("CB")).normalized()*5.4 + he.residue(1).xyz("CB") );
		for ( Size j = 1; j <= CHI2.size(); ++j ) {
			he.set_chi(2,1,CHI2[j]);
			hd.set_chi(2,1,CHI2[j]);
			chi2cen(1,i,j) = stub.global2local( hd.residue(1).xyz("HD1"));
			chi2ori(1,i,j) = stub.global2local( hd.residue(1).xyz("ND1"));
			chi2cen(2,i,j) = stub.global2local( he.residue(1).xyz("HE2"));
			chi2ori(2,i,j) = stub.global2local( he.residue(1).xyz("NE2"));
		}
	}
	// Pose tmp_pose = init_pose;
	//set residues to scan
	vector1<Size> residues;
	if ( option[willmatch::residues].user() ) {
		residues = option[willmatch::residues]();
	} else {
		for ( Size i = 1; i <= nres; ++i ) residues.push_back(i);
	}
	string fname = utility::file_basename(infile)+"_znhis_matches.pdb";
	if ( option[willmatch::splitwork].user() ) {
		fname = utility::file_basename(infile)+"_znhis_matches_part"+string_of(option[willmatch::splitwork]()[1])+".pdb";
	}
	utility::io::ozstream out(fname);
	out << "MODEL BASE" << endl;
	fa_pose.dump_pdb(out);
	out << "ENDMDL" << endl;

	std::pair<vector1<Size>,vector1<Size> > splitwork = makesplitwork(init_pose.size(),init_pose.size());
	vector1<Size> IRES = splitwork.first;
	vector1<Size> JRES = splitwork.second;
	assert(IRES.size()==JRES.size());

	Pose pose = fa_pose;

	Size count = 0;
	for ( Size iwork = 1; iwork <= IRES.size(); ++iwork ) {
		Size irsd = IRES[iwork];
		Size jrsd = JRES[iwork];
		// TR << "HIS HIS genmatch " << irsd << " " << jrsd << std::endl;
		// continue;
		// for(Size irsd = 1; irsd <= init_pose.size(); ++irsd) {
		if ( std::find(residues.begin(),residues.end(),irsd) == residues.end() ) continue; // should just loop?
		if ( std::find(residues.begin(),residues.end(),jrsd) == residues.end() ) continue; // should just loop?
		Stub s1(init_pose.xyz(AtomID(5,irsd)),init_pose.xyz(AtomID(2,irsd)),init_pose.xyz(AtomID(1,irsd)));
		Stub s2(init_pose.xyz(AtomID(5,jrsd)),init_pose.xyz(AtomID(2,jrsd)),init_pose.xyz(AtomID(1,jrsd)));
		// for(Size jrsd = irsd+1; jrsd <= init_pose.size(); ++jrsd) {
		vector1<Vec> foundcen;
		vector1<Vec> foundori;
		if ( init_pose.xyz(AtomID(5,irsd)).distance(init_pose.xyz(AtomID(5,jrsd))) > 11.1 ) continue;
		TR << irsd << " " << jrsd << " found " << count << std::endl;
		for ( Size ich1 = 1; ich1 <= CHI1.size(); ++ich1 ) {
			for ( Size jch1 = 1; jch1 <= CHI1.size(); ++jch1 ) {
				Vec mcen1d = s1.local2global(mcentroid(1,ich1));
				Vec mcen2d = s2.local2global(mcentroid(1,jch1));
				Vec mcen1e = s1.local2global(mcentroid(2,ich1));
				Vec mcen2e = s2.local2global(mcentroid(2,jch1));
				if ( mcen1d.distance_squared(mcen2d) > 16.5*16.5 &&
						mcen1d.distance_squared(mcen2e) > 15.0*15.0 &&
						mcen1e.distance_squared(mcen2d) > 15.0*15.0 &&
						mcen1e.distance_squared(mcen2e) > 14.0*14.0
						) continue;
				for ( Size ich2 = 1; ich2 <= CHI2.size(); ++ich2 ) {
					for ( Size jch2 = 1; jch2 <= CHI2.size(); ++jch2 ) {
						for ( Size ide = 1; ide <= 2; ide++ ) {
							Vec const ceni = s1.local2global(chi2cen(ide,ich1,ich2));
							for ( Size jde = 1; jde <= 2; jde++ ) {
								Vec const cenj = s2.local2global(chi2cen(jde,jch1,jch2));
								Real const d2 = ceni.distance_squared(cenj);
								if ( d2 > MXDSMTL ) {
									// TR << "dist fail " << sqrt(d2) << std::endl;
									continue;
								}
								Vec cen = (ceni+cenj)/2.0;
								Vec const orii = (s1.local2global(chi2ori(ide,ich1,ich2)) - cen).normalized();
								Vec const orij = (s2.local2global(chi2ori(jde,jch1,jch2)) - cen).normalized();
								Real angle = numeric::conversions::degrees(acos(orii.dot(orij)));
								// TR << "ORI " << angle << " " << orii << " " << orij << std::endl;
								if ( /*angle < 90.0-MXAGMTL || (90+MXAGMTL < angle && */(angle < 109.5-MXAGMTL) || 109.5+MXAGMTL < angle ) {
									// TR << "angle fail" << angle << " " << orii << " " << orij << std::endl;
									continue;
								}
								if ( ide==1 ) pose.replace_residue(irsd,hd.residue(1),true);
								else       pose.replace_residue(irsd,he.residue(1),true);
								if ( jde==1 ) pose.replace_residue(jrsd,hd.residue(1),true);
								else       pose.replace_residue(jrsd,he.residue(1),true);
								pose.set_chi(1,irsd,CHI1[ich1]);
								pose.set_chi(2,irsd,CHI2[ich2]);
								pose.set_chi(1,jrsd,CHI1[jch1]);
								pose.set_chi(2,jrsd,CHI2[jch2]);
								bool clash = false;
								for ( Size i = 6; i <= pose.residue(irsd).nheavyatoms(); ++i ) {
									if ( ide==1 && i == pose.residue(irsd).atom_index("ND1") ) continue;
									if ( ide==2 && i == pose.residue(irsd).atom_index("NE2") ) continue;
									for ( Size j = 6; j <= pose.residue(jrsd).nheavyatoms(); ++j ) {
										if ( jde==1 && j == pose.residue(jrsd).atom_index("ND1") ) continue;
										if ( jde==2 && j == pose.residue(jrsd).atom_index("NE2") ) continue;
										if ( pose.xyz(AtomID(i,irsd)).distance_squared(pose.xyz(AtomID(j,jrsd))) < 3.0*3.0 ) {
											clash=true;
											// TR << "HIS CLASH " << i << " " << j << std::endl;
											break;
											// tmp_pose.dump_pdb("test.pdb");
											// std::exit(-1);
										}
									}
									if ( clash ) break;
								}
								if ( clash ) continue;

								// core::id::AtomID_Mask mask;
								// core::pose::initialize_atomid_map(mask,pose,false);
								// mask.fill_with(irsd,true);
								// mask.fill_with(jrsd,true);
								// core::pack::optimize_H_and_notify(pose,mask);

								for ( Size iatm = 6; iatm <= pose.residue(irsd).nheavyatoms(); ++iatm ) {
									if ( !clashcheck.clash_check(pose.residue(irsd).xyz(iatm),irsd) ) clash=true;
								}
								for ( Size jatm = 6; jatm <= pose.residue(jrsd).nheavyatoms(); ++jatm ) {
									if ( !clashcheck.clash_check(pose.residue(jrsd).xyz(jatm),jrsd) ) clash=true;
								}
								if ( clash ) continue;
								// if(!clashcheck.clash_check(pose.residue(irsd).xyz("CE1"))) continue;
								// if(!clashcheck.clash_check(pose.residue(irsd).xyz("ND1"))) continue;
								// if(!clashcheck.clash_check(pose.residue(irsd).xyz("NE2"))) continue;
								// if(!clashcheck.clash_check(pose.residue(jrsd).xyz("CE1"))) continue;
								// if(!clashcheck.clash_check(pose.residue(jrsd).xyz("ND1"))) continue; // always clash w/CB
								// if(!clashcheck.clash_check(pose.residue(jrsd).xyz("NE2"))) continue;

								// check room for other ligands
								if ( !clashcheck.clash_check((ceni-orii*2)       ) ) { continue; }
								if ( !clashcheck.clash_check((cenj-orij*2)       ) ) { continue; }
								if ( !clashcheck.clash_check((ceni-orii*2-orij*2)) ) { continue; }
								if ( !clashcheck.clash_check((ceni-orii*4)       ) ) { continue; }
								if ( !clashcheck.clash_check((cenj-orij*4)       ) ) { continue; }
								if ( !clashcheck.clash_check((ceni-orii*4-orij*2)) ) { continue; }
								if ( !clashcheck.clash_check((cenj-orij*4-orii*2)) ) { continue; }
								if ( !clashcheck.clash_check((ceni-orii*4+orij*2)) ) { continue; }
								if ( !clashcheck.clash_check((cenj-orij*4+orii*2)) ) { continue; }

								myoptH(pose,sf);

								Vec ori = (orii+orij).normalized();
								bool overlap = false;
								for ( Size i = 1; i <= foundori.size(); ++i ) {
									if ( cen.distance_squared(foundcen[i]) < option[willmatch::match_overlap_dis]() &&
											ori.dot(foundori[i])              > MATCH_OVERLAP_DOT                      ) {
										TR << "overlap!" << std::endl;
										overlap = true;
										break;
									}
								}
								if ( overlap ) continue;
								count++;
								foundcen.push_back(cen);
								foundori.push_back(ori);

								string geom = "sqp";
								if ( angle > 100.25 ) geom = "tet";
								string tag = utility::file_basename(infile)+" "+geom+" "+" "+string_of(irsd)+" "+string_of(CHI1[ich1])+" "+string_of(CHI2[ich2])+" "+string_of(jrsd)+" "+string_of(CHI1[jch1])+" "+string_of(CHI2[jch2]);
								TR << "MATCH " << tag << std::endl;
								out << "MODEL HIS HIS "/*+string_of(ceni.distance(cenj))+" "+string_of(angle)+" "*/+tag << std::endl;
								Size one = 1, two = pose.residue(irsd).natoms()+1;
								core::conformation::Residue r1(pose.residue(irsd));
								core::conformation::Residue r2(pose.residue(jrsd));
								r1.chain(1);
								r2.chain(1);
								r1.seqpos(irsd);
								r2.seqpos(jrsd);
								core::io::pdb::dump_pdb_residue(r1,one,out);
								core::io::pdb::dump_pdb_residue(r2,two,out);
								// out << "HETATM" << I(5,9997) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,995) << "    " << F(8,3, ceni.x()) << F(8,3, ceni.y()) << F(8,3, ceni.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
								// out << "HETATM" << I(5,9998) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,996) << "    " << F(8,3, cenj.x()) << F(8,3, cenj.y()) << F(8,3, cenj.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
								// out << "HETATM" << I(5,9997) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,997) << "    " << F(8,3, (ceni+orii).x()) << F(8,3, (ceni+orii).y()) << F(8,3, (ceni+orii).z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
								// out << "HETATM" << I(5,9998) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,998) << "    " << F(8,3, (cenj+orij).x()) << F(8,3, (cenj+orij).y()) << F(8,3, (cenj+orij).z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';

								Vec viz;
								// viz = ceni-orii*2       ; out<<"HETATM"<<I(5,9991)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// viz = cenj-orij*2       ; out<<"HETATM"<<I(5,9992)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// viz = cenj-orii*2-orij*2; out<<"HETATM"<<I(5,9992)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// viz = ceni-orii*4       ; out<<"HETATM"<<I(5,9993)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// viz = cenj-orij*4       ; out<<"HETATM"<<I(5,9994)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// viz = ceni-orii*4-orij*2; out<<"HETATM"<<I(5,9995)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// viz = cenj-orij*4-orii*2; out<<"HETATM"<<I(5,9996)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// viz = ceni-orii*4+orij*2; out<<"HETATM"<<I(5,9997)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// viz = cenj-orij*4+orii*2; out<<"HETATM"<<I(5,9998)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
								// keep this one!
								viz = cen               ; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';


								out << "ENDMDL" << std::endl;
							}
						}
					}
				}
			}
		}

	}
	out.close();

}

// find 2 HIS that can chelate a tetrahedral metal
void run_tyr_his() {
	using namespace basic::options::OptionKeys;
	using namespace core::id;

	core::chemical::ResidueTypeSetCAP cen_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
	core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Real const MNDISHB = option[willmatch::min_dis_hb]();
	Real const MXDISHB = option[willmatch::max_dis_hb]();
	// Real const MXANGHB = option[willmatch::max_ang_metal]();
	Real const MATCH_OVERLAP_DOT = cos(numeric::conversions::radians(option[willmatch::match_overlap_ang]()));
	TR << "MATCH_OVERLAP_DOT: " << MATCH_OVERLAP_DOT << std::endl;
	string infile = option[in::file::s]()[1];
	Pose in_cen,in_fa;
	pose_from_file(in_fa, *fa_residue_set,infile, core::import_pose::PDB_file);
	Pose native = in_fa;
	pose_from_file(in_cen,*cen_residue_set,infile, core::import_pose::PDB_file);
	Size nres = in_cen.size();
	core::chemical::ResidueType const & ala  ( in_cen.residue(1).residue_type_set().name_map("ALA") );
	core::chemical::ResidueType const & alafa( in_fa.residue(1).residue_type_set().name_map("ALA") );
	// core::chemical::ResidueType const & hise ( in_fa.residue(1).residue_type_set().name_map("HIS") );
	// core::chemical::ResidueType const & hisd ( in_fa.residue(1).residue_type_set().name_map("HIS_D") );
	// core::chemical::ResidueType const & tyr  ( in_fa.residue(1).residue_type_set().name_map("TYR") );
	for ( Size i = 1; i <= nres; ++i ) {
		core::pose::replace_pose_residue_copying_existing_coordinates(in_cen,i,ala);
		core::pose::replace_pose_residue_copying_existing_coordinates(in_fa,i,alafa);
	}
	Pose const init_pose = in_cen;
	Pose const fa_pose = in_fa;
	ImplicitFastClashCheck clashcheck(init_pose,basic::options::option[basic::options::OptionKeys::willmatch::clash_dis]());

	// compute elegible residues
	Real BPR = 2.5;
	numeric::xyzVector<Real> metal_cen(option[willmatch::xyz_of_interest]()[1],option[willmatch::xyz_of_interest]()[2],option[willmatch::xyz_of_interest]()[3]);
	TR << "metal cen: " << metal_cen.x() << " " << metal_cen.y() << " " << metal_cen.z() << std::endl;
	utility::vector1<Size> forbid_res = option[willmatch::forbid_residues]();
	core::id::AtomID_Map<Real> bb_sasa = compute_bb_sasa(in_fa,BPR);
	utility::vector1<Size> allowed_res;
	for ( Size i = 1; i <= init_pose.size(); ++i ) {
		if ( native.residue(i).name3() == "GLY" ) continue;
		if ( native.residue(i).name3() == "CYD" ) continue;
		if ( native.residue(i).name3() == "CYS" ) continue;
		if ( native.residue(i).name3() == "PRO" ) continue;
		if ( bb_sasa[AtomID(5,i)]/12.56637/(2.0+BPR)/(2.0+BPR) > 0.10 ) continue;
		if ( std::find(forbid_res.begin(),forbid_res.end(),i) != forbid_res.end() ) continue;
		allowed_res.push_back(i);
	}
	utility::vector1<Size> allowed_res_tyr;
	for ( Size i = 1; i <= allowed_res.size(); ++i ) {
		// no tyr right next to forbidden (metal bound res)
		// if( std::find(forbid_res.begin(),forbid_res.end(),allowed_res[i]-1) != forbid_res.end() ) continue;
		// if( std::find(forbid_res.begin(),forbid_res.end(),allowed_res[i]+1) != forbid_res.end() ) continue;
		numeric::xyzVector<Real> cbxyz = init_pose.residue(allowed_res[i]).xyz("CB");
		if ( metal_cen.distance(cbxyz) <  5.0 ) continue;
		if ( metal_cen.distance(cbxyz) > 22.0 ) continue;
		bool place_for_re = false;
		for ( Size j = (Size)numeric::max(1,(int)allowed_res[i]-7); j <= (Size)numeric::min((int)init_pose.size(),(int)allowed_res[i]+7); ++j ) {
			if ( allowed_res[i] == j ) continue;
			if ( native.residue(j).name3() == "GLY" ) continue;
			if ( native.residue(j).name3() == "CYD" ) continue;
			if ( native.residue(j).name3() == "CYS" ) continue;
			if ( native.residue(j).name3() == "PRO" ) continue;
			if ( bb_sasa[AtomID(5,j)] < 0.2 ) continue;
			numeric::xyzVector<Real> cb2xyz = init_pose.residue(j).xyz("CB");
			if ( metal_cen.distance(cb2xyz) < 10.0 ) continue;
			if ( metal_cen.distance(cb2xyz) > 30.0 ) continue;
			if ( metal_cen.distance(cb2xyz) < 1.1*metal_cen.distance(cbxyz) ) continue;
			if ( cbxyz.distance(cb2xyz) < 6.0 ) continue;
			place_for_re = true;
		}
		if ( !place_for_re ) continue;
		allowed_res_tyr.push_back(allowed_res[i]);
	}
	TR << "ALLOW HIS:     resi " << allowed_res[1];
	for ( Size i = 2; i <= allowed_res.size(); ++i ) TR << "+" << allowed_res[i];
	TR << std::endl;
	TR << "ALLOW TYR:     resi " << allowed_res_tyr[1];
	for ( Size i = 2; i <= allowed_res_tyr.size(); ++i ) TR << "+" << allowed_res_tyr[i];
	TR << std::endl;


	ScoreFunctionOP sf = core::scoring::get_score_function();

	Real chi1incr = option[willmatch::chi1_increment]();
	Real chi2incr = option[willmatch::chi2_increment]();
	vector1<Real> CHI1,CHI2;
	for ( Real i = 0; i < 360; i+= chi1incr ) CHI1.push_back(i);
	for ( Real i = 0; i < 360; i+= chi2incr ) CHI2.push_back(i);

	// setup HIS residues for checking
	Pose he,hd,ty,ph;
	core::pose::make_pose_from_sequence(he,"H[HIS]"  ,*fa_residue_set,false);
	core::pose::make_pose_from_sequence(hd,"H[HIS_D]",*fa_residue_set,false);
	core::pose::make_pose_from_sequence(ty,"Y"       ,*fa_residue_set,false);
	core::pose::make_pose_from_sequence(ph,"F"       ,*fa_residue_set,false);

	ObjexxFCL::FArray3D<Vec> chi2cen(2,CHI1.size(),CHI2.size()),chi2ori(2,CHI1.size(),CHI2.size());
	ObjexxFCL::FArray2D<Vec> chi2centy(CHI1.size(),CHI2.size()),chi2ority(CHI1.size(),CHI2.size());
	ObjexxFCL::FArray2D<Vec> mcentroid(2,CHI1.size());
	Stub stub(he.xyz(AtomID(5,1)),he.xyz(AtomID(2,1)),he.xyz(AtomID(1,1)));
	// precompute cen and ori for hisd and hise
	for ( Size i = 1; i <= CHI1.size(); ++i ) {
		hd.set_chi(1,1,CHI1[i]);
		he.set_chi(1,1,CHI1[i]);
		ty.set_chi(1,1,CHI1[i]);
		mcentroid(1,i) = stub.global2local(  he.residue(1).xyz("CG") );
		mcentroid(2,i) = stub.global2local( (he.residue(1).xyz("CG")-he.residue(1).xyz("CB")).normalized()*5.4 + he.residue(1).xyz("CB") );
		for ( Size j = 1; j <= CHI2.size(); ++j ) {
			he.set_chi(2,1,CHI2[j]);
			hd.set_chi(2,1,CHI2[j]);
			ty.set_chi(2,1,CHI2[j]);
			chi2cen(1,i,j) = stub.global2local( hd.residue(1).xyz("ND1"));
			chi2ori(1,i,j) = stub.global2local( hd.residue(1).xyz("HD1"));
			chi2cen(2,i,j) = stub.global2local( he.residue(1).xyz("NE2"));
			chi2ori(2,i,j) = stub.global2local( he.residue(1).xyz("HE2"));
			chi2centy(i,j) = stub.global2local( ty.residue(1).xyz("OH"));
			chi2ority(i,j) = stub.global2local( ty.residue(1).xyz("CZ"));
		}
	}
	// Pose tmp_pose = init_pose;
	//set residues to scan
	vector1<Size> residues;
	if ( option[willmatch::residues].user() ) {
		residues = option[willmatch::residues]();
	} else {
		for ( Size i = 1; i <= nres; ++i ) residues.push_back(i);
	}
	string fname = utility::file_basename(infile)+"_tyrhis_matches.pdb";
	if ( option[willmatch::splitwork].user() ) {
		fname = utility::file_basename(infile)+"_tyrhis_matches_part"+string_of(option[willmatch::splitwork]()[1])+".pdb";
	}
	utility::io::ozstream out(fname);
	out << "MODEL BASE" << endl;
	fa_pose.dump_pdb(out);
	out << "ENDMDL" << endl;

	std::pair<vector1<Size>,vector1<Size> > splitwork = makesplitwork(allowed_res,allowed_res_tyr);
	vector1<Size> IRES = splitwork.first;
	vector1<Size> JRES = splitwork.second;
	assert(IRES.size()==JRES.size());

	Pose pose = fa_pose;

	Size count = 0;
	for ( Size iwork = 1; iwork <= IRES.size(); ++iwork ) {
		Size irsd = IRES[iwork];
		Size jrsd = JRES[iwork];
		if ( irsd==jrsd ) continue;
		// TR << "HIS HIS genmatch " << irsd << " " << jrsd << std::endl;
		// continue;
		// for(Size irsd = 1; irsd <= init_pose.size(); ++irsd) {
		if ( std::find(residues.begin(),residues.end(),irsd) == residues.end() ) continue; // should just loop?
		if ( std::find(residues.begin(),residues.end(),jrsd) == residues.end() ) continue; // should just loop?
		Stub s1(init_pose.xyz(AtomID(5,irsd)),init_pose.xyz(AtomID(2,irsd)),init_pose.xyz(AtomID(1,irsd)));
		Stub s2(init_pose.xyz(AtomID(5,jrsd)),init_pose.xyz(AtomID(2,jrsd)),init_pose.xyz(AtomID(1,jrsd)));
		// for(Size jrsd = irsd+1; jrsd <= init_pose.size(); ++jrsd) {
		vector1<Vec> foundcen;
		vector1<Vec> foundori;
		if ( init_pose.xyz(AtomID(2,irsd)).distance(init_pose.xyz(AtomID(2,jrsd))) > 20.0 ) continue;
		TR << irsd << " " << jrsd << " found " << count << std::endl;
		for ( Size ich1 = 1; ich1 <= CHI1.size(); ++ich1 ) {
			for ( Size jch1 = 1; jch1 <= CHI1.size(); ++jch1 ) {
				// skip for now
				// Vec mcen1d = s1.local2global(mcentroid(1,ich1));
				// Vec mcen2d = s2.local2global(mcentroid(1,jch1));
				// Vec mcen1e = s1.local2global(mcentroid(2,ich1));
				// Vec mcen2e = s2.local2global(mcentroid(2,jch1));
				// if( mcen1d.distance_squared(mcen2d) > 16.5*16.5 &&
				//     mcen1d.distance_squared(mcen2e) > 15.0*15.0 &&
				//     mcen1e.distance_squared(mcen2d) > 15.0*15.0 &&
				//     mcen1e.distance_squared(mcen2e) > 14.0*14.0
				// ) continue;
				for ( Size ich2 = 1; ich2 <= CHI2.size(); ++ich2 ) {
					for ( Size jch2 = 1; jch2 <= CHI2.size(); jch2 += 9 ) {
						for ( Size ide = 1; ide <= 2; ide++ ) {
							Vec const ceni = s1.local2global(chi2cen(ide,ich1,ich2));
							Vec const cenj = s2.local2global(chi2centy(  jch1,jch2));
							Real const d2 = ceni.distance(cenj);
							if ( MNDISHB > d2 || d2 > MXDISHB ) {
								// TR << "dist fail " << d2 << std::endl;
								continue;
							}
							// TR << "dist " << d2 << std::endl;
							Vec const orii = (ceni - s1.local2global(chi2ori(ide,ich1,ich2))).normalized();
							Vec const orij = (s2.local2global(chi2ority(  jch1,jch2)) - cenj).normalized();
							Real angle = numeric::conversions::degrees(acos(orii.dot(orij)));
							// TR << "ORI " << angle << " " << orii << " " << orij << std::endl;
							// if( angle < 109.5-MXANGHB || 109.5+MXANGHB < angle ) {
							if ( angle < 90 || 130 < angle ) {
								// TR << "angle fail " << angle << " " << orii << " " << orij << std::endl;
								continue;
							}
							// TR << "angle " << angle << std::endl;

							Vec p0 = cenj;
							Vec p1 = ceni;
							Vec p2 = s1.local2global(chi2ori(ide,ich1,ich2));
							if ( (p1-p0).dot(p1-p2) < 0.0 ) continue;
							Real dtoline = ((p0-p1).cross(p0-p2)).length() / (p2-p1).length();
							// TR << "dtoline " << dtoline << std::endl;
							if ( dtoline > option[willmatch::max_dis_hb_colinear]() ) continue;

							if ( ide==1 ) pose.replace_residue(irsd,hd.residue(1),true);
							else       pose.replace_residue(irsd,he.residue(1),true);
							pose.replace_residue(jrsd,ty.residue(1),true);
							pose.set_chi(1,irsd,CHI1[ich1]);
							pose.set_chi(2,irsd,CHI2[ich2]);
							pose.set_chi(1,jrsd,CHI1[jch1]);
							pose.set_chi(2,jrsd,CHI2[jch2]);

							// clash check Tyr/His
							bool clash = false;
							for ( Size i = 6; i <= pose.residue(irsd).nheavyatoms(); ++i ) {
								for ( Size j = 6; j <= pose.residue(jrsd).nheavyatoms(); ++j ) {
									if ( (i == pose.residue(jrsd).atom_index("OH")) && (
											(ide==1 && i == pose.residue(irsd).atom_index("ND1")) ||
											(ide==2 && i == pose.residue(irsd).atom_index("NE2")) ) ) continue;
									if ( pose.xyz(AtomID(i,irsd)).distance_squared(pose.xyz(AtomID(j,jrsd))) < 3.0*3.0 ) {
										clash=true;
										// TR << "HIS CLASH " << i << " " << j << std::endl;
										break;
										// tmp_pose.dump_pdb("test.pdb");
										// std::exit(-1);
									}
								}
								if ( clash ) break;
							}
							if ( clash ) continue;

							// core::id::AtomID_Mask mask;
							// core::pose::initialize_atomid_map(mask,pose,false);
							// mask.fill_with(irsd,true);
							// mask.fill_with(jrsd,true);
							// core::pack::optimize_H_and_notify(pose,mask);

							for ( Size iatm = 6; iatm <= pose.residue(irsd).nheavyatoms(); ++iatm ) {
								if ( !clashcheck.clash_check(pose.residue(irsd).xyz(iatm),irsd) ) clash=true;
							}
							for ( Size jatm = 6; jatm <= pose.residue(jrsd).nheavyatoms(); ++jatm ) {
								if ( !clashcheck.clash_check(pose.residue(jrsd).xyz(jatm),jrsd) ) clash=true;
							}
							if ( clash ) continue;
							// if(!clashcheck.clash_check(pose.residue(irsd).xyz("CE1"))) continue;
							// if(!clashcheck.clash_check(pose.residue(irsd).xyz("ND1"))) continue;
							// if(!clashcheck.clash_check(pose.residue(irsd).xyz("NE2"))) continue;
							// if(!clashcheck.clash_check(pose.residue(jrsd).xyz("CE1"))) continue;
							// if(!clashcheck.clash_check(pose.residue(jrsd).xyz("ND1"))) continue; // always clash w/CB
							// if(!clashcheck.clash_check(pose.residue(jrsd).xyz("NE2"))) continue;

							// TR << "CANDIDATE" << std::endl;
							// pose.dump_pdb("test.pdb");
							// std::exit(-1);


							Vec ori = (orii+orij).normalized();
							bool overlap = false;
							for ( Size i = 1; i <= foundori.size(); ++i ) {
								if ( cenj.distance_squared(foundcen[i]) < option[willmatch::match_overlap_dis]() &&
										ori.dot(foundori[i])               > MATCH_OVERLAP_DOT                      ) {
									TR << "overlap!" << std::endl;
									overlap = true;
									break;
								}
							}
							if ( overlap ) continue;

							myoptH(pose,sf);


							TR.Info << "check if Tyr OH/CZ exposed" << std::endl;
							Pose tmppose = native;
							tmppose.replace_residue(jrsd,ty.residue(1),true);
							tmppose.set_chi(1,jrsd,CHI1[jch1]);
							tmppose.set_chi(2,jrsd,CHI2[jch2]);
							utility::vector1<Real> rsd_sasa(tmppose.size(),0.0);
							core::id::AtomID_Map<Real> atom_sasa;
							core::id::AtomID_Map<bool> atom_mask;
							core::pose::initialize_atomid_map(atom_sasa,tmppose,0.0);
							core::pose::initialize_atomid_map(atom_mask,tmppose,false);
							for ( Size i = 1; i <= tmppose.size(); i++ ) {
								for ( Size j = 1; j <= tmppose.residue(i).nheavyatoms(); j++ ) atom_mask[AtomID(j,i)] = true;
							}
							core::scoring::calc_per_atom_sasa( tmppose, atom_sasa, rsd_sasa, 1.8, true, atom_mask );
							if ( atom_sasa[AtomID(tmppose.residue(jrsd).atom_index("OH"),jrsd)] > 0.0 ||
									/**/atom_sasa[AtomID(tmppose.residue(jrsd).atom_index("CZ"),jrsd)] > 0.0 ) {
								tmppose.dump_pdb("tyr_sasa_fail_"+string_of(count)+".pdb");
								continue;
							}


							TR << "remembering..." << std::endl;
							count++;
							foundcen.push_back(cenj);
							foundori.push_back(ori);

							string tag = utility::file_basename(infile)+" "+" "+string_of(irsd)+" "+string_of(CHI1[ich1])+" "+string_of(CHI2[ich2])+" "+string_of(jrsd)+" "+string_of(CHI1[jch1])+" "+string_of(CHI2[jch2]);
							TR << "MATCH " << tag << std::endl;
							out << "MODEL HIS HIS "/*+string_of(ceni.distance(cenj))+" "+string_of(angle)+" "*/+tag << std::endl;
							Size one = 1, two = pose.residue(irsd).natoms()+1;
							core::conformation::Residue r1(pose.residue(irsd));
							core::conformation::Residue r2(pose.residue(jrsd));
							r1.chain(1);
							r2.chain(1);
							r1.seqpos(1);
							r2.seqpos(2);
							core::io::pdb::dump_pdb_residue(r1,one,out);
							core::io::pdb::dump_pdb_residue(r2,two,out);
							// out << "HETATM" << I(5,9997) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,995) << "    " << F(8,3, ceni.x()) << F(8,3, ceni.y()) << F(8,3, ceni.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
							// out << "HETATM" << I(5,9998) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,996) << "    " << F(8,3, cenj.x()) << F(8,3, cenj.y()) << F(8,3, cenj.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
							// out << "HETATM" << I(5,9997) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,997) << "    " << F(8,3, (ceni+orii).x()) << F(8,3, (ceni+orii).y()) << F(8,3, (ceni+orii).z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
							// out << "HETATM" << I(5,9998) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,998) << "    " << F(8,3, (cenj+orij).x()) << F(8,3, (cenj+orij).y()) << F(8,3, (cenj+orij).z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';

							// Vec viz;
							// viz = ceni-orii*2       ; out<<"HETATM"<<I(5,9991)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = cenj-orij*2       ; out<<"HETATM"<<I(5,9992)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = cenj-orii*2-orij*2; out<<"HETATM"<<I(5,9992)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = ceni-orii*4       ; out<<"HETATM"<<I(5,9993)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = cenj-orij*4       ; out<<"HETATM"<<I(5,9994)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = ceni-orii*4-orij*2; out<<"HETATM"<<I(5,9995)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = cenj-orij*4-orii*2; out<<"HETATM"<<I(5,9996)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = ceni-orii*4+orij*2; out<<"HETATM"<<I(5,9997)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = cenj-orij*4+orii*2; out<<"HETATM"<<I(5,9998)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

							out << "ENDMDL" << std::endl;


							// short design step
							if ( true ) {
								TR << "desiging" << std::endl;
								Pose dpose = native;
								if ( ide==1 ) dpose.replace_residue(irsd,he.residue(1),true);
								else       dpose.replace_residue(irsd,hd.residue(1),true); // swapped....
								dpose.replace_residue(jrsd,ty.residue(1),true);
								// dpose.replace_residue(jrsd,ph.residue(1),true);
								dpose.set_chi(1,irsd,CHI1[ich1]);
								dpose.set_chi(2,irsd,CHI2[ich2]);
								dpose.set_chi(1,jrsd,CHI1[jch1]);
								dpose.set_chi(2,jrsd,CHI2[jch2]);
								// dpose.dump_pdb("1_init.pdb");

								// design_dyad(dpose,irsd,jrsd,sf,false);
								// dpose.dump_pdb("2_des1.pdb");
								// if( dpose.residue(irsd).name3()!="HIS" || dpose.residue(jrsd).name3()!="TYR" ) {
								//  utility_exit_with_message("tyr or his in dyad replaced!!!");
								// }
								//
								// // check burial
								// core::id::AtomID_Mask mask;
								// core::pose::initialize_atomid_map(mask,dpose,false);
								// for(Size i = 1; i <= dpose.size(); ++i) for(Size j = 1; j <= dpose.residue(i).nheavyatoms(); ++j ) mask[AtomID(j,i)] = true;
								// core::id::AtomID_Map<Real> atom_sasa; utility::vector1<Real> sasa;
								// TR << "computing sasa...." << std::endl;
								// core::scoring::calc_per_atom_sasa( dpose, atom_sasa, sasa, 1.8, false, mask);
								// for(Size i = 1; i <= sasa.size(); ++i) if( atom_sasa.n_atom(i) > 4 ) sasa[i] = atom_sasa[AtomID(5,i)];
								// if( /*atom_sasa[AtomID(dpose.residue(irsd).atom_index("CE1"),irsd)] > 0 ||*/
								//     /*atom_sasa[AtomID(dpose.residue(irsd).atom_index("CD2"),jrsd)] > 0 ||*/
								//     // atom_sasa[AtomID(dpose.residue(jrsd).atom_index("CD1"),jrsd)] > 0 ||
								//     // atom_sasa[AtomID(dpose.residue(jrsd).atom_index("CD2"),jrsd)] > 0 ||
								//     // atom_sasa[AtomID(dpose.residue(jrsd).atom_index("CE1"),jrsd)] > 0 ||
								//     // atom_sasa[AtomID(dpose.residue(jrsd).atom_index("CE2"),jrsd)] > 0 ||
								//     atom_sasa[AtomID(dpose.residue(jrsd).atom_index( "CZ"),jrsd)] > 0 //||
								//                  atom_sasa[AtomID(dpose.residue(jrsd).atom_index( "OH"),jrsd)] > 0 ) {
								//  TR << "rejecting non-buried TYR-HIS dyad" << std::endl;
								//  continue;
								// }

								// dpose.dump_pdb("test.pdb");
								// utility_exit_with_message("DBG!");

								TR << "desiging more" << std::endl;
								design_dyad(dpose,irsd,jrsd,sf,true);
								if ( dpose.residue(irsd).name3()!="HIS" || dpose.residue(jrsd).name3()!="TYR" ) {
									utility_exit_with_message("tyr or his in dyad replaced!!!");
								}

								// not needed if use tyr in first place....
								// // put OH back on tyr
								// dpose.replace_residue(jrsd,ty.residue(1),true);
								// dpose.set_chi(1,jrsd,CHI1[jch1]);
								// dpose.set_chi(2,jrsd,CHI2[jch2]);

								string fn = utility::file_basename(infile)+"_"+"_"+string_of(irsd)+"_"+string_of(CHI1[ich1])+"_"+string_of(CHI2[ich2])+"_"+string_of(jrsd)+"_"+string_of(CHI1[jch1])+"_"+string_of(CHI2[jch2])+".pdb";
								dpose.dump_pdb(fn);
								// std::exit(-1);

							}
						}
					}
				}
			}
		}

	}
	out.close();

}

void align_carboxyl_diiron_OLD(core::pose::Pose & pose, Size irsd, Size jrsd, Vec cen) {
	Vec cgi = pose.residue(irsd).xyz("CG");
	Vec cdi = pose.residue(irsd).xyz("CD");
	Vec o1i = pose.residue(irsd).xyz("OE1") - cdi;
	Vec o2i = pose.residue(irsd).xyz("OE2") - cdi;
	Vec cgj = pose.residue(jrsd).xyz("CG");
	Vec cdj = pose.residue(jrsd).xyz("CD");
	Vec o1j = pose.residue(jrsd).xyz("OE1") - cdj;
	Vec o2j = pose.residue(jrsd).xyz("OE2") - cdj;
	Vec axi = cdi-cgi;
	Vec axj = cdj-cgj;
	Mat Ri = rotation_matrix_degrees(axi,5.0);
	Mat Rj = rotation_matrix_degrees(axj,5.0);
	Vec mo1i(0,0,0),mo2i(0,0,0),mo1j(0,0,0),mo2j(0,0,0),ZERO(0,0,0);
	Real mx = 0.0;
	for ( Size i = 1; i <= 72; ++i ) {
		o1i = Ri * o1i;
		o2i = Ri * o2i;
		for ( Size j = 1; j <= 72; ++j ) {
			o1j = Rj * o1j;
			o2j = Rj * o2j;
			Real a1 = angle_degrees(o1i+cdi,cen,o1j+cdj);
			Real a2 = angle_degrees(o2i+cdi,cen,o1j+cdj);
			Real a3 = angle_degrees(o1i+cdi,cen,o2j+cdj);
			Real a4 = angle_degrees(o2i+cdi,cen,o2j+cdj);
			Real o = fabs(a1) - fabs(fabs(a2)-90.0) - fabs(fabs(a2)-90.0) - fabs(fabs(a2)-90.0);
			if ( o > mx ) {
				mo1i = o1i;
				mo2i = o2i;
				mo1j = o1j;
				mo2j = o2j;
				mx = o;
			}
		}
	}
	pose.set_xyz( AtomID(pose.residue(irsd).atom_index("OE1"),irsd), mo1i+cdi );
	pose.set_xyz( AtomID(pose.residue(irsd).atom_index("OE2"),irsd), mo2i+cdi );
	pose.set_xyz( AtomID(pose.residue(jrsd).atom_index("OE1"),jrsd), mo1j+cdj );
	pose.set_xyz( AtomID(pose.residue(jrsd).atom_index("OE2"),jrsd), mo2j+cdj );

}

void align_carboxyl_diiron(core::pose::Pose & pose, Size irsd, Size jrsd, Vec cen, Vec mori, bool posneg = true) {
	Real DIR = 1.0;
	if ( !posneg ) DIR = -1.0;
	Vec a = pose.residue(irsd).xyz("CD")-pose.residue(jrsd).xyz("CD");
	a -= projection_matrix(mori)*a;
	a = rotation_matrix_degrees(mori,DIR*45.0) * a;
	a = 2.0 * a.normalized();
	Vec b = rotation_matrix_degrees(mori,-DIR*90.0) * a;
	b = rotation_matrix_degrees(a,DIR*45.0) * b;
	b = 2*b.normalized();
	Real mn = 9e9, rmn = 0.0;
	for ( Real r = 0; r < 360.0; ++r ) {
		pose.set_chi(3,irsd,r);
		Real o = pose.residue(irsd).xyz("OE1").distance_squared(cen+a) + pose.residue(irsd).xyz("OE2").distance_squared(cen+b);
		if ( o < mn ) {
			rmn = r;
			mn = o;
		}
	}
	pose.set_chi(3,irsd,rmn);
	Vec a2 = -a;
	Vec b2 = rotation_matrix_degrees(a,DIR*90.0) * b;
	mn = 9e9; rmn = 0.0;
	for ( Real r = 0; r < 360.0; ++r ) {
		pose.set_chi(3,jrsd,r);
		Real o = pose.residue(jrsd).xyz("OE1").distance_squared(cen+a2) + pose.residue(jrsd).xyz("OE2").distance_squared(cen+b2);
		if ( o < mn ) {
			rmn = r;
			mn = o;
		}
	}
	pose.set_chi(3,jrsd,rmn);
	//
	//
	// utility::io::ozstream out("test.pdb");
	// Vec viz;
	// viz = cen;
	// out<<"HETATM"<<I(5,9997)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	// viz = cen+4*mori;
	// out<<"HETATM"<<I(5,9997)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"W"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	// viz = cen + a;
	// out<<"HETATM"<<I(5,9997)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"X"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	// viz = cen + b;
	// out<<"HETATM"<<I(5,9998)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"Y"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	// viz = cen + a2;
	// out<<"HETATM"<<I(5,9997)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"V"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	// viz = cen + b2;
	// out<<"HETATM"<<I(5,9998)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"U"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	// pose.dump_pdb(out);
	// out.close();
	//
	// utility_exit_with_message("debug...");
}

void run_diiron_glu() {
	using namespace basic::options::OptionKeys;
	using namespace core::id;

	core::chemical::ResidueTypeSetCAP cen_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
	core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Real tmpdis = option[willmatch::max_dis_metal]();
	Real const MXDSMTL = tmpdis*tmpdis;
	Real const MXAGMTL = option[willmatch::max_ang_metal]();
	Real const MATCH_OVERLAP_DOT = cos(numeric::conversions::radians(option[willmatch::match_overlap_ang]()));
	TR << "MATCH_OVERLAP_DOT: " << MATCH_OVERLAP_DOT << std::endl;
	string infile = option[in::file::s]()[1];
	Pose in_cen,in_fa;
	pose_from_file(in_fa, *fa_residue_set,infile, core::import_pose::PDB_file);
	pose_from_file(in_cen,*cen_residue_set,infile, core::import_pose::PDB_file);
	Size nres = in_cen.size();
	core::chemical::ResidueType const & ala( in_cen.residue(1).residue_type_set().name_map("ALA") );
	core::chemical::ResidueType const & alafa( in_fa.residue(1).residue_type_set().name_map("ALA") );
	// core::chemical::ResidueType const & hise( in_fa.residue(1).residue_type_set().name_map("HIS") );
	// core::chemical::ResidueType const & hisd( in_fa.residue(1).residue_type_set().name_map("HIS_D") );
	for ( Size i = 1; i <= nres; ++i ) {
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
	for ( Real i = 0; i < 360; i+= chi1incr ) CHI1.push_back(i);
	for ( Real i = 0; i < 360; i+= chi2incr ) CHI2.push_back(i);

	// setup HIS residues for checking
	Pose glu;
	core::pose::make_pose_from_sequence(glu,"E",*fa_residue_set,false);
	// glu.dump_pdb("glu.pdb");
	ObjexxFCL::FArray2D<Vec> chi2cen(CHI1.size(),CHI2.size()), chi2ori(CHI1.size(),CHI2.size());
	// ObjexxFCL::FArray2D<Vec> mcentroid(2,CHI1.size());
	Stub stub(glu.xyz(AtomID(5,1)),glu.xyz(AtomID(2,1)),glu.xyz(AtomID(1,1)));
	// precompute cen and ori for hisd and hise
	for ( Size i = 1; i <= CHI1.size(); ++i ) {
		glu.set_chi(1,1,CHI1[i]);
		//mcentroid(1,i) = stub.global2local(  he.residue(1).xyz("CG") );
		for ( Size j = 1; j <= CHI2.size(); ++j ) {
			glu.set_chi(2,1,CHI2[j]);
			Vec cd = glu.residue(1).xyz("CD");
			Vec cg = glu.residue(1).xyz("CG");
			Vec ori = (cd-cg).normalized();
			Vec cen = cd + 2.362957*ori;
			ori = cen + ori;
			chi2cen(i,j) = stub.global2local( cen );
			chi2ori(i,j) = stub.global2local( ori );
		}
	}
	// Pose tmp_pose = init_pose;
	//set residues to scan
	vector1<Size> residues;
	if ( option[willmatch::residues].user() ) {
		residues = option[willmatch::residues]();
	} else {
		for ( Size i = 1; i <= nres; ++i ) residues.push_back(i);
	}
	string fname = utility::file_basename(infile)+"_znhis_matches.pdb";
	if ( option[willmatch::splitwork].user() ) {
		fname = utility::file_basename(infile)+"_znhis_matches_part"+string_of(option[willmatch::splitwork]()[1])+".pdb";
	}
	utility::io::ozstream out(fname);
	out << "MODEL BASE" << endl;
	fa_pose.dump_pdb(out);
	out << "ENDMDL" << endl;

	std::pair<vector1<Size>,vector1<Size> > splitwork = makesplitwork(init_pose.size(),init_pose.size());
	vector1<Size> IRES = splitwork.first;
	vector1<Size> JRES = splitwork.second;
	assert(IRES.size()==JRES.size());

	ScoreFunctionOP cstsf = new core::scoring::ScoreFunction;
	cstsf->set_weight(core::scoring::atom_pair_constraint,1.0);


	Pose pose = fa_pose;

	Size count = 0;
	for ( Size iwork = 1; iwork <= IRES.size(); ++iwork ) {
		Size irsd = IRES[iwork];
		Size jrsd = JRES[iwork];
		vector1<Vec> foundcen_coarse;
		vector1<Vec> foundori_coarse;
		// TR << "HIS HIS genmatch " << irsd << " " << jrsd << std::endl;
		// continue;
		// for(Size irsd = 1; irsd <= init_pose.size(); ++irsd) {
		if ( std::find(residues.begin(),residues.end(),irsd) == residues.end() ) continue; // should just loop?
		if ( std::find(residues.begin(),residues.end(),jrsd) == residues.end() ) continue; // should just loop?
		Stub s1(init_pose.xyz(AtomID(5,irsd)),init_pose.xyz(AtomID(2,irsd)),init_pose.xyz(AtomID(1,irsd)));
		Stub s2(init_pose.xyz(AtomID(5,jrsd)),init_pose.xyz(AtomID(2,jrsd)),init_pose.xyz(AtomID(1,jrsd)));
		// for(Size jrsd = irsd+1; jrsd <= init_pose.size(); ++jrsd) {
		if ( init_pose.xyz(AtomID(5,irsd)).distance_squared(init_pose.xyz(AtomID(5,jrsd))) > 100.0 ) continue;
		TR << irsd << " " << jrsd << " found " << count << std::endl;
		for ( Size ich1 = 1; ich1 <= CHI1.size(); ++ich1 ) {
			for ( Size jch1 = 1; jch1 <= CHI1.size(); ++jch1 ) {
				// Vec mcen1d = s1.local2global(mcentroid(1,ich1));
				// Vec mcen2d = s2.local2global(mcentroid(1,jch1));
				// Vec mcen1e = s1.local2global(mcentroid(2,ich1));
				// Vec mcen2e = s2.local2global(mcentroid(2,jch1));
				// if( mcen1d.distance_squared(mcen2d) > 16.5*16.5 &&
				//     mcen1d.distance_squared(mcen2e) > 15.0*15.0 &&
				//     mcen1e.distance_squared(mcen2d) > 15.0*15.0 &&
				//     mcen1e.distance_squared(mcen2e) > 14.0*14.0
				// ) continue;
				for ( Size ich2 = 1; ich2 <= CHI2.size(); ++ich2 ) {
					for ( Size jch2 = 1; jch2 <= CHI2.size(); ++jch2 ) {
						Vec const ceni = s1.local2global(chi2cen(ich1,ich2));
						Vec const cenj = s2.local2global(chi2cen(jch1,jch2));
						Real const d2 = ceni.distance_squared(cenj);
						if ( d2 > MXDSMTL ) {
							// TR << "dist fail " << sqrt(d2) << std::endl;
							continue;
						}
						Vec cen = (ceni+cenj)/2.0;
						Vec const orii = s1.local2global(chi2ori(ich1,ich2))-ceni;
						Vec const orij = s2.local2global(chi2ori(jch1,jch2))-cenj;
						Real angle = numeric::conversions::degrees(acos(orii.dot(orij)));
						// TR << "ORI " << angle << " " << orii << " " << orij << std::endl;
						if ( angle < 120.0-MXAGMTL || 120.0+MXAGMTL < angle ) {
							// TR << "angle fail" << angle << " " << orii << " " << orij << std::endl;
							continue;
						}
						pose.replace_residue(irsd,glu.residue(1),true);
						pose.replace_residue(jrsd,glu.residue(1),true);
						pose.set_chi(1,irsd,CHI1[ich1]);
						pose.set_chi(2,irsd,CHI2[ich2]);
						pose.set_chi(1,jrsd,CHI1[jch1]);
						pose.set_chi(2,jrsd,CHI2[jch2]);
						bool clash = false;
						for ( Size i = 6; i <= pose.residue(irsd).nheavyatoms() - 2; ++i ) { // excludes OE1 and OE2
							for ( Size j = 6; j <= pose.residue(jrsd).nheavyatoms() - 2; ++j ) {
								if ( pose.xyz(AtomID(i,irsd)).distance_squared(pose.xyz(AtomID(j,jrsd))) < 3.0*3.0 ) {
									clash=true;
									// TR << "HIS CLASH " << i << " " << j << std::endl;
									break;
									// tmp_pose.dump_pdb("test.pdb");
									// std::exit(-1);
								}
							}
							if ( clash ) break;
						}
						if ( clash ) continue;

						// core::id::AtomID_Mask mask;
						// core::pose::initialize_atomid_map(mask,pose,false);
						// mask.fill_with(irsd,true);
						// mask.fill_with(jrsd,true);
						// core::pack::optimize_H_and_notify(pose,mask);


						for ( Size iatm = 6; iatm <= pose.residue(irsd).nheavyatoms()-2; ++iatm ) {
							if ( !clashcheck.clash_check(pose.residue(irsd).xyz(iatm),irsd) ) clash=true;
						}
						for ( Size jatm = 6; jatm <= pose.residue(jrsd).nheavyatoms()-2; ++jatm ) {
							if ( !clashcheck.clash_check(pose.residue(jrsd).xyz(jatm),jrsd) ) clash=true;
						}
						if ( clash ) continue;

						Vec mori = (orii+orij).normalized();
						// check room for other ligands
						if ( !clashcheck.clash_check((cen+0.0*mori)       ) ) { continue; }
						if ( !clashcheck.clash_check((cen+1.0*mori)       ) ) { continue; }
						if ( !clashcheck.clash_check((cen+2.0*mori)       ) ) { continue; }
						if ( !clashcheck.clash_check((cen+3.0*mori)       ) ) { continue; }
						if ( !clashcheck.clash_check((cen+4.0*mori)       ) ) { continue; }
						if ( !clashcheck.clash_check((cen+5.0*mori)       ) ) { continue; }
						if ( !clashcheck.clash_check((cen+6.0*mori)       ) ) { continue; }
						// if( !clashcheck.clash_check((cenj-orij*2)       ) ) { continue; }
						// if( !clashcheck.clash_check((ceni-orii*2-orij*2)) ) { continue; }
						// if( !clashcheck.clash_check((ceni-orii*4)       ) ) { continue; }
						// if( !clashcheck.clash_check((cenj-orij*4)       ) ) { continue; }
						// if( !clashcheck.clash_check((ceni-orii*4-orij*2)) ) { continue; }
						// if( !clashcheck.clash_check((cenj-orij*4-orii*2)) ) { continue; }
						// if( !clashcheck.clash_check((ceni-orii*4+orij*2)) ) { continue; }
						// if( !clashcheck.clash_check((cenj-orij*4+orii*2)) ) { continue; }

						// myoptH(pose,sf);

						align_carboxyl_diiron(pose,irsd,jrsd,cen,mori);
						//now clash check OE1/2
						for ( Size iatm = pose.residue(irsd).nheavyatoms()-1; iatm <= pose.residue(irsd).nheavyatoms(); ++iatm ) {
							if ( !clashcheck.clash_check(pose.residue(irsd).xyz(iatm),irsd) ) clash=true;
						}
						for ( Size jatm = pose.residue(jrsd).nheavyatoms()-1; jatm <= pose.residue(jrsd).nheavyatoms(); ++jatm ) {
							if ( !clashcheck.clash_check(pose.residue(jrsd).xyz(jatm),jrsd) ) clash=true;
						}
						if ( clash ) {
							clash = false;
							// align again, reverse way
							align_carboxyl_diiron(pose,irsd,jrsd,cen,mori,false);
							//now clash check OE1/2
							for ( Size iatm = pose.residue(irsd).nheavyatoms()-1; iatm <= pose.residue(irsd).nheavyatoms(); ++iatm ) {
								if ( !clashcheck.clash_check(pose.residue(irsd).xyz(iatm),irsd) ) clash=true;
							}
							for ( Size jatm = pose.residue(jrsd).nheavyatoms()-1; jatm <= pose.residue(jrsd).nheavyatoms(); ++jatm ) {
								if ( !clashcheck.clash_check(pose.residue(jrsd).xyz(jatm),jrsd) ) clash=true;
							}
						}
						if ( clash ) continue;

						for ( Size iatm = pose.residue(irsd).nheavyatoms()-1; iatm <= pose.residue(irsd).nheavyatoms(); ++iatm ) {
							for ( Size jatm = pose.residue(jrsd).nheavyatoms()-1; jatm <= pose.residue(jrsd).nheavyatoms(); ++jatm ) {
								if ( pose.residue(jrsd).xyz(jatm).distance_squared(pose.residue(irsd).xyz(iatm)) < 9.0 ) clash=true;
							}
						}
						if ( clash ) continue;


						// {
						//  // is cross better than what was in HIS/HIS?
						//  Vec ori = (orii.cross(orij)).normalized();
						//  bool overlap = false;
						//  for(Size i = 1; i <= foundori_coarse.size(); ++i) {
						//    if( cen.distance_squared(foundcen_coarse[i]) < option[willmatch::match_overlap_dis]()/2.0 &&
						//        ori .dot(foundori_coarse [i])            > MATCH_OVERLAP_DOT/2.0                      ){
						//      // TR << "overlap!" << std::endl;
						//      overlap = true;
						//      break;
						//    }
						//  }
						//  if(overlap) continue;
						//  foundcen_coarse.push_back(cen);
						//  foundori_coarse.push_back(ori);
						// }

						vector1<Vec> foundcen;
						vector1<Vec> foundori;
						vector1<Vec> foundori2;
						for ( Size krsd = 1; krsd <= init_pose.size(); ++krsd ) {
							if ( init_pose.xyz(AtomID(5,irsd)).distance_squared(init_pose.xyz(AtomID(5,krsd))) > 144.0 ) continue;
							if ( init_pose.xyz(AtomID(5,jrsd)).distance_squared(init_pose.xyz(AtomID(5,krsd))) > 144.0 ) continue;
							Stub s3(init_pose.xyz(AtomID(5,krsd)),init_pose.xyz(AtomID(2,krsd)),init_pose.xyz(AtomID(1,krsd)));
							for ( Size kch1 = 1; kch1 <= CHI1.size(); ++kch1 ) {
								for ( Size kch1 = 1; kch1 <= CHI1.size(); ++kch1 ) {
									for ( Size kch2 = 1; kch2 <= CHI2.size(); ++kch2 ) {
										Vec const cenk = s3.local2global(chi2cen(kch1,kch2));
										Vec const cenm2 = cen+2*mori;
										Real const d2 = cenm2.distance_squared(cenk);
										if ( d2 > MXDSMTL ) {
											// TR << "dist fail " << sqrt(d2) << std::endl;
											continue;
										}

										// sub in K glu
										pose.replace_residue(krsd,glu.residue(1),true);
										pose.set_chi(1,krsd,CHI1[kch1]);
										pose.set_chi(2,krsd,CHI2[kch2]);
										// Vec cen = (ceni+cenj)/2.0;
										Vec const orik = s3.local2global(chi2ori(kch1,kch2))-cenk;
										Real angle = angle_degrees(pose.residue(krsd).xyz("CD"),cenm2,cen);
										// TR << "ORI " << angle << " " << orii << " " << orij << std::endl;
										if ( angle < 90.0-MXAGMTL || 90.0+MXAGMTL < angle ) {
											// TR << "angle fail" << angle << " " << orii << " " << orij << std::endl;
											continue;
										}


										Real dangle = dihedral_degrees( pose.residue(krsd).xyz("CD"),cenm2,cen,pose.residue(irsd).xyz("CD") );
										bool kori = true;
										if ( dangle < 0 ) {
											kori = false;
											dangle += 180.0;
										}
										if ( dangle < 120.0-2.0*MXAGMTL || 120.0+1.0*MXAGMTL < dangle ) {
											// TR << "angle fail" << angle << " " << orii << " " << orij << std::endl;
											continue;
										}

										bool clash = false;
										for ( Size ij = 6; ij <= pose.residue(irsd).nheavyatoms() - 2; ++ij ) { // excludes OE1 and OE2
											for ( Size k = 6; k <= pose.residue(krsd).nheavyatoms() - 2; ++k ) {
												if ( pose.xyz(AtomID(ij,irsd)).distance_squared(pose.xyz(AtomID(k,krsd))) < 3.0*3.0 ) {
													clash=true;
													break;
												}
												if ( pose.xyz(AtomID(ij,jrsd)).distance_squared(pose.xyz(AtomID(k,krsd))) < 3.0*3.0 ) {
													clash=true;
													break;
												}
											}
											if ( clash ) break;
										}
										if ( clash ) continue;

										// core::id::AtomID_Mask mask;
										// core::pose::initialize_atomid_map(mask,pose,false);
										// mask.fill_with(irsd,true);
										// mask.fill_with(jrsd,true);
										// core::pack::optimize_H_and_notify(pose,mask);


										for ( Size katm = 6; katm <= pose.residue(krsd).nheavyatoms()-2; ++katm ) {
											if ( !clashcheck.clash_check(pose.residue(krsd).xyz(katm),krsd) ) clash=true;
										}
										if ( clash ) continue;

										// align OE1/2 in krsd
										{
											Real mn = 9e9, rmn = 0.0;
											for ( Real r = 0; r < 360.0; ++r ) {
												pose.set_chi(3,krsd,r);
												Real o = pose.residue(krsd).xyz("OE1").distance_squared(cen) + pose.residue(krsd).xyz("OE2").distance_squared(cen+4*mori);
												if ( o < mn ) {
													rmn = r;
													mn = o;
												}
											}
											pose.set_chi(3,krsd,rmn);
										}
										for ( Size katm = pose.residue(krsd).nheavyatoms()-1; katm <= pose.residue(krsd).nheavyatoms(); ++katm ) {
											if ( !clashcheck.clash_check(pose.residue(krsd).xyz(katm),krsd) ) clash=true;
										}
										if ( clash ) continue;

										{
											// utility::io::ozstream out("test.pdb");
											// pose.dump_pdb(out);
											// Vec viz;
											Vec axis = pose.residue(krsd).xyz("CD") - cenm2;
											axis -= projection_matrix(mori)*axis;
											axis = rotation_matrix_degrees(mori,90.0) * axis;
											Mat R = rotation_matrix_degrees(axis,180.0);
											for ( Size atm = pose.residue(krsd).nheavyatoms()-3; atm <= pose.residue(krsd).nheavyatoms()-2; ++atm ) {
												// viz = R*(pose.xyz(AtomID(atm,irsd))-cenm2)+cenm2; out<<"HETATM"<<I(5,9991)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"B"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
												// viz = R*(pose.xyz(AtomID(atm,jrsd))-cenm2)+cenm2; out<<"HETATM"<<I(5,9991)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"C"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
												// viz = R*(pose.xyz(AtomID(atm,krsd))-cenm2)+cenm2; out<<"HETATM"<<I(5,9991)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"D"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
												if ( !clashcheck.clash_check( R*(pose.xyz(AtomID(atm,irsd))-cenm2)+cenm2) ) clash=true;
												if ( !clashcheck.clash_check( R*(pose.xyz(AtomID(atm,jrsd))-cenm2)+cenm2) ) clash=true;
												if ( !clashcheck.clash_check( R*(pose.xyz(AtomID(atm,krsd))-cenm2)+cenm2) ) clash=true;
											}
											// viz = cen+2.0*mori      ; out<<"HETATM"<<I(5,9997)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"B"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
											// viz = cen+4.0*mori      ; out<<"HETATM"<<I(5,9998)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"C"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
											// viz = cen               ; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' ' <<  " ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
											// out.close();
											// utility_exit_with_message("testing");
											if ( clash ) continue;
										}


										// is cross better than what was in HIS/HIS?
										Vec ori = (orii.cross(orij)).normalized();
										bool overlap = false;
										for ( Size i = 1; i <= foundori.size(); ++i ) {
											if ( cen.distance_squared(foundcen[i]) < option[willmatch::match_overlap_dis]() &&
													ori .dot(foundori [i])            > MATCH_OVERLAP_DOT                      &&
													orik.dot(foundori2[i])            > MATCH_OVERLAP_DOT                      ) {
												// TR << "overlap!" << std::endl;
												overlap = true;
												break;
											}
										}
										if ( overlap ) continue;
										count++;
										foundcen.push_back(cen);
										foundori.push_back(ori);
										foundori2.push_back(orik);


										string geom = "EE4";
										if ( kori ) geom += "_POS";
										else     geom += "_NEG";
										string tag = utility::file_basename(infile)+" "+geom+" "+" "+string_of(irsd)+" "+string_of(CHI1[ich1])+" "+string_of(CHI2[ich2])+" "+string_of(jrsd)+" "+string_of(CHI1[jch1])+" "+string_of(CHI2[jch2])+" "+string_of(krsd)+" "+string_of(CHI1[kch1])+" "+string_of(CHI2[kch2])+" END";
										TR << "MATCH " << tag << std::endl;
										out << "MODEL HIS HIS "/*+string_of(ceni.distance(cenj))+" "+string_of(angle)+" "*/+tag << std::endl;
										Size one = 1, two = pose.residue(irsd).natoms()+1, three = 2*pose.residue(irsd).natoms()+1;
										core::conformation::Residue r1(pose.residue(irsd));
										core::conformation::Residue r2(pose.residue(jrsd));
										core::conformation::Residue r3(pose.residue(krsd));
										r1.chain(1);
										r2.chain(1);
										r3.chain(1);
										r1.seqpos(irsd);
										r2.seqpos(jrsd);
										r3.seqpos(krsd);
										core::io::pdb::dump_pdb_residue(r1,one  ,out);
										core::io::pdb::dump_pdb_residue(r2,two  ,out);
										core::io::pdb::dump_pdb_residue(r3,three,out);
										// out << "HETATM" << I(5,9997) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,995) << "    " << F(8,3, ceni.x()) << F(8,3, ceni.y()) << F(8,3, ceni.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
										// out << "HETATM" << I(5,9998) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,996) << "    " << F(8,3, cenj.x()) << F(8,3, cenj.y()) << F(8,3, cenj.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
										// out << "HETATM" << I(5,9997) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,997) << "    " << F(8,3, (ceni+orii).x()) << F(8,3, (ceni+orii).y()) << F(8,3, (ceni+orii).z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
										// out << "HETATM" << I(5,9998) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,998) << "    " << F(8,3, (cenj+orij).x()) << F(8,3, (cenj+orij).y()) << F(8,3, (cenj+orij).z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';

										Vec viz;
										// viz = ceni-orii*2       ; out<<"HETATM"<<I(5,9991)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = cenj-orij*2       ; out<<"HETATM"<<I(5,9992)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = cenj-orii*2-orij*2; out<<"HETATM"<<I(5,9992)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = ceni-orii*4       ; out<<"HETATM"<<I(5,9993)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = cenj-orij*4       ; out<<"HETATM"<<I(5,9994)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = ceni-orii*4-orij*2; out<<"HETATM"<<I(5,9995)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = cenj-orij*4-orii*2; out<<"HETATM"<<I(5,9996)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = ceni-orii*4+orij*2; out<<"HETATM"<<I(5,9997)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = cenj-orij*4+orii*2; out<<"HETATM"<<I(5,9998)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// keep this one!
										// viz = ceni           ; out<<"HETATM"<<I(5,9994)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"B"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = cenj           ; out<<"HETATM"<<I(5,9995)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"C"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = ceni+orii      ; out<<"HETATM"<<I(5,9996)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"D"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										// viz = cenj+orij      ; out<<"HETATM"<<I(5,9997)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"E"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										viz = cen+2.0*mori      ; out<<"HETATM"<<I(5,9997)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"B"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										viz = cen+4.0*mori      ; out<<"HETATM"<<I(5,9998)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"C"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
										viz = cen               ; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';


										out << "ENDMDL" << std::endl;

									}
								}
							}
						}
					}
				}
			}
		}

	}
	out.close();

}


int main (int argc, char *argv[]) {

	try {


		using namespace basic::options::OptionKeys;
		using namespace core::id;
		devel::init(argc,argv);


		run_zn2his();
		// run_diiron_glu();
		// run_3bpy();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


