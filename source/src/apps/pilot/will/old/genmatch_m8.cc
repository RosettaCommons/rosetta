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
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
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
	protocols::minimization_packing::symmetry::SymPackRotamersMover repack( sf, task );
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

	protocols::minimization_packing::symmetry::SymPackRotamersMover repack( sf, task );
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
	protocols::minimization_packing::PackRotamersMover repack( sf, task );
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

	protocols::minimization_packing::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );

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

void align_carboxyl_m8(core::pose::Pose & pose, Size irsd, Size jrsd, Vec cen) {
	Vec cdi = pose.residue(irsd).xyz("CD");
	Vec cdj = pose.residue(jrsd).xyz("CD");
	Vec cgi = pose.residue(irsd).xyz("CG");
	Vec cgj = pose.residue(jrsd).xyz("CG");
	Vec odi = pose.residue(irsd).xyz("OE1");
	Vec odj = pose.residue(jrsd).xyz("OE1");
	Vec a = (cdi-cen).cross(cdj-cen);
	Real angi = dihedral_degrees( odi,cdi,cgi,cen+a );
	Real angj = dihedral_degrees( odj,cdj,cgj,cen+a );
	// pose.dump_pdb("test1.pdb");
	pose.set_chi( 3,irsd, pose.chi(3,irsd)-angi );
	pose.set_chi( 3,jrsd, pose.chi(3,jrsd)-angj );
	// pose.dump_pdb("test2.pdb");
	// utility_exit_with_message("debug");
}


void run_m8() {
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
	string fname = utility::file_basename(infile)+"_glu_m8_matches.pdb";
	if ( option[willmatch::splitwork].user() ) {
		fname = utility::file_basename(infile)+"_glu_m8_matches_part"+string_of(option[willmatch::splitwork]()[1])+".pdb";
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
		vector1<Vec> foundcen;
		vector1<Vec> foundori;
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
						if ( angle < 90.0-MXAGMTL || 90.0+MXAGMTL < angle ) {
							// TR << "angle fail" << angle << " " << orii << " " << orij << std::endl;
							continue;
						}
						pose = fa_pose;
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


						align_carboxyl_m8(pose,irsd,jrsd,cen);
						//now clash check OE1/2
						for ( Size iatm = pose.residue(irsd).nheavyatoms()-1; iatm <= pose.residue(irsd).nheavyatoms(); ++iatm ) {
							if ( !clashcheck.clash_check(pose.residue(irsd).xyz(iatm),irsd) ) clash=true;
						}
						for ( Size jatm = pose.residue(jrsd).nheavyatoms()-1; jatm <= pose.residue(jrsd).nheavyatoms(); ++jatm ) {
							if ( !clashcheck.clash_check(pose.residue(jrsd).xyz(jatm),jrsd) ) clash=true;
						}
						if ( clash ) continue;

						for ( Size iatm = pose.residue(irsd).nheavyatoms()-1; iatm <= pose.residue(irsd).nheavyatoms(); ++iatm ) {
							for ( Size jatm = pose.residue(jrsd).nheavyatoms()-1; jatm <= pose.residue(jrsd).nheavyatoms(); ++jatm ) {
								if ( pose.residue(jrsd).xyz(jatm).distance_squared(pose.residue(irsd).xyz(iatm)) < 7.0 ) clash=true;
							}
						}
						if ( clash ) continue;


						// is cross better than what was in HIS/HIS?
						Vec ori = (orii.cross(orij)).normalized();
						bool overlap = false;
						for ( Size i = 1; i <= foundori.size(); ++i ) {
							if ( cen.distance_squared(foundcen[i]) < option[willmatch::match_overlap_dis]() &&
									ori .dot(foundori [i])            > MATCH_OVERLAP_DOT                      ) {
								// TR << "overlap!" << std::endl;
								overlap = true;
								break;
							}
						}
						if ( overlap ) continue;
						count++;
						foundcen.push_back(cen);
						foundori.push_back(ori);

						string geom = "EE4";
						string tag = utility::file_basename(infile)+"_"+geom+"_"+string_of(irsd)+"_"+string_of(CHI1[ich1])+"_"+string_of(CHI2[ich2])+"_"+string_of(jrsd)+"_"+string_of(CHI1[jch1])+"_"+string_of(CHI2[jch2]);
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

						Vec viz;
						viz = cen+2.0*mori      ; out<<"HETATM"<<I(5,9997)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"B"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						viz = cen+4.0*mori      ; out<<"HETATM"<<I(5,9998)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"C"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						viz = cen               ; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' ' << " ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

						out << "ENDMDL" << std::endl;

						Vec axis1 = ((cen-pose.residue(irsd).xyz("CD"))  +   (cen-pose.residue(jrsd).xyz("CD"))).normalized();
						Vec axis2 = ((cen-pose.residue(irsd).xyz("CD")).cross(cen-pose.residue(jrsd).xyz("CD"))).normalized();
						axis1 = projperp(axis2,axis1).normalized();

						for ( Size ii = 1; ii <= 4; ii++ ) {
							Pose pose2 = pose;
							rot_pose(pose2,axis2,180,cen);

							for ( Size ir = 1; ir <= pose.size(); ++ir ) {
								for ( Size ia = 1; ia <= pose.residue(ir).nheavyatoms()-2; ++ia ) {
									for ( Size jr = 1; jr <= pose.size(); ++jr ) {
										for ( Size ja = 1; ja <= pose.residue(jr).nheavyatoms()-2; ++ja ) {
											if ( pose.residue(ir).xyz(ia).distance_squared(pose2.residue(jr).xyz(ja)) < 9.0 ) clash=true;
										}
									}
								}
							}
							if ( clash ) continue;

							utility::io::ozstream out(tag+"_homodimer_"+string_of(ii)+".pdb");
							pose.dump_pdb(out);
							for ( Size jj = 1; jj <= pose2.size(); ++jj ) {
								core::conformation::Residue r(pose2.residue(jj));
								r.chain(2);
								pose2.replace_residue(jj,r,false);
							}
							pose2.dump_pdb(out);
							viz = cen;         out<<"HETATM"<<I(5,9999)<<' '<<"CA  "<<' ' << " CA"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = cen+2*axis1; out<<"HETATM"<<I(5,9999)<<' '<<"CA  "<<' ' << " CA"<<' '<<"Y"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = cen+2*axis2; out<<"HETATM"<<I(5,9999)<<' '<<"CA  "<<' ' << " CA"<<' '<<"X"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							out.close();
							axis2 = rotation_matrix_degrees(axis1,45.0)*axis2;
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

		run_m8();


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


