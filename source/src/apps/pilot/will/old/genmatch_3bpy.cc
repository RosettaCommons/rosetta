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

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
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
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static basic::Tracer TR( "gentetra" );
static core::io::silent::SilentFileData sfd;


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

Real iface_check_c3(Pose & pose, Size nres) {
	Mat const R1 = rotation_matrix_degrees(Vec(0,0,1), 120.0);
	Mat const R2 = rotation_matrix_degrees(Vec(0,0,1),-120.0);
	Real num = 0;
	for ( Size ir = 1; ir <= nres; ir++ ) {
		if ( !pose.residue(ir).is_protein() ) continue;
		for ( Size jr = 1; jr <= nres; jr++ ) {
			if ( !pose.residue(jr).is_protein() ) continue;
			num += sigmoidish_neighbor(pose.residue(ir).xyz(5).distance_squared( R1 * pose.residue(jr).xyz(5) ));
			num += sigmoidish_neighbor(pose.residue(ir).xyz(5).distance_squared( R2 * pose.residue(jr).xyz(5) ));
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
			if ( !pose.residue(i).is_protein() ) continue;
			if ( pose.residue(i).nheavyatoms() < 5 ) aid.atomno() = 2;
			for ( Size j = nres+1; j <= 3*nres; ++j ) {
				AtomID aid2(5,j);
				if ( !pose.residue(j).is_protein() ) continue;
				if ( pose.residue(j).nheavyatoms() < 5 ) aid.atomno() = 2;
				if ( pose.xyz(aid).distance_squared(pose.xyz(aid2)) < 49.0 ) {
					interface.push_back(i);
				}
			}
		}
		for ( Size i = 1; i <= nres; ++i ) {
			AtomID aid(5,i);
			if ( !pose.residue(i).is_protein() ) continue;
			if ( pose.residue(i).nheavyatoms() < 5 ) aid.atomno() = 2;
			Vec xyz = pose.xyz(aid);
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
	sf->score(pose);

	Real worig = sf->get_weight(core::scoring::res_type_constraint);
	if ( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);

	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}

vector1<Vec> line_cone_intersection(Vec p, Vec d, Vec v, Vec a, Real t) {
	using namespace numeric;
	vector1<Vec> sol;
	t = conversions::radians(t);
	Mat M = outer_product(a,a) - cos(t)*cos(t)*Mat::identity();
	Vec D = p-v;
	Real c2 =   d.dot(M*d);
	Real c1 = 2*d.dot(M*D);
	Real c0 =   D.dot(M*D);
	Real disc = c1*c1 - 4*c0*c2;
	if ( disc == 0 ) sol.push_back( p + (-c1)/(2.0*c2)*d );
	else if ( disc > 0 ) {
		disc = sqrt(disc);
		sol.push_back(p+(-c1+disc)/(2.0*c2)*d);
		sol.push_back(p+(-c1-disc)/(2.0*c2)*d);
	}
	return sol;
}


vector1<std::pair<Vec,Vec> > intersecting_disulfide_axes(Vec CB, Vec SG, Vec symmaxis, Vec symmcen = Vec(0,0,0)) {
	vector1<std::pair<Vec,Vec> > sol;
	Vec p = symmcen;
	Vec d = symmaxis.normalized();
	Vec v = SG;
	Vec a = (SG-CB).normalized();
	Real t = 45.0; // dihedral should be 90 = 2*45
	vector1<Vec> X = line_cone_intersection(p,d,v,a,t);
	for ( Size i = 1; i <= X.size(); ++i ) {
		Vec x = X[i];
		Vec o = projperp(a,x-v);
		Real L = o.length();
		o = o.normalized();
		Real ang = 90.0 - numeric::conversions::degrees( atan(1.0/L) );
		Vec o1 = rotation_matrix_degrees(a, ang) * o;
		Vec o2 = rotation_matrix_degrees(a,-ang) * o;
		sol.push_back(std::pair<Vec,Vec>(x,v+o1));
		sol.push_back(std::pair<Vec,Vec>(x,v+o2));
	}
	return sol;
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


void run() {

	using namespace basic::options::OptionKeys;

	core::chemical::ResidueTypeSetCAP  rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Pose bpy1,bpy2;
	core::import_pose::pose_from_file(bpy1,*rs,"input/py_aa_align_trimer_in.pdb", core::import_pose::PDB_file);
	core::import_pose::pose_from_file(bpy2,*rs,"input/py_aa_align_trimer_in.pdb", core::import_pose::PDB_file);
	for ( int i = 1; i <= 2; ++i ) {
		Pose & bpy( i==1 ? bpy1 : bpy2 );
		core::pose::remove_lower_terminus_type_from_pose_residue(bpy,1);
		core::pose::remove_upper_terminus_type_from_pose_residue(bpy,1);
		bpy.append_residue_by_jump( *core::conformation::ResidueFactory::create_residue(rs->name_map("VRT")),1,"NN1","ORIG");
		core::kinematics::FoldTree ft = bpy.fold_tree();
		ft.reorder(2);
		bpy.fold_tree(ft);
	}

	for ( Size ifile = 1; ifile <= option[in::file::s]().size(); ++ifile ) {
		string infile = utility::file_basename(option[in::file::s]()[ifile]);
		Pose init;
		core::import_pose::pose_from_file(init,*rs,option[in::file::s]()[ifile], core::import_pose::PDB_file);


		Size nres = init.size();
		// ScoreFunctionOP sf = core::scoring::get_score_function();
		ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction(core::scoring::get_score_function());

		//    vector1<int> fres = option[willmatch::forbid_residues]();
		vector1<numeric::xyzTriple<Real> > chis;
		for ( Size i = 2; i <= init.size()-1; ++i ) {
			//if( std::find(fres.begin(),fres.end(),i) != fres.end() ) continue;
			for ( Real chi1 = -180.0; chi1 < 180.0; chi1+=option[willmatch::chi1_increment]() ) {
				for ( Real chi2 = -180.0; chi2 < 180.0; chi2+=option[willmatch::chi2_increment]() ) {
					chis.push_back(numeric::xyzTriple<Real>(i,chi1,chi2));
				}
			}
		}
		std::pair<Size,Size> split = makesplitwork_3bpy(chis.size());
		Pose pose = init;

		Mat const R1 = rotation_matrix_degrees(Vec(0,0,1), 120.0);
		Mat const R2 = rotation_matrix_degrees(Vec(0,0,1),-120.0);


		Size lastrsd = 0;
		for ( Size ichi = max(1lu,split.first); ichi <= min(chis.size(),split.second); ichi++ ) {
			Size irsd = chis[ichi].x();
			if ( pose.residue(irsd).name3()=="GLY" || pose.residue(irsd).name3()=="PRO" || pose.residue(irsd).name3()=="DPR" ) continue;
			Real chi1 = chis[ichi].y();
			Real chi2 = chis[ichi].z();
			// for helix BPY
			// if( chi1 < -80 || -40 < chi1 && chi1 < 160 || 200 < chi1 ) continue;
			// if( chi2 < 70 || 110 < chi2 ) continue;
			bpy.set_phi  (1,pose.phi  (irsd));
			bpy.set_psi  (1,pose.psi  (irsd));
			bpy.set_omega(1,pose.omega(irsd));
			bool newrsd = lastrsd!=irsd;
			lastrsd = irsd;
			if ( newrsd ) {
				TR << "res " << irsd << std::endl;
				pose = init;
				pose.replace_residue(irsd,bpy.residue(1),true);
				core::id::AtomID_Map< core::id::AtomID > atom_map;
				core::pose::initialize_atomid_map( atom_map, pose, core::id::AtomID::BOGUS_ATOM_ID() ); // maps every atomid to bogus atom
				for ( Size i = 6; i <= bpy.residue(1).nheavyatoms(); ++i ) atom_map[AtomID(i,irsd)] = AtomID(i,1);
				core::scoring::superimpose_pose( pose, bpy, atom_map );
				pose.append_residue_by_jump( *core::conformation::ResidueFactory::create_residue(rs->name_map("VRT")),irsd,"NN1","ORIG");
				core::kinematics::FoldTree ft = pose.fold_tree();
				ft.reorder(pose.size());
				pose.fold_tree(ft);
			}
			pose.set_chi(1,irsd,chi1);
			pose.set_chi(2,irsd,chi2);

			// core::id::AtomID_Map< core::id::AtomID > atom_map;
			// core::pose::initialize_atomid_map( atom_map, pose, core::id::AtomID::BOGUS_ATOM_ID() ); // maps every atomid to bogus atom
			// atom_map[AtomID(2,irsd)] = AtomID(2,1);
			// atom_map[AtomID(3,irsd)] = AtomID(3,1);
			// atom_map[AtomID(4,irsd)] = AtomID(4,1);
			// atom_map[AtomID(5,irsd)] = AtomID(5,1);
			// core::scoring::superimpose_pose( pose, bpy, atom_map );

			// Pose tmp = pose;
			// core::pose::symmetry::make_symmetric_pose( tmp );
			// tmp.dump_pdb("test"+string_of(ichi)+".pdb");
			// if(ichi==10) utility_exit_with_message("debug");

			ImplicitFastClashCheck ifc(pose,2.0);
			for ( Size ia = 6; ia <= pose.residue(irsd).nheavyatoms(); ++ia ) if ( !ifc.clash_check(    pose.xyz(AtomID(ia,irsd)), irsd ) ) goto cont0;
			//for(Size ia = 1; ia <= pose.residue(irsd).nheavyatoms(); ++ia) if( !ifc.clash_check( R1*pose.xyz(AtomID(ia,irsd)), irsd ) ) goto cont0;
			//for(Size ia = 1; ia <= pose.residue(irsd).nheavyatoms(); ++ia) if( !ifc.clash_check( R2*pose.xyz(AtomID(ia,irsd)), irsd ) ) goto cont0;
			goto done0; cont0:

			pose.dump_pdb("test.pdb");
			utility_exit_with_message("bpy clash 2.0");

			continue; done0:

			;
			//TR << "no bpy clash" << std::endl;

			// TR << "iface check start" << std::endl;
			// Real iface = iface_check_c3(pose,nres);
			// TR << "iface check DONE" << std::endl;
			// if( iface < 50 ) continue;

			// // clash check
			// for(Size ir = 1; ir <= nres; ++ir) {
			//   for(Size jr = 1; jr <= nres; ++jr) {
			//     for(Size ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia) {
			//       for(Size ja = 1; ja <= pose.residue(jr).nheavyatoms(); ++ja) {
			//         if( pose.xyz(AtomID(ia,ir)).distance_squared( R1*pose.xyz(AtomID(ja,jr)) ) < 9.0) if(pose.residue(ir).name3()!="BPY" || pose.residue(jr).name3()!="BPY") goto cont1;
			//         if( pose.xyz(AtomID(ia,ir)).distance_squared( R2*pose.xyz(AtomID(ja,jr)) ) < 9.0) if(pose.residue(ir).name3()!="BPY" || pose.residue(jr).name3()!="BPY") goto cont1;
			//       }
			//     }
			//   }
			// }
			// goto done1; cont1: continue; done1:

			// string fname = infile+"_"+lead_zero_string_of(irsd,3)+"_"+lead_zero_string_of((Size)(chi1+180.0),3)+"_"+lead_zero_string_of((Size)(chi2+180.0),3)+".pdb";
			// TR << "found nonclashing trimer: " << fname << std::endl;
			// // Pose symm = pose;
			// // core::pose::symmetry::make_symmetric_pose( symm );

			// // design(symm,nres+1,sf);

			// // sf->score(symm);
			// // core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
			// // ss_out->fill_struct(symm,fname);
			// // sfd.write_silent_struct( *ss_out, option[ basic::options::OptionKeys::out::file::silent ]() );
			// // symm.dump_pdb(fname);

		}
	}
}

int main (int argc, char *argv[]) {

	try {


		// Vec p(-6.456746, 5.922204, -0.982538);
		// Vec d(0.393718,  0.677101,  0.621707);
		// Vec v(0.000000,  0.000000,  0.000000);
		// Vec a(0.998233, -0.003844, -0.059301);
		// Real t = 45;
		// vector1<Vec> X = line_cone_intersection(p,d,v,a,t);
		// for(Size i = 1; i <= X.size(); ++i) {
		//  TR << "X " << i << " " << X[i] << std::endl;
		// }
		// utility_exit_with_message("debug line_cone_intersection");

		devel::init(argc,argv);
		run();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


