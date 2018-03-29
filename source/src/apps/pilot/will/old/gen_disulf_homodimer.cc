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
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <utility/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <protocols/symmetry/SymDockingInitialPerturbation.hh>
#include <sstream>
#include <utility/io/izstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <protocols/moves/MoverStatistics.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>

#include <numeric/NumericTraits.hh>

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static basic::Tracer TR( "genmatch" );


inline Real sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
	if ( sqdist > 9.*9. ) {
		return 0.0;
	} else if ( sqdist<6.*6. ) {
		return 1.0;
	} else {
		Real dist=sqrt( sqdist );
		return sqr(1.0  - sqr( (dist - 6.) / (9. - 6.) ) );
	}
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
	PackerTaskOP task=TaskFactory::create_packer_task(pose);
	task->initialize_extra_rotamer_flags_from_command_line();
	for ( Size i=1; i<=nres; ++i ) {
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
	for ( Size ir=1; ir<=pose.size(); ++ir ) {
		atom_map.set(AtomID(2,ir) , true );
		atom_map.set(AtomID(3,ir) , true );
		atom_map.set(AtomID(5,ir) , true );
	}
	core::id::AtomID_Map<Real> atom_sasa; utility::vector1<Real> sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, sasa, 2.3, false, atom_map );
	for ( Size i=1; i<=sasa.size(); ++i ) if ( atom_sasa.n_atom(i) > 4 ) sasa[i]=atom_sasa[AtomID(5,i)];

	using namespace core::pack::task;
	PackerTaskOP task=TaskFactory::create_packer_task(pose);
	vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys]=false;
	aas[core::chemical::aa_his]=false;
	// aas[core::chemical::aa_met]=false;
	aas[core::chemical::aa_pro]=false;
	aas[core::chemical::aa_gly]=false;
	aas[core::chemical::aa_gly]=false;
	if ( option[basic::options::OptionKeys::willmatch::exclude_ala]() ) aas[core::chemical::aa_ala]=false;
	if ( option[basic::options::OptionKeys::smhybrid::design_hydrophobic]() ) {
		aas[core::chemical::aa_ser]=false;
		aas[core::chemical::aa_thr]=false;
		aas[core::chemical::aa_asp]=false;
		aas[core::chemical::aa_glu]=false;
		aas[core::chemical::aa_lys]=false;
		aas[core::chemical::aa_arg]=false;
		aas[core::chemical::aa_asn]=false;
		aas[core::chemical::aa_gln]=false;
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
		for ( Size i=1; i<=nres; ++i ) {
			if ( sasa.size() >= i && sasa[i] > 15.0 ) continue;
			AtomID aid(5,i);
			if ( pose.residue(i).nheavyatoms()<5 ) aid.atomno()=2;
			for ( Size j=nres+1; j<=3*nres; ++j ) {
				AtomID aid2(5,j);
				if ( pose.residue(j).nheavyatoms()<5 ) aid.atomno()=2;
				if ( pose.xyz(aid).distance_squared(pose.xyz(aid2))<49 ) {
					interface.push_back(i);
				}
			}
		}
		for ( Size i=1; i<=nres; ++i ) {
			Vec xyz=pose.xyz(AtomID(5,i));
			xyz.z()=0;
			if ( xyz.length()<8.0 ) interface.push_back(i);
		}
	}
	for ( Size i=1; i<=nres; ++i ) {
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
	for ( Size i=nres+1; i<=pose.size(); ++i ) {
		task->nonconst_residue_task(i).prevent_repacking();
	}
	TR << *task << std::endl;

	protocols::minimization_packing::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}


void minimize(Pose & pose, Size nres, Size , ScoreFunctionOP sf, int bb=0) {
	core::kinematics::MoveMapOP movemap=new core::kinematics::MoveMap;
	// core::pose::symmetry::make_symmetric_movemap(pose,*movemap);
	movemap->set_chi(true);
	movemap->set_bb(false);
	movemap->set_jump(false);
	if ( bb==1 ) for ( Size i=1; i<=nres; ++i ) if ( pose.secstruct(i)=='L' ) movemap->set_bb(i,true);
	if ( bb>=2 ) for ( Size i=1; i<=nres; ++i ) movemap->set_bb(i,true);
	// movemap->set_chi(bpyres,false);

	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );

	protocols::minimization_packing::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );

	m.apply(pose);

}

struct Hit {
	Real score,ch1,ch2;
	Size icys;
	bool idsf;
	Mat rotn;
	Vec rcen;
	Hit(Real sc, Size ic, Real i1, Real i2, bool id, Mat r, Vec c) : score(sc),ch1(i1),ch2(i2),icys(ic),idsf(id),rotn(r),rcen(c) {}
};
struct HitCmp {
	bool operator() (Hit const & i, Hit const & j) { return i.score < j.score; }
};


core::pack::rotamer_set::RotamerSetOP get_rotset(Pose & pose, Size icys) {
	core::pack::rotamer_set::RotamerSetOP rotset;
	core::scoring::ScoreFunction dummy_sfxn;
	dummy_sfxn( pose );
	core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
	dummy_task->initialize_from_command_line();
	dummy_task->nonconst_residue_task( icys ).and_extrachi_cutoff(0);
	dummy_task->nonconst_residue_task( icys ).restrict_to_repacking();
	dummy_task->nonconst_residue_task( icys ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
	dummy_task->nonconst_residue_task( icys ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
	utility::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
	core::pack::rotamer_set::RotamerSetFactory rsf;
	rotset = rsf.create_rotamer_set( pose.residue( icys ) );
	rotset->set_resid( icys );
	rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
	return rotset;
}

void run(std::string fname) {
	using namespace std;
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace chemical;
	using namespace id;
	using namespace pose;
	using namespace scoring;

	// Size ANGLE_INCR=30;

	const Real PI = numeric::NumericTraits<Real>::pi();
	// setup stuff
	ResidueTypeSetCAP frs=ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	ResidueTypeSetCAP crs=ChemicalManager::get_instance()->residue_type_set( CENTROID );
	ScoreFunctionOP sfstd=get_score_function();
	ScoreFunctionOP sfcen=ScoreFunctionFactory::create_score_function("score3");


	// read pose info
	Pose cenp,cysp,natp,cys,ile,homo;
	make_pose_from_sequence(cys,"C",core::chemical::FA_STANDARD,false);
	if ( cys.residue(1).is_lower_terminus() ) remove_lower_terminus_type_from_pose_residue(cenp,1);
	if ( cys.residue(1).is_upper_terminus() ) remove_upper_terminus_type_from_pose_residue(cenp,1);
	make_pose_from_sequence(ile,"I",core::chemical::CENTROID   ,false);
	if ( ile.residue(1).is_lower_terminus() ) remove_lower_terminus_type_from_pose_residue(cenp,1);
	if ( ile.residue(1).is_upper_terminus() ) remove_upper_terminus_type_from_pose_residue(cenp,1);
	pose_from_file(cenp,*crs,fname, core::import_pose::PDB_file);
	pose_from_file(natp,*frs,fname, core::import_pose::PDB_file);
	Size nres=cenp.size();
	for ( Size ir=1; ir<=nres; ++ir ) {
		if ( cenp.residue(ir).is_lower_terminus() ) remove_lower_terminus_type_from_pose_residue(cenp,ir);
		if ( natp.residue(ir).is_lower_terminus() ) remove_lower_terminus_type_from_pose_residue(natp,ir);
		if ( cenp.residue(ir).is_upper_terminus() ) remove_upper_terminus_type_from_pose_residue(cenp,ir);
		if ( natp.residue(ir).is_upper_terminus() ) remove_upper_terminus_type_from_pose_residue(natp,ir);
	}
	core::scoring::dssp::Dssp dssp(natp);
	dssp.insert_ss_into_pose(natp);
	dssp.insert_ss_into_pose(cenp);

	cysp=natp;
	homo=cenp;
	homo.append_residue_by_jump(cenp.residue(1),1);
	for ( Size ir=2; ir<=nres; ++ir ) homo.append_residue_by_bond(cenp.residue(ir));
	for ( Size ir=1; ir<=2*nres; ++ir ) homo.replace_residue(ir,ile.residue(1),true);

	for ( Size ir=1; ir<=nres; ++ir ) {
		cysp.replace_residue(ir,cys.residue(1),true);
		// set HG as other SG for disulf
		cysp.set_dof( DOF_ID(AtomID(11,ir),PHI  ), PI);
		cysp.set_dof( DOF_ID(AtomID(11,ir),THETA), 1.368997);
		cysp.set_dof( DOF_ID(AtomID(11,ir),D    ), 2.020   );
	}
	ImplicitFastClashCheck ifc(cysp,3.0);


	// get optional res list
	vector1<Size> ifres; {
		vector1<Size> tmp=read_res_list(option[willmatch::exclude_res1]());
		for ( Size i=1; i<=nres; ++i ) if ( find(tmp.begin(),tmp.end(),i)==tmp.end() ) ifres.push_back(i);
	}

	vector1<Hit> hits;

	for ( Size icys=3; icys<=nres-2; ++icys ) {
		if ( std::find(ifres.begin(),ifres.end(),icys)==ifres.end() ) continue;
		if ( cysp.secstruct(icys)!='H' ) continue;
		TR << "icys: " << icys << " " << hits.size() << std::endl;
		pack::rotamer_set::RotamerSetOP rotset = get_rotset(cysp,icys);
		for ( Size crot = 1; crot <= rotset->num_rotamers(); ++crot ) {
			cysp.set_chi(1,icys,rotset->rotamer(crot)->chi(1));
			cysp.set_chi(2,icys,rotset->rotamer(crot)->chi(2));
			if ( !ifc.clash_check( cysp.xyz(AtomID(6 ,icys)), icys ) ) continue;
			if ( !ifc.clash_check( cysp.xyz(AtomID(11,icys)), icys ) ) continue;
			for ( Size idsf=0; idsf<2; ++idsf ) {
				Vec axs=projperp( cysp.xyz(AtomID(11,icys))-cysp.xyz(AtomID(6,icys)), cysp.xyz(AtomID(6,icys))-cysp.xyz(AtomID(5,icys))).normalized();
				axs=rotation_matrix_degrees( cysp.xyz(AtomID(11,icys))-cysp.xyz(AtomID(6,icys)), (idsf)?45.0:-45.0 ) * axs;
				Mat rotn=rotation_matrix_degrees(axs,180.0);
				Vec rcen=(cysp.xyz(AtomID(6,icys))+cysp.xyz(AtomID(11,icys)))/2.0; // cen of disulf bond
				bool clash=false;
				for ( Size jr=1; jr<=nres; ++jr ) {
					for ( Size ja=1; ja<=5; ja++ ) {
						Vec const xyz(rotn*(cysp.xyz(AtomID(ja,jr))-rcen)+rcen);
						if ( !ifc.clash_check(xyz) ) clash=true;
					}
					if ( clash ) break;
				}
				if ( clash ) continue;
				for ( Size ir = 1; ir <= nres; ++ir ) {
					for ( Size ia = 1; ia <= homo.residue_type(ir).natoms(); ++ia ) {
						homo.set_xyz( AtomID(ia,ir+nres), rotn * (homo.xyz(AtomID(ia,ir))-rcen) + rcen );
					}
				}
				hits.push_back(Hit(sfcen->score(homo),icys,rotset->rotamer(crot)->chi(1),rotset->rotamer(crot)->chi(2),idsf,rotn,rcen));
			} // idsf
		} // crot
	} // icys

	std::sort(hits.begin(),hits.end(),HitCmp());

	homo=natp;
	homo.append_residue_by_jump(natp.residue(1),1);
	for ( Size ir=2; ir<=nres; ++ir ) homo.append_residue_by_bond(natp.residue(ir));

	for ( Size ih = 1; ih <= 20; ++ih ) {
		Hit const & h(hits[ih]);
		for ( Size ir = 1; ir <= nres; ++ir ) {
			for ( Size ia = 1; ia <= homo.residue_type(ir).natoms(); ++ia ) {
				homo.set_xyz( AtomID(ia,ir+nres), h.rotn * (homo.xyz(AtomID(ia,ir))-h.rcen) + h.rcen );
			}
		}
		homo.replace_residue(h.icys     ,cys.residue(1),true);
		homo.set_chi(1,h.icys     ,h.ch1);
		homo.set_chi(2,h.icys     ,h.ch2);
		homo.replace_residue(h.icys+nres,cys.residue(1),true);
		homo.set_chi(1,h.icys+nres,h.ch1);
		homo.set_chi(2,h.icys+nres,h.ch2);
		homo.dump_pdb("hit"+lzs(ih,4)+".pdb");
	}


}


int main (int argc, char *argv[]) {

	try {


		devel::init(argc,argv);

		using basic::options::option;
		using namespace basic::options::OptionKeys;

		run(option[in::file::s]()[1]);

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


