// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#define sym_Clo 5
#define sym_Chi 5
#define CONTACT_D2 20.25
#define CONTACT_TH 5
#define MIN_HELEX_RES 20
#define MAX_CYS_RES 3
#define MAX_NRES 200
#define MIN_NBR_COUNT 7
#define MAX_ANG_METAL 15.0
#define MAX_DIS_METAL 0.5
#define MIN_DIS_METAL_SEP 10.0

#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <utility/graph/Graph.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>b
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/will_util.ihh>

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

using core::id::AtomID;
using basic::options::option;
using core::pose::Pose;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::rotation_matrix_degrees;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using std::pair;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_file;
using core::kinematics::Stub;

static basic::Tracer tr( "hczndimer" );
static core::io::silent::SilentFileData sfd;

vector1<pair<Real,Real> >
makerots(Pose const & p, Size ir, Pose const & rsd) {
	Pose pose(p);
	pose.replace_residue(ir,rsd.residue(1),true);
	core::scoring::ScoreFunction dummy_sfxn;
	dummy_sfxn.set_weight(core::scoring::fa_rep,0.44);
	dummy_sfxn( pose );
	core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
	dummy_task->nonconst_residue_task(ir).restrict_to_repacking();
	dummy_task->nonconst_residue_task(ir).or_fix_his_tautomer(false);
	dummy_task->nonconst_residue_task(ir).and_extrachi_cutoff(1);
	dummy_task->nonconst_residue_task(ir).or_ex1(true);
	dummy_task->nonconst_residue_task(ir).or_ex2(true);
	// NO_EXtrA_CHI_SAMPLES = 0,      //0
	// EX_ONE_STDDEV,                 //1
	// EX_ONE_HALF_STEP_STDDEV,       //2
	// EX_TWO_FULL_STEP_STDDEVS,      //3
	// EX_TWO_HALF_STEP_STDDEVS,      //4
	// EX_FOUR_HALF_STEP_STDDEVS,     //5
	// EX_THREE_THIRD_STEP_STDDEVS,   //6
	// EX_SIX_QUARTER_STEP_STDDEVS,   //7
	dummy_task->nonconst_residue_task(ir).or_ex1_sample_level(core::pack::task::EX_TWO_FULL_STEP_STDDEVS);
	dummy_task->nonconst_residue_task(ir).or_ex2_sample_level(core::pack::task::EX_TWO_FULL_STEP_STDDEVS);
	utility::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
	core::pack::rotamer_set::RotamerSetFactory rsf;
	core::pack::rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( pose.residue(ir) );
	rotset->set_resid(ir);
	rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
	vector1<pair<Real,Real> > rots;
	for ( Size irot = 1; irot <= rotset->num_rotamers(); ++irot ) {
		Real r1 = rotset->rotamer(irot)->chi(1);
		Real r2 = rotset->rotamer(irot)->chi(2);
		bool seenit = false;
		for ( vector1<pair<Real,Real> >::iterator i = rots.begin(); i != rots.end(); ++i ) {
			if ( fabs(r1-i->first) < 5.0 && fabs(r2-i->second) < 5.0 ) seenit=true;
		}
		if ( !seenit ) rots.push_back(pair<Real,Real>(r1,r2));
	}
	//tr << rotset->num_rotamers() << " " << rots.size() << std::endl;
	return rots;
}

bool clashcheck(Pose const & p, Vec v) {
	for ( Size ir = 1; ir <= p.size(); ++ir ) {
		if ( p.xyz(AtomID(1,ir)).distance_squared(v) < 9.0 ) return false;
		if ( p.xyz(AtomID(2,ir)).distance_squared(v) < 9.0 ) return false;
		if ( p.xyz(AtomID(3,ir)).distance_squared(v) < 9.0 ) return false;
		if ( p.xyz(AtomID(4,ir)).distance_squared(v) < 9.0 ) return false;
		if ( p.residue(ir).name3()=="GLY" ) continue;
		if ( p.xyz(AtomID(5,ir)).distance_squared(v) < 9.0 ) return false;
	}
	return true;
}
bool clashcheckhalf(Pose const & p, Vec v) {
	for ( Size ir = 1; ir <= p.size(); ++ir ) {
		if ( p.xyz(AtomID(1,ir)).distance_squared(v) < 4.0 ) return false;
		if ( p.xyz(AtomID(2,ir)).distance_squared(v) < 4.0 ) return false;
		if ( p.xyz(AtomID(3,ir)).distance_squared(v) < 4.0 ) return false;
		if ( p.xyz(AtomID(4,ir)).distance_squared(v) < 4.0 ) return false;
		if ( p.residue(ir).name3()=="GLY" ) continue;
		if ( p.xyz(AtomID(5,ir)).distance_squared(v) < 4.0 ) return false;
	}
	return true;
}

enum RTYPE {
	CYS = 1,
	HIS1,
	HIS2,
	ASP1,
	ASP2,
	ASP3,
	NRTYPES = ASP3
};

struct Hit {
	Vec cen,axs,ori;
	Size ir,jr,itype,jtype;
	Real ich1,ich2,jch1,jch2;
	Hit(){}
	Hit(Vec c, Vec a, Vec o, Size irs, Size jrs, Size _itype, Size _jtype, Real ic1, Real ic2, Real jc1, Real jc2)
	: cen(c),axs(a.normalized()),ori(o.normalized()),ir(irs),jr(jrs),itype(_itype),jtype(_jtype),ich1(ic1),ich2(ic2),jch1(jc1),jch2(jc2) {}
};
std::ostream & operator<<(std::ostream & o, Hit const & h) {
	o << h.cen  << " " << h.axs  << " " << h.ori   << " ";
	o << h.ir   << " " << h.jr   << " " << h.itype << " " << h.jtype << " ";
	o << h.ich1 << " " << h.ich2 << " " << h.jch1  << " " << h.jch2;
	return o;
}
std::istream & operator>>(std::istream & i, Hit & h) {
	i >> h.cen  >> h.axs  >> h.ori;
	i >> h.ir   >> h.jr   >> h.itype >> h.jtype;
	i >> h.ich1 >> h.ich2 >> h.jch1  >> h.jch2;
	return i;
}

void refine(Pose & pose, ScoreFunctionOP sf, Size r1, Size r2, Size r3, Size r4 ) {
	core::conformation::symmetry::SymmetryInfoCOP syminfo = core::pose::symmetry::symmetry_info(pose);
	Size N = syminfo->num_independent_residues();

	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
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
	utility_exit_with_message("NOT DONE CODING!!!!!!!");
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

	protocols::minimization_packing::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);

	// cleanup 2
	pose.remove_constraints( res_cst );
	sf->set_weight(core::scoring::res_type_constraint,worig);


	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(true);
	movemap->set_bb(true);
	movemap->set_chi(true);
	protocols::minimization_packing::symmetry::SymMinMover( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);


}


Pose make_single_res_pose(string rt) {
	Pose tmp;
	make_pose_from_sequence(tmp,rt,core::chemical::FA_STANDARD,false);
	remove_lower_terminus_type_from_pose_residue(tmp,1);
	remove_upper_terminus_type_from_pose_residue(tmp,1);
	return tmp;
}
Real sqr(Real x) { return x*x; }

void dock(Pose init, std::string const & fn) {
	ScoreFunctionOP sf = core::scoring::get_score_function();
	ScoreFunctionOP sfnosym = new core::scoring::ScoreFunction(*sf);
	if ( sf->get_weight(core::scoring::atom_pair_constraint) < 0.0001 || sf->get_weight(core::scoring::atom_pair_constraint) < 0.0001 ) {
		utility_exit_with_message("no constraint weight on pair if angle!! use -score:weights");
	}
	// ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction;
	// sf->set_weight(core::scoring::atom_pair_constraint,1.0);
	// sf->set_weight(core::scoring::angle_constraint,1.0);
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(false); movemap->set_bb(true); movemap->set_chi(true);
	protocols::minimization_packing::MinMover( movemap, sfnosym, "lbfgs_armijo_nonmonotone", 1e-3, true, false, false ).apply(init);
	/*///////////////////////////////////////////////////////////////////////////////////*/ tr << "make mbcount" << endl; /*//////////////////////*/
	vector1<Size> nbcount(init.size(),0);
	for ( Size ir = 1; ir <= init.size(); ++ir ) {
		for ( Size jr = 1; jr <= init.size(); ++jr ) {
			if ( init.xyz(AtomID(2,ir)).distance_squared(init.xyz(AtomID(2,jr))) < 100.0 ) nbcount[ir]++;
		}
	}
	/*////////////////////////////////////////////////////////////////////////////////////*/ tr << "make poses" << endl; /*//////////////////////*/
	Pose ala = make_single_res_pose("A");
	vector1<Pose> res(NRTYPES);
	res[CYS ] = make_single_res_pose("C[CYS_M]" );
	res[HIS1] = make_single_res_pose("H[HIS_M1]");
	res[HIS2] = make_single_res_pose("H[HIS_M2]");
	res[ASP1] = make_single_res_pose("D[ASP_M1]");
	res[ASP2] = make_single_res_pose("D[ASP_M2]");
	res[ASP3] = make_single_res_pose("D[ASP_M3]");
	vector1<Size> matom(NRTYPES),batom(NRTYPES);
	matom[CYS ] = res[CYS ].residue(1).atom_index("ZN"); batom[CYS ] = res[CYS ].residue(1).atom_index("SG" );
	matom[HIS1] = res[HIS1].residue(1).atom_index("ZN"); batom[HIS1] = res[HIS1].residue(1).atom_index("ND1");
	matom[HIS2] = res[HIS2].residue(1).atom_index("ZN"); batom[HIS2] = res[HIS2].residue(1).atom_index("NE2");
	matom[ASP1] = res[ASP1].residue(1).atom_index("ZN"); batom[ASP1] = res[ASP1].residue(1).atom_index("OD1");
	matom[ASP2] = res[ASP2].residue(1).atom_index("ZN"); batom[ASP2] = res[ASP2].residue(1).atom_index("OD1");
	matom[ASP3] = res[ASP3].residue(1).atom_index("ZN"); batom[ASP3] = res[ASP3].residue(1).atom_index("CG" );
	vector1<Real> discbm(NRTYPES);
	for ( Size i = 1; i <= NRTYPES; ++i ) discbm[i] = res[i].xyz(AtomID(5,1)).distance(res[i].xyz(AtomID(matom[i],1)));
	/*////////////////////////////////////////////////////////////////////////////////*/ tr << "make rotamers" << endl; /*//////////////////////*/
	vector1<vector1<vector1<pair<Real,Real> > > > allrots(NRTYPES);
	Pose p(init);
	//for(Size ir = 1; ir <= p.size(); ++ir) p.replace_residue(ir,ala.residue(1),true);
	// search pairs
	vector1<Hit> hits;

	if ( false ) { /////////////////////////////////////////////////////////////////////////////////////////


		for ( Size ir = 1; ir <= init.size(); ++ir ) {
			allrots[CYS ].push_back(makerots(init,ir,res[CYS ]));
			allrots[HIS1].push_back(makerots(init,ir,res[HIS1]));
			allrots[HIS2] = allrots[HIS1];
			allrots[ASP1].push_back(makerots(init,ir,res[ASP1]));
			allrots[ASP2] = allrots[ASP1];
			allrots[ASP3] = allrots[ASP1];
		}
		/*////////////////////////////////////////////////////////////////////////////////*/ tr << "find pairs" << endl; /*//////////////////////*/
		for ( Size ir = 1; ir <= p.size(); ++ir ) {
			if ( nbcount[ir] < MIN_NBR_COUNT ) continue;
			if ( p.residue(ir).name3()=="GLY"||p.residue(ir).name3()=="PRO" ) continue;
			core::conformation::Residue itmp(p.residue(ir)); // remember replaced res
			tr << ir << " " << hits.size() << endl;
			for ( Size jr = ir+1; jr <= p.size(); ++jr ) {
				if ( nbcount[jr] < MIN_NBR_COUNT ) continue;
				if ( p.xyz(AtomID(5,ir)).distance_squared(p.xyz(AtomID(5,jr))) > 141.0 ) continue;
				if ( p.residue(jr).name3()=="GLY"||p.residue(jr).name3()=="PRO" ) continue;
				core::conformation::Residue jtmp(p.residue(jr)); // remember replaced res
				for ( Size itype = 1; itype <= NRTYPES; ++itype ) {
					if ( p.xyz(AtomID(5,ir)).distance_squared(p.xyz(AtomID(5,jr))) > sqr(discbm[itype]+5.7+0.5) ) continue;
					p.replace_residue(ir,res[itype].residue(1),true);
					vector1<pair<Real,Real> > const & irots( allrots[itype][ir] );
					for ( Size jtype = 1; jtype <= NRTYPES; ++jtype ) {
						if ( p.xyz(AtomID(5,ir)).distance_squared(p.xyz(AtomID(5,jr))) > sqr(discbm[itype]+discbm[jtype]+0.5) ) continue;
						p.replace_residue(jr,res[jtype].residue(1),true);
						vector1<pair<Real,Real> > const & jrots( allrots[jtype][jr] );
						for ( Size irot = 1; irot <= irots.size(); ++irot ) {
							p.set_chi(1,ir,irots[irot].first);
							p.set_chi(2,ir,irots[irot].second);
							Vec const ix = p.xyz(AtomID(matom[itype],ir));
							for ( Size jrot = 1; jrot <= jrots.size(); ++jrot ) {
								p.set_chi(1,jr,jrots[jrot].first);
								p.set_chi(2,jr,jrots[jrot].second);
								Vec const jx = p.xyz(AtomID(matom[jtype],jr));
								if ( ix.distance_squared(jx) > MAX_DIS_METAL*MAX_DIS_METAL ) continue; // dist check
								Vec const m = (ix+jx)/2.0;
								Vec const ib = p.xyz(AtomID(batom[itype],ir));
								Vec const jb = p.xyz(AtomID(batom[jtype],jr));
								if ( fabs(109.471313101-angle_degrees(ib,m,jb)) > MAX_ANG_METAL ) continue; // ang check
								bool clash = false;
								for ( Size ia = 5; ia <= p.residue(ir).nheavyatoms(); ++ia ) {
									if ( p.residue(ir).is_virtual(ia) ) continue;
									for ( Size ja = 5; ja <= p.residue(jr).nheavyatoms(); ++ja ) {
										if ( p.residue(jr).is_virtual(ja) ) continue;
										if ( p.xyz(AtomID(ia,ir)).distance_squared(p.xyz(AtomID(ja,jr))) < 9.0 ) clash=true;
									}
								}
								if ( clash ) continue; // clash check
								Vec a = (((m-ib).normalized()+(m-jb).normalized())/2.0).normalized();
								Vec L3 = rotation_matrix_degrees(a,90.0)*((m+2.0*(m-ib).normalized())-m)+m;
								Vec L4 = rotation_matrix_degrees(a,90.0)*((m+2.0*(m-jb).normalized())-m)+m;
								if ( !clashcheck(p,L3) ) continue;
								if ( !clashcheck(p,L4) ) continue;
								// Pose tmp(ala);
								// tmp.replace_residue(1,p.residue(ir),false);
								// tmp.set_xyz(AtomID(1,1),L3);
								// tmp.dump_pdb("/tmp/work/HIT_"+lzs(ir,3)+"_"+lzs(jr,3)+"_"+lzs(irot,3)+"_"+lzs(jrot,3)+"_H.pdb");
								// tmp.replace_residue(1,p.residue(jr),false);
								// tmp.set_xyz(AtomID(1,1),L4);
								// tmp.dump_pdb("/tmp/work/HIT_"+lzs(ir,3)+"_"+lzs(jr,3)+"_"+lzs(irot,3)+"_"+lzs(jrot,3)+"_C.pdb");
								// tr << "HALFHIT "+lzs(ir,3)+"_"+lzs(jr,3)+"_"+lzs(irot,3)+"_"+lzs(jrot,3) << endl;
								Hit h(m,a,a.cross(L3-m),ir,jr,itype,jtype,p.chi(1,ir),p.chi(2,ir),p.chi(1,jr),p.chi(2,jr));
								hits.push_back(h);
								cout << std::setprecision(10);
								cout << h << endl;
							} // jrot
						} // irot
					} // jtype
				} // itype
				p.replace_residue(jr,jtmp,false);
			} // jr
			p.replace_residue(ir,itmp,false);
		} // ir
	}///////////////////////////////////////////////////////////////////////////////////////////////////////

	tr << "reading from file!!!" << endl;
	std::ifstream in("hits.dat");
	Hit h;
	while ( in >> h ) hits.push_back(h);
	tr << "read " << hits.size() << " hits" << endl;


	Pose & p2(p);
	for ( Size ih = 230; ih <= hits.size(); ++ih ) {
		Hit & hi(hits[ih]);
		tr << ih << "/" << hits.size() << endl;
		for ( Size jh = ih+1; jh <= hits.size(); ++jh ) {
			Hit & hj(hits[jh]);
			if ( hi.ir==hj.ir||hi.jr==hj.ir||hi.ir==hj.jr||hi.jr==hj.jr ) continue;
			Size nh=0,nc=0,nd=0;
			if ( hi.itype==CYS ) nc++;
			if ( hi.jtype==CYS ) nc++;
			if ( hj.itype==CYS ) nc++;
			if ( hj.jtype==CYS ) nc++;
			if ( hi.itype==HIS1 ) nh++;
			if ( hi.jtype==HIS1 ) nh++;
			if ( hj.itype==HIS1 ) nh++;
			if ( hj.jtype==HIS1 ) nh++;
			if ( hi.itype==HIS2 ) nh++;
			if ( hi.jtype==HIS2 ) nh++;
			if ( hj.itype==HIS2 ) nh++;
			if ( hj.jtype==HIS2 ) nh++;
			if ( hi.itype==ASP1 ) nd++;
			if ( hi.jtype==ASP1 ) nd++;
			if ( hj.itype==ASP1 ) nd++;
			if ( hj.jtype==ASP1 ) nd++;
			if ( hi.itype==ASP2 ) nd++;
			if ( hi.jtype==ASP2 ) nd++;
			if ( hj.itype==ASP2 ) nd++;
			if ( hj.jtype==ASP2 ) nd++;
			if ( hi.itype==ASP3 ) nd++;
			if ( hi.jtype==ASP3 ) nd++;
			if ( hj.itype==ASP3 ) nd++;
			if ( hj.jtype==ASP3 ) nd++;
			// if(nc+nd!=2) continue;
			// if(nd > 2) continue;
			// if(nh < 2 || nh > 3) continue;
			if ( nh!=2 ) continue;
			if ( hi.cen.distance_squared(hj.cen) < MIN_DIS_METAL_SEP*MIN_DIS_METAL_SEP ) continue;
			Vec c2cen = (hi.cen+hj.cen) / 2.0;
			if ( !clashcheckhalf(p2,c2cen) ) continue;
			Vec c2ori = (hi.cen-hj.cen).normalized();
			for ( Size iaxs = 0; iaxs < 360; iaxs++ ) {
				Vec c2axs = rotation_matrix_degrees(c2ori,(Real)iaxs)*c2ori.cross(Vec(1,0,0));
				Mat c2rot = rotation_matrix_degrees(c2axs,180.0);
				// axes must be opposite ~180°
				if ( hi.axs.dot(c2rot*hj.axs) > -0.984807753012208 ) continue;// cos(10°)
				// oris must be 90° rotated
				if ( fabs(hi.ori.dot(c2rot*hj.ori)) > 0.17364817766693041 ) continue; //cos(80°)
				bool clash = false;
				for ( Size ir = 1; ir <= p.size(); ++ir ) {
					for ( Size ia = 2; ia <= 2; ++ia ) {
						if ( !clashcheck(p2,c2rot*(p2.xyz(AtomID(ia,ir))-c2cen)+c2cen) ) clash=true;
					} if ( clash ) break;
				}
				if ( clash ) continue;
				tr << "FULLHIT!!!!!!" << endl;
				/*////////////////////////////////////////////////////////////////////////////////*/ tr << "pose setup" << endl; /*//////////////////////*/
				Pose q(p2);
				q.replace_residue(hi.ir,res[hi.itype].residue(1),true);
				q.replace_residue(hi.jr,res[hi.jtype].residue(1),true);
				q.replace_residue(hj.ir,res[hj.itype].residue(1),true);
				q.replace_residue(hj.jr,res[hj.jtype].residue(1),true);
				q.set_chi(1,hi.ir,hi.ich1);
				q.set_chi(2,hi.ir,hi.ich2);
				q.set_chi(1,hi.jr,hi.jch1);
				q.set_chi(2,hi.jr,hi.jch2);
				q.set_chi(1,hj.ir,hj.ich1);
				q.set_chi(2,hj.ir,hj.ich2);
				q.set_chi(1,hj.jr,hj.jch1);
				q.set_chi(2,hj.jr,hj.jch2);
				trans_pose(q,-c2cen);
				alignaxis(q,Vec(0,0,1),c2axs);
				core::pose::symmetry::make_symmetric_pose(q);
				/*////////////////////////////////////////////////////////////////////////////////*/ tr << "csts" << endl; /*//////////////////////*/
				using namespace core::scoring::constraints;
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.itype],hi.ir)              , AtomID(matom[hi.jtype],hi.jr)              , new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hj.itype],hj.ir)              , AtomID(matom[hj.jtype],hj.jr)              , new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.itype],hi.ir              ), AtomID(matom[hj.itype],hj.ir+p.size()), new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.itype],hi.ir+p.size()), AtomID(matom[hj.itype],hj.ir              ), new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.jtype],hi.jr              ), AtomID(matom[hj.itype],hj.ir+p.size()), new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.jtype],hi.jr+p.size()), AtomID(matom[hj.itype],hj.ir              ), new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.itype],hi.ir              ), AtomID(matom[hj.jtype],hj.jr+p.size()), new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.itype],hi.ir+p.size()), AtomID(matom[hj.jtype],hj.jr              ), new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.jtype],hi.jr              ), AtomID(matom[hj.jtype],hj.jr+p.size()), new HarmonicFunc(0.0,0.1) ) );
				q.add_constraint( new AtomPairConstraint( AtomID(matom[hi.jtype],hi.jr+p.size()), AtomID(matom[hj.jtype],hj.jr              ), new HarmonicFunc(0.0,0.1) ) );

				q.add_constraint( new AngleConstraint( AtomID(batom[hi.itype],hi.ir              ), AtomID(matom[hi.itype],hi.ir              ), AtomID(batom[hi.jtype],hi.jr              ), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hj.itype],hj.ir              ), AtomID(matom[hj.itype],hj.ir              ), AtomID(batom[hj.jtype],hj.jr              ), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hi.itype],hi.ir              ), AtomID(matom[hi.itype],hi.ir              ), AtomID(batom[hj.itype],hj.ir+p.size()), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hi.itype],hi.ir+p.size()), AtomID(matom[hi.itype],hi.ir+p.size()), AtomID(batom[hj.itype],hj.ir              ), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hi.jtype],hi.jr              ), AtomID(matom[hi.jtype],hi.jr              ), AtomID(batom[hj.itype],hj.ir+p.size()), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hi.jtype],hi.jr+p.size()), AtomID(matom[hi.jtype],hi.jr+p.size()), AtomID(batom[hj.itype],hj.ir              ), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hi.itype],hi.ir              ), AtomID(matom[hi.itype],hi.ir              ), AtomID(batom[hj.jtype],hj.jr+p.size()), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hi.itype],hi.ir+p.size()), AtomID(matom[hi.itype],hi.ir+p.size()), AtomID(batom[hj.jtype],hj.jr              ), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hi.jtype],hi.jr              ), AtomID(matom[hi.jtype],hi.jr              ), AtomID(batom[hj.jtype],hj.jr+p.size()), new CircularHarmonicFunc(1.91063485009,0.2) ) );
				q.add_constraint( new AngleConstraint( AtomID(batom[hi.jtype],hi.jr+p.size()), AtomID(matom[hi.jtype],hi.jr+p.size()), AtomID(batom[hj.jtype],hj.jr              ), new CircularHarmonicFunc(1.91063485009,0.2) ) );

				//sf->show(q);
				//q.dump_pdb("test0.pdb");
				protocols::minimization_packing::symmetry::SymMinMover( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-3, true, false, false ).apply(q);
				sf->show(q);
				//q.dump_pdb("test1.pdb");
				//utility_exit_with_message("aorsitn");

				q.dump_pdb(utility::file_basename(fn)+"_"+lzs(ih,3)+lzs(jh,3)+lzs(iaxs,3)+str(nc)+"C"+str(nd)+"D"+str(nh)+"H"+"_A.pdb");
				//utility_exit_with_message("aosnrt");
			}
		}
	}

}


int main(int argc, char *argv[]) {

	try {

		devel::init(argc,argv);
		using namespace basic::options::OptionKeys;
		// loop over input files, do some checks, call dock
		for ( Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn ) {
			string fn = option[in::file::s]()[ifn];
			Pose pnat;
			tr << "checking " << fn << std::endl;
			core::import_pose::pose_from_file(pnat,fn, core::import_pose::PDB_file);
			if ( pnat.size() < 20 ) continue;
			if ( pnat.size() > 250 ) continue;
			core::scoring::dssp::Dssp dssp(pnat);
			dssp.insert_ss_into_pose(pnat);
			Size cyscnt=0, nhelix=0;
			if ( pnat.size() > MAX_NRES ) goto cont1;
			for ( Size ir = 2; ir <= pnat.size()-1; ++ir ) {
				if ( pnat.secstruct(ir) == 'H' ) nhelix++;
				//if(!pnat.residue(ir).is_protein()) goto cont1;
				if ( pnat.residue(ir).is_lower_terminus() ) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
				if ( pnat.residue(ir).is_upper_terminus() ) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
				if ( pnat.residue(ir).name3()=="CYS" ) { if ( ++cyscnt > MAX_CYS_RES ) goto cont1; }
			} goto done1; cont1: tr << "skipping " << fn << std::endl; continue; done1:
			dock(pnat,fn);
		}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

