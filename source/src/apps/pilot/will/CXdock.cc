// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#define CONTACT_TH 0
#define MIN_HELEX_RES 20
#define MAX_CYS_RES 3

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <numeric/xyzVector.hh>
#include <protocols/sic_dock/SICFast.hh>
#include <protocols/sic_dock/RigidScore.hh>
#include <protocols/sic_dock/util.hh>
// #include <protocols/sic_dock/designability_score.hh>

#include <apps/pilot/will/will_util.ihh>

#include <numeric/xyzTransform.hh>

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
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;

static THREAD_LOCAL basic::Tracer TR( "CXdock" );
static core::io::silent::SilentFileData sfd;

OPT_1GRP_KEY( FileVector, cxdock, bench2 )
OPT_1GRP_KEY( FileVector, cxdock, bench3 )
OPT_1GRP_KEY( FileVector, cxdock, bench4 )
OPT_1GRP_KEY( FileVector, cxdock, bench5 )
OPT_1GRP_KEY( FileVector, cxdock, bench6 )
OPT_1GRP_KEY( IntegerVector, cxdock, nfold )
OPT_1GRP_KEY( Integer, cxdock, max_res )
OPT_1GRP_KEY( Real, cxdock, contact_dis )
OPT_1GRP_KEY( Real, cxdock, clash_dis )
OPT_1GRP_KEY( Real, cxdock, sample )
OPT_1GRP_KEY( Real, cxdock, redundancy_angle_cut )
OPT_1GRP_KEY( Real, cxdock, redundancy_dist_cut )

typedef numeric::xyzTransform<core::Real> Xform;

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( cxdock::bench2 , "benchmark dimers (must be aligned-z)" , "" );
	NEW_OPT( cxdock::bench3 , "benchmark trimers (must be aligned-z)" , "" );
	NEW_OPT( cxdock::bench4 , "benchmark tetra (must be aligned-z)" , "" );
	NEW_OPT( cxdock::bench5 , "benchmark penta (must be aligned-z)" , "" );
	NEW_OPT( cxdock::bench6 , "benchmark hexa (must be aligned-z)" , "" );
	NEW_OPT( cxdock::nfold  , "C syms to test", 1 );
	NEW_OPT( cxdock::max_res  , "", 99999 );
	NEW_OPT( cxdock::contact_dis , "", 10.0 );
	NEW_OPT( cxdock::clash_dis   , "",  4.0 );
	NEW_OPT( cxdock::sample      , "",  5.0 );
	NEW_OPT( cxdock::redundancy_angle_cut, "", 10.0 );
	NEW_OPT( cxdock::redundancy_dist_cut, "",   4.0 );
}

#ifdef USE_OPENMP
#include <omp.h>
#endif
int num_threads() {
	#ifdef USE_OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
}
int thread_num() {
	#ifdef USE_OPENMP
		return omp_get_thread_num();
	#else
		return 0;
	#endif
}

// info about a "hit"
struct Hit {
	int iss,irt,sym;
	Real rscore,rmsd;
	Xform x1,x2;
	Hit(int is, int ir, Real cb, int sm) : iss(is),irt(ir),sym(sm),rscore(cb) {}
};
bool cmpscore (Hit i,Hit j) { return i.rscore    > j.rscore   ; }
bool cmprmsd(Hit i,Hit j) { return i.rmsd   < j.rmsd; }


// core::Real
// compute_xform_score(
// 	protocols::sic_dock::XfoxmScore const & xfs,
// 	Xform const & x1,
// 	Xform const & x2,
// 	utility::vector1<Xform> const & stubs,
// 	utility::vector1<char> const & ss
// ){
// 	assert(ss.size()==stubs.size());
// 	core::Real score = 0.0;
// 	utility::vector1<char>::const_iterator ssi = ss.begin();
// 	for(utility::vector1<Xform>::const_iterator i = stubs.begin(); i != stubs.end(); ++i,++ssi){
// 		utility::vector1<char>::const_iterator ssj = ss.begin();
// 		for(utility::vector1<Xform>::const_iterator j = stubs.begin(); j != stubs.end(); ++j,++ssj){
// 			Xform const a( x1.R * i->M, x1.R * i->v + x1.t );
// 			Xform const b( x2.R * j->M, x2.R * j->v + x2.t );
// 			if( (a.t-b.t).length_squared() > 64.0 ) continue;
// 			score += xfs.score(a,b,*ssi,*ssj);
// 		}
// 	}
// 	return score;
// }

// core::Real
// compute_xform_score(
// 	protocols::sic_dock::XfoxmScore const & xfs,
// 	core::pose::Pose const & pose1,
// 	core::pose::Pose const & pose2
// ){
// 	float tot_score = 0.0;
// 	for(core::Size ir = 1; ir <= pose1.n_residue(); ++ir){
// 		if(!pose1.residue(ir).is_protein()) continue;
// 		if(!pose1.residue(ir).has("CB")) continue;
// 		Vec CBi = pose1.residue(ir).xyz("CB");
// 		Vec CAi = pose1.residue(ir).xyz("CA");
// 		Vec  Ni = pose1.residue(ir).xyz( "N");
// 		Xform sir(CBi,CAi,Ni);
// 		for(core::Size jr = 1; jr <= pose2.n_residue(); ++jr){
// 			if(!pose2.residue(jr).is_protein()) continue;
// 			if(!pose2.residue(jr).has("CB")) continue;
// 			Vec CBj = pose2.residue(jr).xyz("CB");
// 			Vec CAj = pose2.residue(jr).xyz("CA");
// 			Vec  Nj = pose2.residue(jr).xyz( "N");
// 			if( CBi.distance_squared(CBj) > 64.0 ) continue;
// 			Xform sjr(CBj,CAj,Nj);
// 			core::Real tmp = xfs.score(sir,sjr,pose1.secstruct(ir),pose2.secstruct(jr));
// 			tot_score += tmp;
// 		}
// 	}
// 	return tot_score;
// }

void
make_native_olig(
	core::pose::Pose const & pose,
	core::pose::Pose       & olig,
	int nfold
){
	olig = pose;
	Pose tmp(pose);
	for(int i = 2; i <= nfold; ++i){
		rot_pose(tmp,Vec(0,0,1),360.0/nfold);
		for(Size i = 1; i <= tmp.n_residue(); ++i){
			if(olig.residue(i).is_lower_terminus()||olig.residue(i).is_ligand()){
				olig.append_residue_by_jump(tmp.residue(i),1);
			} else {
				olig.append_residue_by_bond(tmp.residue(i));
			}
		}
	}
}

// inline
// Xform
// get_cx_stub(Hit h){
// 	Real halfside = (h.x1.t-h.x2.t).length()/2.0;
// 	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
// 	Xform x(h.x1);
// 	x.t = x.t + Vec(0,r,halfside);
// 	{
// 		Real ang = dihedral_degrees(Vec(0,0,1),Vec(0,0,0),Vec(1,0,0), x*(Vec(0,0,0)) );
// 		Mat R = numeric::x_rotation_matrix_degrees( -ang );
// 		x.R = R * x.R;
// 		x.t = R * x.t;
// 	}
// 	return Xform(x.R,x.t);
// }

inline Vec get_rot_center(Xform const & xsym, int sym){
	Vec cen(0,0,0),pt(0,0,0);
	for(int i = 1; i <= sym; ++i){ cen += pt; pt = xsym*pt; }
	cen /= Real(sym);
	return cen;
}
inline Vec get_rot_center(Xform const & x1, Xform const & x2, int sym){
	return get_rot_center( x2 * ~x1, sym );
}

inline
Xform
get_cx_xform(Hit const & h){
	// Real halfside = (h.x1.t-h.x2.t).length()/2.0;
	// Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
	// cout << halfside << " " << r << endl;
	// Xform x(h.x2);
	// x.t = x.t + Vec(0,r,halfside);
	// Mat R = rotation_matrix_degrees(Vec(1,0,1),180.0);
	// return Xform( R*x.R, R*x.t );
	Xform x1(h.x1.R,h.x1.t),x2(h.x2.R,h.x2.t);
	Vec cen = get_rot_center(x1,x2, h.sym);
	// cout << cen << endl;
	return Xform( h.x2.R, -cen );
}

void
make_dock_olig(
	core::pose::Pose const & pose,
	core::pose::Pose       & olig,
	Hit h
){
	// Pose tmp1(pose);
	// Real halfside = (h.x1.t-h.x2.t).length()/2.0;
	// Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
	// xform_pose(tmp1,h.x1);
	// trans_pose(tmp1,Vec(0,0,halfside));
	// trans_pose(tmp1,Vec(0,r,0));
	// rot_pose(tmp1,Vec(1,0,1),180.0);
	// olig = tmp1;
	// tmp1.dump_pdb(outdir+tag);
	// make_native_olig(tmp1,olig,h.sym);

	Xform x = get_cx_xform(h);
	Mat swapXZ = rotation_matrix_degrees(Vec(1,0,1),180.0); // swap axis
	x = swapXZ * x;
	olig = pose;
	xform_pose(olig,x);
}


inline
Real
get_rmsd(
	utility::vector1<Vec> const & native_ca,
	utility::vector1<Vec> const & init_ca,
	Hit const & h//,
	// Pose pose
){
	Real halfside = (h.x1.t-h.x2.t).length()/2.0;
	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
	Xform x(h.x1);
	x.t = x.t + Vec(0,r,halfside);
	{
		Real ang = dihedral_degrees(Vec(0,0,1),Vec(0,0,0),Vec(1,0,0), x*Vec(0,0,0) );
		Mat R = numeric::x_rotation_matrix_degrees( -ang );
		x.R = R * x.R;
		x.t = R * x.t;
	}

	Mat Rz = numeric::z_rotation_matrix_degrees(180.0);
	Real rmsd=0.0, rev_rmsd=0.0;
	utility::vector1<Vec>::const_iterator i=native_ca.begin(), j=init_ca.begin();
	for(; i != native_ca.end(); ++i,++j){
		Vec v = x*(*j);
		rmsd     += i->distance_squared(      v );
		rev_rmsd += i->distance_squared( Rz * v );
	}
	if(rev_rmsd < rmsd){
		rmsd = rev_rmsd;
	}
	rmsd = sqrt( rmsd/(Real)native_ca.size() );

	return rmsd;
}

inline
Real
get_rmsd_debug(
	utility::vector1<Vec> const & native_ca,
	utility::vector1<Vec> const & init_ca,
	Hit const & h,
	Pose pose
){
	Real halfside = (h.x1.t-h.x2.t).length()/2.0;
	Real r = halfside / tan(numeric::constants::d::pi/(Real)h.sym); // tan(a/2) * edge_len/2
	Xform x(h.x1);
	x.t = x.t + Vec(0,r,halfside);
	{
		Real ang = dihedral_degrees(Vec(0,0,1),Vec(0,0,0),Vec(1,0,0), x*(Vec(0,0,0)) );
		Mat R = numeric::x_rotation_matrix_degrees( -ang );
		x.R = R * x.R;
		x.t = R * x.t;
	}

	xform_pose(pose,x);

	Mat Rz = numeric::z_rotation_matrix_degrees(180.0);
	Real rmsd=0.0, rev_rmsd=0.0;
	utility::vector1<Vec>::const_iterator i=native_ca.begin(), j=init_ca.begin();
	for(; i != native_ca.end(); ++i,++j){
		Vec v = x*(*j);
		rmsd     += i->distance_squared(      v );
		rev_rmsd += i->distance_squared( Rz * v );
	}
	if(rev_rmsd < rmsd){
		rmsd = rev_rmsd;
		rot_pose(pose,Vec(0,0,1),180.0);
	}
	rmsd = sqrt( rmsd/(Real)native_ca.size() );


	pose.dump_pdb("align_rmsd.pdb");
	// utility_exit_with_message("foo");

	return rmsd;
}

// dock pose against rotated version of itself
// rotate all different ways by looping over rotation axes on unit sphere point
// samples in ssamp, then loopint over rotation magnitudes from 1-180°
// picks cases with most CB contacts in helices
void
dock(
	Pose const & init_pose,
	std::string const & fn,
	vector1<Vec> const & ssamp,
	Pose const & /*native_olig*/,
	int native_nfold = 1,
	vector1<Vec> const & native_ca = vector1<Vec>()
){
	TR << "dock " << fn << endl;

	using namespace basic::options;
	using namespace OptionKeys;

	bool bench = native_nfold != 1;

	// set up one SIC per thread
	protocols::sic_dock::SICFast sic(option[cxdock::clash_dis]());
	sic.init(init_pose);

	protocols::sic_dock::JointScoreOP rigidsfxn = new protocols::sic_dock::JointScore;
	protocols::sic_dock::CBScoreCOP cbscore = new protocols::sic_dock::CBScore(init_pose,init_pose,option[cxdock::clash_dis](),option[cxdock::contact_dis]());
	//lnscore_ = new protocols::sic_dock::LinkerScore(init_pose,init_pose,option[tcdock::max_linker_len](),option[tcdock::linker_lookup_radius]());
	rigidsfxn->add_score(cbscore,1.0);
	//rigidsfxn->add_score(lnscore_,1.0);

	utility::vector1<Vec> init_ca;
	for(Size ir = 1; ir <= init_pose.n_residue(); ++ir){
		if(init_pose.residue(ir).has("CA")){
			init_ca.push_back(init_pose.residue(ir).xyz("CA"));
		}
	}

	//store initial coords & angle samples 1-180 (other distribution of samps would be better...)
	Real sampang = option[cxdock::sample]();
	vector1<core::Real> asamp; for(Real i = sampang/2.0; i < 180.0; i += sampang) asamp.push_back(i);

	vector1<int> syms;
	if(!bench) syms = option[cxdock::nfold]();
	else syms.push_back(native_nfold);

	// setup rotation matricies
	vector1<Mat> Rsym(syms.size());
	for(int ic = 1; ic <= (int)syms.size(); ic++) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/Real(syms[ic]));
	vector1<vector1<Hit> > hits(syms.size());

	// utility::vector1<Xform> stubs;
	// utility::vector1<char> ss;
	// for(Size ir = 1; ir <= init_pose.n_residue(); ++ir){
	// 	if( !init_pose.residue(ir).is_protein() ) continue;
	// 	if( !init_pose.residue(ir).has("CB") ) continue;
	// 	if( init_pose.secstruct(ir) == 'L' ) continue;
	// 	stubs.push_back( Xform(init_pose.residue(ir).xyz("CB"),init_pose.residue(ir).xyz("CA"),init_pose.residue(ir).xyz("N")));
	// 	ss.push_back(init_pose.secstruct(ir));
	// }

	Size totsamp = ssamp.size()*asamp.size();
	Size outinterval = min((Size)50000,totsamp/10);
	// loop over axes on sphere points & rotaton amounts 1-180
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for(int iss = 1; iss <= (int)ssamp.size(); ++iss){
		for(int irt = 1; irt <= (int)asamp.size(); ++irt) {
			Size progress = (iss-1)*asamp.size()+irt;
			if(progress%outinterval==0) TR << progress << " of " << totsamp << std::endl;
			Mat R = rotation_matrix_degrees(ssamp[iss],(Real)asamp[irt]);
			Xform const x1( R, Vec(0,0,0));
			// loop over symmitreis to check
			for(int ic = 1; ic <= (int)syms.size(); ic++){
				Xform const x2( Rsym[ic]*R, Vec(0,0,0) );
				// Real t = sic.slide_into_contact(x1,x2,Vec(0,0,1));
				Real rscore = 0.0;
				Xform tmp(x1);
				Real const t = slide_into_contact_and_score(sic,*rigidsfxn,tmp,x2,Vec(0,0,1),rscore);

				Hit h(iss,irt,rscore,syms[ic]);
				h.x1 = x1;
				h.x2 = x2;
				h.x1.t += t*Vec(0,0,1);
				// h.xscore = compute_xform_score(xfs,h.x1,h.x2,stubs,ss);
				if(bench) h.rmsd = get_rmsd( native_ca, init_ca, h );
				// Xform x = get_cx_stub(h);

				// cout << "XFORM " << x.R << " " << x.t << endl;

				if(rscore >= CONTACT_TH){
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					{
						hits[ic].push_back(h);

						// Pose olig;
						// make_dock_olig(init_pose,olig,h);
						// Real actualrmsd = core::scoring::CA_rmsd(olig,native_olig);

							// std::cout << "RMSD "
							//           // << actualrmsd << " "
							//           << h.rmsd << " "
							//           << h.rscore << " "
							//           << h.xscore << " "
							//           << std::endl;

						// if( fabs(actualrmsd - h.rmsd) > 5.0 && h.rmsd < 5.0 ){
						// 	get_rmsd_debug( native_ca, init_ca, h, init_pose );
						// 	// core::scoring::calpha_superimpose_pose(olig,native_olig);
						// 	olig.dump_pdb("olig.pdb");
						// 	native_olig.dump_pdb("native.pdb");
						// 	utility_exit_with_message("test rmsd");
						// }


					}
				}
			}
		}
	}

	Size nstruct = basic::options::option[basic::options::OptionKeys::out::nstruct]();
	std::string outdir = "./";
	if(basic::options::option[basic::options::OptionKeys::out::file::o].user()){
		outdir = basic::options::option[basic::options::OptionKeys::out::file::o].user()+"/";
	}

	// // test
	// for(int ic = 1; ic <= (int)syms.size(); ic++) {
	// 	std::sort(hits[ic].begin(),hits[ic].end(),cmpscore);
	// 	for(int i = 1; i <= (int)hits[ic].size(); ++i) {
	// 		Hit & h(hits[ic][i]);

	// 		// Pose tmp1(init_pose),tmp2(init_pose);
	// 		// xform_pose(tmp1,h.x1);
	// 		// xform_pose(tmp2,h.x2);
	// 		// tmp1.dump_pdb("test1.pdb");
	// 		// tmp2.dump_pdb("test2.pdb");

	// 		// Xform const a( h.x1.R * stubs[1].R, h.x1.R * stubs[1].t + h.x1.t );
	// 		// Xform const b( h.x2.R * stubs[1].R, h.x2.R * stubs[1].t + h.x2.t );
	// 		// std::cout << Xform(tmp1.residue(1).xyz("CB"),tmp1.residue(1).xyz("CA"),tmp1.residue(1).xyz("N")) << std::endl;
	// 		// std::cout << a << std::endl;
	// 		// std::cout << Xform(tmp2.residue(1).xyz("CB"),tmp2.residue(1).xyz("CA"),tmp2.residue(1).xyz("N")) << std::endl;
	// 		// std::cout << b << std::endl;

	// 		// std::cout << compute_xform_score(xfs,tmp1,tmp2) << " == "
	// 		//           << compute_xform_score(xfs,h.x1,h.x2,stubs,ss) << " "
	// 		//           << h.rscore << " " << h.xscore << std::endl;;

	// 		// utility_exit_with_message("test");

	// 		// std::cout <<"SCORES " << h.rscore <<" "<< compute_xform_score(xfs,h.x1,h.x2,stubs,ss) << " " << h.xscore << std::endl;;

	// 	}
	// }


	// dump top results
	vector1<Hit> dumpedit;
	for(int ic = 1; ic <= (int)syms.size(); ic++) {
		TR << "for sym C" << syms[ic] << " found " << hits[ic].size() << " hits" << endl;
		std::sort(hits[ic].begin(),hits[ic].end(),cmpscore);
		Size nout=0,i=0;
		while( nout < nstruct && ++i <= hits[ic].size() ){
		// for( i = 1; i <= min(nstruct,); ++i) {
			Hit & h(hits[ic][i]);

			// for debugging
			// Pose tmp1(init_pose),tmp2(init_pose);
			// xform_pose(tmp1,h.x1);
			// xform_pose(tmp2,h.x2);
			// tmp1.dump_pdb("tmp1.pdb");
			// tmp2.dump_pdb("tmp2.pdb");
			// utility_exit_with_message("test");

			// redundency check
			bool toosimilar = false;
			for(vector1<Hit>::const_iterator j = dumpedit.begin(); j != dumpedit.end(); ++j){
				Xform xrel = h.x1 * ~(j->x1);
				Real ang=0.0;
				numeric::rotation_axis(xrel.R,ang);
				ang = numeric::conversions::degrees(ang);
				Real dis = xrel.t.length();
				if( ang < option[cxdock::redundancy_angle_cut]() || dis < option[cxdock::redundancy_dist_cut]() ){
					toosimilar = true;
					// debug
					// std::cerr << dis << " " << ang << std::endl;
					// if( 1 ){
					// 	Pose tmp1,tmp2;
					// 	make_dock_olig(init_pose,tmp1, h);
					// 	make_dock_olig(init_pose,tmp2,*j);
					// 	tmp1.dump_pdb("tosim1.pdb");
					// 	tmp2.dump_pdb("tosim2.pdb");
					// 	utility_exit_with_message("redundancy_cut test");
					// }
				}
			}
			if(toosimilar) continue;
			dumpedit.push_back(h);

			std::string tag = utility::file_basename(fn);
			if(tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
			if(tag.substr(tag.size()-4)==".pdb") tag = tag.substr(0,tag.size()-4);
			tag += "_C"+ObjexxFCL::string_of(h.sym)+"_score"+ObjexxFCL::string_of(i)+".pdb";
			cout << "RESULT " << h.sym << " " << h.iss << " " << ssamp.size()  << " " << h.irt << " " << h.rscore << " " << tag << endl;
			Pose tmp;
			make_dock_olig(init_pose,tmp,h);
			// TR << "dumping top rscore to " << outdir+tag << endl;
			tmp.dump_pdb(outdir+tag);
			++nout;
			// utility_exit_with_message("test");
		}
	}
	// // dump top results
	// for(int ic = 1; ic <= (int)syms.size(); ic++) {
	// 	std::sort(hits[ic].begin(),hits[ic].end(),cmpxsc);
	// 	for(int i = 1; i <= (int)min(nstruct,hits[ic].size()); ++i) {
	// 		Hit & h(hits[ic][i]);
	// 		std::string tag = utility::file_basename(fn);
	// 		if(tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
	// 		if(tag.substr(tag.size()-4)==".pdb") tag = tag.substr(0,tag.size()-4);
	// 		tag += "_C"+ObjexxFCL::string_of(h.sym)+"_xsc"+ObjexxFCL::string_of(i)+".pdb";
	// 		cout << "RESULT " << h.sym << " " << h.iss << " " << ssamp.size()  << " " << h.irt << " " << h.rscore << " " << tag << endl;
	// 		Pose tmp;
	// 		make_dock_olig(init_pose,tmp,h);
	// 		tmp.dump_pdb( outdir+tag);
	// 		// utility_exit_with_message("test");
	// 	}
	// }

	if(bench){
		// dump top results
		for(int ic = 1; ic <= (int)syms.size(); ic++) {
			std::sort(hits[ic].begin(),hits[ic].end(),cmprmsd);
			for(int i = 1; i <= (int)min(nstruct,hits[ic].size()); ++i) {
				Hit & h(hits[ic][i]);
				std::string tag = utility::file_basename(fn);
				if(tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
				if(tag.substr(tag.size()-4)==".pdb") tag = tag.substr(0,tag.size()-4);
				tag += "_C"+ObjexxFCL::string_of(h.sym)+"_rmsd"+ObjexxFCL::string_of(i)+".pdb";
				cout << "RESULT " << h.sym << " " << h.iss << " " << ssamp.size()  << " " << h.irt << " " << h.rscore << " " << tag << endl;
				Pose tmp;
				make_dock_olig(init_pose,tmp,h);
				TR << "dumping top rmsd to " << outdir+tag << endl;
				tmp.dump_pdb(outdir+tag);
				// utility_exit_with_message("test");
			}
		}
	}

	if(bench){ // benchmark
		// Real rms_goodscore=999.0, rms_goodxsc=999.0;
		std::cout << "RESULT RMSD for best by score";
		std::sort(hits[1].begin(),hits[1].end(),cmpscore);
		for(int i = 1; i <= (int)min(20, int(hits[1].size()) ); ++i) std::cout << " " << hits[1][i].rmsd;
		std::cout << endl;

		// std::cout << "RESULT RMSD for best by xform score  ";
		// std::sort(hits[1].begin(),hits[1].end(),cmpxsc);
		// for(int i = 1; i <= (int)min(20ul,hits[1].size()); ++i) std::cout << " " << hits[1][i].rmsd;
		// std::cout << endl;

		Real rms_goodscore=999.0; //, rms_goodxsc=999.0;
		std::sort(hits[1].begin(),hits[1].end(),cmpscore);
		for(int i = 1; i <= (int)min(20, int(hits[1].size()) ); ++i) rms_goodscore = min(rms_goodscore,hits[1][i].rmsd);
		Real best_score = hits[1][1].rmsd;
		// std::sort(hits[1].begin(),hits[1].end(),cmpxsc);
		// for(int i = 1; i <= (int)min(20ul,hits[1].size()); ++i) rms_goodxsc = min(rms_goodxsc,hits[1][i].rmsd);
		// Real best_xsc = hits[1][1].rmsd;
		// rms_goodscore /= 10.0;
		// rms_goodxsc /= 10.0;
		std::cout << "Best RMSD in top 20 by num. contacts: " << rms_goodscore << " top " << best_score << std::endl;
		// std::cout << "Best RMSD in top 20 by xform score:   " << rms_goodxsc << " top " << best_xsc << std::endl;

	}


}

utility::vector1<Vec>
align_native_state(
	core::pose::Pose & pose,
	int nfold
){
	// align so that axis is X and COM is along Z
	// for convenience of RMSD calc
	if(nfold<2) return utility::vector1<Vec>();
	Vec cen(0,0,0);
	Real ncen = 0.0;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		if(pose.residue(ir).has("CA")){
			cen += pose.residue(ir).xyz("CA");
			ncen += 1.0;
		}
	}
	cen /= ncen;
	rot_pose( pose, Vec(0,0,1), -dihedral_degrees(Vec(1,0,0),Vec(0,0,0),Vec(0,0,1),cen) );
	rot_pose( pose, Vec(1,0,1), 180.0 ); // symZ to symX
	utility::vector1<Vec> CAs;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		if(pose.residue(ir).has("CA")){
			CAs.push_back(pose.residue(ir).xyz("CA"));
		}
	}
	return CAs;
}

void
read_sphere(
	vector1<Vec> & ssamp
){
	// 5.0°       sphere_672.dat.gz
	// 4.0°       sphere_1032.dat.gz
	// 3.0°       sphere_1812.dat.gz
	// 2.0°       sphere_4092.dat.gz
	// 1.4°       sphere_8192.dat.gz
	// 1.0°       sphere_17282.dat.gz
	Real ang = min(5.0,max(1.0,round(basic::options::option[basic::options::OptionKeys::cxdock::sample]())));
	Size nss = 672;
	if(     4.5 <= ang && ang < 5.5 ){ ang=5.0; nss =   672; }
	else if(3.5 <= ang && ang < 4.5 ){ ang=4.0; nss =  1032; }
	else if(2.5 <= ang && ang < 3.5 ){ ang=3.0; nss =  1812; }
	else if(1.7 <= ang && ang < 2.5 ){ ang=2.0; nss =  4092; }
	else if(1.2 <= ang && ang < 1.7 ){ ang=1.4; nss =  8192; }
	else if(0.8 <= ang && ang < 1.2 ){ ang=1.0; nss = 17282; }
	else utility_exit_with_message("sampling level unsupported, 1.0 - 5.0 degrees ( currently 1.0, 1.4, 2.0, 3.0, 4.0, 5.0 )");
	TR << "sphere sampling resolution: " << ang << std::endl;
	ssamp.resize(nss);
	izstream is;
	if(!basic::database::open(is,"sampling/spheres/sphere_"+str(nss)+".dat"))
		utility_exit_with_message("can't open sphere data: sampling/spheres/sphere_"+str(nss)+".dat");
	for(Size i = 1; i <= nss; ++i) {
		core::Real x,y,z;
		is >> x >> y >> z;
		ssamp[i] = Vec(x,y,z);
	}
	is.close();
}

void
get_tasks_from_command_line(
	utility::vector1<std::pair<string,int> > & tasks
){
	using namespace basic::options;
	using namespace OptionKeys;
	if(option[in::file::s].user()){
		for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[in::file::s]()[ifn],1));
		}
	}
	if(option[cxdock::bench2].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench2]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench2]()[ifn],2));
		}
	}
	if(option[cxdock::bench3].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench3]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench3]()[ifn],3));
		}
	}
	if(option[cxdock::bench4].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench4]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench4]()[ifn],4));
		}
	}
	if(option[cxdock::bench5].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench5]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench5]()[ifn],5));
		}
	}
	if(option[cxdock::bench6].user()){
		for(Size ifn = 1; ifn <= option[cxdock::bench6]().size(); ++ifn){
			tasks.push_back(std::make_pair(option[cxdock::bench6]()[ifn],6));
		}
	}
}

int main(int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;

	vector1<Vec> ssamp; read_sphere(ssamp);

	// protocols::sic_dock::XfoxmScore const xfs("/Users/sheffler/project/designability_stats/results/");

	vector1<std::pair<string,int> > tasks;
	get_tasks_from_command_line(tasks);

	for(vector1<std::pair<string,int> >::const_iterator i = tasks.begin(); i != tasks.end(); ++i){
		string fn = i->first;
		int nfold = i->second;
		TR << "searching " << fn << " " << nfold << std::endl;

		Pose pnat;
		core::import_pose::pose_from_pdb(pnat,fn);
		core::scoring::dssp::Dssp dssp(pnat);
		dssp.insert_ss_into_pose(pnat);

		// center on z
		trans_pose(pnat,Vec(0,0,-center_of_geom(pnat,1,pnat.n_residue()).z()));

		Pose olig; make_native_olig(pnat,olig,nfold);

		// sym axis now along X, rotated around X so CA com in XZ plane,
		// this makes computing the RMSD of docked states easier
		vector1<Vec> native_ca = align_native_state(pnat,nfold); // align CA com
		// pnat.dump_pdb("native_align.pdb");

		Vec cen = center_of_geom(pnat,1,pnat.n_residue());
		trans_pose(pnat,-Vec(0,cen.y(),cen.z()));

		// inpuc checks
		if( pnat.n_residue() > (Size)option[cxdock::max_res]() ) continue;
		// Size cyscnt=0, nhelix=0;
		// for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
		// 	if(pnat.secstruct(ir) == 'H') nhelix++;
		// 	//if(!pnat.residue(ir).is_protein()) goto cont1;
		// 	if(pnat.residue(ir).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
		// 	if(pnat.residue(ir).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
		// 	if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > MAX_CYS_RES) goto cont1; }
		// } goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		// if( nhelix < MIN_HELEX_RES ) continue;

		dock(pnat,fn,ssamp,olig,nfold,native_ca);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

