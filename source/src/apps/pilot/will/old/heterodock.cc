// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

#include <apps/pilot/will/pick_normal_sampling.hh>

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
using core::import_pose::pose_from_file;

static basic::Tracer TR( "heterodock" );
static core::io::silent::SilentFileData sfd;

OPT_1GRP_KEY( FileVector, heterodock, dockpairs     )
OPT_1GRP_KEY( FileVector, heterodock, allbyall      )
OPT_1GRP_KEY( FileVector, heterodock, somebyothers1 )
OPT_1GRP_KEY( FileVector, heterodock, somebyothers2 )
OPT_1GRP_KEY( Integer, heterodock, max_res )
OPT_1GRP_KEY( Real, heterodock, score_thresh )
OPT_1GRP_KEY( Real, heterodock, contact_dis )
OPT_1GRP_KEY( Real, heterodock, clash_dis )
OPT_1GRP_KEY( Real, heterodock, sample_angle )
OPT_1GRP_KEY( Real, heterodock, redundancy_angle_cut )
OPT_1GRP_KEY( Real, heterodock, redundancy_dist_cut )

typedef numeric::xyzTransform<core::Real> Xform;

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( heterodock::dockpairs         , "" , "" );
	NEW_OPT( heterodock::allbyall      , "" , "" );
	NEW_OPT( heterodock::somebyothers1 , "" , "" );
	NEW_OPT( heterodock::somebyothers2 , "" , "" );
	NEW_OPT( heterodock::max_res  , "", 99999 );
	NEW_OPT( heterodock::score_thresh , "", 0.0 );
	NEW_OPT( heterodock::contact_dis , "", 10.0 );
	NEW_OPT( heterodock::clash_dis   , "",  4.0 );
	NEW_OPT( heterodock::sample_angle      , "",  5.0 );
	NEW_OPT( heterodock::redundancy_angle_cut, "", 10.0 );
	NEW_OPT( heterodock::redundancy_dist_cut, "",   4.0 );
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


void
make_olig(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	Xform const & x1,
	Xform const & x2,
	core::pose::Pose       & olig
){
	olig = pose1;
	olig.append_residue_by_jump(pose2.residue(1),1,"","",true);
	for ( Size i = 2; i <= pose2.size(); ++i ) {
		if ( olig.residue(i).is_lower_terminus()||olig.residue(i).is_ligand() ) {
			olig.append_residue_by_jump(pose2.residue(i),1);
		} else {
			olig.append_residue_by_bond(pose2.residue(i));
		}
	}
	xform_pose(olig,x1,1,pose1.size());
	xform_pose(olig,x2,pose1.size()+1,pose1.size()+pose2.size());
}


// info about a "hit"
struct Hit {
	Real iori,irot,islide;
	Real sicscore,rmsd;
	Xform x1,x2;
	Hit(Real _iori, Real _irot, Real _islide, Real _score) : iori(_iori),irot(_irot),islide(_islide),sicscore(_score) {}
};
bool cmpscore (Hit i,Hit j) { return i.sicscore > j.sicscore   ; }
bool cmprmsd(Hit i,Hit j) { return i.rmsd   < j.rmsd; }


Real myxdis(Xform const & x1, Xform const & x2){
	return sqrt(
		x1.R.col_x().distance_squared(x2.R.col_x())*900.0 +
		x1.R.col_y().distance_squared(x2.R.col_y())*900.0 +
		x1.R.col_z().distance_squared(x2.R.col_z())*900.0 +
		x1.t        .distance_squared(x2.t        )
	);
}
// dock pose against rotated version of itself
// rotate all different ways by looping over rotation axes on unit sphere point
// samples in normals, then loopint over rotation magnitudes from 1-180°
// picks cases with most CB contacts in helices
void
dock(
	Pose const & init_pose1,
	Pose const & init_pose2,
	std::string const & tag0,
	vector1<Vec> const & normals,
	Xform const & xnative
){
	using namespace basic::options;
	using namespace OptionKeys;

	//store initial coords & angle samples 1-180 (other distribution of samps would be better...)
	Real Nang = sqrt(2*(Real)normals.size())     * 2.0;
	// vector1<core::Real> asamp; for(Real i = 1.0; i < 2.0*Nang; i += 2.0) asamp.push_back(i/2.0/Nang*180.0);
	Size totsamp = normals.size()*normals.size()*Nang;
	Size outinterval = min((Size)50000,totsamp/10);

	TR << "dock " << tag0 << " " << normals.size() << " " << Nang << " total samples: " << totsamp/1000 << "K  resolution: ~" << 180.0/Nang << " degrees" << endl;

	// set up one SIC per thread -- now thread safe!
	protocols::sic_dock::SICFast sic(option[heterodock::clash_dis]());
	sic.init(init_pose1,init_pose2);

	protocols::sic_dock::JointScoreOP rigidsfxn = new protocols::sic_dock::JointScore;
	protocols::sic_dock::CBScoreCOP cbscore = new protocols::sic_dock::CBScore(init_pose1,init_pose2,option[heterodock::clash_dis](),option[heterodock::contact_dis]());
	//lnscore_ = new protocols::sic_dock::LinkerScore(init_pose,init_pose,option[tcdock::max_linker_len](),option[tcdock::linker_lookup_radius]());
	rigidsfxn->add_score(cbscore,1.0);
	//rigidsfxn->add_score(lnscore_,1.0);

	// setup rotation matricies
	vector1<Hit> hits;

	Xform const x2; // null -- keep pose2 fixed

	// loop over axes on sphere points & rotaton amounts 1-180
#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
#endif
	for ( int iori = 1; iori <= (int)normals.size(); ++iori ) {

		for ( int iang = 1; iang <= Nang; ++iang ) {

			Real ang = uniform()*3.14159;
			while ( (ang-sin(2.0*ang)/2.0)/3.14159/2.0 < uniform() ) ang = uniform()*3.14159;

			Xform const xori( rotation_matrix_degrees(normals[iori],ang), Vec(0,0,0) );

			// Xform const xori;

			// Real angerr = xori.R.col_x().distance_squared(xnative.R.col_x())*900.0 +xori.R.col_y().distance_squared(xnative.R.col_y())*900.0 + xori.R.col_z().distance_squared(xnative.R.col_z())*900.0;
			// if( angerr > 0.1 ){
			//  // cout << angerr << endl;;
			//              continue;
			//       }

			for ( int islide = 1; islide <= (int)normals.size(); ++islide ) {

				Size progress = ((iori-1)*Nang+iang)*normals.size() + islide;
				if ( progress%outinterval==0 ) {
#ifdef USE_OPENMP
					#pragma omp critical
#endif
					TR << progress/1000 << "K of " << totsamp/1000 << "K " << hits.size() << std::endl;
				}

				Vec const & slide_dir = normals[islide];

				Real sicscore = 0.0;
				Xform x1(xori);

				/*Real const t = */slide_into_contact_and_score(sic,*rigidsfxn,x1,x2,slide_dir,sicscore);

				// for(int ibackoff = 1; ibackoff <= 10; ++ibackoff){

				Hit h(ang,ang,islide,sicscore);
				h.x1 = x1;
				// h.x1.t -= ((Real)ibackoff-1)/5.0*slide_dir;
				h.x2 = x2;
				h.rmsd = myxdis(h.x1,xnative); // LEI this is a terrible hack, not a real rmsd...

				if ( sicscore >= option[heterodock::score_thresh]() ) {
#ifdef USE_OPENMP
					#pragma omp critical
#endif
					{
						hits.push_back(h);
					}
				}
				// }

			}
		}
	}

	Size nstruct = basic::options::option[basic::options::OptionKeys::out::nstruct]();
	std::string outdir = "./";
	if ( basic::options::option[basic::options::OptionKeys::out::file::o].user() ) {
		outdir = basic::options::option[basic::options::OptionKeys::out::file::o]()+"/";
	}

	// dump top results
	vector1<Hit> dumpedit;
	TR << "NUM_HITS " << tag0 << " " << hits.size() << " hits" << endl;
	// std::sort(hits.begin(),hits.end(),cmpscore);
	std::sort(hits.begin(),hits.end(),cmprmsd);
	Size nout=0,i=0;
	while ( nout < nstruct && ++i <= hits.size() ) {
		Hit & h(hits[i]);

		// // for debugging
		// Pose tmp1(init_pose1),tmp2(init_pose2); xform_pose(tmp1,h.x1); xform_pose(tmp2,h.x2); tmp1.dump_pdb("tmp1.pdb"); tmp2.dump_pdb("tmp2.pdb"); utility_exit_with_message("test");

		// redundency check
		bool toosimilar = false;
		for ( vector1<Hit>::const_iterator j = dumpedit.begin(); j != dumpedit.end(); ++j ) {
			Xform xrel = h.x1 * ~(j->x1);
			Real ang=0.0;
			numeric::rotation_axis(xrel.R,ang);
			ang = numeric::conversions::degrees(ang);
			Real dis = xrel.t.length();
			if ( ang < option[heterodock::redundancy_angle_cut]() || dis < option[heterodock::redundancy_dist_cut]() ) {
				toosimilar = true;
				// debug
				// std::cerr << dis << " " << ang << std::endl;
				// if( 1 ){
				//  Pose tmp1,tmp2;
				//  make_dock_olig(init_pose,tmp1, h);
				//  make_dock_olig(init_pose,tmp2,*j);
				//  tmp1.dump_pdb("tosim1.pdb");
				//  tmp2.dump_pdb("tosim2.pdb");
				//  utility_exit_with_message("redundancy_cut test");
				// }
			}
		}
		if ( toosimilar ) continue;
		dumpedit.push_back(h);

		// cout << xnative << endl;
		// cout << h.x1 << endl;
		// cout << h.x2 << endl;
		// cout << h.x1.t-xnative.t << endl;
		// cout << normals[h.islide] << endl;
		// cout << h.x1.R.col_x().distance(xnative.R.col_x())*30 << endl;
		// cout << h.x1.R.col_y().distance(xnative.R.col_y())*30 << endl;
		// cout << h.x1.R.col_z().distance(xnative.R.col_z())*30 << endl;
		// cout << h.x1.t.distance(xnative.t) << endl;

		string tag = tag0;
		tag += "_heterodock_top"+ObjexxFCL::string_of(i)+".pdb";
		cout << "TOPHITS nnorm: "  << normals.size() << " " << h.iori << " " << " " << h.irot << " " << h.islide << " " << h.sicscore << " " << tag << " " << h.rmsd << endl;
		// TR << "dumping top sicscore to " << outdir+tag << endl;
		Pose olig;
		make_olig(init_pose1,init_pose2,h.x1,Xform(),olig);
		olig.dump_pdb(outdir+tag);
		++nout;

		// Pose pnat;
		// make_olig(init_pose1,init_pose2,xnative,Xform(),pnat);
		// pnat.dump_pdb("native.pdb");

	}

}

void
read_sphere(
	vector1<Vec> & normals
){
	Real ang = basic::options::option[basic::options::OptionKeys::heterodock::sample_angle]();
	int nss = pick_normal_sampling( ang );
	TR << "angular sampling resolution: requested: " << ang << " actual: " << 180.0/sqrt(2*(Real)nss) << std::endl;
	normals.resize(nss);
	izstream is;
	if ( !basic::database::open(is,"sampling/spheres/sphere_"+str(nss)+".dat") ) {
		utility_exit_with_message("can't open sphere data: sampling/spheres/sphere_"+str(nss)+".dat");
	}
	for ( int i = 1; i <= nss; ++i ) {
		core::Real x,y,z;
		is >> x >> y >> z;
		normals[i] = Vec(x,y,z);
	}
	is.close();
}

void
get_tasks_from_command_line(
	utility::vector1<std::pair<string,string> > & tasks
){
	using namespace basic::options;
	using namespace OptionKeys;
	if ( option[heterodock::dockpairs].user() ) {
		vector1<string> files = option[heterodock::dockpairs]();
		if ( files.size()%2==1 ) utility_exit_with_message("heterodock::dockpairs should be even number!");
		for ( Size ifn = 1; ifn <= files.size(); ifn+=2 ) {
			tasks.push_back(std::make_pair(files[ifn],files[ifn+1]));
		}
	}
	if ( option[heterodock::allbyall     ].user() ) cerr << "all-by-all not implemented yet";
	if ( option[heterodock::somebyothers1].user() ) cerr << "somebyothers1 not implemented yet";
	if ( option[heterodock::somebyothers2].user() ) cerr << "somebyothers2 not implemented yet";
	if ( tasks.size()==0 ) utility_exit_with_message("no heterodock tasks!");
}

Pose get_centered_pose(string fname, Vec & cen){
	Pose todock;
	core::import_pose::pose_from_file(todock,fname, core::import_pose::PDB_file);
	core::scoring::dssp::Dssp dssp(todock);
	dssp.insert_ss_into_pose(todock);
	cen = center_of_geom(todock,1,todock.size());
	trans_pose(todock,-cen);
	return todock;
}

int main(int argc, char *argv[]) {
	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;

	vector1<Vec> normals; read_sphere(normals);

	vector1<std::pair<string,string> > tasks; get_tasks_from_command_line(tasks);

	for ( vector1<std::pair<string,string> >::const_iterator i = tasks.begin(); i != tasks.end(); ++i ) {
		string fn1 = i->first;
		string fn2 = i->second;
		string tag = utility::file_basename(fn1);
		if ( tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
		if ( tag.substr(tag.size()-4)==".pdb" ) tag = tag.substr(0,tag.size()-4);
		tag += "_"+utility::file_basename(fn2);
		if ( tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
		if ( tag.substr(tag.size()-4)==".pdb" ) tag = tag.substr(0,tag.size()-4);

		TR << "searching " << tag << std::endl;

		Vec cen1,cen2;
		Pose todock1 = get_centered_pose(fn1,cen1);
		Pose todock2 = get_centered_pose(fn2,cen2);

		Xform xnative; xnative.t = cen1-cen2;
		// Pose olig;
		// make_olig(todock1,todock2,xnative,Xform(),olig);
		// olig.dump_pdb("native.pdb");
		// utility_exit_with_message("arst");

		if ( todock1.size() > (Size)option[heterodock::max_res]() ) continue;
		if ( todock2.size() > (Size)option[heterodock::max_res]() ) continue;

		dock(todock1,todock2,tag,normals,xnative);
	}
}

