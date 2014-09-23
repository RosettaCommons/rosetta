#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>
	#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
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
	#include <core/pose/symmetry/util.hh>
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
	#include <protocols/sic_dock/xyzStripeHashPose.hh>
	#include <protocols/sic_dock/scores/MotifHashRigidScore.hh>
	#include <protocols/sic_dock/scores/TrisBpyScore.hh>
	#include <protocols/motif_hash/motif_hash_stuff.hh>
	#include <numeric/xyzTransform.hh>
	#include <basic/sampling/orientations/QuaternionGrid.hh>
	#include <core/scoring/constraints/ConstraintSet.hh>
	#include <core/scoring/constraints/AtomPairConstraint.hh>

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
	using numeric::x_rotation_matrix_degrees;
	using numeric::y_rotation_matrix_degrees;
	using numeric::z_rotation_matrix_degrees;
	using numeric::x_rotation_matrix;
	using numeric::y_rotation_matrix;
	using numeric::z_rotation_matrix;
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
	using namespace protocols::sic_dock;

static thread_local basic::Tracer TR( "sicdock" );
static core::io::silent::SilentFileData sfd;

OPT_1GRP_KEY( FileVector, sicdock, bench2 )
	OPT_1GRP_KEY( FileVector, sicdock, bench3 )
	OPT_1GRP_KEY( FileVector, sicdock, bench4 )
	OPT_1GRP_KEY( FileVector, sicdock, bench5 )
	OPT_1GRP_KEY( FileVector, sicdock, bench6 )
	OPT_1GRP_KEY( IntegerVector, sicdock, nfold )
	OPT_1GRP_KEY( Boolean, sicdock, dihedral )
	OPT_1GRP_KEY( Boolean, sicdock, rotate_x_only )
	OPT_1GRP_KEY( Boolean, sicdock, no_recenter )
	OPT_1GRP_KEY( Integer, sicdock, max_res )
	// OPT_1GRP_KEY( Real, sicdock, contact_dis )
	// OPT_1GRP_KEY( Real, sicdock, clash_dis )
	OPT_1GRP_KEY( Real, sicdock, sample_radius )
	OPT_1GRP_KEY( Real, sicdock, sample_number )
	OPT_1GRP_KEY( Real, sicdock, redundancy_angle_cut )
	OPT_1GRP_KEY( Real, sicdock, redundancy_dist_cut )
	OPT_1GRP_KEY( Real, sicdock, motif_hash_weight )
	OPT_1GRP_KEY( Real, sicdock, sample_thickness )
	OPT_1GRP_KEY( Real, sicdock, min_score )
	OPT_1GRP_KEY( Real, sicdock, min_contact_score_cut )
	OPT_1GRP_KEY( Integer, sicdock, nscore )
	OPT_1GRP_KEY( IntegerVector, sicdock, iori )
	OPT_1GRP_KEY( Real, sicdock, trisbpy )
	OPT_1GRP_KEY( Real, sicdock, trisbpy_min_contacts )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( sicdock::bench2               , "benchmark dimers (must be aligned-z)" , "" );
	NEW_OPT( sicdock::bench3               , "benchmark trimers (must be aligned-z)" , "" );
	NEW_OPT( sicdock::bench4               , "benchmark tetra (must be aligned-z)" , "" );
	NEW_OPT( sicdock::bench5               , "benchmark penta (must be aligned-z)" , "" );
	NEW_OPT( sicdock::bench6               , "benchmark hexa (must be aligned-z)" , "" );
	NEW_OPT( sicdock::nfold                , "C syms to test", 1 );
	NEW_OPT( sicdock::max_res              , "", 99999 );
	NEW_OPT( sicdock::dihedral             , "", "false" );
	// NEW_OPT( sicdock::contact_dis          , "", 10.0 );
	// NEW_OPT( sicdock::clash_dis            , "",  3.5 );
	NEW_OPT( sicdock::sample_radius        , "",  5.0 );
	NEW_OPT( sicdock::sample_number        , "",  70728 );
	NEW_OPT( sicdock::redundancy_angle_cut , "",  3 );
	NEW_OPT( sicdock::redundancy_dist_cut, "",   1 );
	NEW_OPT( sicdock::nscore, "",   100 );
	NEW_OPT( sicdock::motif_hash_weight, "",  10.0 );
	NEW_OPT( sicdock::sample_thickness, "",  0.0 );
	NEW_OPT( sicdock::min_score, "",  200.0 );
	NEW_OPT( sicdock::min_contact_score_cut, "",  20.0 );
	NEW_OPT( sicdock::rotate_x_only, "",  false );
	NEW_OPT( sicdock::no_recenter, "",  false );
	NEW_OPT( sicdock::iori, "",  -1 );
	NEW_OPT( sicdock::trisbpy, "", 0.1 );
	NEW_OPT( sicdock::trisbpy_min_contacts, "", 15 );
}


struct Hit {
	int iqg,sym;
	Real rscore,rmsd,dang,ddis,backoff;
	Xform x1,x2,x3,x4;
	Hit(int iq, Real cb, int sm) : iqg(iq),sym(sm),rscore(cb),rmsd(0),dang(0),ddis(0),backoff(0) {}
	Hit() : iqg(-1),sym(-1),rscore(-9e9),rmsd(0),dang(0),ddis(0),backoff(0) {}
};
bool cmpscore (Hit i,Hit j) { return i.rscore    > j.rscore   ; }
bool cmprmsd(Hit i,Hit j) { return i.rmsd   < j.rmsd; }


inline Vec get_rot_center(Xform const & xsym, int sym){
	Vec cen(0,0,0),pt(0,0,0);
	for(int i = 1; i <= sym; ++i){ cen += pt; pt = xsym*pt; }
	cen /= Real(sym);
	return cen;
}
inline Vec get_rot_center(Xform const & x1, Xform const & x2, int sym){
	return get_rot_center( x2 * ~x1, sym );
}
inline Xform get_cx_xform(Xform const & x1, Xform const & x2, int sym){
	Vec cen = get_rot_center(x1,x2,sym);
	return Xform( x2.R, -cen );
}



inline
Real
get_rmsd(
	utility::vector1<Vec> const & native_ca,
	utility::vector1<Vec> const & init_ca,
	Hit const & h//,
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

	Mat Rz = z_rotation_matrix_degrees(180.0);
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


void
dock(
	Pose const & init_pose,
	std::string const & fn,
	basic::sampling::orientations::QuaternionGridCOP qgrid,
	Pose const & /*native_olig*/,
	int native_nfold = 1,
	Vec cen=Vec(0,0,0)
){
	TR << "dock " << fn << endl;

	using namespace basic::options;
	using namespace OptionKeys;
	using namespace core::scoring::constraints;
	using core::id::AtomID;

	bool bench = native_nfold != 1;


	core::id::AtomID_Map<Real> clashmap;
	core::pose::initialize_atomid_map(clashmap,init_pose,-1.0);
	for(Size ir = 1; ir <= init_pose.n_residue(); ++ir){
		if(!init_pose.residue(ir).is_protein()) continue;
		clashmap[AtomID(1,ir)] = init_pose.residue(ir).atom_type(1).lj_radius();
		clashmap[AtomID(2,ir)] = init_pose.residue(ir).atom_type(2).lj_radius();
		clashmap[AtomID(3,ir)] = init_pose.residue(ir).atom_type(3).lj_radius();
		// clashmap[AtomID(4,ir)] = init_pose.residue(ir).atom_type(4).lj_radius();
		clashmap[AtomID(5,ir)] = init_pose.residue(ir).atom_type(5).lj_radius();
	}
	protocols::sic_dock::SICFast sic(option[sicdock::clash_dis]());
	sic.init(init_pose,clashmap);

	protocols::sic_dock::scores::MotifHashRigidScoreOP mhscore_(NULL);
	protocols::sic_dock::scores::TrisBpyScoreCOP bpy3_sc_(NULL);
	protocols::sic_dock::JointScoreOP rigidsfxn = new protocols::sic_dock::JointScore;
	protocols::sic_dock::RigidScoreCOP cbscore_  = new protocols::sic_dock::CBScore(init_pose,init_pose,option[sicdock::clash_dis](),option[sicdock::contact_dis]());
	// rigidsfxn->add_score(cbscore_,1.0);

	if(option[basic::options::OptionKeys::sicdock::trisbpy].user()){
		bpy3_sc_ = new protocols::sic_dock::scores::TrisBpyScore(init_pose,option[sicdock::trisbpy](),option[sicdock::trisbpy_min_contacts]());
		rigidsfxn->add_score(bpy3_sc_,1.0);
	}

	if(option[basic::options::OptionKeys::mh::xform_score_data].user()){
		mhscore_ = new protocols::sic_dock::scores::MotifHashRigidScore(init_pose,init_pose);
		rigidsfxn->add_score(mhscore_, option[sicdock::motif_hash_weight]() );
	}

	utility::vector1<Vec> init_ca;
	for(Size ir = 1; ir <= init_pose.n_residue(); ++ir){
		if(init_pose.residue(ir).has("CA")){
			init_ca.push_back(init_pose.residue(ir).xyz("CA"));
		}
	}

	vector1<int> syms;
	if(!bench) syms = option[sicdock::nfold]();
	else syms.push_back(native_nfold);
	std::sort(syms.begin(),syms.end());
	vector1<Mat> Rsym(syms.size());
	for(int ic = 1; ic <= (int)syms.size(); ic++) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/Real(syms[ic]));
	vector1<Hit> hits;
	Hit tophit;

	Size totsamp = qgrid->num_samples();
	Size outinterval = min((Size)100000,totsamp/10);

	Size nsamp = 0, nsamp_motif=0, nsamp_slide=0;
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for(int iqg = 1; iqg <= (int)totsamp; ++iqg){
		if(option[sicdock::iori].user() && std::find( option[sicdock::iori]().begin(), option[sicdock::iori]().end(),iqg) ==option[sicdock::iori]().end()) continue;

		if(iqg%outinterval==0){
			#ifdef USE_OPENMP
				#pragma omp critical
			#endif
			{
				cout << F(5,1,(Real)iqg/(Real)totsamp*100.0) << "%" << ", nsamp: " << F(7,4,(Real)nsamp/1000000.0) << "M, Nsic: " << F(7,4,(Real)nsamp_slide/1000000.0) << "M ";
				if(hits.size()){ cout<<"   topscore: " <<F(8,2,tophit.rscore); rigidsfxn->show(cout,tophit.x1,tophit.x2); }
				cout << F(7,3,(Real)mhscore_->nhashlookups()/1000000.0) <<  endl;
			}
		}

		Mat R = qgrid->quaternion(iqg).rotation_matrix();

		Real dummy;
		if( option[sicdock::rotate_x_only]() && fabs(rotation_axis(R,dummy).x()) < 0.999 ) continue;

		Xform const x1_0( R, Vec(0,0,0));
		for(int ic = 1; ic <= (int)syms.size(); ic++){
			Xform const x2( Rsym[ic]*R, Vec(0,0,0) );

			Xform x1(x1_0);
			Real const t = sic.slide_into_contact(x1,x2,Vec(0,0,1));
			++nsamp_slide;
			x1.t.z() += t;

			Real const backoff_delta = min(1.0,max(0.1547933,qgrid->maxrad() / 180.0 * 3.14159 * fabs(t/4.0)));
			Size const nang = (Size)(360.0/syms[ic]/1.5/qgrid->maxrad()+1);
			Real const ang_delta = 360.0/(Real)syms[ic]/(Real)nang;

			for(Real sdelta = 0; sdelta <= option[sicdock::sample_thickness](); sdelta += backoff_delta){
				++nsamp;
				x1.t.z() += backoff_delta;

				Hit h(iqg,0.0,syms[ic]);
				h.backoff = sdelta;
				h.x1=x1; h.x2=x2;
				h.x1 = rotation_matrix_degrees(Vec(0,1,0),90.0) * get_cx_xform(x1,x2,h.sym);
				h.x1 = z_rotation_matrix( atan2(h.x1.t.x(),h.x1.t.y()) ) * h.x1;
				h.x2 = h.x1 * ~x1 * x2;

				Real rscore = ((protocols::sic_dock::RigidScoreOP)rigidsfxn)->score(h.x1,h.x2);
				if( rscore < option[sicdock::min_contact_score_cut]() ) break;
				if( rscore < option[sicdock::min_score            ]() ) continue;
				h.rscore = rscore;

				if(!option[sicdock::dihedral]()){
					#ifdef USE_OPENMP
						#pragma omp critical
					#endif
					hits.push_back(h);
					if(h.rscore>tophit.rscore) tophit = h;
					continue;
				}

				Hit const h0 = h;
				for(Real dang = 180.0/(Real)syms[ic]; dang < 540.0/(Real)syms[ic]-0.0001; dang += ang_delta){
					h = h0;

					Xforms x1s(syms[ic]),x2s(syms[ic]);
					x1s[1] = z_rotation_matrix_degrees(dang/2.0) * h0.x1;
					x1s[2] = z_rotation_matrix_degrees(dang/2.0) * h0.x2;
					x2s[1] = y_rotation_matrix_degrees(   180.0) * x1s[1];
					x2s[2] = y_rotation_matrix_degrees(   180.0) * x1s[2];
					// for(Size isym = 1; isym <= syms[ic]; ++isym){
					// 	x1s[sym] = z_rotation_matrix_degrees(360.0/(Real)syms[ic]*(Real)(isym-1)+dang/2.0) * h0.x1;
					// 	x2s[sym] = y_rotation_matrix_degrees(180.0) * x1s[isym];
					// }

					Real td = sic.slide_into_contact(x1s,x2s,Vec(0,0,1));
					++nsamp_slide;
					x1s[1].t.z() += td/2.0; x1s[2].t.z() += td/2.0;
					x2s[1].t.z() -= td/2.0; x2s[2].t.z() -= td/2.0;

					for(Real ddelta = 0; ddelta <= option[sicdock::sample_thickness](); ddelta += backoff_delta){
						x1s[1].t.z() += backoff_delta/2.0; x1s[2].t.z() += backoff_delta/2.0;
						x2s[1].t.z() -= backoff_delta/2.0; x2s[2].t.z() -= backoff_delta/2.0;

						rscore = h0.rscore + rigidsfxn->score(x1s[1],x2s);
						if( rscore < 2.0*option[sicdock::min_contact_score_cut]() ) break;
						if( rscore < 2.0*option[sicdock::min_score            ]() ) continue;

						++nsamp;
						h.rscore = rscore;
						h.dang = dang;
						h.ddis = ddelta;
						h.x1=x1s[1]; h.x2=x1s[2]; h.x3=x2s[1]; h.x4=x2s[2];

						#ifdef USE_OPENMP
							#pragma omp critical
						#endif
						if( h.rscore >= option[sicdock::min_score]() ){
							hits.push_back(h);
							if(h.rscore>tophit.rscore) tophit = h;
							// if(uniform() <= 0.001){
							// 	// cout << td <<' '<< t13 <<' '<< t14 <<' '<< t23 <<' '<< t24 << endl;
							// 	// cout << h.x1 << endl;
							// 	// cout << h.x2 << endl;
							// 	// cout << h.x3 << endl;
							// 	// cout << h.x4 << endl;
							// 	Pose tmp1(init_pose),tmp2(init_pose),tmp3(init_pose),tmp4(init_pose);
							// 	protocols::sic_dock::xform_pose(tmp1,h.x1);
							// 	protocols::sic_dock::xform_pose(tmp2,h.x2);
							// 	protocols::sic_dock::xform_pose(tmp3,h.x3);
							// 	protocols::sic_dock::xform_pose(tmp4,h.x4);
							// 	tmp1.dump_pdb("tmp1.pdb");
							// 	tmp2.dump_pdb("tmp2.pdb");
							// 	tmp3.dump_pdb("tmp3.pdb");
							// 	tmp4.dump_pdb("tmp4.pdb");
							// 	utility_exit_with_message("rastrast");
							// }
						}
					}
				}
			}
		}
	}

	cout << "num samples:       " << (Real)nsamp/1000000.0 << "M" << endl;
	cout << "num samples motif: " << (Real)nsamp_motif/1000.0 << "K" << endl;

	Size nstruct = basic::options::option[basic::options::OptionKeys::out::nstruct]();
	std::string outdir = "./";
	if(basic::options::option[basic::options::OptionKeys::out::file::o].user()){
		outdir = basic::options::option[basic::options::OptionKeys::out::file::o]()+"/";
	}
		bool NAT=false;
		vector1<vector1<Hit> > dumpedit(syms.back());
		vector1<Size> order(syms.back(),0);
			TR << " found " << hits.size() << " hits" << endl;
			Size nout=0,i=0;
			cout << "  file                       nf    iori    dang   order   rscore"; if(NAT) cout << "     rms"; cout << "   rd1   rd2";
			cbscore_->show(cout); cout << " ";
			if(bpy3_sc_) bpy3_sc_->show(cout);
			if(mhscore_) mhscore_->show(cout);
			cout << endl;

			Pose native;
			if(NAT){ // native stuff
				Xform x1(Mat::identity(),cen);
				Xform x2(z_rotation_matrix_degrees(360.0/(Real)syms.front()));
				x2.t = x2.R * cen;
				Hit h(-1,rigidsfxn->score(x1,x2),syms.front());
				h.x1 = x1;
				h.x2 = x2;
				cout << "| " << ObjexxFCL::format::LJ(40,"NATIVE") <<' '<< I(3,h.sym) <<' '<< I(7,h.iqg) <<' '<< I(7,qgrid->num_samples())  << "      0 " << F(9,3,h.rscore) << "   0.000   0.000";
				rigidsfxn->show(cout,h.x1,h.x2);
				Pose tmp1(init_pose),tmp2(init_pose);
				protocols::sic_dock::xform_pose(tmp1,h.x1);
				protocols::sic_dock::xform_pose(tmp2,h.x2);
				// if(mhscore_){
				// 	utility::io::ozstream out(outdir+"/native_olig_motifs.pdb");
				// 	int ndump = mhscore_->dump_matching_motifs(tmp1,tmp2,out);
				// 	out.close();
				// }
				make_Cx(tmp1,syms.front());
				native = tmp1;
				// tmp1.dump_pdb(outdir+"/native_olig.pdb");
			}

			std::sort(hits.begin(),hits.end(),cmpscore);
			while( nout < (Size)option[sicdock::nscore]() && ++i <= hits.size() ){
				Hit & h(hits[i]);
				Xforms x1s0,x2s0;
				x1s0.push_back(h.x1);
				x2s0.push_back(h.x2);
				if(option[sicdock::dihedral]()){ x2s0.push_back(h.x3); x2s0.push_back(h.x4); }
				Xforms const &x1s(x1s0),&x2s(x2s0);
				order[h.sym]++;
				bool toosimilar = false;
				for(vector1<Hit>::const_iterator j = dumpedit[h.sym].begin(); j != dumpedit[h.sym].end(); ++j){
					Xform xrel = h.x1 * ~(j->x1);
					Real ang=0.0;
					numeric::rotation_axis(xrel.R,ang);
					ang = numeric::conversions::degrees(ang);
					Real dis = xrel.t.length();
					if( ang < option[sicdock::redundancy_angle_cut]() || dis < option[sicdock::redundancy_dist_cut]() ){
						toosimilar = true;
					}
				}
				if(toosimilar) continue;
				dumpedit[h.sym].push_back(h);

				std::string tag = utility::file_basename(fn);
				if(tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
				if(tag.substr(tag.size()-4)==".pdb") tag = tag.substr(0,tag.size()-4);
				tag += "_C"+ObjexxFCL::string_of(h.sym)+"_"+ObjexxFCL::format::I(4,4,nout+1);
				Pose olig(init_pose);
				protocols::sic_dock::xform_pose(olig,h.x1);
				protocols::sic_dock::make_Cx(olig,h.sym);
				if(option[sicdock::dihedral]()) protocols::sic_dock::make_Cx(olig,2,Vec(0,1,0));
				Real rms = -1;
				if(NAT) rms = core::scoring::CA_rmsd(olig,native);

				using namespace protocols::sic_dock;
				xyzStripeHashPoseCOP clash_checker = NULL;
				if(nout < nstruct){
					if(mhscore_){
						utility::io::ozstream out(outdir+tag+"_motifs.pdb");
						clash_checker = new xyzStripeHashPose(olig,NCO,2.5);
						mhscore_->dump_matching_motifs(x1s,x2s,out,clash_checker);
						out.close();
						// protocols::motif_hash::MotifRotamerSetOperationCOP rso = mhscore_->motif_hash()->get_matching_motifs_rotsetop(olig,NULL,0.6,NULL);
						// rso->dump_pdb("rso.pdb");
					}
					if(bpy3_sc_) bpy3_sc_->dump_trisbpy(x1s,x2s,outdir+tag+"_bpy.pdb");
					olig.dump_pdb(outdir+tag+".pdb");
					// Pose tmp1(init_pose),tmp2(init_pose),tmp3(init_pose),tmp4(init_pose);
					// protocols::sic_dock::xform_pose(tmp1,x1s[1]);
					// protocols::sic_dock::xform_pose(tmp2,x2s[1]);
					// protocols::sic_dock::xform_pose(tmp3,x2s[2]);
					// protocols::sic_dock::xform_pose(tmp4,x2s[3]);
					// tmp1.dump_pdb("tmp1.pdb");
					// tmp2.dump_pdb("tmp2.pdb");
					// tmp3.dump_pdb("tmp3.pdb");
					// tmp4.dump_pdb("tmp4.pdb");
					// utility_exit_with_message("rastrast");
				}

				cout << "| " << ObjexxFCL::format::LJ(25,tag) <<' ';
				cout << I(3,h.sym) <<' ';
				cout << I(7,h.iqg) <<' ';
				cout << F(7,3,h.dang) <<' ';
				cout << I(6,order[h.sym]) <<' ';
				cout << F(9,3,h.rscore) <<' ';
				if(0) cout << F(7,3,rms) <<' ';
				cout << F(5,3,h.backoff) <<' ';
				cout << F(5,3,h.ddis);;
				cbscore_->show(cout,x1s,x2s); cout << " ";
				if(bpy3_sc_) bpy3_sc_->show(cout,x1s,x2s);
				if(mhscore_) mhscore_->show(cout,x1s,x2s,clash_checker,10);
				cout << endl;
				++nout;
			}
}

utility::vector1<Vec>
align_native_state(
	core::pose::Pose & pose,
	int nfold
){
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
	protocols::sic_dock::rot_pose( pose, Vec(0,0,1), -dihedral_degrees(Vec(1,0,0),Vec(0,0,0),Vec(0,0,1),cen) );
	protocols::sic_dock::rot_pose( pose, Vec(1,0,1), 180.0 ); // symZ to symX
	utility::vector1<Vec> CAs;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		if(pose.residue(ir).has("CA")){
			CAs.push_back(pose.residue(ir).xyz("CA"));
		}
	}
	return CAs;
}

void
get_tasks_from_command_line(
	utility::vector1<std::pair<string,int> > & tasks
){
	using namespace basic::options;
	using namespace OptionKeys;
	if(option[in::file::s   ].user()) for(Size ifn = 1; ifn <= option[in::file::s   ]().size(); ++ifn) tasks.push_back(std::make_pair(option[in::file::s   ]()[ifn],1));
	if(option[sicdock::bench2].user()) for(Size ifn = 1; ifn <= option[sicdock::bench2]().size(); ++ifn) tasks.push_back(std::make_pair(option[sicdock::bench2]()[ifn],2));
	if(option[sicdock::bench3].user()) for(Size ifn = 1; ifn <= option[sicdock::bench3]().size(); ++ifn) tasks.push_back(std::make_pair(option[sicdock::bench3]()[ifn],3));
	if(option[sicdock::bench4].user()) for(Size ifn = 1; ifn <= option[sicdock::bench4]().size(); ++ifn) tasks.push_back(std::make_pair(option[sicdock::bench4]()[ifn],4));
	if(option[sicdock::bench5].user()) for(Size ifn = 1; ifn <= option[sicdock::bench5]().size(); ++ifn) tasks.push_back(std::make_pair(option[sicdock::bench5]()[ifn],5));
	if(option[sicdock::bench6].user()) for(Size ifn = 1; ifn <= option[sicdock::bench6]().size(); ++ifn) tasks.push_back(std::make_pair(option[sicdock::bench6]()[ifn],6));
}

int main(int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;

	std::cout << "reading quaternion grid data" << endl;
	basic::sampling::orientations::QuaternionGridCOP qgrid; {
		basic::sampling::orientations::QuaternionGridManager *m = basic::sampling::orientations::QuaternionGridManager::get_instance();
		if( option[sicdock::sample_number].user() && option[sicdock::sample_radius].user() ){
			utility_exit_with_message("cant specify both -sicdock::sample_radius and -sicdock::sample_number");
		} else if( option[sicdock::sample_number].user() ){
			qgrid = m->request_by_size(option[sicdock::sample_number]());
		} else if( option[sicdock::sample_radius].user() ){
			qgrid = m->request_by_radius(option[sicdock::sample_radius]());
		} else {
			qgrid = m->request_by_radius(option[sicdock::sample_radius]());
		}
		std::cout << "CXdock sampling: " <<  *qgrid << endl;
	}

	vector1<std::pair<string,int> > tasks;
	get_tasks_from_command_line(tasks);
	cout << "tasks: " << tasks.size() << endl;

	for(vector1<std::pair<string,int> >::const_iterator i = tasks.begin(); i != tasks.end(); ++i){
		string fn = i->first;
		int nfold = i->second;
		TR << "searching " << fn <<' '<< nfold << std::endl;

		Pose pnat;
		core::import_pose::pose_from_pdb(pnat,fn);
		if(pnat.n_residue()>300){ cout << fn << " is too big" << endl; continue; }
		core::scoring::dssp::Dssp dssp(pnat);
		dssp.insert_ss_into_pose(pnat);

		Vec const cen = protocols::sic_dock::center_of_geom(pnat,1,pnat.n_residue());
		if(!option[sicdock::no_recenter]()) protocols::sic_dock::trans_pose(pnat,Vec(0,0,-cen.z()));

		Pose olig(pnat); make_Cx(olig,nfold);
		// olig.dump_pdb("native_olig.pdb");
		// vector1<Vec> native_ca = align_native_state(pnat,nfold); // align CA com

		if(!option[sicdock::no_recenter]()) protocols::sic_dock::trans_pose(pnat,-Vec(cen.x(),cen.y(),0));
		if( pnat.n_residue() > (Size)option[sicdock::max_res]() ) continue;

		dock(pnat,fn,qgrid,olig,nfold,cen);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	cout << "SICDOCK_DONE" << endl;
}
