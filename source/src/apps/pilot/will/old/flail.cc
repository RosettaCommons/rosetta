// includes ///////////////////////////////////
	#include <devel/init.hh>
	#include <numeric/xyzTransform.io.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/random/random.hh>
	#include <protocols/sic_dock/Rose.hh>
	#include <core/pose/Pose.hh>
	#include <core/conformation/Residue.hh>
	#include <core/scoring/Ramachandran.hh>
	#include <sstream>
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/option.hh>
	#include <basic/options/option_macros.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/pose/annotated_sequence.hh>
	#include <platform/types.hh>
	#include <utility/io/ozstream.hh>
// NAMEPACES
	#ifdef USE_OPENMP
	#include <omp.h>
	#endif
	using std::string;
	using utility::vector1;
	using ObjexxFCL::format::I;
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::RJ;
	using ObjexxFCL::string_of;
	using numeric::min;
	using numeric::max;
	using std::cout;
	using std::cerr;
	using std::endl;
	using namespace numeric;
	using namespace numeric::random;
	using core::pose::Pose;
	using core::pose::PoseOP;
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using protocols::sic_dock::Rose;
	typedef xyzVector   <core::Real> V;
	typedef xyzMatrix   <core::Real> M;
	typedef xyzTransform<core::Real> X;
	typedef core::id::AtomID AID;
// Math & constants ////////////////////////////////////
	core::Real const eps = 0.000001;
	static inline float
	fastpow2 (float p)
	{
	  float offset = (p < 0) ? 1.0f : 0.0f;
	  float clipp = (p < -126) ? -126.0f : p;
	  int w = (int)clipp;
	  float z = clipp - w + offset;
	  union { uint32_t i; float f; } v = { (1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z) };
	  return v.f;
	}

	static inline float
	fastexp(float p)
	{
	  return fastpow2 (1.442695040f * p);
	}
// OPTS
	OPT_1GRP_KEY( StringVector  , flail, linker_seqs )
	OPT_1GRP_KEY( IntegerVector , flail, linker_lens )
	OPT_1GRP_KEY( IntegerVector , flail, target_residues )
	// OPT_1GRP_KEY( core::Real          , flail, magnitude )
	OPT_1GRP_KEY( Real          , flail, temperature )
	OPT_1GRP_KEY( Real          , flail, output_frac )
	OPT_1GRP_KEY( Real          , flail, move_sd )
	OPT_1GRP_KEY( Real          , flail, walltime_limit )

	void register_options() {
		NEW_OPT( flail::linker_seqs, "linker seqs in order", "" );
		NEW_OPT( flail::linker_lens, "linker lengths, assume GGS in order", 0 );
		NEW_OPT( flail::target_residues, "target residues specifying dimer interaction", 0 );
		// NEW_OPT( flail::magnitude, "perturbation magnitude", 10.0 );
		NEW_OPT( flail::temperature, "temperature", 10.0 );
		NEW_OPT( flail::output_frac, "output_frac", 0.001 );
		NEW_OPT( flail::move_sd, "move_sd", 60.0 );
		NEW_OPT( flail::walltime_limit, "time limit", 9e9 );
	}
// classes
	struct Link {
		typedef utility::vector1<V> PTS;
		Link( std::string const & seq ) : N(seq.size()*3),last_move_res(2),last_move_phi(-90.0),last_move_psi(90.0) {
			P.reserve(N);
			Pose p;
			make_pose_from_sequence(p,seq,core::chemical::CENTROID);
			for(core::Size ir = 1; ir <= p.size(); ++ir){
				p.set_phi(ir,-90.0);
				p.set_psi(ir, 90.0);
				p.set_omega(ir,180.0);
			}
			for(core::Size ir = 1; ir <= p.size(); ++ir){
				P.push_back( p.residue(ir).xyz( "N") );
				P.push_back( p.residue(ir).xyz("CA") );
				P.push_back( p.residue(ir).xyz( "C") );
				aas.push_back(p.residue(ir).aa());
			}
		}
		void undo_most_recent_move(){
			set_phi(last_move_res,last_move_phi);
			set_psi(last_move_res,last_move_psi);
			last_move_res = 2;
			last_move_phi = -90.0;
			last_move_psi =  90.0;
		}
		void set_phi_psi(core::Size i, core::Real phi, core::Real psi){
			last_move_res = i;
			last_move_phi = phi;
			last_move_psi = psi;
			set_phi(i,phi);
			set_psi(i,psi);
		}
		void dump_pdb(std::ostream & out, char chain='L') const {
			std::string na[3] = {"  N "," CA ","  C "};
			for(core::Size i = 1; i <= N; ++i){
				V const & v(P[i]);
				core::Size ir = (i-1)/3+1;
				out<<"ATOM  "<<I(5,i)<<' '<<na[(i-1)%3]<<' '<<"GLY"<<' '<<chain<<I(4,ir)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			}
		}
		void dump_pdb(std::string const & fn, char chain='L') const {
			utility::io::ozstream out(fn);
			dump_pdb(out,chain);
			out.close();
		}
		X const n_anchor(){
			return X( P[ 1 ], P[ 2 ], P[3] );
		}
		X const c_anchor(){
			return X( P[N-2], P[N-1], P[N] );
		}
		void align_n(X const & anchor){
			X const x = anchor * ~n_anchor();
			for(PTS::iterator i = P.begin(); i != P.end(); ++i) *i = x * *i;
			assert(is_connected());
		}
		void align_c(X const & anchor){
			X const x = anchor * ~c_anchor();
			for(PTS::iterator i = P.begin(); i != P.end(); ++i) *i = x * *i;
			assert(is_connected());
		}
		bool clashes() const {
			// for(PTS::const_iterator i = P.begin()+3; i != P.end(); ++i) if( i->distance_squared(*(i-3)) < 6.0 ) return true;
			// for(PTS::const_iterator i = P.begin()+4; i != P.end(); ++i) if( i->distance_squared(*(i-4)) < 7.0 ) return true;
			for(PTS::const_iterator i = P.begin(); i != P.end()-12; ++i)
				for(PTS::const_iterator j = i+11; j != P.end(); ++j)
					if( i->distance_squared(*j) < 7.0 )
						return true;
			return false;
		}
		bool clashes(Link const & o) const {
			for(PTS::const_iterator i = o.P.begin()+5; i != o.P.end()-5; ++i){
				for(PTS::const_iterator j = P.begin()+5; j != P.end()-5; ++j){
					if( i->distance_squared(*j) < 7.0 )
						return true;
				}
			}
			return false;
		}
		bool clashes(Rose const & rose) const {
			if( N > 16 ){
				// TODO would be more efficient to xform once
				X xinv = ~rose.x;
				for(PTS::const_iterator i = P.begin()+12; i != P.end()-12; ++i){
					if( rose.h->clash(xinv * *i) ) return true;
				}
			}
			return false;
		}
		core::Size nres() const { return N/3; }
		bool is_connected() const {
			for(PTS::const_iterator i = P.begin()+1; i != P.end(); ++i){
				core::Real d2 = i->distance_squared(*(i-1));
				if( fabs(d2-1.458001*1.458001) > 0.00001 &&
					fabs(d2-1.523259*1.523259) > 0.00001 &&
					fabs(d2-1.328685*1.328685) > 0.00001
				) return false;
			}
			return true;
		}
		Real  get_phi(core::Size i){ return dihedral_degrees(P[3*i-3],P[3*i-2],P[3*i-1],P[3*i]); }
		Real  get_psi(core::Size i){ return dihedral_degrees(P[3*i-2],P[3*i-1],P[3*i],P[3*i+1]); }
		void move_phi(core::Size i, core::Real deg){ rotate(3*i-2,deg); }
		void move_psi(core::Size i, core::Real deg){ rotate(3*i-1,deg); }
	private:
		void rotate(core::Size const & ibond, core::Real const & deg){
			if( ibond < 1 || ibond+1 > N ) utility_exit_with_message("ibond out of bounds!!");
			V const cen(P[ibond]);
			V const axs(P[ibond+1]-P[ibond]);
			X const rot( X::rot_deg(axs,deg,cen) );
			for(PTS::iterator i = P.begin()+ibond+1; i != P.end(); ++i) *i = rot * *i;
			assert(is_connected());
		}
		void  set_phi(core::Size i, core::Real deg){ rotate(3*i-2,deg-dihedral_degrees(P[3*i-3],P[3*i-2],P[3*i-1],P[3*i])); }
		void  set_psi(core::Size i, core::Real deg){ rotate(3*i-1,deg-dihedral_degrees(P[3*i-2],P[3*i-1],P[3*i],P[3*i+1])); }
		core::Size N;
		PTS P;
		core::Size last_move_res;
		core::Real last_move_phi,last_move_psi;
	public:
		vector1<core::chemical::AA> aas;
	};
	struct Bola {
		vector1<Rose> roses;
		vector1<Link> links;
		vector1<Link>::iterator last_moved_link;
		X target,target_stub1,target_stub2;
		core::Size target_rose1,target_rose2;
		core::Size ntrials,nclashes,naccepts;
		V T11,T12,T13,T21,T22,T23;
		core::scoring::RamachandranCOP rama;
		utility::vector1<Size> clash_intralink_,clash_crosslink_,clash_roselink_,clash_roserose_;
		Real move_sd_;
		bool started_;
		vector1<xyzVector<Real> > tgtatoms_;
		Bola(
			vector1<PoseOP> const & p,
			vector1<string> const & l,
			vector1<core::Size> const & target_residues,
			core::scoring::RamachandranCOP _rama
		):
		   ntrials(0),
		   nclashes(0),
		   naccepts(0),
		   rama(_rama),
		   clash_intralink_(l.size(),0),clash_crosslink_(l.size(),0),clash_roselink_(p.size(),0),clash_roserose_(p.size(),0),
		   move_sd_(option[flail::move_sd]()),
		   started_(false)
		{
			tgtatoms_.push_back(p[1]->xyz(core::id::AtomID(2,24)));
			tgtatoms_.push_back(p[1]->xyz(core::id::AtomID(2,36)));
			tgtatoms_.push_back(p[1]->xyz(core::id::AtomID(2,53)));
			tgtatoms_.push_back(p[1]->xyz(core::id::AtomID(2,59)));
			tgtatoms_.push_back(p[1]->xyz(core::id::AtomID(2,64)));
			tgtatoms_.push_back(p[1]->xyz(core::id::AtomID(2,82)));
			tgtatoms_.push_back(p[1]->xyz(core::id::AtomID(2,88)));
			for(vector1<string>::const_iterator i = l.begin(); i != l.end(); ++i) links.push_back(Link(*i));
			for(vector1<PoseOP>::const_iterator i = p.begin(); i != p.end(); ++i) roses.push_back(Rose(*i));
			last_moved_link = links.end();
			if(p.size()!=l.size()+1) utility_exit_with_message("number of poses and linkers unmatches");
			if(p.size()!=target_residues.size()) utility_exit_with_message("number of targets and linkers unmatches");
			target_rose1=0,target_rose2=0;
			for(core::Size i = 1; i <= target_residues.size(); ++i){
				if(target_residues[i]==0) continue;
				if(target_rose2!=0) utility_exit_with_message("two targets only!");
				if(target_rose1==0){
					core::conformation::Residue const & res(p[i]->residue(target_residues[i]));
					target_rose1 = i;
					target_stub1 = X( res.xyz("CB"),res.xyz("CA"),res.xyz("N") );
					T11 = res.xyz("CB");
					T12 = res.xyz("CA");
					T13 = res.xyz("N");
				} else {
					core::conformation::Residue const & res(p[i]->residue(target_residues[i]));
					target_rose2 = i;
					target_stub2 = X( res.xyz("CB"),res.xyz("CA"),res.xyz("N") );
					T21 = res.xyz("CB");
					T22 = res.xyz("CA");
					T23 = res.xyz("N");
				}
			}
			if(target_rose1==0||target_rose2==0) utility_exit_with_message("not enough targets!");
			target = get_target_xform();
			refold();
			// cout << "TARGETS " << target_rose1 << " " << target_rose2 << " start dis: " << get_target_distance() << endl;
		}
		inline X get_target_xform(){
			return (roses[target_rose2].x*target_stub2) - (roses[target_rose1].x*target_stub1);
		}
		// inline X get_target_xform_safe(){
			// 	X const & x1(roses[target_rose1].x);
			// 	X const & x2(roses[target_rose2].x);
			// 	V xT11 = x1*T11;
			// 	V xT12 = x1*T12;
			// 	V xT13 = x1*T13;
			// 	V xT21 = x2*T21;
			// 	V xT22 = x2*T22;
			// 	V xT23 = x2*T23;
			// 	X xT1(xT11,xT12,xT13);
			// 	X xT2(xT21,xT22,xT23);
			// 	if( xT1.distance_squared( x1*target_stub1 ) > 0.001 ) utility_exit_with_message("bad xform-based 1");
			// 	if( xT2.distance_squared( x2*target_stub2 ) > 0.001 ) utility_exit_with_message("bad xform-based 2");

			// 	V translation = xT1.R.transposed() * ( xT2.t - xT1.t );
			// 	M rotation = xT1.R.transposed() * xT2.R;

			// 	X x3(rotation,translation);
			// 	if( x3.distance(get_target_xform()) > 0.00001 ) utility_exit_with_message("get_target_xform vs safe fail");

			// 	return  x3;
			// }
		inline core::Real get_target_distance(){
			Real dist2 = 0.0;
			dist2 += (roses[1].x * tgtatoms_[1]).distance_squared( tgtatoms_[1] );
			dist2 += (roses[1].x * tgtatoms_[2]).distance_squared( tgtatoms_[2] );
			dist2 += (roses[1].x * tgtatoms_[3]).distance_squared( tgtatoms_[3] );
			dist2 += (roses[1].x * tgtatoms_[4]).distance_squared( tgtatoms_[4] );
			dist2 += (roses[1].x * tgtatoms_[5]).distance_squared( tgtatoms_[5] );
			dist2 += (roses[1].x * tgtatoms_[6]).distance_squared( tgtatoms_[6] );
			dist2 += (roses[1].x * tgtatoms_[7]).distance_squared( tgtatoms_[7] );
			return sqrt(dist2/7.0);
			// return 10.0*(get_target_xform().distance( target ) - 0.9*get_target_xform().t.distance( target.t ));
		}
		void random_close_phipsi_from_rama(Real sd, core::chemical::AA const aa, Real const phi0, Real const psi0, Real & phi, Real & psi){
			for(int i = 0; i < 1000; ++i){
				rama->random_phipsi_from_rama(aa,phi,psi);
				Real const dphi = min(fabs(phi-phi0),360.0-fabs(phi-phi0));
				Real const dpsi = min(fabs(psi-psi0),360.0-fabs(psi-psi0));
				Real const P = fastexp(-dphi/sd*dphi/sd) * fastexp(-dpsi/sd*dpsi/sd);
				if( uniform() < P ){
					// cout << "random_close_phipsi_from_rama " << sd << " " << P << " " << phi0 << " " << psi0 << " " << phi << " " << psi << " " << dphi << " " << dpsi << endl;
					return;
				}
			}
			cout << "random_close_phipsi_from_rama FAIL" << endl;
		}
		void random_move(core::chemical::AA & moved_aa, Real & moved_phi, Real & moved_psi){
			vector1<std::pair<vector1<Link>::iterator,core::Size> > choices;
			for(vector1<Link>::iterator i=links.begin(); i!=links.end(); ++i){
				for(core::Size j = 2; j < i->nres(); ++j) choices.push_back(std::make_pair(i,j));
			}
			std::pair<vector1<Link>::iterator,core::Size> choice = numeric::random::random_element(choices);
			vector1<Link>::iterator link = choice.first;
			core::Size ir = choice.second;
			// core::Real phi=0,psi=0;
			// // rama->random_phipsi_from_rama(link->aas[ir],phi,psi);
			// random_close_phipsi_from_rama(started_?move_sd_:9999999.9,link->aas[ir],phi0,psi0,phi,psi);
			bool phiorpsi = uniform() < 0.5;
			Real gauss = (uniform()+uniform()+uniform()+uniform()+uniform()+uniform()+uniform()+uniform()+uniform()-4.5)/3.0;
			if(phiorpsi){
				link->move_phi(ir,gauss*move_sd_);
			} else {
				link->move_psi(ir,gauss*move_sd_);
			}
			moved_aa  = link->aas[ir];
			moved_phi = link->get_phi(ir);
			moved_psi = link->get_psi(ir);
			last_moved_link = link;
			++ntrials;
		}

		void random_nonclashing_move(core::chemical::AA & moved_aa, Real & moved_phi, Real & moved_psi){
			// for(int i=1;i<=clash_intralink_.size();++i) clash_intralink_[i]=0;
			// for(int i=1;i<=clash_crosslink_.size();++i) clash_crosslink_[i]=0;
			// for(int i=1;i<=clash_roselink_ .size();++i) clash_roselink_ [i]=0;
			// for(int i=1;i<=clash_roserose_ .size();++i) clash_roserose_ [i]=0;
			for(core::Size i = 0; i < started_?1000:1000000; ++i){
				random_move(moved_aa,moved_phi,moved_psi);
				if(!clashes()){
					started_ = true;
					// if(i < 100 && move_sd_ < 360.0){
					// 	move_sd_ *= 1.01;
					// }
					// cout << "made nonclash move" << endl;
					// dump_pdb("nonclash_move.pdb");
					// utility_exit_with_message("test");

					return;
				}
				++nclashes;
				undo_most_recent_move();
			}
			// dump_pdb("stuck_clash.pdb");
			// cerr << "clash_intralink_: "; for(int i=1;i<=clash_intralink_.size();++i) << clash_intralink_[i] << " ";
			// cerr << "clash_crosslink_: "; for(int i=1;i<=clash_crosslink_.size();++i) << clash_crosslink_[i] << " ";
			// cerr << "clash_roselink_:  "; for(int i=1;i<=clash_roselink_ .size();++i) << clash_roselink_ [i] << " ";
			// cerr << "clash_roserose_:  "; for(int i=1;i<=clash_roserose_ .size();++i) << clash_roserose_ [i] << endl;
			// utility_exit_with_message("no non-clashing after 1000000 tries!");
			// if(move_sd_>5.0){
			// 	move_sd_ *= 0.9;
			// 	cout << "failed to find non-clashing move in 10K trials, reducing move_sd " << move_sd_ << " " << naccepts << endl;
			// } else {
			// 	cout << "failed to find non-clashing move in 10K trials, move_sd at min " << move_sd_ << " " << naccepts << endl;
			// }
			if(started_) dump_pdb(option[out::file::o]()+"/nonclash_move_fail.pdb");
			utility_exit_with_message("failed to find non-clashing move in 1000 trials");

		}
		void undo_most_recent_move(){
			if(last_moved_link != links.end()){
				last_moved_link->undo_most_recent_move();
			} else {
				utility_exit_with_message("UNDO ERROR!!!!!");
			}
			last_moved_link = links.end();
		}
		void refold(){
			for(core::Size i = 1; i <= links.size(); ++i){
				links[links.size()-i+1].align_c( roses[links.size()-i+2].n_anchor() );
				roses[links.size()-i+1].align_c( links[links.size()-i+1].n_anchor() );
			}
		}
		void dump_pdb(std::string const & fn){
			utility::io::ozstream out(fn);
			refold();
			for(core::Size i=1; i<=links.size(); ++i){
				roses[i].dump_minimal_pdb(out,char(63+2*i));
				links[i].dump_pdb(out,char(63+2*i+1));
			}
			// roses.back().dump_minimal_pdb(out,char(63+2*links.size()+2));
			V v;
			// v=roses[target_rose2].x*T23; out<<"ATOM  "<<I(5,10006)<<' '<<"  N "<<' '<<"TS2"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			// v=roses[target_rose2].x*T22; out<<"ATOM  "<<I(5,10005)<<' '<<" CA "<<' '<<"TS2"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			// v=roses[target_rose2].x*T21; out<<"ATOM  "<<I(5,10004)<<' '<<" CB "<<' '<<"TS2"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			// v=roses[target_rose1].x*T13; out<<"ATOM  "<<I(5,10003)<<' '<<"  N "<<' '<<"TS1"<<' '<<'T'<<I(4,998)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			// v=roses[target_rose1].x*T12; out<<"ATOM  "<<I(5,10002)<<' '<<" CA "<<' '<<"TS1"<<' '<<'T'<<I(4,998)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			// v=roses[target_rose1].x*T11; out<<"ATOM  "<<I(5,10001)<<' '<<" CB "<<' '<<"TS1"<<' '<<'T'<<I(4,998)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=           tgtatoms_[1]; out<<"ATOM  "<<I(5,10001)<<' '<<"TGT0"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=           tgtatoms_[2]; out<<"ATOM  "<<I(5,10002)<<' '<<"TGT0"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=           tgtatoms_[3]; out<<"ATOM  "<<I(5,10003)<<' '<<"TGT0"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=           tgtatoms_[4]; out<<"ATOM  "<<I(5,10004)<<' '<<"TGT0"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=           tgtatoms_[5]; out<<"ATOM  "<<I(5,10005)<<' '<<"TGT0"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=           tgtatoms_[6]; out<<"ATOM  "<<I(5,10006)<<' '<<"TGT0"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=           tgtatoms_[7]; out<<"ATOM  "<<I(5,10007)<<' '<<"TGT0"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=roses[1].x*tgtatoms_[1]; out<<"ATOM  "<<I(5,10001)<<' '<<"TGT1"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=roses[1].x*tgtatoms_[2]; out<<"ATOM  "<<I(5,10002)<<' '<<"TGT1"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=roses[1].x*tgtatoms_[3]; out<<"ATOM  "<<I(5,10003)<<' '<<"TGT1"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=roses[1].x*tgtatoms_[4]; out<<"ATOM  "<<I(5,10004)<<' '<<"TGT1"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=roses[1].x*tgtatoms_[5]; out<<"ATOM  "<<I(5,10005)<<' '<<"TGT1"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=roses[1].x*tgtatoms_[6]; out<<"ATOM  "<<I(5,10006)<<' '<<"TGT1"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			v=roses[1].x*tgtatoms_[7]; out<<"ATOM  "<<I(5,10007)<<' '<<"TGT1"<<' '<<"TGT"<<' '<<'T'<<I(4,999)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
			out.close();
			// roses[1].pose()->dump_pdb("results/test.pdb");
			// utility_exit_with_message("DBG");
		}
		bool clashes(){
			bool clash = false;
			int ii=1;
			for(vector1<Link>::const_iterator i=links.begin(); i!=links.end(); ++i,++ii) if(i->clashes()){
				++clash_intralink_[ii];
				return true; //clash = true;
			}
			refold();
			ii=1;
			for(vector1<Link>::const_iterator i=links.begin(); i!=links.end(); ++i,++ii){
				for(vector1<Link>::const_iterator j=i+1; j!=links.end(); ++j){
					if( j->clashes(*i) ){
						++clash_crosslink_[ii];
						return true; //clash = true;
					}
				}
			}
			ii=1;
			for(vector1<Rose>::const_iterator i = roses.begin(); i != roses.end(); ++i,++ii){
				for(vector1<Link>::const_iterator j = links.begin(); j != links.end(); ++j){
					if( j->clashes(*i) ){
						++clash_roselink_[ii];
						return true; //clash = true;
					}
				}
				for(vector1<Rose>::const_iterator j = i+1; j != roses.end(); ++j){
					if( j->clashes(*i) ){
						++clash_roserose_[ii];
						return true; //clash = true;
					}
				}
			}
			return clash;
		}
		bool is_connected(){
			for(vector1<Link>::const_iterator i=links.begin(); i!=links.end(); ++i) if(!i->is_connected()) return false;
			return true;
		}
	};
int main(int argv, char **argc){
  try {

	register_options();
	devel::init(argv,argc);

	core::Real temperature = option[flail::temperature]();
	vector1<core::Size> hist(1000,0);

	core::scoring::RamachandranCOP rama = new core::scoring::Ramachandran;

	vector1<std::string> linker_seqs;
	if( option[flail::linker_seqs].user() ){
		linker_seqs = option[flail::linker_seqs]();
	} else if( option[flail::linker_lens].user() ){
		for(vector1<int>::const_iterator i = option[flail::linker_lens]().begin(); i != option[flail::linker_lens]().end(); ++i){
			linker_seqs.push_back(string("GGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGS").substr(0,*i));
		}
	} else {
		utility_exit_with_message("no linker params specified!!");
	}
	vector1<core::Size> target_residues = option[flail::target_residues]();
	vector1<PoseOP> poses = core::import_pose::poseOPs_from_files(option[in::file::s](),"centroid", core::import_pose::PDB_file);

	Bola bola(poses,linker_seqs,target_residues,rama);
	// bola.random_nonclashing_move();
	// assert(bola.is_connected());

	std::string tag = string_of(option[out::nstruct]())+"_";
	for(vector1<string>::const_iterator i = linker_seqs.begin(); i != linker_seqs.end(); ++i) tag += string_of(i->size())+"_";
	tag += string_of(temperature)+"_";
	tag += string_of(bola.move_sd_)+"_";
	tag += string_of(numeric::random::uniform()).substr(2,11);

	time_t t_start,t_now;
	time (&t_start);

	cout << "starting..." << endl;
	Real topscore = 9e9;
	for(Size sampling_round = 1; sampling_round < 99999999999; ++sampling_round){
		string sampling_round_str = F(8,6,((Real)sampling_round)/1000000).substr(2);
		core::Real last_score = bola.get_target_distance();
		for(core::Size i = 1; i <= (core::Size)option[out::nstruct](); ++i){
			core::chemical::AA moved_aa;
			Real moved_phi,moved_psi;
			bola.random_nonclashing_move(moved_aa,moved_phi,moved_psi);
			core::Real score = bola.get_target_distance();
			if( numeric::random::uniform() < fastexp((last_score-score)/temperature) ){
				++bola.naccepts;
				++hist[(int)score+1];
				last_score = score;
				if(score < 8.0) {
					cout << tag<<"_"<<sampling_round_str<< " GOODSCORE " << score << " " << option[out::file::o]()+"/good_score_"+tag+"_"+sampling_round_str+".pdb.gz" << endl;
					bola.dump_pdb(option[out::file::o]()+"/good_score_"+tag+"_"+sampling_round_str+".pdb.gz");
				}
				if(score<topscore){
					topscore = score;
					cout << tag<<"_"<<sampling_round_str<< " BESTSCORE " << score << " " << option[out::file::o]()+"/best_score_"+F(7,3,score)+"_"+tag+"_"+sampling_round_str+".pdb.gz" << endl;
					if(score<10.0) bola.dump_pdb(option[out::file::o]()+"/best_score_"+F(7,3,score)+"_"+tag+"_"+sampling_round_str+".pdb.gz");
				}
			} else {
				bola.undo_most_recent_move();
			}
		}

		// hist crap
		Size tot=0;
		for(vector1<core::Size>::const_iterator i = hist.begin(); i != hist.end(); ++i) tot += *i;
		vector1<core::Real> hist2(hist.size(),0.0);
		for(core::Size i = 1; i <= hist.size(); ++i) hist2[i] = (core::Real)hist[i]*exp(((core::Real)i-0.5)/temperature);
		core::Real sum=0.0;
		for(core::Size i = 1; i <= hist2.size(); ++i) sum += hist2[i];
		for(core::Size i = 1; i <= hist2.size(); ++i) hist2[i] = hist2[i]/sum;

		time(&t_now);
		Real timedif = difftime(t_now,t_start);

		// cout << tag<<"_"<<sampling_round_str<<" ";
		// cout << "clash_intralink_: "; for(Size i=1;i<=bola.clash_intralink_.size();++i) cout << bola.clash_intralink_[i] << " ";
		// cout << "clash_crosslink_: "; for(Size i=1;i<=bola.clash_crosslink_.size();++i) cout << bola.clash_crosslink_[i] << " ";
		// cout << "clash_roselink_:  "; for(Size i=1;i<=bola.clash_roselink_ .size();++i) cout << bola.clash_roselink_ [i] << " ";
		// cout << "clash_roserose_:  "; for(Size i=1;i<=bola.clash_roserose_ .size();++i) cout << bola.clash_roserose_ [i] << " ";
		// cout << endl;
		cout << tag<<"_"<<sampling_round_str<<" " << "trials:           " << bola.ntrials << endl;
		// cout << tag<<"_"<<sampling_round_str<<" " << "nclash trials:    " << bola.ntrials-bola.nclashes << endl;
		cout << tag<<"_"<<sampling_round_str<<" " << "naccepts:         " << bola.naccepts << endl;
		cout << tag<<"_"<<sampling_round_str<<" " << "clash/ntrial:     " << (core::Real)bola.nclashes/(core::Real)bola.ntrials << endl;
		cout << tag<<"_"<<sampling_round_str<<" " << "naccepts/ntrials: " << (core::Real)bola.naccepts/(core::Real)bola.ntrials << endl;
		cout << tag<<"_"<<sampling_round_str<<" " << "profiling: " << timedif << " " << ((Real)bola.naccepts)/timedif << " " << ((Real)bola.ntrials)/timedif << " " << bola.move_sd_ << endl;

		cout << F(7,3,temperature)<<" "<<F(7,3,bola.move_sd_)<<" "<< I(12,tot) <<" HIST ";
		for(int i=1;i<=20;++i) cout << " "<<I(3,hist[i]);
		for(int i=1;i<=20;++i) cout << " "<<hist2[i];
		cout << endl;

		utility::io::ozstream out(option[out::file::o]()+"/hist_"+tag+"_"+sampling_round_str+".dat");
		for(vector1<core::Size>::const_iterator i = hist.begin(); i != hist.end(); ++i) out << *i << " ";
		out << endl;
		out.close();
		cout << tag<<"_"<<sampling_round_str<< " SCORE " << bola.get_target_distance() << " " << option[out::file::o]()+"/lastaccept_"+tag+"_"+sampling_round_str+".pdb.gz" << endl;
		bola.dump_pdb(option[out::file::o]()+"/lastaccept_"+tag+"_"+sampling_round_str+".pdb.gz");

		if( option[flail::walltime_limit]() < timedif) break;
	}

  } catch ( utility::excn::EXCN_Base const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }

  return 0;

}
