// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/scoring/motif/xfrags.hh>

#include <core/scoring/motif/motif_hash_stuff.hh>

	#include <ObjexxFCL/FArray2D.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <basic/Tracer.hh>
	#include <basic/database/open.hh>
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>
	#include <basic/options/option_macros.hh>
	#include <basic/pymol_chains.hh>
	#include <core/chemical/AtomType.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/conformation/symmetry/util.hh>
	#include <core/conformation/Residue.hh>
	#include <core/conformation/ResidueFactory.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/io/silent/SilentFileData.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/annotated_sequence.hh>
	#include <core/pose/util.hh>
	#include <core/pose/symmetry/util.hh>
	#include <core/conformation/symmetry/SymmetryInfo.hh>
	#include <core/pose/xyzStripeHashPose.hh>
	#include <core/io/pdb/pose_io.hh>
	#include <core/scoring/Energies.hh>
	#include <core/scoring/EnergyGraph.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/ScoreTypeManager.hh>
	#include <core/scoring/dssp/Dssp.hh>
	#include <core/scoring/dssp/StrandPairing.hh>
	#include <core/scoring/hbonds/HBondOptions.hh>
	#include <core/scoring/methods/EnergyMethodOptions.hh>
	#include <core/scoring/packing/compute_holes_score.hh>
	#include <core/scoring/rms_util.hh>
	#include <core/scoring/sasa.hh>
	#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
	#include <numeric/conversions.hh>
	#include <numeric/model_quality/rms.hh>
	#include <numeric/random/random.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <numeric/xyzVector.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/file/file_sys_util.hh>
	#include <utility/fixedsizearray1.hh>
	#include <numeric/xyzTransform.hh>

	#include <numeric/geometry/hashing/SixDHasher.hh>

	#include <boost/unordered_set.hpp>

	#include <boost/foreach.hpp>
	#include <bitset>

    #ifndef _WIN32
	#include <pthread.h>
    #endif


	// #include <core/pack/task/PackerTask.hh>
	// #include <core/pack/packer_neighbors.hh>
	// #include <core/pack/rotamer_set/RotamerSet.hh>
	// #include <core/pack/rotamer_set/RotamerSetFactory.hh>
	// #include <core/pack/task/TaskFactory.hh>

	#define MAX_UINT16 65535
	#define MAX_UINT8    255
	#define XFORM_SCORE_FILE_VERSION 1

namespace core {
namespace scoring {
namespace motif {

	using numeric::Xforms;


/************************************************* types ************************************************/
	static basic::Tracer TR("core.scoring.motif");

	using core::pose::PoseCoordPickMode_BB;
	using core::pose::PoseCoordPickMode_N_C_O;

	typedef utility::fixedsizearray1<float,20> float20;
	typedef utility::fixedsizearray1<float,9> float9;
	using std::make_pair;
	using core::chemical::AA;
	using core::id::AtomID;
	using basic::options::option;
	namespace mh = basic::options::OptionKeys::mh;
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
	using namespace ObjexxFCL::format;
	using ObjexxFCL::string_of;
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::ostream;
	using std::string;
	using utility::io::izstream;
	using utility::io::ozstream;
	using utility::file_basename;
	using utility::vector1;
	using std::endl;
	using core::import_pose::pose_from_pdb;
	using numeric::geometry::hashing::Real3;
	using numeric::geometry::hashing::Real6;
	using core::pose::xyzStripeHashPoseCOP;

/************************************************ helpers ***********************************************/
	// static inline double sqr(double const & x) { return x*x; }
	// static inline float  sqr(float  const & x) { return x*x; }

	// static inline float
	// fastpow2 (float p)
	// {
	//   float offset = (p < 0) ? 1.0f : 0.0f;
	//   float clipp = (p < -126) ? -126.0f : p;
	//   int w = clipp;
	//   float z = clipp - w + offset;
	//   union { uint32_t i; float f; } v = { (1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z) };

	//   return v.f;
	//  }
	// static inline float
	// fastexp (float p)
	// {
	//   return fastpow2 (1.442695040f * p);
	//  }
	// static inline float
	// fastsinh (float p)
	// {
	//   return 0.5f * (fastexp (p) - fastexp (-p));
	//  }
	// static void xform_pose( core::pose::Pose & pose, Xform const & s, Size sres=1, Size eres=0 ) {
	//   if(eres==0) eres = pose.n_residue();
	//   for(Size ir = sres; ir <= eres; ++ir) {
	//     for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
	//       core::id::AtomID const aid(core::id::AtomID(ia,ir));
	//       pose.set_xyz( aid, s*pose.xyz(aid) );
	//     }
	//   }
	//  }

	static utility::fixedsizearray1<Real,6> get_bins(Real c, Real a){
		utility::fixedsizearray1<Real,6> bins(c);
		bins[4] = a;
		bins[5] = a;
		bins[6] = a;
		return bins;
	 }

 	static inline Real uint16_to_real(uint16_t const & val, Real const & lb, Real const & ub){
		return min(ub,max(lb,((Real)val)/((Real)MAX_UINT16) * (ub-lb) + lb));
	 }
	static inline uint16_t real_to_uint16(Real const & val, Real const & lb, Real const & ub){
		return (uint16_t)( min(1.0,max(0.0,(val-lb)/(ub-lb))) * (Real)MAX_UINT16 );
	 }

	static inline Real uint8_to_real(uint8_t const & val, Real const & lb, Real const & ub){
		return min(ub,max(lb,((Real)val)/((Real)MAX_UINT8) * (ub-lb) + lb));
	 }
	static inline uint8_t real_to_uint8(Real const & val, Real const & lb, Real const & ub){
		return (uint16_t)( min(1.0,max(0.0,(val-lb)/(ub-lb))) * (Real)MAX_UINT8 );
	 }


		Xfres::Xfres(char _aa, char _ss, Real _phi, Real _psi, Real _omg, Real _chi1, Real _chi2, Real _chi3, Real _chi4){
			aa  (_aa  );
			ss  (_ss  );
			phi (_phi );
			psi (_psi );
			omg (_omg );
			chi1(_chi1);
			chi2(_chi2);
			chi3(_chi3);
			chi4(_chi4);
		}
		Real Xfres::phi () const { return uint16_to_real(phi_  ,-180.0,180.0); }
		Real Xfres::psi () const { return uint16_to_real(psi_  ,-180.0,180.0); }
		Real Xfres::omg () const { return uint16_to_real(omg_  ,-180.0,180.0); }
		Real Xfres::chi1() const { return uint8_to_real(chi_[1],-180.0,180.0); }
		Real Xfres::chi2() const { return uint8_to_real(chi_[2],-180.0,180.0); }
		Real Xfres::chi3() const { return uint8_to_real(chi_[3],-180.0,180.0); }
		Real Xfres::chi4() const { return uint8_to_real(chi_[4],-180.0,180.0); }
		void Xfres::phi (Real _phi) { phi_    = real_to_uint16(_phi,-180.0,180.0); }
		void Xfres::psi (Real _psi) { psi_    = real_to_uint16(_psi,-180.0,180.0); }
		void Xfres::omg (Real _omg) { omg_    = real_to_uint16(_omg,-180.0,180.0); }
		void Xfres::chi1(Real _chi1){ chi_[1] = real_to_uint8(_chi1,-180.0,180.0); }
		void Xfres::chi2(Real _chi2){ chi_[2] = real_to_uint8(_chi2,-180.0,180.0); }
		void Xfres::chi3(Real _chi3){ chi_[3] = real_to_uint8(_chi3,-180.0,180.0); }
		void Xfres::chi4(Real _chi4){ chi_[4] = real_to_uint8(_chi4,-180.0,180.0); }
		void Xfres::place_sidechain_in_pose(Pose & pose, int ir) const{
			string const & name( core::chemical::name_from_aa(core::chemical::aa_from_oneletter_code(aa())) );
			core::conformation::ResidueOP newres = core::conformation::ResidueFactory::create_residue(pose.residue(ir).residue_type_set().name_map(name));
			pose.replace_residue(ir,*newres,true);
			if(pose.residue(ir).nchi() >= 1) pose.set_chi(1,ir,chi1());
			if(pose.residue(ir).nchi() >= 2) pose.set_chi(2,ir,chi2());
			if(pose.residue(ir).nchi() >= 3) pose.set_chi(3,ir,chi3());
			if(pose.residue(ir).nchi() >= 4) pose.set_chi(4,ir,chi4());
		}
		ostream & operator << (ostream & out, Xfres const & x){
			out<<"Xfres "<<x.aa()<<" "<<x.ss()<<" "<<F(8,3,x.phi())<<" "<<F(8,3,x.psi())<<" "<<F(8,3,x.omg())<<" "<<F(8,3,x.chi1())<<" "<<F(8,3,x.chi2())<<" "<<F(8,3,x.chi3())<<" "<<F(8,3,x.chi4());
			return out;
		}
		std::istream & operator >> (std::istream & in , Xfres & x){
			string tag; in >> tag; if("Xfres"!=tag) utility_exit_with_message("bad tag for Xfres");
			char _aa,_ss; Real _phi,_psi,_omg,_chi1,_chi2,_chi3,_chi4;
			in>>_aa>>_ss>>_phi>>_psi>>_omg>>_chi1>>_chi2>>_chi3>>_chi4;
			x.aa  (_aa  );
			x.ss  (_ss  );
			x.phi (_phi );
			x.psi (_psi );
			x.omg (_omg );
			x.chi1(_chi1);
			x.chi2(_chi2);
			x.chi3(_chi3);
			x.chi4(_chi4);
			return in;
		}
		bool read_xfres_binary( string  const & fname , vector1<Xfres> & xfres){
			utility::io::izstream in(fname,std::ios::in|std::ios::binary);
			if(!in.good()){
				TR.Error << "read_xfres_binary(fname): problem opening xfres input file " << fname << endl;
				return false;
			}
			if(!read_xfres_binary(in,xfres)){
				TR.Error << "read_xfres_binary(fname): problem while reading xfres file " << fname << endl;
				return false;
			}
			return true;
		}
		bool read_xfres_binary( std::istream & in, vector1<Xfres> & xfres){
			uint64_t n0=xfres.size(),n=0;
			in.read((char*)&n,sizeof(uint64_t));
			xfres.resize(xfres.size()+n);
			for(Size i = 1; i <= n; ++i){
				in.read((char*)&xfres[n0+i],sizeof(Xfres));
			}
			return true;
		}
		bool write_xfres_binary( ostream & out , vector1<Xfres> const & xfres){
			uint64_t n = xfres.size();
			out.write((char*)&n,sizeof(uint64_t));
			cout << "write_xfres_binary " << n << endl;
			for(vector1<Xfres>::const_iterator i = xfres.begin(); i != xfres.end(); ++i){
				Xfres const & sm( *i );
				out.write((char*)&sm,sizeof(Xfres));
				if(!out.good()) return false;
			}
			return true;
		}
		bool write_xfres_binary( string  const & fname , vector1<Xfres> const & xfres){
			utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
			if(!out.good()){
				TR.Error << "write_xfres_binary(fname): problem opening xfres output file " << fname << endl;
				return false;
			}
			if(!write_xfres_binary(out,xfres)){
				TR.Error << "write_xfres_binary(fname): problem while writing xfres file " << fname << endl;
				return false;
			}
			if(!out.good()){
				TR.Error << "write_xfres_binary(fname): problem after writing xfres file " << fname << endl;
				out.close();
				return false;
			}
			out.close();
			return true;
		}
		bool read_xfres_binary( vector1<string> const & fnames, vector1<Xfres> & xfres){
			Size count = 0;
			TR << "Nfiles: " << fnames.size() << endl;
			for(vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i){
				if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
				if(!read_xfres_binary(*i,xfres)){
					TR.Error << "read_xfres_binary(fnames): error reading file "+*i << endl;
					return false;
				}
			}
			TR << endl;
			return true;
		}
		Xfrag::Xfrag(Real6 _rt, Size _position, uint8_t _size, char _sscomp, Real _bfac, Real _burial){
			rt6(_rt);
			position(_position);
			sscomp(_sscomp);
			size(_size);
			bfac(_bfac);
			burial(_burial);
		}
		Real6 Xfrag::rt6() const {
			Real6 rt6;
			for(Size i = 1; i <= 3; ++i) rt6[i] = uint16_to_real(rt6_[i],
				-MOTIF_HASH_CART_SIZE,
				 MOTIF_HASH_CART_SIZE);
			for(Size i = 4; i <= 5; ++i) rt6[i] = uint16_to_real(rt6_[i], 0.0, 360.0 );
			for(Size i = 6; i <= 6; ++i) rt6[i] = uint16_to_real(rt6_[i], 0.0, 180.0 );
			return rt6;
		}
		void  Xfrag::rt6(Real6 const & rt_in) {
			Real6 rt6 = rt_in;
			rt6[4] = fmod(rt_in[4],360.0);
			rt6[5] = fmod(rt_in[5],360.0);
			rt6[6] = fmod(rt_in[6],360.0);
			rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
			rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
			rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
			for(Size i = 1; i <= 3; ++i) rt6_[i] = real_to_uint16(rt_in[i],
				-MOTIF_HASH_CART_SIZE,
				 MOTIF_HASH_CART_SIZE);
			for(Size i = 4; i <= 5; ++i) rt6_[i] = real_to_uint16(rt_in[i], 0.0, 360.0 );
			for(Size i = 6; i <= 6; ++i) rt6_[i] = real_to_uint16(rt_in[i], 0.0, 180.0 );
		}
		Size    Xfrag::position() const { return position_; }
		Size    Xfrag::size    () const { return size_; }
		char    Xfrag::sscomp  () const { return sscomp_; }
		Real    Xfrag::bfac    () const { return uint8_to_real(bfac_,  0.0, 2.55); }
		Real    Xfrag::burial  () const { return uint8_to_real(burial_,0.0, 1.0 ); }
		Real    Xfrag::ex1     () const { return uint8_to_real(ex1_,   0.0, 1.0 ); }
		Real    Xfrag::ex2     () const { return uint8_to_real(ex2_,   0.0, 1.0 ); }
		Real    Xfrag::ex3     () const { return uint8_to_real(ex3_,   0.0, 1.0 ); }
		Real    Xfrag::ex4     () const { return uint8_to_real(ex4_,   0.0, 1.0 ); }
		void    Xfrag::position(Size    _position) { position_ = _position; }
		void    Xfrag::size    (Size    _size    ) { size_     = _size; }
		void    Xfrag::sscomp  (char    _sscomp  ) { sscomp_   = _sscomp; }
		void    Xfrag::bfac    (Real    _bfac    ) { bfac_     = real_to_uint8(_bfac,  0.0, 2.55); }
		void    Xfrag::burial  (Real    _burial  ) { burial_   = real_to_uint8(_burial,0.0, 1.0 ); }
		void    Xfrag::ex1     (Real    _ex1     ) { ex1_      = real_to_uint8(_ex1,   0.0, 1.0 ); }
		void    Xfrag::ex2     (Real    _ex2     ) { ex2_      = real_to_uint8(_ex2,   0.0, 1.0 ); }
		void    Xfrag::ex3     (Real    _ex3     ) { ex3_      = real_to_uint8(_ex3,   0.0, 1.0 ); }
		void    Xfrag::ex4     (Real    _ex4     ) { ex4_      = real_to_uint8(_ex4,   0.0, 1.0 ); }
		void    Xfrag::insert (Pose & pose, vector1<Xfres> const & xfres, int lowres, int highres) const {
			int cutres=0;
			for(int ic=1; ic <= pose.fold_tree().num_cutpoint(); ++ic){
				int c = pose.fold_tree().cutpoint(ic);
				if(c < lowres || c >= highres) continue;
				if(cutres) utility_exit_with_message("more than one cutpoint in insertion");
				cutres = c;
			}
			if(!cutres) utility_exit_with_message("must be a cutpoint in inertion");
			core::kinematics::FoldTree ft( pose.fold_tree() );
			for(int ij=1; ij <= (int)pose.num_jump(); ++ij){
				int ur = pose.fold_tree().upstream_jump_residue(ij);
				int dr = pose.fold_tree().downstream_jump_residue(ij);
				if( lowres < ur && ur <= cutres  ) { ft.slide_jump(ij, lowres, dr);      cout << "ft.slide_jump("<<ij<<", "<<lowres   <<", "<<dr<<");" << endl; }
				if( lowres < dr && dr <= cutres  ) { ft.slide_jump(ij, ur, lowres);      cout << "ft.slide_jump("<<ij<<", "<<ur       <<", "<<lowres<<");" << endl; }
				if( cutres < ur && ur <= highres ) { ft.slide_jump(ij, highres+1, dr);   cout << "ft.slide_jump("<<ij<<", "<<highres+1<<", "<<dr<<");" << endl; }
				if( cutres < dr && dr <= highres ) { ft.slide_jump(ij, ur, highres+1);   cout << "ft.slide_jump("<<ij<<", "<<ur       <<", "<<highres+1<<");" << endl; }
				cout << ft;
				for(int i=1;i<=(int)ft.nres();++i) cout << i << " " << ft.is_jump_point(i) << endl;
			}
			ft.slide_cutpoint(cutres, highres);
			ft.reorder(1);
			cout << ft;
			for(int i=1;i<=(int)ft.nres();++i) cout << i << " " << ft.is_jump_point(i) << endl;
			cout << "OLD " << pose.fold_tree();
			pose.fold_tree(ft);
			cout << "NEW " << pose.fold_tree();

			for(int ir=lowres+1; ir <= highres; ++ir){
				cout << pose.fold_tree();
				cout << "delete " << lowres+1 << endl;
				pose.delete_polymer_residue(lowres+1);
			}
			core::conformation::ResidueOP dummyres = core::conformation::ResidueFactory::create_residue(pose.residue(1).residue_type_set().name_map("GLY"));
			for(int ir=2; ir <= (int)size(); ++ir) pose.append_polymer_residue_after_seqpos(*dummyres,lowres-2+ir,true);
			vector1<Xfres>::const_iterator ifrag = xfres.begin() + position();
			for(int ir=1; ir <= (int)size(); ++ir,++ifrag){
				int resno = lowres+ir-1;
				pose.set_phi  (resno,ifrag->phi());
				pose.set_psi  (resno,ifrag->psi());
				pose.set_omega(resno,ifrag->omg());
				ifrag->place_sidechain_in_pose(pose,resno);
			}
		}

		ostream & operator << (ostream & out, Xfrag const & x){
			out<<"Xfrag "<<I(9,x.position())<<" "<<I(2,x.size())<<" "<<x.sscomp()<<" "<<F(7,3,x.bfac())<<" "<<F(7,3,x.burial())<<" "<<F(7,3,x.ex1())<<" "<<F(7,3,x.ex2())<<" "<<F(7,3,x.ex3())<<" "<<F(7,3,x.ex4());
			out<<F(7,3,x.rt6()[1])<<" "<<F(7,3,x.rt6()[2])<<" "<<F(7,3,x.rt6()[3])<<" "<<F(7,3,x.rt6()[4])<<" "<<F(7,3,x.rt6()[5])<<" "<<F(7,3,x.rt6()[6]);
			return out;
		}
		std::istream & operator >> (std::istream & in , Xfrag & x){
			string tag; in >> tag; if("Xfrag"!=tag) utility_exit_with_message("bad tag for Xfrag");
			char _position,_size,_sscomp; Real _bfac,_burial,_ex1,_ex2,_ex3,_ex4;
			in>>_position>>_size>>_sscomp>>_bfac>>_burial>>_ex1>>_ex2>>_ex3>>_ex4;
			x.position  (_position  );
			x.size  (_size  );
			x.sscomp (_sscomp );
			x.bfac (_bfac );
			x.burial (_burial );
			x.ex1(_ex1);
			x.ex2(_ex2);
			x.ex3(_ex3);
			x.ex4(_ex4);
			return in;
		}
		bool read_xfrag_binary( string  const & fname , vector1<Xfrag> & xfrag){
			utility::io::izstream in(fname,std::ios::in|std::ios::binary);
			if(!in.good()){
				TR.Error << "read_xfrag_binary(fname): problem opening xfrag input file " << fname << endl;
				return false;
			}
			if(!read_xfrag_binary(in,xfrag)){
				TR.Error << "read_xfrag_binary(fname): problem while reading xfrag file " << fname << endl;
				return false;
			}
			return true;
		}
		bool read_xfrag_binary( std::istream & in, vector1<Xfrag> & xfrag){
			uint64_t n0=xfrag.size(),n=0;
			in.read((char*)&n,sizeof(uint64_t));
			xfrag.resize(xfrag.size()+n);
			for(Size i = 1; i <= n; ++i){
				in.read((char*)&xfrag[n0+i],sizeof(Xfrag));
			}
			return true;
		}
		bool write_xfrag_binary( ostream & out , vector1<Xfrag> const & xfrag){
			uint64_t n = xfrag.size();
			out.write((char*)&n,sizeof(uint64_t));
			cout << "write_xfrag_binary " << n << endl;
			for(vector1<Xfrag>::const_iterator i = xfrag.begin(); i != xfrag.end(); ++i){
				Xfrag const & sm( *i );
				out.write((char*)&sm,sizeof(Xfrag));
				if(!out.good()) return false;
			}
			return true;
		}
		bool write_xfrag_binary( string  const & fname , vector1<Xfrag> const & xfrag){
			utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
			if(!out.good()){
				TR.Error << "write_xfrag_binary(fname): problem opening xfrag output file " << fname << endl;
				return false;
			}
			if(!write_xfrag_binary(out,xfrag)){
				TR.Error << "write_xfrag_binary(fname): problem while writing xfrag file " << fname << endl;
				return false;
			}
			if(!out.good()){
				TR.Error << "write_xfrag_binary(fname): problem after writing xfrag file " << fname << endl;
				out.close();
				return false;
			}
			out.close();
			return true;
		}
		bool read_xfrag_binary( vector1<string> const & fnames, vector1<Xfrag> & xfrag){
			Size count = 0;
			TR << "Nfiles: " << fnames.size() << endl;
			for(vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i){
				if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
				if(!read_xfrag_binary(*i,xfrag)){
					TR.Error << "read_xfrag_binary(fnames): error reading file "+*i << endl;
					return false;
				}
			}
			TR << endl;
			return true;
		}




		bool read_xfrag_binary (        string  const & fname , vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres){
			utility::io::izstream in(fname,std::ios::in|std::ios::binary);
			if(!in.good()){
				TR.Error << "read_xfrag_binary(fname): problem opening xfrag input file " << fname << endl;
				return false;
			}
			if(!read_xfrag_binary(in,xfrag,xfres)){
				TR.Error << "read_xfrag_binary(fname): problem while reading xfrag file " << fname << endl;
				return false;
			}
			return true;

		}
		bool read_xfrag_binary (  std::istream        & in    , vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres){
			if(!read_xfres_binary(in,xfres)) return false;
			if(!read_xfrag_binary(in,xfrag)) return false;
			return true;
		}
		bool write_xfrag_binary(  ostream        & out   , vector1<Xfrag> const & xfrag, vector1<Xfres> const & xfres){
			if(!write_xfres_binary(out,xfres)) return false;
			if(!write_xfrag_binary(out,xfrag)) return false;
			return true;
		}
		bool write_xfrag_binary(        string  const & fname , vector1<Xfrag> const & xfrag, vector1<Xfres> const & xfres){
			utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
			if(!out.good()){
				TR.Error << "write_xfrag_binary(fname): problem opening xfrag output file " << fname << endl;
				return false;
			}
			if(!write_xfrag_binary(out,xfrag,xfres)){
				TR.Error << "write_xfrag_binary(fname): problem while writing xfrag file " << fname << endl;
				return false;
			}
			if(!out.good()){
				TR.Error << "write_xfrag_binary(fname): problem after writing xfrag file " << fname << endl;
				out.close();
				return false;
			}
			out.close();
			return true;
		}
		bool read_xfrag_binary (vector1<string> const & fnames, vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres){
			Size count = 0;
			TR << "Nfiles: " << fnames.size() << endl;
			for(vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i){
				if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
				if(!read_xfrag_binary(*i,xfrag,xfres)){
					TR.Error << "read_xfrag_binary(fnames): error reading file "+*i << endl;
					return false;
				}
			}
			TR << endl;
			return true;

		}


		XfragSet::XfragSet(/*core::Real cartsize,*/ core::Real cartresl, core::Real angleresl)
		 : /*cart_size_(cartsize),*/ cart_resl_(cartresl), angle_resl_(angleresl),
				hasher_(
				numeric::geometry::BoundingBox<Vector>(
					Vector(-MOTIF_HASH_CART_SIZE),
					Vector( MOTIF_HASH_CART_SIZE)
				),
				utility::fixedsizearray1<Size,3>(0),
				get_bins(cart_resl_,angle_resl_)
			)
	 {}

}
}
}

