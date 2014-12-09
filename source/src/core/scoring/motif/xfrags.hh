// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// includes
#ifndef INCLUDED_core_scoring_motif_xfrags_hh
#define INCLUDED_core_scoring_motif_xfrags_hh

// #include <core/scoring/motif/xfrags.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/xyzStripeHashPose.fwd.hh>
#include <numeric/xyzVector.hh>
#include <utility/fixedsizearray1.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/HomogeneousTransform.hh>
#include <boost/unordered_map.hpp>

namespace core {
namespace scoring {
namespace motif {

//types

	using core::id::AtomID;
	using core::pose::Pose;
	using core::pose::PoseCOP;
	using core::pose::PoseCAP;
	using core::Real;
	using core::scoring::ScoreFunctionOP;
	using core::Size;
	using std::string;
	using utility::vector1;
	using numeric::geometry::hashing::Real3;
	using numeric::geometry::hashing::Real6;
	using core::pose::xyzStripeHashPose;
	using core::pose::xyzStripeHashPoseCAP;

	typedef utility::vector1<Real> Reals;
	typedef utility::vector1<Size> Sizes;
	typedef utility::vector1<int>  Ints;
	typedef utility::vector1<float> Floats;
	typedef utility::vector1<bool> Bools;
	typedef numeric::xyzVector<core::Real> Vec;
	typedef numeric::xyzMatrix<core::Real> Mat;
	typedef numeric::xyzTransform<core::Real> Xform;

	class Xfres {
		/*  2 */ char aa_,ss_;
		/*  6 */ uint16_t phi_,psi_,omg_;
		/*  4 */ utility::fixedsizearray1<uint8_t,4> chi_;
		public:
		Xfres() {}
		Xfres(char aa, char ss, Real phi, Real psi, Real omg, Real chi1=-180.0, Real chi2=-180.0, Real chi3=-180.0, Real chi4=-180.0);
		char aa  () const { return aa_; }
		char ss  () const { return ss_; }
		Real phi () const;
		Real psi () const;
		Real omg () const;
		Real chi1() const;
		Real chi2() const;
		Real chi3() const;
		Real chi4() const;
		void aa  (char aa  ) { aa_ = aa; }
		void ss  (char ss  ) { ss_ = ss; }
		void phi (Real _phi );
		void psi (Real _psi );
		void omg (Real _omg );
		void chi1(Real chi1);
		void chi2(Real chi2);
		void chi3(Real chi3);
		void chi4(Real chi4);
		void place_sidechain_in_pose(Pose & pose, int ir) const;
		friend std::ostream & operator << (std::ostream & out, Xfres const & x);
	};
	std::ostream & operator << (std::ostream & out, Xfres const & x);
	std::istream & operator >> (std::istream & in , Xfres       & x);
	bool read_xfres_binary (        string const & fname , vector1<Xfres>       & Xfres);
	bool read_xfres_binary (  std::istream        & in    , vector1<Xfres>       & Xfres);
	bool write_xfres_binary(  std::ostream        & out   , vector1<Xfres> const & Xfres);
	bool write_xfres_binary(        string const & fname , vector1<Xfres> const & Xfres);
	bool read_xfres_binary (vector1<string> const & fnames, vector1<Xfres>       & Xfres);

	class Xfrag {
		/* 12 */ utility::fixedsizearray1<uint16_t,6> rt6_;
		/*  4 */ uint32_t position_;
		char sscomp_;
		/*  4 */ uint8_t  size_,bfac_,burial_;
		/*  4 */ uint8_t  ex1_,ex2_,ex3_,ex4_;
		public:
		Xfrag(){}
		Xfrag(Real6 _rt, Size _position, uint8_t _size, char _sscomp, Real _bfac, Real _burial);
		Real6   rt6     () const;
		Size    position() const;
		Size    size    () const;
		char    sscomp  () const;
		Real    bfac    () const;
		Real    burial  () const;
		Real    ex1     () const;
		Real    ex2     () const;
		Real    ex3     () const;
		Real    ex4     () const;
		void    rt6     (Real6 const & _rt6);
		void    position(Size    _position);
		void    size    (Size    _size    );
		void    sscomp  (char    _sscomp  );
		void    bfac    (Real    _bfac    );
		void    burial  (Real    _burial  );
		void    ex1     (Real    _ex1     );
		void    ex2     (Real    _ex2     );
		void    ex3     (Real    _ex3     );
		void    ex4     (Real    _ex4     );
		void insert(Pose & pose,  vector1<Xfres> const & xfres, int lowres, int highres) const;
	};
	std::ostream & operator << (std::ostream & out, Xfrag const & x);
	std::istream & operator >> (std::istream & in , Xfrag       & x);
	bool read_xfrag_binary (        string const & fname , vector1<Xfrag>       & xfrag);
	bool read_xfrag_binary (  std::istream        & in    , vector1<Xfrag>       & xfrag);
	bool write_xfrag_binary(  std::ostream        & out   , vector1<Xfrag> const & xfrag);
	bool write_xfrag_binary(        string const & fname , vector1<Xfrag> const & xfrag);
	bool read_xfrag_binary (vector1<string> const & fnames, vector1<Xfrag>       & xfrag);


	bool read_xfrag_binary (        string const & fname , vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres);
	bool read_xfrag_binary (  std::istream        & in    , vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres);
	bool write_xfrag_binary(  std::ostream        & out   , vector1<Xfrag> const & xfrag, vector1<Xfres> const & xfres);
	bool write_xfrag_binary(        string const & fname , vector1<Xfrag> const & xfrag, vector1<Xfres> const & xfres);
	bool read_xfrag_binary (vector1<string> const & fnames, vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres);

	class XfragSet : public utility::pointer::ReferenceCount {
		typedef boost::uint64_t Key;
		typedef numeric::geometry::hashing::SixDCoordinateBinner SixDCoordinateBinner;
		typedef numeric::geometry::hashing::bin_index_hasher bin_index_hasher;
		typedef boost::unordered_multimap< Key, Xfrag, bin_index_hasher > XfragMap;
		double const /*cart_size_,*/cart_resl_,angle_resl_;
		SixDCoordinateBinner hasher_;
		XfragMap xfragmap_;
		utility::vector1<Xfres> xfres_;
		public:
		XfragSet(/*core::Real cartsize=10.0,*/ core::Real cartresl=1.0, core::Real angleresl=15.0);
		// Undefined, commenting out to fix PyRosetta build  void write_binary(std::string fname) const;
		// Undefined, commenting out to fix PyRosetta build  void read_binary(std::string fname);
		// Undefined, commenting out to fix PyRosetta build  void add_frags(utility::vector1<Xfrag> const & frags, utility::vector1<Xfres> const & res);
		// Undefined, commenting out to fix PyRosetta build  void lookup_frags(Real6 const & rt6, utility::vector1<Xfrag> & frags) const;
	};






}
}
}

#endif


