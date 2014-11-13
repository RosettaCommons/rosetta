// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_numeric_geometry_hashing_xyzStripeHash_hh
#define INCLUDED_numeric_geometry_hashing_xyzStripeHash_hh

#include <numeric/geometry/hashing/xyzStripeHash.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzTransform.hh>
#include <ObjexxFCL/format.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <boost/foreach.hpp>

#include <map>

namespace numeric {
namespace geometry {
namespace hashing {

struct Ball {
	float x_,y_,z_;
	uint8_t radius_,atomno_;
	uint16_t resid_;
	float & x() { return x_; }
	float & y() { return y_; }
	float & z() { return z_; }
	float const & x() const { return x_; }
	float const & y() const { return y_; }
	float const & z() const { return z_; }
	xyzVector<Real> xyz_real() const { return xyzVector<Real>(x_,y_,z_); }
	xyzVector<float> const & xyz() const { return *((xyzVector_float const*)this); }
	float radius() const { return ((float)radius_)/20.0; }
	void  radius(float radius) { radius_ = (uint8_t)(20.0*radius); }
	Size  resi() const { return resid_; }
	Size atomno() const { return atomno_; }
};

struct Counter {
	int count;
	Counter():count(0){}
	void visit( numeric::xyzVector<float> const &, numeric::xyzVector<float> const & ){ ++count; }
 };

class xyzStripeHash : public utility::pointer::ReferenceCount {
	short  short_min( short const a,  short const b) { return (a < b) ? a : b; }
	short  short_max( short const a,  short const b) { return (a > b) ? a : b; }
	short ushort_min(unsigned short const a, unsigned short const b) { return (a < b) ? a : b; }
	short ushort_max(unsigned short const a, unsigned short const b) { return (a > b) ? a : b; }
 public:
	typedef unsigned short ushort;
	typedef struct { unsigned short x,y; } ushort2;
	typedef numeric::xyzVector<float> Vec;

	// iterators:
	template<class C>
	struct iter_base : public std::iterator<std::input_iterator_tag,float> {
		iter_base(Ball const *p) : p_(p) {}
		C & operator=(C const & r) { p_ = r.p_; return *this; }
		C & operator++() { ++p_; return static_cast<C &>(*this); }
		bool operator!=(C const & r) const { return (p_ != r.p_); }
		bool operator==(C const & r) const { return (p_ == r.p_); }
	protected:
		Ball const *p_;
	};
	struct const_iterator : public iter_base<const_iterator> {
		const_iterator(Ball const *p) : iter_base<const_iterator>(p) {}
		Vec const & operator *() { return *((Vec const *)(this->p_)); }
		Vec const * operator->() { return  ((Vec const *)(this->p_)); }
		float radius() { return this->p_->radius(); }
	};

 public:
	xyzStripeHash(
		float grid_size = 0.0,
	    utility::vector1<Ball> const & balls=utility::vector1<Ball>()
	);

	void init( utility::vector1<Ball> const & balls );

	virtual ~xyzStripeHash() {
		if(grid_balls_)  delete grid_balls_;
		if(grid_stripe_) delete grid_stripe_;
	 }

	const_iterator begin() const { return const_iterator(grid_balls_       ) ; }
	const_iterator end()   const { return const_iterator(grid_balls_+nballs_) ; }

	bool sanity_check() const;

	std::string debug_pdb(Xform const & x=Xform::identity()) const;

	int nbcount( Vec const & v_in ) const;
	int nbcount_raw( Vec const & v ) const;
	bool clash( Vec const & v_in ) const;
	bool clash_not_resid( Vec const & v_in, int const & resid, int const & resid2 = 0) const;
	bool clash_raw( Vec const & v ) const;
	float clash_amount( Vec const & v_in ) const;

	// @brief Check if dist(b_in, hb) < (b_in.radius + hb.radius) for any hb in hash.
	//
	// Input ball provided in global coordinate frame, clashes are evaluated in global frame.
	//
	// Returns 1-based ball index if clash is found, 0 otherwise.
	int clash_check_ball( Ball const & b ) const;

	// @brief Generate residue mapping (r_t, r) where:
	// 	any(dist(b_t, b) < (b_t.radius + b.radius)) for 
	// 		{b_t in target_balls | b_t.resi = r_t}, {b in this | b.resi == r}
	//
	// Populated residue_pairs with the first identified clash per r_t, meaning that
	// all clashing resi in test_balls are included in mapping, however not all clashing
	// resi in this are guarenteed to be included.
	//
	// If residue_pairs already contains an entry for resi r_t the residue will not
	// be checked, and the mapping will not be updated.
	//
	// @returns true if any clashing pair was identified, false otherwise.
	bool clash_check_residue_pairs(
			utility::vector1<Ball> const & test_balls,
			std::map<Size, Size> & residue_pairs
	) const;

	void fill_pairs(
		xyzVector_float const & v,
		int const & ir,
		utility::vector1<std::pair<int,int> > & pairs,
		float maxd2=0.0
	 ) const;

	template<typename Visitor> void	visit( Vec const & v_in, Visitor & visitor ) const {
		Vec const v = v_in+translation_;
		float x = v.x(); float y = v.y(); float z = v.z();
		if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return; // worth it iff
		if( x > xmx_        || y > ymx_        || z > zmx_        ) return;                      // worth it iff
		int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,(int)(x/grid_size_));
		int const iy0  = (y<0) ? 0 : y/grid_size_;
		int const iz0  = (z<0) ? 0 : z/grid_size_;
		int const iyl = numeric::max(0,iy0-1);
		int const izl = numeric::max(0,iz0-1);
		int const iyu = numeric::min((int)ydim_,iy0+2);
		int const izu = numeric::min((int)zdim_,(int)iz0+2);
		for(int iy = iyl; iy < iyu; ++iy) {
			for(int iz = izl; iz < izu; ++iz) {
				int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
				assert(ig < xdim_*ydim_*zdim_);
				assert(ix < xdim_);
				assert(iy < ydim_);
				assert(iz < zdim_);
				int const & igl = grid_stripe_[ig].x;
				int const & igu = grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					Vec const & c = *((Vec*)(grid_balls_+i));
					float const d2 = (x-c.x())*(x-c.x()) + (y-c.y())*(y-c.y()) + (z-c.z())*(z-c.z());
					if( d2 <= grid_size2_ ) {
						visitor.visit(v,c,d2);
					}
				}
			}
		}
	 }
	template<typename Visitor> void	visit_lax( Vec const & v_in, float const vr, Visitor & visitor ) const {
		Vec const v = v_in+translation_;
		float x = v.x(); float y = v.y(); float z = v.z();
		if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return; // worth it iff
		if( x > xmx_        || y > ymx_        || z > zmx_        ) return;                      // worth it iff
		int const ix   = (x<0) ? 0 : numeric::min(xdim_-1, static_cast<int>(x/grid_size_));
		int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
		int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
		int const iyl = numeric::max(0,iy0-1);
		int const izl = numeric::max(0,iz0-1);
		int const iyu = numeric::min(static_cast<int>(ydim_), iy0+2);
		int const izu = numeric::min(static_cast<int>(zdim_), static_cast<int>(iz0+2));
		for(int iy = iyl; iy < iyu; ++iy) {
			for(int iz = izl; iz < izu; ++iz) {
				int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
				assert(ig < xdim_*ydim_*zdim_);
				assert(ix < xdim_);
				assert(iy < ydim_);
				assert(iz < zdim_);
				int const & igl = grid_stripe_[ig].x;
				int const & igu = grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					Vec const & c = *((Vec*)(grid_balls_+i));
					float cr = grid_balls_[i].radius();
					visitor.visit(v,c,vr,cr);
				}
			}
		}
	 }

	Ball const * grid_atoms() const { return grid_balls_; }
	Size size() const { return nballs_; }
	int natom() const { return nballs_; }
	int xdim () const { return  xdim_; }
	int ydim () const { return  ydim_; }
	int zdim () const { return  zdim_; }
	float  grid_size() const { return  grid_size_; }
	float  grid_size2() const { return  grid_size2_; }
	xyzVector_float const & translation() const { return translation_; }
	xyzVector<Real> translation_real() const { return xyzVector<Real>(translation_.x(),translation_.y(),translation_.z()); }

	ushort2 const * grid_stripe() const { return grid_stripe_; }

	Ball const &   ball(Size const & ib) const { assert(1ul<=ib&&ib<=(Size)nballs_); return grid_balls_[ib-1]; }
	xyzVector_float xyz(Size const & ib) const { assert(1ul<=ib&&ib<=(Size)nballs_); return grid_balls_[ib-1].xyz()-translation_; }
	Size           resi(Size const & ib) const { assert(1ul<=ib&&ib<=(Size)nballs_); return grid_balls_[ib-1].resi(); }

 private:
	float grid_size_,grid_size2_;
	int nballs_;
	Ball  const * grid_balls_;
	ushort2 const * grid_stripe_;
	int xdim_,ydim_,zdim_;
	float xmx_,ymx_,zmx_;
	//numeric::xyzMatrix<Real> rotation_;
	numeric::xyzVector<float> translation_;
	// neighbor_iterator neighbor_end_;
 };



} // namespace hashing
} // namespace geometry
} // namespace numeric

#endif
