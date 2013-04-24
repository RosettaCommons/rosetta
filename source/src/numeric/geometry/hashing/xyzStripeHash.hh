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
#include <ObjexxFCL/format.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <boost/foreach.hpp>

namespace numeric {
namespace geometry {
namespace hashing {

struct Ball {
	float x_,y_,z_;
	uint8_t radius_,atomno_;
	uint16_t resid_;
	inline float & x() { return x_; }
	inline float & y() { return y_; }
	inline float & z() { return z_; }
	inline float const & x() const { return x_; }
	inline float const & y() const { return y_; }
	inline float const & z() const { return z_; }
	inline float radius() const { return ((float)radius_)/20.0; }
	inline void  radius(float radius) { radius_ = (uint8_t)(20.0*radius); }
};

struct Counter {
	int count;
	Counter():count(0){}
	void visit( numeric::xyzVector<float> const &, numeric::xyzVector<float> const & ){ ++count; }
};

class xyzStripeHash : public utility::pointer::ReferenceCount {
	inline short  short_min( short const a,  short const b) { return (a < b) ? a : b; }
	inline short  short_max( short const a,  short const b) { return (a > b) ? a : b; }
	inline short ushort_min(unsigned short const a, unsigned short const b) { return (a < b) ? a : b; }
	inline short ushort_max(unsigned short const a, unsigned short const b) { return (a > b) ? a : b; }
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
	void
	init( utility::vector1<Ball> const & balls );
	virtual ~xyzStripeHash() {
		if(grid_atoms_)  delete grid_atoms_;
		if(grid_stripe_) delete grid_stripe_;
	}

	inline const_iterator begin() const { return const_iterator(grid_atoms_       ) ; }
	inline const_iterator end()   const { return const_iterator(grid_atoms_+natom_) ; }

	bool sanity_check() const;

	int
	nbcount( Vec const & v_in ) const;
	int
	nbcount_raw( Vec const & v ) const;
	bool
	clash( Vec const & v_in ) const;
	bool
	clash_not_resid( Vec const & v_in, int resid ) const;
	bool
	clash_raw( Vec const & v ) const;

	void fill_pairs(
		Vec const & v,
		int const & ir,
		utility::vector1<std::pair<int,int> > & pairs,
		float maxd2=0.0
	) const;

	template<typename V>
	void
	visit( Vec const & v_in, V & visitor ) const {
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
					Vec const & c = *((Vec*)(grid_atoms_+i));
					float const d2 = (x-c.x())*(x-c.x()) + (y-c.y())*(y-c.y()) + (z-c.z())*(z-c.z());
					if( d2 <= grid_size2_ ) {
						visitor.visit(v,c,d2);
					}
				}
			}
		}
	}

	template<typename V>
	void
	visit_lax( Vec const & v_in, float const vr, V & visitor ) const {
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
					Vec const & c = *((Vec*)(grid_atoms_+i));
					float cr = grid_atoms_[i].radius();
					visitor.visit(v,c,vr,cr);
				}
			}
		}
	}


	inline Ball  const * grid_atoms() const { return grid_atoms_; }
	inline ushort2 const * grid_stripe() const { return grid_stripe_; }
	inline int natom() const { return natom_; }
	inline int xdim () const { return  xdim_; }
	inline int ydim () const { return  ydim_; }
	inline int zdim () const { return  zdim_; }
	inline float  grid_size() const { return  grid_size_; }
	inline float  grid_size2() const { return  grid_size2_; }
	inline numeric::xyzVector<float> const & translation() const { return translation_; }

private:
	float grid_size_,grid_size2_;
	int natom_;
	Ball  const * grid_atoms_;
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
