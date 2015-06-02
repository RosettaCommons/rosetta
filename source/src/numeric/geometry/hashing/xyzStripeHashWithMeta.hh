// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_numeric_geometry_hashing_xyzStripeHashWithMeta_hh
#define INCLUDED_numeric_geometry_hashing_xyzStripeHashWithMeta_hh

#include <numeric/geometry/hashing/xyzStripeHashWithMeta.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/format.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace numeric {
namespace geometry {
namespace hashing {

template<typename T>
class xyzStripeHashWithMeta : public utility::pointer::ReferenceCount {

	inline short  short_min( short const a,  short const b) { return (a < b) ? a : b; }
	inline short  short_max( short const a,  short const b) { return (a > b) ? a : b; }
	inline short ushort_min(unsigned short const a, unsigned short const b) { return (a < b) ? a : b; }
	inline short ushort_max(unsigned short const a, unsigned short const b) { return (a > b) ? a : b; }

public:
	typedef std::pair<numeric::xyzVector<T>,T> VecandVal;
	typedef struct { T x,y,z,w; } float4;
	//typedef unsigned int uint;
	typedef unsigned short ushort;
	typedef struct { unsigned short x,y; } ushort2;
	typedef numeric::xyzVector<T> Vec;

	// iterators:
	template<class C>
	struct iter_base : public std::iterator<std::input_iterator_tag,T> {
		iter_base(float4 const *p) : p_(p) {}
		C & operator=(C const & r) { p_ = r.p_; return *this; }
		C & operator++() { ++p_; return static_cast<C &>(*this); }
		bool operator!=(C const & r) const { return (p_ != r.p_); }
		bool operator==(C const & r) const { return (p_ == r.p_); }
	protected:
		float4 const *p_;
	};
	struct xyzmeta_iterator : public iter_base<xyzmeta_iterator> {
		xyzmeta_iterator(float4 const *p) : iter_base<xyzmeta_iterator>(p) {}
		const VecandVal & operator*() { return *((VecandVal const *)(this->p_)); }
	};
	struct xyz_iterator : public iter_base<xyz_iterator> {
		xyz_iterator(float4 const *p) : iter_base<xyz_iterator>(p) {}
		const xyzVector<T> & operator*() { return ((VecandVal const *)(this->p_))->first; }
	};
	struct meta_iterator : public iter_base<meta_iterator> {
		meta_iterator(float4 const *p) : iter_base<meta_iterator>(p) {}
		const T & operator*() { return *((VecandVal const *)(this->p_)); }
	};

public:

	xyzStripeHashWithMeta(
		T grid_size
	):
		grid_size_(grid_size),
		grid_size2_(grid_size*grid_size),
		grid_atoms_(NULL),
		grid_stripe_(NULL)//,
		//	neighbor_end_(*this)
	{}
	xyzStripeHashWithMeta( T grid_size,
	               utility::vector1<numeric::xyzVector<T> > const & atoms,
	               utility::vector1<T> const & meta
	):
		grid_size_(grid_size),
		grid_size2_(grid_size*grid_size),
		grid_atoms_(NULL),
		grid_stripe_(NULL)//,
		//neighbor_end_(*this)
	{
		init(atoms,meta);
	}

	void init(
		utility::vector1<numeric::xyzVector<T> > const & atoms,
		utility::vector1<T> const & meta
	){
		// if( sizeof(T) < sizeof(M) ) utility_exit_with_message("octree metadata must fit in sizeof(T)!");
		if( meta.size() > 0 && atoms.size() != meta.size() ) utility_exit_with_message("must be metadata for each point!");
		if( atoms.size() > 65535 ) utility_exit_with_message("xyzStripeHashWithMeta con only handle < 65535 atoms!");

// #define FUDGE 0.0f

		natom_ = atoms.size();
		//neighbor_end_.end();

		T xmn= 9e9,ymn= 9e9,zmn= 9e9;
		T xmx=-9e9,ymx=-9e9,zmx=-9e9;
		for(int i = 1; i <= natom_; ++i) {
			xmn = numeric::min(xmn,atoms[i].x());
			ymn = numeric::min(ymn,atoms[i].y());
			zmn = numeric::min(zmn,atoms[i].z());
			xmx = numeric::max(xmx,atoms[i].x());
			ymx = numeric::max(ymx,atoms[i].y());
			zmx = numeric::max(zmx,atoms[i].z());
		}
		//TR<<xmx-xmn<<" "<<ymx-ymn<<" "<<zmx-zmn<<std::endl;
		// for(int i = 0; i < natom_; ++i) {
		//   atoms[i].x -= xmn/*-FUDGE*/;
		//   atoms[i].y -= ymn/*-FUDGE*/;
		//   atoms[i].z -= zmn/*-FUDGE*/;
		// }

		//std::cout << "xyzStripeHashWithMeta: " << xmn << " "  << ymn << " "  << zmn << " " << xmx << " "  << ymx << " "  << zmx << std::endl;

		xdim_ = (int)((xmx-xmn+0.0001)/grid_size_+0.999999);
		ydim_ = (int)((ymx-ymn+0.0001)/grid_size_+0.999999);
		zdim_ = (int)((zmx-zmn+0.0001)/grid_size_+0.999999);
		assert(xdim_ < 9999); assert(ydim_ < 9999); assert(zdim_ < 9999);
		int const gsize = xdim_*ydim_*zdim_;
		ushort2 *gindex  = new ushort2[gsize];
		ushort2 *gstripe = new ushort2[gsize];
		for(int i = 0; i < gsize; ++i) { gindex[i].y = 0; gindex[i].x = 0; }
		//TR<<"atom "<<natom_<<" grid1 "<<xdim_*ydim_*zdim_<<" "<<xdim_<<" "<<ydim_<<" "<<zdim_<<std::endl;

		for(int i = 1; i <= natom_; ++i) {
			int ix = (int)((atoms[i].x()-xmn/*+FUDGE*/)/grid_size_);
			int iy = (int)((atoms[i].y()-ymn/*+FUDGE*/)/grid_size_);
			int iz = (int)((atoms[i].z()-zmn/*+FUDGE*/)/grid_size_);
			assert(ix >= 0); assert(iy >= 0); assert(iz >= 0); assert(ix < xdim_); assert(iy < ydim_); assert(iz < zdim_);
			int ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig>=0);assert(ig<9999999);
			++(gindex[ig].y);
		}
		for(int i = 1; i < gsize; ++i) gindex[i].x = gindex[i-1].x + gindex[i-1].y;
		for(int i = 1; i < gsize; ++i) gindex[i].y = gindex[i  ].x + gindex[i  ].y;
		for( int iz = 0; iz < zdim_; ++iz) for( int iy = 0; iy < ydim_; ++iy) for( int ix = 0; ix < xdim_; ++ix) {
					int const ixl = (int)numeric::max(      0 ,(int)ix-1 );
					int const ixu =       numeric::min(xdim_-1u,     ix+1u);
					int const ig0 = xdim_*iy+xdim_*ydim_*iz;
					gstripe[ix+ig0].x = gindex[ixl+ig0].x;
					gstripe[ix+ig0].y = gindex[ixu+ig0].y;
				}
		grid_stripe_ = gstripe;
		// for(int iz = 0; iz < zdim_; ++iz) for(int iy = 0; iy < ydim_; ++iy) for(int ix = 0; ix < xdim_; ++ix) {
		//       int i = ix+xdim_*iy+xdim_*ydim_*iz;
		//       TR<<ix<<" "<<iy<<" "<<iz<<" "<<I(3,gindex[i].x)<<" "<<I(3,gindex[i].y) <<" "<<I(3,grid_stripe_[i].x)<<" "<<I(3,grid_stripe_[i].y)<<std::endl;
		//     }
		float4 *gatom = new float4[natom_+4]; // space for 4 overflow atoms
		for(int i=0;i<4;++i) {gatom[natom_+i].x=9e9;gatom[natom_+i].y=9e9;gatom[natom_+i].z=9e9;gatom[natom_+i].w=9e9;}
		ushort *gridc = new ushort[gsize];
		for(int i = 0; i < gsize; ++i) gridc[i] = 0;
		for(int i = 1; i <= natom_; ++i) {
			int const ix = (int)((atoms[i].x()-xmn/*+FUDGE*/)/grid_size_);
			int const iy = (int)((atoms[i].y()-ymn/*+FUDGE*/)/grid_size_);
			int const iz = (int)((atoms[i].z()-zmn/*+FUDGE*/)/grid_size_);
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			int const idx = gindex[ig].x + gridc[ig];
			gatom[ idx ].x = atoms[i].x()-xmn/*+FUDGE*/;
			gatom[ idx ].y = atoms[i].y()-ymn/*+FUDGE*/;
			gatom[ idx ].z = atoms[i].z()-zmn/*+FUDGE*/;
			gatom[ idx ].w = meta.size()>0 ? meta[i] : i;
			++(gridc[ig]);
		}
		grid_atoms_ = gatom;
		translation_.x() =/* FUDGE*/ - xmn;
		translation_.y() =/* FUDGE*/ - ymn;
		translation_.z() =/* FUDGE*/ - zmn;
		xmx_ = xmx-xmn/*+FUDGE*/+grid_size_;
		ymx_ = ymx-ymn/*+FUDGE*/+grid_size_;
		zmx_ = zmx-zmn/*+FUDGE*/+grid_size_;
		// for(int iz = 0; iz < zdim(); ++iz) for(int iy = 0; iy < ydim(); ++iy) for(int ix = 0; ix < xdim(); ++ix) {
		//       int i = ix+xdim_*iy+xdim_*ydim_*iz;
		//       TR<<"GRID CELL "<<ix<<" "<<iy<<" "<<iz<<std::endl;
		//       for(int ig = gindex[i].x; ig < gindex[i].y; ++ig) {
		//       TR<<F(7,3,gatom[ig].x)<<" "<<F(7,3,gatom[ig].y)<<" "<<F(7,3,gatom[ig].z)<<std::endl;
		//     }
		//   }
		delete [] gridc;
		delete [] gindex;
	}
	virtual ~xyzStripeHashWithMeta() {
		if(grid_atoms_)  delete grid_atoms_;
		if(grid_stripe_) delete grid_stripe_;
	}

	    xyz_iterator     xyz_begin() const { return     xyz_iterator(grid_atoms_       ) ; }
	    xyz_iterator     xyz_end()   const { return     xyz_iterator(grid_atoms_+natom_) ; }
	   meta_iterator    meta_begin() const { return    meta_iterator(grid_atoms_       ) ; }
	   meta_iterator    meta_end()   const { return    meta_iterator(grid_atoms_+natom_) ; }
	xyzmeta_iterator xyzmeta_begin() const { return xyzmeta_iterator(grid_atoms_       ) ; }
	xyzmeta_iterator xyzmeta_end()   const { return xyzmeta_iterator(grid_atoms_+natom_) ; }

	// neighbor_iterator neighbor_begin( xyzVector<T> v ) const { return neighbor_iterator(*this,v); }
	// neighbor_iterator const & neighbor_end() const { return neighbor_end_; }

	bool sanity_check() const {
		using namespace ObjexxFCL::format;
		for(int ix = 0; ix < xdim_; ++ix) {
			for(int iy = 0; iy < ydim_; ++iy) {
				for(int iz = 0; iz < zdim_; ++iz) {
					//std::cout << ix << " " << iy << " " << iz << endl;
					ushort const ig  = ix+xdim_*iy+ydim_*xdim_*iz;
					ushort const igl = grid_stripe_[ig].x;
					ushort const igu = grid_stripe_[ig].y;
					for(int i = igl; i < igu; ++i) {
						// float const & x(grid_atoms_[i].x);
						float const & y(grid_atoms_[i].y);
						float const & z(grid_atoms_[i].z);
					 // if(i==igl) std::cout << endl;
						// bool xc = grid_size_*(float)ix <= x && x <= grid_size_*(float)(ix+1);
						bool yc = grid_size_*(float)iy <= y && y <= grid_size_*(float)(iy+1);
						bool zc = grid_size_*(float)iz <= z && z <= grid_size_*(float)(iz+1);
						if(/*!xc||*/!yc||!zc) utility_exit_with_message("INSANE!");
						//std::cout<<I(2,ix)<<" "<<I(2,iy)<<" "<<I(2,iz)<<" "<<F(8,3,x)<<" "<<F(8,3,y)<<" "<<F(8,3,z)<<" "<<xc<<" "<<yc<<" "<<zc<<std::endl;
					}
				}
				return true;
			}
		}
		return true;
	}

	inline
	int
	nbcount(
		xyzVector<T> const & v_in
	) const {
		Vec const v = v_in+translation_;
		T x = v.x(); T y = v.y(); T z = v.z();
		if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return 0; // worth it iff
		if( x > xmx_ || y > ymx_ || z > zmx_ ) return 0;                      // worth it iff
		int count = 0;
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
				int const igl = grid_stripe_[ig].x;
				int const igu = grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					float4 const a2 = grid_atoms_[i];
					float const d2 = (x-a2.x)*(x-a2.x) + (y-a2.y)*(y-a2.y) + (z-a2.z)*(z-a2.z);
					if( d2 <= grid_size2_ ) {
						++count;
					}
				}
			}
		}
		return count;
	}

	template<typename V>
	inline
	void
	visit(
		Vec const & v_in,
		T const & m_in,
		V & visitor
	) const {
		Vec const v = v_in+translation_;
		T x = v.x(); T y = v.y(); T z = v.z();
		if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return; // worth it iff
		if( x > xmx_        || y > ymx_        || z > zmx_        ) return;                      // worth it iff
		int const ix   = (x<0) ? 0 : numeric::min(xdim_-1, static_cast<int>(x/grid_size_));
		int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
		int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
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
				int const igl = grid_stripe_[ig].x;
				int const igu = grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					Vec const & c = *((Vec*)(grid_atoms_+i));
					T   const & m = (grid_atoms_+i)->w;
					float const d2 = (x-c.x())*(x-c.x()) + (y-c.y())*(y-c.y()) + (z-c.z())*(z-c.z());
					if( d2 <= grid_size2_ ) {
						visitor.visit(v,m_in,c,m,d2);
					}
				}
			}
		}
	}

	template<typename V>
	inline
	void
	visit_lax(
		Vec const & v_in,
		T const & m_in,
		V & visitor
	) const {
		Vec const v = v_in+translation_;
		T x = v.x(); T y = v.y(); T z = v.z();
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
				int const igl = grid_stripe_[ig].x;
				int const igu = grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					Vec const & c = *((Vec*)(grid_atoms_+i));
					T   const & m = (grid_atoms_+i)->w;
					visitor.visit(v,m_in,c,m);
				}
			}
		}
	}


	inline float4  const * grid_atoms() const { return grid_atoms_; }
	inline ushort2 const * grid_stripe() const { return grid_stripe_; }
	inline int natom() const { return natom_; }
	inline int xdim () const { return  xdim_; }
	inline int ydim () const { return  ydim_; }
	inline int zdim () const { return  zdim_; }
	inline float  grid_size() const { return  grid_size_; }
	inline const numeric::xyzVector<T> translation() const { return translation_; }

private:
	xyzStripeHashWithMeta();
	T const grid_size_,grid_size2_;
	int natom_;
	float4  const * grid_atoms_;
	ushort2 const * grid_stripe_;
	int xdim_,ydim_,zdim_;
	T xmx_,ymx_,zmx_;
	//numeric::xyzMatrix<T> rotation_;
	numeric::xyzVector<T> translation_;
	// neighbor_iterator neighbor_end_;
};

class xyzStripeHashWithMeta_float : public xyzStripeHashWithMeta<float>
{

};

} // namespace hashing
} // namespace geometry
} // namespace numeric

#endif


	// struct neighbor_iterator : public std::iterator<std::input_iterator_tag,T> {
	// 	neighbor_iterator( xyzStripeHashWithMeta<T> const & h ): h_(h) {}
	// 	neighbor_iterator( xyzStripeHashWithMeta<T> const & h, numeric::xyzVector<T> const & v_in):
	// 		h_(h),
	// 		x  (v_in.x()+h_.translation_.x()),
	// 		y  (v_in.y()+h_.translation_.y()),
	// 		z  (v_in.z()+h_.translation_.z()),
	// 		ix ( (x<0) ? 0 : numeric::min(h_.xdim_-1,(int)(x/h_.grid_size_)) ),
	// 		iy0( (y<0) ? 0 : y/h_.grid_size_ ),
	// 		iz0( (z<0) ? 0 : z/h_.grid_size_ ),
	// 		iyl( numeric::max(0,iy0-1) ),
	// 		izl( numeric::max(0,iz0-1) ),
	// 		iyu( numeric::min((int)h_.ydim_,iy0+2) ),
	// 		izu( numeric::min((int)h_.zdim_,(int)iz0+2) )
	// 	{
	// 		// if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return 0; // worth it iff
	// 		// if( x > xmx_ || y > ymx_ || z > zmx_ ) return 0;                      // worth it iff
	// 		iy = iyl-1;
	// 		iz=9999999,igu=0,i=9999999;
	// 		// return count;
	// 		++(*this);
	// 	}
	// 	neighbor_iterator & operator++() {
	// 		using ObjexxFCL::format::I;
	// 		while(iy < iyu){
	// 			if(iz >= izu){
	// 				++iy;
	// 				iz = izl-1;
	// 				i=9999999,igu=0;
	// 			} else {
	// 				if( i >= igu ){
	// 					++iz;
	// 					ig = ix+h_.xdim_*iy+h_.xdim_*h_.ydim_*iz;
	// 					igl = h_.grid_stripe_[ig].x;
	// 					igu = h_.grid_stripe_[ig].y;
	// 					i = igl-1;
	// 				} else {
	// 					++i;
	// 					float4 const & a2 = h_.grid_atoms_[i];
	// 					float const d2 = (x-a2.x)*(x-a2.x) + (y-a2.y)*(y-a2.y) + (z-a2.z)*(z-a2.z);
	// 					if( d2 <= h_.grid_size2_ ){
	// 						std::cout << "iter "
	// 						          << I(3,iyl) << " " << I(3,iy) << " " << I(3,iyu) << "   "
	// 						          << I(3,izl) << " " << I(3,iz) << " " << I(3,izu) << "   "
	// 						          << I(5,igl) << " " << I(5,i-1) << " " << I(5,igu) << "   "
	// 						          << std::endl;
	// 						return *this;
	// 					}
	// 				}
	// 			}
	// 		}
	// 		end();
	// 		return *this;
	// 	}
	// 	const VecandVal & operator*() { return *((VecandVal const *)(grid_atoms_+i)); }
	// 	// neighbor_iterator & operator=(const neighbor_iterator& r) { p_ = r.p_; return *this; }
	// 	bool operator!=(neighbor_iterator const & r) const { return (i != r.i); }
	// 	bool operator==(neighbor_iterator const & r) const { return (i == r.i); }
	// 	void end() { i = h_.natom(); }
	// private:
	// 	xyzStripeHashWithMeta const & h_;
	// 	T x,y,z;
	// 	int ix,iy0,iz0,iyl,izl,iyu,izu;
	// 	int iy,iz,ig,igl,igu,i;
	// };
