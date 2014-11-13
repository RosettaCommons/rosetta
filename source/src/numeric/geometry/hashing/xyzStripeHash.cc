// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <numeric/geometry/hashing/xyzStripeHash.hh>

namespace numeric {
namespace geometry {
namespace hashing {

using std::cout;
using std::cerr;
using std::endl;

xyzStripeHash::xyzStripeHash(
	float grid_size,
    utility::vector1<Ball> const & balls
):
	grid_size_(grid_size),
	grid_size2_(grid_size*grid_size),
	grid_balls_(NULL),
	grid_stripe_(NULL)//,
	// neighbor_end_(*this)
{
	if(balls.size()) init(balls);
 }

void
xyzStripeHash::init(
	utility::vector1<Ball> const & balls
){
	if( balls.size() < 1 || balls.size() > 65535 ){
		// std::cout << std::endl << "nballs: " << balls.size() << std::endl;
		// utility_exit_with_message("xyzStripeHash requires 0 < N < 65535 balls!");
		nballs_ = 0;
		grid_size_ = 0;
		grid_size2_ = 0;
		xmx_ = 0;
		ymx_ = 0;
		zmx_ = 0;
		return;
	}
	if(grid_size_ == 0.0 ){
		BOOST_FOREACH(Ball const & b, balls) grid_size_ = max(grid_size_,2.0f*b.radius());
		grid_size2_ = grid_size_*grid_size_;
	}
	nballs_ = balls.size();
	// neighbor_end_.end();

	if( grid_size_ <= 0.0 ) utility_exit_with_message("grid_size_ <= 0");

	float xmn= 9e9,ymn= 9e9,zmn= 9e9;
	float xmx=-9e9,ymx=-9e9,zmx=-9e9;
	for(int i = 1; i <= nballs_; ++i) {
		xmn = numeric::min(xmn,balls[i].x());
		ymn = numeric::min(ymn,balls[i].y());
		zmn = numeric::min(zmn,balls[i].z());
		xmx = numeric::max(xmx,balls[i].x());
		ymx = numeric::max(ymx,balls[i].y());
		zmx = numeric::max(zmx,balls[i].z());
	}

	xdim_ = static_cast<int>((xmx-xmn+0.0001)/grid_size_+0.999999);
	ydim_ = static_cast<int>((ymx-ymn+0.0001)/grid_size_+0.999999);
	zdim_ = static_cast<int>((zmx-zmn+0.0001)/grid_size_+0.999999);
	assert(xdim_ < 9999); assert(ydim_ < 9999); assert(zdim_ < 9999);
	int const gsize = xdim_*ydim_*zdim_;
	ushort2 *gindex  = new ushort2[gsize];
	ushort2 *gstripe = new ushort2[gsize];
	for(int i = 0; i < gsize; ++i) { gindex[i].y = 0; gindex[i].x = 0; }
	//TR<<"atom "<<nballs_<<" grid1 "<<xdim_*ydim_*zdim_<<" "<<xdim_<<" "<<ydim_<<" "<<zdim_<<std::endl;

	for(int i = 1; i <= nballs_; ++i) {
		int ix = static_cast<int>((balls[i].x()-xmn/*+FUDGE*/)/grid_size_);
		int iy = static_cast<int>((balls[i].y()-ymn/*+FUDGE*/)/grid_size_);
		int iz = static_cast<int>((balls[i].z()-zmn/*+FUDGE*/)/grid_size_);
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
	Ball *gatom = new Ball[nballs_+4]; // space for 4 overflow balls
	for(int i=0;i<4;++i) {gatom[nballs_+i].x()=9e9;gatom[nballs_+i].y()=9e9;gatom[nballs_+i].z()=9e9;}
	ushort *gridc = new ushort[gsize];
	for(int i = 0; i < gsize; ++i) gridc[i] = 0;
	for(int i = 1; i <= nballs_; ++i) {
		int const ix = static_cast<int>((balls[i].x()-xmn/*+FUDGE*/)/grid_size_);
		int const iy = static_cast<int>((balls[i].y()-ymn/*+FUDGE*/)/grid_size_);
		int const iz = static_cast<int>((balls[i].z()-zmn/*+FUDGE*/)/grid_size_);
		int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
		int const idx = gindex[ig].x + gridc[ig];
		gatom[ idx ].x() = balls[i].x()-xmn/*+FUDGE*/;
		gatom[ idx ].y() = balls[i].y()-ymn/*+FUDGE*/;
		gatom[ idx ].z() = balls[i].z()-zmn/*+FUDGE*/;
		gatom[ idx ].radius( balls[i].radius() );
		gatom[ idx ].resid_ = balls[i].resid_;
		gatom[ idx ].atomno_ = balls[i].atomno_;
		++(gridc[ig]);
	}
	grid_balls_ = gatom;
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
	delete gridc;
	delete gindex;
 }

bool xyzStripeHash::sanity_check() const {
	using namespace ObjexxFCL::format;
	for(int ix = 0; ix < xdim_; ++ix) {
		for(int iy = 0; iy < ydim_; ++iy) {
			for(int iz = 0; iz < zdim_; ++iz) {
				//std::cout << ix << " " << iy << " " << iz << endl;
				ushort const ig  = ix+xdim_*iy+ydim_*xdim_*iz;
				ushort const igl = grid_stripe_[ig].x;
				ushort const igu = grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					// float const & x(grid_balls_[i].x);
					float const & y(grid_balls_[i].y());
					float const & z(grid_balls_[i].z());
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

int
xyzStripeHash::nbcount( Vec const & v_in ) const {
	Vec const v = v_in+translation_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return 0; // worth it iff
	if( x > xmx_ || y > ymx_ || z > zmx_ ) return 0;                      // worth it iff
	int count = 0;
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
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
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 <= grid_size2_ ) {
					++count;
				}
			}
		}
	}
	return count;
 }

int
xyzStripeHash::nbcount_raw( Vec const & v ) const {
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return 0; // worth it iff
	if( x > xmx_ || y > ymx_ || z > zmx_ ) return 0;                      // worth it iff
	int count = 0;
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
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 <= grid_size2_ ) {
					++count;
				}
			}
		}
	}
	return count;
 }

bool
xyzStripeHash::clash( Vec const & v_in ) const {
	Vec const v = v_in+translation_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
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
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 < grid_size2_ ) {
					return true;
				}
			}
		}
	}
	return false;
 }

float
xyzStripeHash::clash_amount( Vec const & v_in ) const {
	Vec const v = v_in+translation_;
	float clash_amount = grid_size2_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
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
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 < grid_size2_ ) {
					cout << grid_size_-sqrt(d2) << " " << grid_size2_ << " " << a2.resi() << " " << a2.atomno() << " " <<  endl;
					clash_amount = numeric::min(d2,clash_amount);
			}
		}
	}
 }
	clash_amount = grid_size_ - sqrt(clash_amount);
	return clash_amount;
 }

bool
xyzStripeHash::clash_not_resid( Vec const & v_in, int const & resid, int const & resid2 ) const {
	Vec const v = v_in+translation_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
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
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 < grid_size2_ && (int)a2.resid_ != resid && (int)a2.resid_ != resid2 ) {
					return true;
				}
			}
		}
	}
	return false;
 }

bool
xyzStripeHash::clash_raw( Vec const & v ) const {
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1, static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
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
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 < grid_size2_ ) {
					return true;
				}
			}
		}
	}
	return false;
 }

int
xyzStripeHash::clash_check_ball( Ball const & b_in ) const {
	// Input ball coordinates are in global coordinate space,
	// translate ball into hash frame.
	Vec b_in_hash_space = b_in.xyz() + translation_;
	float x = b_in_hash_space.x(); float y = b_in_hash_space.y(); float z = b_in_hash_space.z();

	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1, static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
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
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				
				if( d2 < (a2.radius() + b_in.radius()) ) {
					// Return type is 1-based ball number.
					return i + 1;
				}
			}
		}
	}

	return 0;
 }

bool xyzStripeHash::clash_check_residue_pairs(
		utility::vector1<Ball> const & test_balls,
		std::map<Size, Size> & residue_pairs
) const
{
	bool found_clash = false;

	for (Size i = 1; i <= test_balls.size(); i++)
	{
		if (residue_pairs.count(test_balls[i].resi()))
		{
			continue;
		}

		int clashing_ball = clash_check_ball(test_balls[i]);

		if (clashing_ball)
		{
			residue_pairs[test_balls[i].resi()] = ball(clashing_ball).resi();
			found_clash = true;
		}
	}

	return found_clash;
}

void
xyzStripeHash::fill_pairs(
	Vec const & v_in,
	int const & ir,
	utility::vector1<std::pair<int,int> > & pairs,
	float maxd2
) const {
	if(0.0==maxd2) maxd2=grid_size2_;
	Vec const v = v_in+translation_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
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
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 <= maxd2 ) {
					pairs.push_back(std::make_pair(ir,a2.resid_));
				}
			}
		}
	}
 }

std::string xyzStripeHash::debug_pdb(Xform const & x) const {
	using namespace ObjexxFCL::format;
	std::ostringstream out;
	int atomno = 0;
	for(int i = 0; i < nballs_; ++i){
		Ball const b = grid_balls_[i];
		Vec v = x*(b.xyz() - translation_);
		++atomno;
		out<<"ATOM  "<<I(5,atomno)<<' '<<" XSH"<<' '<<"XSH"<<' '<<'~'<<I(4,atomno/100)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
	}
	return out.str();

}



 } // namespace hashing
 } // namespace geometry
 } // namespace numeric

