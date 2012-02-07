#ifndef INCLUDED_apps_pilot_will_xyzStripeHash_hh
#define INCLUDED_apps_pilot_will_xyzStripeHash_hh

#include <utility/vector1.hh>
#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/format.hh>

template<typename T>
class xyzStripeHash {
public:
  typedef struct { T x,y,z,w; } float4;
	//typedef unsigned int uint;
	typedef unsigned short ushort;
  typedef struct { ushort x,y; } ushort2;
public:
  T const grid_size_,grid_size2_;
  int natom_;
  float4  const * grid_atoms_;
  ushort2 const * grid_stripe_;
  int xdim_,ydim_,zdim_;
	T xmx_,ymx_,zmx_;
  //numeric::xyzMatrix<Real> rotation_;
  numeric::xyzVector<numeric::Real> translation_;
public:

  xyzStripeHash( T grid_size ) : grid_size_(grid_size), grid_size2_(grid_size*grid_size), grid_atoms_(NULL), grid_stripe_(NULL)//,
                           /*rotation_(numeric::xyzMatrix<Real>::identity()), translation_(T(0),T(0),T(0))*/ {}
  xyzStripeHash( T grid_size,
           utility::vector1<numeric::xyzVector<T> > const & atoms,
           utility::vector1<T> const & meta
           ) : grid_size_(grid_size), grid_atoms_(NULL), grid_stripe_(NULL)//,
               /*rotation_(numeric::xyzMatrix<Real>::identity()), translation_(T(0),T(0),T(0))*/
  {
    init(atoms,meta);
  }

  void init(
            utility::vector1<numeric::xyzVector<T> > const & atoms,
            utility::vector1<T> const & meta
            )
  {
    // if( sizeof(T) < sizeof(M) ) utility_exit_with_message("octree metadata must fit in sizeof(T)!");
    if( atoms.size() != meta.size() ) utility_exit_with_message("must be metadata for each point!");
    if( atoms.size() > 65535 ) utility_exit_with_message("xyzStripeHash con only handle < 65535 atoms!");

#define FUDGE 0.0f

    natom_ = atoms.size();

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
    //   atoms[i].x -= xmn-FUDGE;
    //   atoms[i].y -= ymn-FUDGE;
    //   atoms[i].z -= zmn-FUDGE;
    // }

		//std::cout << "xyzStripeHash: " << xmn << " "  << ymn << " "  << zmn << " " << xmx << " "  << ymx << " "  << zmx << std::endl;

    xdim_ = ceil((xmx-xmn+0.0001)/grid_size_);
    ydim_ = ceil((ymx-ymn+0.0001)/grid_size_);
    zdim_ = ceil((zmx-zmn+0.0001)/grid_size_);
    assert(xdim_ < 9999); assert(ydim_ < 9999); assert(zdim_ < 9999);
    int const gsize = xdim_*ydim_*zdim_;
    ushort2 *gindex  = new ushort2[gsize];
    ushort2 *gstripe = new ushort2[gsize];
    for(int i = 0; i < gsize; ++i) { gindex[i].y = 0; gindex[i].x = 0; }
    //TR<<"atom "<<natom_<<" grid1 "<<xdim_*ydim_*zdim_<<" "<<xdim_<<" "<<ydim_<<" "<<zdim_<<std::endl;

    for(int i = 1; i <= natom_; ++i) {
      int ix = (atoms[i].x()-xmn+FUDGE)/grid_size_;
      int iy = (atoms[i].y()-ymn+FUDGE)/grid_size_;
      int iz = (atoms[i].z()-zmn+FUDGE)/grid_size_;
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
      int const ix = (atoms[i].x()-xmn+FUDGE)/grid_size_;
      int const iy = (atoms[i].y()-ymn+FUDGE)/grid_size_;
      int const iz = (atoms[i].z()-zmn+FUDGE)/grid_size_;
      int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
      int const idx = gindex[ig].x + gridc[ig];
      gatom[ idx ].x = atoms[i].x()-xmn+FUDGE;
      gatom[ idx ].y = atoms[i].y()-ymn+FUDGE;
      gatom[ idx ].z = atoms[i].z()-zmn+FUDGE;
      gatom[ idx ].w = meta[i];
      ++(gridc[ig]);
    }
    grid_atoms_ = gatom;
    translation_.x() = FUDGE - xmn;
    translation_.y() = FUDGE - ymn;
    translation_.z() = FUDGE - zmn;
		xmx_ = xmx-xmn+FUDGE+grid_size_;
		ymx_ = ymx-ymn+FUDGE+grid_size_;
		zmx_ = zmx-zmn+FUDGE+grid_size_;				
    // for(int iz = 0; iz < zdim(); ++iz) for(int iy = 0; iy < ydim(); ++iy) for(int ix = 0; ix < xdim(); ++ix) {
    //       int i = ix+xdim_*iy+xdim_*ydim_*iz;
    //       TR<<"GRID CELL "<<ix<<" "<<iy<<" "<<iz<<std::endl;
    //       for(int ig = gindex[i].x; ig < gindex[i].y; ++ig) {
    //       TR<<F(7,3,gatom[ig].x)<<" "<<F(7,3,gatom[ig].y)<<" "<<F(7,3,gatom[ig].z)<<std::endl;
    //     }
    //   }
    delete gridc,gindex;
  }
  virtual ~xyzStripeHash() {
    if(grid_atoms_)  delete grid_atoms_;
    if(grid_stripe_) delete grid_stripe_;
  }

  bool sanity_check() const {
		using namespace ObjexxFCL::fmt;
    for(int ix = 0; ix < xdim_; ++ix) {
      for(int iy = 0; iy < ydim_; ++iy) {
        for(int iz = 0; iz < zdim_; ++iz) {
          //std::cout << ix << " " << iy << " " << iz << endl;
          ushort const ig  = ix+xdim_*iy+ydim_*xdim_*iz;
          ushort const igl = grid_stripe_[ig].x;
          ushort const igu = grid_stripe_[ig].y;
          for(int i = igl; i < igu; ++i) {
            float const & x(grid_atoms_[i].x);
            float const & y(grid_atoms_[i].y);
            float const & z(grid_atoms_[i].z);
           // if(i==igl) std::cout << endl;
            bool xc = grid_size_*(float)ix <= x && x <= grid_size_*(float)(ix+1);
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

  inline int const nbcount( float x, float y, float z ) const {
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
  inline float const ljenergy( float x, float y, float z ) {
    float e = 0.0f;;
    int const ix  = (x<0) ? 0u : numeric::min(xdim_-1,(int)(x/grid_size_));
    int const iy0  = (y<0) ? 0u : y/grid_size_;
    int const iz0  = (z<0) ? 0u : z/grid_size_;
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
          if( d2 <= grid_size_*grid_size_ ) {
            float const r = 4.0 / sqrt(d2);
            float const r2 = r*r;
            float const r3 = r2*r;
            float const r6 = r3*r2;
            float const ljtmp = r6*r6 - r6;
            if( d2 >= 4.0f ) e += ljtmp;
          }
        }
      }
    }
    return e;
  }

	inline float4  const * grid_atoms() const { return grid_atoms_; }
	inline ushort2 const * grid_stripe() const { return grid_stripe_; }
	inline int const natom() const { return natom_; }
	inline int const xdim () const { return  xdim_; }
	inline int const ydim () const { return  ydim_; }
	inline int const zdim () const { return  zdim_; }
	inline float const & grid_size() const { return  grid_size_; }
	inline const numeric::xyzVector<numeric::Real> translation() const { return translation_; }

private:
	xyzStripeHash();
};



inline short const  short_min( short const a,  short const b) { return (a < b) ? a : b; }
inline short const  short_max( short const a,  short const b) { return (a > b) ? a : b; }
inline short const ushort_min(ushort const a, ushort const b) { return (a < b) ? a : b; }
inline short const ushort_max(ushort const a, ushort const b) { return (a > b) ? a : b; }



#endif
