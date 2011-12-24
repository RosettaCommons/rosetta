#include <apps/pilot/will/pca_simple.hh>
#include <apps/pilot/will/gpu_bit_utils.hh>

template<typename T, typename M>
class xyzStripeHash {
public:
  //typedef struct { T x,y,z,w; } float4;
  //typedef struct { ushort x,y; } ushort2;
private:
  T const grid_size_;
  uint natom_;
  float4  const * grid_atoms_;
  ushort2 const * grid_stripe_;
  uint xdim_,ydim_,zdim_;
  //numeric::xyzMatrix<Real> rotation_;
  //numeric::xyzVector<Real> translation_;
public:
  xyzStripeHash( T grid_size ) : grid_size_(grid_size), grid_atoms_(NULL), grid_stripe_(NULL)//,
                           /*rotation_(numeric::xyzMatrix<Real>::identity()), translation_(T(0),T(0),T(0))*/ {}
  xyzStripeHash( T grid_size,
           vector1<numeric::xyzVector<T> > const & atoms,
           vector1<M> const & meta
           ) : grid_size_(grid_size), grid_atoms_(NULL), grid_stripe_(NULL)//,
               /*rotation_(numeric::xyzMatrix<Real>::identity()), translation_(T(0),T(0),T(0))*/
  {
    init(atoms,meta);
  }

  void init(
            vector1<numeric::xyzVector<T> > const & atoms,
            vector1<M> const & meta
            )
  {
    if( sizeof(T) < sizeof(M) ) utility_exit_with_message("octree metadata must fit in sizeof(T)!");
    if( atoms.size() != meta.size() ) utility_exit_with_message("must be metadata for each point!");
    if( atoms.size() > 65535 ) utility_exit_with_message("xyzStripeHash con only handle < 65535 atoms!");

#define FUDGE 0.0f

    natom_ = atoms.size();

    T xmn= 9e9,ymn= 9e9,zmn= 9e9;
    T xmx=-9e9,ymx=-9e9,zmx=-9e9;
    for(uint i = 1; i <= natom_; ++i) {
      xmn = numeric::min(xmn,atoms[i].x());
      ymn = numeric::min(ymn,atoms[i].y());
      zmn = numeric::min(zmn,atoms[i].z());
      xmx = numeric::max(xmx,atoms[i].x());
      ymx = numeric::max(ymx,atoms[i].y());
      zmx = numeric::max(zmx,atoms[i].z());
    }
    //TR<<xmx-xmn<<" "<<ymx-ymn<<" "<<zmx-zmn<<std::endl;
    // for(uint i = 0; i < natom_; ++i) {
    //   atoms[i].x -= xmn-FUDGE;
    //   atoms[i].y -= ymn-FUDGE;
    //   atoms[i].z -= zmn-FUDGE;
    // }
    xdim_ = ceil((xmx-xmn+0.0001)/grid_size_);
    ydim_ = ceil((ymx-ymn+0.0001)/grid_size_);
    zdim_ = ceil((zmx-zmn+0.0001)/grid_size_);
    assert(xdim_ < 9999); assert(ydim_ < 9999); assert(zdim_ < 9999);
    uint const gsize = xdim_*ydim_*zdim_;
    ushort2 *gindex  = new ushort2[gsize];
    ushort2 *gstripe = new ushort2[gsize];
    for(Size i = 0; i < gsize; ++i) { gindex[i].y = 0; gindex[i].x = 0; }
    //TR<<"atom "<<natom_<<" grid1 "<<xdim_*ydim_*zdim_<<" "<<xdim_<<" "<<ydim_<<" "<<zdim_<<std::endl;

    for(Size i = 1; i <= natom_; ++i) {
      int ix = (atoms[i].x()-xmn+FUDGE)/grid_size_;
      int iy = (atoms[i].y()-ymn+FUDGE)/grid_size_;
      int iz = (atoms[i].z()-zmn+FUDGE)/grid_size_;
      assert(ix >= 0)1; assert(iy >= 0); assert(iz >= 0); assert(ix < xdim_); assert(iy < ydim_); assert(iz < zdim_);
      int ig = ix+xdim_*iy+xdim_*ydim_*iz;
      assert(ig>=0);assert(ig<9999999);
      ++(gindex[ig].y);
    }
    for(uint i = 1; i < gsize; ++i) gindex[i].x = gindex[i-1].x + gindex[i-1].y;
    for(uint i = 1; i < gsize; ++i) gindex[i].y = gindex[i  ].x + gindex[i  ].y;
    for( int iz = 0; iz < zdim_; ++iz) for( int iy = 0; iy < ydim_; ++iy) for( int ix = 0; ix < xdim_; ++ix) {
          uint const ixl = (uint)numeric::max(     0 ,(int)ix-1 );
          uint const ixu =       numeric::min(xdim_-1u,     ix+1u);
          uint const ig0 = xdim_*iy+xdim_*ydim_*iz;
          gstripe[ix+ig0].x = gindex[ixl+ig0].x;
          gstripe[ix+ig0].y = gindex[ixu+ig0].y;
        }
    grid_stripe_ = gstripe;
    // for(uint iz = 0; iz < zdim_; ++iz) for(uint iy = 0; iy < ydim_; ++iy) for(uint ix = 0; ix < xdim_; ++ix) {
    //       uint i = ix+xdim_*iy+xdim_*ydim_*iz;
    //       TR<<ix<<" "<<iy<<" "<<iz<<" "<<I(3,gindex[i].x)<<" "<<I(3,gindex[i].y) <<" "<<I(3,grid_stripe_[i].x)<<" "<<I(3,grid_stripe_[i].y)<<std::endl;
    //     }
    float4 *gatom = new float4[natom_+4]; // space for 4 overflow atoms
    for(uint i=0;i<4;++i) {gatom[natom_+i].x=9e9;gatom[natom_+i].y=9e9;gatom[natom_+i].z=9e9;gatom[natom_+i].w=9e9;}
    ushort *gridc = new ushort[gsize];
    for(Size i = 0; i < gsize; ++i) gridc[i] = 0;
    for(Size i = 1; i <= natom_; ++i) {
      uint const ix = (atoms[i].x()-xmn+FUDGE)/grid_size_;
      uint const iy = (atoms[i].y()-ymn+FUDGE)/grid_size_;
      uint const iz = (atoms[i].z()-zmn+FUDGE)/grid_size_;
      uint const ig = ix+xdim_*iy+xdim_*ydim_*iz;
      uint const idx = gindex[ig].x + gridc[ig];
      gatom[ idx ].x = atoms[i].x()-xmn+FUDGE;
      gatom[ idx ].y = atoms[i].y()-ymn+FUDGE;
      gatom[ idx ].z = atoms[i].z()-zmn+FUDGE;
      gatom[ idx ].w = *((float*)(&meta[i]));
      ++(gridc[ig]);
    }
    grid_atoms_ = gatom;
    // translation_.x() = FUDGE - xmn;
    // translation_.y() = FUDGE - ymn;
    // translation_.z() = FUDGE - zmn;
    // for(uint iz = 0; iz < zdim(); ++iz) for(uint iy = 0; iy < ydim(); ++iy) for(uint ix = 0; ix < xdim(); ++ix) {
    //       uint i = ix+xdim_*iy+xdim_*ydim_*iz;
    //       TR<<"GRID CELL "<<ix<<" "<<iy<<" "<<iz<<std::endl;
    //       for(Size ig = gindex[i].x; ig < gindex[i].y; ++ig) {
    //       TR<<F(7,3,gatom[ig].x)<<" "<<F(7,3,gatom[ig].y)<<" "<<F(7,3,gatom[ig].z)<<std::endl;
    //     }
    //   }
    delete gridc,gindex;
  }
  virtual ~xyzStripeHash() {
    if(grid_atoms_)  delete grid_atoms_;
    if(grid_stripe_) delete grid_stripe_;
  }

  bool sanity_check() {
    for(int ix = 0; ix < xdim_; ++ix) {
      for(int iy = 0; iy < ydim_; ++iy) {
        for(int iz = 0; iz < zdim_; ++iz) {
          //cout << ix << " " << iy << " " << iz << endl;
          ushort const ig  = ix+xdim_*iy+ydim_*xdim_*iz;
          ushort const igl = grid_stripe_[ig].x;
          ushort const igu = grid_stripe_[ig].y;
          for(int i = igl; i < igu; ++i) {
            float const & x(grid_atoms_[i].x);
            float const & y(grid_atoms_[i].y);
            float const & z(grid_atoms_[i].z);
            if(i==igl) cout << endl;
            bool xc = grid_size_*(float)ix <= x && x <= grid_size_*(float)(ix+1);
            bool yc = grid_size_*(float)iy <= y && y <= grid_size_*(float)(iy+1);
            bool zc = grid_size_*(float)iz <= z && z <= grid_size_*(float)(iz+1);
            if(/*!xc||*/!yc||!zc) utility_exit_with_message("INSANE!");
            cout<<I(2,ix)<<" "<<I(2,iy)<<" "<<I(2,iz)<<" "<<F(8,3,x)<<" "<<F(8,3,y)<<" "<<F(8,3,z)<<" "<<xc<<" "<<yc<<" "<<zc<<std::endl;
          }
        }
        return true;
      }
    }
    return true;
  }

  inline Size const nbcount( float x, float y, float z ) {
    Size count = 0u;
    uint const ix  = (x<0) ? 0u : numeric::min(xdim_,(uint)(x/grid_size_));
    int const iy0  = (y<0) ? 0u : y/grid_size_;
    int const iz0  = (z<0) ? 0u : z/grid_size_;
    uint const iyl = numeric::max(0,iy0-1);
    uint const izl = numeric::max(0,iz0-1);
    uint const iyu = numeric::min((int)ydim_,iy0+2);
    uint const izu = numeric::min((int)zdim_,(int)iz0+2);
    for(uint iy = iyl; iy < iyu; ++iy) {
      for(uint iz = izl; iz < izu; ++iz) {
        uint const ig = ix+xdim_*iy+xdim_*ydim_*iz;
        uint const igl = grid_stripe_[ig].x;
        uint const igu = grid_stripe_[ig].y;
        for(uint i = igl; i < igu; ++i) {
          float4 const a2 = grid_atoms_[i];
          float const d2 = sqr(x-a2.x) + sqr(y-a2.y) + sqr(z-a2.z);
          if( d2 <= grid_size_*grid_size_ ) {
            ++count;
          }
        }
      }
    }
    return count;
  }
  inline float const ljenergy( float x, float y, float z ) {
    float e = 0.0f;;
    uint const ix  = (x<0) ? 0u : numeric::min(xdim_,(uint)(x/grid_size_));
    int const iy0  = (y<0) ? 0u : y/grid_size_;
    int const iz0  = (z<0) ? 0u : z/grid_size_;
    uint const iyl = numeric::max(0,iy0-1);
    uint const izl = numeric::max(0,iz0-1);
    uint const iyu = numeric::min((int)ydim_,iy0+2);
    uint const izu = numeric::min((int)zdim_,(int)iz0+2);
    for(uint iy = iyl; iy < iyu; ++iy) {
      for(uint iz = izl; iz < izu; ++iz) {
        uint const ig = ix+xdim_*iy+xdim_*ydim_*iz;
        uint const igl = grid_stripe_[ig].x;
        uint const igu = grid_stripe_[ig].y;
        for(uint i = igl; i < igu; ++i) {
          float4 const a2 = grid_atoms_[i];
          float const d2 = sqr(x-a2.x) + sqr(y-a2.y) + sqr(z-a2.z);
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
  inline uint const natom() const { return natom_; }
  inline uint const xdim () const { return  xdim_; }
  inline uint const ydim () const { return  ydim_; }
  inline uint const zdim () const { return  zdim_; }
  inline float const & grid_size() const { return  grid_size_; }

};



inline short const  short_min( short const a,  short const b) { return (a < b) ? a : b; }
inline short const  short_max( short const a,  short const b) { return (a > b) ? a : b; }
inline short const ushort_min(ushort const a, ushort const b) { return (a < b) ? a : b; }
inline short const ushort_max(ushort const a, ushort const b) { return (a > b) ? a : b; }

float
octree_test( float4  const * gatom,
             ushort2 const * gstripe,
             float   const * gsize_in,
             uint8   const * gdim,
             uint tmp
             //        float         * output
             ){
  float e = 0.0f;
  ushort const xdim = gdim[0].x;
  ushort const ydim = gdim[0].y;
  ushort const zdim = gdim[0].z;
  float const gsize = gsize_in[0];
  float const gsize2 = gsize*gsize;

  for(ushort ia = ((ushort)0); ia < gdim[0].s4; ++ia) {
    float tmp = 0.0f;
    float4 const a = gatom[ia];// + (10.0f*((float)get_global_id(0))/((float)get_global_size(0)));
    short const ix   = ((a.x < 0.0000001f) ? ((short)0) : short_min(xdim,(ushort)(a.x/gsize)));
    short const iy0  = ((a.y < 0.0000001f) ? ((short)0) : a.y/gsize);
    short const iz0  = ((a.z < 0.0000001f) ? ((short)0) : a.z/gsize);
    ushort const iyl = short_max((int)0,iy0-(int)1);
    ushort const izl = short_max((int)0,iz0-(int)1);
    ushort const iyu = short_min(ydim,(uint)iy0+((ushort)2));
    ushort const izu = short_min(zdim,(ushort)iz0+((ushort)2));
    for(ushort iy = iyl; iy < iyu; ++iy) {
      for(ushort iz = izl; iz < izu; ++iz) {
        ushort const ig = ix+xdim*iy+ydim*xdim*iz;
        ushort const igl = short_max(((ushort)gstripe[ig].x),ia+((ushort)1u));
        ushort const igu = short_max(((ushort)gstripe[ig].y),ia+((ushort)1u));
        if( igu == ia ) continue;
        ushort i = igl;
        while(i < igu) {
          float4 const a0 = gatom[i+0];
          float const d2a = mad(a.x-a0.x,a.x-a0.x,mad(a.y-a0.y,a.y-a0.y,sqr(a.z-a0.z)));
          float const ra = (d2a < 4.00000f || d2a > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2a));
          float const ra2 = ra*ra;
          float const ra3 = ra*ra2;
          float const ra6 = ra2*ra3;
          float const sa = mad(ra6,ra6,-ra6) + exp(-d2a);
          tmp += sa;
          i = i+((ushort)1);
        }
      }
    }
    e += tmp;
  }
  return e;
}
