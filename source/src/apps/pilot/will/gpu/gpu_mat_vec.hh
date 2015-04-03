#include <numeric/IOTraits.hh>
#include <iostream>
#include <iomanip>

typedef cl_float16 float16;
typedef cl_float8  float8;
typedef cl_float4  float4;
typedef cl_float3  float3;
typedef cl_float2  float2;
typedef cl_ushort2 ushort2;
typedef cl_uint8   uint8;

inline float native_sin(float x) { return std::sin(x); }
inline float native_cos(float x) { return std::cos(x); }
inline float native_divide(float a, float b) { return a/b; }
inline float native_sqrt(float a) { return sqrt(a); }
inline float mad(float a, float b, float c) { return a*b+c; }

inline float sqr(float const & x) {
  return x*x;
}


/////////////////////// VEC ///////////////////////

struct VEC {
  float x,y,z;
#ifdef __cplusplus
  VEC() {}
  VEC(float _x, float _y, float _z) { x=_x; y=_y; z=_z; }
  VEC(Vec v) {
    x = v.x();
    y = v.y();
    z = v.z();
  }
  Vec xyzVector() const { return Vec(x,y,z); }
#endif
};
inline struct VEC vec(float x, float y, float z) { VEC v; v.x=x; v.y=y; v.z=z; return v; }

/////////////////////// MAT ///////////////////////

struct MAT {
  float xx,xy,xz,
    yx,yy,yz,
    zx,zy,zz;
#ifdef __cplusplus
  MAT() {}
  MAT(Mat m) {
    xx = m.xx();
    yx = m.yx();
    zx = m.zx();
    xy = m.xy();
    yy = m.yy();
    zy = m.zy();
    xz = m.xz();
    yz = m.yz();
    zz = m.zz();
  }
  Mat xyzMatrix() const { return Mat::rows(xx,xy,xz,yx,yy,yz,zx,zy,zz); };
#endif
};
inline struct MAT rows(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz) {
  struct MAT m;
  m.xx=xx; m.xy=xy; m.xz=xz;
  m.yx=yx; m.yy=yy; m.yz=yz;
  m.zx=zx; m.zy=zy; m.zz=zz;
  return m;
}
inline struct MAT cols(float xx, float yx, float zx, float xy, float yy, float zy, float xz, float yz, float zz) {
  struct MAT m;
  m.xx=xx; m.xy=xy; m.xz=xz;
  m.yx=yx; m.yy=yy; m.yz=yz;
  m.zx=zx; m.zy=zy; m.zz=zz;
  return m;
}
inline struct MAT rows(struct VEC const rx, struct VEC const ry, struct VEC const rz) {
  struct MAT m;
  m.xx=rx.x; m.xy=rx.y; m.xz=rx.z;
  m.yx=ry.x; m.yy=ry.y; m.yz=ry.z;
  m.zx=rz.x; m.zy=rz.y; m.zz=rz.z;
  return m;
}
inline struct MAT cols(struct VEC const cx, struct VEC const cy, struct VEC const cz) {
  struct MAT m;
  m.xx=cx.x; m.xy=cy.x; m.xz=cz.x;
  m.yx=cx.y; m.yy=cy.y; m.yz=cz.y;
  m.zx=cx.z; m.zy=cy.z; m.zz=cz.z;
  return m;
}

/////////////////////// XFORM ///////////////////////

struct XFORM {
  struct MAT R;
  struct VEC t;
#ifdef __cplusplus
  // XFORM(core::pose::Pose const & p, Size rsd) {    
  // }
  core::kinematics::Stub stub() const { return core::kinematics::Stub(R.xyzMatrix(),t.xyzVector()); }
  void apply(core::pose::Pose & p) { xform_pose(p,stub());  }
  void aprev(core::pose::Pose & p) { xform_pose_rev(p,stub());  }
#endif
};
inline struct XFORM xform(struct MAT const R, struct VEC const t) { struct XFORM x; x.R = R;  x.t = t;  return x; }


////////////////////// operations //////////////////////


inline struct MAT multmm(struct MAT const a, struct MAT const b) {
  struct MAT c;
  c.xx = mad(a.xx,b.xx,mad(a.xy,b.yx,a.xz*b.zx));
  c.xy = mad(a.xx,b.xy,mad(a.xy,b.yy,a.xz*b.zy));
  c.xz = mad(a.xx,b.xz,mad(a.xy,b.yz,a.xz*b.zz));
  c.yx = mad(a.yx,b.xx,mad(a.yy,b.yx,a.yz*b.zx));
  c.yy = mad(a.yx,b.xy,mad(a.yy,b.yy,a.yz*b.zy));
  c.yz = mad(a.yx,b.xz,mad(a.yy,b.yz,a.yz*b.zz));
  c.zx = mad(a.zx,b.xx,mad(a.zy,b.yx,a.zz*b.zx));
  c.zy = mad(a.zx,b.xy,mad(a.zy,b.yy,a.zz*b.zy));
  c.zz = mad(a.zx,b.xz,mad(a.zy,b.yz,a.zz*b.zz));
  return c;
}
inline struct VEC multmv(struct MAT const a, struct VEC const b) {
  struct VEC c;
  c.x = mad(a.xx,b.x,mad(a.xy,b.y,a.xz*b.z));
  c.y = mad(a.yx,b.x,mad(a.yy,b.y,a.yz*b.z));
  c.z = mad(a.zx,b.x,mad(a.zy,b.y,a.zz*b.z));
  return c;
}
inline struct VEC multfv(float const a, struct VEC const v) {
  struct VEC r;
  r.x = a*v.x;
  r.y = a*v.y;
  r.z = a*v.z;
  return r;
}
inline struct MAT multfm(float const a, struct MAT const m) {
  struct MAT r;
  r.xx= a*m.xx;  r.xy= a*m.xy;  r.xz= a*m.xz;
  r.yx= a*m.yx;  r.yy= a*m.yy;  r.yz= a*m.yz;
  r.zx= a*m.zx;  r.zy= a*m.zy;  r.zz= a*m.zz;
  return r;
}
inline struct VEC crossvv(struct VEC const a, struct VEC const b) {
  struct VEC r;
  r.x = mad(a.y,b.z,-a.z*b.y);
  r.y = mad(a.z,b.x,-a.x*b.z);
  r.z = mad(a.x,b.y,-a.y*b.x);
  return r;
}
inline struct VEC addvv(struct VEC const a, struct VEC const b) {
  struct VEC r;
  r.x = a.x+b.x;
  r.y = a.y+b.y;
  r.z = a.z+b.z;
  return r;
}
inline struct VEC subvv(struct VEC const a, struct VEC const b) {
  struct VEC r;
  r.x = a.x-b.x;
  r.y = a.y-b.y;
  r.z = a.z-b.z;
  return r;
}
inline float dotvv(struct VEC const a, struct VEC const b) {
  return mad(a.x,b.x,mad(a.y,b.y,a.z*b.z));
}
inline float length2v(struct VEC v) {
  return mad(v.x,v.x,mad(v.y,v.y,v.z*v.z));
}
inline float length2f(float const x, float const y, float const z) {
  return mad(x,x,mad(y,y,z*z));
}
inline float lengthv(struct VEC v) {
  return native_sqrt(mad(v.x,v.x,mad(v.y,v.y,v.z*v.z)));
}
inline float lengthf(float const x, float const y, float const z) {
  return native_sqrt(mad(x,x,mad(y,y,z*z)));
}
inline struct VEC normalizedv(struct VEC v) {
  return multfv( native_divide(1.0,lengthv(v)) , v );
}
inline struct VEC normalizedf(float const x, float const y, float const z) {
  struct VEC r;
  float l = 1.0f / lengthf(x,y,z);
  r.x = x*l;
  r.y = y*l;
  r.z = z*l;
  return r;
}
inline struct VEC proj(struct VEC const a, struct VEC const v) {
  float d = dotvv(a,v) / length2v(a);
  struct VEC r;
  r.x = d*a.x;
  r.y = d*a.y;
  r.z = d*a.z;
  return r;
}
inline struct VEC pproj(struct VEC const a, struct VEC const v) {
  float d = native_divide( dotvv(a,v), length2v(a) );
  struct VEC r;
  r.x = v.x-d*a.x;
  r.y = v.y-d*a.y;
  r.z = v.z-d*a.z;
  return r;
}
inline void rotx(struct VEC *v, float const sin, float const cos) {
  float const tmp = cos * v->y - sin * v->z;
  ;          v->z = sin * v->y + cos * v->z;
  v->y = tmp;
}
inline void roty(struct VEC *v, float const sin, float const cos) {
  float const tmp = cos * v->z - sin * v->x;
  ;          v->x = sin * v->z + cos * v->x;
  v->z = tmp;
}
inline void rotz(struct VEC *v, float const sin, float const cos) {
  float const tmp = cos * v->x - sin * v->y;
  ;          v->y = sin * v->x + cos * v->y;
  v->x = tmp;
}
inline struct MAT projectionMAT(struct VEC const a) {
  struct MAT P;
  float l2 = 1.0f/length2v(a);
  P.xx=l2*a.x*a.x; P.xy=l2*a.x*a.y; P.xz=l2*a.x*a.z;
  P.yx=l2*a.y*a.x; P.yy=l2*a.y*a.y; P.yz=l2*a.y*a.z;
  P.zx=l2*a.z*a.x; P.zy=l2*a.z*a.y; P.zz=l2*a.z*a.z;
  return P;
}
inline struct MAT projectionMAT(float const x, float const y, float const z) {
  struct MAT P;
  float l2 = native_divide(1.0f,length2f(x,y,z));
  P.xx=l2*x*x; P.xy=l2*x*y; P.xz=l2*x*z;
  P.yx=l2*y*x; P.yy=l2*y*y; P.yz=l2*y*z;
  P.zx=l2*z*x; P.zy=l2*z*y; P.zz=l2*z*z;
  return P;
}
inline struct MAT rotationMAT(struct VEC const a, float const t) {
  struct VEC const n = normalizedv(a);
  float sin_t = native_sin(t);
  float cos_t = native_cos(t);
  struct MAT R = multfm(1.0f-cos_t,projectionMAT(n));
  R.xx += cos_t;       R.xy -= sin_t * n.z; R.xz += sin_t * n.y;
  R.yx += sin_t * n.z; R.yy += cos_t;       R.yz -= sin_t * n.x;
  R.zx -= sin_t * n.y; R.zy += sin_t * n.x; R.zz += cos_t;
  return R;
}
inline struct MAT rotationMAT(float const x, float const y, float const z, float const t) {
  struct VEC const n = normalizedf(x,y,z);
  float sin_t = native_sin( t );
  float cos_t = native_cos( t );
  struct MAT R = multfm(1.0f-cos_t,projectionMAT(n));
  R.xx += cos_t;       R.xy -= sin_t * n.z; R.xz += sin_t * n.y;
  R.yx += sin_t * n.z; R.yy += cos_t;       R.yz -= sin_t * n.x;
  R.zx -= sin_t * n.y; R.zy += sin_t * n.x; R.zz += cos_t;
  return R;
}
inline struct MAT transposed(struct MAT const m) {
  return cols(m.xx,m.xy,m.xz,m.yx,m.yy,m.yz,m.zx,m.zy,m.zz);
}

// XFORM
inline struct XFORM multxx(struct XFORM const x2, struct XFORM const x1) {
  //x2.R*( x1.R*v + x1.t )+x2.t
  //x2.R*x1.R*v + x2.R*x1.t + x2.t
  struct XFORM x;
  x.R = multmm(x2.R,x1.R);
  x.t = addvv(multmv(x2.R,x1.t),x2.t);
  return x;
}
inline struct XFORM xrev(struct XFORM const x){
  struct XFORM r;
  r.R = transposed(x.R);
  r.t = multmv(r.R,multfv(-1.0f,x.t));
  return r;
}
inline struct VEC multxv(struct XFORM const x, struct VEC const v) {
  return addvv(multmv(x.R,v),x.t);
}
inline struct XFORM const
vvcxform(struct VEC const _x1, struct VEC const _x2,
         struct VEC const _y1, struct VEC const _y2,
         struct VEC const _c1, struct VEC const _c2)
{
  struct VEC const x1 = normalizedv(_x1);
  struct VEC const y1 = normalizedv(pproj(_x1,_y1));
  struct VEC const z1 = crossvv(x1,y1);
  struct VEC const x2 = normalizedv(_x2);
  struct VEC const y2 = normalizedv(pproj(_x2,_y2));
  struct VEC const z2 = crossvv(x2,y2);
  struct MAT const Rto = cols(x2,y2,z2);
  struct MAT const Rfr = rows(x1,y1,z1);
  struct XFORM x;
  x.R = multmm(Rto,Rfr);
  x.t = subvv(_c2,multmv(x.R,_c1));
  return x;
}
inline struct XFORM const
stub(struct VEC const a, struct VEC const b, struct VEC const c)
{
  struct VEC const x = normalizedv(subvv(a,b));
  struct VEC const z = normalizedv(crossvv(x,subvv(c,b)));
  struct VEC const y = crossvv(z,x);
  struct XFORM s;
  s.R = cols(x,y,z);
  s.t = a;
  return s;
}
inline struct XFORM const
stubrev(struct VEC const a, struct VEC const b, struct VEC const c)
{
  struct VEC const x = normalizedv(subvv(a,b));
  struct VEC const z = normalizedv(crossvv(x,subvv(c,b)));
  struct VEC const y = crossvv(z,x);
  struct XFORM s;
  s.R = rows(x,y,z);
  s.t = multmv(s.R,vec(-a.x,-a.y,-a.z));
  return s;
}
inline struct XFORM const
stubc(struct VEC const cen, struct VEC const a, struct VEC const b, struct VEC const c)
{
  struct VEC const x = normalizedv(subvv(a,b));
  struct VEC const z = normalizedv(crossvv(x,subvv(c,b)));
  struct VEC const y = crossvv(z,x);
  struct XFORM s;
  s.R = cols(x,y,z);
  s.t = cen;
  return s;
}
inline struct XFORM const
stubcrev(struct VEC const cen, struct VEC const a, struct VEC const b, struct VEC const c)
{
  struct VEC const x = normalizedv(subvv(a,b));
  struct VEC const z = normalizedv(crossvv(x,subvv(c,b)));
  struct VEC const y = crossvv(z,x);
  struct XFORM s;
  s.R = rows(x,y,z);
  s.t = multmv(s.R,vec(-cen.x,-cen.y,-cen.z));
  return s;
}
#ifdef __cplusplus
XFORM const stub(core::pose::Pose const & p, Size const rsd) {
  return stub   (p.residue(rsd).xyz("CA"),p.residue(rsd).xyz("N"),p.residue(rsd).xyz("C"));
}
XFORM const stubrev(core::pose::Pose const & p, Size const rsd) {
  return stubrev(p.residue(rsd).xyz("CA"),p.residue(rsd).xyz("N"),p.residue(rsd).xyz("C"));
}
#endif


///////////////////////////////////
inline bool eq(float u, float v) {
  return fabs(u-v) < 0.0001;
}
inline bool eq(struct VEC u, struct VEC v) {
  return ( eq(u.x,v.x) && eq(u.y,v.y) && eq(u.z,v.z) );
}
inline bool eq(struct MAT  u, struct MAT  v) {
  return ( eq(u.xx,v.xx) && eq(u.yx,v.yx) && eq(u.zx,v.zx) &&
           eq(u.xy,v.xy) && eq(u.yy,v.yy) && eq(u.zy,v.zy) &&
           eq(u.xz,v.xz) && eq(u.yz,v.yz) && eq(u.zz,v.zz) );
}
inline bool eq(Vec u, struct VEC  v) { eq(VEC(u),v); }
inline bool eq(struct VEC  u, Vec v) { eq(u,VEC(v)); }
inline bool eq(Mat u, struct MAT  v) { eq(MAT(u),v); }
inline bool eq(struct MAT  u, Mat v) { eq(u,MAT(v)); }
std::ostream & operator<<(std::ostream & out, struct MAT m) {
  out << m.xx << " " << m.xy << " " << m.xz << std::endl;
  out << m.yx << " " << m.yy << " " << m.yz << std::endl;
  out << m.zx << " " << m.zy << " " << m.zz << std::endl;
  return out;
}
std::ostream & operator<<(std::ostream & out, struct VEC v) {
  using std::setw;
  typedef numeric::IOTraits<float> Traits;
  std::ios_base::fmtflags const old_flags = out.flags();
  int const old_precision = out.precision( Traits::precision() );
  out << std::right << std::showpoint << std::uppercase;
  int const w = Traits::width();
  out << setw(w) << v.x << ' ' << setw(w) << v.y << ' ' << setw(w) << v.z;
  out.precision( old_precision );
  out.flags( old_flags );
  return out;
}
void myasserteq(float u, float v, string s) {
  if( !eq(u,v) ) {
    cout << u << std::endl;
    cout << v << std::endl;
    utility_exit_with_message(s);
  }
}
void myasserteq(struct VEC u, struct VEC v, string s) {
  if( !eq(u,v) ) {
    std::cerr << u.xyzVector() << std::endl;
    std::cerr << v.xyzVector() << std::endl;
    utility_exit_with_message(s);
  }
}
inline bool myasserteq(Vec u, struct VEC  v, string s) { myasserteq(VEC(u),v,s); }
inline bool myasserteq(struct VEC  u, Vec v, string s) { myasserteq(u,VEC(v),s); }

void myasserteq(struct MAT u, struct MAT v, string s) {
  if( !eq(u,v) ) {
    std::cerr << u.xyzMatrix() << std::endl;
    std::cerr << v.xyzMatrix() << std::endl;
    utility_exit_with_message(s);
  }
}
inline bool myasserteq(Mat u, struct MAT  v, string s) { myasserteq(MAT(u),v,s); }
inline bool myasserteq(struct MAT  u, Mat v, string s) { myasserteq(u,MAT(v),s); }

void test_MAT_VEC() {
  xyzMatrix<float> m(xyzMatrix<float>::cols(10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f));
  xyzMatrix<float> n(xyzMatrix<float>::cols(10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f));
  xyzVector<float> u(xyzVector<float>(10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f));
  xyzVector<float> v(xyzVector<float>(10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f));
  xyzVector<float> w(xyzVector<float>(10.0f*uniform()-5.0f,10.0f*uniform()-5.0f,10.0f*uniform()-5.0f));
  struct MAT M(m),N(n);
  struct VEC U(u),V(v),W(w);

  myasserteq(  M                    ,       m                 , "conversion" );
  myasserteq(  M.xyzMatrix()        ,       m                 , "conversion" );
  myasserteq(  N                    ,       n                 , "conversion" );
  myasserteq(  U                    ,       u                 , "conversion" );
  myasserteq(  V                    ,       v                 , "conversion" );
  myasserteq(  length2v(V)    ,       dotvv(V,V)          , "dot/len2" );
  myasserteq(  dotvv(U,V)    ,       u.dot(v)          , "dot" );
  myasserteq(  multmm(M,N)            ,       m*n               , "mult MM" );
  //myasserteq(  mult(M,V)            ,       m*v               , "mult MV" );
  myasserteq(   proj(U,U)           ,       projection_matrix(u)*u     , "proj" );
  myasserteq(  pproj(U,V)           ,     projperp(u,v)     , "projperp" );
  myasserteq(  crossvv(U,V)           ,     u.cross(v)     , "cross" );
  myasserteq(  normalizedv(V)        ,     v.normalized()     , "normalize" );
  myasserteq(  rows(U,V,W)       ,     Mat::rows(u,v,w)     , "rows" );
  myasserteq(  cols(U,V,W)       ,     Mat::cols(u,v,w)     , "cols" );
  myasserteq(  rotationMAT(U,123.0),   rotation_matrix(Vec(u),123.0)     , "rotation" );
  myasserteq(  projectionMAT(U),   projection_matrix(Vec(u))     , "projection" );

  //  cout << "PASS!" << std::endl;
}


void test_chi_xform() {
  {
    core::pose::Pose p;
    core::pose::make_pose_from_sequence(p,"N","fa_standard",false);
    remove_lower_terminus_type_from_pose_residue(p,1);
    remove_upper_terminus_type_from_pose_residue(p,1);

    float CHI1 = numeric::conversions::radians(25.0);
    float CHI2 = numeric::conversions::radians(35.0);

    core::pose::Pose tgt(p);
    core::kinematics::Stub s( p.residue(1).xyz("CB"), p.residue(1).xyz("CG"), p.residue(1).xyz("CB"), p.residue(1).xyz("CA") );
    xform_pose_rev(tgt,s);
    Vec  n0 = tgt.residue(1).xyz( "N");
    Vec ca0 = tgt.residue(1).xyz("CA");
    tgt.dump_pdb("out/tgt.pdb");
    rot_pose(tgt, Vec(1,0,0), -35.0);
    tgt.set_chi(2,1,35);
    Vec ca = rotation_matrix_degrees(Vec(1,0,0),-35.0) * ca0;
    Vec  n = rotation_matrix_degrees(Vec(1,0,0),-35.0) *  n0;
    Vec uca = ca.normalized();
    rot_pose(tgt,uca,25.0);
    tgt.set_chi(1,1,25);
    n = rotation_matrix_degrees(uca,25.0) * n0;
    tgt.dump_pdb("out/tgt_chi.pdb");

    core::pose::Pose src(p);
    trans_pose(src,-src.residue(1).xyz("CB"));
    src.set_chi(1,1,25.0);
    src.set_chi(2,1,35.0);
    src.dump_pdb("out/src.pdb");

    VEC CA0(ca0);
    VEC  N0( n0);
    VEC CAl =                              multmv( rotationMAT(1.0f,0.0f,0.0f,-CHI2), CA0 )  ;
    VEC Nl  = multmv( rotationMAT(CAl,CHI1), multmv( rotationMAT(1.0f,0.0f,0.0f,-CHI2),  N0 ) );
    VEC  Ng( src.residue(1).xyz( "N") );
    VEC CAg( src.residue(1).xyz("CA") );

    VEC TX( normalizedv(Nl) );
    VEC TY( normalizedv(pproj(TX,CAl)) );
    VEC TZ (crossvv(TX,TY) );
    VEC SX( normalizedv(Ng) );
    VEC SY( normalizedv(pproj(SX,CAg)) );
    VEC SZ( crossvv(SX,SY) );
    MAT R = multmm(cols(TX,TY,TZ),rows(SX,SY,SZ));

    Vec tx( tgt.residue(1).xyz("N").normalized() );
    Vec ty( projperp(tx,tgt.residue(1).xyz("CA")).normalized() );
    Vec tz( tx.cross(ty));
    Vec sx( src.residue(1).xyz("N").normalized() );
    Vec sy( projperp(sx,src.residue(1).xyz("CA")).normalized() );
    Vec sz( sx.cross(sy));
    Mat tr(Mat::cols(tx,ty,tz));
    Mat sr(Mat::rows(sx,sy,sz));
    Mat r = tr*sr;

    myasserteq( tx, TX, "TX" );
    myasserteq( ty, TY, "TY" );
    myasserteq( tz, TZ, "TZ" );
    myasserteq( sx, SX, "SX" );
    myasserteq( sy, SY, "SY" );
    myasserteq( sz, SZ, "SZ" );
    myasserteq( R,  r , "TR*SR" );

    //rot_pose(src,r);
    rot_pose(src,R.xyzMatrix());

    src.dump_pdb("out/test.pdb");
    utility_exit_with_message("TESTING");

  }
}


// // going from float16 to bb coords
// inline struct VEC  N(float16 const r) {
//   struct VEC v;
//   v.x = r.s0;
//   v.y = r.s1;
//   v.z = r.s2;
// }
// inline struct VEC CA(float16 const r) {
//   struct VEC v;
//   v.x = r.s3;
//   v.y = r.s4;
//   v.z = r.s5;
// }
// inline struct VEC  C(float16 const r) {
//   struct VEC v;
//   v.x = r.s6;
//   v.y = r.s7;
//   v.z = r.s8;
// }
// inline struct VEC  O(float16 const r) {
//   struct VEC v;
//   v.x = r.s9;
//   v.y = r.sa;
//   v.z = r.sb;
// }
// inline struct VEC CB(float16 const r) {
//   struct VEC v;
//   v.x = r.sc;
//   v.y = r.sd;
//   v.z = r.se;
// }

// xform between stubs


// ////  code for packing xform+32bits metadata into float8
// struct TRANS {
//   float xx,xy,xz,yx,yy,yz,zx,zy,zz,x,y,z;
//   uint meta;
// };


// inline float8 packtrans(struct TRANS const t) {
//   float8 f;
//   f.s0 = t.xx;
//   f.s1 = t.yx;
//   f.s2 = t.xy;
//   f.s3 = t.yy;
//   f.s4 = t.x;
//   f.s5 = t.y;
//   f.s6 = t.z;
//   uint tmp = ((t.meta & 4294967293u) | ((t.zx<0.0f)?0u:1u)) | ((t.zy<0.0f)?0u:2u);
//   f.s7 = *((float*)(&tmp));
//   return f;
// }

// inline struct TRANS unpacktrans(float8 const f) {
//   struct TRANS t;
//   t.meta = *((uint*)&f.s7);
//   t.xx = f.s0;
//   t.yx = f.s1;
//   t.xy = f.s2;
//   t.yy = f.s3;
//   t.zx = native_sqrt(-mad(t.xx,t.xx,mad(t.yx,t.yx,-1.0f))) * ((t.meta&1u) ? 1.0f : -1.0f);
//   t.zy = native_sqrt(-mad(t.xy,t.xy,mad(t.yy,t.yy,-1.0f))) * ((t.meta&2u) ? 1.0f : -1.0f);
//   t.xz = mad(t.xy,t.yz,-t.xz*t.yy);
//   t.yz = mad(t.xz,t.yx,-t.xx*t.yz);
//   t.zz = mad(t.xx,t.yy,-t.xy*t.yx);
//   t.x = f.s4;
//   t.y = f.s5;
//   t.z = f.s6;
//   return t;
// }

// inline float8 packmvm(Mat const R, Vec const t, uint const meta) { // first 2 bits of meta cleared
//   float8 f;
//   f.s0 = R.xx();
//   f.s1 = R.yx();
//   f.s2 = R.xy();
//   f.s3 = R.yy();
//   f.s4 = t.x();
//   f.s5 = t.y();
//   f.s6 = t.z();
//   uint tmp = ((meta & 4294967293u) | ((R.zx()<0.0f)?0u:1u)) | ((R.zy()<0.0f)?0u:2u);
//   f.s7 = *((float*)(&tmp));
//   return f;
// }
// inline uint unpackmvm(float8 const f, Mat &R, Vec &t) {
//   uint meta = *((uint*)&f.s7);
//   R.xx() = f.s0;
//   R.yx() = f.s1;
//   R.xy() = f.s2;
//   R.yy() = f.s3;
//   R.zx() = native_sqrt(-mad(R.xx(),R.xx(),mad(R.yx(),R.yx(),-1.0f))) * ((meta&1u) ? 1.0f : -1.0f);
//   R.zy() = native_sqrt(-mad(R.xy(),R.xy(),mad(R.yy(),R.yy(),-1.0f))) * ((meta&2u) ? 1.0f : -1.0f);
//   R.xz() = mad(R.yx(),R.zy(),-R.zx()*R.yy());
//   R.yz() = mad(R.zx(),R.xy(),-R.xx()*R.zy());
//   R.zz() = mad(R.xx(),R.yy(),-R.yx()*R.xy());
//   t.x() = f.s4;
//   t.y() = f.s5;
//   t.z() = f.s6;
//   return meta;
// }

