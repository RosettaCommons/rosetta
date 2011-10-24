
struct VEC;
struct MAT;
struct MAT rowsf(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz) ;
struct MAT colsf(float xx, float yx, float zx, float xy, float yy, float zy, float xz, float yz, float zz) ;
struct MAT rows(struct VEC const rx, struct VEC const ry, struct VEC const rz) ;
struct MAT cols(struct VEC const cx, struct VEC const cy, struct VEC const cz) ;
struct MAT multmm(struct MAT const a, struct MAT const b) ;
struct VEC multmv(struct MAT const a, struct VEC const b) ;
struct VEC multfv(float const a, struct VEC const v) ;
struct MAT multfm(float const a, struct MAT const m) ;
struct VEC crossvv(struct VEC const a, struct VEC const b) ;
struct VEC addvv(struct VEC const a, struct VEC const b) ;
struct VEC subvv(struct VEC const a, struct VEC const b) ;
float dotvv(struct VEC const a, struct VEC const b) ;
float length2v(struct VEC v) ;
float length2f(float const x, float const y, float const z) ;
float lengthv(struct VEC v) ;
float lengthf(float const x, float const y, float const z) ;
struct VEC normalizedv(struct VEC v) ;
struct VEC normalizedf(float const x, float const y, float const z) ;
struct VEC proj(struct VEC const a, struct VEC const v) ;
struct VEC pproj(struct VEC const a, struct VEC const v) ;
struct MAT projectionMAT(struct VEC const a) ;
struct MAT projectionMATf(float const x, float const y, float const z) ;
struct MAT rotationMAT(struct VEC const a, float const t) ;
struct MAT rotationMATf(float const x, float const y, float const z, float const t) ;




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
  Vec xyzVector() { return Vec(x,y,z); }
#endif
};
inline struct VEC vec(float x, float y, float z) { struct VEC v; v.x=x; v.y=y; v.z=z; return v; }
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
  Mat xyzMatrix() { return Mat::rows(xx,xy,xz,yx,yy,yz,zx,zy,zz); };
#endif
};
inline struct MAT rowsf(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz) {
  struct MAT m;
  m.xx=xx; m.xy=xy; m.xz=xz;
  m.yx=yx; m.yy=yy; m.yz=yz;
  m.zx=zx; m.zy=zy; m.zz=zz;
  return m;
}
inline struct MAT colsf(float xx, float yx, float zx, float xy, float yy, float zy, float xz, float yz, float zz) {
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
  return multfv( native_recip(lengthv(v)) , v );
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
inline struct MAT projectionMAT(struct VEC const a) {
  struct MAT P;
  float l2 = 1.0f/length2v(a);
  P.xx=l2*a.x*a.x; P.xy=l2*a.x*a.y; P.xz=l2*a.x*a.z;
  P.yx=l2*a.y*a.x; P.yy=l2*a.y*a.y; P.yz=l2*a.y*a.z;
  P.zx=l2*a.z*a.x; P.zy=l2*a.z*a.y; P.zz=l2*a.z*a.z;
  return P;
}
inline struct MAT projectionMATf(float const x, float const y, float const z) {
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
inline struct MAT rotationMATf(float const x, float const y, float const z, float const t) {
  struct VEC const n = normalizedf(x,y,z);
  float sin_t = native_sin( t );
  float cos_t = native_cos( t );
  struct MAT R = multfm(1.0f-cos_t,projectionMAT(n));
  R.xx += cos_t;       R.xy -= sin_t * n.z; R.xz += sin_t * n.y;
  R.yx += sin_t * n.z; R.yy += cos_t;       R.yz -= sin_t * n.x;
  R.zx -= sin_t * n.y; R.zy += sin_t * n.x; R.zz += cos_t;
  return R;
}
