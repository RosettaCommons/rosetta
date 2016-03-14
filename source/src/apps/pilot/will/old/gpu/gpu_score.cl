bool count_pair(uint i, uint j);
float sqr(float x);

bool count_pair(uint i, uint j) {
  return i != j;
}

inline float sqr(float x) { return x*x; }

__kernel void test_score(
                         __global float4 const * input,
                         __global float        * output
                         ){
  float4 ma = input[get_global_id(0)];
  float tscore = 0.0;

  for(uint i = 0; i < get_global_size(0); i += 4) {
    //    if( i > get_global_id(0) ) continue;
    float4 const aa = input[i+0];
    float4 const ab = input[i+1];
    float4 const ac = input[i+2];
    float4 const ad = input[i+3];
    float const d2a = mad(ma.x-aa.x,ma.x-aa.x,mad( ma.y-aa.y,ma.y-aa.y,sqr(ma.z-aa.z) ));
    float const d2b = mad(ma.x-ab.x,ma.x-ab.x,mad( ma.y-ab.y,ma.y-ab.y,sqr(ma.z-ab.z) ));
    float const d2c = mad(ma.x-ac.x,ma.x-ac.x,mad( ma.y-ac.y,ma.y-ac.y,sqr(ma.z-ac.z) ));
    float const d2d = mad(ma.x-ad.x,ma.x-ad.x,mad( ma.y-ad.y,ma.y-ad.y,sqr(ma.z-ad.z) ));
    //float const d2m = min(min(d2a,d2b),min(d2c,d2d));
    //if( d2m > 100.0 ) continue;
    float const ra = /*d2a <= 100.0 ?*/ 4.0*native_recip(native_sqrt(d2a)); //: 0.0f;
    float const rb = /*d2b <= 100.0 ?*/ 4.0*native_recip(native_sqrt(d2b)); //: 0.0f;
    float const rc = /*d2c <= 100.0 ?*/ 4.0*native_recip(native_sqrt(d2c)); //: 0.0f;
    float const rd = /*d2d <= 100.0 ?*/ 4.0*native_recip(native_sqrt(d2d)); //: 0.0f;
    float const ra2 = ra*ra;
    float const rb2 = rb*rb;
    float const rc2 = rc*rc;
    float const rd2 = rd*rd;
    float const ra3 = ra2*ra;
    float const rb3 = rb2*rb;
    float const rc3 = rc2*rc;
    float const rd3 = rd2*rd;
    float const ra6 = ra2*ra3;
    float const rb6 = rb2*rb3;
    float const rc6 = rc2*rc3;
    float const rd6 = rd2*rd3;
    float const sa = mad(ra6,ra6,-ra6);
    float const sb = mad(rb6,rb6,-rb6);
    float const sc = mad(rc6,rc6,-rc6);
    float const sd = mad(rd6,rd6,-rd6);
    /*                                      */ tscore += count_pair(get_global_id(0),i+0) ? sa : 0.0f;
    /*if( i+1 > get_global_id(0) ) continue;*/ tscore += count_pair(get_global_id(0),i+1) ? sb : 0.0f;
    /*if( i+2 > get_global_id(0) ) continue;*/ tscore += count_pair(get_global_id(0),i+2) ? sc : 0.0f;
    /*if( i+3 > get_global_id(0) ) continue;*/ tscore += count_pair(get_global_id(0),i+3) ? sd : 0.0f;
  }

  __local float lsum[256];
  lsum[get_local_id(0)] = tscore;//get_global_id(0);
  for(uint c=get_local_size(0)/2;c>0;c/=2) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if(c>get_local_id(0)) lsum[get_local_id(0)] += lsum[get_local_id(0)+c];
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  if( get_local_id(0)==0 ) output[get_global_id(0)/get_local_size(0)] = lsum[0];
}

// correct, but slower than cpu
// __kernel void
// minmax1(__global float4 const * input,
//         __global float8       * tmp)
// {
//   float4 const mypt = input[get_global_id(0)];
//   __local float8 lminmax[128];
//   lminmax[ get_local_id(0) ] = (float8){mypt,-mypt};
//   for(uint c=get_local_size(0)/2;c>0;c/=2) {
//     barrier(CLK_LOCAL_MEM_FENCE);
//     if(c>get_local_id(0)) lminmax[get_local_id(0)] = min( lminmax[get_local_id(0)], lminmax[get_local_id(0)+c] );
//   }
//   barrier(CLK_LOCAL_MEM_FENCE);
//   if( get_local_id(0)==0 ) tmp[get_global_id(0)/get_local_size(0)] = lminmax[0];
// }


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


inline short const  short_min( short const a,  short const b) { return (a < b) ? a : b; }
inline short const  short_max( short const a,  short const b) { return (a > b) ? a : b; }
inline short const ushort_min(ushort const a, ushort const b) { return (a < b) ? a : b; }
inline short const ushort_max(ushort const a, ushort const b) { return (a > b) ? a : b; }


ushort nbcount( __constant float4  const * gatom,
                __constant ushort2 const * gstripe,
                float  const gsize,
                uint8  const gdim,
                ushort const ia)
{
  ushort count = 0u;
  ushort const xdim = gdim.s0;
  ushort const ydim = gdim.s1;
  ushort const zdim = gdim.s2;
  float  const gsize2 = gsize*gsize;

  float4 const a = gatom[ia];
  short const ix   = ((a.x < 0.0000001f) ? ((short)0) : ushort_min(xdim,(ushort)(a.x/gsize)));
  short const iy0  = ((a.y < 0.0000001f) ? ((short)0) : a.y/gsize);
  short const iz0  = ((a.z < 0.0000001f) ? ((short)0) : a.z/gsize);
  ushort const iyl = (ushort)short_max(((short)0),((short)iy0)-((short)1));
  ushort const izl = (ushort)short_max(((short)0),((short)iz0)-((short)1));
  ushort const iyu = ushort_min(ydim,(ushort)iy0+((ushort)2));
  ushort const izu = ushort_min(zdim,(ushort)iz0+((ushort)2));
  for(ushort iy = iyl; iy < iyu; ++iy) {
    for(ushort iz = izl; iz < izu; ++iz) {
      ushort const ig = ix+xdim*iy+ydim*xdim*iz;
      ushort const igl = ushort_max(((ushort)gstripe[ig].x),ia+((ushort)1u));
      ushort const igu = ushort_max(((ushort)gstripe[ig].y),ia+((ushort)1u));
      if( igu == ia ) continue;
      ushort i = igl;
      while(i < igu) {
        float4 const a0 = gatom[i+0];
        float4 const a1 = gatom[i+1];
        float4 const a2 = gatom[i+2];
        float4 const a3 = gatom[i+3];
        ushort c0 =                                      (mad(a.x-a0.x,a.x-a0.x,mad(a.y-a0.y,a.y-a0.y,sqr(a.z-a0.z))) <= gsize2);
        ushort c1 = (i+((ushort)1)>=igu) ? ((ushort)0) : (mad(a.x-a1.x,a.x-a1.x,mad(a.y-a1.y,a.y-a1.y,sqr(a.z-a1.z))) <= gsize2);
        ushort c2 = (i+((ushort)2)>=igu) ? ((ushort)0) : (mad(a.x-a2.x,a.x-a2.x,mad(a.y-a2.y,a.y-a2.y,sqr(a.z-a2.z))) <= gsize2);
        ushort c3 = (i+((ushort)3)>=igu) ? ((ushort)0) : (mad(a.x-a3.x,a.x-a3.x,mad(a.y-a3.y,a.y-a3.y,sqr(a.z-a3.z))) <= gsize2);
        count += c0+c1+c2+c3;

        // i+=4;
        /* float const d2a = mad(a.x-a0.x,a.x-a0.x,mad(a.y-a0.y,a.y-a0.y,sqr(a.z-a0.z))); */
        /* float const d2b = mad(a.x-a1.x,a.x-a1.x,mad(a.y-a1.y,a.y-a1.y,sqr(a.z-a1.z))); */
        /* float const d2c = mad(a.x-a2.x,a.x-a2.x,mad(a.y-a2.y,a.y-a2.y,sqr(a.z-a2.z))); */
        /* //float const d2d = mad(a.x-a3.x,a.x-a3.x,mad(a.y-a3.y,a.y-a3.y,sqr(a.z-a3.z))); */
        /* float const ra = (d2a < 4.00000f || d2a > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2a)); */
        /* float const rb = (d2b < 4.00000f || d2b > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2b)); */
        /* float const rc = (d2c < 4.00000f || d2c > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2c)); */
        /* //float const rd = (d2d < 4.00000f || d2d > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2d)); */
        /* float const ra2 = ra*ra; */
        /* float const rb2 = rb*rb; */
        /* float const rc2 = rc*rc; */
        /* //float const rd2 = rd*rd; */
        /* float const ra3 = ra*ra2; */
        /* float const rb3 = rb*rb2; */
        /* float const rc3 = rc*rc2; */
        /* //float const rd3 = rd*rd2; */
        /* float const ra6 = ra2*ra3; */
        /* float const rb6 = rb2*rb3; */
        /* float const rc6 = rc2*rc3; */
        /* //float const rd6 = rd2*rd3; */
        /* float const sa = mad(ra6,ra6,-ra6); */
        /* float const sb = mad(rb6,rb6,-rb6); */
        /* float const sc = mad(rc6,rc6,-rc6); */
        /* //float const sd = mad(rd6,rd6,-rd6); */
        /* gpue += sa+sb+sc;//+sd; */

        i += ((ushort)4);
      }
    }
  }

  return count;
}


float score_test( __constant float4  const * gatom,
                  __constant ushort2 const * gstripe,
                  float  const gsize,
                  uint8  const gdim,
                  ushort const ia)
{
  float gpue = 0.0f;
  ushort const xdim = gdim.s0;
  ushort const ydim = gdim.s1;
  ushort const zdim = gdim.s2;
  float  const gsize2 = gsize*gsize;

  float4 const a = gatom[ia];
  short const ix   = ((a.x < 0.0000001f) ? ((short)0) : ushort_min(xdim,(ushort)(a.x/gsize)));
  short const iy0  = ((a.y < 0.0000001f) ? ((short)0) : a.y/gsize);
  short const iz0  = ((a.z < 0.0000001f) ? ((short)0) : a.z/gsize);
  ushort const iyl = (ushort)short_max(((short)0),((short)iy0)-((short)1));
  ushort const izl = (ushort)short_max(((short)0),((short)iz0)-((short)1));
  ushort const iyu = ushort_min(ydim,(ushort)iy0+((ushort)2));
  ushort const izu = ushort_min(zdim,(ushort)iz0+((ushort)2));
  for(ushort iy = iyl; iy < iyu; ++iy) {
    for(ushort iz = izl; iz < izu; ++iz) {
      ushort const ig = ix+xdim*iy+ydim*xdim*iz;
      ushort const igl = ushort_max(((ushort)gstripe[ig].x),ia+((ushort)1u));
      ushort const igu = ushort_max(((ushort)gstripe[ig].y),ia+((ushort)1u));
      if( igu == ia ) continue;
      ushort i = igl;
      while(i < igu) {
        float4 const a0 = gatom[i+0];
        float4 const a1 = gatom[i+1];
        float4 const a2 = gatom[i+2];
        float4 const a3 = gatom[i+3];
        float4 const a4 = gatom[i+4];
        float const d2a = mad(a.x-a0.x,a.x-a0.x,mad(a.y-a0.y,a.y-a0.y,sqr(a.z-a0.z)));
        float const d2b = mad(a.x-a1.x,a.x-a1.x,mad(a.y-a1.y,a.y-a1.y,sqr(a.z-a1.z)));
        float const d2c = mad(a.x-a2.x,a.x-a2.x,mad(a.y-a2.y,a.y-a2.y,sqr(a.z-a2.z)));
        float const d2d = mad(a.x-a3.x,a.x-a3.x,mad(a.y-a3.y,a.y-a3.y,sqr(a.z-a3.z)));
        float const d2e = mad(a.x-a4.x,a.x-a4.x,mad(a.y-a4.y,a.y-a4.y,sqr(a.z-a4.z)));
        float const ra = (d2a < 4.00000f || d2a > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2a));
        float const rb = (d2b < 4.00000f || d2b > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2b));
        float const rc = (d2c < 4.00000f || d2c > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2c));
        float const rd = (d2d < 4.00000f || d2d > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2d));
        float const re = (d2e < 4.00000f || d2e > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2e));
        float const ra2 = ra*ra;
        float const rb2 = rb*rb;
        float const rc2 = rc*rc;
        float const rd2 = rd*rd;
        float const re2 = re*re;
        float const ra3 = ra*ra2;
        float const rb3 = rb*rb2;
        float const rc3 = rc*rc2;
        float const rd3 = rd*rd2;
        float const re3 = re*re2;
        float const ra6 = ra2*ra3;
        float const rb6 = rb2*rb3;
        float const rc6 = rc2*rc3;
        float const rd6 = rd2*rd3;
        float const re6 = re2*re3;
        float const sa =                              mad(ra6,ra6,-ra6) + native_exp(-d2a);
        float const sb = (i+((ushort)1)>=igu) ? 0.f : mad(rb6,rb6,-rb6) + native_exp(-d2b);
        float const sc = (i+((ushort)2)>=igu) ? 0.f : mad(rc6,rc6,-rc6) + native_exp(-d2c);
        float const sd = (i+((ushort)3)>=igu) ? 0.f : mad(rd6,rd6,-rd6) + native_exp(-d2d);
        float const se = (i+((ushort)4)>=igu) ? 0.f : mad(re6,re6,-re6) + native_exp(-d2e);
        gpue += sa+sb+sc+sd+se;

        i += ((ushort)5);
      }
    }
  }

  return gpue;
}


float score_test_nooctree( __constant float4  const * gatom,
                           __constant ushort2 const * gstripe,
                           float  const gsize,
                           uint8  const gdim,
                           ushort const ia)
{
  float gpue = 0.0f;
  float const gsize2 = gsize*gsize;

  float4 const a = gatom[ia];
  for(ushort i = ((ushort)get_local_id(0)); i < gdim.s4; i += 5) {
    if(i+(ushort)4 <= ia) continue;
    float4 const a0 = gatom[i+0];
    float4 const a1 = gatom[i+1];
    float4 const a2 = gatom[i+2];
    float4 const a3 = gatom[i+3];
    float4 const a4 = gatom[i+4];
    float const d2a = mad(a.x-a0.x,a.x-a0.x,mad(a.y-a0.y,a.y-a0.y,sqr(a.z-a0.z)));
    float const d2b = mad(a.x-a1.x,a.x-a1.x,mad(a.y-a1.y,a.y-a1.y,sqr(a.z-a1.z)));
    float const d2c = mad(a.x-a2.x,a.x-a2.x,mad(a.y-a2.y,a.y-a2.y,sqr(a.z-a2.z)));
    float const d2d = mad(a.x-a3.x,a.x-a3.x,mad(a.y-a3.y,a.y-a3.y,sqr(a.z-a3.z)));
    float const d2e = mad(a.x-a4.x,a.x-a4.x,mad(a.y-a4.y,a.y-a4.y,sqr(a.z-a4.z)));
    float const ra = (d2a < 4.00000f || d2a > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2a));
    float const rb = (d2b < 4.00000f || d2b > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2b));
    float const rc = (d2c < 4.00000f || d2c > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2c));
    float const rd = (d2d < 4.00000f || d2d > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2d));
    float const re = (d2e < 4.00000f || d2e > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2e));
    float const ra2 = ra*ra;
    float const rb2 = rb*rb;
    float const rc2 = rc*rc;
    float const rd2 = rd*rd;
    float const re2 = re*re;
    float const ra3 = ra*ra2;
    float const rb3 = rb*rb2;
    float const rc3 = rc*rc2;
    float const rd3 = rd*rd2;
    float const re3 = re*re2;
    float const ra6 = ra2*ra3;
    float const rb6 = rb2*rb3;
    float const rc6 = rc2*rc3;
    float const rd6 = rd2*rd3;
    float const re6 = re2*re3;
    float const sa = (i+((ushort)1)>=gdim.s4 || ia >= i+((ushort)0) ) ? 0.f : mad(ra6,ra6,-ra6) + native_exp(-d2a);
    float const sb = (i+((ushort)1)>=gdim.s4 || ia >= i+((ushort)1) ) ? 0.f : mad(rb6,rb6,-rb6) + native_exp(-d2b);
    float const sc = (i+((ushort)2)>=gdim.s4 || ia >= i+((ushort)2) ) ? 0.f : mad(rc6,rc6,-rc6) + native_exp(-d2c);
    float const sd = (i+((ushort)3)>=gdim.s4 || ia >= i+((ushort)3) ) ? 0.f : mad(rd6,rd6,-rd6) + native_exp(-d2d);
    float const se = (i+((ushort)4)>=gdim.s4 || ia >= i+((ushort)4) ) ? 0.f : mad(re6,re6,-re6) + native_exp(-d2e);
    gpue += sa+sb+sc+sd+se;
  }

  return gpue;
}


__kernel void
octree( __constant float4  const * gatom,
        __constant ushort2 const * gstripe,
        __constant float   const * gsize_in,
        __constant uint8   const * gdim,
        __global   float         * output )
{
  float gpue = 0.0f;
  //uint count = 0u;
  ushort const xdim = gdim[0].s0;
  ushort const ydim = gdim[0].s1;
  ushort const zdim = gdim[0].s2;
  float const gsize = gsize_in[0];
  float const gsize2 = gsize*gsize;

  for(ushort ia = ((ushort)get_local_id(0)); ia < gdim[0].s4; ia += get_local_size(0)) {
    gpue += score_test(gatom,gstripe,gsize,gdim[0],ia);
  }

  //  if( get_local_id(0)==0 ) output[get_global_id(0)/get_local_size(0)] = gpue;//lsum[0];
  /* //__local uint lsum[256]; */
  //lsum[get_local_id(0)] = count;
  __local float lsum[256];
  lsum[get_local_id(0)] = gpue;
  for(uint c=get_local_size(0)/2;c>0;c/=2) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if(c>get_local_id(0)) lsum[get_local_id(0)] += lsum[get_local_id(0)+c];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if( get_local_id(0)==0 ) output[get_global_id(0)/get_local_size(0)] = lsum[0];

}


