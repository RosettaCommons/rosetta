float
octree_test( xyzStripeHashWithMeta<float,float>::float4  const * gatom,
             xyzStripeHashWithMeta<float,float>::ushort2 const * gstripe,
             float   const * gsize_in,
             xyzStripeHashWithMeta<float,float>::uint8   const * gdim,
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
//          float const d2a = mad(a.x-a0.x,a.x-a0.x,mad(a.y-a0.y,a.y-a0.y, (a.z-a0.z)*(a.z-a0.z) ));!!!!!!!!
          float const d2a = (a.x-a0.x)*(a.x-a0.x) + (a.y-a0.y)*(a.y-a0.y) + (a.z-a0.z)*(a.z-a0.z) ;
//          float const ra = (d2a < 4.00000f || d2a > gsize2) ? 0.0f : 4.0*native_recip(native_sqrt(d2a));
          float const ra = (d2a < 4.00000f || d2a > gsize2) ? 0.0f : 4.0/sqrt(d2a);
          float const ra2 = ra*ra;
          float const ra3 = ra*ra2;
          float const ra6 = ra2*ra3;
          float const sa = ra6*ra6-ra6 + exp(-d2a);
          tmp += sa;
          i = i+((ushort)1);
        }
      }
    }
    e += tmp;
  }
  return e;
}

