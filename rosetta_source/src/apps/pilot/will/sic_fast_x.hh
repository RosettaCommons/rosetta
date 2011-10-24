typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

double
sic_fast_x(
  vector1<Vec> & pa , vector1<Vec> & pb ,
  vector1<Vec> & cba, vector1<Vec> & cbb,
  int & cbcount,
  Real const contact_dis = 4.5,
  Real const   clash_dis = 3.5
){
  double const contact_dis2 = contact_dis*contact_dis;
  double const clash_dis2 = clash_dis*clash_dis;
  double const BIN = clash_dis/2.0;

  // get bounds for plane hashes
  double ymx1=-9e9,ymn1=9e9,zmx1=-9e9,zmn1=9e9,ymx=-9e9,ymn=9e9,zmx=-9e9,zmn=9e9;
  for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
    ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
    zmx1 = max(zmx1,ia->z()); zmn1 = min(zmn1,ia->z());
  }
  for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
    ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
    zmx = max(zmx,ib->z()); zmn = min(zmn,ib->z());
  }
  ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);
  zmx = min(zmx,zmx1); zmn = max(zmn,zmn1);

  int ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // one extra on each side for correctness,
  int zlb = (int)floor(zmn/BIN)-2; int zub = (int)ceil(zmx/BIN)+2; // and one extra for outside atoms

  // insert points into hashes
  int const ysize = yub-ylb+1;
  int const zsize = zub-zlb+1;
  ObjexxFCL::FArray2D<Vec> ha(ysize,zsize,Vec(0,0,-9e9)),hb(ysize,zsize,Vec(0,0,9e9));
  for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
    int const iy = (int)ceil(ia->y()/BIN)-ylb;
    int const iz = (int)ceil(ia->z()/BIN)-zlb;
    if( iy < 1 || iy > ysize || iz < 1 || iz > zsize ) continue;
    if( ha(iy,iz).z() < ia->z() ) ha(iy,iz) = *ia;
  }
  for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
    int const iy = (int)ceil(ib->y()/BIN)-ylb;
    int const iz = (int)ceil(ib->z()/BIN)-zlb;
    if( iy < 1 || iy > ysize || iz < 1 || iz > zsize ) continue;
    if( hb(iy,iz).z() > ib->z() ) hb(iy,iz) = *ib;
  }

  // check hashes for min dis
  int imna=0,jmna=0,imnb=0,jmnb=0;
  double mindis = 9e9;
  for(int i = 1; i <= ysize; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
    for(int j = 1; j <= zsize; ++j) {
      for(int k = -2; k <= 2; ++k) {
        if(i+k < 1 || i+k > ysize) continue;
        for(int l = -2; l <= 2; ++l) {
          if(j+l < 1 || j+l > zsize) continue;
          double const ya = ha(i  ,j  ).y();
          double const za = ha(i  ,j  ).z();
          double const yb = hb(i+k,j+l).y();
          double const zb = hb(i+k,j+l).z();
          double const d2 = (ya-yb)*(ya-yb) + (za-zb)*(za-zb);
          if( d2 < clash_dis2 ) {
            double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(BIN*BIN*4.0-d2);
            mindis = min(mindis,dz);
          }
        }
      }
    }
  }

  cbcount = 0;
  for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
    for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
      if( ib->distance_squared( (*ia) + (mindis*Vec(0,0,1)) ) < contact_dis2 ) {
        cbcount++;
      }
    }
  }

 return mindis;
}
