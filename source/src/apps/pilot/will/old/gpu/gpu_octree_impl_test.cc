  cl_float4 *atoms = new cl_float4[natom];
  uint count = 0;
  for(int ir = 1; ir <= p.n_residue(); ++ir) {
    for(int ia = 1; ia <= p.residue(ir).nheavyatoms(); ++ia) {
      cl_float4 & xyz( atoms[count++] );
      xyz.x = p.residue(ir).atom(ia).xyz().x();
      xyz.y = p.residue(ir).atom(ia).xyz().y();
      xyz.z = p.residue(ir).atom(ia).xyz().z();
      xyz.w = p.residue(ir).atom_type_index(ia);
    }
  }

  float xmn= 9e9,ymn= 9e9,zmn= 9e9;
  float xmx=-9e9,ymx=-9e9,zmx=-9e9;
  for(uint i = 0; i < natom; ++i) {
    xmn = min(xmn,atoms[i].x);
    ymn = min(ymn,atoms[i].y);
    zmn = min(zmn,atoms[i].z);
    xmx = max(xmx,atoms[i].x);
    ymx = max(ymx,atoms[i].y);
    zmx = max(zmx,atoms[i].z);
  }
  // TR << xmx-xmn << " " << ymx-ymn << " " << zmx-zmn << std::endl;
  // utility_exit_with_message("");
  for(uint i = 0; i < natom; ++i) {
    atoms[i].x -= xmn-0.01;
    atoms[i].y -= ymn-0.01;
    atoms[i].z -= zmn-0.01;
  }
  float gsize = GRADIUS;
  uint xdim = ceil((xmx-xmn+0.02)/gsize);
  uint ydim = ceil((ymx-ymn+0.02)/gsize);
  uint zdim = ceil((zmx-zmn+0.02)/gsize);
  cl_ushort2 *grid1 = new cl_ushort2[xdim*ydim*zdim];
  cl_ushort2 *grid3 = new cl_ushort2[xdim*ydim*zdim];
  for(Size i = 0; i < xdim*ydim*zdim; ++i) { grid1[i].x = 0; grid1[i].y = 0; }
  // TR << "atom " << natom << " grid1 " << xdim*ydim*zdim << " " << xdim << " " << ydim << " " << zdim << std::endl;

  for(Size i = 0; i < natom; ++i) {
    uint ix = atoms[i].x/gsize;
    uint iy = atoms[i].y/gsize;
    uint iz = atoms[i].z/gsize;
    ++(grid1[ix+xdim*iy+xdim*ydim*iz].y);
  }
  for(uint i = 1; i < xdim*ydim*zdim; ++i) grid1[i].x = grid1[i-1].x + grid1[i-1].y;
  for(uint i = 1; i < xdim*ydim*zdim; ++i) grid1[i].y = grid1[i  ].x + grid1[i  ].y;
  //compute grid3
  for( int iz = 0; iz < zdim; ++iz) for( int iy = 0; iy < ydim; ++iy) for( int ix = 0; ix < xdim; ++ix) {
        uint const ixl = (uint)max(     0 ,(int)ix-1 );
        uint const ixu =       min(xdim-1u,     ix+1u);
        grid3[ix+xdim*iy+xdim*ydim*iz].x = grid1[ixl+xdim*iy+xdim*ydim*iz].x;
        grid3[ix+xdim*iy+xdim*ydim*iz].y = grid1[ixu+xdim*iy+xdim*ydim*iz].y;
      }

  // for(uint iz = 0; iz < zdim; ++iz) for(uint iy = 0; iy < ydim; ++iy) for(uint ix = 0; ix < xdim; ++ix) {
  //       uint i = ix+xdim*iy+xdim*ydim*iz;
  //       TR<<ix<<" "<<iy<<" "<<iz<<" "<<I(3,grid1[i].x)<<" "<<I(3,grid1[i].y) <<" "<<I(3,grid3[i].x)<<" "<<I(3,grid3[i].y)<<std::endl;
  //     }
  cl_float4 *gatom = new cl_float4[natom];
  cl_ushort *gridc = new cl_ushort[xdim*ydim*zdim];
  for(Size i = 0; i < xdim*ydim*zdim; ++i) gridc[i] = 0;
  for(Size i = 0; i < natom; ++i) {
    uint const ix = atoms[i].x/gsize;
    uint const iy = atoms[i].y/gsize;
    uint const iz = atoms[i].z/gsize;
    uint const ig = ix+xdim*iy+xdim*ydim*iz;
    gatom[ grid1[ig].x + gridc[ig] ] = atoms[i];
    ++(gridc[ig]);
  }

  // for(uint iz = 0; iz < zdim; ++iz) for(uint iy = 0; iy < ydim; ++iy) for(uint ix = 0; ix < xdim; ++ix) {
  //       uint i = ix+xdim*iy+xdim*ydim*iz;
  //       TR << "GRID CELL " << ix << " " << iy << " " << iz << std::endl;
  //       for(Size ig = grid1[i].x; ig < grid1[i].y; ++ig) {
  //         TR << F(7,3,gatom[ig].x) << " " << F(7,3,gatom[ig].y) << " " << F(7,3,gatom[ig].z) << std::endl;
  //       }
  //     }
