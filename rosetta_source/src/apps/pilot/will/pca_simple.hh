template<typename T>
inline void swap(T & a, T & b) {
  T tmp = a;
  a = b;
  b = tmp;
}

Size pca_align_BROKEN(core::pose::Pose & p){
  Vec u;
  Size natom = 0;
  for(int ir = 1; ir <= p.n_residue(); ++ir) {
    for(int ia = 1; ia <= p.residue(ir).nheavyatoms(); ++ia) {
      u += p.residue(ir).atom(ia).xyz();
      ++natom;
    }
  }
  u /= natom;
  //trans_pose(p,-u);
  Mat cov(0.0);
  for(int ir = 1; ir <= p.n_residue(); ++ir) {
    for(int ia = 1; ia <= p.residue(ir).nheavyatoms(); ++ia) {
      cov += outer_product( p.xyz(AtomID(ia,ir)), p.xyz(AtomID(ia,ir)) );
    }
  }
  cov /= natom;
  Mat evec;
  Vec eval = eigenvector_jacobi(cov,0.001,evec);
  if(    eval. x()<eval. y()) {
    swap(eval. x(),eval. y());
    swap(evec.xx(),evec.xy());
    swap(evec.yx(),evec.yy());
    swap(evec.zx(),evec.zy());
  }
  if(    eval. x()<eval. z()) {
    swap(eval. x(),eval. z());
    swap(evec.xx(),evec.xz());
    swap(evec.yx(),evec.yz());
    swap(evec.zx(),evec.zz());
  }
  if(    eval. y()<eval. z()) {
    swap(eval. y(),eval. z());
    swap(evec.xy(),evec.xz());
    swap(evec.yy(),evec.yz());
    swap(evec.zy(),evec.zz());
  }
  // cerr<<eval<<std::endl;
  // cerr<<evec<<std::endl;
  //rot_pose(p,evec.transposed());
  return natom;
}
