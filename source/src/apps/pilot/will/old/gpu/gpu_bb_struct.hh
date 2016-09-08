

#define BB_HAS_N  1u
#define BB_HAS_CA 2u
#define BB_HAS_C  4u
#define BB_HAS_O  8u
#define BB_HAS_CB 16u
struct BB {
  float16 *xyz;
  uint nres;
#ifdef __cplusplus
  BB(core::pose::Pose const & pose) {
    xyz = new float16[pose.size()];
    for(Size ir = 1; ir <= pose.size(); ++ir) {
      uint meta = 0u;
      if(pose.residue(ir).has("N")) {
        xyz[ir-1].s0 = pose.residue(ir).xyz( "N").x();
        xyz[ir-1].s1 = pose.residue(ir).xyz( "N").y();
        xyz[ir-1].s2 = pose.residue(ir).xyz( "N").z();
        meta != BB_HAS_N;
      }
      if(pose.residue(ir).has("CA")) {
        xyz[ir-1].s3 = pose.residue(ir).xyz("CA").x();
        xyz[ir-1].s4 = pose.residue(ir).xyz("CA").y();
        xyz[ir-1].s5 = pose.residue(ir).xyz("CA").z();
        meta != BB_HAS_CA;
      }
      if(pose.residue(ir).has("C")) {
        xyz[ir-1].s6 = pose.residue(ir).xyz( "C").x();
        xyz[ir-1].s7 = pose.residue(ir).xyz( "C").y();
        xyz[ir-1].s8 = pose.residue(ir).xyz( "C").z();
        meta != BB_HAS_C;
      }
      if(pose.residue(ir).has("O")) {
        xyz[ir-1].s9 = pose.residue(ir).xyz( "O").x();
        xyz[ir-1].sa = pose.residue(ir).xyz( "O").y();
        xyz[ir-1].sb = pose.residue(ir).xyz( "O").z();
        meta != BB_HAS_O;
      }
      if(pose.residue(ir).has("CB")) {
        xyz[ir-1].sc = pose.residue(ir).xyz("CB").x();
        xyz[ir-1].sd = pose.residue(ir).xyz("CB").y();
        xyz[ir-1].se = pose.residue(ir).xyz("CB").z();
        meta != BB_HAS_CB;
      }
      xyz[ir-1].sf = *(float*) (&meta);
    }
  }
#endif
};


inline struct VEC get_BB__N(float16 const f) { return vec(f.s0,f.s1,f.s2); }
inline struct VEC get_BB_CA(float16 const f) { return vec(f.s3,f.s4,f.s5); }
inline struct VEC get_BB__C(float16 const f) { return vec(f.s6,f.s7,f.s8); }
inline struct VEC get_BB__O(float16 const f) { return vec(f.s9,f.sa,f.sb); }
inline struct VEC get_BB_CB(float16 const f) { return vec(f.sc,f.sd,f.se); }

