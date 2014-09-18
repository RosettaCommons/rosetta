#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/gpu.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/packing/surf_vol.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/gpu/CL.hh>
#include <apps/pilot/will/gpu/gpu_mat_vec.hh>
#include <apps/pilot/will/xyzStripeHashWithMeta.hh>

static thread_local basic::Tracer TR( "gpu_cbdock" );


struct xyzStripeHashWithMeta {
  float4  const * gatom  ;
  ushort2 const * gstripe;
  float   const   gsize  ;
  uint8   const   gdim   ;
  float   const   thresh2;
};

inline bool clash( xyzStripeHashWithMeta const * h, float4 const a )
{
  ushort const xdim = h->gdim.s0;
  ushort const ydim = h->gdim.s1;
  ushort const zdim = h->gdim.s2;
  float  const gsize2 = h->gsize*h->gsize;
  short  const ix   = ((a.x < 0.0000001f) ? ((short)0) : ushort_min(xdim,(ushort)(a.x/h->gsize)));
  short  const iy0  = ((a.y < 0.0000001f) ? ((short)0) :                          a.y/h->gsize  );
  short  const iz0  = ((a.z < 0.0000001f) ? ((short)0) :                          a.z/h->gsize  );
  ushort const iyl = (ushort)short_max(((short)0),((short)iy0)-((short)1));
  ushort const izl = (ushort)short_max(((short)0),((short)iz0)-((short)1));
  ushort const iyu = ushort_min(ydim,(ushort)iy0+((ushort)2));
  ushort const izu = ushort_min(zdim,(ushort)iz0+((ushort)2));
  for(ushort iy = iyl; iy < iyu; ++iy) {
    for(ushort iz = izl; iz < izu; ++iz) {
      ushort const ig = ix+xdim*iy+ydim*xdim*iz;
      ushort const igl = h->gstripe[ig].x;
      ushort const igu = h->gstripe[ig].y;
      ushort i = igl;
      while(i < igu) {
        float4 const a0 = h->gatom[i+0];
        float4 const a1 = h->gatom[i+1];
        float4 const a2 = h->gatom[i+2];
        float4 const a3 = h->gatom[i+3];
        float const d2a = mad(a.x-a0.x,a.x-a0.x,mad(a.y-a0.y,a.y-a0.y,sqr(a.z-a0.z)));
        float const d2b = mad(a.x-a1.x,a.x-a1.x,mad(a.y-a1.y,a.y-a1.y,sqr(a.z-a1.z)));
        float const d2c = mad(a.x-a2.x,a.x-a2.x,mad(a.y-a2.y,a.y-a2.y,sqr(a.z-a2.z)));
        float const d2d = mad(a.x-a3.x,a.x-a3.x,mad(a.y-a3.y,a.y-a3.y,sqr(a.z-a3.z)));
        bool cond = d2a <= h->thresh2 || d2b <= h->thresh2 || d2c <= h->thresh2 || d2d <= h->thresh2;
        if(cond) return true;
        i += ((ushort)4);
      }
    }
  }
  return false;
}

inline uint contacts( xyzStripeHashWithMeta const * h, float4 const a )
{
  unsigned char ctcnt = (unsigned char)0u;
  ushort const xdim   = h->gdim.s0;
  ushort const ydim   = h->gdim.s1;
  ushort const zdim   = h->gdim.s2;
  float  const gsize2 = h->gsize*h->gsize;
  short  const ix  = ((a.x < 0.0000001f) ? ((short)0) : ushort_min(xdim,(ushort)(a.x/h->gsize)));
  short  const iy0 = ((a.y < 0.0000001f) ? ((short)0) :                          a.y/h->gsize  );
  short  const iz0 = ((a.z < 0.0000001f) ? ((short)0) :                          a.z/h->gsize  );
  ushort const iyl = (ushort)short_max(((short)0),((short)iy0)-((short)1));
  ushort const izl = (ushort)short_max(((short)0),((short)iz0)-((short)1));
  ushort const iyu = ushort_min(ydim,(ushort)iy0+((ushort)2));
  ushort const izu = ushort_min(zdim,(ushort)iz0+((ushort)2));
  for(ushort iy = iyl; iy < iyu; ++iy) {
    for(ushort iz = izl; iz < izu; ++iz) {
      ushort const ig  = ix+xdim*iy+ydim*xdim*iz;
      ushort const igl = h->gstripe[ig].x;
      ushort const igu = h->gstripe[ig].y;
      ushort i = igl;
      while(i < igu) {
        float4 const a0 = h->gatom[i+0];
        float4 const a1 = h->gatom[i+1];
        float4 const a2 = h->gatom[i+2];
        float4 const a3 = h->gatom[i+3];
        float const d2a = mad(a.x-a0.x,a.x-a0.x,mad(a.y-a0.y,a.y-a0.y,sqr(a.z-a0.z)));
        float const d2b = mad(a.x-a1.x,a.x-a1.x,mad(a.y-a1.y,a.y-a1.y,sqr(a.z-a1.z)));
        float const d2c = mad(a.x-a2.x,a.x-a2.x,mad(a.y-a2.y,a.y-a2.y,sqr(a.z-a2.z)));
        float const d2d = mad(a.x-a3.x,a.x-a3.x,mad(a.y-a3.y,a.y-a3.y,sqr(a.z-a3.z)));
        if(d2a <= h->thresh2) ++ctcnt;
        if(d2b <= h->thresh2) ++ctcnt;
        if(d2c <= h->thresh2) ++ctcnt;
        if(d2d <= h->thresh2) ++ctcnt;
        TR << d2a << " " << h->thresh2 << " " << i << std::endl;
        //        ctcnt += (d2a <= h->thresh2) + (d2b <= h->thresh2) + (d2c <= h->thresh2) + (d2d <= h->thresh2);
        i += ((ushort)4);
      }
    }
  }
  return (uint)ctcnt;
}

int main(int argc, char *argv[]) {
  devel::init(argc,argv);
  using namespace basic::options;
  using namespace core::scoring::packing;

  Pose p; core::import_pose::pose_from_pdb(p,option[OptionKeys::in::file::s]()[1]);

  PoseHash clash(3.4,p);
  PoseHash ctact(5.0,p);
  clash.sanity_check();
  // for(int i = 0; i < 20; ++i) {
  //   TR << I(3,i) << " " << F(8,3,clash.grid_atoms()[i].x) << " " << F(8,3,clash.grid_atoms()[i].y) << " " << F(8,3,clash.grid_atoms()[i].z) << std::endl;
  // }
  utility_exit_with_message("otnir");


  // uint8 dimcl; dimcl.s0=clash.xdim(); dimcl.s1=clash.ydim(); dimcl.s2=clash.zdim(); dimcl.s3=dimcl.s0*dimcl.s1*dimcl.s2; dimcl.s4=clash.natom();
  // uint8 dimct; dimct.s0=ctact.xdim(); dimct.s1=ctact.ydim(); dimct.s2=ctact.zdim(); dimct.s3=dimct.s0*dimct.s1*dimct.s2; dimct.s4=ctact.natom();
  // xyzStripeHashWithMeta hcl = { clash.grid_atoms(), clash.grid_stripe(), clash.grid_size(), dimcl, 3.4f*3.4f };
  // xyzStripeHashWithMeta hct = { ctact.grid_atoms(), ctact.grid_stripe(), ctact.grid_size(), dimcl, 5.0f*5.0f };

  // core::id::AtomID_Map<Real> surf = get_surf_vol(p,1.4).surf;

  // uint testi = 100;


  // uint count = 0;
  // for(int i = 0; i < ctact.natom(); ++i) {
  //   float const dx = ctact.grid_atoms()[i].x-ctact.grid_atoms()[testi].x;
  //   float const dy = ctact.grid_atoms()[i].y-ctact.grid_atoms()[testi].y;
  //   float const dz = ctact.grid_atoms()[i].z-ctact.grid_atoms()[testi].z;
  //   if( dx*dx+dy*dy+dz+dz <= 25.0 ) count++;
  // }
  // TR << contacts(&hct,ctact.grid_atoms()[testi]) << " " << count << endl;

  // TR << true + true + true << endl;
  // TR << true + false + true << endl;

  // TR << sizeof(true) << endl;
}
