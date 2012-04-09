// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:




#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

struct xyzStripeHashPose;

namespace protocols {
  namespace sic_dock {

    struct SICFast {
      typedef numeric::xyzVector<core::Real> Vec;

      double xmx1,xmn1,ymx1,ymn1,xmx,xmn,ymx,ymn;
      double CTD,CLD,CTD2,CLD2,BIN;
      int xlb,ylb,xub,yub;
      xyzStripeHashPose *xh1_bb_,*xh2_bb_,*xh1_cb_,*xh2_cb_;
      utility::vector1<double> w1_,w2_;
      ObjexxFCL::FArray2D<Vec> ha,hb;



    public:
      SICFast();

      void init(
		core::pose::Pose const & cmp1in,
		utility::vector1<Vec> cmp1cbs, 
		utility::vector1<double> const & cmp1wts, 
		Vec cmp1axs,
		core::pose::Pose const & cmp2in,
		utility::vector1<Vec> cmp2cbs,
		utility::vector1<double> const & cmp2wts,
		Vec cmp2axs
		);

      double slide_into_contact(
				utility::vector1<Vec> pa,
				utility::vector1<Vec> pb,
				utility::vector1<Vec> const & cba, 
				utility::vector1<Vec> const & cbb, 
				Vec ori, 
				double & score,
				double anga, 
				double angb, 
				Vec axsa, 
				Vec axsb,
				bool cmp1or2_a, 
				bool cmp1or2_b
				);


    private:
      void rotate_points(utility::vector1<Vec> & pa, utility::vector1<Vec> & pb, Vec ori);

      void get_bounds(utility::vector1<Vec> & pa, utility::vector1<Vec> & pb);

      double get_score(utility::vector1<Vec> const & cba, utility::vector1<Vec> const & cbb, Vec ori, double mindis,
		       double anga, double angb, Vec axsa, Vec axsb, bool cmp1or2_a, bool cmp1or2_b);

      double refine_mindis_with_xyzHash(utility::vector1<Vec> const & pa, utility::vector1<Vec> const & pb, double mindis, Vec ori,
					double anga, double angb, Vec axsa, Vec axsb, bool cmp1or2_a, bool cmp1or2_b);

      void fill_plane_hash(utility::vector1<Vec> & pa, utility::vector1<Vec> & pb);

      double get_mindis_with_plane_hashes();

    };
    int flood_fill3D(int i, int j, int k, ObjexxFCL::FArray3D<double> & grid, double t) {
      if( grid(i,j,k) <= t ) return 0;
      grid(i,j,k) = t;
      int nmark = 1;
      if(i>1           ) nmark += flood_fill3D(i-1,j  ,k  ,grid,t);
      if(i<grid.size1()) nmark += flood_fill3D(i+1,j  ,k  ,grid,t);
      if(j>1           ) nmark += flood_fill3D(i  ,j-1,k  ,grid,t);
      if(j<grid.size2()) nmark += flood_fill3D(i  ,j+1,k  ,grid,t);
      if(k>1           ) nmark += flood_fill3D(i  ,j  ,k-1,grid,t);
      if(k<grid.size3()) nmark += flood_fill3D(i  ,j  ,k+1,grid,t);
      return nmark;
    }


  } // namespace sic_dock
} // namespace protocols
