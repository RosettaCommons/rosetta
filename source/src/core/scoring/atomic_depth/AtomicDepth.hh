// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/atomic_depth/AtomicDepth.hh
/// @brief  Calculates depth of atoms from the edge of the Sasa surface
/// @author Dong Xu and Yang Zhang
/// @author Brian Coventry (port to Rosetta)

// This class is a port of "EDTSurf: Quick and accurate construction of macromolecular surfaces" by Dong Xu and Yang Zhang
//   https://zhanglab.ccmb.med.umich.edu/EDTSurf/
// References:
//   D. Xu, Y. Zhang (2009) Generating Triangulated Macromolecular Surfaces by Euclidean Distance Transform. PLoS ONE 4(12): e8140.
//   D. Xu, H. Li, Y. Zhang (2013) Protein Depth Calculation and the Use for Improving Accuracy of Protein Fold Recognition.
//         Journal of Computational Biology 20(10):805-816.

// The following is the license originally packaged with EDTSurf
//
/*////////////////////////////////////////////////////////////////
Permission to use, copy, modify, and distribute this program for
any purpose, with or without fee, is hereby granted, provided that
the notices on the head, the reference information, and this
copyright notice appear in all copies or substantial portions of
the Software. It is provided "as is" without express or implied
warranty.
*////////////////////////////////////////////////////////////////

// EDTSurf has been modified here to remove everything besides what is needed for the atomic depth calculation

#ifndef INCLUDED_core_scoring_atomic_depth_AtomicDepth_hh
#define INCLUDED_core_scoring_atomic_depth_AtomicDepth_hh

#include <core/scoring/atomic_depth/AtomicDepth.fwd.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace atomic_depth {


const static char nb[26][3]={{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1},
{1,1,0}, {1,-1,0}, {-1,1,0}, {-1,-1,0}, {1,0,1}, {1,0,-1}, {-1,0,1}, {-1,0,-1}, {0,1,1}, {0,1,-1}, {0,-1,1}, {0,-1,-1},
{1,1,1}, {1,1,-1}, {1,-1,1}, {-1,1,1}, {1,-1,-1}, {-1,-1,1}, {-1,1,-1}, {-1,-1,-1}};

const int fast_one_shell_cut[3] = {6, 12, 9};
const int fast_one_shell_lowb[3] = {0, 6, 18};
const int fast_one_shell_highb[3] = {6, 18, 26};

typedef struct point3d
{
	double x,y,z;
}point3d;
typedef struct volumepixel
{
	float distance;
	bool inout;
	bool isbound;
	bool isdone;
}volumepixel;
typedef struct voxel
{
	int ix,iy,iz;
}voxel;
typedef struct voxel2
{
	short int ix,iy,iz;
}voxel2;

class AtomicDepth : public utility::pointer::ReferenceCount
{
public:
	AtomicDepth( pose::Pose const & pose, Real probe_radius=1.4, bool poly_leu_depth=false, Real resolution=0.25);

	utility::vector1<Real> calcdepth( utility::vector1<conformation::Atom> const & atoms, chemical::AtomTypeSet const & type_set ) const;
	Real calcdepth( conformation::Atom const & atom, chemical::AtomTypeSet const & type_set ) const;

	void visualize_at_depth( Real depth, std::string const & fname, Real fraction=1.0 ) const;

private:
	void boundbox(pose::Pose const & pose,point3d & minp,point3d & maxp);
	void boundingatom( pose::Pose const & pose );
	void initpara( pose::Pose const & pose);
	void fill_vp();
	void fillvoxels( pose::Pose const & pose );
	void fillatom(conformation::Atom const & atom);
	void fastoneshell(int innum,int & allocout,std::vector<std::vector<std::vector<voxel2> > > & boundpoint,int & outnum, int & elimi);
	void fastdistancemap(int type);
	void buildboundary();

private:
	point3d ptran_;
	int boxlength_;
	bool flagradius_;
	double proberadius_;
	double fixsf_;
	double scalefactor_;
	point3d pmin_,pmax_;
	int pheight_,pwidth_,plength_;
	std::vector<int> widxz_;
	std::vector<std::vector<int> > depty_;
	double cutradis_;
	std::vector<std::vector<std::vector<volumepixel > > > vp_;
	std::shared_ptr<std::vector<voxel2> > inarray_;
	std::shared_ptr<std::vector<voxel2> > outarray_;
	int totalsurfacevox_;
	int totalinnervox_;

};




} // namespace atomic_depth
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_atomic_depth_AtomicDepth_fwd_hh
