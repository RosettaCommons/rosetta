// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/PocketGrid.hh
/// @brief  protocols::pockets::PocketGrid header
/// @author David Johnson
/// @author Ragul Gowthaman

#ifndef INCLUDED_protocols_pockets_PocketGrid_hh
#define INCLUDED_protocols_pockets_PocketGrid_hh

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/pockets/PocketGrid.fwd.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>

// Numeric Headers
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>

#include <utility/vector1_bool.hh>
#include <list>
#include <string>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace pockets {

///@ Cluster of Pocket points
class PCluster
{

	friend class PocketGrid;
	friend class TargetPocketGrid;
	friend class EggshellGrid;
	friend class ComparisonGrid;
	friend class PClusterSet;

public:

  PCluster(core::Size x, core::Size y, core::Size z, core::Real step_);
  PCluster(const PCluster& old);
  ~PCluster() {};
  int size() const {return points_.size();};
  bool testNeighbor(PCluster & c2);
  bool isClose(PCluster const & c2) const;
  bool isTarget(core::Size numTargets = 2) const {return target_ && (numTargets==1 || subtarget_);};
  bool isSolventExposed() const {return solventExposed_;};
  void add (core::Size x, core::Size y, core::Size z);
  typedef struct {
    core::Size x;
    core::Size y;
    core::Size z;
  } Cxyz;
  std::list < Cxyz > points_;

private:
  int count_;
  bool target_, subtarget_, solventExposed_;
  core::Size maxX, minX, maxY, minY, maxZ, minZ;
  core::Real step;

}; //PCluster

///@ Set of cluster of pocket points
class PClusterSet
{

	friend class PocketGrid;
	friend class TargetPocketGrid;
	friend class EggshellGrid;
	friend class ComparisonGrid;

public:
  PClusterSet();
  PClusterSet& operator= (const PClusterSet& old);
  void clear();
  std::list<PCluster>::iterator add (core::Size x, core::Size y, core::Size z, core::Real step);
  void findClusters();
  void join(std::list<PCluster>::iterator c1, std::list<PCluster>::iterator c2);
  core::Real getLargestClusterSize( core::Real const & stepSize, core::Real const & minClusterSize, core::Size const & numTargets, bool ignoreBuried, bool ignoreSurface );
  core::Real getNetClusterSize( core::Real const & stepSize, core::Real const & minClusterSize, core::Size const & numTargets, bool ignoreBuried, bool ignoreSurface );
  core::Size size() { return clusters_.size(); }

private:
  std::list <PCluster> clusters_;

}; //PClusterSet

	///@ Cluster of exemplar points
	class CCluster
	{
		
		friend class PocketGrid;
		friend class TargetPocketGrid;
		friend class EggshellGrid;
		friend class ComparisonGrid;
		friend class CClusterSet;
		
	public:
		
		CCluster(core::Size x, core::Size y, core::Size z, std::string atype, core::Real step_, core::Real absX=0., core::Real absY=0., core::Real absZ=0.);
		CCluster(const CCluster& old);
		~CCluster() {};
		int size() const {return points_.size();};
		bool testNeighbor(CCluster & c2);
		bool isClose(CCluster const & c2) const;
		bool isTarget(core::Size numTargets = 2) const {return target_ && (numTargets==1 || subtarget_);};
		bool isSolventExposed() const {return solventExposed_;};
		void add (core::Size x, core::Size y, core::Size z, std::string atype, core::Real absX=0., core::Real absY=0., core::Real absZ=0.);
		typedef struct {
			core::Size x;
			core::Size y;
			core::Size z;
			core::Real absX;
			core::Real absY;
			core::Real absZ;
			std::string atom_type;
		} Cxyz;
		std::list < Cxyz > points_;
		
	private:
		int count_;
		bool target_, subtarget_, solventExposed_;
		core::Size maxX, minX, maxY, minY, maxZ, minZ;
		core::Real step;
		
	}; //CCluster
	
	///@ Set of clusters of exemplar points
	class CClusterSet
	{
		
		friend class PocketGrid;
		friend class TargetPocketGrid;
		friend class EggshellGrid;
		friend class ComparisonGrid;
		
	public:
		CClusterSet();
		CClusterSet& operator= (const CClusterSet& old);
		void clear();
		std::list<CCluster>::iterator add (core::Size x, core::Size y, core::Size z, std::string aname, core::Real step, core::Real absX=0., core::Real absY=0., core::Real absZ=0.);
		void findClusters();
		void join(std::list<CCluster>::iterator c1, std::list<CCluster>::iterator c2);
		core::Real getLargestClusterSize( core::Real const & stepSize, core::Real const & minClusterSize, core::Size const & numTargets, bool ignoreBuried, bool ignoreSurface );
		core::Real getNetClusterSize( core::Real const & stepSize, core::Real const & minClusterSize, core::Size const & numTargets, bool ignoreBuried, bool ignoreSurface );
		core::Size size() { return clusters_.size(); }
		
	private:
		std::list <CCluster> clusters_;
		
	}; //CClusterSet
	
///@
class PocketGrid : public utility::pointer::ReferenceCount
{

	friend class EggshellGrid;
	friend class TargetPocketGrid;
	friend class ComparisonGrid;

protected:
  enum PtType {EMPTY, PROTEIN, TARGET, SUBTARGET, HSURFACE, PSURFACE, POCKET, PO_SURF, PO_BURIED,T_SURFACE, ST_SURFACE, TP_POCKET, TP_SURF, TP_BURIED, PO_EDGE, TP_EDGE,PSP, EGGSHELL, EGGSHELL_SURROUNDING };
	std::vector < std::vector < std::vector <PtType> > > grid_;
	std::vector < std::vector < std::vector <PtType> > > c_grid_;
	std::vector < std::vector < std::vector <core::Size> > > pockets_;
  //std::vector < conformation::Atom > atoms_;
  //core::Size numAtoms_;
  core::Size xdim_, ydim_, zdim_;
  core::Real xcorn_;
  core::Real ycorn_;
  core::Real zcorn_;
  core::Real stepSize_;
  core::Real maxLen_;
  std::string tag_;
  void init();
  void setup_default_options();
  void newSearch(core::Size thr1, core::Size thr2, core::Size thr3, core::Size max1, core::Size max2, core::Size max3, bool psp=false, bool sps=false);
  bool fill(core::Size x, core::Size y,core::Size z);
  bool touchesSolvent(core::Size x, core::Size y,core::Size z) const;
  bool touchesSS(core::Size x, core::Size y,core::Size z) const;
  bool touchesPS(core::Size x, core::Size y,core::Size z) const;
  bool touchesSurface(core::Size x, core::Size y,core::Size z, bool polar, bool either=false) const;
  bool isSurfacePoint(core::Size x, core::Size y, core::Size z) const;
	bool isTargetPocketPoint(core::Size x, core::Size y, core::Size z) const {return grid_[x][y][z]==TP_EDGE || grid_[x][y][z]==TP_POCKET || grid_[x][y][z]==TP_BURIED || grid_[x][y][z]==TP_SURF;};
	bool isDeepTargetPocketPoint(core::Size x, core::Size y, core::Size z) const {return grid_[x][y][z]==TP_POCKET || grid_[x][y][z]==TP_BURIED;};
	bool isProteinPoint(core::Size x, core::Size y, core::Size z) const {return grid_[x][y][z]==PROTEIN || grid_[x][y][z]==TARGET|| grid_[x][y][z]==SUBTARGET;};
  numeric::xyzVector<core::Real> rotatePoint(core::Real x, core::Real y, core::Real z);

  core::Size pdbno_;
  core::Size expdbno_;
  core::Size numTargets_;
  static const core::Size MAX_TARGETS=2;
	PClusterSet clusters_;
	CClusterSet c_clusters_;

  core::Real size_x_;
  core::Real size_y_;
  core::Real size_z_;
  core::Real spacing_;
  core::Real limit_x_;
  core::Real limit_y_;
  core::Real limit_z_;
  bool restrictSize_;
  bool ignoreBuriedPockets_;
  bool ignoreExposedPockets_;
  core::Real probe_rad_;
  core::Real surf_score_;
  core::Real surf_dist_;
  core::Real bur_score_;
  core::Real bur_dist_;
  bool side_chains_only_;
  bool markpsp_;
  bool marksps_;
  bool exemplarRestriction_;
  bool dumpExemplars_;
  bool search13_;
  core::Real minPockSize_;
  core::Real maxPockSize_;
  numeric::xyzMatrix<core::Real> rot_mat_;
  bool static_grid_;

public:
  PocketGrid();
	virtual ~PocketGrid();
  PocketGrid(const PocketGrid& gr);
  PocketGrid& operator= (const PocketGrid& gr);

  PocketGrid( core::conformation::Residue const & central_rsd );
  PocketGrid( std::vector< core::conformation::ResidueOP > const & central_rsd );
  PocketGrid( core::conformation::Residue const & central_rsd, core::Real x, core::Real y, core::Real z );
  PocketGrid( std::vector< core::conformation::ResidueOP > const & central_rsd, core::Real x, core::Real y, core::Real z );
  PocketGrid( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real y, core::Real z, core::Real const & stepSize );
  PocketGrid( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real y, core::Real z );
  PocketGrid( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real y, core::Real z, core::Real const & stepSize, bool psp, bool sps);

  PocketGrid( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real const & stepSize=1, bool psp=false, bool sps=false);
  PocketGrid( core::Vector const & center, core::Real x, core::Real const & stepSize=1, bool psp=false, bool sps=false);
  PocketGrid( core::Vector const & center, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize=1, bool psp=false, bool sps=false);

  void initialize( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real const size_x, core::Real const size_y, core::Real const size_z, core::Real const & stepSize=1, bool psp=false, bool sps=false);

  void initialize( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real const & stepSize=1, bool psp=false, bool sps=false);
  void initialize( core::Vector const & center, core::Real x, core::Real const & stepSize=1, bool psp=false, bool sps=false);
  void initialize( core::Vector const & center, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize=1, bool psp=false, bool sps=false);
  void initialize( core::conformation::Residue const & central_rsd, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize=1, bool psp=false, bool sps=false);
  void initialize( std::vector< core::conformation::ResidueOP > const & central_rsd, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize=1, bool psp=false, bool sps=false);
  void clear();

  void recenter(core::Real const & xc, core::Real const & yc, core::Real const & zc);
	void recenter( core::conformation::Residue const & central_rsd );
	void recenter( std::vector< core::conformation::ResidueOP > const & central_rsds );
  void recenter( core::Vector const & center);

  void findPockets(core::Size thr, core::Real max);
  void findPSP(core::Size thr, core::Real max);
  void findSPS(core::Size thr, core::Real max);

  void mark (core::Real x, core::Real y, core::Real z, core::Real const & vdWd, core::Real const & buffer, bool polar=false, bool targetResi=false);
  void mark (core::Real x, core::Real y, core::Real z, core::Real const & vdWd, core::Real const & buffer, bool polar=false, int targetResi=0);
  void mark (core::Vector const & center, core::Real const & vdWd, core::Real const & buffer, bool polar=false, bool targetResi=false);
  void mark (core::Vector const & center, core::Real const & vdWd, core::Real const & buffer, bool polar=false, int targetResi=0);

  void clearSmallPockets(core::Size minsize);
	void findClusters();
	void findExemplars(core::pose::Pose const & inPose, Size const total_residues);
	void findClustersByExemplars();
	void linkExemplarsThroughSolvent();

  void dumpGridToFile();
  void dumpGridToFile( std::string const & output_filename );
  void dumpExemplarToFile();
  void dumpExemplarToFile( std::string const & output_filename );
  void dumpTargetPocketsToFile( std::string const & output_filename );
  void dumpTargetPocketsToFile( std::string const & output_filename, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 );
  void dumpTargetPocketsToPDB( std::string const & output_filename, bool minipock=false );
  void dumpTargetPocketsToPDB( std::string const & output_filename, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 );
  void fillTargetPockets();
  void print() const;
  core::Real targetPocketVolume(core::Real const & surf_sc, core::Real const & bur_sc) const ;
  core::Real largestTargetPocketVolume();
  void markPocketDepth(core::Real const & surf_d, core::Real const & bur_d);
  void markEdgeDepth(core::Real const & surf_d, core::Real const & bur_d);
  core::Real netTargetPocketVolume();
  void markDepth(core::Size x, core::Size y, core::Size z, core::Real const & surf_d, core::Real const & bur_d);
  bool markOneEdgeDepth(core::Size x, core::Size y, core::Size z, core::Real const & surf_d, core::Real const & bur_d, bool isTarget);
	void clusterPockets();
	void clusterCPockets();
  core::Real targetPocketSolventSurface() const ;
  core::Real targetPocketProteinSurface() const ;
  core::Real targetPocketHydrophobicProteinSurface() const ;
  core::Real targetPocketPolarProteinSurface() const ;
  core::Real targetPocketHeuristicScore() const ;
  bool isTooSmall() const;
  core::Vector whatIsTooSmall() const;

	bool autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::scoring::constraints::XYZ_Func const & xyz_func, Size const total_residues, bool center_target=true, core::Real x=0.0, core::Real y=0.0, core::Real z=0.0 );

	bool autoexpanding_pocket_eval( std::vector< core::conformation::ResidueOP > const & central_rsd, core::scoring::constraints::XYZ_Func const & xyz_func, Size const total_residues, bool center_target=true, core::Real x=0.0, core::Real y=0.0, core::Real z=0.0 );

	bool autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::pose::Pose const & inPose, bool center_target=true, core::Real x=0.0, core::Real y=0.0, core::Real z=0.0);

	bool autoexpanding_pocket_eval( std::vector< core::conformation::ResidueOP > const & central_rsd, core::pose::Pose const & inPose, bool center_target=true, core::Real x=0.0, core::Real y=0.0, core::Real z=0.0);

	void move_pose_to_standard_orie( core::Size const & central_seqpos, core::pose::Pose & pose );
  static std::vector< core::conformation::ResidueOP > getRelaxResidues( core::pose::Pose const & input_pose, std::string const & resids );
  void randomAngle();
  void zeroAngle();

	core::Real get_pocket_distance( PocketGrid const & template_pocket ) const { return get_pocket_distance( template_pocket, ""); };
	core::Real get_pocket_distance( PocketGrid const & template_pocket, std::string const & comparison_pdbname ) const;

  utility::vector1<core::Real> getBounds();

}; //class PocketGrid


///@
class TargetPocketGrid : public PocketGrid
{

public:

  TargetPocketGrid( const PocketGrid& gr );
  TargetPocketGrid( std::string const & fname );

  void findClusters();
	void dumpTargetPocketsToPDB( std::string const & output_filename, bool minipock = false );
  //void dump_gridfile( std::string const & fname ) const;

	//core::Real get_pocket_distance( TargetPocketGrid const & template_pocket ) const { return get_pocket_distance( template_pocket, ""); };
	//core::Real get_pocket_distance( TargetPocketGrid const & template_pocket, std::string const & comparison_pdbname ) const;

}; //class TargetPocketGrid

///@
class EggshellGrid : public PocketGrid
{

	friend class NonPlaidFingerprint;

private:

	std::list< numeric::xyzVector<core::Real> > eggshell_coord_list_;
	std::list< numeric::xyzVector<core::Real> > extra_coord_list_;
	numeric::xyzVector<core::Real> eggshell_CoM_;

public:

  EggshellGrid( const PocketGrid& gr );
	EggshellGrid( const PocketGrid& grd, core::pose::Pose const & ligand_pose );
  EggshellGrid( const PocketGrid& ext_grd, std::list< numeric::xyzVector<core::Real> > const & eggshell_coord_list);
  EggshellGrid( std::string const & fname );

	std::list< numeric::xyzVector<core::Real> > const & eggshell_coord_list() const { return eggshell_coord_list_; };
	std::list< numeric::xyzVector<core::Real> > const & extra_coord_list() const { return extra_coord_list_; };

  void dump_eggshell( std::string const & fname ) const;
  void dump_eggshell( std::string const & fname, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 ) const;
  void write_eggshell_to_pdb( std::string const & output_pdbname ) const;
  void write_eggshell_to_pdb( std::string const & output_pdbname, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 ) const;
  void write_grid_to_pdb( const PocketGrid& gr, std::string const & output_pdbname ) const;

	numeric::xyzVector<core::Real> calculate_center_xyz (std::list< numeric::xyzVector<core::Real> > const & pocketshell_coord_list );

	core::Real get_eggshell_distance( EggshellGrid const & template_eggshell ) const { return get_eggshell_distance( template_eggshell, ""); };
	core::Real get_eggshell_distance( EggshellGrid const & template_eggshell, std::string const & comparison_pdbname ) const;

}; //class EggshellGrid

///@
class ComparisonGrid {

private:
	
	enum PtType {EMPTY, PROTEIN, TARGET, SUBTARGET, HSURFACE, PSURFACE, POCKET, PO_SURF, PO_BURIED,T_SURFACE, ST_SURFACE, TP_POCKET, TP_SURF, TP_BURIED, PO_EDGE, TP_EDGE,PSP, EGGSHELL, EGGSHELL_SURROUNDING, LIGAND };
	std::vector < std::vector < std::vector <PtType> > > grid_;
	core::Size xdim_, ydim_, zdim_;
	core::Real xcorn_;
	core::Real ycorn_;
	core::Real zcorn_;
	core::Real stepSize_;

public:
	
	ComparisonGrid(const PocketGrid& gr);
	core::Real mark (const protocols::pockets::PocketGrid& gr, core::Real x, core::Real y, core::Real z, core::Real const & vdWd, core::Real const & penalty);
	core::Real compareCoverage(const PocketGrid& gr);

}; //class ComparisonGrid

}//pockets
}//protocols

#endif
