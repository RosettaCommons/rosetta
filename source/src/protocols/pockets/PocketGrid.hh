// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pockets/PocketGrid.hh
/// @brief  protocols::pockets::PocketGrid header
/// @author David Johnson
/// @author Ragul Gowthaman

#ifndef INCLUDED_protocols_pockets_PocketGrid_hh
#define INCLUDED_protocols_pockets_PocketGrid_hh

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/pockets/PocketGrid.fwd.hh>
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>

// Numeric Headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1_bool.hh>
#include <list>
#include <string>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace pockets {

/// @ Cluster of Pocket points
class PCluster
{

	friend class PocketGrid;
	friend class TargetPocketGrid;
	friend class EggshellGrid;
	friend class ElectrostaticpotentialGrid;
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
	struct Cxyz {
		core::Size x;
		core::Size y;
		core::Size z;
#ifdef    SERIALIZATION
		template< class Archive > void save( Archive & arc ) const;
		template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

	};
	std::list < Cxyz > points_;

private:
	int count_;
	bool target_, subtarget_, solventExposed_;
	core::Size maxX, minX, maxY, minY, maxZ, minZ;
	core::Real step;

#ifdef    SERIALIZATION
public:
	PCluster();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //PCluster

/// @ Set of cluster of pocket points
class PClusterSet
{

	friend class PocketGrid;
	friend class TargetPocketGrid;
	friend class EggshellGrid;
	friend class ElectrostaticpotentialGrid;
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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //PClusterSet

/// @ Cluster of exemplar points
class CCluster
{

	friend class PocketGrid;
	friend class TargetPocketGrid;
	friend class EggshellGrid;
	friend class ElectrostaticpotentialGrid;
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
	struct Cxyz {
		core::Size x;
		core::Size y;
		core::Size z;
		core::Real absX;
		core::Real absY;
		core::Real absZ;
		std::string atom_type;

#ifdef    SERIALIZATION
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
	};
	std::list < Cxyz > points_;

private:
	int count_;
	bool target_, subtarget_, solventExposed_;
	core::Size maxX, minX, maxY, minY, maxZ, minZ;
	core::Real step;

#ifdef    SERIALIZATION
public:
	CCluster();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //CCluster

/// @ Set of clusters of exemplar points
class CClusterSet
{

	friend class PocketGrid;
	friend class TargetPocketGrid;
	friend class EggshellGrid;
	friend class ElectrostaticpotentialGrid;
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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //CClusterSet

/// @
class PocketGrid : public utility::pointer::ReferenceCount
{

	friend class EggshellGrid;
	friend class ElectrostaticpotentialGrid;
	friend class TargetPocketGrid;
	friend class ComparisonGrid;
	friend class NonPlaidFingerprint;

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
	core::Real spacing_;
	core::Real stepSize_;
	numeric::xyzVector<core::Real> mid_;
	numeric::xyzVector<core::Size> dim_;
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
	numeric::xyzVector<core::Real> dim() const { return dim_; };
	core::Real spacing() const { return spacing_; };
	PocketGrid( core::conformation::Residue const & central_rsd );
	PocketGrid( std::vector< core::conformation::ResidueCOP > const & central_rsd );
	PocketGrid( core::conformation::Residue const & central_rsd, core::Real x, core::Real y, core::Real z );
	PocketGrid( std::vector< core::conformation::ResidueCOP > const & central_rsd, core::Real x, core::Real y, core::Real z );
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
	void initialize( std::vector< core::conformation::ResidueCOP > const & central_rsd, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize=1, bool psp=false, bool sps=false);
	void clear();

	void recenter(core::Real const & xc, core::Real const & yc, core::Real const & zc);
	void recenter( core::conformation::Residue const & central_rsd );
	void recenter( std::vector< core::conformation::ResidueCOP > const & central_rsds );
	void recenter( core::Vector const & center);

	void findPockets(core::Size thr, core::Real max);
	void findPSP(core::Size thr, core::Real max);
	void findSPS(core::Size thr, core::Real max);

	void mark (core::Real x, core::Real y, core::Real z, core::Real const & vdWd, core::Real const & buffer, bool polar=false, int targetResi=0);
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

	bool autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::scoring::func::XYZ_Func const & xyz_func, Size const total_residues, bool center_target=true, core::Real x=0.0, core::Real y=0.0, core::Real z=0.0 );

	bool autoexpanding_pocket_eval( std::vector< core::conformation::ResidueCOP > const & central_rsd, core::scoring::func::XYZ_Func const & xyz_func, Size const total_residues, bool center_target=true, core::Real x=0.0, core::Real y=0.0, core::Real z=0.0 );

	bool autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::pose::Pose const & inPose, bool center_target=true, core::Real x=0.0, core::Real y=0.0, core::Real z=0.0);

	bool autoexpanding_pocket_eval( std::vector< core::conformation::ResidueCOP > const & central_rsd, core::pose::Pose const & inPose, bool center_target=true, core::Real x=0.0, core::Real y=0.0, core::Real z=0.0);
	void move_pose_to_standard_orie( core::Size const & central_seqpos, core::pose::Pose & pose );
	static std::vector< core::conformation::ResidueCOP > getRelaxResidues( core::pose::Pose const & input_pose, std::string const & resids );
	void randomAngle();
	void zeroAngle();
	core::Real get_pocket_distance( PocketGrid const & template_pocket ) const { return get_pocket_distance( template_pocket, ""); };
	core::Real get_pocket_distance( PocketGrid const & template_pocket, std::string const & comparison_pdbname ) const;
	utility::vector1<core::Real> getBounds();

	void write_pocketGrid_to_pdb( std::string const & output_filename );
	void alter_espGrid( std::string const & espGrid_filename );
	void alter_espGrid_with_bound_ligand( std::string const & espGrid_filename, core::pose::Pose const & protein_pose );
	std::list< numeric::xyzVector<core::Real> > get_connolly_surfacePoints( core::pose::Pose const & protein_pose ) const;
	std::list< numeric::xyzVector<core::Real> > get_connolly_surfacePoints_within_grid( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list);

	virtual bool operator == ( PocketGrid const & other ) const;
	virtual bool same_type_as_me( PocketGrid const & other ) const;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class PocketGrid


/// @
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

	virtual bool same_type_as_me( PocketGrid const & other ) const;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	TargetPocketGrid();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class TargetPocketGrid

/// @
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
	EggshellGrid( core::pose::Pose const & protein_pose, const PocketGrid& pocket_grid, const PocketGrid& grid_for_extshell );
	EggshellGrid( std::string const & fname );
	std::list< numeric::xyzVector<core::Real> > const & eggshell_coord_list() const { return eggshell_coord_list_; };
	std::list< numeric::xyzVector<core::Real> > const & extra_coord_list() const { return extra_coord_list_; };
	void get_connolly_eggshell( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list, const PocketGrid& eggGrid, const PocketGrid& extGrid );
	void get_connolly_eggshell_on_grid( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list, const PocketGrid& eggGrid, const PocketGrid& extGrid );
	void dump_eggshell( std::string const & fname ) const;
	void dump_eggshell( std::string const & fname, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 ) const;
	void write_eggshell_to_pdb( std::string const & output_pdbname ) const;
	void write_eggshell_to_pdb( std::string const & output_pdbname, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 ) const;
	void write_Grid_to_pdb( const PocketGrid& gr, std::string const & output_pdbname ) const;
	numeric::xyzVector<core::Real> calculate_center_xyz (std::list< numeric::xyzVector<core::Real> > const & pocketshell_coord_list );
	core::Real get_eggshell_distance( EggshellGrid const & template_eggshell ) const { return get_eggshell_distance( template_eggshell, ""); };
	core::Real get_eggshell_distance( EggshellGrid const & template_eggshell, std::string const & comparison_pdbname ) const;

	virtual bool operator == ( PocketGrid const & rhs ) const;
	virtual bool same_type_as_me( PocketGrid const & other ) const;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	EggshellGrid();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class EggshellGrid

class ElectrostaticpotentialGrid : public PocketGrid
{

	friend class NonPlaidFingerprint;
	friend class PocketGrid;

private:

	numeric::xyzVector<core::Size> espGrid_dim_;
	numeric::xyzVector<core::Real> espGrid_mid_;
	core::Real espGrid_spacing_;
	std::list< utility::vector1<core::Real> > espGrid_points_list_;

public:
	std::vector < std::vector < std::vector <core::Real> > > espGrid_;
	enum PtType {PROTEIN, SURFACE, SOLVENT, PSEUDO_PROTEIN};
	std::vector < std::vector < std::vector <PtType> > > typGrid_;

	ElectrostaticpotentialGrid(){};//empty default constructor
	void  initialize_espGrid();
	void  initialize_typGrid();
	ElectrostaticpotentialGrid( std::string const & input_espGrid_filename, core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, bool delphi=false );
	void mark_protein_espGrid_points( core::pose::Pose const & protein_pose );
	void mark_atom_espGrid_points( numeric::xyzVector<core::Real> atom_coord, core::Real const & vdWd );
	std::list< utility::vector1<core::Real> > const & espGrid_point_list() const { return espGrid_points_list_; };
	void write_espGrid_to_pdb( std::string const & output_pdbname ) const;
	void write_connollySurface_to_pdb( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list,  std::string const & output_pdbname ) const;
	void get_oegrid_dimensions( std::string const & input_OEgrid_filename );
	void fill_espGrid_values( std::string const & input_OEgrid_filename );
	void fill_espGrid_values_with_type( std::string const & input_OEgrid_filename );
	void print_to_oegrid( std::string const & output_OEgrid_filename );
	void resize_espGrid_to_match_pocketGrid( std::string const & espGrid_filename, PocketGrid const & pg );
	void assign_esp_for_protein_grid_points( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list );
	void assign_esp_for_atom_grid_points( numeric::xyzVector<core::Real> atom_coord, core::Real const & inflated_radius);
	void assign_esp_for_surface_grid_points( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list );
	void assign_esp_for_surface_grid_points_by_nn( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list );
	core::Real get_atom_dist_from_connolly_surface( numeric::xyzVector<core::Real> atom_coord, std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list );
	void mark_buried_solvent_points( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list );
	void get_ZAP_espGrid_values( std::string const & input_espGrid_filename );
	void get_ZAP_espGrid_values_with_type( std::string const & input_espGrid_filename );
	void get_DELPHI_espGrid_values( std::string const & input_espGrid_filename );
	void change_espGrid_for_darc( std::string const & espGrid_filename, core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid );
	void set_surface_esp_to_zero();
	void set_protein_esp_to_zero();
	void cap_espGrid();

	virtual bool operator == ( PocketGrid const & other ) const;
	virtual bool same_type_as_me( PocketGrid const & other ) const;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class ElectrostaticpotentialGrid

/// @
class ComparisonGrid {

private:

	enum PtType {EMPTY, PROTEIN, TARGET, SUBTARGET, HSURFACE, PSURFACE, POCKET, PO_SURF, PO_BURIED,T_SURFACE, ST_SURFACE, TP_POCKET, TP_SURF, TP_BURIED, PO_EDGE, TP_EDGE,PSP, EGGSHELL, EGGSHELL_SURROUNDING, LIGAND };
	std::vector < std::vector < std::vector <PtType> > > grid_;
	core::Size xdim_, ydim_, zdim_;
	//core::Size midx_, midy_, midz_;
	core::Real xcorn_;
	core::Real ycorn_;
	core::Real zcorn_;
	core::Real stepSize_;

public:

	ComparisonGrid(const PocketGrid& gr);
	core::Real mark (const protocols::pockets::PocketGrid& gr, core::Real x, core::Real y, core::Real z, core::Real const & vdWd, core::Real const & penalty);
	core::Real compareCoverage(const PocketGrid& gr);

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ComparisonGrid();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class ComparisonGrid

inline void convert_cartesian_to_grid( numeric::xyzVector<core::Real> const & cart_coord, numeric::xyzVector<core::Real> const & mid, numeric::xyzVector<core::Size> const & dim, core::Real const & spacing, numeric::xyzVector<core::Real> & grid_coord ) {
	grid_coord.x() = (cart_coord.x() - ( mid.x() - ((static_cast<core::Real>(dim.x()-1)/2) * spacing) ))/spacing;
	grid_coord.y() = (cart_coord.y() - ( mid.y() - ((static_cast<core::Real>(dim.y()-1)/2) * spacing) ))/spacing;
	grid_coord.z() = (cart_coord.z() - ( mid.z() - ((static_cast<core::Real>(dim.z()-1)/2) * spacing) ))/spacing;
}

inline void convert_grid_to_cartesian( numeric::xyzVector<core::Real> & grid_coord, numeric::xyzVector<core::Real> const & mid, numeric::xyzVector<core::Size> const & dim, core::Real const & spacing, numeric::xyzVector<core::Real> & cart_coord ) {
	cart_coord.x() = grid_coord.x() * spacing + ( mid.x() - ((static_cast<core::Real>(dim.x()-1)/2) * spacing) );
	cart_coord.y() = grid_coord.y() * spacing + ( mid.y() - ((static_cast<core::Real>(dim.y()-1)/2) * spacing) );
	cart_coord.z() = grid_coord.z() * spacing + ( mid.z() - ((static_cast<core::Real>(dim.z()-1)/2) * spacing) );
}

}//pockets
}//protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_pockets_PocketGrid )
#endif // SERIALIZATION


#endif
