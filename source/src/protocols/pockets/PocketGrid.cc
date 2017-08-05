// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pockets/PocketGrid.cc
/// @brief protocols::pockets::PocketGrid functions
/// @author David Johnson
/// @author Ragul Gowthaman

#include <protocols/pockets/PocketGrid.hh>

// Core Headers
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/Stub.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/conversions.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <core/pack/task/TaskFactory.hh>

// Utility Headers
#include <string>
#include <algorithm>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <cmath>
#include <map>

// Eigen header for SVD support
#include <Eigen/SVD>

#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <numeric/random/random.fwd.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>


namespace numeric {

template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzVector< core::Size > const & xyz
)
{
	arch( xyz.x(), xyz.y(), xyz.z() );
}


template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzVector< core::Size > & xyz
)
{
	arch( xyz.x(), xyz.y(), xyz.z() );
}

}
#endif

namespace protocols {
namespace pockets {

static THREAD_LOCAL basic::Tracer TR("core.grid.Pockets.PocketGrid");

// Visual Studio is more particular about implicit upconversion than gcc
#ifdef _WIN32
core::Real floor(core::Size x) {
	return std::floor((double) x);
}

core::Real pow(core::Size x, core::Size y) {
	return std::pow((double) x, (double) y);
}

core::Real sqrt(core::Size x) {
	return std::sqrt((double) x);
}
#endif

PCluster::PCluster(core::Size x, core::Size y, core::Size z, core::Real step_){
	Cxyz point;
	point.x=x;
	minX=x;
	maxX=x;
	point.y=y;
	minY=y;
	maxY=y;
	point.z=z;
	minZ=z;
	maxZ=z;
	points_.push_back(point);
	target_=false;
	subtarget_=false;
	solventExposed_=false;
	count_=1;
	step=step_;
}

void PCluster::add(core::Size x, core::Size y, core::Size z){
	Cxyz point;
	point.x=x;
	point.y=y;
	point.z=z;
	points_.push_back(point);
}


PCluster::PCluster(const PCluster& old){
	count_ = old.count_;
	points_ = old.points_;
	target_ = old.target_;
	subtarget_ = old.subtarget_;
	solventExposed_ = old.solventExposed_;
	maxX = old.maxX; maxY = old.maxY; maxZ = old.maxZ;
	minX = old.minX; minY = old.minY; minZ = old.minZ;
	step = old.step;
}


bool PCluster::isClose(PCluster const & c2) const{
	if ( minX>c2.maxX+1 ) return false;
	if ( (int)maxX<(int)c2.minX-1 ) return false;
	if ( minY>c2.maxY+1 ) return false;
	if ( (int)maxY<(int)c2.minY-1 ) return false;
	if ( minZ>c2.maxZ+1 ) return false;
	if ( (int)maxZ<(int)c2.minZ-1 ) return false;
	return true;
}

bool PCluster::testNeighbor(PCluster & c2){
	if ( c2.points_.size()>100 &&points_.size()>100 ) {
		auto i = points_.begin();
		int icount=0;
		while ( i!=points_.end() ) {
			auto j = c2.points_.begin();
			int ictend=std::min(500,(int)points_.size()-icount);
			for ( int ict=0; ict<ictend; ++i, ++ict ) {
				int jcount=0;
				while ( j!=c2.points_.end() ) {
					int jctend=std::min(500,(int)c2.points_.size()-jcount);
					for ( int jct=0; jct<jctend; ++j, ++jct ) {
						if ( std::abs((int) i->x - (int) j->x)<=1 ) {
							if ( std::abs((int) i->y - (int) j->y)<=1 ) {
								if ( std::abs((int) i->z - (int) j->z)<=1 ) {
									points_.splice(points_.end(), c2.points_);
									minX=std::min(minX, c2.minX);
									minY=std::min(minY, c2.minY);
									minZ=std::min(minZ, c2.minZ);
									maxX=std::max(maxX, c2.maxX);
									maxY=std::max(maxY, c2.maxY);
									maxZ=std::max(maxZ, c2.maxZ);
									return true;
								}
							}
						}
					}
					jcount+=jctend;
				}
			}
			icount+=ictend;
		}

	} else {

		for ( auto i = points_.begin(); i!=points_.end(); ++i ) {
			for ( auto j = c2.points_.begin(); j!=c2.points_.end(); ++j ) {
				if ( std::abs((int) i->x - (int) j->x)<=1 ) {
					if ( std::abs((int) i->y - (int) j->y)<=1 ) {
						if ( std::abs((int) i->z - (int) j->z)<=1 ) {

							minX=std::min(minX, c2.minX);
							minY=std::min(minY, c2.minY);
							minZ=std::min(minZ, c2.minZ);
							maxX=std::max(maxX, c2.maxX);
							maxY=std::max(maxY, c2.maxY);
							maxZ=std::max(maxZ, c2.maxZ);
							points_.splice(points_.end(), c2.points_);
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}

PClusterSet::PClusterSet(){
	clusters_.clear();
}

PClusterSet& PClusterSet::operator= (const PClusterSet& old){
	if ( this != &old ) {
		clusters_ = old.clusters_;
	}
	return *this;
}

void PClusterSet::clear(){
	clusters_.clear();
}

std::list<PCluster>::iterator PClusterSet::add(core::Size x, core::Size y, core::Size z, core::Real step){
	PCluster tmp(x,y,z, step);
	clusters_.push_front(tmp);
	return clusters_.begin();
}

void PClusterSet::join(std::list<PCluster>::iterator c1, std::list<PCluster>::iterator c2){
	c1->points_.splice(c1->points_.end(), c2->points_);
	c1->target_= c1->target_ || c2->target_;
	c1->subtarget_= c1->subtarget_ || c2->subtarget_;
	c1->solventExposed_= c1->solventExposed_ || c2->solventExposed_;
	clusters_.erase(c2);
}

void PClusterSet::findClusters(){
	for ( auto i = clusters_.begin(); i!=clusters_.end(); ++i ) {
		for ( auto j = i; j!=clusters_.end(); ++j ) {
			if ( i==j ) continue;
			if ( j->isClose(*i) ) {
				if ( j->testNeighbor(*i) ) {
					auto oldI = i;
					--i;
					clusters_.erase(oldI);
					break;
				}
			}
		}
	}
}

core::Real PClusterSet::getLargestClusterSize( core::Real const & stepSize, core::Real const & minClusterSize, core::Size const & numTargets, bool ignoreBuried, bool ignoreSurface){
	core::Real largest_size=0.;
	for ( auto & cluster : clusters_ ) {
		if ( !cluster.isTarget(numTargets) ) continue;
		if ( ignoreBuried && !cluster.isSolventExposed() ) continue;
		if ( ignoreSurface && cluster.isSolventExposed() ) continue;
		if ( cluster.size()>largest_size ) {
			largest_size=cluster.size();
		}
	}
	largest_size = largest_size*pow(stepSize, 3)-minClusterSize;
	if ( largest_size > 0 ) return largest_size;

	return 0;
}

core::Real PClusterSet::getNetClusterSize( core::Real const & stepSize, core::Real const & minClusterSize, core::Size const & numTargets, bool ignoreBuried, bool ignoreSurface ){
	core::Real total_size=0.;
	for ( auto & cluster : clusters_ ) {
		if ( !cluster.isTarget(numTargets) ) continue;
		if ( ignoreBuried && !cluster.isSolventExposed() ) continue;
		if ( ignoreSurface && cluster.isSolventExposed() ) continue;
		if ( cluster.size()*pow(stepSize, 3) >minClusterSize ) {
			total_size+=cluster.size()*pow(stepSize, 3)-minClusterSize;
		}
	}
	if ( total_size>0 ) return total_size;

	return 0;
}


CCluster::CCluster(core::Size x, core::Size y, core::Size z, std::string atype, core::Real step_, core::Real absX, core::Real absY, core::Real absZ){
	using namespace basic::options;
	Cxyz point;
	point.x=x;
	minX=x;
	maxX=x;
	point.y=y;
	minY=y;
	maxY=y;
	point.z=z;
	minZ=z;
	maxZ=z;
	point.atom_type = atype;
	point.absX = absX;
	point.absY = absY;
	point.absZ = absZ;
	points_.push_back(point);
	target_=false;
	subtarget_=false;
	solventExposed_=false;
	count_=1;
	step=step_;
	filterExemplars_=option[ OptionKeys::pocket_grid::pocket_limit_exemplar_color ]();
}

void CCluster::add(core::Size x, core::Size y, core::Size z, std::string atype, core::Real absX, core::Real absY, core::Real absZ){
	Cxyz point;
	point.x=x;
	point.y=y;
	point.z=z;
	point.atom_type = atype;
	point.absX = absX;
	point.absY = absY;
	point.absZ = absZ;
	points_.push_back(point);
}


CCluster::CCluster(const CCluster& old){
	count_ = old.count_;
	points_ = old.points_;
	target_ = old.target_;
	subtarget_ = old.subtarget_;
	solventExposed_ = old.solventExposed_;
	maxX = old.maxX; maxY = old.maxY; maxZ = old.maxZ;
	minX = old.minX; minY = old.minY; minZ = old.minZ;
	step = old.step;
	filterExemplars_=old.filterExemplars_;
}


bool CCluster::isClose(CCluster const & c2) const{
	if ( minX>c2.maxX+7 ) return false;
	if ( (int)maxX<(int)c2.minX-7 ) return false;
	if ( minY>c2.maxY+7 ) return false;
	if ( (int)maxY<(int)c2.minY-7 ) return false;
	if ( minZ>c2.maxZ+7 ) return false;
	if ( (int)maxZ<(int)c2.minZ-7 ) return false;
	return true;
}

bool CCluster::testNeighbor(CCluster & c2){
	for ( auto i = points_.begin(); i!=points_.end(); ++i ) {
		for ( auto j = c2.points_.begin(); j!=c2.points_.end(); ++j ) {

			//Do not allow H-bond atoms to be neighbors with other H-bond atoms
			if ( filterExemplars_ && (i->atom_type.compare("C") != 0 ) &&(j->atom_type.compare("C") != 0 ) ) {
				continue;
			}

			if ( sqrt(pow((double)i->x - (double)j->x,2) + pow((double)i->y - (double)j->y,2) + pow((double)i->z - (double)j->z,2) ) <= 5./step ) {
				if ( ((i->atom_type.compare("C") == 0  || i->atom_type.compare("BBe") == 0 || i->atom_type.compare("BNe") == 0 ) && (j->atom_type.compare("C") == 0  || j->atom_type.compare("BBe") == 0 || j->atom_type.compare("BNe") == 0)) || !filterExemplars_ || (sqrt(pow((double)i->x - (double)j->x,2) + pow((double)i->y - (double)j->y,2) + pow((double)i->z - (double)j->z,2))  <= 3./step ) ) {
					if ( ((i->atom_type.compare("C") != 0) && (j->atom_type.compare("C") != 0)) && filterExemplars_ ) continue;

					minX=std::min(minX, c2.minX);
					minY=std::min(minY, c2.minY);
					minZ=std::min(minZ, c2.minZ);
					maxX=std::max(maxX, c2.maxX);
					maxY=std::max(maxY, c2.maxY);
					maxZ=std::max(maxZ, c2.maxZ);
					if ( c2.target_ ) target_=true;
					if ( c2.subtarget_ ) subtarget_=true;
					if ( c2.solventExposed_ ) solventExposed_=true;
					points_.splice(points_.end(), c2.points_);
					return true;
				}
			}
		}
	}
	return false;
}

CClusterSet::CClusterSet(){
	clusters_.clear();
}

CClusterSet& CClusterSet::operator= (const CClusterSet& old){
	if ( this != &old ) {
		clusters_ = old.clusters_;
	}
	return *this;
}

void CClusterSet::clear(){
	clusters_.clear();
}

std::list<CCluster>::iterator CClusterSet::add(core::Size x, core::Size y, core::Size z, std::string aname, core::Real step, core::Real absX, core::Real absY, core::Real absZ){
	CCluster tmp(x,y,z, aname, step, absX, absY, absZ);
	clusters_.push_front(tmp);
	return clusters_.begin();
}

void CClusterSet::join(std::list<CCluster>::iterator c1, std::list<CCluster>::iterator c2){
	c1->points_.splice(c1->points_.end(), c2->points_);
	clusters_.erase(c2);
}

void CClusterSet::findClusters(){
	for ( auto i = clusters_.begin(); i!=clusters_.end(); ++i ) {
		for ( auto j = i; j!=clusters_.end(); ++j ) {
			if ( i==j ) continue;
			//    if (j->isClose(*i)) {
			if ( j->testNeighbor(*i) ) {
				auto oldI = i;
				--i;
				clusters_.erase(oldI);
				break;
			}
			//    }
		}
	}
}

core::Real CClusterSet::getLargestClusterSize( core::Real const & stepSize, core::Real const & minClusterSize, core::Size const & numTargets, bool ignoreBuried, bool ignoreSurface){
	core::Real largest_size=0.;
	for ( auto & cluster : clusters_ ) {
		if ( !cluster.isTarget(numTargets) ) continue;
		if ( ignoreBuried && !cluster.isSolventExposed() ) continue;
		if ( ignoreSurface && cluster.isSolventExposed() ) continue;
		if ( cluster.size()>largest_size ) {
			largest_size=cluster.size();
		}
	}
	largest_size = largest_size*pow(stepSize, 3)-minClusterSize;
	if ( largest_size > 0 ) return largest_size;

	return 0;
}

core::Real CClusterSet::getNetClusterSize( core::Real const & stepSize, core::Real const & minClusterSize, core::Size const & numTargets, bool ignoreBuried, bool ignoreSurface ){
	core::Real total_size=0.;
	for ( auto & cluster : clusters_ ) {
		if ( !cluster.isTarget(numTargets) ) continue;
		if ( ignoreBuried && !cluster.isSolventExposed() ) continue;
		if ( ignoreSurface && cluster.isSolventExposed() ) continue;
		if ( cluster.size()*pow(stepSize, 3) >minClusterSize ) {
			total_size+=cluster.size()*pow(stepSize, 3)-minClusterSize;
		}
	}
	if ( total_size>0 ) return total_size;

	return 0;
}


void PocketGrid::clear(){
	grid_.clear();
	c_grid_.clear();
	pockets_.clear();
	clusters_.clear();
}

void PocketGrid::init(){
	clear();
	grid_.resize(xdim_);
	c_grid_.resize(xdim_);
	pockets_.resize(xdim_);
	for ( core::Size tx=0; tx<xdim_; tx++ ) {
		grid_[tx].resize(ydim_);
		c_grid_[tx].resize(ydim_);
		pockets_[tx].resize(ydim_);
		for ( core::Size ty=0; ty<ydim_; ty++ ) {
			grid_[tx][ty].resize(zdim_, EMPTY);
			c_grid_[tx][ty].resize(zdim_, EMPTY);
			pockets_[tx][ty].resize(zdim_, EMPTY);
		}
	}
}


void PocketGrid::setup_default_options(){

	using namespace basic::options;
	spacing_=option[ OptionKeys::pocket_grid::pocket_grid_spacing ]();
	if ( option[ OptionKeys::pocket_grid::pocket_grid_size ]()>0 ) {
		size_x_ = option[ OptionKeys::pocket_grid::pocket_grid_size ]();
		size_y_ = option[ OptionKeys::pocket_grid::pocket_grid_size ]();
		size_z_ = option[ OptionKeys::pocket_grid::pocket_grid_size ]();
	} else {
		size_x_ = option[ OptionKeys::pocket_grid::pocket_grid_size_x ]();
		size_y_ = option[ OptionKeys::pocket_grid::pocket_grid_size_y ]();
		size_z_ = option[ OptionKeys::pocket_grid::pocket_grid_size_z ]();
	}
	restrictSize_ = option[ OptionKeys::pocket_grid::pocket_restrict_size ]();
	maxLen_ = option[ OptionKeys::pocket_grid::pocket_max_spacing ]();
	limit_x_ = size_x_*2;
	limit_y_ = size_y_*2;
	limit_z_ = size_z_*2;
	probe_rad_=option[ OptionKeys::pocket_grid::pocket_probe_radius ]();
	side_chains_only_=option[ OptionKeys::pocket_grid::pocket_side ]();
	exemplarRestriction_=option[ OptionKeys::pocket_grid::pocket_filter_by_exemplar ]();
	dumpExemplars_=option[ OptionKeys::pocket_grid::pocket_dump_exemplars ]();
	filterExemplars_=option[ OptionKeys::pocket_grid::pocket_limit_exemplar_color ]();
	limitExemplars_=option[ OptionKeys::pocket_grid::pocket_limit_small_exemplars ]();
	ignoreBuriedPockets_=option[ OptionKeys::pocket_grid::pocket_ignore_buried ]();
	ignoreExposedPockets_=option[ OptionKeys::pocket_grid::pocket_only_buried ]();
	if ( ignoreExposedPockets_ ) ignoreBuriedPockets_=false;
	markpsp_=option[ OptionKeys::pocket_grid::pocket_psp ]();
	marksps_=option[ OptionKeys::pocket_grid::pocket_sps ]();
	search13_=option[ OptionKeys::pocket_grid::pocket_search13 ]();
	surf_score_=option[ OptionKeys::pocket_grid::pocket_surface_score ]();
	surf_dist_=option[ OptionKeys::pocket_grid::pocket_surface_dist ]();
	bur_score_=option[ OptionKeys::pocket_grid::pocket_buried_score ]();
	bur_dist_=option[ OptionKeys::pocket_grid::pocket_buried_dist ]();
	minPockSize_=option[ OptionKeys::pocket_grid::pocket_min_size ]();
	maxPockSize_=option[ OptionKeys::pocket_grid::pocket_max_size ]();
	rot_mat_ = numeric::x_rotation_matrix_degrees((core::Real)0);
	static_grid_=option[ OptionKeys::pocket_grid::pocket_static_grid ]();
}


PocketGrid::PocketGrid(const PocketGrid& gr) :
	utility::pointer::ReferenceCount()
{

	setup_default_options();

	xdim_ = gr.xdim_; ydim_ = gr.ydim_; zdim_ = gr.zdim_;
	xcorn_ = gr.xcorn_; ycorn_ = gr.ycorn_; zcorn_ = gr.zcorn_;
	stepSize_ = gr.stepSize_;
	markpsp_ = gr.markpsp_;
	marksps_ = gr.marksps_;
	exemplarRestriction_ = gr.exemplarRestriction_;
	dumpExemplars_ = gr.dumpExemplars_;
	filterExemplars_ = gr.filterExemplars_;
	limitExemplars_ = gr.limitExemplars_;
	search13_ = gr.search13_;
	tag_=gr.tag_;
	pdbno_=gr.pdbno_;
	expdbno_=gr.expdbno_;
	minPockSize_ = gr.minPockSize_;
	static_grid_ = gr.static_grid_;

	grid_.resize(xdim_);
	c_grid_.resize(xdim_);
	pockets_.resize(xdim_);
	for ( core::Size tx=0; tx<xdim_; tx++ ) {
		grid_[tx].resize(ydim_);
		c_grid_[tx].resize(ydim_);
		pockets_[tx].resize(ydim_);
		for ( core::Size ty=0; ty<ydim_; ty++ ) {
			grid_[tx][ty].resize(zdim_);
			c_grid_[tx][ty].resize(zdim_);
			pockets_[tx][ty].resize(zdim_);
			for ( core::Size tz=0; tz<zdim_; tz++ ) {
				grid_[tx][ty][tz] = gr.grid_[tx][ty][tz];
				c_grid_[tx][ty][tz] = gr.c_grid_[tx][ty][tz];
				pockets_[tx][ty][tz] = gr.pockets_[tx][ty][tz];
			}
		}
	}

}

PocketGrid::~PocketGrid() = default;

PocketGrid& PocketGrid::operator=(const PocketGrid& gr){
	if ( this != &gr ) {

		xdim_ = gr.xdim_; ydim_ = gr.ydim_; zdim_ = gr.zdim_;
		xcorn_ = gr.xcorn_; ycorn_ = gr.ycorn_; zcorn_ = gr.zcorn_;
		stepSize_ = gr.stepSize_;
		markpsp_ = gr.markpsp_;
		marksps_ = gr.marksps_;
		exemplarRestriction_ = gr.exemplarRestriction_;
		dumpExemplars_ = gr.dumpExemplars_;
		filterExemplars_ = gr.filterExemplars_;
		limitExemplars_ = gr.limitExemplars_;
		search13_ = gr.search13_;
		tag_=gr.tag_;
		pdbno_=gr.pdbno_;
		expdbno_=gr.expdbno_;
		minPockSize_ = gr.minPockSize_;

		grid_.resize(xdim_);
		c_grid_.resize(xdim_);
		pockets_.resize(xdim_);
		for ( core::Size tx=0; tx<xdim_; tx++ ) {
			grid_[tx].resize(ydim_);
			c_grid_[tx].resize(ydim_);
			pockets_[tx].resize(ydim_);
			for ( core::Size ty=0; ty<ydim_; ty++ ) {
				grid_[tx][ty].resize(zdim_);
				c_grid_[tx][ty].resize(zdim_);
				pockets_[tx][ty].resize(zdim_);
				for ( core::Size tz=0; tz<zdim_; tz++ ) {
					grid_[tx][ty][tz] = gr.grid_[tx][ty][tz];
					c_grid_[tx][ty][tz] = gr.c_grid_[tx][ty][tz];
					pockets_[tx][ty][tz] = gr.pockets_[tx][ty][tz];
				}
			}
		}
	}
	return *this;
}


PocketGrid::PocketGrid(){
	if ( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
	}
	setup_default_options();
}

PocketGrid::PocketGrid( core::conformation::Residue const & central_rsd ) {
	setup_default_options();
	initialize(central_rsd, size_x_, size_y_, size_z_, spacing_, markpsp_, marksps_);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
	if ( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
	}
}

PocketGrid::PocketGrid( std::vector< core::conformation::ResidueCOP > const & central_rsds ) {
	setup_default_options();
	initialize(central_rsds, size_x_, size_y_, size_z_, spacing_, markpsp_, marksps_);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
	if ( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
	}
}

PocketGrid::PocketGrid( core::conformation::Residue const & central_rsd, core::Real x, core::Real y, core::Real z ) {
	setup_default_options();
	initialize(central_rsd, x, y, z, spacing_, markpsp_, marksps_);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
	if ( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
	}
}

PocketGrid::PocketGrid( std::vector< core::conformation::ResidueCOP > const & central_rsds, core::Real x, core::Real y, core::Real z ) {
	setup_default_options();
	initialize(central_rsds, x, y, z, spacing_, markpsp_, marksps_);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
	if ( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
	}
}

void PocketGrid::randomAngle(){
	core::Real x,y,z;
	x = (int) (numeric::random::uniform() *89 +1);
	y = (int) (numeric::random::uniform() *89 +1);
	z = (int) (numeric::random::uniform() *89 +1);
	rot_mat_ = numeric::x_rotation_matrix_degrees(x);
	rot_mat_ = rot_mat_ * numeric::y_rotation_matrix_degrees(y);
	rot_mat_ = rot_mat_ * numeric::z_rotation_matrix_degrees(z);
}

void PocketGrid::zeroAngle(){
	rot_mat_ = numeric::x_rotation_matrix_degrees((core::Real)0);
}

numeric::xyzVector<core::Real> PocketGrid::rotatePoint(core::Real x, core::Real y, core::Real z){
	numeric::xyzVector<core::Real> coord(x,y,z);
	coord = rot_mat_ * coord;
	return coord;
}


bool PocketGrid::isSurfacePoint(core::Size x, core::Size y, core::Size z) const{
	return (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE || grid_[x][y][z]==ST_SURFACE);
}

void PocketGrid::initialize( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real y, core::Real z, core::Real const & stepSize, bool psp, bool sps){
	markpsp_=psp;
	marksps_=sps;
	stepSize_=stepSize;
	xdim_=2*(core::Size)ceil(x/stepSize_)+1;
	ydim_=2*(core::Size)ceil(y/stepSize_)+1;
	zdim_=2*(core::Size)ceil(z/stepSize_)+1;
	numTargets_=1;
	init();
	recenter(xc,yc,zc);
	pdbno_=0;
	expdbno_=0;
}

PocketGrid::PocketGrid( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real y, core::Real z, core::Real const & stepSize ){
	setup_default_options();
	initialize(xc, yc, zc, x, y, z, stepSize, markpsp_, marksps_);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
}

PocketGrid::PocketGrid( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real y, core::Real z ){
	setup_default_options();
	initialize(xc, yc, zc, x, y, z, spacing_, markpsp_, marksps_);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
}

PocketGrid::PocketGrid( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real y, core::Real z, core::Real const & stepSize, bool psp, bool sps){
	setup_default_options();
	initialize(xc, yc, zc, x, y, z, stepSize, psp, sps);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
}

void PocketGrid::initialize (core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real const & stepSize, bool psp, bool sps){
	markpsp_=psp;
	marksps_=sps;
	stepSize_=stepSize;
	xdim_=2*(core::Size)ceil(x/stepSize_)+1;
	ydim_=2*(core::Size)ceil(x/stepSize_)+1;
	zdim_=2*(core::Size)ceil(x/stepSize_)+1;
	numTargets_=1;
	init();
	recenter(xc,yc,zc);
	pdbno_=0;
	expdbno_=0;
}

PocketGrid::PocketGrid (core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real const & stepSize, bool psp, bool sps){
	setup_default_options();
	initialize(xc, yc, zc, x, stepSize, psp, sps);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
}

void PocketGrid::initialize (core::Vector const & center, core::Real x, core::Real const & stepSize, bool psp, bool sps){
	markpsp_=psp;
	marksps_=sps;
	stepSize_=stepSize;
	core::Real xc=center(1);
	core::Real yc=center(2);
	core::Real zc=center(3);
	xdim_=2*(core::Size)ceil(x/stepSize_)+1;
	ydim_=2*(core::Size)ceil(x/stepSize_)+1;
	zdim_=2*(core::Size)ceil(x/stepSize_)+1;
	numTargets_=1;
	init();
	recenter(xc,yc,zc);
	pdbno_=0;
	expdbno_=0;
}

PocketGrid::PocketGrid (core::Vector const & center, core::Real x, core::Real const & stepSize, bool psp, bool sps){
	setup_default_options();
	initialize(center, x, stepSize, psp, sps);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
}

void PocketGrid::initialize (core::Vector const & center, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize, bool psp, bool sps){
	markpsp_=psp;
	marksps_=sps;
	stepSize_=stepSize;
	xdim_=2*(core::Size)ceil(x/stepSize_)+1;
	ydim_=2*(core::Size)ceil(y/stepSize_)+1;
	zdim_=2*(core::Size)ceil(z/stepSize_)+1;
	numTargets_=1;
	init();
	recenter(center);
	pdbno_=0;
	expdbno_=0;
}

void PocketGrid::initialize (core::conformation::Residue const & central_rsd, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize, bool psp, bool sps){
	markpsp_=psp;
	marksps_=sps;
	stepSize_=stepSize;
	xdim_=2*(core::Size)ceil(x/stepSize_)+1;
	ydim_=2*(core::Size)ceil(y/stepSize_)+1;
	zdim_=2*(core::Size)ceil(z/stepSize_)+1;
	numTargets_=1;
	init();
	recenter(central_rsd);
	pdbno_=0;
	expdbno_=0;
}

void PocketGrid::initialize (std::vector< core::conformation::ResidueCOP > const & central_rsds, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize, bool psp, bool sps){
	markpsp_=psp;
	marksps_=sps;
	stepSize_=stepSize;
	xdim_=2*(core::Size)ceil(x/stepSize_)+1;
	ydim_=2*(core::Size)ceil(y/stepSize_)+1;
	zdim_=2*(core::Size)ceil(z/stepSize_)+1;
	numTargets_=central_rsds.size();
	init();
	recenter(central_rsds);
	pdbno_=0;
	expdbno_=0;
}

PocketGrid::PocketGrid (core::Vector const & center, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize, bool psp, bool sps){
	setup_default_options();
	initialize(center, x, y, z, stepSize, psp, sps);
	using namespace basic::options;
	tag_=option[ OptionKeys::out::output_tag ]();
}

void PocketGrid::recenter(core::Real const & xc, core::Real const & yc, core::Real const & zc){
	init();
	clusters_.clear();
	numeric::xyzVector<core::Real> rpoint;
	using namespace basic::options;
	if ( option[ OptionKeys::pocket_grid::round_pocketGrid_center ]() ) {
		rpoint = rotatePoint(xc, yc, zc);
		xcorn_=(core::Real)std::floor(rpoint.x()-((xdim_-1)/2*stepSize_) + 0.5);
		ycorn_=(core::Real)std::floor(rpoint.y()-((ydim_-1)/2*stepSize_) + 0.5);
		zcorn_=(core::Real)std::floor(rpoint.z()-((zdim_-1)/2*stepSize_) + 0.5);
	} else {
		rpoint = rotatePoint(xc, yc, zc);
		xcorn_=(core::Real)rpoint.x()-((xdim_-1)/2*stepSize_);
		ycorn_=(core::Real)rpoint.y()-((ydim_-1)/2*stepSize_);
		zcorn_=(core::Real)rpoint.z()-((zdim_-1)/2*stepSize_);
	}

	dim_.x() = xdim_;
	dim_.y() = ydim_;
	dim_.z() = zdim_;
	mid_.x() = xcorn_+ ((xdim_-1)/2 * stepSize_);
	mid_.y() = ycorn_+ ((ydim_-1)/2 * stepSize_);
	mid_.z() = zcorn_+ ((zdim_-1)/2 * stepSize_);

}

void PocketGrid::recenter( core::conformation::Residue const & central_rsd ){

	debug_assert( central_rsd.is_protein() );
	core::Vector center;
	center(1)=0.0;
	center(2)=0.0;
	center(3)=0.0;
	core::Size count=0;
	for ( Size i = 1, i_end = central_rsd.nheavyatoms(); i <= i_end; ++i ) {
		//if (central_rsd.atom(i).type()<18||central_rsd.atom(i).type()>21){
		if ( !central_rsd.atom_is_backbone(i) && !central_rsd.atom_is_hydrogen(i) ) {
			center(1)+=central_rsd.atom(i).xyz()(1);
			center(2)+=central_rsd.atom(i).xyz()(2);
			center(3)+=central_rsd.atom(i).xyz()(3);
			count++;
		}
	}
	//Use CA if glycine
	if ( !count ) {
		center(1)=central_rsd.xyz("CA")(1);
		center(2)=central_rsd.xyz("CA")(2);
		center(3)=central_rsd.xyz("CA")(3);
		count++;
	}
	if ( count ) {
		center(1)/=count;
		center(2)/=count;
		center(3)/=count;
	}
	recenter(center(1), center(2), center(3));

}

void PocketGrid::recenter( std::vector< core::conformation::ResidueCOP > const & central_rsds ){
	core::Vector center;
	center(1)=0.0;
	center(2)=0.0;
	center(3)=0.0;
	core::Vector left;
	left(1)=100000.0;
	left(2)=100000.0;
	left(3)=100000.0;
	core::Vector right;
	right(1)=-100000.0;
	right(2)=-100000.0;
	right(3)=-100000.0;
	core::Size count=0;
	int sz = central_rsds.size();
	debug_assert (sz>0);
	//this should restrict the recentering to those atoms that define the grid.
	if ( sz > (int)MAX_TARGETS ) sz = (int)MAX_TARGETS;
	for ( int rnum = 0; rnum < sz; ++rnum ) {
		core::conformation::ResidueCOP central_rsd = central_rsds[rnum];
		core::Vector local_center;
		local_center(1)=0.0;
		local_center(2)=0.0;
		local_center(3)=0.0;
		core::Size local_count=0;

		debug_assert( central_rsd->is_protein() );
		int rcount=0;
		for ( Size i = 1, i_end = central_rsd->nheavyatoms(); i <= i_end; ++i ) {
			//if (central_rsd->atom(i).type()<18||central_rsd->atom(i).type()>21){
			if ( !central_rsd->atom_is_backbone(i) && !central_rsd->atom_is_hydrogen(i) ) {
				local_center(1)+=central_rsd->atom(i).xyz()(1);
				local_center(2)+=central_rsd->atom(i).xyz()(2);
				local_center(3)+=central_rsd->atom(i).xyz()(3);
				local_count++;
				rcount++;
			}
		}
		//Use CA if glycine
		if ( !rcount ) {
			local_center(1)=central_rsd->xyz("CA")(1);
			local_center(2)=central_rsd->xyz("CA")(2);
			local_center(3)=central_rsd->xyz("CA")(3);
			local_count++;
		}

		if ( local_count ) {
			local_center(1)/=local_count;
			local_center(2)/=local_count;
			local_center(3)/=local_count;
			center(1)+=local_center(1);
			center(2)+=local_center(2);
			center(3)+=local_center(3);
			count++;

			if ( left(1)==100000.0 ) {
				left(1)=local_center(1);
				left(2)=local_center(2);
				left(3)=local_center(3);
				right(1)=local_center(1);
				right(2)=local_center(2);
				right(3)=local_center(3);
			} else {
				if ( left(1) > local_center(1) ) left(1) = local_center(1);
				if ( left(2) > local_center(2) ) left(2) = local_center(2);
				if ( left(3) > local_center(3) ) left(3) = local_center(3);
				if ( right(1) < local_center(1) ) right(1) = local_center(1);
				if ( right(2) < local_center(2) ) right(2) = local_center(2);
				if ( right(3) < local_center(3) ) right(3) = local_center(3);
			}
		}
	}
	if ( count ) {
		center(1)/=count;
		center(2)/=count;
		center(3)/=count;
	}
	int buffer=14;
	if ( right(1)-left(1) + buffer > (xdim_-1)*stepSize_ ) {
		core::Size delta = (core::Size)ceil((right(1)-left(1) + buffer - (xdim_-1)*stepSize_)/stepSize_);
		if ( delta % 2 ) delta++;
		TR << "Expanding PocketGrid in the X direction to fit target residues"<<std::endl;
		xdim_ += delta;
	}
	if ( right(2)-left(2) + buffer > (ydim_-1)*stepSize_ ) {
		core::Size delta = (core::Size)ceil((right(2)-left(2) + buffer - (ydim_-1)*stepSize_)/stepSize_);
		if ( delta % 2 ) delta++;
		TR << "Expanding PocketGrid in the Y direction to fit target residues"<<std::endl;
		ydim_ += delta;
	}
	if ( right(3)-left(3) + buffer > (zdim_-1)*stepSize_ ) {
		core::Size delta = (core::Size)ceil((right(3)-left(3) + buffer - (zdim_-1)*stepSize_)/stepSize_);
		if ( delta % 2 ) delta++;
		TR << "Expanding PocketGrid in the Z direction to fit target residues"<<std::endl;
		zdim_ += delta;
	}
	recenter(center(1), center(2), center(3));

}
void PocketGrid::recenter(core::Vector const & center){
	recenter(center(1), center(2), center(3));
}

std::vector< core::conformation::ResidueCOP > PocketGrid::getRelaxResidues( core::pose::Pose const & input_pose, std::string const & resids ) {
	std::vector< core::conformation::ResidueCOP > residues;

	const std::string & delimiters = ",";
	// Skip delimiters at beginning.
	std::string::size_type lastPos = resids.find_first_not_of(delimiters, 0);

	// Find first "non-delimiter".
	std::string::size_type pos = resids.find_first_of(delimiters, lastPos);

	while ( std::string::npos != pos || std::string::npos != lastPos ) {
		std::string const & resid = resids.substr( lastPos, pos - lastPos );
		lastPos = resids.find_first_not_of(delimiters, pos);
		pos = resids.find_first_of(delimiters, lastPos);

		int  central_relax_pdb_number;
		char chain = ' ';
		std::size_t fpos( resid.find(':') );
		if ( fpos != std::string::npos ) {
			central_relax_pdb_number = ObjexxFCL::int_of( resid.substr(0,fpos) );
			if ( fpos != resid.size()-1 ) {
				chain = resid[ fpos+1 ];
			}
		} else {
			central_relax_pdb_number = ObjexxFCL::int_of( resid );
		}
		for ( int j = 1, resnum = input_pose.size(); j <= resnum; ++j ) {
			if ( input_pose.pdb_info()->number(j) == central_relax_pdb_number ) {
				if ( chain != ' ' ) {
					if ( input_pose.pdb_info()->chain(j) == chain ) {
						residues.push_back(input_pose.conformation().residue(j).clone());
						continue;
					}
				} else {
					residues.push_back(input_pose.conformation().residue(j).clone());
					continue;
				}
			}
		}
	}
	if ( residues.size() > MAX_TARGETS ) {
		std::cerr<<"PocketGrid warning: "<<residues.size()<<" target residues specified, but only a maximum of "<<MAX_TARGETS<<" are supported. Ignoring "<<residues.size()-MAX_TARGETS<<" residue(s)."<<std::endl;
	}
	return residues;
}

void PocketGrid::findPockets(core::Size thr, core::Real max){
	thr=(core::Size)std::floor((core::Real)(thr)/stepSize_+.5);
	core::Size thr2=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(3.)+.5);
	core::Size thr3=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(2.)+.5);
	core::Size max1=(core::Size)std::floor(max/stepSize_+.5);
	core::Size max2=(core::Size)std::floor((max/stepSize_)/sqrt(3.)+.5);
	core::Size max3=(core::Size)std::floor((max/stepSize_)/sqrt(2.)+.5);
	newSearch(thr, thr2, thr3, max1, max2, max3);
}

void PocketGrid::findPSP(core::Size thr, core::Real max){
	thr=(core::Size)std::floor((core::Real)(thr)/stepSize_+.5);
	core::Size thr2=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(3.)+.5);
	core::Size thr3=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(2.)+.5);
	core::Size max1=(core::Size)std::floor(max/stepSize_+.5);
	core::Size max2=(core::Size)std::floor((max/stepSize_)/sqrt(3.)+.5);
	core::Size max3=(core::Size)std::floor((max/stepSize_)/sqrt(2.)+.5);
	newSearch(thr, thr2, thr3, max1, max2, max3, true);
}

void PocketGrid::findSPS(core::Size thr, core::Real max){
	thr=(core::Size)std::floor((core::Real)(thr)/stepSize_+.5);
	core::Size thr2=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(3.)+.5);
	core::Size thr3=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(2.)+.5);
	core::Size max1=(core::Size)std::floor(max/stepSize_+.5);
	core::Size max2=(core::Size)std::floor((max/stepSize_)/sqrt(3.)+.5);
	core::Size max3=(core::Size)std::floor((max/stepSize_)/sqrt(2.)+.5);
	newSearch(thr, thr2, thr3, max1, max2, max3, false, true);
}

bool PocketGrid::fill(core::Size x, core::Size y,core::Size z){
	if ( grid_[x][y][z]!=TP_POCKET && grid_[x][y][z]!=POCKET  ) {
		return false;
	}
	if ( pockets_[x][y][z]!=EMPTY ) return false;

	grid_[x][y][z]=TP_POCKET;
	pockets_[x][y][z]=TP_POCKET;
	bool touching=false;
	core::Size x1,y1,z1;
	if ( x==0 ) x1=0;
	else x1=x-1;
	if ( y==0 ) y1=0;
	else y1=y-1;
	if ( z==0 ) z1=0;
	else z1=z-1;
	for ( core::Size x2=x1; x2<std::min(xdim_, x+2); ++x2 ) {
		for ( core::Size y2=y1; y2<std::min(ydim_, y+2); ++y2 ) {
			for ( core::Size z2=z1; z2<std::min(zdim_, z+2); ++z2 ) {
				touching=touching || fill(x2,y2,z2);
			}
		}
	}
	return touching;
}

bool PocketGrid::touchesSolvent(core::Size x, core::Size y,core::Size z) const{
	bool touching=false;
	if ( (x==0) || (y==0)||(z==0) ) return true;
	if ( (x==xdim_-1) || (y==ydim_-1)||(z==zdim_-1) ) return true;
	if ( grid_[x-1][y][z]==EMPTY ) return true;
	if ( grid_[x+1][y][z]==EMPTY ) return true;
	if ( grid_[x][y-1][z]==EMPTY ) return true;
	if ( grid_[x][y+1][z]==EMPTY ) return true;
	if ( grid_[x][y][z-1]==EMPTY ) return true;
	if ( grid_[x][y][z+1]==EMPTY ) return true;

	return touching;
}

bool PocketGrid::touchesSS(core::Size x, core::Size y,core::Size z) const {
	if ( x!=0 ) if ( grid_[x-1][y][z]==EMPTY ) return true;
	if ( y!=0 ) if ( grid_[x][y-1][z]==EMPTY ) return true;
	if ( z!=0 ) if ( grid_[x][y][z-1]==EMPTY ) return true;
	if ( (x!=xdim_-1) ) if ( grid_[x+1][y][z]==EMPTY ) return true;
	if ( (y!=ydim_-1) ) if ( grid_[x][y+1][z]==EMPTY ) return true;
	if ( (z!=zdim_-1) ) if ( grid_[x][y][z+1]==EMPTY ) return true;
	return false;
}

bool PocketGrid::touchesPS(core::Size x, core::Size y,core::Size z) const {
	if ( x!=0 ) {
		if ( isSurfacePoint(x-1,y,z) ) return true;
	}
	if ( y!=0 ) {
		if ( isSurfacePoint(x,y-1,z) ) return true;
	}
	if ( z!=0 ) {
		if ( isSurfacePoint(x,y,z-1) ) return true;
	}
	if ( (x!=xdim_-1) ) {
		if ( isSurfacePoint(x+1,y,z) ) return true;
	}
	if ( (y!=ydim_-1) ) {
		if ( isSurfacePoint(x,y+1,z) ) return true;
	}
	if ( (z!=zdim_-1) ) {
		if ( isSurfacePoint(x,y,z+1) ) return true;
	}
	return false;
}


bool PocketGrid::touchesSurface(core::Size x, core::Size y,core::Size z, bool polar, bool either) const {
	bool touching=false;
	if ( polar ) {
		if ( x!=0 ) if ( grid_[x-1][y][z]==PSURFACE ) return true;
		if ( x!=xdim_-1 ) if ( grid_[x+1][y][z]==PSURFACE ) return true;
		if ( y!=0 ) if ( grid_[x][y-1][z]==PSURFACE ) return true;
		if ( y!=ydim_-1 ) if ( grid_[x][y+1][z]==PSURFACE ) return true;
		if ( z!=0 ) if ( grid_[x][y][z-1]==PSURFACE ) return true;
		if ( z!=zdim_-1 ) if ( grid_[x][y][z+1]==PSURFACE ) return true;
	} else {
		if ( x!=0 ) if ( grid_[x-1][y][z]==HSURFACE ) return true;
		if ( x!=xdim_-1 ) if ( grid_[x+1][y][z]==HSURFACE ) return true;
		if ( y!=0 ) if ( grid_[x][y-1][z]==HSURFACE ) return true;
		if ( y!=ydim_-1 ) if ( grid_[x][y+1][z]==HSURFACE ) return true;
		if ( z!=0 ) if ( grid_[x][y][z-1]==HSURFACE ) return true;
		if ( z!=zdim_-1 ) if ( grid_[x][y][z+1]==HSURFACE ) return true;
	}
	if ( either ) {
		if ( x!=0 ) if ( grid_[x-1][y][z]==PSURFACE ) return true;
		if ( x!=xdim_-1 ) if ( grid_[x+1][y][z]==PSURFACE ) return true;
		if ( y!=0 ) if ( grid_[x][y-1][z]==PSURFACE ) return true;
		if ( y!=ydim_-1 ) if ( grid_[x][y+1][z]==PSURFACE ) return true;
		if ( z!=0 ) if ( grid_[x][y][z-1]==PSURFACE ) return true;
		if ( z!=zdim_-1 ) if ( grid_[x][y][z+1]==PSURFACE ) return true;
		if ( x!=0 ) if ( grid_[x-1][y][z]==HSURFACE ) return true;
		if ( x!=xdim_-1 ) if ( grid_[x+1][y][z]==HSURFACE ) return true;
		if ( y!=0 ) if ( grid_[x][y-1][z]==HSURFACE ) return true;
		if ( y!=ydim_-1 ) if ( grid_[x][y+1][z]==HSURFACE ) return true;
		if ( z!=0 ) if ( grid_[x][y][z-1]==HSURFACE ) return true;
		if ( z!=zdim_-1 ) if ( grid_[x][y][z+1]==HSURFACE ) return true;
	}
	return touching;
}

void PocketGrid::fillTargetPockets(){
	for ( core::Size x=0; x<xdim_; ++x ) {
		for ( core::Size y=0; y<ydim_; ++y ) {
			for ( core::Size z=0; z<zdim_; ++z ) {
				pockets_[x][y][z]=EMPTY;
			}
		}
	}
	for ( core::Size x=0; x<xdim_; ++x ) {
		for ( core::Size y=0; y<ydim_; ++y ) {
			for ( core::Size z=0; z<zdim_; ++z ) {
				if ( grid_[x][y][z]==TP_POCKET ) {
					pockets_[x][y][z]=TP_POCKET;
					core::Size x1,y1,z1;
					if ( x==0 ) x1=0;
					else x1=x-1;
					if ( y==0 ) y1=0;
					else y1=y-1;
					if ( z==0 ) z1=0;
					else z1=z-1;
					for ( core::Size x2=x1; x2<std::min(xdim_, x+2); ++x2 ) {
						for ( core::Size y2=y1; y2<std::min(ydim_, y+2); ++y2 ) {
							for ( core::Size z2=z1; z2<std::min(zdim_, z+2); ++z2 ) {
								fill(x2,y2,z2);
							}
						}
					}
				}
			}
		}
	}
}


void PocketGrid::newSearch(core::Size thr1, core::Size thr2, core::Size thr3, core::Size max1, core::Size max2, core::Size max3, bool psp, bool sps){
	//vector of deltas to add to x,y,z.  Index % 3 = 0 for x deltas, 1 for y, 2 for z
	std::vector<int> deltas;
	deltas.push_back(0);deltas.push_back(0);deltas.push_back(1);
	deltas.push_back(0);deltas.push_back(1);deltas.push_back(0);
	deltas.push_back(1);deltas.push_back(0);deltas.push_back(0);
	deltas.push_back(1);deltas.push_back(1);deltas.push_back(1);
	deltas.push_back(-1);deltas.push_back(1);deltas.push_back(1);
	deltas.push_back(-1);deltas.push_back(-1);deltas.push_back(1);
	deltas.push_back(1);deltas.push_back(-1);deltas.push_back(1);
	if ( search13_ ) {
		deltas.push_back(1);deltas.push_back(0);deltas.push_back(1);
		deltas.push_back(-1);deltas.push_back(0);deltas.push_back(1);
		deltas.push_back(1);deltas.push_back(1);deltas.push_back(0);
		deltas.push_back(-1);deltas.push_back(1);deltas.push_back(0);
		deltas.push_back(0);deltas.push_back(1);deltas.push_back(1);
		deltas.push_back(0);deltas.push_back(-1);deltas.push_back(1);
	}
	int dirs=deltas.size();
	dirs /=3;

	core::Size thr;
	core::Size max;
	//iterate over the entire grid
	for ( core::Size cx=0; cx<xdim_; cx++ ) {
		for ( core::Size cy=0; cy<ydim_; cy++ ) {
			for ( core::Size cz=0; cz<zdim_; cz++ ) {
				//bool tar_surf=false;

				//if a point is a surface, or pocket and doind psp searching, or solvent and doing a sps search, search in all directions from that point.
				if ( (!sps&&(isSurfacePoint(cx, cy, cz)||(grid_[cx][cy][cz]==POCKET&&psp)||(grid_[cx][cy][cz]==TP_POCKET&&psp)||(grid_[cx][cy][cz]==POCKET&&sps)))||(grid_[cx][cy][cz]==EMPTY&&sps) ) {
					bool tar_surf=false;
					if ( grid_[cx][cy][cz]==T_SURFACE||grid_[cx][cy][cz]==ST_SURFACE||(grid_[cx][cy][cz]==TP_POCKET&&psp) ) tar_surf=true;
					//iterate over the searching directions, get the index associated with x
					for ( int i=0; i<dirs*3; i+=3 ) {
						if ( i<=6 ) {
							thr=thr1;
							max=max1;
						} else if ( i<=18 ) {
							thr=thr2;
							max=max2;
						} else {
							thr=thr3;
							max=max3;
						}
						bool t_surf=tar_surf;
						int x;
						int y;
						int z;
						int count;
						bool marked=true;
						//starting at that point, search outward in the direction being searched until max distance
						for ( x=(int)cx+deltas[i],y=(int)cy+deltas[i+1],z=(int)cz+deltas[i+2], count=0; count<(int)max; x+=deltas[i],y+=deltas[i+1],z+=deltas[i+2], count++ ) {
							if ( x<0 || x>=(int)xdim_ ) break;
							if ( y<0 || y>=(int)ydim_ ) break;
							if ( z<0 || z>=(int)zdim_ ) break;

							//Just keep going for EMPTY (solvent)
							if ( grid_[x][y][z]==EMPTY&&!sps ) {
								marked=false;
								continue;
							} else if ( (grid_[x][y][z]==POCKET&&!psp) || (grid_[x][y][z]==TP_POCKET&&!psp) ) {
								//if not doing psp, treat pockets as EMPTY coz they are still solvent
								if ( sps ) {
									marked=false;
								}
								continue;
							} else if ( sps && isSurfacePoint(x,y,z) ) break;
							//if it's a surface, we want to go back to the starting point and fill the solvents in as pockets
							else if ( (!sps && isSurfacePoint(x,y,z)) ||
									(grid_[x][y][z]==POCKET&&psp) ||
									(grid_[x][y][z]==TP_POCKET&&psp) ||
									(grid_[x][y][z]==EMPTY&&sps) ) {

								//if there is no EMPTY, then everything has been marked already.  No need to go further.
								if ( marked ) break;


								if ( grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==ST_SURFACE||(grid_[x][y][z]==TP_POCKET&&psp) ) t_surf=true;
								if ( (count>(int)thr && !psp)||psp||sps ) {
									for ( int c=count; c>0; c-- ) {
										if ( isSurfacePoint(x-c*deltas[i],y-c*deltas[i+1],z-c*deltas[i+2]) ) TR.Fatal<<"MAJOR ERROR, overwriting surface with pocket\n";
										if ( sps ) {
											grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]=EMPTY;
										} else if ( markpsp_||!psp ) {
											if ( t_surf ) {
												grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]=TP_POCKET;
											} else if ( grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]!=TP_POCKET )  grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]=POCKET;
										} else {
											if ( grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]==EMPTY ) grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]=PSP;
										}
									}
								}
								break;
							} else break;
							//if it's not a surface, pocket, or solvent then can't be a pocket.  Stop searching this direction
						}

					}

				}
			}
		}
	}

}

void PocketGrid::mark(core::Vector const & center, core::Real const & vdWd, core::Real const & buffer, bool polar, int targetResi) {
	mark(center(1),center(2),center(3), vdWd, buffer, polar, targetResi);
}

void PocketGrid::mark(core::Real x, core::Real y, core::Real z, core::Real const & vdWd, core::Real const & buffer, bool polar, int targetResi) {
	x-=xcorn_;
	y-=ycorn_;
	z-=zcorn_;
	int xcen=(int)floor(x/stepSize_+.5);
	int ycen=(int)floor(y/stepSize_+.5);
	int zcen=(int)floor(z/stepSize_+.5);
	int minX,maxX;
	int minY,maxY;
	int minZ,maxZ;
	core::Real radius=(vdWd+buffer)/stepSize_;
	core::Real vdW=vdWd/stepSize_;

	minX=int(std::max(ceil(((core::Real)(xcen)-radius)-.5), 0.));
	maxX=int(std::min(floor(((core::Real)(xcen)+radius)+.5), (core::Real)(xdim_)-1.));
	minY=int(std::max(ceil(((core::Real)(ycen)-radius)-.5), 0.));
	maxY=int(std::min(floor(((core::Real)(ycen)+radius)+.5), (core::Real)(ydim_)-1.));
	minZ=int(std::max(ceil(((core::Real)(zcen)-radius)-.5), 0.));
	maxZ=int(std::min(floor(((core::Real)(zcen)+radius)+.5), (core::Real)(zdim_)-1.));

	for ( int xIter=minX; xIter<=maxX; ++xIter ) {
		for ( int yIter=minY; yIter<=maxY; ++yIter ) {
			int centerZ=std::max(minZ, zcen);
			for ( int zIter=centerZ; zIter <= maxZ; ++zIter ) {
				if ( pow((xIter-xcen),2)+pow((yIter-ycen),2)+pow((zIter-zcen),2)>pow(radius, 2) ) continue;
				if ( pow((xIter-xcen),2)+pow((yIter-ycen),2)+pow((zIter-zcen),2)>pow(vdW, 2) ) {
					if ( grid_[xIter][yIter][zIter]!=PROTEIN && grid_[xIter][yIter][zIter]!=TARGET && grid_[xIter][yIter][zIter]!=SUBTARGET ) {
						if ( targetResi==1 ) {
							grid_[xIter][yIter][zIter]=T_SURFACE;
						} else if ( targetResi==2 ) {
							grid_[xIter][yIter][zIter]=ST_SURFACE;
						} else if ( polar ) { if ( grid_[xIter][yIter][zIter]!=HSURFACE ) grid_[xIter][yIter][zIter]=PSURFACE;
						} else { grid_[xIter][yIter][zIter]=HSURFACE;
						}
					} else if ( targetResi==1 ) {
						grid_[xIter][yIter][zIter]=TARGET;
					} else if ( targetResi==2 ) grid_[xIter][yIter][zIter]=SUBTARGET;
				} else {
					if ( targetResi == 1 ) {
						grid_[xIter][yIter][zIter]=TARGET;
					} else if ( targetResi==2 ) {
						grid_[xIter][yIter][zIter]=SUBTARGET;
					} else {
						if ( grid_[xIter][yIter][zIter]!=TARGET && grid_[xIter][yIter][zIter]!=SUBTARGET ) grid_[xIter][yIter][zIter]=PROTEIN;
					}
				}
			}

			centerZ=std::min(maxZ, zcen-1);
			for ( int zIter=centerZ; zIter >= minZ; --zIter ) {
				if ( pow(xIter-xcen,2)+pow(yIter-ycen,2)+pow(zIter-zcen,2)>pow(radius, 2) ) continue;
				if ( pow(xIter-xcen,2)+pow(yIter-ycen,2)+pow(zIter-zcen,2)>pow(vdW, 2) ) {
					if ( grid_[xIter][yIter][zIter]!=PROTEIN && grid_[xIter][yIter][zIter]!=TARGET && grid_[xIter][yIter][zIter]!=SUBTARGET ) {
						if ( targetResi == 1 ) {
							grid_[xIter][yIter][zIter]=T_SURFACE;
						} else if ( targetResi==2 ) {
							grid_[xIter][yIter][zIter]=ST_SURFACE;
						} else if ( polar ) { if ( grid_[xIter][yIter][zIter]!=HSURFACE ) grid_[xIter][yIter][zIter]=PSURFACE;
						} else { grid_[xIter][yIter][zIter]=HSURFACE;
						}
					} else if ( targetResi==1 ) {
						grid_[xIter][yIter][zIter]=TARGET;
					} else if ( targetResi==2 ) grid_[xIter][yIter][zIter]=SUBTARGET;
				} else {
					if ( targetResi == 1 ) {
						grid_[xIter][yIter][zIter]=TARGET;
					} else if ( targetResi==2 ) {
						grid_[xIter][yIter][zIter]=SUBTARGET;
					} else {
						grid_[xIter][yIter][zIter]=PROTEIN;
					}
				}
			}
		}
	}

}

void PocketGrid::clearSmallPockets(core::Size minsize){
	core::Size tmp=minsize;
	tmp++;
}

void PocketGrid::print() const{
	for ( core::Size z=0; z<zdim_; z++ ) {
		for ( core::Size x=0; x<xdim_; x++ ) {
			for ( core::Size y=0; y<ydim_; y++ ) {
				TR << pockets_[x][y][z];
			}
			TR<<std::endl;
		}
		TR<<std::endl;
	}
}

void PocketGrid::dumpGridToFile() {
	std::stringstream filename;
	filename<<tag_<<"pocket"<<pdbno_<<".pdb";
	pdbno_++;
	dumpGridToFile( filename.str() );
}


void PocketGrid::dumpGridToFile( std::string const & output_filename ) {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_filename, std::ios::out);
	int counter=1;
	int counter2=9;

	outPDB_stream<<"ATOM      1  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<0*stepSize_+xcorn_<<std::setw(8)<<0*stepSize_+ycorn_<<std::setw(8)<<0*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      2  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<0*stepSize_+xcorn_<<std::setw(8)<<(ydim_-1)*stepSize_+ycorn_<<std::setw(8)<<0*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      3  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<0*stepSize_+xcorn_<<std::setw(8)<<(ydim_-1)*stepSize_+ycorn_<<std::setw(8)<<(zdim_-1)*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      4  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<0*stepSize_+xcorn_<<std::setw(8)<<0*stepSize_+ycorn_<<std::setw(8)<<(zdim_-1)*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      5  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<(xdim_-1)*stepSize_+xcorn_<<std::setw(8)<<0*stepSize_+ycorn_<<std::setw(8)<<0*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      6  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<(xdim_-1)*stepSize_+xcorn_<<std::setw(8)<<(ydim_-1)*stepSize_+ycorn_<<std::setw(8)<<0*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      7  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<(xdim_-1)*stepSize_+xcorn_<<std::setw(8)<<0*stepSize_+ycorn_<<std::setw(8)<<(zdim_-1)*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      8  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<(xdim_-1)*stepSize_+xcorn_<<std::setw(8)<<(ydim_-1)*stepSize_+ycorn_<<std::setw(8)<<(zdim_-1)*stepSize_+zcorn_<<std::endl;


	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				std::string concatenated_pdb_info;
				concatenated_pdb_info += "ATOM  ";
				std::stringstream  tmp;
				tmp<<counter2;
				if ( counter2<10 ) concatenated_pdb_info += "    ";
				else if ( counter2<100 ) concatenated_pdb_info += "   ";
				else if ( counter2<1000 ) concatenated_pdb_info += "  ";
				else if ( counter2<10000 ) concatenated_pdb_info += " ";
				concatenated_pdb_info += tmp.str()+"  ";
				if ( grid_[x][y][z]==EMPTY ) {
					continue;
				}
				if ( grid_[x][y][z]==POCKET )  {
					continue;
				}
				if ( grid_[x][y][z]==PO_SURF ) {
					continue;
				}
				if ( grid_[x][y][z]==PO_BURIED )  {
					continue;
				}
				if ( grid_[x][y][z]==HSURFACE ) {
					continue;
				}
				if ( grid_[x][y][z]==PSURFACE ) {
					continue;
				}
				if ( grid_[x][y][z]==TP_POCKET )  {
					continue;
				}
				if ( grid_[x][y][z]==TP_SURF )  {
					concatenated_pdb_info += "TPS TPS";
				}
				if ( grid_[x][y][z]==TP_BURIED )  {
					continue;
				}
				if ( grid_[x][y][z]==PO_EDGE )  {
					continue;
				}
				if ( grid_[x][y][z]==TP_EDGE )  {
					concatenated_pdb_info += "TPE TPE";
				}
				if ( grid_[x][y][z]==T_SURFACE )  {
					concatenated_pdb_info += "TS  TS ";
				}
				if ( grid_[x][y][z]==ST_SURFACE )  {
					concatenated_pdb_info += "STS STS";
				}
				if ( grid_[x][y][z]==TARGET )  {
					continue;
				}
				if ( grid_[x][y][z]==SUBTARGET )  {
					continue;
				}
				if ( grid_[x][y][z]==PROTEIN ) {
					concatenated_pdb_info += "PR  PR ";
					//continue;
				}
				if ( grid_[x][y][z]==PSP ) {
					continue;
				}

				tmp.str(std::string());
				tmp<<"          "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x*stepSize_+xcorn_<<std::setw(8)<<y*stepSize_+ycorn_<<std::setw(8)<<z*stepSize_+zcorn_<<std::endl;
				concatenated_pdb_info += tmp.str();
				counter++;
				outPDB_stream<<concatenated_pdb_info;
			}
		}
	}

	int clustNo=1;
	bool smallPocket;
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) smallPocket=true;
		else smallPocket=false;
		for ( auto & point : cluster.points_ ) {
			std::string concatenated_pdb_info;
			concatenated_pdb_info += "ATOM  ";
			std::stringstream  tmp;
			tmp<<counter2;
			if ( counter2<10 ) concatenated_pdb_info += "    ";
			else if ( counter2<100 ) concatenated_pdb_info += "   ";
			else if ( counter2<1000 ) concatenated_pdb_info += "  ";
			else if ( counter2<10000 ) concatenated_pdb_info += " ";
			else concatenated_pdb_info += "";
			concatenated_pdb_info += tmp.str()+"  ";
			if ( smallPocket ) {
				if ( grid_[point.x][point.y][point.z]==TP_POCKET ) concatenated_pdb_info += "STP STP  ";
				if ( grid_[point.x][point.y][point.z]==TP_SURF ) concatenated_pdb_info += "STS STS  ";
				if ( grid_[point.x][point.y][point.z]==TP_BURIED ) concatenated_pdb_info += "STB STB  ";
				if ( grid_[point.x][point.y][point.z]==TP_EDGE ) concatenated_pdb_info += "STE STE  ";
			} else {
				if ( grid_[point.x][point.y][point.z]==TP_POCKET ) concatenated_pdb_info += "TP  TP   ";
				if ( grid_[point.x][point.y][point.z]==TP_SURF ) concatenated_pdb_info += "TPS TPS  ";
				if ( grid_[point.x][point.y][point.z]==TP_BURIED ) concatenated_pdb_info += "TPB TPB  ";
				if ( grid_[point.x][point.y][point.z]==TP_EDGE ) concatenated_pdb_info += "TPE TPE  ";
			}
			if ( grid_[point.x][point.y][point.z]==POCKET ) concatenated_pdb_info += "PC  PC   ";
			if ( grid_[point.x][point.y][point.z]==PO_SURF ) concatenated_pdb_info += "PCS PCS  ";
			if ( grid_[point.x][point.y][point.z]==PO_BURIED ) concatenated_pdb_info += "PCB PCB  ";
			if ( grid_[point.x][point.y][point.z]==PO_EDGE ) concatenated_pdb_info += "PCE PCE  ";

			tmp.str(std::string());
			tmp<<clustNo;
			if ( clustNo<10 ) concatenated_pdb_info += "   ";
			else if ( clustNo<100 ) concatenated_pdb_info += "  ";
			else if ( clustNo<1000 ) concatenated_pdb_info += " ";
			concatenated_pdb_info += tmp.str()+"  ";
			tmp.str(std::string());
			tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.x*stepSize_+xcorn_<<std::setw(8)<<point.y*stepSize_+ycorn_<<std::setw(8)<<point.z*stepSize_+zcorn_<<std::endl;
			concatenated_pdb_info += tmp.str();
			counter2++;
			outPDB_stream<<concatenated_pdb_info;

		}
		clustNo++;
	}

	outPDB_stream.close();
	outPDB_stream.clear();

}

void PocketGrid::dumpExemplarToFile() {
	std::stringstream filename;
	filename<<tag_<<"exemplar"<<expdbno_<<".pdb";
	expdbno_++;
	dumpExemplarToFile( filename.str() );
}


void PocketGrid::dumpExemplarToFile( std::string const & output_filename ) {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_filename, std::ios::out);
	//int counter=1; bazzoli: unused
	int counter2=1;
	int clustNo=1;
	//bool smallPocket;
	int count=0;
	for ( auto & cluster : c_clusters_.clusters_ ) {
		count++;

		if ( !cluster.isTarget(numTargets_) ) {
			continue;
		}

		//Require at least 5 hydrophobic carbon atoms to dump exemplar to a file
		int ccount = 0;
		for ( auto & point : cluster.points_ ) {
			if ( point.atom_type.compare("C") == 0 ) {
				ccount++;
			}
		}
		if ( ccount <5 && limitExemplars_ ) {
			continue;
		}

		for ( auto & point : cluster.points_ ) {
			std::string concatenated_pdb_info;
			concatenated_pdb_info += "HETATM";
			std::stringstream  tmp;
			tmp<<counter2;
			if ( counter2<10 ) concatenated_pdb_info += "    ";
			else if ( counter2<100 ) concatenated_pdb_info += "   ";
			else if ( counter2<1000 ) concatenated_pdb_info += "  ";
			else if ( counter2<10000 ) concatenated_pdb_info += " ";
			else concatenated_pdb_info += "";
			if ( point.atom_type.compare("BBe") == 0 ) {
				tmp << "   Be";
			} else if ( point.atom_type.compare("BNe") == 0 ) {
				tmp << "   Ne";
			} else tmp << "   "<<point.atom_type;
			concatenated_pdb_info += tmp.str()+" ";
			if ( point.atom_type.length() ==1 ) concatenated_pdb_info += " ";
			concatenated_pdb_info += "TMP A";

			tmp.str(std::string());
			tmp<<clustNo;
			if ( clustNo<10 ) concatenated_pdb_info += "   ";
			else if ( clustNo<100 ) concatenated_pdb_info += "  ";
			else if ( clustNo<1000 ) concatenated_pdb_info += " ";
			concatenated_pdb_info += tmp.str()+"  ";
			tmp.str(std::string());
			if ( point.atom_type.compare("C") == 0 ) {
				tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.x*stepSize_+xcorn_<<std::setw(8)<<point.y*stepSize_+ycorn_<<std::setw(8)<<point.z*stepSize_+zcorn_<<"  1.00  2.03           "<<point.atom_type<<std::endl;
			} else {
				if ( point.atom_type.compare("BBe") == 0 ) {
					tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.absX<<std::setw(8)<<point.absY<<std::setw(8)<<point.absZ<<"  1.00  2.03          Be"<<std::endl;
				} else if ( point.atom_type.compare("BNe") == 0 ) {
					tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.absX<<std::setw(8)<<point.absY<<std::setw(8)<<point.absZ<<"  1.00  2.03          Ne"<<std::endl;
				} else {
					tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.absX<<std::setw(8)<<point.absY<<std::setw(8)<<point.absZ<<"  1.00  2.03          "<<point.atom_type<<std::endl;
				}
			}

			concatenated_pdb_info += tmp.str();
			counter2++;
			outPDB_stream<<concatenated_pdb_info;

		}
		clustNo++;
	}


	outPDB_stream.close();
	outPDB_stream.clear();

}


void PocketGrid::dumpTargetPocketsToPDB( std::string const & output_filename, bool minipock ){

	//this sets up a temporary grid I can tick off pocket points from for minimal pockets
	std::vector < std::vector < std::vector <PtType> > > tmpgrid_;
	if ( minipock ) {
		tmpgrid_.resize(xdim_);
		for ( core::Size tx=0; tx<xdim_; tx++ ) {
			tmpgrid_[tx].resize(ydim_);
			for ( core::Size ty=0; ty<ydim_; ty++ ) {
				tmpgrid_[tx][ty].resize(zdim_, EMPTY);
			}
		}
		for ( auto & cluster : clusters_.clusters_ ) {
			if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
			if ( !cluster.isTarget(numTargets_) ) continue;
			for ( auto & point : cluster.points_ ) {
				tmpgrid_[point.x][point.y][point.z] = grid_[point.x][point.y][point.z];
			}
		}
		/*    for (core::Size tx=0;tx<xdim_;tx++)
		for (core::Size ty=0;ty<ydim_;ty++)
		for (core::Size tz=0;tz<zdim_;tz++)
		if (grid_[tx][ty][tz] == TP_SURF) tmpgrid_[tx][ty][tz] = TP_SURF;
		*/
	}

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_filename, std::ios::out);
	//int counter=1;
	int counter2=1;

	int clustNo=1;
	//bool smallPocket;
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			if ( minipock ) {
				if ( tmpgrid_[point.x][point.y][point.z] == EMPTY ) continue;
				if ( point.x < 2 || point.y < 2|| point.z < 2 || point.x > xdim_-3 || point.y > ydim_-3|| point.z > zdim_-3 ) continue;
				if ( tmpgrid_[point.x+2][point.y][point.z] != EMPTY && tmpgrid_[point.x-2][point.y][point.z] != EMPTY && tmpgrid_[point.x][point.y+2][point.z] != EMPTY && tmpgrid_[point.x][point.y-2][point.z] != EMPTY &&tmpgrid_[point.x][point.y][point.z+2] != EMPTY && tmpgrid_[point.x][point.y][point.z-2] != EMPTY ) {
					bool pockpt=true;
					for ( int i = -1; i <2; i++ ) {
						for ( int j = -1; j <2; j++ ) {
							for ( int k = -1; k <2; k++ ) {
								if ( i==0 && j==0 && k==0 ) continue;
								if ( tmpgrid_[point.x+i][point.y+j][point.z+k] == EMPTY ) {
									pockpt=false;
								}
							}
						}
					}
					if ( pockpt ) {
						for ( int i = -1; i <2; i++ ) {
							for ( int j = -1; j <2; j++ ) {
								for ( int k = -1; k <2; k++ ) {
									tmpgrid_[point.x+i][point.y+j][point.z+k] = EMPTY;
								}
							}
						}
						tmpgrid_[point.x-2][point.y][point.z] = EMPTY;
						tmpgrid_[point.x][point.y-2][point.z] = EMPTY;
						tmpgrid_[point.x][point.y][point.z-2] = EMPTY;
					} else continue;
				} else continue;
			}
			std::string concatenated_pdb_info;
			if ( minipock ) {
				concatenated_pdb_info += "HETATM";
			} else {
				concatenated_pdb_info += "ATOM  ";
			}
			std::stringstream  tmp;
			tmp<<counter2;
			if ( counter2<10 ) concatenated_pdb_info += "    ";
			else if ( counter2<100 ) concatenated_pdb_info += "   ";
			else if ( counter2<1000 ) concatenated_pdb_info += "  ";
			else if ( counter2<10000 ) concatenated_pdb_info += " ";
			else concatenated_pdb_info += "";
			concatenated_pdb_info += tmp.str()+"  ";
			if ( minipock ) {
				concatenated_pdb_info += " C  TMP A";
			} else {
				if ( grid_[point.x][point.y][point.z]==TP_POCKET ) concatenated_pdb_info += "TP  TP   ";
				if ( grid_[point.x][point.y][point.z]==TP_SURF ) concatenated_pdb_info += "TPS TPS  ";
				if ( grid_[point.x][point.y][point.z]==TP_BURIED ) concatenated_pdb_info += "TPB TPB  ";
				if ( grid_[point.x][point.y][point.z]==TP_EDGE ) concatenated_pdb_info += "TPE TPE  ";
			}
			tmp.str(std::string());
			tmp<<clustNo;
			if ( clustNo<10 ) concatenated_pdb_info += "   ";
			else if ( clustNo<100 ) concatenated_pdb_info += "  ";
			else if ( clustNo<1000 ) concatenated_pdb_info += " ";
			concatenated_pdb_info += tmp.str()+"  ";
			tmp.str(std::string());
			tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.x*stepSize_+xcorn_<<std::setw(8)<<point.y*stepSize_+ycorn_<<std::setw(8)<<point.z*stepSize_+zcorn_;
			if ( minipock ) tmp << "  1.00  2.03           C";
			tmp << std::endl;
			concatenated_pdb_info += tmp.str();
			counter2++;
			outPDB_stream<<concatenated_pdb_info;

		}
		clustNo++;
	}

	if ( minipock ) {
		for ( core::Size tx=2; tx<xdim_-3; tx++ ) {
			for ( core::Size ty=2; ty<ydim_-3; ty++ ) {
				for ( core::Size tz=2; tz<zdim_-3; tz++ ) {
					if ( tmpgrid_[tx][ty][tz] == EMPTY ) continue;
					if ( tmpgrid_[tx+2][ty][tz] != EMPTY && tmpgrid_[tx-2][ty][tz] != EMPTY && tmpgrid_[tx][ty+2][tz] != EMPTY && tmpgrid_[tx][ty-2][tz] != EMPTY &&tmpgrid_[tx][ty][tz+2] != EMPTY && tmpgrid_[tx][ty][tz-2] != EMPTY ) {
						bool pockpt=true;
						for ( int i = -1; i <2; i++ ) {
							for ( int j = -1; j <2; j++ ) {
								for ( int k = -1; k <2; k++ ) {
									if ( i==0 && j==0 && k==0 ) continue;
									if ( tmpgrid_[tx+i][ty+j][tz+k] == EMPTY ) {
										pockpt=false;
									}
								}
							}
						}
						if ( pockpt ) {
							for ( int i = -1; i <2; i++ ) {
								for ( int j = -1; j <2; j++ ) {
									for ( int k = -1; k <2; k++ ) {
										tmpgrid_[tx+i][ty+j][tz+k] = EMPTY;
									}
								}
							}
							tmpgrid_[tx-2][ty][tz] = EMPTY;
							tmpgrid_[tx][ty-2][tz] = EMPTY;
							tmpgrid_[tx][ty][tz-2] = EMPTY;
						} else continue;
					} else continue;
					std::string concatenated_pdb_info;
					concatenated_pdb_info += "HETATM";
					std::stringstream  tmp;
					tmp<<counter2;
					if ( counter2<10 ) concatenated_pdb_info += "    ";
					else if ( counter2<100 ) concatenated_pdb_info += "   ";
					else if ( counter2<1000 ) concatenated_pdb_info += "  ";
					else if ( counter2<10000 ) concatenated_pdb_info += " ";
					else concatenated_pdb_info += "";
					concatenated_pdb_info += tmp.str()+"  ";
					concatenated_pdb_info += " C  TMP A";
					tmp.str(std::string());
					tmp<<clustNo;
					if ( clustNo<10 ) concatenated_pdb_info += "   ";
					else if ( clustNo<100 ) concatenated_pdb_info += "  ";
					else if ( clustNo<1000 ) concatenated_pdb_info += " ";
					concatenated_pdb_info += tmp.str()+"  ";
					tmp.str(std::string());
					tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<tx*stepSize_+xcorn_<<std::setw(8)<<ty*stepSize_+ycorn_<<std::setw(8)<<tz*stepSize_+zcorn_;
					if ( minipock ) tmp << "  1.00  2.03           C";
					tmp << std::endl;
					concatenated_pdb_info += tmp.str();
					counter2++;
					outPDB_stream<<concatenated_pdb_info;


					clustNo++;
				}
			}
		}
	}


	outPDB_stream.close();
	outPDB_stream.clear();

}
void PocketGrid::dumpTargetPocketsToPDB( std::string const & output_filename, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 ){

	core::Real xmin = 1,xmax=0,ymin=1,ymax=0,zmin=1,zmax=0;
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			numeric::xyzVector<core::Real> coord(point.x,point.y,point.z);
			coord = rot1 * coord;
			coord = rot2 * coord;
			coord = rot3 * coord;
			if ( xmin>xmax ) {
				xmin = coord(1);
				xmax = coord(1);
				ymin = coord(2);
				ymax = coord(2);
				zmin = coord(3);
				zmax = coord(3);
			} else {
				if ( xmin > coord(1) ) xmin = coord(1);
				if ( xmax < coord(1) ) xmax = coord(1);
				if ( ymin > coord(2) ) ymin = coord(2);
				if ( ymax < coord(2) ) ymax = coord(2);
				if ( zmin > coord(3) ) zmin = coord(3);
				if ( zmax < coord(3) ) zmax = coord(3);
			}
		}
	}

	//core::Real new_xcorn=(core::Real)std::floor(xmin + 0.5);
	//core::Real new_ycorn=(core::Real)std::floor(ymin + 0.5);
	//core::Real new_zcorn=(core::Real)std::floor(zmin + 0.5);

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_filename, std::ios::out);
	//int counter=1;
	int counter2=1;

	int clustNo=1;
	//bool smallPocket;
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			std::string concatenated_pdb_info;
			concatenated_pdb_info += "ATOM  ";
			std::stringstream  tmp;
			tmp<<counter2;
			if ( counter2<10 ) concatenated_pdb_info += "    ";
			else if ( counter2<100 ) concatenated_pdb_info += "   ";
			else if ( counter2<1000 ) concatenated_pdb_info += "  ";
			else if ( counter2<10000 ) concatenated_pdb_info += " ";
			else concatenated_pdb_info += "";
			concatenated_pdb_info += tmp.str()+"  ";
			if ( grid_[point.x][point.y][point.z]==TP_POCKET ) concatenated_pdb_info += "TP  TP   ";
			if ( grid_[point.x][point.y][point.z]==TP_SURF ) concatenated_pdb_info += "TPS TPS  ";
			if ( grid_[point.x][point.y][point.z]==TP_BURIED ) concatenated_pdb_info += "TPB TPB  ";
			if ( grid_[point.x][point.y][point.z]==TP_EDGE ) concatenated_pdb_info += "TPE TPE  ";

			tmp.str(std::string());
			tmp<<clustNo;
			if ( clustNo<10 ) concatenated_pdb_info += "   ";
			else if ( clustNo<100 ) concatenated_pdb_info += "  ";
			else if ( clustNo<1000 ) concatenated_pdb_info += " ";
			concatenated_pdb_info += tmp.str()+"  ";
			tmp.str(std::string());
			numeric::xyzVector<core::Real> coord(point.x*stepSize_+xcorn_,point.y*stepSize_+ycorn_,point.z*stepSize_+zcorn_);
			coord = rot1 * coord;
			coord = rot2 * coord;
			coord = rot3 * coord;

			tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<coord(1)<<std::setw(8)<<coord(2)<<std::setw(8)<<coord(3)<<std::endl;
			concatenated_pdb_info += tmp.str();
			counter2++;
			outPDB_stream<<concatenated_pdb_info;

		}
		clustNo++;
	}
	outPDB_stream.close();
	outPDB_stream.clear();

}

void PocketGrid::dumpTargetPocketsToFile( std::string const & output_filename ){


	// NOTE: THIS FUNCTION DOESN'T WRITE THE WHOLE GRID_
	//         HERE WE ONLY PROVIDE MINIMAL FUNCTIONALITY TO DO POCKET VS POCKET COMPARISONS....

	utility::io::ozstream outstream;
	outstream.open(output_filename, std::ios::out);

	// Print stepSize_ xdim_ ydim_ zdim_ xcorn_ ycorn_ zcorn_
	outstream << stepSize_ << " " << xdim_ << " " << ydim_ << " " << zdim_ << " " << xcorn_ << " " << ycorn_ << " " << zcorn_ << std::endl;

	// Need to print each grid point if it's of type EGGSHELL along with the corresponding cartesian coors ( for eggshell_coord_list_ )
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			if ( grid_[point.x][point.y][point.z]==TP_POCKET ) {
				outstream << point.x << " " << point.y << " " << point.z << " TP_POCKET" << std::endl;
			}
			if ( grid_[point.x][point.y][point.z]==TP_BURIED ) {
				outstream << point.x << " " << point.y << " " << point.z << " TP_BURIED" << std::endl;
			}

		}

	}

	outstream.close();
	outstream.clear();

}

void PocketGrid::dumpTargetPocketsToFile( std::string const & output_filename, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 ){

	core::Real xmin = 1,xmax=0,ymin=1,ymax=0,zmin=1,zmax=0;
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			numeric::xyzVector<core::Real> coord(xcorn_ + stepSize_ * point.x,ycorn_ + stepSize_ * point.y,zcorn_ + stepSize_ * point.z);
			coord = rot1 * coord;
			coord = rot2 * coord;
			coord = rot3 * coord;
			if ( xmin>xmax ) {
				xmin = coord(1);
				xmax = coord(1);
				ymin = coord(2);
				ymax = coord(2);
				zmin = coord(3);
				zmax = coord(3);
			} else {
				if ( xmin > coord(1) ) xmin = coord(1);
				if ( xmax < coord(1) ) xmax = coord(1);
				if ( ymin > coord(2) ) ymin = coord(2);
				if ( ymax < coord(2) ) ymax = coord(2);
				if ( zmin > coord(3) ) zmin = coord(3);
				if ( zmax < coord(3) ) zmax = coord(3);
			}
		}
	}

	core::Real new_xcorn=(core::Real)std::floor(xmin);
	core::Real new_ycorn=(core::Real)std::floor(ymin);
	core::Real new_zcorn=(core::Real)std::floor(zmin);
	core::Real x_far_corn=(core::Real)std::ceil(xmax);
	core::Real y_far_corn=(core::Real)std::ceil(ymax);
	core::Real z_far_corn=(core::Real)std::ceil(zmax);
	core::Size new_xdim = (core::Size)std::ceil((x_far_corn-new_xcorn)/stepSize_)+1;
	core::Size new_ydim = (core::Size)std::ceil((y_far_corn-new_ycorn)/stepSize_)+1;
	core::Size new_zdim = (core::Size)std::ceil((z_far_corn-new_zcorn)/stepSize_)+1;

	// NOTE: THIS FUNCTION DOESN'T WRITE THE WHOLE GRID_
	//         HERE WE ONLY PROVIDE MINIMAL FUNCTIONALITY TO DO POCKET VS POCKET COMPARISONS....

	utility::io::ozstream outstream;
	outstream.open(output_filename, std::ios::out);

	// Print stepSize_ xdim_ ydim_ zdim_ xcorn_ ycorn_ zcorn_
	outstream << stepSize_ << " " << new_xdim << " " << new_ydim << " " << new_zdim << " " << new_xcorn << " " << new_ycorn << " " << new_zcorn << std::endl;

	// Need to print each grid point if it's of type EGGSHELL along with the corresponding cartesian coors ( for eggshell_coord_list_ )
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			numeric::xyzVector<core::Real> coord(xcorn_ + point.x * stepSize_, ycorn_ + point.y * stepSize_, zcorn_ + point.z * stepSize_);
			coord = rot1 * coord;
			coord = rot2 * coord;
			coord = rot3 * coord;
			numeric::xyzVector<core::Size> newCoord(0,0,0);
			newCoord(1) = (core::Size)std::floor((coord(1) - new_xcorn)/stepSize_ + .5);
			newCoord(2) = (core::Size)std::floor((coord(2) - new_ycorn)/stepSize_ + .5);
			newCoord(3) = (core::Size)std::floor((coord(3) - new_zcorn)/stepSize_ + .5);

			if (  grid_[point.x][point.y][point.z]==TP_POCKET )  {
				outstream << newCoord(1) << " " << newCoord(2) << " " << newCoord(3) << " TP_POCKET" << std::endl;
			}
			if ( grid_[point.x][point.y][point.z]==TP_BURIED )  {
				outstream << newCoord(1) << " " << newCoord(2) << " " << newCoord(3) << " TP_BURIED" << std::endl;
			}

		}
	}


	outstream.close();
	outstream.clear();

}


void PocketGrid::clusterPockets(){

	//vector of deltas to add to x,y,z.  Index % 3 = 0 for x deltas, 1 for y, 2 for z
	std::vector<int> deltas;
	int dirs=13;
	deltas.push_back(0);deltas.push_back(-1);deltas.push_back(-1);
	deltas.push_back(0);deltas.push_back(-1);deltas.push_back(0);
	deltas.push_back(0);deltas.push_back(-1);deltas.push_back(1);
	deltas.push_back(0);deltas.push_back(0);deltas.push_back(-1);
	//deltas.push_back(0);deltas.push_back(1);deltas.push_back(-1);
	deltas.push_back(-1);deltas.push_back(-1);deltas.push_back(-1);
	deltas.push_back(-1);deltas.push_back(-1);deltas.push_back(0);
	deltas.push_back(-1);deltas.push_back(-1);deltas.push_back(1);
	deltas.push_back(-1);deltas.push_back(0);deltas.push_back(-1);
	deltas.push_back(-1);deltas.push_back(0);deltas.push_back(0);
	deltas.push_back(-1);deltas.push_back(0);deltas.push_back(1);
	deltas.push_back(-1);deltas.push_back(1);deltas.push_back(-1);
	deltas.push_back(-1);deltas.push_back(1);deltas.push_back(0);
	deltas.push_back(-1);deltas.push_back(1);deltas.push_back(1);

	// a grid of pointers back to clusters.  Non-clustered points point to end()
	std::vector < std::vector < std::vector <std::list<PCluster>::iterator> > > cluster_grid;
	cluster_grid.resize(xdim_);
	for ( core::Size tx=0; tx<xdim_; tx++ ) {
		cluster_grid[tx].resize(ydim_);
		for ( core::Size ty=0; ty<ydim_; ty++ ) {
			cluster_grid[tx][ty].resize(zdim_, clusters_.clusters_.end());
		}
	}

	// Search the grid for pockets
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				// If you find a pocket, check if its neighbors are already part of a cluster
				if ( grid_[x][y][z]==TP_POCKET || grid_[x][y][z]==POCKET || grid_[x][y][z]==PO_EDGE || grid_[x][y][z]==TP_EDGE || grid_[x][y][z]==PO_BURIED || grid_[x][y][z]==TP_BURIED ) {
					int found=0;
					for ( int i=0; i<dirs*3; i+=3 ) {
						int cx = (int)x+deltas[i];
						int cy = (int)y+deltas[i+1];
						int cz = (int)z+deltas[i+2];
						if ( cx<0 || cy<0 || cz<0 || cx>=(int)xdim_ || cy>=(int)ydim_ || cz>=(int)zdim_ ) continue;

						// If a neighbor is found...
						if ( cluster_grid[cx][cy][cz]!=clusters_.clusters_.end() ) {
							// If this is the first neighbor found, add it to that cluster
							if ( !found ) {
								found=1;
								cluster_grid[x][y][z]=cluster_grid[cx][cy][cz];
								cluster_grid[x][y][z]->add(x,y,z);
							} else {
								// If this cluster differs from the one that you just added, merge them
								if ( cluster_grid[cx][cy][cz]!=cluster_grid[x][y][z] ) {
									auto old = cluster_grid[cx][cy][cz];
									for ( auto & point : old->points_ ) {
										cluster_grid[point.x][point.y][point.z]=cluster_grid[x][y][z];
									}
									clusters_.join(cluster_grid[x][y][z], old);
								}
							}
						}
					}
					// If none of the neighbors are in a cluster, create a new cluster
					if ( !found ) {
						cluster_grid[x][y][z]=clusters_.add(x,y,z, stepSize_);
					}
				}
			}
		}
	}
}

void PocketGrid::clusterCPockets(){

	//vector of deltas to add to x,y,z.  Index % 3 = 0 for x deltas, 1 for y, 2 for z
	std::vector<int> deltas;
	int dirs=13;
	deltas.push_back(0);deltas.push_back(-1);deltas.push_back(-1);
	deltas.push_back(0);deltas.push_back(-1);deltas.push_back(0);
	deltas.push_back(0);deltas.push_back(-1);deltas.push_back(1);
	deltas.push_back(0);deltas.push_back(0);deltas.push_back(-1);
	//deltas.push_back(0);deltas.push_back(1);deltas.push_back(-1);
	deltas.push_back(-1);deltas.push_back(-1);deltas.push_back(-1);
	deltas.push_back(-1);deltas.push_back(-1);deltas.push_back(0);
	deltas.push_back(-1);deltas.push_back(-1);deltas.push_back(1);
	deltas.push_back(-1);deltas.push_back(0);deltas.push_back(-1);
	deltas.push_back(-1);deltas.push_back(0);deltas.push_back(0);
	deltas.push_back(-1);deltas.push_back(0);deltas.push_back(1);
	deltas.push_back(-1);deltas.push_back(1);deltas.push_back(-1);
	deltas.push_back(-1);deltas.push_back(1);deltas.push_back(0);
	deltas.push_back(-1);deltas.push_back(1);deltas.push_back(1);

	// a grid of pointers back to clusters.  Non-clustered points point to end()
	std::vector < std::vector < std::vector <std::list<PCluster>::iterator> > > cluster_grid;
	cluster_grid.resize(xdim_);
	for ( core::Size tx=0; tx<xdim_; tx++ ) {
		cluster_grid[tx].resize(ydim_);
		for ( core::Size ty=0; ty<ydim_; ty++ ) {
			cluster_grid[tx][ty].resize(zdim_, clusters_.clusters_.end());
		}
	}

	// Search the grid for pockets
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				// If you find a pocket, check if its neighbors are already part of a cluster
				if ( c_grid_[x][y][z]==TP_POCKET || c_grid_[x][y][z]==POCKET || c_grid_[x][y][z]==PO_EDGE || c_grid_[x][y][z]==TP_EDGE || c_grid_[x][y][z]==PO_BURIED || c_grid_[x][y][z]==TP_BURIED ) {
					int found=0;
					for ( int i=0; i<dirs*3; i+=3 ) {
						int cx = (int)x+deltas[i];
						int cy = (int)y+deltas[i+1];
						int cz = (int)z+deltas[i+2];
						if ( cx<0 || cy<0 || cz<0 || cx>=(int)xdim_ || cy>=(int)ydim_ || cz>=(int)zdim_ ) continue;

						// If a neighbor is found...
						if ( cluster_grid[cx][cy][cz]!=clusters_.clusters_.end() ) {
							// If this is the first neighbor found, add it to that cluster
							if ( !found ) {
								found=1;
								cluster_grid[x][y][z]=cluster_grid[cx][cy][cz];
								cluster_grid[x][y][z]->add(x,y,z);
								cluster_grid[x][y][z]->target_=true;
								cluster_grid[x][y][z]->subtarget_=true;
							} else {
								// If this cluster differs from the one that you just added, merge them
								if ( cluster_grid[cx][cy][cz]!=cluster_grid[x][y][z] ) {
									auto old = cluster_grid[cx][cy][cz];
									for ( auto & point : old->points_ ) {
										cluster_grid[point.x][point.y][point.z]=cluster_grid[x][y][z];
									}
									clusters_.join(cluster_grid[x][y][z], old);
								}
							}
						}
					}
					// If none of the neighbors are in a cluster, create a new cluster
					if ( !found ) {
						cluster_grid[x][y][z]=clusters_.add(x,y,z, stepSize_);
						cluster_grid[x][y][z]->target_=true;
						cluster_grid[x][y][z]->subtarget_=true;
					}
				}
			}
		}
	}
	for ( auto & cluster : clusters_.clusters_ ) {
		cluster.target_=true;
		cluster.subtarget_=true;
		cluster.solventExposed_ = true;
	}
}


void PocketGrid::markPocketDepth(core::Real const & surf_d, core::Real const & bur_d){
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				if ( grid_[x][y][z]==TP_POCKET || grid_[x][y][z]==POCKET ) markDepth(x,y,z, surf_d, bur_d);
			}
		}
	}
}

void PocketGrid::markEdgeDepth(core::Real const & surf_d, core::Real const & bur_d ){
	for ( auto & cluster : clusters_.clusters_ ) {
		for ( auto pit=cluster.points_.begin(); pit != cluster.points_.end(); ++pit ) {
			if ( grid_[pit->x][pit->y][pit->z]==TP_EDGE || grid_[pit->x][pit->y][pit->z]==PO_EDGE ) {
				bool surf=markOneEdgeDepth(pit->x,pit->y,pit->z, surf_d, bur_d, cluster.isTarget(numTargets_));
				if ( surf ) {
					auto tmp = pit;
					--pit;
					cluster.points_.erase(tmp);
				}
			}
		}
	}
}

void PocketGrid::markDepth(core::Size x, core::Size y, core::Size z, core::Real const & surf_d, core::Real const & bur_d){
	core::Size ssteps=(core::Size)ceil(surf_d/stepSize_);
	core::Size st;
	core::Size en;
	if ( ssteps>x ) {
		st=0;
		if ( grid_[x][y][z]==POCKET ) {
			grid_[x][y][z]=PO_EDGE;
		} else grid_[x][y][z]=TP_EDGE;
		//clusters_.add(x,y,z, stepSize_);
		return;
	} else st=x-ssteps;
	if ( x+ssteps>=xdim_ ) {
		en=xdim_;
		if ( grid_[x][y][z]==POCKET ) {
			grid_[x][y][z]=PO_EDGE;
		} else grid_[x][y][z]=TP_EDGE;
		//clusters_.add(x,y,z, stepSize_);
		return ;
	} else en=x+ssteps+1;
	for ( core::Size i=st; i<en; i++ ) {
		if ( grid_[i][y][z]==EMPTY ) {
			if ( grid_[x][y][z]==POCKET ) {
				grid_[x][y][z]=PO_SURF;
			} else grid_[x][y][z]=TP_SURF;
			return ;
		}
	}
	if ( ssteps>y ) {
		st=0;
		if ( grid_[x][y][z]==POCKET ) {
			grid_[x][y][z]=PO_EDGE;
		} else grid_[x][y][z]=TP_EDGE;
		//clusters_.add(x,y,z, stepSize_);
		return ;
	} else st=y-ssteps;
	if ( y+ssteps>=ydim_ ) {
		en=ydim_;
		if ( grid_[x][y][z]==POCKET ) {
			grid_[x][y][z]=PO_EDGE;
		} else grid_[x][y][z]=TP_EDGE;
		//clusters_.add(x,y,z, stepSize_);
		return ;
	} else en=y+ssteps+1;
	for ( core::Size i=st; i<en; i++ ) {
		if ( grid_[x][i][z]==EMPTY ) {
			if ( grid_[x][y][z]==POCKET ) {
				grid_[x][y][z]=PO_SURF;
			} else grid_[x][y][z]=TP_SURF;
			return ;
		}
	}
	if ( ssteps>z ) {
		st=0;
		if ( grid_[x][y][z]==POCKET ) {
			grid_[x][y][z]=PO_EDGE;
		} else grid_[x][y][z]=TP_EDGE;
		//clusters_.add(x,y,z, stepSize_);
		return ;
	} else st=z-ssteps;
	if ( z+ssteps>=zdim_ ) {
		en=zdim_;
		if ( grid_[x][y][z]==POCKET ) {
			grid_[x][y][z]=PO_EDGE;
		} else grid_[x][y][z]=TP_EDGE;
		//clusters_.add(x,y,z, stepSize_);
		return ;
	} else en=z+ssteps+1;
	for ( core::Size i=st; i<en; i++ ) {
		if ( grid_[x][y][i]==EMPTY ) {
			if ( grid_[x][y][z]==POCKET ) {
				grid_[x][y][z]=PO_SURF;
			} else grid_[x][y][z]=TP_SURF;
			return ;
		}
	}

	ssteps=(core::Size)ceil(bur_d/stepSize_);
	if ( ssteps>x ) st=0;
	else st=x-ssteps;
	if ( x+ssteps>=xdim_ ) en=xdim_;
	else en=x+ssteps+1;
	for ( core::Size i=st; i<en; i++ ) {
		if ( isSurfacePoint(i,y,z) ) {
			if ( grid_[x][y][z]==POCKET ) {
				grid_[x][y][z]=PO_BURIED;
			} else grid_[x][y][z]=TP_BURIED;
			//clusters_.add(x,y,z, stepSize_);
			return ;
		}
	}
	if ( ssteps>y ) st=0;
	else st=y-ssteps;
	if ( y+ssteps>=ydim_ ) en=ydim_;
	else en=y+ssteps+1;
	for ( core::Size i=st; i<en; i++ ) {
		if ( isSurfacePoint(x,i,z) ) {
			if ( grid_[x][y][z]==POCKET ) {
				grid_[x][y][z]=PO_BURIED;
			} else grid_[x][y][z]=TP_BURIED;
			//clusters_.add(x,y,z, stepSize_);
			return ;
		}
	}
	if ( ssteps>z ) st=0;
	else st=z-ssteps;
	if ( z+ssteps>=zdim_ ) en=zdim_;
	else en=z+ssteps+1;
	for ( core::Size i=st; i<en; i++ ) {
		if ( isSurfacePoint(x,y,i) ) {
			if ( grid_[x][y][z]==POCKET ) {
				grid_[x][y][z]=PO_BURIED;
			} else grid_[x][y][z]=TP_BURIED;
			//clusters_.add(x,y,z, stepSize_);
			return ;
		}
	}
	//clusters_.add(x,y,z, stepSize_);
}

bool PocketGrid::markOneEdgeDepth(core::Size x, core::Size y, core::Size z, core::Real const & surf_d, core::Real const & bur_d, bool isTarget){
	core::Size xst;
	core::Size xen;
	core::Size yst;
	core::Size yen;
	core::Size zst;
	core::Size zen;
	core::Size ssteps=(core::Size)ceil(surf_d/stepSize_);
	if ( ssteps>x ) {
		xst=0;
	} else xst=x-ssteps;
	if ( x+ssteps>=xdim_ ) {
		xen=xdim_;
	} else xen=x+ssteps+1;
	if ( ssteps>y ) {
		yst=0;
	} else yst=y-ssteps;
	if ( y+ssteps>=ydim_ ) {
		yen=ydim_;
	} else yen=y+ssteps+1;
	if ( ssteps>z ) {
		zst=0;
	} else zst=z-ssteps;
	if ( z+ssteps>=zdim_ ) {
		zen=zdim_;
	} else zen=z+ssteps+1;
	for ( core::Size i=xst; i<xen; i++ ) {
		if ( grid_[i][y][z]==EMPTY ) {
			if ( grid_[x][y][z]==PO_EDGE||!isTarget ) {
				grid_[x][y][z]=PO_SURF;
			} else grid_[x][y][z]=TP_SURF;
			return true;
		}
	}
	for ( core::Size i=yst; i<yen; i++ ) {
		if ( grid_[x][i][z]==EMPTY ) {
			if ( grid_[x][y][z]==PO_EDGE||!isTarget ) {
				grid_[x][y][z]=PO_SURF;
			} else grid_[x][y][z]=TP_SURF;
			return true;
		}
	}
	for ( core::Size i=zst; i<zen; i++ ) {
		if ( grid_[x][y][i]==EMPTY ) {
			if ( grid_[x][y][z]==PO_EDGE||!isTarget ) {
				grid_[x][y][z]=PO_SURF;
			} else grid_[x][y][z]=TP_SURF;
			return true;
		}
	}

	ssteps=(core::Size)ceil(bur_d/stepSize_);
	if ( ssteps>x ) {
		xst=0;
	} else xst=x-ssteps;
	if ( x+ssteps>=xdim_ ) {
		xen=xdim_;
	} else xen=x+ssteps+1;
	if ( ssteps>y ) {
		yst=0;
	} else yst=y-ssteps;
	if ( y+ssteps>=ydim_ ) {
		yen=ydim_;
	} else yen=y+ssteps+1;
	if ( ssteps>z ) {
		zst=0;
	} else zst=z-ssteps;
	if ( z+ssteps>=zdim_ ) {
		zen=zdim_;
	} else zen=z+ssteps+1;
	for ( core::Size i=xst; i<xen; i++ ) {
		if ( grid_[i][y][z]==HSURFACE || grid_[i][y][z]==PSURFACE ) {
			if ( grid_[x][y][z]==POCKET||!isTarget ) {
				grid_[x][y][z]=PO_BURIED;
			} else grid_[x][y][z]=TP_BURIED;
			return false;
		}
	}
	for ( core::Size i=yst; i<yen; i++ ) {
		if ( grid_[x][i][z]==HSURFACE || grid_[x][i][z]==PSURFACE ) {
			if ( grid_[x][y][z]==POCKET||!isTarget ) {
				grid_[x][y][z]=PO_BURIED;
			} else grid_[x][y][z]=TP_BURIED;
			return false;
		}
	}

	for ( core::Size i=zst; i<zen; i++ ) {
		if ( grid_[x][y][i]==HSURFACE || grid_[x][y][i]==PSURFACE ) {
			if ( grid_[x][y][z]==POCKET||!isTarget ) {
				grid_[x][y][z]=PO_BURIED;
			} else grid_[x][y][z]=TP_BURIED;
			return false;
		}
	}

	if ( grid_[x][y][z]==PO_EDGE||!isTarget ) {
		grid_[x][y][z]=POCKET;
	} else grid_[x][y][z]=TP_POCKET;
	return false;

}

core::Real PocketGrid::targetPocketVolume(core::Real const & surf_sc, core::Real const & bur_sc) const {
	core::Real vol=0;
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				if ( grid_[x][y][z]==TP_POCKET ) vol+=1;
				if ( grid_[x][y][z]==TP_SURF ) vol+=surf_sc;
				if ( grid_[x][y][z]==TP_BURIED ) vol+=bur_sc;
			}
		}
	}
	return vol*pow(stepSize_,3);
}

core::Real PocketGrid::targetPocketSolventSurface() const {
	core::Real vol=0;
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				if ( grid_[x][y][z]==TP_POCKET ) {
					if ( touchesSolvent(x,y,z) ) vol+=1;
				}
			}
		}
	}
	vol*=pow(stepSize_,2);
	return vol;


}

core::Real PocketGrid::targetPocketProteinSurface() const {
	core::Real vol=0;
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				if ( grid_[x][y][z]==TP_POCKET ) {
					if ( touchesSurface(x,y,z, true, true) ) vol+=1;
				}
			}
		}
	}
	vol*=pow(stepSize_,2);
	return vol;


}

core::Real PocketGrid::targetPocketHydrophobicProteinSurface() const {
	core::Real vol=0;
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				if ( grid_[x][y][z]==TP_POCKET ) {
					if ( touchesSurface(x,y,z, false) ) vol+=1;
				}
			}
		}
	}
	vol*=pow(stepSize_,2);
	return vol;


}

core::Real PocketGrid::targetPocketPolarProteinSurface() const {
	core::Real vol=0;
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				if ( grid_[x][y][z]==TP_POCKET ) {
					if ( touchesSurface(x,y,z, true) ) vol+=1;
				}
			}
		}
	}
	vol*=pow(stepSize_,2);
	return vol;


}

core::Real PocketGrid::targetPocketHeuristicScore() const {
	core::Real vol=0;
	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				if ( grid_[x][y][z]==TP_POCKET ) {
					if ( !touchesSS(x,y,z) ) {
						core::Size protTouch=touchesPS(x,y,z);
						if ( protTouch>0 ) {
							vol+=3*protTouch;
						} else {
							vol+=1;
						}
					}
				}
			}
		}
	}
	vol*=pow(stepSize_,2);
	return vol;


}

void PocketGrid::findClusters(){
	//      clusters_.findClusters();
	clusterPockets();
	//mark target clusters
	for ( auto & cluster : clusters_.clusters_ ) {
		for ( auto pit=cluster.points_.begin(); pit != cluster.points_.end(); ++pit ) {
			for ( int x = -1; x <2; x++ ) {
				for ( int y = -1; y <2; y++ ) {
					for ( int z = -1; z <2; z++ ) {
						if ( x==-1 ) if ( pit->x == 0 ) continue;
						if ( y==-1 ) if ( pit->y == 0 ) continue;
						if ( z==-1 ) if ( pit->z == 0 ) continue;
						if ( x==1 ) if ( pit->x == xdim_-1 ) continue;
						if ( y==1 ) if ( pit->y == ydim_-1 ) continue;
						if ( z==1 ) if ( pit->z == zdim_-1 ) continue;
						if ( grid_[pit->x+x][pit->y+y][pit->z+z]==T_SURFACE ) {
							cluster.target_=true;
						} else if ( grid_[pit->x+x][pit->y+y][pit->z+z]==ST_SURFACE ) {
							cluster.subtarget_=true;
						}
						if ( grid_[pit->x+x][pit->y+y][pit->z+z]==EMPTY || grid_[pit->x+x][pit->y+y][pit->z+z]==PO_SURF || grid_[pit->x+x][pit->y+y][pit->z+z]==TP_SURF ) {
							cluster.solventExposed_=true;
						}
						if ( cluster.target_ && (cluster.subtarget_ || numTargets_ == 1) && cluster.solventExposed_ ) {
							x=2;y=2;z=2;
							pit=cluster.points_.end();
							--pit;
						}
					}
				}
			}
		}
	}

	//Change target clusters to target types on grid
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( !cluster.isTarget(numTargets_) ) {
			for ( auto & point : cluster.points_ ) {
				if ( grid_[point.x][point.y][point.z]==TP_POCKET ) grid_[point.x][point.y][point.z]=POCKET;
				if ( grid_[point.x][point.y][point.z]==TP_SURF ) grid_[point.x][point.y][point.z]=PO_SURF;
				if ( grid_[point.x][point.y][point.z]==TP_BURIED ) grid_[point.x][point.y][point.z]=PO_BURIED;
				if ( grid_[point.x][point.y][point.z]==TP_EDGE ) grid_[point.x][point.y][point.z]=PO_EDGE;
			}
			continue;
		}

		for ( auto & point : cluster.points_ ) {
			if ( grid_[point.x][point.y][point.z]==POCKET ) grid_[point.x][point.y][point.z]=TP_POCKET;
			if ( grid_[point.x][point.y][point.z]==PO_SURF ) grid_[point.x][point.y][point.z]=TP_SURF;
			if ( grid_[point.x][point.y][point.z]==PO_BURIED ) grid_[point.x][point.y][point.z]=TP_BURIED;
			if ( grid_[point.x][point.y][point.z]==PO_EDGE ) grid_[point.x][point.y][point.z]=TP_EDGE;
		}
	}
}

void PocketGrid::findExemplars(core::pose::Pose const & inPose, Size const total_residues){
	//ideal H-bond distance
	//core::Real const opt_distance( 2.75 );
	//core::Real const distance( 3.0 );
	core::Real const opt_distance( 2.75 );
	//core::Real const distance( 3.5 );
	using namespace core::chemical;
	using namespace core::kinematics;

	//this sets up a temporary grid I can tick off pocket points from for minimal pockets
	std::vector < std::vector < std::vector <PtType> > > tmpgrid_;
	tmpgrid_.resize(xdim_);
	for ( core::Size tx=0; tx<xdim_; tx++ ) {
		tmpgrid_[tx].resize(ydim_);
		for ( core::Size ty=0; ty<ydim_; ty++ ) {
			tmpgrid_[tx][ty].resize(zdim_, EMPTY);
		}
	}

	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			tmpgrid_[point.x][point.y][point.z] = grid_[point.x][point.y][point.z];
		}
	}

	c_clusters_.clear();

	// We only want to find donors/acceptors that are solvent accessible
	basic::MetricValue< utility::vector1< core::Real > > resisasa;
	inPose.metric( "sasa", "residue_sasa", resisasa );
	basic::MetricValue< core::id::AtomID_Map< core::Real> > atomsasa;
	inPose.metric( "sasa", "atom_sasa", atomsasa );
	core::id::AtomID_Map< core::Real > atom_sasas = atomsasa.value();


	//Go through residues near the the pocket and look for hydrogen donors and acceptors
	for ( Size j = 1, resnum = total_residues; j <= resnum; ++j ) {
		core::conformation::Residue const & rsd( inPose.conformation().residue(j) );

		// Skip buried residues
		//   if (resisasa.value()[j] < 10) {
		//    std::cout<<rsd.seqpos()-3<<" "<<rsd.name()<<" SASA "<<resisasa.value()[j]<<" being ignored"<<std::endl;
		//    continue;
		//   }

		//int target=0;
		//core::Size total_atoms(0);
		using namespace basic::options;
		/*if (option[ OptionKeys::fingerprint::include_hydrogens ]()){
		total_atoms = rsd.natoms();
		} else {
		total_atoms = rsd.nheavyatoms();
		}*/

		//fill in points that are ideal for a hydrogen acceptor with an O
		for ( auto
				hnum  = rsd.Hpos_polar().begin(),
				hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {

			Size const hatm( *hnum );

			// Skip buried residues
			int offset = 9;
			if ( rsd.seqpos() < 72 ) offset=3;
			offset=0;
			if ( atom_sasas(j, hatm) < 0.1 && atom_sasas(j, rsd.atom_base(hatm)) < 0.1 ) {
				std::cout<<rsd.seqpos()+offset<<" Donor "<<rsd.name()<<" "<<rsd.atom_name(rsd.atom_base(hatm))<<" H SASA "<<atom_sasas(j, hatm)<<" Base SASA "<<atom_sasas(j, rsd.atom_base(hatm))<<" being ignored"<<std::endl;
				continue;
			}
			std::cout<<rsd.seqpos()+offset<<" Donor "<<rsd.name()<<" "<<rsd.atom_name(rsd.atom_base(hatm))<<" H SASA "<<atom_sasas(j, hatm)<<" Base SASA "<<atom_sasas(j, rsd.atom_base(hatm))<<std::endl;

			numeric::xyzVector<core::Real> const & hatm_xyz( rsd.xyz( hatm ) );
			numeric::xyzVector<core::Real> const & datm_xyz( rsd.xyz( rsd.atom_base( hatm ) ) );

			for ( int step = 0; step<4; step++ ) {
				numeric::xyzVector<core::Real> const o1(datm_xyz + (opt_distance + step * stepSize_) * ( hatm_xyz - datm_xyz ).normalized() );
				numeric::xyzVector<core::Real> const ro1(datm_xyz + opt_distance * ( hatm_xyz - datm_xyz ).normalized());
				numeric::xyzVector<core::Real> rpoint = rotatePoint(o1.x(),o1.y(),o1.z());
				numeric::xyzVector<core::Real> rrpoint = rotatePoint(ro1.x(),ro1.y(),ro1.z());
				if ( step==0 ) std::cout<<"Optimal acceptor at "<<ro1.x()<<" "<<ro1.y()<<" "<<ro1.z()<<std::endl;
				core::Size x,y,z;
				if ( rpoint.x()<xcorn_ || rpoint.y()<ycorn_ || rpoint.z()<zcorn_ ) continue;
				x=(int)floor((rpoint.x()-xcorn_)/stepSize_+.5);
				y=(int)floor((rpoint.y()-ycorn_)/stepSize_+.5);
				z=(int)floor((rpoint.z()-zcorn_)/stepSize_+.5);
				if ( x<2 || x>xdim_ -3 || y<2 || y>ydim_ -3 || z<2 || z>zdim_ -3 ) continue;
				//    if (!isProteinPoint(x,y,z)) {
				//    if (isTargetPocketPoint(x,y,z)) {
				if ( true ) {

					if ( !isTargetPocketPoint(x,y,z) && !isTargetPocketPoint(x+2,y,z) && !isTargetPocketPoint(x-2,y,z) && !isTargetPocketPoint(x,y+2,z) && !isTargetPocketPoint(x,y-2,z) &&
							!isTargetPocketPoint(x,y,z+2) &&  !isTargetPocketPoint(x,y,z-2) ) continue;
					/*      std::cout<<x<<" "<<y<<" "<<z<<" "<<grid_[x][y][z]<<" "<<rrpoint.x()<<" "<<rrpoint.y()<<" "<<rrpoint.z()<<" "<<std::endl;
					std::cout<<x+2<<" "<<y<<" "<<z<<" "<<grid_[x+2][y][z]<<std::endl;
					std::cout<<x-2<<" "<<y<<" "<<z<<" "<<grid_[x-2][y][z]<<std::endl;
					std::cout<<x<<" "<<y+2<<" "<<z<<" "<<grid_[x][y+2][z]<<std::endl;
					std::cout<<x<<" "<<y-2<<" "<<z<<" "<<grid_[x][y-2][z]<<std::endl;
					std::cout<<x<<" "<<y<<" "<<z+2<<" "<<grid_[x][y][z+2]<<std::endl;
					std::cout<<x<<" "<<y<<" "<<z-2<<" "<<grid_[x][y][z-2]<<std::endl;
					for (int i = -1; i <2; i++){
					for (int j = -1; j <2; j++){
					for (int k = -1; k <2; k++){
					if (i==0 && j==0 && k==0) continue;
					std::cout<<x+i<<" "<<y+j<<" "<<z+k<<" "<<grid_[x+i][y+j][z+k]<<std::endl;
					}
					}
					}
					std::cout<<std::endl;
					*/

					if ( !isProteinPoint(x,y,z) && !isProteinPoint(x+2,y,z) && !isProteinPoint(x-2,y,z) && !isProteinPoint(x,y+2,z) && !isProteinPoint(x,y-2,z) && !isProteinPoint(x,y,z+2) &&
							!isProteinPoint(x,y,z-2) ) {
						//std::cout<<x<<" "<<y<<" "<<z<<std::endl;
						bool pockpt=true;
						//bool tpt=false ;
						bool tpt=true;
						int buriedness=0;
						for ( int i = -1; i <2; i++ ) {
							for ( int j = -1; j <2; j++ ) {
								for ( int k = -1; k <2; k++ ) {
									//std::cout<<x+i<<" "<<y+j<<" "<<z+k<<" "<<grid_[x+i][y+j][z+k]<<std::endl;
									if ( isDeepTargetPocketPoint(x+i,y+j,z+k) ) {
										tpt=true;
										buriedness++;
									}
									if ( i==0 && j==0 && k==0 ) continue;
									if ( isProteinPoint(x+i,y+j,z+k) ) {
										pockpt=false;
									}
								}
							}
						}
						//std::cout<<std::endl;
						if ( pockpt && tpt ) {
							/*        for (int i = -1; i <2; i++){
							for (int j = -1; j <2; j++){
							for (int k = -1; k <2; k++){
							tmpgrid_[x+i][y+j][z+k] = EMPTY;
							}
							}
							}
							//tmpgrid_[x-2][y][z] = EMPTY;
							//tmpgrid_[x][y-2][z] = EMPTY;
							//tmpgrid_[x][y][z-2] = EMPTY;
							*/
							core::Real offset = (step - 1)*stepSize_;
							if ( offset < 0 ) offset=0;
							numeric::xyzVector<core::Real> const offro1(datm_xyz + (opt_distance+offset) * ( hatm_xyz - datm_xyz ).normalized());
							numeric::xyzVector<core::Real> orrpoint = rotatePoint(offro1.x(),offro1.y(),offro1.z());
							std::cout<<rrpoint.x()<<" "<<rrpoint.y()<<" "<<rrpoint.z()<<", "<<orrpoint.x()<<" "<<orrpoint.y()<<" "<<orrpoint.z()<<" "<<x<<" "<<y<<" "<<z<<", Buriedness: "<<buriedness<<"/27"<<std::endl;

							//Treat more than half buried H-bonders the same as fully buried hydrophobics for special cases
							if ( buriedness < 14 || !filterExemplars_ ) {
								c_clusters_.add(x,y,z, "Ne", stepSize_, orrpoint.x(), orrpoint.y(), orrpoint.z());
							} else c_clusters_.add(x,y,z, "BNe", stepSize_, orrpoint.x(), orrpoint.y(), orrpoint.z());
							break;
						}
					}
				}
			}
		}

		//fill in points that are ideal for a hydrogen donor with an N
		for ( auto
				anum  = rsd.accpt_pos().begin(),
				anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
			Size const aatm( *anum );
			// Skip buried residues
			int offset = 9;
			if ( rsd.seqpos() < 72 ) offset=3;
			offset=0;
			if ( atom_sasas(j, aatm) < 0.1 ) {
				std::cout<<rsd.seqpos()+offset<<" Acceptor "<<rsd.name()<<" "<<rsd.atom_name(aatm)<<" SASA "<<atom_sasas(j, aatm)<<" being ignored"<<std::endl;
				continue;
			}
			std::cout<<rsd.seqpos()+offset<<" Acceptor "<<rsd.name()<<" "<<rsd.atom_name(aatm)<<" SASA "<<atom_sasas(j, aatm)<<std::endl;

			numeric::xyzVector<core::Real> const & aatm_xyz( rsd.xyz( aatm ) );
			numeric::xyzVector<core::Real> aatm_base_xyz( rsd.xyz( rsd.atom_base( aatm ) ) );
			numeric::xyzVector<core::Real> const & aatm_base2_xyz( rsd.xyz( rsd.abase2( aatm ) ) );
			Hybridization const & hybrid( rsd.atom_type(aatm).hybridization() );

			core::Real theta(0.0); // step_size(0.0);
			utility::vector1< core::Real > phi_list, phi_steps;
			phi_steps.push_back(  0 );
			switch( hybrid ) {
			case SP2_HYBRID :
				theta = 180.0 - 120.0;
				//step_size = 15.0;
				phi_list.push_back(   0.0 );
				phi_list.push_back( 180.0 );
				break;
			case SP3_HYBRID :
				theta = 180.0 - 109.0;
				//step_size = 10.0;
				phi_list.push_back( 120.0 );
				phi_list.push_back( 240.0 );
				break;
			case RING_HYBRID :
				{
				numeric::xyzVector<core::Real> const & avg_base_xyz (0.5 * ( aatm_base_xyz + aatm_base2_xyz ));
				aatm_base_xyz(1)=avg_base_xyz(1);
				aatm_base_xyz(2)=avg_base_xyz(2);
				aatm_base_xyz(3)=avg_base_xyz(3);
				theta = 0.0;
				phi_steps.clear();
				phi_steps.push_back( 0.0 );
				phi_list.push_back( 0.0 ); // doesnt matter
				//step_size = 0.0; // doesnt matter
				break;
			}
			default :
				std::cerr << "Bad hybridization type for acceptor " << hybrid << '\n';
				exit(1000);
			}
			Stub stub( aatm_xyz, aatm_base_xyz, aatm_base2_xyz );
			for ( Size i=1; i<= phi_list.size(); ++i ) {
				for ( int step = 0; step<4; step++ ) {
					numeric::xyzVector<core::Real> const o1(stub.spherical( numeric::conversions::radians( phi_list[i]), numeric::conversions::radians( theta ), opt_distance + step*stepSize_));
					numeric::xyzVector<core::Real> rpoint = rotatePoint(o1.x(),o1.y(),o1.z());
					numeric::xyzVector<core::Real> const ro1(stub.spherical( numeric::conversions::radians( phi_list[i]), numeric::conversions::radians( theta ), opt_distance));
					numeric::xyzVector<core::Real> rrpoint = rotatePoint(ro1.x(),ro1.y(),ro1.z());
					if ( step==0 ) std::cout<<"Optimal donor at "<<ro1.x()<<" "<<ro1.y()<<" "<<ro1.z()<<std::endl;
					core::Size x,y,z;
					if ( rpoint.x()<xcorn_ || rpoint.y()<ycorn_ || rpoint.z()<zcorn_ ) continue;
					x=(int)floor((rpoint.x()-xcorn_)/stepSize_+.5);
					y=(int)floor((rpoint.y()-ycorn_)/stepSize_+.5);
					z=(int)floor((rpoint.z()-zcorn_)/stepSize_+.5);
					if ( x<2 || x>xdim_ -3 || y<2 || y>ydim_ -3 || z<2 || z>zdim_ -3 ) continue;
					//    if (!isProteinPoint(x,y,z)) {
					if ( rrpoint.x() < 5 &&  rrpoint.x()>3 && rrpoint.y()>18 ) {
						std::cout<<x<<" "<<y<<" "<<z<<" "<<grid_[x][y][z]<<" "<<rrpoint.x()<<" "<<rrpoint.y()<<" "<<rrpoint.z()<<" "<<std::endl;
					}
					//     if (isTargetPocketPoint(x,y,z)) {
					if ( true ) {
						if ( !isTargetPocketPoint(x,y,z) && !isTargetPocketPoint(x+2,y,z) && !isTargetPocketPoint(x-2,y,z) && !isTargetPocketPoint(x,y+2,z) && !isTargetPocketPoint(x,y-2,z) &&
								!isTargetPocketPoint(x,y,z+2) && !isTargetPocketPoint(x,y,z-2) ) continue;

						if ( rrpoint.x() < 5 &&  rrpoint.x()>2 && rrpoint.y()>18 ) {
							std::cout<<x<<" "<<y<<" "<<z<<" "<<grid_[x][y][z]<<" "<<rrpoint.x()<<" "<<rrpoint.y()<<" "<<rrpoint.z()<<" "<<std::endl;
							std::cout<<x+2<<" "<<y<<" "<<z<<" "<<grid_[x+2][y][z]<<std::endl;
							std::cout<<x-2<<" "<<y<<" "<<z<<" "<<grid_[x-2][y][z]<<std::endl;
							std::cout<<x<<" "<<y+2<<" "<<z<<" "<<grid_[x][y+2][z]<<std::endl;
							std::cout<<x<<" "<<y-2<<" "<<z<<" "<<grid_[x][y-2][z]<<std::endl;
							std::cout<<x<<" "<<y<<" "<<z+2<<" "<<grid_[x][y][z+2]<<std::endl;
							std::cout<<x<<" "<<y<<" "<<z-2<<" "<<grid_[x][y][z-2]<<std::endl;
							for ( int i = -1; i <2; i++ ) {
								for ( int j = -1; j <2; j++ ) {
									for ( int k = -1; k <2; k++ ) {
										if ( i==0 && j==0 && k==0 ) continue;
										std::cout<<x+i<<" "<<y+j<<" "<<z+k<<" "<<grid_[x+i][y+j][z+k]<<std::endl;
									}
								}
							}
							std::cout<<std::endl;

						}

						if ( !isProteinPoint(x,y,z) && !isProteinPoint(x+2,y,z) && !isProteinPoint(x-2,y,z) && !isProteinPoint(x,y+2,z) && !isProteinPoint(x,y-2,z) && !isProteinPoint(x,y,z+2) &&
								!isProteinPoint(x,y,z-2) ) {
							//std::cout<<x<<" "<<y<<" "<<z<<std::endl;
							bool pockpt=true;
							//bool tpt=false;
							bool tpt=true;
							int buriedness=0;
							for ( int c = -1; c <2; c++ ) {
								for ( int j = -1; j <2; j++ ) {
									for ( int k = -1; k <2; k++ ) {
										//std::cout<<x+i<<" "<<y+j<<" "<<z+k<<" "<<grid_[x+i][y+j][z+k]<<std::endl;
										if ( isDeepTargetPocketPoint(x+c,y+j,z+k) ) {
											tpt=true;
											buriedness++;
										}
										if ( c==0 && j==0 && k==0 ) continue;
										if ( isProteinPoint(x+c,y+j,z+k) ) {
											pockpt=false;
										}
									}
								}
							}
							//std::cout<<std::endl;
							if ( pockpt && tpt ) {
								/*         for (int c = -1; c <2; c++){
								for (int j = -1; j <2; j++){
								for (int k = -1; k <2; k++){
								tmpgrid_[x+c][y+j][z+k] = EMPTY;
								}
								}
								}
								//tmpgrid_[x-2][y][z] = EMPTY;
								//tmpgrid_[x][y-2][z] = EMPTY;
								//tmpgrid_[x][y][z-2] = EMPTY;
								*/
								core::Real offset = (step - 1)*stepSize_;
								if ( offset < 0 ) offset=0;
								numeric::xyzVector<core::Real> const offro1(stub.spherical( numeric::conversions::radians( phi_list[i]), numeric::conversions::radians( theta ), opt_distance+offset));
								numeric::xyzVector<core::Real> orrpoint = rotatePoint(offro1.x(),offro1.y(),offro1.z());
								std::cout<<rrpoint.x()<<" "<<rrpoint.y()<<" "<<rrpoint.z()<<" "<<x<<" "<<y<<" "<<z<<", Buriedness: "<<buriedness<<"/27"<<std::endl;

								//Treat more than half buried H-bonders the same as fully buried hydrophobics for special cases
								if ( buriedness < 14 || !filterExemplars_ ) {
									c_clusters_.add(x,y,z, "Be", stepSize_, orrpoint.x(), orrpoint.y(), orrpoint.z());
								} else c_clusters_.add(x,y,z, "BBe", stepSize_, orrpoint.x(), orrpoint.y(), orrpoint.z());
								break;
							}
						}
					}
				}
			}

		}


	}

	std::cout<<std::endl<<"After Donors/Acceptors"<<std::endl;
	for ( auto & cluster : c_clusters_.clusters_ ) {
		for ( auto & point : cluster.points_ ) {
			if ( point.atom_type.compare("C") == 0 ) {
				std::cout<<point.atom_type<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.x*stepSize_+xcorn_<<std::setw(8)<<point.y*stepSize_+ycorn_<<std::setw(8)<<point.z*stepSize_+zcorn_<<std::endl;
			} else {
				std::cout<<point.atom_type<<" "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.absX<<std::setw(8)<<point.absY<<std::setw(8)<<point.absZ<<std::endl;
			}
			for ( int c = -1; c <2; c++ ) {
				for ( int j = -1; j <2; j++ ) {
					for ( int k = -1; k <2; k++ ) {
						tmpgrid_[point.x+c][point.y+j][point.z+k] = EMPTY;
					}
				}
			}


		}
	}


	//fill in remaining potential exemplar points with C.
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			if ( tmpgrid_[point.x][point.y][point.z] == EMPTY ) continue;
			if ( point.x < 2 || point.y < 2|| point.z < 2 || point.x > xdim_-3 || point.y > ydim_-3|| point.z > zdim_-3 ) continue;
			if ( (tmpgrid_[point.x+2][point.y][point.z] != EMPTY && tmpgrid_[point.x-2][point.y][point.z] != EMPTY) &&
					(tmpgrid_[point.x][point.y+2][point.z] != EMPTY && tmpgrid_[point.x][point.y-2][point.z] != EMPTY) &&
					(tmpgrid_[point.x][point.y][point.z+2] != EMPTY && tmpgrid_[point.x][point.y][point.z-2] != EMPTY) ) {
				bool pockpt=true;
				for ( int i = -1; i <2; i++ ) {
					for ( int j = -1; j <2; j++ ) {
						for ( int k = -1; k <2; k++ ) {
							if ( i==0 && j==0 && k==0 ) continue;
							if ( tmpgrid_[point.x+i][point.y+j][point.z+k] == EMPTY ) {
								pockpt=false;
							}
						}
					}
				}
				if ( pockpt ) {
					for ( int i = -1; i <2; i++ ) {
						for ( int j = -1; j <2; j++ ) {
							for ( int k = -1; k <2; k++ ) {
								tmpgrid_[point.x+i][point.y+j][point.z+k] = EMPTY;
							}
						}
					}
					//       tmpgrid_[pit->x-2][pit->y][pit->z] = EMPTY;
					//       tmpgrid_[pit->x][pit->y-2][pit->z] = EMPTY;
					//       tmpgrid_[pit->x][pit->y][pit->z-2] = EMPTY;
					c_clusters_.add(point.x,point.y,point.z, "C", stepSize_);
				} else continue;
			} else continue;
		}
	}

	std::cout<<std::endl<<"After shape"<<std::endl;
	for ( auto & cluster : c_clusters_.clusters_ ) {
		for ( auto & point : cluster.points_ ) {
			if ( point.atom_type.compare("C") == 0 ) {
				std::cout<<point.atom_type<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.x*stepSize_+xcorn_<<std::setw(8)<<point.y*stepSize_+ycorn_<<std::setw(8)<<point.z*stepSize_+zcorn_<<std::endl;
			} else {
				std::cout<<point.atom_type<<" "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.absX<<std::setw(8)<<point.absY<<std::setw(8)<<point.absZ<<std::endl;
			}
		}
	}

	c_clusters_.findClusters();


	std::cout<<std::endl<<"After clustering"<<std::endl;
	for ( auto & cluster : c_clusters_.clusters_ ) {
		for ( auto & point : cluster.points_ ) {
			if ( point.atom_type.compare("C") == 0 ) {
				std::cout<<point.atom_type<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.x*stepSize_+xcorn_<<std::setw(8)<<point.y*stepSize_+ycorn_<<std::setw(8)<<point.z*stepSize_+zcorn_<<std::endl;
			} else {
				std::cout<<point.atom_type<<" "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.absX<<std::setw(8)<<point.absY<<std::setw(8)<<point.absZ<<std::endl;
			}
		}
	}

	//int test=1;
	for ( auto & cluster : c_clusters_.clusters_ ) {
		for ( auto pit=cluster.points_.begin(); pit != cluster.points_.end(); ++pit ) {
			for ( int x = -1* (int) (5/stepSize_)-1; x <= (int) (5/stepSize_)+1; x++ ) {
				for ( int y = -1* (int) (5/stepSize_)-1; y <= (int) (5/stepSize_)+1; y++ ) {
					for ( int z = -1* (int) (5/stepSize_)-1; z <= (int) (5/stepSize_)+1; z++ ) {
						if ( (int) pit->x + x <0 ) continue;
						if ( (int) pit->y + y <0 ) continue;
						if ( (int) pit->z + z <0 ) continue;
						if ( (int) pit->x + x > (int) xdim_-1 ) continue;
						if ( (int) pit->y + y > (int) ydim_-1 ) continue;
						if ( (int) pit->z + z > (int) zdim_-1 ) continue;
						if ( sqrt(pow((double)x,2)+pow((double)y,2)+pow((double)z,2)) > 5./stepSize_ ) continue;
						if ( grid_[pit->x+x][pit->y+y][pit->z+z]==T_SURFACE ) {
							cluster.target_=true;
						} else if ( grid_[pit->x+x][pit->y+y][pit->z+z]==ST_SURFACE ) {
							cluster.subtarget_=true;
						}
						if ( grid_[pit->x+x][pit->y+y][pit->z+z]==EMPTY || grid_[pit->x+x][pit->y+y][pit->z+z]==PO_SURF || grid_[pit->x+x][pit->y+y][pit->z+z]==TP_SURF ) {
							cluster.solventExposed_=true;
						}
						if ( cluster.target_ && (cluster.subtarget_ || numTargets_ == 1) && cluster.solventExposed_ ) {
							x=200;y=200;z=200;
							pit=cluster.points_.end();
							--pit;
						}
					}
				}
			}

		}
	}
	std::cout<<std::endl<<"End of findExemplars"<<std::endl;
	for ( auto & cluster : c_clusters_.clusters_ ) {
		for ( auto & point : cluster.points_ ) {
			if ( point.atom_type.compare("C") == 0 ) {
				std::cout<<point.atom_type<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.x*stepSize_+xcorn_<<std::setw(8)<<point.y*stepSize_+ycorn_<<std::setw(8)<<point.z*stepSize_+zcorn_<<std::endl;
			} else {
				std::cout<<point.atom_type<<" "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.absX<<std::setw(8)<<point.absY<<std::setw(8)<<point.absZ<<std::endl;
			}
		}
	}

}

void PocketGrid::findClustersByExemplars(){
	//clear clusters
	clusters_.clear();
	//go through c clusters
	for ( auto & cluster : c_clusters_.clusters_ ) {
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			for ( int x = -1* (int) (3/stepSize_)-1; x <= (int) (3/stepSize_)+1; x++ ) {
				for ( int y = -1* (int) (3/stepSize_)-1; y <= (int) (3/stepSize_)+1; y++ ) {
					for ( int z = -1* (int) (3/stepSize_)-1; z <= (int) (3/stepSize_)+1; z++ ) {
						if ( (int) point.x + x <0 ) continue;
						if ( (int) point.y + y <0 ) continue;
						if ( (int) point.z + z <0 ) continue;
						if ( (int) point.x + x > (int) xdim_-1 ) continue;
						if ( (int) point.y + y > (int) ydim_-1 ) continue;
						if ( (int) point.z + z > (int) zdim_-1 ) continue;
						//add each c cluster point
						if ( grid_[point.x+x][point.y+y][point.z+z]==TP_BURIED || grid_[point.x+x][point.y+y][point.z+z]==TP_POCKET || grid_[point.x+x][point.y+y][point.z+z]==TP_EDGE ) {
							c_grid_[point.x+x][point.y+y][point.z+z] = TP_BURIED;
						}
					}
				}
			}
		}
	}

	//cluster those
	clusterCPockets();

}

void PocketGrid::linkExemplarsThroughSolvent(){

	//go through c clusters
	for ( auto cit=c_clusters_.clusters_.begin(); cit != c_clusters_.clusters_.end(); ++cit ) {
		if ( !cit->isTarget(numTargets_) ) continue;
		for ( auto cit2=cit; cit2 != c_clusters_.clusters_.end(); ++cit2 ) {
			if ( cit == cit2 ) continue;
			if ( !cit2->isTarget(numTargets_) ) continue;
			core::Real closest=101.;
			std::list<CCluster::Cxyz>::iterator p1;
			std::list<CCluster::Cxyz>::iterator p2;
			for ( auto pit=cit->points_.begin(); pit != cit->points_.end(); ++pit ) {
				for ( auto pit2=cit2->points_.begin(); pit2 != cit2->points_.end(); ++pit2 ) {
					core::Real dist=pow(pit->x-pit2->x,2)+pow(pit->z-pit2->z,2)+pow(pit->z-pit2->z,2);
					if ( dist<closest ) {
						closest=dist;
						p1=pit;
						p2=pit2;
					}
				}
			}
			if ( closest<101. ) {
				int txdim = (int) ceil(p1->x-p2->x/stepSize_)+(int)ceil(3./stepSize_);
				core::Real txcorn=(core::Real)std::floor(std::min(p1->x,p2->x)-1.0);
				int txoff = (int)floor((txcorn-xcorn_)/stepSize_ + 0.5);
				int tydim = (int) ceil(p1->y-p2->y/stepSize_)+(int)ceil(3./stepSize_);
				core::Real tycorn=(core::Real)std::floor(std::min(p1->y,p2->y)-1.0);
				int tyoff = (int)floor((tycorn-ycorn_)/stepSize_ + 0.5);
				int tzdim = (int) ceil(p1->z-p2->z/stepSize_)+(int)ceil(3./stepSize_);
				core::Real tzcorn=(core::Real)std::floor(std::min(p1->z,p2->z)-1.0);
				int tzoff = (int)floor((tzcorn-zcorn_)/stepSize_ + 0.5);
				std::vector < std::vector < std::vector <PtType> > > tmpgrid_;
				tmpgrid_.resize(txdim);
				for ( int tx=0; tx<txdim; tx++ ) {
					tmpgrid_[tx].resize(tydim);
					for ( int ty=0; ty<tydim; ty++ ) {
						tmpgrid_[tx][ty].resize(tzdim, PROTEIN);
					}
				}

				//fill in smaller grid
				for ( int tx=0; tx<txdim; tx++ ) {
					for ( int ty=0; ty<tydim; ty++ ) {
						for ( int tz=0; tz<tzdim; tz++ ) {
							if ( grid_[tx+txoff][ty+tyoff][tz+tzoff] == TP_BURIED ) tmpgrid_[tx][ty][tz]=TP_POCKET;
							if ( grid_[tx+txoff][ty+tyoff][tz+tzoff] == TP_SURF ) tmpgrid_[tx][ty][tz]=TP_POCKET;
							if ( grid_[tx+txoff][ty+tyoff][tz+tzoff] == TP_EDGE ) tmpgrid_[tx][ty][tz]=TP_POCKET;
							if ( grid_[tx+txoff][ty+tyoff][tz+tzoff] == TP_POCKET ) tmpgrid_[tx][ty][tz]=TP_POCKET;
							if ( grid_[tx+txoff][ty+tyoff][tz+tzoff] == EMPTY ) tmpgrid_[tx][ty][tz]=EMPTY;
						}
					}
				}

				//take away already exemplar points
				int x,y,z;
				x=(int)floor((p1->x - txcorn)/stepSize_ + 0.5);
				y=(int)floor((p1->y - tycorn)/stepSize_ + 0.5);
				z=(int)floor((p1->z - tzcorn)/stepSize_ + 0.5);
				for ( int tx=x-1; tx<x+1; tx++ ) {
					if ( tx<0 ) continue;
					for ( int ty=y-1; ty<y+1; ty++ ) {
						if ( ty<0 ) continue;
						for ( int tz=z-1; tz<z+1; tz++ ) {
							if ( tz<0 ) continue;
							tmpgrid_[tx][ty][tz]=POCKET;
						}
					}
				}
				x=(int)floor((p2->x - txcorn)/stepSize_ + 0.5);
				y=(int)floor((p2->y - tycorn)/stepSize_ + 0.5);
				z=(int)floor((p2->z - tzcorn)/stepSize_ + 0.5);
				for ( int tx=x-1; tx<x+1; tx++ ) {
					if ( tx<0 ) continue;
					for ( int ty=y-1; ty<y+1; ty++ ) {
						if ( ty<0 ) continue;
						for ( int tz=z-1; tz<z+1; tz++ ) {
							if ( tz<0 ) continue;
							tmpgrid_[tx][ty][tz]=POCKET;
						}
					}
				}


				for ( auto & cluster : clusters_.clusters_ ) {
					if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
					if ( !cluster.isTarget(numTargets_) ) continue;
					for ( auto & point : cluster.points_ ) {
						tmpgrid_[point.x][point.y][point.z] = grid_[point.x][point.y][point.z];
					}
				}

			}
		}
	}

	//cluster those
	clusterCPockets();

}


core::Real PocketGrid::largestTargetPocketVolume() {
	return clusters_.getLargestClusterSize( stepSize_, minPockSize_, numTargets_, ignoreBuriedPockets_, ignoreExposedPockets_);
}

core::Real PocketGrid::netTargetPocketVolume() {
	if ( maxPockSize_ ) {
		return std::min(maxPockSize_, clusters_.getNetClusterSize( stepSize_ , minPockSize_, numTargets_, ignoreBuriedPockets_, ignoreExposedPockets_));
	} else return clusters_.getNetClusterSize( stepSize_ , minPockSize_, numTargets_, ignoreBuriedPockets_, ignoreExposedPockets_);
}


bool PocketGrid::isTooSmall() const {
	core::Size x,y,z;
	if ( xdim_==0||ydim_==0||zdim_==0 ) return false;
	if ( grid_[0][0][0] == TP_POCKET ) return true;
	if ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_POCKET ) return true;
	if ( grid_[0][0][0] == TP_SURF ) return true;
	if ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_SURF ) return true;
	if ( grid_[0][0][0] == TP_BURIED ) return true;
	if ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_BURIED ) return true;
	if ( grid_[0][0][0] == TP_EDGE ) return true;
	if ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_EDGE ) return true;

	for ( x=0, y=0, z=0; x<xdim_ || y<ydim_ || z<zdim_; x++,y++,z++ ) {
		if ( x==xdim_ ) x--;
		if ( y==ydim_ ) y--;
		if ( z==zdim_ ) z--;
		if ( grid_[0][y][z] == TP_POCKET ) return true;
		if ( grid_[x][0][z] == TP_POCKET ) return true;
		if ( grid_[x][y][0] == TP_POCKET ) return true;
		if ( grid_[xdim_-1][y][z] == TP_POCKET ) return true;
		if ( grid_[x][ydim_-1][z] == TP_POCKET ) return true;
		if ( grid_[x][y][zdim_-1] == TP_POCKET ) return true;
		if ( grid_[0][y][z] == TP_SURF ) return true;
		if ( grid_[x][0][z] == TP_SURF ) return true;
		if ( grid_[x][y][0] == TP_SURF ) return true;
		if ( grid_[xdim_-1][y][z] == TP_SURF ) return true;
		if ( grid_[x][ydim_-1][z] == TP_SURF ) return true;
		if ( grid_[x][y][zdim_-1] == TP_SURF ) return true;
		if ( grid_[0][y][z] == TP_BURIED ) return true;
		if ( grid_[x][0][z] == TP_BURIED ) return true;
		if ( grid_[x][y][0] == TP_BURIED ) return true;
		if ( grid_[xdim_-1][y][z] == TP_BURIED ) return true;
		if ( grid_[x][ydim_-1][z] == TP_BURIED ) return true;
		if ( grid_[x][y][zdim_-1] == TP_BURIED ) return true;
		if ( grid_[0][y][z] == TP_EDGE ) return true;
		if ( grid_[x][0][z] == TP_EDGE ) return true;
		if ( grid_[x][y][0] == TP_EDGE ) return true;
		if ( grid_[xdim_-1][y][z] == TP_EDGE ) return true;
		if ( grid_[x][ydim_-1][z] == TP_EDGE ) return true;
		if ( grid_[x][y][zdim_-1] == TP_EDGE ) return true;
	}


	return false;
}

core::Vector PocketGrid::whatIsTooSmall() const {
	core::Vector dims;
	dims(1)=0;
	dims(2)=0;
	dims(3)=0;
	core::Size x,y,z;
	if ( xdim_==0||ydim_==0||zdim_==0 ) return dims;
	if ( ( grid_[0][0][0] == TP_POCKET) || ( grid_[0][0][0] == TP_SURF) || ( grid_[0][0][0] == TP_BURIED) || ( grid_[0][0][0] == TP_EDGE) || ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_POCKET) || ( grid_[xdim_-1][ydim_-1][zdim_-1]) || ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_BURIED) || ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_EDGE) ) {
		dims(1)=1;
		dims(2)=1;
		dims(3)=1;
		return dims;
	}

	for ( x=0, y=0, z=0; x<xdim_ || y<ydim_ || z<zdim_; x++,y++,z++ ) {
		if ( x==xdim_ ) x--;
		if ( y==ydim_ ) y--;
		if ( z==zdim_ ) z--;
		if ( ( grid_[0][y][z] == TP_POCKET) || ( grid_[0][y][z] == TP_BURIED) || ( grid_[0][y][z] == TP_EDGE) || ( grid_[xdim_-1][y][z] == TP_POCKET) || ( grid_[xdim_-1][y][z] == TP_SURF) || ( grid_[xdim_-1][y][z] == TP_BURIED) || ( grid_[xdim_-1][y][z] == TP_EDGE) ) {
			dims(1)=1;
		}
		if ( ( grid_[x][0][z] == TP_POCKET) || ( grid_[x][0][z] == TP_SURF) || ( grid_[x][0][z] == TP_BURIED) || ( grid_[x][0][z] == TP_EDGE) || ( grid_[x][ydim_-1][z] == TP_POCKET) || ( grid_[x][ydim_-1][z] == TP_SURF) ||  ( grid_[x][ydim_-1][z] == TP_BURIED) || ( grid_[x][ydim_-1][z] == TP_EDGE) ) {
			dims(2)=1;
		}
		if ( ( grid_[x][y][0] == TP_POCKET) || ( grid_[x][y][0] == TP_SURF) || ( grid_[x][y][0] == TP_BURIED) || ( grid_[x][y][0] == TP_EDGE) || ( grid_[x][y][zdim_-1] == TP_POCKET) || ( grid_[x][y][zdim_-1] == TP_SURF) || ( grid_[x][y][zdim_-1] == TP_BURIED) || ( grid_[x][y][zdim_-1] == TP_EDGE) ) {
			dims(3)=1;
		}
	}


	return dims;
}

utility::vector1<core::Real> PocketGrid::getBounds(){
	utility::vector1<core::Real> tmp(6);
	bool first=true;
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			if ( grid_[point.x][point.y][point.z]==TP_POCKET || grid_[point.x][point.y][point.z]==TP_BURIED ||
					grid_[point.x][point.y][point.z]==TP_EDGE ) {
				if ( first ) {
					tmp[1]=xcorn_+stepSize_*point.x;
					tmp[2]=xcorn_+stepSize_*point.x;
					tmp[3]=ycorn_+stepSize_*point.y;
					tmp[4]=ycorn_+stepSize_*point.y;
					tmp[5]=zcorn_+stepSize_*point.z;
					tmp[6]=zcorn_+stepSize_*point.z;
					first=false;
				} else {
					if ( tmp[1]>xcorn_+stepSize_*point.x ) tmp[1]=xcorn_+stepSize_*point.x;
					if ( tmp[2]<xcorn_+stepSize_*point.x ) tmp[2]=xcorn_+stepSize_*point.x;
					if ( tmp[3]>ycorn_+stepSize_*point.y ) tmp[3]=ycorn_+stepSize_*point.y;
					if ( tmp[4]<ycorn_+stepSize_*point.y ) tmp[4]=ycorn_+stepSize_*point.y;
					if ( tmp[5]>zcorn_+stepSize_*point.z ) tmp[5]=zcorn_+stepSize_*point.z;
					if ( tmp[6]<zcorn_+stepSize_*point.z ) tmp[6]=zcorn_+stepSize_*point.z;
				}
			}
		}
	}
	return tmp;
}

void PocketGrid::write_pocketGrid_to_pdb( std::string const & output_filename ) {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_filename, std::ios::out);
	int counter=1;
	int counter2=9;

	outPDB_stream<<"ATOM      1  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<0*stepSize_+xcorn_<<std::setw(8)<<0*stepSize_+ycorn_<<std::setw(8)<<0*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      2  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<0*stepSize_+xcorn_<<std::setw(8)<<(ydim_-1)*stepSize_+ycorn_<<std::setw(8)<<0*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      3  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<0*stepSize_+xcorn_<<std::setw(8)<<(ydim_-1)*stepSize_+ycorn_<<std::setw(8)<<(zdim_-1)*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      4  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<0*stepSize_+xcorn_<<std::setw(8)<<0*stepSize_+ycorn_<<std::setw(8)<<(zdim_-1)*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      5  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<(xdim_-1)*stepSize_+xcorn_<<std::setw(8)<<0*stepSize_+ycorn_<<std::setw(8)<<0*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      6  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<(xdim_-1)*stepSize_+xcorn_<<std::setw(8)<<(ydim_-1)*stepSize_+ycorn_<<std::setw(8)<<0*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      7  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<(xdim_-1)*stepSize_+xcorn_<<std::setw(8)<<0*stepSize_+ycorn_<<std::setw(8)<<(zdim_-1)*stepSize_+zcorn_<<std::endl;
	outPDB_stream<<"ATOM      8  C   C            "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<(xdim_-1)*stepSize_+xcorn_<<std::setw(8)<<(ydim_-1)*stepSize_+ycorn_<<std::setw(8)<<(zdim_-1)*stepSize_+zcorn_<<std::endl;

	core::Size x,y,z;
	for ( x=0; x<(xdim_); x++ ) {
		for ( y=0; y<(ydim_); y++ ) {
			for ( z=0; z<(zdim_); z++ ) {
				std::string concatenated_pdb_info;
				concatenated_pdb_info += "ATOM  ";
				std::stringstream  tmp;
				tmp<<counter2;
				if ( counter2<10 ) concatenated_pdb_info += "    ";
				else if ( counter2<100 ) concatenated_pdb_info += "   ";
				else if ( counter2<1000 ) concatenated_pdb_info += "  ";
				else if ( counter2<10000 ) concatenated_pdb_info += " ";
				concatenated_pdb_info += tmp.str()+"  ";
				if ( grid_[x][y][z]==EMPTY ) {
					concatenated_pdb_info += "EMP EMP";
				}
				if ( grid_[x][y][z]==POCKET )  {
					concatenated_pdb_info += "PKT PKT";
				}
				if ( grid_[x][y][z]==PO_SURF ) {
					concatenated_pdb_info += "POS POS";
				}
				if ( grid_[x][y][z]==PO_BURIED )  {
					concatenated_pdb_info += "BUR BUR";
				}
				if ( grid_[x][y][z]==HSURFACE ) {
					concatenated_pdb_info += "HSU HSU";
				}
				if ( grid_[x][y][z]==PSURFACE ) {
					concatenated_pdb_info += "PSU PSU";
				}
				if ( grid_[x][y][z]==TP_POCKET )  {
					concatenated_pdb_info += "TPP TPP";
				}
				if ( grid_[x][y][z]==TP_SURF )  {
					concatenated_pdb_info += "TPS TPS";
				}
				if ( grid_[x][y][z]==TP_BURIED )  {
					concatenated_pdb_info += "TPB TPB";
				}
				if ( grid_[x][y][z]==PO_EDGE )  {
					concatenated_pdb_info += "POE POE";
				}
				if ( grid_[x][y][z]==TP_EDGE )  {
					concatenated_pdb_info += "TPE TPE";
				}
				if ( grid_[x][y][z]==T_SURFACE )  {
					concatenated_pdb_info += "TS  TS ";
				}
				if ( grid_[x][y][z]==ST_SURFACE )  {
					concatenated_pdb_info += "STS STS";
				}
				if ( grid_[x][y][z]==TARGET )  {
					concatenated_pdb_info += "TGT TGT";
				}
				if ( grid_[x][y][z]==SUBTARGET )  {
					concatenated_pdb_info += "STG STG";
				}
				if ( grid_[x][y][z]==PROTEIN ) {
					concatenated_pdb_info += "PR  PR ";
				}
				if ( grid_[x][y][z]==PSP ) {
					concatenated_pdb_info += "PSP PSP";
				}

				tmp.str(std::string());
				tmp<<"          "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x*stepSize_+xcorn_<<std::setw(8)<<y*stepSize_+ycorn_<<std::setw(8)<<z*stepSize_+zcorn_<<std::endl;
				concatenated_pdb_info += tmp.str();
				counter++;
				outPDB_stream<<concatenated_pdb_info;
			}
		}
	}

	outPDB_stream.close();
	outPDB_stream.clear();

}

void  PocketGrid::alter_espGrid_with_bound_ligand( std::string const & espGrid_filename, core::pose::Pose const & protein_pose ) {
	ElectrostaticpotentialGrid oldG;
	oldG.get_oegrid_dimensions(espGrid_filename);
	oldG.fill_espGrid_values(espGrid_filename);

	debug_assert (stepSize_ == oldG.espGrid_spacing_);
	debug_assert (spacing_  == oldG.espGrid_spacing_);

	ElectrostaticpotentialGrid newG;
	newG.espGrid_.clear();
	newG.espGrid_spacing_ = oldG.espGrid_spacing_;
	newG.espGrid_dim_.x() = xdim_;
	newG.espGrid_dim_.y() = ydim_;
	newG.espGrid_dim_.z() = zdim_;

	newG.espGrid_mid_.x() = xcorn_+ ((xdim_-1)/2*newG.espGrid_spacing_);
	newG.espGrid_mid_.y() = ycorn_+ ((ydim_-1)/2*newG.espGrid_spacing_);
	newG.espGrid_mid_.z() = zcorn_+ ((zdim_-1)/2*newG.espGrid_spacing_);

	core::Real espGrid_xcorn = oldG.espGrid_mid_.x() - ((oldG.espGrid_dim_.x()-1)/2*oldG.espGrid_spacing_);
	core::Real espGrid_ycorn = oldG.espGrid_mid_.y() - ((oldG.espGrid_dim_.y()-1)/2*oldG.espGrid_spacing_);
	core::Real espGrid_zcorn = oldG.espGrid_mid_.z() - ((oldG.espGrid_dim_.z()-1)/2*oldG.espGrid_spacing_);

	core::Real xskip = (xcorn_-espGrid_xcorn)/stepSize_;
	core::Real yskip = (ycorn_-espGrid_ycorn)/stepSize_;
	core::Real zskip = (zcorn_-espGrid_zcorn)/stepSize_;

	std::cout<<"spacing : "<<spacing_<<" "<<stepSize_<<" "<<newG.espGrid_spacing_<<" "<<oldG.espGrid_spacing_<<std::endl;
	std::cout<<"new dimensions : "<<newG.espGrid_dim_.x()<<" "<<newG.espGrid_dim_.y()<<" "<<newG.espGrid_dim_.z()<<std::endl;
	std::cout<<"old dimensions : "<<xdim_<<" "<<ydim_<<" "<<zdim_<<std::endl;
	std::cout<<"old dimensions : "<<oldG.espGrid_dim_.x()<<" "<<oldG.espGrid_dim_.y()<<" "<<oldG.espGrid_dim_.z()<<std::endl;
	std::cout<<"new center : "<<newG.espGrid_mid_.x()<<" "<<newG.espGrid_mid_.y()<<" "<<newG.espGrid_mid_.z()<<std::endl;
	std::cout<<"old center : "<<oldG.espGrid_mid_.x()<<" "<<oldG.espGrid_mid_.y()<<" "<<oldG.espGrid_mid_.z()<<std::endl;
	std::cout<<"corners : "<<espGrid_xcorn<<" "<<espGrid_ycorn<<" "<<espGrid_zcorn<<std::endl;
	std::cout<<"new corners : "<<xcorn_<<" "<<ycorn_<<" "<<zcorn_<<std::endl;
	std::cout<<"SKIP : "<<xskip<<" "<<yskip<<" "<<zskip<<std::endl;

	newG.espGrid_.resize(xdim_);
	for ( core::Size tx=0; tx<xdim_; tx++ ) {
		newG.espGrid_[tx].resize(ydim_);
		for ( core::Size ty=0; ty<ydim_; ty++ ) {
			newG.espGrid_[tx][ty].resize(zdim_, 0);
		}
	}


	//core::Real const min_esp(5);
	//core::Real const max_esp(-5);

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open("highesp.pdb", std::ios::out);
	Size count = 1;
	core::Real prot_grid_dist;
	for ( int x=(int)(0+xskip),newx=0; x<(int)(xskip+xdim_); ++x,++newx ) {
		for ( int y=(int)(0+yskip),newy=0; y<(int)(yskip+ydim_); ++y,++newy ) {
			for ( int z=(int)(0+zskip),newz=0; z<(int)(zskip+zdim_); ++z,++newz ) {

				numeric::xyzVector<core::Real> grid_coord, prot_coord;
				grid_coord.x() = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
				grid_coord.y() = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
				grid_coord.z() = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );

				bool found_protein = false;
				for ( Size j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
					core::conformation::Residue const & curr_rsd = protein_pose.residue(j);
					if ( !curr_rsd.is_protein() ) continue;
					for ( Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i ) {
						if ( curr_rsd.atom_is_hydrogen(i) ) { prot_grid_dist = 0.53; } else { prot_grid_dist = 1.6; }
						prot_coord.x() = curr_rsd.atom(i).xyz()(1);
						prot_coord.y() = curr_rsd.atom(i).xyz()(2);
						prot_coord.z() = curr_rsd.atom(i).xyz()(3);
						// std::cout<<"HETATM"<<std::setw(5)<<count<<"  C   PRT A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< prot_coord.x() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< prot_coord.y() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< prot_coord.z() <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
						if ( grid_coord.distance(prot_coord) < prot_grid_dist ) {
							newG.espGrid_[newx][newy][newz] = 1000;
							found_protein = true;
							//std::cout<<"Found prot" <<std::endl;
							outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   PRO A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< grid_coord.x() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< grid_coord.y() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< grid_coord.z() <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
						}
						if ( found_protein ) break;
					}
					if ( found_protein ) break;
				}

				if ( !found_protein ) {
					newG.espGrid_[newx][newy][newz] = oldG.espGrid_[x][y][z];
					outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   SOL A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< grid_coord.x() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< grid_coord.y() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< grid_coord.z() <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
				}
				/*
				//Assign zero potentials grid points inside protein
				if ( grid_[newx][newy][newz]==PROTEIN || grid_[newx][newy][newz]==TARGET ) { //|| grid_[newx][newy][newz]==EMPTY ) {
				newG.espGrid_[newx][newy][newz] = 1000;
				core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
				core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
				core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
				outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   PRO A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
				}
				//Assign potentials for surface points
				else if( (grid_[newx][newy][newz] == HSURFACE) || (grid_[newx][newy][newz] == PSURFACE) || (grid_[newx][newy][newz] == T_SURFACE) || (grid_[newx][newy][newz] == ST_SURFACE) ){
				newG.espGrid_[newx][newy][newz] = 0.;
				core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
				core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
				core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
				outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   SUR A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
				}
				//Assign potentials for other points
				else {
				//check whether the potentials exceed the max & min bounds
				if ( oldG.espGrid_[x][y][z] > max_esp ) {
				newG.espGrid_[newx][newy][newz] = max_esp;
				core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
				core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
				core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
				outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   MAX A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
				}
				else if ( oldG.espGrid_[x][y][z] < min_esp ) {
				newG.espGrid_[newx][newy][newz] = min_esp;
				core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
				core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
				core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
				outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   MIN A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
				}
				else {
				newG.espGrid_[newx][newy][newz] = oldG.espGrid_[x][y][z];
				core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
				core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
				core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
				outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   NML A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
				}
				}
				*/
			}
		}
	}

	std::list< numeric::xyzVector<core::Real> > surfacePoints_list;
	surfacePoints_list = get_connolly_surfacePoints( protein_pose );

	bool found_surface = false;
	for ( Size iz=0; iz<newG.espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<newG.espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<newG.espGrid_dim_.x(); ++ix ) {
				if ( newG.espGrid_[ix][iy][iz] > 999. ) continue; // skip points that are already assigned as protein
				numeric::xyzVector<core::Real> grid_coord;
				grid_coord.x() = ix * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
				grid_coord.y() = iy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
				grid_coord.z() = iz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
				found_surface = false;
				for ( std::list< numeric::xyzVector<core::Real> >::const_iterator sf = surfacePoints_list.begin(); sf != surfacePoints_list.end(); ++sf ) {
					if ( grid_coord.distance(*sf) < 0.5 ) {
						newG.espGrid_[ix][iy][iz] = 0.;
						found_surface = true;
					}
					if ( found_surface ) break;
				}
			}
		}
	}
	outPDB_stream.close();
	outPDB_stream.clear();

	std::string output_espGrid_filename = "_" + espGrid_filename;
	newG.print_to_oegrid(output_espGrid_filename);
	newG.write_espGrid_to_pdb("espGrid.pdb");

	return;

}

std::list< numeric::xyzVector<core::Real> > PocketGrid::get_connolly_surfacePoints( core::pose::Pose const & protein_pose ) const {

	std::list< numeric::xyzVector<core::Real> > surfacePoints_list;
	surfacePoints_list.clear();
	using namespace core::scoring::sc;
	core::scoring::sc::MolecularSurfaceCalculator calculator;
	calculator.Init();
	calculator.Calc(protein_pose);
	std::vector<DOT> surface_dots = calculator.GetDots(0);
	//std::cout << "Generated surface dots: " << surface_dots.size() << std::endl;
	//core::Size zeronormal_dots = 0;
	//core::Size skipped_dots = 0;
	//core::Size selected_dots = 0;
	core::pack::task::PackerTaskOP task;
	task = core::pack::task::TaskFactory::create_packer_task( protein_pose );
	for ( auto & surface_dot : surface_dots ) {
		if ( task->pack_residue(surface_dot.atom->nresidue) ) {
			numeric::xyzVector<core::Real> surf_coord;
			surf_coord.x() = surface_dot.coor.x();
			surf_coord.y() = surface_dot.coor.y();
			surf_coord.z() = surface_dot.coor.z();
			surfacePoints_list.push_back(surf_coord);
		}
	}

	return surfacePoints_list;
}

std::list< numeric::xyzVector<core::Real> > PocketGrid::get_connolly_surfacePoints_within_grid( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list) {

	std::list< numeric::xyzVector<core::Real> > surfacePoints_within_grid;
	surfacePoints_within_grid.clear();
	numeric::xyzVector<core::Real> grid_coord;
	for ( auto const & sf : surfacePoints_list ) {
		convert_cartesian_to_grid(sf, mid_, dim_, spacing_, grid_coord);
		//ignore if the surface point goes out of the grid
		if ( ((grid_coord.x() >= dim_.x())||(grid_coord.x() < 0)) || ((grid_coord.y() >= dim_.y())||(grid_coord.y() < 0)) || ((grid_coord.z() >= dim_.z())||(grid_coord.z() < 0)) ) {
			continue;
		} else {
			surfacePoints_within_grid.push_back(sf);
		}
	}

	return surfacePoints_within_grid;

}


void PocketGrid::alter_espGrid( std::string const & espGrid_filename ) {

	ElectrostaticpotentialGrid oldG;
	oldG.get_oegrid_dimensions(espGrid_filename);
	oldG.fill_espGrid_values(espGrid_filename);

	debug_assert (stepSize_ == oldG.espGrid_spacing_);
	debug_assert (spacing_  == oldG.espGrid_spacing_);

	ElectrostaticpotentialGrid newG;
	newG.espGrid_.clear();
	newG.espGrid_spacing_ = oldG.espGrid_spacing_;
	newG.espGrid_dim_.x() = xdim_;
	newG.espGrid_dim_.y() = ydim_;
	newG.espGrid_dim_.z() = zdim_;

	newG.espGrid_mid_.x() = xcorn_+ ((xdim_-1)/2*newG.espGrid_spacing_);
	newG.espGrid_mid_.y() = ycorn_+ ((ydim_-1)/2*newG.espGrid_spacing_);
	newG.espGrid_mid_.z() = zcorn_+ ((zdim_-1)/2*newG.espGrid_spacing_);

	core::Real espGrid_xcorn = oldG.espGrid_mid_.x() - ((oldG.espGrid_dim_.x()-1)/2*oldG.espGrid_spacing_);
	core::Real espGrid_ycorn = oldG.espGrid_mid_.y() - ((oldG.espGrid_dim_.y()-1)/2*oldG.espGrid_spacing_);
	core::Real espGrid_zcorn = oldG.espGrid_mid_.z() - ((oldG.espGrid_dim_.z()-1)/2*oldG.espGrid_spacing_);

	core::Real xskip = (xcorn_-espGrid_xcorn)/stepSize_;
	core::Real yskip = (ycorn_-espGrid_ycorn)/stepSize_;
	core::Real zskip = (zcorn_-espGrid_zcorn)/stepSize_;

	std::cout<<"spacing : "<<spacing_<<" "<<stepSize_<<" "<<newG.espGrid_spacing_<<" "<<oldG.espGrid_spacing_<<std::endl;
	std::cout<<"new dimensions : "<<newG.espGrid_dim_.x()<<" "<<newG.espGrid_dim_.y()<<" "<<newG.espGrid_dim_.z()<<std::endl;
	std::cout<<"old dimensions : "<<xdim_<<" "<<ydim_<<" "<<zdim_<<std::endl;
	std::cout<<"old dimensions : "<<oldG.espGrid_dim_.x()<<" "<<oldG.espGrid_dim_.y()<<" "<<oldG.espGrid_dim_.z()<<std::endl;
	std::cout<<"new center : "<<newG.espGrid_mid_.x()<<" "<<newG.espGrid_mid_.y()<<" "<<newG.espGrid_mid_.z()<<std::endl;
	std::cout<<"old center : "<<oldG.espGrid_mid_.x()<<" "<<oldG.espGrid_mid_.y()<<" "<<oldG.espGrid_mid_.z()<<std::endl;
	std::cout<<"corners : "<<espGrid_xcorn<<" "<<espGrid_ycorn<<" "<<espGrid_zcorn<<std::endl;
	std::cout<<"new corners : "<<xcorn_<<" "<<ycorn_<<" "<<zcorn_<<std::endl;
	std::cout<<"SKIP : "<<xskip<<" "<<yskip<<" "<<zskip<<std::endl;

	newG.espGrid_.resize(xdim_);
	for ( core::Size tx=0; tx<xdim_; tx++ ) {
		newG.espGrid_[tx].resize(ydim_);
		for ( core::Size ty=0; ty<ydim_; ty++ ) {
			newG.espGrid_[tx][ty].resize(zdim_, 0);
		}
	}

	using namespace basic::options;
	core::Real const min_esp(50);
	core::Real const max_esp(-50);

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open("highesp.pdb", std::ios::out);
	Size count = 1;
	for ( int x=(int)(0+xskip),newx=0; x<(int)(xskip+xdim_); ++x,++newx ) {
		for ( int y=(int)(0+yskip),newy=0; y<(int)(yskip+ydim_); ++y,++newy ) {
			for ( int z=(int)(0+zskip),newz=0; z<(int)(zskip+zdim_); ++z,++newz ) {
				//Assign potentials only for pocket points, all other points are set to zeroa
				if ( grid_[newx][newy][newz]==POCKET || grid_[newx][newy][newz]==PO_SURF || grid_[newx][newy][newz]==PO_BURIED || grid_[newx][newy][newz]==T_SURFACE || grid_[newx][newy][newz]==ST_SURFACE || grid_[newx][newy][newz]==TP_POCKET || grid_[newx][newy][newz]==TP_SURF || grid_[newx][newy][newz]==TP_BURIED || grid_[newx][newy][newz]==PO_EDGE || grid_[newx][newy][newz]==TP_EDGE ) {

					//check whether the potentials exceed the max & min bounds
					if ( oldG.espGrid_[x][y][z] > max_esp ) {
						newG.espGrid_[newx][newy][newz] = max_esp;
						core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
						core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
						core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
						outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   MAX A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;

					} else if ( oldG.espGrid_[x][y][z] < min_esp ) {
						newG.espGrid_[newx][newy][newz] = min_esp;
						core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
						core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
						core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
						outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  N   MIN A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;

					} else {
						newG.espGrid_[newx][newy][newz] = oldG.espGrid_[x][y][z];
						core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
						core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
						core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
						outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  O   NML A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;

					}

					//Also include surface points that touches the pocket
					core::Size searchxmin, searchxmax, searchymin, searchymax, searchzmin, searchzmax, xx,yy,zz;
					if ( newx != 0 ) {
						searchxmin=newx-1;
					} else {
						searchxmin=newx;
					}
					if ( ((unsigned long)newx) != newG.espGrid_dim_.x()-1 ) {
						searchxmax=newx+1;
					} else {
						searchxmax=newx;
					}
					if ( newy != 0 ) {
						searchymin=newy-1;
					} else {
						searchymin=newy;
					}
					if ( ((unsigned long)newy) != newG.espGrid_dim_.y()-1 ) {
						searchymax=newy+1;
					} else {
						searchymax=newy;
					}
					if ( newz != 0 ) {
						searchzmin=newz-1;
					} else {
						searchzmin=z;
					}
					if ( ((unsigned long)newz) != newG.espGrid_dim_.z()-1 ) {
						searchzmax=newz+1;
					} else {
						searchzmax=newz;
					}
					bool found = false;
					for ( xx = searchxmin; xx<=searchxmax; ++xx ) {
						for ( yy = searchymin; yy<=searchymax; ++yy ) {
							for ( zz = searchzmin; zz<=searchzmax; ++zz ) {
								if ( xx==Size(newx) && yy==Size(newy) && zz==Size(newz) ) continue;
								if ( (grid_[xx][yy][zz] == HSURFACE) || (grid_[xx][yy][zz] == PSURFACE) || (grid_[xx][yy][zz] == T_SURFACE) || (grid_[xx][yy][zz] == ST_SURFACE) ) {
									//check whether the potentials exceed the max & min bounds
									if ( oldG.espGrid_[xx][yy][zz] > max_esp ) {
										newG.espGrid_[xx][yy][zz] = max_esp;
										core::Real xcoord = xx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
										core::Real ycoord = yy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
										core::Real zcoord = zz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
										outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  C   SMX A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[xx][yy][zz]<<std::endl;

									} else if ( oldG.espGrid_[xx][yy][zz] < min_esp ) {
										newG.espGrid_[xx][yy][zz] = min_esp;
										core::Real xcoord = xx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
										core::Real ycoord = yy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
										core::Real zcoord = zz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
										outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  N   SMI A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[xx][yy][zz]<<std::endl;

									} else {
										newG.espGrid_[xx][yy][zz] = oldG.espGrid_[xx][yy][zz];
										core::Real xcoord = xx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
										core::Real ycoord = yy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
										core::Real zcoord = zz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
										outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  O   SML A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[xx][yy][zz]<<std::endl;

									}
									found = true;
								}
							}
							if ( found ) break;
						}
						if ( found ) break;
					}
				} else {
					newG.espGrid_[newx][newy][newz] = 0;
					core::Real xcoord = newx * newG.espGrid_spacing_ + ( newG.espGrid_mid_.x() - ((static_cast<core::Real>(newG.espGrid_dim_.x()-1)/2) * newG.espGrid_spacing_) );
					core::Real ycoord = newy * newG.espGrid_spacing_ + ( newG.espGrid_mid_.y() - ((static_cast<core::Real>(newG.espGrid_dim_.y()-1)/2) * newG.espGrid_spacing_) );
					core::Real zcoord = newz * newG.espGrid_spacing_ + ( newG.espGrid_mid_.z() - ((static_cast<core::Real>(newG.espGrid_dim_.z()-1)/2) * newG.espGrid_spacing_) );
					outPDB_stream<<"HETATM"<<std::setw(5)<<count<<"  S   ZRO A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<newG.espGrid_[newx][newy][newz]<<std::endl;
				}
			}
		}
	}

	outPDB_stream.close();
	outPDB_stream.clear();

	std::string output_espGrid_filename = "new_" + espGrid_filename;
	newG.print_to_oegrid(output_espGrid_filename);

	return;

}


bool PocketGrid::autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::scoring::func::XYZ_Func const & xyz_func, Size const total_residues, bool center_target, core::Real x, core::Real y, core::Real z ) {
	std::vector< core::conformation::ResidueCOP > residues;
	residues.push_back(central_rsd.clone());
	core::pose::Pose tmp_pose;
	int term=1;
	for ( Size j = 1, resnum = total_residues; j <= resnum; ++j ) {
		if ( term ) {
			tmp_pose.append_residue_by_jump(xyz_func.residue(j), 1);
		} else tmp_pose.append_residue_by_bond(xyz_func.residue(j));
		if ( xyz_func.residue(j).is_terminus() ) {
			term=1;
		} else {
			term=0;
		}
	}

	return autoexpanding_pocket_eval(residues, tmp_pose, center_target, x, y, z);
}


bool PocketGrid::autoexpanding_pocket_eval( std::vector< core::conformation::ResidueCOP > const & central_rsds, core::scoring::func::XYZ_Func const & xyz_func, Size const total_residues, bool center_target, core::Real x, core::Real y, core::Real z ) {
	core::pose::Pose tmp_pose;
	int term=1;
	for ( Size j = 1, resnum = total_residues; j <= resnum; ++j ) {
		if ( term ) {
			tmp_pose.append_residue_by_jump(xyz_func.residue(j), 1);
		} else tmp_pose.append_residue_by_bond(xyz_func.residue(j));
		if ( xyz_func.residue(j).is_terminus() ) {
			term=1;
		} else {
			term=0;
		}
	}
	return autoexpanding_pocket_eval(central_rsds, tmp_pose,center_target, x, y, z);
}


bool PocketGrid::autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::pose::Pose const & inPose, bool center_target, core::Real x, core::Real y, core::Real z ) {
	std::vector< core::conformation::ResidueCOP > residues;
	residues.push_back(central_rsd.clone());
	return autoexpanding_pocket_eval(residues, inPose, center_target, x, y, z);
}


bool PocketGrid::autoexpanding_pocket_eval( std::vector< core::conformation::ResidueCOP > const & central_rsds, core::pose::Pose const & inPose, bool center_target, core::Real x, core::Real y, core::Real z  ) {

	using namespace core::chemical;
	bool too_small=true;
	Size total_residues = inPose.size();
	while ( too_small ) {
		if ( center_target ) {
			recenter( central_rsds );
		} else {
			recenter(x,y,z);
		}
		for ( Size j = 1, resnum = total_residues; j <= resnum; ++j ) {
			core::conformation::Residue const & rsd( inPose.conformation().residue(j) );
			//ligands should be ignored and treated like solvent
			if ( !rsd.is_protein() ) continue;
			int target=0;
			int sz = central_rsds.size();
			//this should restrict the recentering to those atoms that define the grid.
			if ( sz > (int)MAX_TARGETS ) sz = (int)MAX_TARGETS;
			for ( int rnum = 0; rnum < sz; ++rnum ) {
				if ( j == central_rsds[rnum]->seqpos() ) target=rnum+1;
			}
			core::Size total_atoms(0);
			using namespace basic::options;
			if ( option[ OptionKeys::fingerprint::include_hydrogens ]() ) {
				total_atoms = rsd.natoms();
			} else {
				total_atoms = rsd.nheavyatoms();
			}
			for ( Size i = 1, i_end = total_atoms; i <= i_end; ++i ) {
				int target_res=target;
				if ( target_res>0 && side_chains_only_ ) {
					if ( (rsd.atom(i).type()>=18)&&(rsd.atom(i).type()<=21) ) {
						target_res=0;
					}
				}

				core::Real vdwOffset=0.;
				//if (rsd.atom_type(i).hybridization() == RING_HYBRID)
				//if (rsd.atom_type(i).is_aromatic())
				// vdwOffset =  stepSize_;
				numeric::xyzVector<core::Real> rpoint = rotatePoint(rsd.atom(i).xyz().x(),rsd.atom(i).xyz().y(),rsd.atom(i).xyz().z());

				mark(rpoint, rsd.atom_type(i).lj_radius()-vdwOffset,probe_rad_, rsd.is_polar(), target_res);
			}
		}

		findPockets(0, maxLen_);
		//core::Size xx,yy,zz;

		if ( markpsp_ ) {
			findPSP(0,maxLen_);
		}

		if ( marksps_ ) {
			findSPS(0,ceil((xdim_+ydim_+zdim_)/stepSize_));
		}
		markPocketDepth(surf_dist_, bur_dist_);

		findClusters();


		if ( static_grid_ == false ) {
			//  Code to check edge of grid goes here
			if ( isTooSmall() ) {
				core::Vector dims=whatIsTooSmall();
				if ( (size_x_ + dims(1) > limit_x_) || (size_y_ + dims(2) > limit_y_) ||(size_z_ + dims(3) > limit_z_) ) {
					if ( dumpExemplars_ || exemplarRestriction_ ) findExemplars(inPose, total_residues);
					if ( exemplarRestriction_ ) findClustersByExemplars();
					std::cerr<<"PocketGrid limit exceded\n";
					return false;
				}
				size_x_ = size_x_ + dims(1);
				size_y_ = size_y_ + dims(2);
				size_z_ = size_z_ + dims(3);
				initialize(central_rsds, size_x_, size_y_, size_z_, spacing_, markpsp_, marksps_);
			} else {
				too_small=false;
			}
		} else if ( static_grid_ == true ) {
			too_small = false;
		}
	}

	if ( dumpExemplars_ || exemplarRestriction_ ) findExemplars(inPose, total_residues);
	if ( exemplarRestriction_ ) findClustersByExemplars();

	/*  core::Size xx,yy,zz;
	for (xx=0;xx<(xdim_); xx++){
	for (yy=0;yy<(ydim_); yy++){
	std::cout << "XYZ: ";
	for (zz=0;zz<(zdim_); zz++){
	std::cout<<grid_[xx][yy][zz]<<" ";
	}
	std::cout<<std::endl;
	}
	std::cout<<"XYZ:"<<std::endl;
	}
	*/

	return true;
}

void PocketGrid::move_pose_to_standard_orie( core::Size const & central_seqpos, core::pose::Pose & pose) {

	// generate the pocket
	autoexpanding_pocket_eval( pose.conformation().residue(central_seqpos), pose );

	// get the pocket CoM and number of pocket points
	core::Size num_pocket_points=0;
	core::Real CoM_x=0.;
	core::Real CoM_y=0.;
	core::Real CoM_z=0.;
	for ( std::list<PCluster>::const_iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit ) {
		if ( cit->points_.size()*pow(stepSize_,3)<minPockSize_ ) continue; // this is a smallpocket
		for ( auto const & point : cit->points_ ) {
			if ( (grid_[point.x][point.y][point.z]==TP_POCKET) || (grid_[point.x][point.y][point.z]==TP_SURF) || (grid_[point.x][point.y][point.z]==TP_EDGE) || (grid_[point.x][point.y][point.z]==TP_BURIED) ) {
				CoM_x += point.x*stepSize_+xcorn_;
				CoM_y += point.y*stepSize_+ycorn_;
				CoM_z += point.z*stepSize_+zcorn_;
				++num_pocket_points;
			}
		}
	}
	CoM_x/=num_pocket_points;
	CoM_y/=num_pocket_points;
	CoM_z/=num_pocket_points;

	Eigen::MatrixXf pocket_coors_matrix(num_pocket_points,3);
	core::Size pocket_inx(0);
	// put pocket points into a matrix for Eigen, subtracting out the CoM so they'll be centered at zero
	for ( std::list<PCluster>::const_iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit ) {
		if ( cit->points_.size()*pow(stepSize_,3)<minPockSize_ ) continue; // this is a smallpocket
		for ( auto const & point : cit->points_ ) {
			if ( (grid_[point.x][point.y][point.z]==TP_POCKET) || (grid_[point.x][point.y][point.z]==TP_SURF) || (grid_[point.x][point.y][point.z]==TP_EDGE) || (grid_[point.x][point.y][point.z]==TP_BURIED) ) {
				// add this pocket point to the Eigen matrix (subtracting out the pocket CoM)
				pocket_coors_matrix(pocket_inx,0) = point.x*stepSize_+xcorn_ - CoM_x;
				pocket_coors_matrix(pocket_inx,1) = point.y*stepSize_+ycorn_ - CoM_y;
				pocket_coors_matrix(pocket_inx,2) = point.z*stepSize_+zcorn_ - CoM_z;
				++pocket_inx;
			}
		}
	}

	// do the SVD calculation
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(pocket_coors_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV );
	svd.computeV();
	Eigen::MatrixXf svd_result(3,3);
	svd_result = svd.matrixV();
	// std::cout << "Matrix V is: " << std::endl << svd.matrixV() << std::endl;
	// MatrixXf transformed = A * svd.matrixV();

	// setup a couple null operations
	numeric::xyzMatrix<core::Real> identity_mat;
	identity_mat.to_identity();
	numeric::xyzVector<core::Real> zero_vec;
	zero_vec.zero();

	// setup the translation vector that will move the pose such that the pocket will be centered at zero
	numeric::xyzVector<core::Real> trans_vec;
	trans_vec.x() = -CoM_x;
	trans_vec.y() = -CoM_y;
	trans_vec.z() = -CoM_z;

	// setup the rotation matrix that will rotate the pose to align the pocket PCA to the axes
	numeric::xyzMatrix<core::Real> rot_mat;
	rot_mat.xx()=svd_result(0,0);
	rot_mat.xy()=svd_result(0,1);
	rot_mat.xz()=svd_result(0,2);
	rot_mat.yx()=svd_result(1,0);
	rot_mat.yy()=svd_result(1,1);
	rot_mat.yz()=svd_result(1,2);
	rot_mat.zx()=svd_result(2,0);
	rot_mat.zy()=svd_result(2,1);
	rot_mat.zz()=svd_result(2,2);

	// first do the translation (with no rotation), then do the rotation (with no translation)
	pose.apply_transform_Rx_plus_v(identity_mat, trans_vec);
	pose.apply_transform_Rx_plus_v(rot_mat, zero_vec);

	// jk DEBUG CODE
	/*
	// transform then print the original pocket using this transform,
	// to ensure PCA works and see how different pocket is when it's recalculated
	utility::io::ozstream outPDB_stream;
	outPDB_stream.open("transformed_starting_grid.pdb", std::ios::out);
	int counter=1;
	int counter2=1;
	for (std::list<PCluster>::const_iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit){
	if (cit->points_.size()*pow(stepSize_,3)<minPockSize_) continue; // this is a smallpocket
	int clustNo=1;
	for (std::list<PCluster::Cxyz>::const_iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
	if ( (grid_[pit->x][pit->y][pit->z]==TP_POCKET) || (grid_[pit->x][pit->y][pit->z]==TP_SURF) || (grid_[pit->x][pit->y][pit->z]==TP_EDGE) || (grid_[pit->x][pit->y][pit->z]==TP_BURIED) ) {
	numeric::xyzVector<core::Real> orig_coors;
	orig_coors.x() = pit->x*stepSize_+xcorn_;
	orig_coors.y() = pit->y*stepSize_+ycorn_;
	orig_coors.z() = pit->z*stepSize_+zcorn_;

	orig_coors += trans_vec;
	numeric::xyzVector<core::Real> new_coors = rot_mat * orig_coors;

	bool smallPocket;
	std::string concatenated_pdb_info;
	concatenated_pdb_info += "ATOM  ";
	std::stringstream  tmp;
	tmp<<counter2;
	if (counter2<10) concatenated_pdb_info += "    ";
	else if (counter2<100) concatenated_pdb_info += "   ";
	else if (counter2<1000) concatenated_pdb_info += "  ";
	else if (counter2<10000) concatenated_pdb_info += " ";
	else concatenated_pdb_info += "";
	concatenated_pdb_info += tmp.str()+"  ";
	if (smallPocket){
	if (grid_[pit->x][pit->y][pit->z]==TP_POCKET) concatenated_pdb_info += "STP STP  ";
	if (grid_[pit->x][pit->y][pit->z]==TP_SURF) concatenated_pdb_info += "STS STS  ";
	if (grid_[pit->x][pit->y][pit->z]==TP_BURIED) concatenated_pdb_info += "STB STB  ";
	if (grid_[pit->x][pit->y][pit->z]==TP_EDGE) concatenated_pdb_info += "STE STE  ";
	}else{
	if (grid_[pit->x][pit->y][pit->z]==TP_POCKET) concatenated_pdb_info += "TP  TP   ";
	if (grid_[pit->x][pit->y][pit->z]==TP_SURF) concatenated_pdb_info += "TPS TPS  ";
	if (grid_[pit->x][pit->y][pit->z]==TP_BURIED) concatenated_pdb_info += "TPB TPB  ";
	if (grid_[pit->x][pit->y][pit->z]==TP_EDGE) concatenated_pdb_info += "TPE TPE  ";
	}
	if (grid_[pit->x][pit->y][pit->z]==POCKET) concatenated_pdb_info += "PC  PC   ";
	if (grid_[pit->x][pit->y][pit->z]==PO_SURF) concatenated_pdb_info += "PCS PCS  ";
	if (grid_[pit->x][pit->y][pit->z]==PO_BURIED) concatenated_pdb_info += "PCB PCB  ";
	if (grid_[pit->x][pit->y][pit->z]==PO_EDGE) concatenated_pdb_info += "PCE PCE  ";
	tmp.str(std::string());
	tmp<<clustNo;
	if (clustNo<10) concatenated_pdb_info += "   ";
	else if (clustNo<100) concatenated_pdb_info += "  ";
	else if (clustNo<1000) concatenated_pdb_info += " ";
	concatenated_pdb_info += tmp.str()+"  ";
	tmp.str(std::string());
	tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coors.x()<<std::setw(8)<<new_coors.y()<<std::setw(8)<<new_coors.z()<<std::endl;
	concatenated_pdb_info += tmp.str();
	counter2++;
	outPDB_stream<<concatenated_pdb_info;
	}
	}
	clustNo++;
	}
	outPDB_stream.close();
	outPDB_stream.clear();
	*/

	// clear the pocket (to be regenerated outside this function)
	init();
	setup_default_options();

	return;
}


core::Real PocketGrid::get_pocket_distance( PocketGrid const & template_pocket, std::string const & comparison_pdbname ) const {


	bool write_comparison_pdb = false;
	// note: string comparison will return non-zero if the strings are NOT equal
	if ( comparison_pdbname.compare("") ) {
		write_comparison_pdb = true;
	}

	utility::io::ozstream comparisonPDB_stream;
	if ( write_comparison_pdb ) comparisonPDB_stream.open(comparison_pdbname, std::ios::out);

	core::Real const full_match_weight = 1.;
	core::Real const partial_match_weight = 0.25;
	core::Real const mismatch_weight = 0.;
	core::Real match_score = 0.;
	core::Size output_res_num = 1;

	// Loop over all points in the template pocket
	core::Size template_num_points = 0;
	core::Size common_points = 0;
	for ( auto cit=template_pocket.clusters_.clusters_.begin(); cit != template_pocket.clusters_.clusters_.end(); ++cit ) {
		if ( cit->points_.size()*pow(template_pocket.stepSize_,3)<template_pocket.minPockSize_ ) continue; // this is a smallpocket
		for ( auto const & point : cit->points_ ) {
			if ( (template_pocket.grid_[point.x][point.y][point.z]==template_pocket.TP_POCKET) || (template_pocket.grid_[point.x][point.y][point.z]==template_pocket.TP_SURF) || (template_pocket.grid_[point.x][point.y][point.z]==template_pocket.TP_EDGE) || (template_pocket.grid_[point.x][point.y][point.z]==template_pocket.TP_BURIED) ) {

				++template_num_points;

				// For this template pocket point, get the cartesian coors and convert to comparison's grid
				core::Real const template_x = point.x*template_pocket.stepSize_+template_pocket.xcorn_;
				core::Real const template_y = point.y*template_pocket.stepSize_+template_pocket.ycorn_;
				core::Real const template_z = point.z*template_pocket.stepSize_+template_pocket.zcorn_;

				core::Size const self_x_index = (core::Size) floor( ( ( template_x - xcorn_ ) / stepSize_ ) + 0.5 );
				core::Size const self_y_index = (core::Size) floor( ( ( template_y - ycorn_ ) / stepSize_ ) + 0.5 );
				core::Size const self_z_index = (core::Size) floor( ( ( template_z - zcorn_ ) / stepSize_ ) + 0.5 );

				//Check to see if template point is within comparison's grid range
				if ( ( self_x_index < xdim_ ) && ( self_y_index < ydim_ ) && ( self_z_index < zdim_ ) ) {
					// note: no need to check that self_index >= 0 since it's of type "Size"
					//    if ( ( self_x_index >= 0 ) && ( self_x_index < xdim_ ) && ( self_y_index >= 0 ) && ( self_y_index < ydim_ ) && ( self_z_index >= 0 ) && ( self_z_index < zdim_ ) ) {

					// jk note: later, we could actually match sub-pocket types - will this give better results??
					//See if template point is pocket for comparison
					if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_POCKET ) || ( grid_[self_x_index][self_y_index][self_z_index] == TP_SURF ) || ( grid_[self_x_index][self_y_index][self_z_index] == TP_EDGE ) || ( grid_[self_x_index][self_y_index][self_z_index] == TP_BURIED ) ) {

						core::Real curr_weight = partial_match_weight;
						++common_points;
						if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_POCKET ) && (template_pocket.grid_[point.x][point.y][point.z]==template_pocket.TP_POCKET) ) {
							curr_weight = full_match_weight;
						} else if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_SURF ) && (template_pocket.grid_[point.x][point.y][point.z]==template_pocket.TP_SURF) ) {
							curr_weight = full_match_weight;
						} else if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_EDGE ) && (template_pocket.grid_[point.x][point.y][point.z]==template_pocket.TP_EDGE) ) {
							curr_weight = full_match_weight;
						} else if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_BURIED ) && (template_pocket.grid_[point.x][point.y][point.z]==template_pocket.TP_BURIED) ) {
							curr_weight = full_match_weight;
						}
						match_score += curr_weight;

						// KK output_res_num needed
						if ( write_comparison_pdb ) comparisonPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   MPO A   "<<std::setw(8)<<output_res_num<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_x<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_y<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_z<<std::endl;
					} else {
						match_score += mismatch_weight;
						// KK output_res_num needed
						if ( write_comparison_pdb ) comparisonPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   NPO A   "<<std::setw(8)<<output_res_num<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_x<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_y<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_z<<std::endl;
					}
				} else {
					//Score as mismatch if it's out of comparison's grid
					match_score += mismatch_weight;
					// KK output_res_num needed
					if ( write_comparison_pdb ) comparisonPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   NPO A   "<<std::setw(8)<<output_res_num<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_x<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_y<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_z<<std::endl;
				}
				output_res_num++;
			}

		}
	}

	if ( write_comparison_pdb ) {
		comparisonPDB_stream.close();
		comparisonPDB_stream.clear();
	}


	core::Size self_num_points = 0;
	for ( auto const & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue; // this is a smallpocket
		for ( const auto & point : cluster.points_ ) {
			if ( (grid_[point.x][point.y][point.z]==TP_POCKET) || (grid_[point.x][point.y][point.z]==TP_SURF) || (grid_[point.x][point.y][point.z]==TP_EDGE) || (grid_[point.x][point.y][point.z]==TP_BURIED) ) {
				++self_num_points;
			}
		}
	}

	// note: non-trivial to average number of points in each pocket, since we don't know how many points in the
	match_score = match_score / ( template_num_points + self_num_points - common_points);
	TR << "PocketGrid match score = " << match_score << std::endl;

	return 1. - sqrt( match_score );
}


bool PocketGrid::operator == ( PocketGrid const & other ) const
{
	if ( ! same_type_as_me( other ) || ! other.same_type_as_me( *this ) ) return false;

	if ( size_x_  != other.size_x_  ) return false;
	if ( size_y_  != other.size_y_  ) return false;
	if ( size_z_  != other.size_z_  ) return false;
	if ( xcorn_   != other.xcorn_   ) return false;
	if ( ycorn_   != other.ycorn_   ) return false;
	if ( zcorn_   != other.zcorn_   ) return false;
	if ( spacing_ != other.spacing_ ) return false;

	return true;
}

bool PocketGrid::same_type_as_me( PocketGrid const & other ) const
{
	PocketGrid const * casted_other = dynamic_cast< PocketGrid const * > (&other);
	return casted_other;
}


TargetPocketGrid::TargetPocketGrid( const PocketGrid& ext_grd): PocketGrid(ext_grd)
{
	init();
	for ( core::Size x=0; x<(ext_grd.xdim_); x++ ) {
		for ( core::Size y=0; y<(ext_grd.ydim_); y++ ) {
			for ( core::Size z=0; z<(ext_grd.zdim_); z++ ) {
				if ( (ext_grd.grid_[x][y][z] == ext_grd.TP_POCKET) || (ext_grd.grid_[x][y][z] == ext_grd.TP_BURIED) || (ext_grd.grid_[x][y][z] == ext_grd.TP_EDGE) ) {
					grid_[x][y][z] = ext_grd.grid_[x][y][z];
				}
			}
		}
	}
	findClusters();
}

TargetPocketGrid::TargetPocketGrid( std::string const & fname ) {
	// NOTE: THIS FUNCTION DOESN'T FULLY SETUP THE GRID
	// HERE WE ONLY PROVIDE MINIMAL FUNCTIONALITY TO DO POCKET VS POCKET COMPARISONS....

	utility::io::izstream instream;
	instream.open(fname, std::ios::in);

	// Read stepSize_ xdim_ ydim_ zdim_ xcorn_ ycorn_ zcorn_
	instream >> stepSize_ >> xdim_ >> ydim_ >> zdim_ >> xcorn_ >> ycorn_ >> zcorn_;

	// Reset the grid for the new dimensions
	init();

	// Read the grid points which should be marked EGGSHELL, use these to push onto eggshell_coord_list_
	// Format is x y z
	core::Size x, y, z;
	std::string type;
	while ( ! instream.eof() ) {
		instream >> x >> y >> z >> type;
		if ( type.compare("TP_POCKET") == 0 ) {
			grid_[x][y][z] = TP_POCKET;
		}
		if ( type.compare("TP_BURIED") == 0 ) {
			grid_[x][y][z] = TP_BURIED;
		}
	}

	instream.close();
	instream.clear();
	findClusters();
}

void TargetPocketGrid::findClusters(){
	//      clusters_.findClusters();
	clusterPockets();
	//mark target clusters
	for ( auto & cluster : clusters_.clusters_ ) {
		cluster.target_=true;
		cluster.subtarget_=true;
		cluster.solventExposed_=true;
	}

}

void TargetPocketGrid::dumpTargetPocketsToPDB( std::string const & output_filename, bool minipock ){

	//this sets up a temporary grid I can tick off pocket points from for minimal pockets
	std::vector < std::vector < std::vector <PtType> > > tmpgrid_;
	if ( minipock ) {
		tmpgrid_.resize(xdim_);
		for ( core::Size tx=0; tx<xdim_; tx++ ) {
			tmpgrid_[tx].resize(ydim_);
			for ( core::Size ty=0; ty<ydim_; ty++ ) {
				tmpgrid_[tx][ty].resize(zdim_, EMPTY);
			}
		}
		for ( auto & cluster : clusters_.clusters_ ) {
			if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
			if ( !cluster.isTarget(numTargets_) ) continue;
			for ( auto & point : cluster.points_ ) {
				tmpgrid_[point.x][point.y][point.z] = grid_[point.x][point.y][point.z];
			}
		}
		/*    for (core::Size tx=1;tx<xdim_-1;tx++)
		for (core::Size ty=1;ty<ydim_-1;ty++)
		for (core::Size tz=1;tz<zdim_-1;tz++){
		if (grid_[tx][ty][tz] == EMPTY) continue;
		for (int dx=-1;dx<2;dx++)
		for (int dy=-1;dy<2;dy++)
		for (int dz=-1;dz<2;dz++){
		if (dz == 0 && dy==0 && dz==0) continue;
		if (grid_[tx+dx][ty+dy][tz+dz] == EMPTY) tmpgrid_[tx+dx][ty+dy][tz+dz] = TP_SURF;
		}
		}
		*/
	}

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_filename, std::ios::out);
	//int counter=1;
	int counter2=1;

	int clustNo=1;
	//bool smallPocket;
	for ( auto & cluster : clusters_.clusters_ ) {
		if ( cluster.points_.size()*pow(stepSize_,3)<minPockSize_ ) continue;
		if ( !cluster.isTarget(numTargets_) ) continue;
		for ( auto & point : cluster.points_ ) {
			if ( minipock ) {
				if ( tmpgrid_[point.x][point.y][point.z] == EMPTY ) continue;
				if ( point.x < 2 || point.y < 2|| point.z < 2 || point.x > xdim_-3 || point.y > ydim_-3|| point.z > zdim_-3 ) continue;
				if ( tmpgrid_[point.x+2][point.y][point.z] != EMPTY && tmpgrid_[point.x-2][point.y][point.z] != EMPTY && tmpgrid_[point.x][point.y+2][point.z] != EMPTY && tmpgrid_[point.x][point.y-2][point.z] != EMPTY &&tmpgrid_[point.x][point.y][point.z+2] != EMPTY && tmpgrid_[point.x][point.y][point.z-2] != EMPTY ) {
					bool pockpt=true;
					for ( int i = -1; i <2; i++ ) {
						for ( int j = -1; j <2; j++ ) {
							for ( int k = -1; k <2; k++ ) {
								if ( i==0 && j==0 && k==0 ) continue;
								if ( tmpgrid_[point.x+i][point.y+j][point.z+k] == EMPTY ) {
									pockpt=false;
								}
							}
						}
					}
					if ( pockpt ) {
						for ( int i = -1; i <2; i++ ) {
							for ( int j = -1; j <2; j++ ) {
								for ( int k = -1; k <2; k++ ) {
									tmpgrid_[point.x+i][point.y+j][point.z+k] = EMPTY;
								}
							}
						}
						tmpgrid_[point.x-2][point.y][point.z] = EMPTY;
						tmpgrid_[point.x][point.y-2][point.z] = EMPTY;
						tmpgrid_[point.x][point.y][point.z-2] = EMPTY;
					} else continue;
				} else continue;
			}
			std::string concatenated_pdb_info;
			concatenated_pdb_info += "HETATM";
			std::stringstream  tmp;
			tmp<<counter2;
			if ( counter2<10 ) concatenated_pdb_info += "    ";
			else if ( counter2<100 ) concatenated_pdb_info += "   ";
			else if ( counter2<1000 ) concatenated_pdb_info += "  ";
			else if ( counter2<10000 ) concatenated_pdb_info += " ";
			else concatenated_pdb_info += "";
			concatenated_pdb_info += tmp.str()+"  ";
			if ( minipock ) {
				concatenated_pdb_info += " C  TMP A";
			} else {
				if ( grid_[point.x][point.y][point.z]==TP_POCKET ) concatenated_pdb_info += "TP  TP   ";
				if ( grid_[point.x][point.y][point.z]==TP_SURF ) concatenated_pdb_info += "TPS TPS  ";
				if ( grid_[point.x][point.y][point.z]==TP_BURIED ) concatenated_pdb_info += "TPB TPB  ";
				if ( grid_[point.x][point.y][point.z]==TP_EDGE ) concatenated_pdb_info += "TPE TPE  ";
			}
			tmp.str(std::string());
			tmp<<clustNo;
			if ( clustNo<10 ) concatenated_pdb_info += "   ";
			else if ( clustNo<100 ) concatenated_pdb_info += "  ";
			else if ( clustNo<1000 ) concatenated_pdb_info += " ";
			concatenated_pdb_info += tmp.str()+"  ";
			tmp.str(std::string());
			tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<point.x*stepSize_+xcorn_<<std::setw(8)<<point.y*stepSize_+ycorn_<<std::setw(8)<<point.z*stepSize_+zcorn_;
			if ( minipock ) tmp << "  1.00  2.03           C";
			tmp << std::endl;
			concatenated_pdb_info += tmp.str();
			counter2++;
			outPDB_stream<<concatenated_pdb_info;

		}
		clustNo++;
	}

	outPDB_stream.close();
	outPDB_stream.clear();

}

void ElectrostaticpotentialGrid::initialize_espGrid(){
	espGrid_.clear();
	espGrid_.resize(espGrid_dim_.x());
	for ( core::Size tx=0; tx<espGrid_dim_.x(); tx++ ) {
		espGrid_[tx].resize(espGrid_dim_.y());
		for ( core::Size ty=0; ty<espGrid_dim_.y(); ty++ ) {
			espGrid_[tx][ty].resize(espGrid_dim_.z(), 0);
		}
	}
}

void ElectrostaticpotentialGrid::initialize_typGrid(){
	typGrid_.clear();
	typGrid_.resize(espGrid_dim_.x());
	for ( core::Size tx=0; tx<espGrid_dim_.x(); tx++ ) {
		typGrid_[tx].resize(espGrid_dim_.y());
		for ( core::Size ty=0; ty<espGrid_dim_.y(); ty++ ) {
			typGrid_[tx][ty].resize(espGrid_dim_.z(), SOLVENT);
		}
	}
}

bool TargetPocketGrid::same_type_as_me( PocketGrid const & other ) const
{
	TargetPocketGrid const * casted_other = dynamic_cast< TargetPocketGrid const * > (&other);
	return casted_other;
}


void ElectrostaticpotentialGrid::set_surface_esp_to_zero() {

	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				if ( typGrid_[ix][iy][iz] == SURFACE ) {
					espGrid_[ix][iy][iz] = 0.;
				}
			}
		}
	}

	return;

}

void ElectrostaticpotentialGrid::set_protein_esp_to_zero() {

	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				if ( typGrid_[ix][iy][iz] == PROTEIN ) {
					espGrid_[ix][iy][iz] = 0.;
				}
			}
		}
	}

	return;

}

void ElectrostaticpotentialGrid::change_espGrid_for_darc( std::string const & espGrid_filename, core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid ) {

	std::list< numeric::xyzVector<core::Real> > surfacePoints_list = pocket_grid.get_connolly_surfacePoints( protein_pose );
	resize_espGrid_to_match_pocketGrid(espGrid_filename, pocket_grid );
	initialize_typGrid();
	assign_esp_for_surface_grid_points( surfacePoints_list );
	//assign_esp_for_surface_grid_points_by_nn( surfacePoints_list );
	assign_esp_for_protein_grid_points( protein_pose, surfacePoints_list );
	mark_buried_solvent_points( protein_pose, surfacePoints_list );
	cap_espGrid();
	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::set_surface_esp_to_zero ]() ) {
		set_surface_esp_to_zero();
	}

	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::set_protein_esp_to_zero ]() ) {
		set_protein_esp_to_zero();
	}

	int pfounddir = espGrid_filename.find_last_of("/\\");
	int pfounddot = espGrid_filename.find_last_of(".");
	std::string tag = espGrid_filename.substr((pfounddir+1),(pfounddot-(pfounddir+1)));
	std::string output_espGrid_filename = "DARC_" + tag + ".agd";
	print_to_oegrid(output_espGrid_filename);
	std::cout<< "Written espGrid file : "<< output_espGrid_filename << std::endl;
	//Print grid points into a PDB file
	using namespace basic::options;
	if ( option[ OptionKeys::pocket_grid::dump_espGrid ]() ) {
		std::string output_espGrid_pdbname = tag + "_espGrid.pdb";
		write_espGrid_to_pdb(output_espGrid_pdbname);
	}
	//Print connolly surface points into a PDB file
	using namespace basic::options;
	if ( option[ OptionKeys::pocket_grid::dump_connollySurface ]() ) {
		std::string output_conSurf_pdbname = tag + "_conSurf.pdb";
		write_connollySurface_to_pdb(surfacePoints_list, output_conSurf_pdbname);
	}

	return;

}

void ElectrostaticpotentialGrid::get_ZAP_espGrid_values_with_type( std::string const & input_espGrid_filename ) {

	std::string input_OEgrid_filename = input_espGrid_filename;
	//get_oegrid_dimensions_from_input_file
	//initialize grid and fill values
	get_oegrid_dimensions(input_OEgrid_filename);
	fill_espGrid_values_with_type(input_OEgrid_filename);

	return;
}

void ElectrostaticpotentialGrid::get_ZAP_espGrid_values( std::string const & input_espGrid_filename ) {

	std::string input_OEgrid_filename = input_espGrid_filename;
	//get_oegrid_dimensions_from_input_file
	//initialize grid and fill values
	get_oegrid_dimensions(input_OEgrid_filename);
	fill_espGrid_values(input_OEgrid_filename);

	return;
}

void ElectrostaticpotentialGrid::cap_espGrid() {

	using namespace basic::options;
	core::Real const cap_esp(option[ OptionKeys::fingerprint::cap_e_potential ]);

	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				if ( espGrid_[ix][iy][iz] > cap_esp ) {
					espGrid_[ix][iy][iz] = cap_esp;
				} else if ( espGrid_[ix][iy][iz] < -(cap_esp) ) {
					espGrid_[ix][iy][iz] = -(cap_esp);
				}
			}
		}
	}

	return;

}

void ElectrostaticpotentialGrid::get_DELPHI_espGrid_values( std::string const & input_espGrid_filename ) {

	std::string input_DELPHI_phimap = input_espGrid_filename;
	std::ifstream inFile(input_DELPHI_phimap.c_str());
	if ( !inFile ) {
		std::cout<< "Can't open input DELPHI phimap file " << input_DELPHI_phimap << std::endl;
		exit(1);
	} else  if ( inFile ) {
		std::cout<< "Input DELPHI phimap file : " << input_DELPHI_phimap << std::endl;
	}
	int linenum = 1;
	std::string line;
	while ( std::getline(inFile, line) ) {
		if ( linenum == 3 ) {
			std::string midX = line.substr(7,14);
			espGrid_mid_.x() = atof(midX.c_str());
			std::string midY = line.substr(21,14);
			espGrid_mid_.y() = atof(midY.c_str());
			std::string midZ = line.substr(35,14);
			espGrid_mid_.z() = atof(midZ.c_str());
		}
		if ( linenum == 4 ) {
			std::string dimX = line.substr(2,3);
			espGrid_dim_.x() = (core::Size)(atof(dimX.c_str()));
		}
		if ( linenum == 5 ) {
			std::string dimY = line.substr(2,3);
			espGrid_dim_.y() = (core::Size)(atof(dimY.c_str()));
		}
		if ( linenum == 6 ) {
			std::string dimZ = line.substr(2,3);
			espGrid_dim_.z() = (core::Size)(atof(dimZ.c_str()));
		}
		if ( linenum == 4 ) {
			std::string spacing = line.substr(7,12);
			espGrid_spacing_ = atof(spacing.c_str());
		}
		linenum++;
	}
	inFile.close();
	//Refer DELPHI manual for these conversions
	core::Real coeff = 0.5291772108;
	core::Real stepsize = espGrid_spacing_ * coeff;
	//core::Real tmp_pt = stepsize * (espGrid_dim_.x() - 1)/2/coeff;
	espGrid_mid_.x() = espGrid_mid_.x() + stepsize * ((espGrid_dim_.x() - 1)/2) * coeff;
	espGrid_mid_.y() = espGrid_mid_.y() + stepsize * ((espGrid_dim_.y() - 1)/2) * coeff;
	espGrid_mid_.z() = espGrid_mid_.z() + stepsize * ((espGrid_dim_.z() - 1)/2) * coeff;
	espGrid_spacing_ = espGrid_spacing_ * coeff;
	//std::cout<<"Midpoint X Y Z : "<<esp_midpoint_x_<<" "<<espGrid_mid_.y()<<" "<<espGrid_mid_.z()<<"\n"<<"Dimension X Y Z : "<<espGrid_dim_.x()<<" "<<espGrid_dim_.y()<<" "<<espGrid_dim_.z()<<"\n"<<"spacing : "<<espGrid_spacing_<<std::endl;
	init();//intialize electrostatic potential grid
	std::ifstream inputFile(input_DELPHI_phimap.c_str());
	std::string newline;
	//skip first 7 lines
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				inputFile >> espGrid_[ix][iy][iz];
			}
		}
	}
	inputFile.close();

	return;
}

ElectrostaticpotentialGrid::ElectrostaticpotentialGrid( std::string const & input_espGrid_filename, core::pose::Pose const &, PocketGrid const &, bool delphi ) {

	if ( delphi==false ) { //Read as OpenEye ZAP grid
		get_ZAP_espGrid_values( input_espGrid_filename );
	} else if ( delphi==true ) { //Read as DELPHI grid
		get_DELPHI_espGrid_values( input_espGrid_filename );
	}

	return;

}

void ElectrostaticpotentialGrid::get_oegrid_dimensions( std::string const & input_OEgrid_filename ) {

	std::ifstream inFile(input_OEgrid_filename.c_str());
	if ( !inFile ) {
		std::cout<< "Can't open input file " << input_OEgrid_filename << std::endl;
		exit(1);
	}

	int linenum = 1;
	std::string line;
	while ( std::getline(inFile, line) ) {
		if ( linenum == 2 ) {
			std::string midX = line.substr(6,12);
			espGrid_mid_.x() = atof(midX.c_str());
			std::string midY = line.substr(19,12);
			espGrid_mid_.y() = atof(midY.c_str());
			std::string midZ = line.substr(32,12);
			espGrid_mid_.z() = atof(midZ.c_str());
		}
		if ( linenum == 3 ) {
			std::string dimX = line.substr(6,6);
			espGrid_dim_.x() = (core::Size)(atof(dimX.c_str()));
			std::string dimY = line.substr(13,6);
			espGrid_dim_.y() = (core::Size)(atof(dimY.c_str()));
			std::string dimZ = line.substr(20,6);
			espGrid_dim_.z() = (core::Size)(atof(dimZ.c_str()));
		}
		if ( linenum == 4 ) {
			std::string spacing = line.substr(10,12);
			espGrid_spacing_ = atof(spacing.c_str());
		}

		linenum++;
	}
	inFile.close();
	//std::cout<<"Midpoint X Y Z : "<<espGrid_mid_.x()<<" "<<espGrid_mid_.y()<<" "<<espGrid_mid_.z()<<"\n"<<"Dimension X Y Z : "<<espGrid_dim_.x()<<" "<<espGrid_dim_.y()<<" "<<espGrid_dim_.z()<<"\n"<<"spacing : "<<espGrid_spacing_<<std::endl;

	return;
}

void ElectrostaticpotentialGrid::print_to_oegrid( std::string const & output_OEgrid_filename ) {

	utility::io::ozstream outFile;
	outFile.open(output_OEgrid_filename, std::ios::out);

	if ( !outFile ) {
		std::cout<< "Can't open output file " << output_OEgrid_filename << std::endl;
		exit(1);
	}

	outFile << "Title: " <<std::endl;
	char buffer[80];
	sprintf(buffer, "Mid: %12.6f %12.6f %12.6f", espGrid_mid_.x(), espGrid_mid_.y(), espGrid_mid_.z());
	outFile << buffer << std::endl;
	sprintf(buffer, "Dim: %6lu %6lu %6lu", (unsigned long) espGrid_dim_.x(), (unsigned long) espGrid_dim_.y(), (unsigned long) espGrid_dim_.z());
	outFile << buffer << std::endl;
	sprintf(buffer, "Spacing: %12.6f", espGrid_spacing_);
	outFile << buffer << std::endl;
	outFile << "Values:" << std::endl;

	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				outFile<<espGrid_[ix][iy][iz]<<'\t'<<typGrid_[ix][iy][iz]<<std::endl;
			}
		}
	}

	outFile.close();

	return;

}

void  ElectrostaticpotentialGrid::resize_espGrid_to_match_pocketGrid( std::string const & espGrid_filename, PocketGrid const & pg ) {

	get_oegrid_dimensions(espGrid_filename);
	fill_espGrid_values(espGrid_filename);
	numeric::xyzVector<core::Size> trim_espGrid_dim = pg.dim_;
	core::Real trim_espGrid_spacing = espGrid_spacing_;
	std::vector < std::vector < std::vector <core::Real> > > trim_espGrid;
	trim_espGrid.clear();
	//assgn mid points for resized esp grid
	numeric::xyzVector<core::Real> trim_espGrid_mid;
	//find the center coordinates of the pocketGrid
	core::Real cen_x = pg.xcorn_+ ((pg.xdim_-1)/2*pg.stepSize_);
	core::Real cen_y = pg.ycorn_+ ((pg.ydim_-1)/2*pg.stepSize_);
	core::Real cen_z = pg.zcorn_+ ((pg.zdim_-1)/2*pg.stepSize_);
	//convert center coordinates of the pocketGrid to grid points(in terms of espGrid)
	core::Real grdX = (cen_x - ( espGrid_mid_.x() - ((static_cast<core::Real>(espGrid_dim_.x()+1)/2) * espGrid_spacing_) ))/espGrid_spacing_;
	core::Real grdY = (cen_y - ( espGrid_mid_.y() - ((static_cast<core::Real>(espGrid_dim_.y()+1)/2) * espGrid_spacing_) ))/espGrid_spacing_;
	core::Real grdZ = (cen_z - ( espGrid_mid_.z() - ((static_cast<core::Real>(espGrid_dim_.z()+1)/2) * espGrid_spacing_) ))/espGrid_spacing_;
	//rounded grid points are the nearest neighbour points
	core::Size nn_grdX = (core::Size) std::floor(grdX+.5);
	core::Size nn_grdY = (core::Size) std::floor(grdY+.5);
	core::Size nn_grdZ = (core::Size) std::floor(grdZ+.5);
	//convert bact to carteesian coordinates and assign as mid points for the new trimmed espGrid
	trim_espGrid_mid.x() = nn_grdX * espGrid_spacing_ + ( espGrid_mid_.x() - ((static_cast<core::Real>(espGrid_dim_.x()+1)/2) * espGrid_spacing_) );
	trim_espGrid_mid.y() = nn_grdY * espGrid_spacing_ + ( espGrid_mid_.y() - ((static_cast<core::Real>(espGrid_dim_.y()+1)/2) * espGrid_spacing_) );
	trim_espGrid_mid.z() = nn_grdZ * espGrid_spacing_ + ( espGrid_mid_.z() - ((static_cast<core::Real>(espGrid_dim_.z()+1)/2) * espGrid_spacing_) );

	//assign midpoint of pocket grid as mid point of esp grid
	//trim_espGrid_mid.x() = pg.xcorn_+ ((pg.xdim_-1)/2*trim_espGrid_spacing);
	//trim_espGrid_mid.y() = pg.ycorn_+ ((pg.ydim_-1)/2*trim_espGrid_spacing);
	//trim_espGrid_mid.z() = pg.zcorn_+ ((pg.zdim_-1)/2*trim_espGrid_spacing);

	core::Real espGrid_xcorn = espGrid_mid_.x() - ((espGrid_dim_.x()-1)/2*espGrid_spacing_);
	core::Real espGrid_ycorn = espGrid_mid_.y() - ((espGrid_dim_.y()-1)/2*espGrid_spacing_);
	core::Real espGrid_zcorn = espGrid_mid_.z() - ((espGrid_dim_.z()-1)/2*espGrid_spacing_);
	core::Real xskip = (pg.xcorn_-espGrid_xcorn)/espGrid_spacing_;
	core::Real yskip = (pg.ycorn_-espGrid_ycorn)/espGrid_spacing_;
	core::Real zskip = (pg.zcorn_-espGrid_zcorn)/espGrid_spacing_;

	trim_espGrid.resize(trim_espGrid_dim.x());
	for ( core::Size tx=0; tx<trim_espGrid_dim.x(); tx++ ) {
		trim_espGrid[tx].resize(trim_espGrid_dim.y());
		for ( core::Size ty=0; ty<trim_espGrid_dim.y(); ty++ ) {
			trim_espGrid[tx][ty].resize(trim_espGrid_dim.z(), 0);
		}
	}
	//std::cout<<" skip :"<<xskip<<" "<<yskip<<" "<<zskip<<std::endl;
	//std::cout<<" dim :"<<trim_espGrid_dim.x()<<" "<<trim_espGrid_dim.y()<<" "<<trim_espGrid_dim.z()<<std::endl;
	//std::cout<<" mid :"<<trim_espGrid_mid.x()<<" "<<trim_espGrid_mid.y()<<" "<<trim_espGrid_mid.z()<<std::endl;
	for ( int x=(int)(0+xskip),newx=0; x<(int)(xskip+trim_espGrid_dim.x()); ++x,++newx ) {
		for ( int y=(int)(0+yskip),newy=0; y<(int)(yskip+trim_espGrid_dim.y()); ++y,++newy ) {
			for ( int z=(int)(0+zskip),newz=0; z<(int)(zskip+trim_espGrid_dim.z()); ++z,++newz ) {
				trim_espGrid[newx][newy][newz] = espGrid_[x][y][z];
			}
		}
	}

	espGrid_.clear();
	espGrid_dim_ = trim_espGrid_dim;
	espGrid_mid_ = trim_espGrid_mid;
	espGrid_spacing_ = trim_espGrid_spacing;
	espGrid_ = trim_espGrid;

	return;

}

void ElectrostaticpotentialGrid::fill_espGrid_values( std::string const & input_OEgrid_filename ) {

	//intialize electrostatic potential grid
	initialize_espGrid();

	std::ifstream inputFile(input_OEgrid_filename.c_str());
	if ( !inputFile ) {
		std::cout<< "Can't open input file " << input_OEgrid_filename << std::endl;
		exit(1);
	}

	std::string newline;
	//skip first 5 lines
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);

	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				inputFile >> espGrid_[ix][iy][iz];
			}
		}
	}
	inputFile.close();

	return;

}

void ElectrostaticpotentialGrid::fill_espGrid_values_with_type( std::string const & input_OEgrid_filename ) {

	//intialize electrostatic potential grid
	initialize_espGrid();
	initialize_typGrid();

	std::ifstream inputFile(input_OEgrid_filename.c_str());
	if ( !inputFile ) {
		std::cout<< "Can't open input file " << input_OEgrid_filename << std::endl;
		exit(1);
	}

	std::string newline;
	//skip first 5 lines
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);
	std::getline(inputFile, newline);

	core::Real a;
	core::Size b;
	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				std::getline(inputFile, newline);
				std::istringstream lin;
				lin.str(newline);
				if ( lin >> a >> b ) {
					espGrid_[ix][iy][iz] = a;
					typGrid_[ix][iy][iz] = PtType(b);
				}
			}
		}
	}

	inputFile.close();

	return;

}

void ElectrostaticpotentialGrid::mark_protein_espGrid_points( core::pose::Pose const & protein_pose ) {

	//Assign ESP=0 for grid points inside protein
	for ( int j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
		core::conformation::Residue const & rsd(protein_pose.conformation().residue(j) );
		core::Size total_atoms(0);
		using namespace basic::options;
		if ( option[ OptionKeys::fingerprint::include_hydrogens ]() ) {
			total_atoms = rsd.natoms();
		} else {
			total_atoms = rsd.nheavyatoms();
		}
		for ( Size i = 1, i_end = total_atoms; i <= i_end; ++i ) {
			mark_atom_espGrid_points(rsd.atom(i).xyz(), rsd.atom_type(i).lj_radius());
		}
	}

	return;

}

void ElectrostaticpotentialGrid::mark_atom_espGrid_points( numeric::xyzVector<core::Real> protein_atom, core::Real const & protein_atom_vdw ) {
	// core::Size count = 1;

	core::Real xcen = (protein_atom.x() - ( espGrid_mid_.x() - ((static_cast<core::Real>(espGrid_dim_.x()+1)/2) * espGrid_spacing_) ))/espGrid_spacing_;
	core::Real ycen = (protein_atom.y() - ( espGrid_mid_.y() - ((static_cast<core::Real>(espGrid_dim_.y()+1)/2) * espGrid_spacing_) ))/espGrid_spacing_;
	core::Real zcen = (protein_atom.z() - ( espGrid_mid_.z() - ((static_cast<core::Real>(espGrid_dim_.z()+1)/2) * espGrid_spacing_) ))/espGrid_spacing_;

	// core::Real xcoord = xcen * espGrid_spacing_ + ( espGrid_mid_.x() - ((static_cast<core::Real>(espGrid_dim_.x()+1)/2) * espGrid_spacing_) );
	// core::Real ycoord = ycen * espGrid_spacing_ + ( espGrid_mid_.y() - ((static_cast<core::Real>(espGrid_dim_.y()+1)/2) * espGrid_spacing_) );
	// core::Real zcoord = zcen * espGrid_spacing_ + ( espGrid_mid_.z() - ((static_cast<core::Real>(espGrid_dim_.z()+1)/2) * espGrid_spacing_) );

	core::Real minX,maxX;
	core::Real minY,maxY;
	core::Real minZ,maxZ;

	using namespace basic::options;
	core::Real const esp_buffer = option[ OptionKeys::pocket_grid::esp_buffer_dist ];//Angstroms

	core::Real radius = (protein_atom_vdw + esp_buffer)/espGrid_spacing_;

	// std::cout<<"RADIUS : "<<protein_atom_vdw<<" "<<esp_buffer<<" "<<radius<<" "<<espGrid_spacing_<<std::endl;

	minX = std::max( ceil( (xcen - radius) - 0.5), 0.);
	minY = std::max( ceil( (ycen - radius) - 0.5), 0.);
	minZ = std::max( ceil( (zcen - radius) - 0.5), 0.);
	maxX = std::min( floor( (xcen + radius) + 0.5 ), espGrid_dim_.x() -1.);
	maxY = std::min( floor( (ycen + radius) + 0.5 ), espGrid_dim_.y() -1.);
	maxZ = std::min( floor( (zcen + radius) + 0.5 ), espGrid_dim_.z() -1.);

	for ( int xIter=(int)minX; xIter<=(int)maxX; ++xIter ) {
		for ( int yIter=(int)minY; yIter<=(int)maxY; ++yIter ) {
			for ( int zIter=(int)minZ; zIter<=(int)maxZ; ++zIter ) {
				//core::Real xcoord = xIter * espGrid_spacing_ + ( espGrid_mid_.x() - ((static_cast<core::Real>(espGrid_dim_.x()+1)/2) * espGrid_spacing_) );
				//core::Real ycoord = yIter * espGrid_spacing_ + ( espGrid_mid_.y() - ((static_cast<core::Real>(espGrid_dim_.y()+1)/2) * espGrid_spacing_) );
				//core::Real zcoord = zIter * espGrid_spacing_ + ( espGrid_mid_.z() - ((static_cast<core::Real>(espGrid_dim_.z()+1)/2) * espGrid_spacing_) );
				if ( pow((xIter-xcen),2)+pow((yIter-ycen),2)+pow((zIter-zcen),2) <= pow(protein_atom_vdw + esp_buffer, 2) ) {
					espGrid_[xIter][yIter][zIter] = 0;
				} else { }
			}
		}
	}

	return;

}

void  ElectrostaticpotentialGrid::assign_esp_for_surface_grid_points( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list ) {

	bool found_surface = false;
	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				numeric::xyzVector<core::Real> grid_cart;
				grid_cart.x() = ix * espGrid_spacing_ + ( espGrid_mid_.x() - ((static_cast<core::Real>(espGrid_dim_.x()-1)/2) * espGrid_spacing_) );
				grid_cart.y() = iy * espGrid_spacing_ + ( espGrid_mid_.y() - ((static_cast<core::Real>(espGrid_dim_.y()-1)/2) * espGrid_spacing_) );
				grid_cart.z() = iz * espGrid_spacing_ + ( espGrid_mid_.z() - ((static_cast<core::Real>(espGrid_dim_.z()-1)/2) * espGrid_spacing_) );
				found_surface = false;
				for ( auto const & sf : surfacePoints_list ) {
					if ( grid_cart.distance(sf) < 0.5 ) {
						typGrid_[ix][iy][iz] = SURFACE;
						found_surface = true;
					}
					if ( found_surface ) break;
				}
			}
		}
	}


	return;
}

void  ElectrostaticpotentialGrid::assign_esp_for_surface_grid_points_by_nn( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list ) {
	//nn - nearest neighbour method
	for ( auto const & sf : surfacePoints_list ) {
		//convert to grid points
		numeric::xyzVector<core::Real> grid_coord;
		convert_cartesian_to_grid(sf, espGrid_mid_, espGrid_dim_, espGrid_spacing_, grid_coord);
		//rounded grid points are the nearest neighbour points
		core::Size nnX = (core::Size) std::floor(grid_coord.x()+.5);
		core::Size nnY = (core::Size) std::floor(grid_coord.y()+.5);
		core::Size nnZ = (core::Size) std::floor(grid_coord.z()+.5);
		// ignore if the surface point goes out of the grid
		if ( ((nnX >= espGrid_dim_.x())||(nnX <= 0)) || ((nnY >= espGrid_dim_.y())||(nnY <= 0)) || ((nnZ >= espGrid_dim_.z())||(nnZ <= 0)) ) {
			continue;
		} else {
			typGrid_[nnX][nnY][nnZ] = SURFACE;
		}
	}

	return;
}

void  ElectrostaticpotentialGrid::assign_esp_for_protein_grid_points( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list ) {

	//Assign ESP for grid points inside protein
	for ( int j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
		core::conformation::Residue const & rsd(protein_pose.conformation().residue(j) );
		core::Size total_atoms(0);
		using namespace basic::options;
		if ( option[ OptionKeys::fingerprint::include_hydrogens ]() ) {
			total_atoms = rsd.natoms();
		} else {
			total_atoms = rsd.nheavyatoms();
		}
		for ( Size i = 1, i_end = total_atoms; i <= i_end; ++i ) {
			numeric::xyzVector<core::Real> grid_coord;
			convert_cartesian_to_grid(rsd.atom(i).xyz(),  espGrid_mid_, espGrid_dim_, espGrid_spacing_, grid_coord);
			if ( ((grid_coord.x() >= espGrid_dim_.x())||(grid_coord.x() < 0)) || ((grid_coord.y() >= espGrid_dim_.y())||(grid_coord.y() < 0)) || ((grid_coord.z() >= espGrid_dim_.z())||(grid_coord.z() < 0)) ) {
				continue;
			} else {
				core::Real inflated_radius = get_atom_dist_from_connolly_surface(rsd.atom(i).xyz(), surfacePoints_list);
				assign_esp_for_atom_grid_points(rsd.atom(i).xyz(), inflated_radius);
			}
		}
	}

	return;

}

void  ElectrostaticpotentialGrid::assign_esp_for_atom_grid_points( numeric::xyzVector<core::Real> atom_coord, core::Real const & inflated_radius ) {

	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				if ( (typGrid_[ix][iy][iz] == PROTEIN) || (typGrid_[ix][iy][iz] == SURFACE) ) continue; // skip points that are already assigned as protein or surface
				numeric::xyzVector<core::Real> grid_cart;
				grid_cart.x() = ix * espGrid_spacing_ + ( espGrid_mid_.x() - ((static_cast<core::Real>(espGrid_dim_.x()-1)/2) * espGrid_spacing_) );
				grid_cart.y() = iy * espGrid_spacing_ + ( espGrid_mid_.y() - ((static_cast<core::Real>(espGrid_dim_.y()-1)/2) * espGrid_spacing_) );
				grid_cart.z() = iz * espGrid_spacing_ + ( espGrid_mid_.z() - ((static_cast<core::Real>(espGrid_dim_.z()-1)/2) * espGrid_spacing_) );
				if ( atom_coord.distance(grid_cart) < inflated_radius ) {
					typGrid_[ix][iy][iz] = PROTEIN;
				}
			}
		}
	}

	return;

}

core::Real  ElectrostaticpotentialGrid::get_atom_dist_from_connolly_surface( numeric::xyzVector<core::Real> atom_coord, std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list ) {

	core::Real min_dist = 999;
	for ( auto const & sf : surfacePoints_list ) {
		core::Real curr_dist = atom_coord.distance(sf);
		if ( curr_dist < min_dist ) {
			min_dist = curr_dist;
		}
	}
	return min_dist;
}

void ElectrostaticpotentialGrid::write_connollySurface_to_pdb( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list, std::string const & output_pdbname) const {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_pdbname, std::ios::out);
	for ( auto const & sf : surfacePoints_list ) {
		outPDB_stream<<"HETATM   "<<std::setw(2)<<3<<"  N   CNY A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf.x() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf.y() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf.z() <<std::endl;
	}
	outPDB_stream.close();
	outPDB_stream.clear();


	return;

}

void  ElectrostaticpotentialGrid::mark_buried_solvent_points( core::pose::Pose const &, std::list< numeric::xyzVector<core::Real> > const & ) {

	bool found_protein;
	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				found_protein = false;
				if ( (typGrid_[ix][iy][iz] == PROTEIN) || (typGrid_[ix][iy][iz] == SURFACE) ) continue; // skip points that are already assigned as protein or surface

				for ( core::Size i = ix+1; i < espGrid_dim_.x(); ++i ) {
					if ( typGrid_[i][iy][iz] == PROTEIN ) {
						typGrid_[ix][iy][iz] = PSEUDO_PROTEIN;
						found_protein = true;
					}
					if ( typGrid_[i][iy][iz] == SURFACE ) break;
					if ( found_protein ) break;
				}
				if ( !found_protein ) {
					for ( core::Size j = iy+1; j < espGrid_dim_.y(); ++j ) {
						if ( typGrid_[ix][j][iz] == PROTEIN ) {
							typGrid_[ix][iy][iz] = PSEUDO_PROTEIN;
							found_protein = true;
						}
						if ( typGrid_[ix][j][iz] == SURFACE ) break;
						if ( found_protein ) break;
					}
				}
				if ( !found_protein ) {
					for ( core::Size k = iz+1; k < espGrid_dim_.z(); ++k ) {
						if ( typGrid_[ix][iy][k] == PROTEIN ) {
							typGrid_[ix][iy][iz] = PSEUDO_PROTEIN;
							found_protein = true;
						}
						if ( typGrid_[ix][iy][k] == SURFACE ) break;
						if ( found_protein ) break;
					}
				}

				if ( !found_protein ) {
					for ( int i = ix-1; i >= 0 ; --i ) {
						if ( typGrid_[i][iy][iz] == PROTEIN ) {
							typGrid_[ix][iy][iz] = PSEUDO_PROTEIN;
							found_protein = true;
						}
						if ( typGrid_[i][iy][iz] == SURFACE ) break;
						if ( found_protein ) break;
					}
				}

				if ( !found_protein ) {
					for ( int j = iy-1; j >= 0; --j ) {
						if ( typGrid_[ix][j][iz] == PROTEIN ) {
							typGrid_[ix][iy][iz] = PSEUDO_PROTEIN;
							found_protein = true;
						}
						if ( typGrid_[ix][j][iz] == SURFACE ) break;
						if ( found_protein ) break;
					}
				}

				if ( !found_protein ) {
					for ( int k = iz-1; k >= 0; --k ) {
						if ( typGrid_[ix][iy][k] == PROTEIN ) {
							typGrid_[ix][iy][iz] = PSEUDO_PROTEIN;
							found_protein = true;
						}
						if ( typGrid_[ix][iy][k] == SURFACE ) break;
						if ( found_protein ) break;
					}
				}

			}
		}
	}

	//Assign PSUDO_PROTIEN points as PROTEIN points
	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {
				if ( typGrid_[ix][iy][iz] == PSEUDO_PROTEIN ) { typGrid_[ix][iy][iz] = PROTEIN;}
			}
		}
	}

	return;

}

void ElectrostaticpotentialGrid::write_espGrid_to_pdb( std::string const & output_pdbname ) const {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_pdbname, std::ios::out);

	for ( Size iz=0; iz<espGrid_dim_.z(); ++iz ) {
		for ( Size iy=0; iy<espGrid_dim_.y(); ++iy ) {
			for ( Size ix=0; ix<espGrid_dim_.x(); ++ix ) {

				core::Real xcoord = ix * espGrid_spacing_ + ( espGrid_mid_.x() - ((static_cast<core::Real>(espGrid_dim_.x()-1)/2) * espGrid_spacing_) );
				core::Real ycoord = iy * espGrid_spacing_ + ( espGrid_mid_.y() - ((static_cast<core::Real>(espGrid_dim_.y()-1)/2) * espGrid_spacing_) );
				core::Real zcoord = iz * espGrid_spacing_ + ( espGrid_mid_.z() - ((static_cast<core::Real>(espGrid_dim_.z()-1)/2) * espGrid_spacing_) );

				if ( (typGrid_[ix][iy][iz] == PROTEIN) || (typGrid_[ix][iy][iz] == PSEUDO_PROTEIN) ) {
					outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   PRO A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<espGrid_[ix][iy][iz]<<std::endl;
				} else if ( typGrid_[ix][iy][iz] == SURFACE ) {
					outPDB_stream<<"HETATM   "<<std::setw(2)<<2<<"  O   SUR A   2    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<espGrid_[ix][iy][iz]<<std::endl;
				} else if ( typGrid_[ix][iy][iz] == SOLVENT ) {
					outPDB_stream<<"HETATM   "<<std::setw(2)<<3<<"  N   SOL A   3    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< xcoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< ycoord <<std::setw(8)<<std::fixed<<std::setprecision(3)<< zcoord <<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<espGrid_[ix][iy][iz]<<std::endl;
				}
			}
		}
	}


	outPDB_stream.close();
	outPDB_stream.clear();

	return;

}

void EggshellGrid::get_connolly_eggshell_on_grid( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list, const PocketGrid& eggGrid, const PocketGrid& extGrid ) {

	bool found_surface = false;
	//assign surface point for eggGrid to get the eggshell points
	eggshell_coord_list_.clear();
	for ( Size iz=0; iz<eggGrid.dim_.z(); ++iz ) {
		for ( Size iy=0; iy<eggGrid.dim_.y(); ++iy ) {
			for ( Size ix=0; ix<eggGrid.dim_.x(); ++ix ) {
				numeric::xyzVector<core::Real> grid_coord(ix,iy,iz);
				numeric::xyzVector<core::Real> grid_cart;
				convert_grid_to_cartesian(grid_coord, eggGrid.mid_, eggGrid.dim_, eggGrid.spacing_, grid_cart);
				found_surface = false;
				for ( auto const & sf : surfacePoints_list ) {
					if ( grid_cart.distance(sf) < 0.5 ) {
						eggshell_coord_list_.push_back(grid_cart);
						found_surface = true;
					}
					if ( found_surface ) break;
				}
			}
		}
	}

	//assign surface point for extGrid (to get extra eggshell points)
	extra_coord_list_.clear();
	for ( Size iz=0; iz<extGrid.zdim_; ++iz ) {
		for ( Size iy=0; iy<extGrid.ydim_; ++iy ) {
			for ( Size ix=0; ix<extGrid.xdim_; ++ix ) {
				numeric::xyzVector<core::Real> grid_coord(ix,iy,iz);
				numeric::xyzVector<core::Real> grid_cart;
				convert_grid_to_cartesian(grid_coord, eggGrid.mid_, eggGrid.dim_, eggGrid.spacing_, grid_cart);
				found_surface = false;
				for ( auto const & sf : surfacePoints_list ) {
					if ( grid_coord.distance(sf) < 0.5 ) {
						extra_coord_list_.push_back(grid_cart);
						found_surface = true;
					}
					if ( found_surface ) break;
				}
			}
		}
	}

	//Now remove 'eggshell point' from extra points
	for ( std::list< numeric::xyzVector<core::Real> >::const_iterator eg = eggshell_coord_list_.begin(); eg != eggshell_coord_list_.end(); ++eg ) {
		for ( auto ex = extra_coord_list_.begin(); ex != extra_coord_list_.end(); ) {
			if ( eg == ex ) {
				ex = extra_coord_list_.erase(ex);
			} else {
				++ex;
			}
		}
	}

	eggshell_CoM_ = calculate_center_xyz(eggshell_coord_list_);

	//write_eggshell_to_pdb( "surf_eggshell.pdb" );

	return;
}

// David, does this look right? Is the right data being compared to make sure two grids are the same?
bool ElectrostaticpotentialGrid::operator == ( PocketGrid const & other ) const {
	if ( ! PocketGrid::operator == ( other ) ) return false;

	ElectrostaticpotentialGrid const & other_downcast( static_cast< ElectrostaticpotentialGrid const & > ( other ));
	if ( espGrid_dim_         != other_downcast.espGrid_dim_         ) return false;
	if ( espGrid_mid_         != other_downcast.espGrid_mid_         ) return false;
	if ( espGrid_spacing_     != other_downcast.espGrid_spacing_     ) return false;
	if ( espGrid_points_list_ != other_downcast.espGrid_points_list_ ) return false;
	return true;
}


bool ElectrostaticpotentialGrid::same_type_as_me( PocketGrid const & other ) const
{
	ElectrostaticpotentialGrid const * casted_other = dynamic_cast< ElectrostaticpotentialGrid const * > (&other);
	return casted_other;
}


void EggshellGrid::get_connolly_eggshell( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list, const PocketGrid& eggGrid, const PocketGrid& extGrid ) {

	numeric::xyzVector<core::Real> grid_coord;

	//assign connolly surface point for eggGrid
	eggshell_coord_list_.clear();
	for ( auto const & sf : surfacePoints_list ) {
		convert_cartesian_to_grid(sf, eggGrid.mid_, eggGrid.dim_, eggGrid.spacing_, grid_coord);
		//ignore if the surface point goes out of the grid
		//std::cout<<" GRID SURF : "<<grid_coord.x()<<" "<<grid_coord.y()<<" "<<grid_coord.z()<<std::endl;
		//std::cout<<" GRID DIM : "<<eggGrid.dim_.x()<<" "<<eggGrid.dim_.y()<<" "<<eggGrid.dim_.z()<<std::endl;
		//std::cout<<" GRID MID : "<<eggGrid.mid_.x()<<" "<<eggGrid.mid_.y()<<" "<<eggGrid.mid_.z()<<std::endl;
		if ( ((grid_coord.x() >= eggGrid.dim_.x())||(grid_coord.x() < 0)) || ((grid_coord.y() >= eggGrid.dim_.y())||(grid_coord.y() < 0)) || ((grid_coord.z() >= eggGrid.dim_.z())||(grid_coord.z() < 0)) ) {

			continue;
		} else {
			eggshell_coord_list_.push_back(sf);
		}
	}

	//assign connolly surface point for extGrid
	extra_coord_list_.clear();
	for ( auto const & sf : surfacePoints_list ) {
		convert_cartesian_to_grid(sf, extGrid.mid_, extGrid.dim_, extGrid.spacing_, grid_coord);
		//ignore if the surface point goes out of the grid
		if ( ((grid_coord.x() >= extGrid.dim_.x())||(grid_coord.x() < 0)) || ((grid_coord.y() >= extGrid.dim_.y())||(grid_coord.y() < 0)) || ((grid_coord.z() >= extGrid.dim_.z())||(grid_coord.z() < 0)) ) {
			continue;
		} else {
			extra_coord_list_.push_back(sf);
		}
	}

	eggshell_CoM_ = calculate_center_xyz(eggshell_coord_list_);

	write_eggshell_to_pdb( "surf_eggshell.pdb" );

	return;

}

EggshellGrid::EggshellGrid( core::pose::Pose const & protein_pose, const PocketGrid& eggGrid, const PocketGrid& extGrid ) {

	std::list< numeric::xyzVector<core::Real> > surfacePoints_list;
	surfacePoints_list = eggGrid.get_connolly_surfacePoints( protein_pose );
	//get_connolly_eggshell_on_grid(surfacePoints_list, eggGrid, extGrid);
	get_connolly_eggshell(surfacePoints_list, eggGrid, extGrid);

	return;

}

EggshellGrid::EggshellGrid( const PocketGrid& ext_grd, std::list< numeric::xyzVector<core::Real> > const & eggshell_coord_list):
	PocketGrid( ext_grd)
{
	init();
	extra_coord_list_.clear();
	using namespace basic::options;
	core::Real const extra_dist = option[ OptionKeys::pocket_grid::extra_eggshell_dist ];
	core::Size searchxmin, searchxmax, searchymin, searchymax, searchzmin, searchzmax, xx,yy,zz;
	numeric::xyzVector<core::Real> extra_coord;
	for ( core::Size x=0; x<(ext_grd.xdim_); x++ ) {
		for ( core::Size y=0; y<(ext_grd.ydim_); y++ ) {
			for ( core::Size z=0; z<(ext_grd.zdim_); z++ ) {
				if ( (ext_grd.grid_[x][y][z] == ext_grd.HSURFACE) || (ext_grd.grid_[x][y][z] == ext_grd.PSURFACE) || (ext_grd.grid_[x][y][z] == ext_grd.T_SURFACE) || (ext_grd.grid_[x][y][z] == ext_grd.ST_SURFACE) ) {
					if ( x != 0 ) {
						searchxmin=x-1;
					} else {
						searchxmin=x;
					}
					if ( x != ext_grd.xdim_-1 ) {
						searchxmax=x+1;
					} else {
						searchxmax=x;
					}
					if ( y != 0 ) {
						searchymin=y-1;
					} else {
						searchymin=y;
					}
					if ( y != ext_grd.ydim_-1 ) {
						searchymax=y+1;
					} else {
						searchymax=y;
					}
					if ( z != 0 ) {
						searchzmin=z-1;
					} else {
						searchzmin=z;
					}
					if ( z != ext_grd.zdim_-1 ) {
						searchzmax=z+1;
					} else {
						searchzmax=z;
					}
					bool found = false;
					for ( xx = searchxmin; xx<=searchxmax; ++xx ) {
						for ( yy = searchymin; yy<=searchymax; ++yy ) {
							for ( zz = searchzmin; zz<=searchzmax; ++zz ) {
								if ( xx==x && yy==y && zz==z ) continue;
								if ( ( ext_grd.grid_[xx][yy][zz]==ext_grd.EMPTY ) || (ext_grd.grid_[xx][yy][zz] == ext_grd.TP_POCKET) || (ext_grd.grid_[xx][yy][zz] == ext_grd.TP_EDGE) || (ext_grd.grid_[xx][yy][zz] == ext_grd.TP_BURIED) || (ext_grd.grid_[xx][yy][zz] == ext_grd.TP_SURF) || (ext_grd.grid_[xx][yy][zz] == ext_grd.POCKET) || (ext_grd.grid_[xx][yy][zz] == ext_grd.PO_BURIED) || (ext_grd.grid_[xx][yy][zz] == ext_grd.PO_SURF) || (ext_grd.grid_[xx][yy][zz] == ext_grd.PO_EDGE) ) {
									extra_coord.x() = (x*ext_grd.stepSize_+ext_grd.xcorn_);
									extra_coord.y() = (y*ext_grd.stepSize_+ext_grd.ycorn_);
									extra_coord.z() = (z*ext_grd.stepSize_+ext_grd.zcorn_);
									found = true;
									bool nearby = false;
									for ( auto const & eg : eggshell_coord_list ) {
										if ( extra_coord == eg ) {
											nearby = false;
											break;
										}
										if ( extra_coord.distance(eg) < extra_dist ) {
											nearby = true;
										}
									}
									if ( nearby ) {
										extra_coord_list_.push_back(extra_coord);
									}
								}
								if ( found ) break;
							}
							if ( found ) break;
						}
						if ( found ) break;
					}
				}
			}
		}
	}

	//Print grid points into a PDB file
	using namespace basic::options;
	if ( option[ OptionKeys::pocket_grid::dump_Grid ]() ) {
		write_Grid_to_pdb(ext_grd, "ext_grd.pdb");
	}
	return;

}

EggshellGrid::EggshellGrid( const PocketGrid& gr ) :
	// initialize member data using PocketGrid
	PocketGrid( gr )
{
	init();

	eggshell_coord_list_.clear();
	// use Grid to fill in Eggshell grid
	numeric::xyzVector<core::Real> grid_coord;
	core::Size searchxmin, searchxmax, searchymin, searchymax, searchzmin, searchzmax, xx,yy,zz;
	for ( core::Size x=0; x<(gr.xdim_); x++ ) {
		for ( core::Size y=0; y<(gr.ydim_); y++ ) {
			for ( core::Size z=0; z<(gr.zdim_); z++ ) {
				if ( (gr.grid_[x][y][z] == gr.HSURFACE) || (gr.grid_[x][y][z] == gr.PSURFACE) || (gr.grid_[x][y][z] == gr.T_SURFACE)| (gr.grid_[x][y][z] == gr.ST_SURFACE) /*||(gr.grid_[x][y][z]==gr.TP_SURF)*/ ) {
					if ( x != 0 ) {
						searchxmin=x-1;
					} else {
						searchxmin=x;
					}
					if ( x != gr.xdim_-1 ) {
						searchxmax=x+1;
					} else {
						searchxmax=x;
					}
					if ( y != 0 ) {
						searchymin=y-1;
					} else {
						searchymin=y;
					}
					if ( y != gr.ydim_-1 ) {
						searchymax=y+1;
					} else {
						searchymax=y;
					}
					if ( z != 0 ) {
						searchzmin=z-1;
					} else {
						searchzmin=z;
					}
					if ( z != gr.zdim_-1 ) {
						searchzmax=z+1;
					} else {
						searchzmax=z;
					}
					bool found=false;
					for ( xx = searchxmin; (xx<=searchxmax) && ! found; xx++ ) {
						for ( yy = searchymin; (yy<=searchymax) && ! found; yy++ ) {
							for ( zz = searchzmin; (zz<=searchzmax) && ! found; zz++ ) {
								if ( xx==x && yy==y && zz==z ) continue;
								if ( (gr.grid_[xx][yy][zz] == gr.TP_POCKET) || (gr.grid_[xx][yy][zz] == gr.TP_EDGE) || (gr.grid_[xx][yy][zz] == gr.TP_BURIED) ) {
									found = true;
									grid_[x][y][z] = EGGSHELL;
									grid_coord.x() = (x*gr.stepSize_+gr.xcorn_);
									grid_coord.y() = (y*gr.stepSize_+gr.ycorn_);
									grid_coord.z() = (z*gr.stepSize_+gr.zcorn_);
									eggshell_coord_list_.push_back(grid_coord);
								}
							}
						}
					}
				}
			}
		}
	}

	//set eggshell CoM from eggshell_coord_list
	eggshell_CoM_ = calculate_center_xyz(eggshell_coord_list_);

	//Print grid points into a PDB file
	using namespace basic::options;
	if ( option[ OptionKeys::pocket_grid::dump_Grid ]() ) {
		write_Grid_to_pdb(gr, "egg_grd.pdb");
	}
	return;

}

EggshellGrid::EggshellGrid( const PocketGrid& grd, core::pose::Pose const & ligand_pose ) :
	// initialize member data using PocketGrid
	PocketGrid( grd )
{
	init();
	using namespace basic::options;
	core::Real const egg_dist = option[ OptionKeys::pocket_grid::eggshell_dist ];
	core::Real const ext_dist = option[ OptionKeys::pocket_grid::extra_eggshell_dist ];

	std::list< numeric::xyzVector<core::Real> > grd_coord_list;
	// use Grid to fill in Eggshell grid
	numeric::xyzVector<core::Real> grd_coord;
	//core::Size searchxmin, searchxmax, searchymin, searchymax, searchzmin, searchzmax, xx,yy,zz;
	//bool found;

	for ( core::Size x=0; x<(grd.xdim_); x++ ) {
		for ( core::Size y=0; y<(grd.ydim_); y++ ) {
			for ( core::Size z=0; z<(grd.zdim_); z++ ) {
				//found = false;
				if ( (grd.grid_[x][y][z] == grd.HSURFACE) || (grd.grid_[x][y][z] == grd.PSURFACE) || (grd.grid_[x][y][z] == grd.T_SURFACE) || (grd.grid_[x][y][z] == grd.ST_SURFACE) ) {
					grd_coord.x() = (x*grd.stepSize_+grd.xcorn_);
					grd_coord.y() = (y*grd.stepSize_+grd.ycorn_);
					grd_coord.z() = (z*grd.stepSize_+grd.zcorn_);
					//found = true;
					grd_coord_list.push_back(grd_coord);
				}
			}
		}
	}

	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = ligand_pose.size(); j <= resnum; ++j ) {
		if ( !ligand_pose.residue(j).is_protein() ) {
			lig_res_num = j;
			break;
		}
	}
	if ( lig_res_num == 0 ) {
		TR.Fatal<<"Error, no ligand to include_eggshell_points_based_on_known_ligand" << std::endl;
		exit(1);
	}

	core::conformation::Residue const & curr_rsd = ligand_pose.conformation().residue(lig_res_num);
	core::Size ligand_total_atoms = curr_rsd.nheavyatoms();
	numeric::xyzVector<core::Real> lig_atom_coord;
	std::list< numeric::xyzVector<core::Real> > lig_atom_coord_list;
	for ( Size i = 1, i_end = ligand_total_atoms; i <= i_end; ++i ) {
		lig_atom_coord.x() = curr_rsd.atom(i).xyz()(1);
		lig_atom_coord.y() = curr_rsd.atom(i).xyz()(2);
		lig_atom_coord.z() = curr_rsd.atom(i).xyz()(3);
		lig_atom_coord_list.push_back(lig_atom_coord);
	}

	numeric::xyzVector<core::Real> xyz_coord;
	for ( std::list< numeric::xyzVector<core::Real> >::const_iterator aa = grd_coord_list.begin(); aa != grd_coord_list.end(); ++aa ) {
		xyz_coord = *aa;
		bool found_egg = false;
		for ( std::list< numeric::xyzVector<core::Real> >::const_iterator bb = lig_atom_coord_list.begin(); bb != lig_atom_coord_list.end(); ++bb ) {
			if ( xyz_coord.distance(*bb) <= egg_dist ) {
				found_egg = true;break;
			}
		}
		if ( found_egg ) eggshell_coord_list_.push_back(xyz_coord);
	}

	for ( std::list< numeric::xyzVector<core::Real> >::const_iterator aa = grd_coord_list.begin(); aa != grd_coord_list.end(); ++aa ) {
		xyz_coord = *aa;
		bool found_ext = false;
		for ( std::list< numeric::xyzVector<core::Real> >::const_iterator bb = eggshell_coord_list_.begin(); bb != eggshell_coord_list_.end(); ++bb ) {
			if ( *bb == *aa ) continue;
			if ( xyz_coord.distance(*bb) <= ext_dist ) {
				found_ext = true;break;
			}
		}
		if ( found_ext ) extra_coord_list_.push_back(xyz_coord);
	}

	//set eggshell CoM from eggshell_coord_list
	eggshell_CoM_ = calculate_center_xyz(eggshell_coord_list_);
	//Print grid points into a PDB file
	using namespace basic::options;
	if ( option[ OptionKeys::pocket_grid::dump_Grid ]() ) {
		write_Grid_to_pdb(grd, "grd.pdb");
	}
	return;
}

EggshellGrid::EggshellGrid( std::string const & fname ) {

	// NOTE: THIS FUNCTION DOESN'T FULLY SETUP THE GRID OR extra_coord_list_
	// HERE WE ONLY PROVIDE MINIMAL FUNCTIONALITY TO DO EGGSHELL VS EGGSHELL COMPARISONS....

	utility::io::izstream instream;
	instream.open(fname, std::ios::in);

	// Read stepSize_ xdim_ ydim_ zdim_ xcorn_ ycorn_ zcorn_
	instream >> stepSize_ >> xdim_ >> ydim_ >> zdim_ >> xcorn_ >> ycorn_ >> zcorn_;

	// Reset the grid for the new dimensions
	init();
	eggshell_coord_list_.clear();
	extra_coord_list_.clear();

	// Read the grid points which should be marked EGGSHELL, use these to push onto eggshell_coord_list_
	// Format is x y z
	core::Size x, y, z;
	while ( ! instream.eof() ) {
		instream >> x >> y >> z;
		grid_[x][y][z] = EGGSHELL;
		numeric::xyzVector<core::Real> grid_coord;
		grid_coord.x() = x*stepSize_+xcorn_;
		grid_coord.y() = y*stepSize_+ycorn_;
		grid_coord.z() = z*stepSize_+zcorn_;
		eggshell_coord_list_.push_back(grid_coord);
	}

	instream.close();
	instream.clear();

}


void EggshellGrid::dump_eggshell( std::string const & fname ) const {

	// NOTE: THIS FUNCTION DOESN'T WRITE THE WHOLE GRID OR extra_coord_list_
	//         HERE WE ONLY PROVIDE MINIMAL FUNCTIONALITY TO DO EGGSHELL VS EGGSHELL COMPARISONS....

	utility::io::ozstream outstream;
	outstream.open(fname, std::ios::out);

	// Print stepSize_ xdim_ ydim_ zdim_ xcorn_ ycorn_ zcorn_
	outstream << stepSize_ << " " << xdim_ << " " << ydim_ << " " << zdim_ << " " << xcorn_ << " " << ycorn_ << " " << zcorn_ << std::endl;

	// Need to print each grid point if it's of type EGGSHELL along with the corresponding cartesian coors ( for eggshell_coord_list_ )
	for ( core::Size x=0; x<xdim_; x++ ) {
		for ( core::Size y=0; y<ydim_; y++ ) {
			for ( core::Size z=0; z<zdim_; z++ ) {
				if ( ( grid_[x][y][z]==EGGSHELL ) ) {
					outstream << x << " " << y << " " << z << std::endl;
					//     outstream << x << " " << y << " " << z << " " << x*stepSize_+xcorn_ << " " << y*stepSize_+ycorn_ << " " << z*stepSize_+zcorn_ << std::endl;
				}
			}
		}
	}

	outstream.close();
	outstream.clear();

}

void EggshellGrid::dump_eggshell( std::string const & fname, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 ) const {

	core::Real xmin = 1,xmax=0,ymin=1,ymax=0,zmin=1,zmax=0;
	for ( auto const & pit : eggshell_coord_list_ ) {
		numeric::xyzVector<core::Real> coord(pit.x(),pit.y(),pit.z());
		coord = rot1 * coord;
		coord = rot2 * coord;
		coord = rot3 * coord;
		if ( xmin>xmax ) {
			xmin = coord(1);
			xmax = coord(1);
			ymin = coord(2);
			ymax = coord(2);
			zmin = coord(3);
			zmax = coord(3);
		} else {
			if ( xmin > coord(1) ) xmin = coord(1);
			if ( xmax < coord(1) ) xmax = coord(1);
			if ( ymin > coord(2) ) ymin = coord(2);
			if ( ymax < coord(2) ) ymax = coord(2);
			if ( zmin > coord(3) ) zmin = coord(3);
			if ( zmax < coord(3) ) zmax = coord(3);
		}
	}

	core::Real new_xcorn=(core::Real)std::floor(xmin);
	core::Real new_ycorn=(core::Real)std::floor(ymin);
	core::Real new_zcorn=(core::Real)std::floor(zmin);
	core::Real x_far_corn=(core::Real)std::ceil(xmax);
	core::Real y_far_corn=(core::Real)std::ceil(ymax);
	core::Real z_far_corn=(core::Real)std::ceil(zmax);
	core::Size new_xdim = (core::Size)std::ceil((x_far_corn-new_xcorn)/stepSize_)+1;
	core::Size new_ydim = (core::Size)std::ceil((y_far_corn-new_ycorn)/stepSize_)+1;
	core::Size new_zdim = (core::Size)std::ceil((z_far_corn-new_zcorn)/stepSize_)+1;

	// NOTE: THIS FUNCTION DOESN'T WRITE THE WHOLE GRID_
	// NOTE: THIS FUNCTION DOESN'T WRITE THE WHOLE GRID OR extra_coord_list_
	//         HERE WE ONLY PROVIDE MINIMAL FUNCTIONALITY TO DO EGGSHELL VS EGGSHELL COMPARISONS....

	utility::io::ozstream outstream;
	outstream.open(fname, std::ios::out);

	// Print stepSize_ xdim_ ydim_ zdim_ xcorn_ ycorn_ zcorn_
	outstream << stepSize_ << " " << new_xdim << " " << new_ydim << " " << new_zdim << " " << new_xcorn << " " << new_ycorn << " " << new_zcorn << std::endl;

	// Need to print each grid point if it's of type EGGSHELL along with the corresponding cartesian coors ( for eggshell_coord_list_ )
	for ( core::Size x=0; x<xdim_; x++ ) {
		for ( core::Size y=0; y<ydim_; y++ ) {
			for ( core::Size z=0; z<zdim_; z++ ) {
				if ( ( grid_[x][y][z]==EGGSHELL ) ) {
					numeric::xyzVector<core::Real> coord(xcorn_ + x * stepSize_, ycorn_ + y * stepSize_, zcorn_ + z * stepSize_);
					coord = rot1 * coord;
					coord = rot2 * coord;
					coord = rot3 * coord;
					numeric::xyzVector<core::Size> newCoord(0,0,0);
					newCoord(1) = (core::Size)std::floor((coord(1) - new_xcorn)/stepSize_ + .5);
					newCoord(2) = (core::Size)std::floor((coord(2) - new_ycorn)/stepSize_ + .5);
					newCoord(3) = (core::Size)std::floor((coord(3) - new_zcorn)/stepSize_ + .5);

					outstream << newCoord(1) << " " << newCoord(2) << " " << newCoord(3) << std::endl;
					//          outstream << x << " " << y << " " << z << " " << x*stepSize_+xcorn_ << " " << y*stepSize_+ycorn_ << " " << z*stepSize_+zcorn_ << std::endl;
				}
			}
		}
	}

	outstream.close();
	outstream.clear();

}

void EggshellGrid::write_eggshell_to_pdb( std::string const & output_pdbname ) const {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_pdbname, std::ios::out);
	for ( auto const & pd : eggshell_coord_list_ ) {
		outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   EGG A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.z()<<std::endl;
	}
	for ( auto const & pd : extra_coord_list_ ) {
		outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  O   EXT B   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.z()<<std::endl;
	}
	outPDB_stream.close();
	outPDB_stream.clear();

}

void EggshellGrid::write_eggshell_to_pdb( std::string const & output_pdbname, numeric::xyzMatrix<core::Real> rot1, numeric::xyzMatrix<core::Real> rot2, numeric::xyzMatrix<core::Real> rot3 ) const {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_pdbname, std::ios::out);
	for ( auto const & pd : eggshell_coord_list_ ) {
		numeric::xyzVector<core::Real> coord(pd.x(),pd.y(),pd.z());
		coord = rot1 * coord;
		coord = rot2 * coord;
		coord = rot3 * coord;
		outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   EGG A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<coord(1)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<coord(2)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<coord(3)<<std::endl;
	}
	for ( auto const & pd : extra_coord_list_ ) {
		numeric::xyzVector<core::Real> coord(pd.x(),pd.y(),pd.z());
		coord = rot1 * coord;
		coord = rot2 * coord;
		coord = rot3 * coord;
		outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  O   EGG B   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<coord(1)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<coord(2)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<coord(3)<<std::endl;
	}
	outPDB_stream.close();
	outPDB_stream.clear();

}

void EggshellGrid::write_Grid_to_pdb( const PocketGrid& gr, std::string const & output_pdbname ) const {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_pdbname, std::ios::out);

	for ( core::Size x=0; x<(gr.xdim_); x++ ) {
		for ( core::Size y=0; y<(gr.ydim_); y++ ) {
			for ( core::Size z=0; z<(gr.zdim_); z++ ) {
				outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x*gr.stepSize_+gr.xcorn_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y*gr.stepSize_+gr.ycorn_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z*gr.stepSize_+gr.zcorn_<<std::endl;
			}
		}
	}
	outPDB_stream.close();
	outPDB_stream.clear();

	return;

}

numeric::xyzVector<core::Real> EggshellGrid::calculate_center_xyz(std::list< numeric::xyzVector<core::Real> > const & pocketshell_coord_list ) {

	numeric::xyzVector<core::Real> pocketshell_CoM(0.);
	for ( auto const & pd : pocketshell_coord_list ) {
		pocketshell_CoM.x() += pd.x();
		pocketshell_CoM.y() += pd.y();
		pocketshell_CoM.z() += pd.z();
	}
	pocketshell_CoM /= pocketshell_coord_list.size();
	return pocketshell_CoM;
}

core::Real EggshellGrid::get_eggshell_distance( EggshellGrid const & template_eggshell, std::string const & comparison_pdbname ) const {

	TR << "My CoM is " << eggshell_CoM_.x() << " "  << eggshell_CoM_.y() << " "  << eggshell_CoM_.z() << std::endl;
	TR << "Template CoM is " << template_eggshell.eggshell_CoM_.x() << " "  << template_eggshell.eggshell_CoM_.y() << " "  << template_eggshell.eggshell_CoM_.z() << std::endl;

	bool write_comparison_pdb = false;
	// note: string comparison will return non-zero if the strings are NOT equal
	if ( comparison_pdbname.compare("") ) {
		write_comparison_pdb = true;
	}

	utility::io::ozstream comparisonPDB_stream;
	if ( write_comparison_pdb ) comparisonPDB_stream.open(comparison_pdbname, std::ios::out);

	core::Real const match_weight = 1.;
	core::Real const mismatch_weight = 0.;
	core::Real match_score = 0.;
	core::Size output_res_num = 1;
	core::Size common_points = 0;

	// Loop over all points in the template eggshell
	for ( auto const & pd : template_eggshell.eggshell_coord_list() ) {

		// For this template eggshell point, get the cartesian coors and convert to comparison's grid
		core::Real const template_x = pd.x();
		core::Real const template_y = pd.y();
		core::Real const template_z = pd.z();
		core::Size const self_x_index = (core::Size) floor( ( ( template_x - xcorn_ ) / stepSize_ ) + 0.5 );
		core::Size const self_y_index = (core::Size) floor( ( ( template_y - ycorn_ ) / stepSize_ ) + 0.5 );
		core::Size const self_z_index = (core::Size) floor( ( ( template_z - zcorn_ ) / stepSize_ ) + 0.5 );

		//Check to see if template point is within comparison's grid range
		if ( ( self_x_index < xdim_ ) && ( self_y_index < ydim_ ) && ( self_z_index < zdim_ ) ) {
			// note: no need to check that self_index >= 0 since it's of type "Size"
			//  if ( ( self_x_index >= 0 ) && ( self_x_index < xdim_ ) && ( self_y_index >= 0 ) && ( self_y_index < ydim_ ) && ( self_z_index >= 0 ) && ( self_z_index < zdim_ ) ) {
			//See if template point is eggshell for comparison
			if ( grid_[self_x_index][self_y_index][self_z_index] == EGGSHELL ) {
				++common_points;
				match_score += match_weight;
				// KK output_res_num needed
				if ( write_comparison_pdb ) comparisonPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   MEG A   "<<std::setw(8)<<output_res_num<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_x<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_y<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_z<<std::endl;

			} else {
				match_score += mismatch_weight;
				// KK output_res_num needed
				if ( write_comparison_pdb ) comparisonPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   NEG A   "<<std::setw(8)<<output_res_num<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_x<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_y<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_z<<std::endl;
			}
		} else {
			//Score as mismatch if it's out of comparison's grid
			match_score += mismatch_weight;
			// KK output_res_num needed
			if ( write_comparison_pdb ) comparisonPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   NEG A   "<<std::setw(8)<<output_res_num<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_x<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_y<<std::setw(8)<<std::fixed<<std::setprecision(3)<<template_z<<std::endl;
		}
		output_res_num++;
	}

	if ( write_comparison_pdb ) {
		comparisonPDB_stream.close();
		comparisonPDB_stream.clear();
	}

	match_score = match_score / ( template_eggshell.eggshell_coord_list().size() + eggshell_coord_list().size() - common_points);
	TR << "Eggshell match score = " << match_score << std::endl;

	return 1. - match_score;
}

// David, does this look right? Is the right data being compared to make sure two grids are the same?
bool EggshellGrid::operator == ( PocketGrid const & other ) const {
	if ( ! PocketGrid::operator == ( other ) ) return false;
	EggshellGrid const & other_downcast( static_cast< EggshellGrid const & > ( other ));
	if ( eggshell_coord_list_ != other_downcast.eggshell_coord_list_ ) return false;
	if ( extra_coord_list_    != other_downcast.extra_coord_list_    ) return false;
	if ( eggshell_CoM_        != other_downcast.eggshell_CoM_        ) return false;
	return true;
}

bool EggshellGrid::same_type_as_me( PocketGrid const & other ) const
{
	EggshellGrid const * casted_other = dynamic_cast< EggshellGrid const * > (&other);
	return casted_other;
}


ComparisonGrid::ComparisonGrid(const PocketGrid& gr){
	xdim_ = gr.xdim_;
	xcorn_ = gr.xcorn_;
	ydim_ = gr.ydim_;
	ycorn_ = gr.ycorn_;
	zdim_ = gr.zdim_;
	zcorn_ = gr.zcorn_;
	stepSize_ = gr.stepSize_;
	grid_.resize(xdim_);
	for ( core::Size tx=0; tx<xdim_; tx++ ) {
		grid_[tx].resize(ydim_);
		for ( core::Size ty=0; ty<ydim_; ty++ ) {
			grid_[tx][ty].resize(zdim_, EMPTY);
		}
	}
}

core::Real ComparisonGrid::mark (const PocketGrid& gr, core::Real x, core::Real y, core::Real z, core::Real const & vdWd, core::Real const & penalty){
	core::Real return_penalty = 0.;

	x-=xcorn_;
	y-=ycorn_;
	z-=zcorn_;
	int xcen=(int)floor(x/stepSize_+.5);
	int ycen=(int)floor(y/stepSize_+.5);
	int zcen=(int)floor(z/stepSize_+.5);
	int minX,maxX;
	int minY,maxY;
	int minZ,maxZ;
	core::Real vdW=vdWd/stepSize_;

	minX=int(std::max(ceil(((core::Real)(xcen)-vdW)-.5), 0.));
	maxX=int(std::min(floor(((core::Real)(xcen)+vdW)+.5), (core::Real)(xdim_)-1.));
	minY=int(std::max(ceil(((core::Real)(ycen)-vdW)-.5), 0.));
	maxY=int(std::min(floor(((core::Real)(ycen)+vdW)+.5), (core::Real)(ydim_)-1.));
	minZ=int(std::max(ceil(((core::Real)(zcen)-vdW)-.5), 0.));
	maxZ=int(std::min(floor(((core::Real)(zcen)+vdW)+.5), (core::Real)(zdim_)-1.));
	core::Real volume=0;
	core::Real pocketVolume=0;

	for ( int xIter=minX; xIter<=maxX; ++xIter ) {
		for ( int yIter=minY; yIter<=maxY; ++yIter ) {
			int centerZ=std::max(minZ, zcen);
			for ( int zIter=centerZ; zIter <= maxZ; ++zIter ) {
				if ( pow((xIter-xcen),2)+pow((yIter-ycen),2)+pow((zIter-zcen),2)<=pow(vdW, 2) ) {
					volume += 1;
					if ( gr.grid_[xIter][yIter][zIter]==gr.PROTEIN || gr.grid_[xIter][yIter][zIter]==gr.TARGET || gr.grid_[xIter][yIter][zIter]==gr.SUBTARGET ) {
						return_penalty += penalty;
					} else {
						if ( gr.grid_[xIter][yIter][zIter]==gr.TP_POCKET ||gr.grid_[xIter][yIter][zIter]==gr.TP_BURIED || gr.grid_[xIter][yIter][zIter]==gr.TP_EDGE ) pocketVolume+=1;
						grid_[xIter][yIter][zIter]=LIGAND;
					}
				}
			}
		}
	}
	//penalize for not being in the pocket
	return_penalty += volume * penalty * (1 - std::min(pocketVolume*2,volume)/volume);
	return return_penalty;

}

core::Real ComparisonGrid::compareCoverage(const PocketGrid& gr){
	core::Real return_penalty = 0;
	for ( auto cit=gr.clusters_.clusters_.begin(); cit != gr.clusters_.clusters_.end(); ++cit ) {
		if ( cit->points_.size()*pow(gr.stepSize_,3)<gr.minPockSize_ ) continue; // this is a smallpocket
		for ( auto const & point : cit->points_ ) {
			if ( (gr.grid_[point.x][point.y][point.z]==gr.TP_POCKET) || (gr.grid_[point.x][point.y][point.z]==gr.TP_BURIED ||
					gr.grid_[point.x][point.y][point.z]==gr.TP_EDGE) ) {
				if ( grid_[point.x][point.y][point.z]==EMPTY ) return_penalty +=1.;
			}
		}
	}
	return return_penalty;
}

} // pockets
} // protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::PCluster::Cxyz::save( Archive & arc ) const {
	arc( CEREAL_NVP( x ) ); // core::Size
	arc( CEREAL_NVP( y ) ); // core::Size
	arc( CEREAL_NVP( z ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::PCluster::Cxyz::load( Archive & arc ) {
	arc( x ); // core::Size
	arc( y ); // core::Size
	arc( z ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::PCluster::Cxyz );

/// @brief Default constructor required by cereal to deserialize this class
protocols::pockets::EggshellGrid::EggshellGrid() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::EggshellGrid::save( Archive & arc ) const {
	arc( cereal::base_class< PocketGrid >( this ) );
	arc( CEREAL_NVP( eggshell_coord_list_ ) ); // std::list<numeric::xyzVector<core::Real> >
	arc( CEREAL_NVP( extra_coord_list_ ) ); // std::list<numeric::xyzVector<core::Real> >
	arc( CEREAL_NVP( eggshell_CoM_ ) ); // numeric::xyzVector<core::Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::EggshellGrid::load( Archive & arc ) {
	arc( cereal::base_class< PocketGrid >( this ) );
	arc( eggshell_coord_list_ ); // std::list<numeric::xyzVector<core::Real> >
	arc( extra_coord_list_ ); // std::list<numeric::xyzVector<core::Real> >
	arc( eggshell_CoM_ ); // numeric::xyzVector<core::Real>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::EggshellGrid );
CEREAL_REGISTER_TYPE( protocols::pockets::EggshellGrid )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::PocketGrid::save( Archive & arc ) const {
	arc( CEREAL_NVP( grid_ ) ); // std::vector<std::vector<std::vector<PtType> > >
	arc( CEREAL_NVP( c_grid_ ) ); // std::vector<std::vector<std::vector<PtType> > >
	arc( CEREAL_NVP( pockets_ ) ); // std::vector<std::vector<std::vector<core::Size> > >
	arc( CEREAL_NVP( xdim_ ) ); // core::Size
	arc( CEREAL_NVP( ydim_ ) ); // core::Size
	arc( CEREAL_NVP( zdim_ ) ); // core::Size
	arc( CEREAL_NVP( xcorn_ ) ); // core::Real
	arc( CEREAL_NVP( ycorn_ ) ); // core::Real
	arc( CEREAL_NVP( zcorn_ ) ); // core::Real
	arc( CEREAL_NVP( spacing_ ) ); // core::Real
	arc( CEREAL_NVP( stepSize_ ) ); // core::Real
	arc( CEREAL_NVP( mid_ ) ); // numeric::xyzVector<core::Real>
	arc( CEREAL_NVP( dim_ ) ); // numeric::xyzVector<core::Size>
	arc( CEREAL_NVP( maxLen_ ) ); // core::Real
	arc( CEREAL_NVP( tag_ ) ); // std::string
	arc( CEREAL_NVP( pdbno_ ) ); // core::Size
	arc( CEREAL_NVP( expdbno_ ) ); // core::Size
	arc( CEREAL_NVP( numTargets_ ) ); // core::Size
	arc( CEREAL_NVP( clusters_ ) ); // class protocols::pockets::PClusterSet
	arc( CEREAL_NVP( c_clusters_ ) ); // class protocols::pockets::CClusterSet
	arc( CEREAL_NVP( size_x_ ) ); // core::Real
	arc( CEREAL_NVP( size_y_ ) ); // core::Real
	arc( CEREAL_NVP( size_z_ ) ); // core::Real
	arc( CEREAL_NVP( limit_x_ ) ); // core::Real
	arc( CEREAL_NVP( limit_y_ ) ); // core::Real
	arc( CEREAL_NVP( limit_z_ ) ); // core::Real
	arc( CEREAL_NVP( restrictSize_ ) ); // _Bool
	arc( CEREAL_NVP( ignoreBuriedPockets_ ) ); // _Bool
	arc( CEREAL_NVP( ignoreExposedPockets_ ) ); // _Bool
	arc( CEREAL_NVP( probe_rad_ ) ); // core::Real
	arc( CEREAL_NVP( surf_score_ ) ); // core::Real
	arc( CEREAL_NVP( surf_dist_ ) ); // core::Real
	arc( CEREAL_NVP( bur_score_ ) ); // core::Real
	arc( CEREAL_NVP( bur_dist_ ) ); // core::Real
	arc( CEREAL_NVP( side_chains_only_ ) ); // _Bool
	arc( CEREAL_NVP( markpsp_ ) ); // _Bool
	arc( CEREAL_NVP( marksps_ ) ); // _Bool
	arc( CEREAL_NVP( exemplarRestriction_ ) ); // _Bool
    arc( CEREAL_NVP( dumpExemplars_ ) ); // _Bool
    arc( CEREAL_NVP( filterExemplars_ ) ); // _Bool
    arc( CEREAL_NVP( limitExemplars_ ) ); // _Bool
	arc( CEREAL_NVP( search13_ ) ); // _Bool
	arc( CEREAL_NVP( minPockSize_ ) ); // core::Real
	arc( CEREAL_NVP( maxPockSize_ ) ); // core::Real
	arc( CEREAL_NVP( rot_mat_ ) ); // numeric::xyzMatrix<core::Real>
	arc( CEREAL_NVP( static_grid_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::PocketGrid::load( Archive & arc ) {
	arc( grid_ ); // std::vector<std::vector<std::vector<PtType> > >
	arc( c_grid_ ); // std::vector<std::vector<std::vector<PtType> > >
	arc( pockets_ ); // std::vector<std::vector<std::vector<core::Size> > >
	arc( xdim_ ); // core::Size
	arc( ydim_ ); // core::Size
	arc( zdim_ ); // core::Size
	arc( xcorn_ ); // core::Real
	arc( ycorn_ ); // core::Real
	arc( zcorn_ ); // core::Real
	arc( spacing_ ); // core::Real
	arc( stepSize_ ); // core::Real
	arc( mid_ ); // numeric::xyzVector<core::Real>
	arc( dim_ ); // numeric::xyzVector<core::Size>
	arc( maxLen_ ); // core::Real
	arc( tag_ ); // std::string
	arc( pdbno_ ); // core::Size
	arc( expdbno_ ); // core::Size
	arc( numTargets_ ); // core::Size
	arc( clusters_ ); // class protocols::pockets::PClusterSet
	arc( c_clusters_ ); // class protocols::pockets::CClusterSet
	arc( size_x_ ); // core::Real
	arc( size_y_ ); // core::Real
	arc( size_z_ ); // core::Real
	arc( limit_x_ ); // core::Real
	arc( limit_y_ ); // core::Real
	arc( limit_z_ ); // core::Real
	arc( restrictSize_ ); // _Bool
	arc( ignoreBuriedPockets_ ); // _Bool
	arc( ignoreExposedPockets_ ); // _Bool
	arc( probe_rad_ ); // core::Real
	arc( surf_score_ ); // core::Real
	arc( surf_dist_ ); // core::Real
	arc( bur_score_ ); // core::Real
	arc( bur_dist_ ); // core::Real
	arc( side_chains_only_ ); // _Bool
	arc( markpsp_ ); // _Bool
	arc( marksps_ ); // _Bool
	arc( exemplarRestriction_ ); // _Bool
    arc( dumpExemplars_ ); // _Bool
    arc( filterExemplars_ ); // _Bool
    arc( limitExemplars_ ); // _Bool
	arc( search13_ ); // _Bool
	arc( minPockSize_ ); // core::Real
	arc( maxPockSize_ ); // core::Real
	arc( rot_mat_ ); // numeric::xyzMatrix<core::Real>
	arc( static_grid_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::PocketGrid );
CEREAL_REGISTER_TYPE( protocols::pockets::PocketGrid )


/// @brief Default constructor required by cereal to deserialize this class
protocols::pockets::CCluster::CCluster() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::CCluster::save( Archive & arc ) const {
	arc( CEREAL_NVP( points_ ) ); // std::list<Cxyz>
	arc( CEREAL_NVP( count_ ) ); // int
	arc( CEREAL_NVP( target_ ) ); // _Bool
	arc( CEREAL_NVP( subtarget_ ) ); // _Bool
	arc( CEREAL_NVP( solventExposed_ ) ); // _Bool
	arc( CEREAL_NVP( maxX ) ); // core::Size
	arc( CEREAL_NVP( minX ) ); // core::Size
	arc( CEREAL_NVP( maxY ) ); // core::Size
	arc( CEREAL_NVP( minY ) ); // core::Size
	arc( CEREAL_NVP( maxZ ) ); // core::Size
	arc( CEREAL_NVP( minZ ) ); // core::Size
	arc( CEREAL_NVP( step ) ); // core::Real
	arc( CEREAL_NVP( filterExemplars_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::CCluster::load( Archive & arc ) {
	arc( points_ ); // std::list<Cxyz>
	arc( count_ ); // int
	arc( target_ ); // _Bool
	arc( subtarget_ ); // _Bool
	arc( solventExposed_ ); // _Bool
	arc( maxX ); // core::Size
	arc( minX ); // core::Size
	arc( maxY ); // core::Size
	arc( minY ); // core::Size
	arc( maxZ ); // core::Size
	arc( minZ ); // core::Size
	arc( step ); // core::Real
	arc( filterExemplars_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::CCluster );

/// @brief Default constructor required by cereal to deserialize this class
protocols::pockets::PCluster::PCluster() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::PCluster::save( Archive & arc ) const {
	arc( CEREAL_NVP( points_ ) ); // std::list<Cxyz>
	arc( CEREAL_NVP( count_ ) ); // int
	arc( CEREAL_NVP( target_ ) ); // _Bool
	arc( CEREAL_NVP( subtarget_ ) ); // _Bool
	arc( CEREAL_NVP( solventExposed_ ) ); // _Bool
	arc( CEREAL_NVP( maxX ) ); // core::Size
	arc( CEREAL_NVP( minX ) ); // core::Size
	arc( CEREAL_NVP( maxY ) ); // core::Size
	arc( CEREAL_NVP( minY ) ); // core::Size
	arc( CEREAL_NVP( maxZ ) ); // core::Size
	arc( CEREAL_NVP( minZ ) ); // core::Size
	arc( CEREAL_NVP( step ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::PCluster::load( Archive & arc ) {
	arc( points_ ); // std::list<Cxyz>
	arc( count_ ); // int
	arc( target_ ); // _Bool
	arc( subtarget_ ); // _Bool
	arc( solventExposed_ ); // _Bool
	arc( maxX ); // core::Size
	arc( minX ); // core::Size
	arc( maxY ); // core::Size
	arc( minY ); // core::Size
	arc( maxZ ); // core::Size
	arc( minZ ); // core::Size
	arc( step ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::PCluster );

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::ElectrostaticpotentialGrid::save( Archive & arc ) const {
	arc( cereal::base_class< PocketGrid >( this ) );
	arc( CEREAL_NVP( espGrid_dim_ ) ); // numeric::xyzVector<core::Size>
	arc( CEREAL_NVP( espGrid_mid_ ) ); // numeric::xyzVector<core::Real>
	arc( CEREAL_NVP( espGrid_spacing_ ) ); // core::Real
	arc( CEREAL_NVP( espGrid_points_list_ ) ); // std::list<utility::vector1<core::Real> >
	arc( CEREAL_NVP( espGrid_ ) ); // std::vector<std::vector<std::vector<core::Real> > >
	arc( CEREAL_NVP( typGrid_ ) ); // std::vector<std::vector<std::vector<PtType> > >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::ElectrostaticpotentialGrid::load( Archive & arc ) {
	arc( cereal::base_class< PocketGrid >( this ) );
	arc( espGrid_dim_ ); // numeric::xyzVector<core::Size>
	arc( espGrid_mid_ ); // numeric::xyzVector<core::Real>
	arc( espGrid_spacing_ ); // core::Real
	arc( espGrid_points_list_ ); // std::list<utility::vector1<core::Real> >
	arc( espGrid_ ); // std::vector<std::vector<std::vector<core::Real> > >
	arc( typGrid_ ); // std::vector<std::vector<std::vector<PtType> > >
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::ElectrostaticpotentialGrid );
CEREAL_REGISTER_TYPE( protocols::pockets::ElectrostaticpotentialGrid )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::PClusterSet::save( Archive & arc ) const {
	arc( CEREAL_NVP( clusters_ ) ); // std::list<PCluster>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::PClusterSet::load( Archive & arc ) {
	arc( clusters_ ); // std::list<PCluster>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::PClusterSet );

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::CClusterSet::save( Archive & arc ) const {
	arc( CEREAL_NVP( clusters_ ) ); // std::list<CCluster>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::CClusterSet::load( Archive & arc ) {
	arc( clusters_ ); // std::list<CCluster>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::CClusterSet );

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::CCluster::Cxyz::save( Archive & arc ) const {
	arc( CEREAL_NVP( x ) ); // core::Size
	arc( CEREAL_NVP( y ) ); // core::Size
	arc( CEREAL_NVP( z ) ); // core::Size
	arc( CEREAL_NVP( absX ) ); // core::Real
	arc( CEREAL_NVP( absY ) ); // core::Real
	arc( CEREAL_NVP( absZ ) ); // core::Real
	arc( CEREAL_NVP( atom_type ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::CCluster::Cxyz::load( Archive & arc ) {
	arc( x ); // core::Size
	arc( y ); // core::Size
	arc( z ); // core::Size
	arc( absX ); // core::Real
	arc( absY ); // core::Real
	arc( absZ ); // core::Real
	arc( atom_type ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::CCluster::Cxyz );

/// @brief Default constructor required by cereal to deserialize this class
protocols::pockets::TargetPocketGrid::TargetPocketGrid() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::TargetPocketGrid::save( Archive & arc ) const {
	arc( cereal::base_class< PocketGrid >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::TargetPocketGrid::load( Archive & arc ) {
	arc( cereal::base_class< PocketGrid >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::TargetPocketGrid );
CEREAL_REGISTER_TYPE( protocols::pockets::TargetPocketGrid )


/// @brief Default constructor required by cereal to deserialize this class
protocols::pockets::ComparisonGrid::ComparisonGrid() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::ComparisonGrid::save( Archive & arc ) const {
	arc( CEREAL_NVP( grid_ ) ); // std::vector<std::vector<std::vector<PtType> > >
	arc( CEREAL_NVP( xdim_ ) ); // core::Size
	arc( CEREAL_NVP( ydim_ ) ); // core::Size
	arc( CEREAL_NVP( zdim_ ) ); // core::Size
	arc( CEREAL_NVP( xcorn_ ) ); // core::Real
	arc( CEREAL_NVP( ycorn_ ) ); // core::Real
	arc( CEREAL_NVP( zcorn_ ) ); // core::Real
	arc( CEREAL_NVP( stepSize_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::ComparisonGrid::load( Archive & arc ) {
	arc( grid_ ); // std::vector<std::vector<std::vector<PtType> > >
	arc( xdim_ ); // core::Size
	arc( ydim_ ); // core::Size
	arc( zdim_ ); // core::Size
	arc( xcorn_ ); // core::Real
	arc( ycorn_ ); // core::Real
	arc( zcorn_ ); // core::Real
	arc( stepSize_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::ComparisonGrid );
CEREAL_REGISTER_DYNAMIC_INIT( protocols_pockets_PocketGrid )
#endif // SERIALIZATION
