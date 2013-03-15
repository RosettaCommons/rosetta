// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/PocketGrid.cc
/// @brief  protocols::pockets::PocketGrid functions
/// @author David Johnson
/// @author Ragul Gowthaman

#include <protocols/pockets/PocketGrid.hh>

// Core Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyz.functions.fwd.hh>

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

namespace protocols {
namespace pockets {

static basic::Tracer TR("core.grid.Pockets.PocketGrid");

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
  target=false;
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
  target = old.target;
  maxX = old.maxX; maxY = old.maxY; maxZ = old.maxZ;
  minX = old.minX; minY = old.minY; minZ = old.minZ;
  step = old.step;
}


bool PCluster::isClose(PCluster const & c2) const{
  if (minX>c2.maxX+1) return false;
  if ((int)maxX<(int)c2.minX-1) return false;
  if (minY>c2.maxY+1) return false;
  if ((int)maxY<(int)c2.minY-1) return false;
  if (minZ>c2.maxZ+1) return false;
  if ((int)maxZ<(int)c2.minZ-1) return false;
  return true;
}

bool PCluster::testNeighbor(PCluster & c2){
  if (c2.points_.size()>100 &&points_.size()>100){
    std::list<Cxyz>::iterator i = points_.begin();
    int icount=0;
    while (i!=points_.end()){
      std::list<Cxyz>::iterator j = c2.points_.begin();
      int ictend=std::min(500,(int)points_.size()-icount);
      for (int ict=0; ict<ictend; ++i, ++ict){
        int jcount=0;
        while (j!=c2.points_.end()){
          int jctend=std::min(500,(int)c2.points_.size()-jcount);
          for (int jct=0; jct<jctend; ++j, ++jct){
            if (abs((int) i->x - (int) j->x)<=1)
              if (abs((int) i->y - (int) j->y)<=1)
                if (abs((int) i->z - (int) j->z)<=1){
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
          jcount+=jctend;
        }
      }
      icount+=ictend;
    }

  }else{

    for (std::list<Cxyz>::iterator i = points_.begin(); i!=points_.end(); ++i){
      for (std::list<Cxyz>::iterator j = c2.points_.begin(); j!=c2.points_.end(); ++j){
        if (abs((int) i->x - (int) j->x)<=1)
          if (abs((int) i->y - (int) j->y)<=1)
            if (abs((int) i->z - (int) j->z)<=1){

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
  return false;
}

PClusterSet::PClusterSet(){
  clusters_.clear();
}

PClusterSet& PClusterSet::operator= (const PClusterSet& old){
  if (this != &old){
    clusters_ = old.clusters_;
  }
  return *this;
}

void PClusterSet::clear(){
  clusters_.clear();
}

std::_List_iterator<protocols::pockets::PCluster> PClusterSet::add(core::Size x, core::Size y, core::Size z, core::Real step){
  PCluster tmp(x,y,z, step);
  clusters_.push_front(tmp);
  return clusters_.begin();
}

void PClusterSet::join(std::list<PCluster>::iterator c1, std::list<PCluster>::iterator c2){
  c1->points_.splice(c1->points_.end(), c2->points_);
  clusters_.erase(c2);
}

void PClusterSet::findClusters(){
  for (std::list<PCluster>::iterator i = clusters_.begin(); i!=clusters_.end(); ++i){
    for (std::list<PCluster>::iterator j = i; j!=clusters_.end(); ++j){
      if (i==j) continue;
      if (j->isClose(*i)) {
        if (j->testNeighbor(*i)) {
          std::list<PCluster>::iterator oldI = i;
          i--;
          clusters_.erase(oldI);
          break;
        }
      }
    }
  }
}

core::Real PClusterSet::getLargestClusterSize( core::Real const & stepSize, core::Real const & minClusterSize ){
  core::Real largest_size=0.;
  for (std::list<PCluster>::iterator i = clusters_.begin(); i!=clusters_.end(); ++i){
    if (!i->isTarget()) continue;
    if (i->size()>largest_size){
      largest_size=i->size();
    }
  }
  largest_size = largest_size*pow(stepSize, 3)-minClusterSize;
  if (largest_size > 0) return largest_size;

  return 0;
}

core::Real PClusterSet::getNetClusterSize( core::Real const & stepSize, core::Real const & minClusterSize ){
  core::Real total_size=0.;
  for (std::list<PCluster>::iterator i = clusters_.begin(); i!=clusters_.end(); ++i){
    if (!i->isTarget()) continue;
    if (i->size()*pow(stepSize, 3) >minClusterSize){
      total_size+=i->size()*pow(stepSize, 3)-minClusterSize;
    }
  }
  if (total_size>0) return total_size;

  return 0;
}

void PocketGrid::clear(){
  grid_.clear();
  pockets_.clear();
  clusters_.clear();
}

void PocketGrid::init(){
  clear();
  grid_.resize(xdim_);
  pockets_.resize(xdim_);
  for (core::Size tx=0;tx<xdim_;tx++){
    grid_[tx].resize(ydim_);
    pockets_[tx].resize(ydim_);
    for (core::Size ty=0;ty<ydim_;ty++){
      grid_[tx][ty].resize(zdim_, EMPTY);
      pockets_[tx][ty].resize(zdim_, EMPTY);
    }
  }
}


void PocketGrid::setup_default_options(){

  using namespace basic::options;
  spacing_=option[ OptionKeys::pocket_grid::pocket_grid_spacing ]();
  if (option[ OptionKeys::pocket_grid::pocket_grid_size ]()>0){
    size_x_ = option[ OptionKeys::pocket_grid::pocket_grid_size ]();
    size_y_ = option[ OptionKeys::pocket_grid::pocket_grid_size ]();
    size_z_ = option[ OptionKeys::pocket_grid::pocket_grid_size ]();
  }else{
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
  markpsp_=option[ OptionKeys::pocket_grid::pocket_psp ]();
  marksps_=option[ OptionKeys::pocket_grid::pocket_sps ]();
  surf_score_=option[ OptionKeys::pocket_grid::pocket_surface_score ]();
  surf_dist_=option[ OptionKeys::pocket_grid::pocket_surface_dist ]();
  bur_score_=option[ OptionKeys::pocket_grid::pocket_buried_score ]();
  bur_dist_=option[ OptionKeys::pocket_grid::pocket_buried_dist ]();
  minPockSize_=option[ OptionKeys::pocket_grid::pocket_min_size ]();
  maxPockSize_=option[ OptionKeys::pocket_grid::pocket_max_size ]();
  rot_mat_ = numeric::x_rotation_matrix_degrees((core::Real)0);

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
  tag_=gr.tag_;
  pdbno_=gr.pdbno_;
  minPockSize_ = gr.minPockSize_;

  grid_.resize(xdim_);
  pockets_.resize(xdim_);
  for (core::Size tx=0;tx<xdim_;tx++){
    grid_[tx].resize(ydim_);
    pockets_[tx].resize(ydim_);
    for (core::Size ty=0;ty<ydim_;ty++){
      grid_[tx][ty].resize(zdim_);
      pockets_[tx][ty].resize(zdim_);
      for (core::Size tz=0;tz<zdim_;tz++){
        grid_[tx][ty][tz] = gr.grid_[tx][ty][tz];
        pockets_[tx][ty][tz] = gr.pockets_[tx][ty][tz];
      }
    }
  }

}

PocketGrid::~PocketGrid() {}

PocketGrid& PocketGrid::operator=(const PocketGrid& gr){
  if (this != &gr){

    xdim_ = gr.xdim_; ydim_ = gr.ydim_; zdim_ = gr.zdim_;
    xcorn_ = gr.xcorn_; ycorn_ = gr.ycorn_; zcorn_ = gr.zcorn_;
    stepSize_ = gr.stepSize_;
    markpsp_ = gr.markpsp_;
    marksps_ = gr.marksps_;
    tag_=gr.tag_;
    pdbno_=gr.pdbno_;
    minPockSize_ = gr.minPockSize_;

    grid_.resize(xdim_);
    pockets_.resize(xdim_);
    for (core::Size tx=0;tx<xdim_;tx++){
      grid_[tx].resize(ydim_);
      pockets_[tx].resize(ydim_);
      for (core::Size ty=0;ty<ydim_;ty++){
        grid_[tx][ty].resize(zdim_);
        pockets_[tx][ty].resize(zdim_);
        for (core::Size tz=0;tz<zdim_;tz++){
          grid_[tx][ty][tz] = gr.grid_[tx][ty][tz];
          pockets_[tx][ty][tz] = gr.pockets_[tx][ty][tz];
        }
      }
    }
  }
  return *this;
}


PocketGrid::PocketGrid(){
  setup_default_options();
}

PocketGrid::PocketGrid( core::conformation::Residue const & central_rsd ) {
  setup_default_options();
  initialize(central_rsd, size_x_, size_y_, size_z_, spacing_, markpsp_, marksps_);
  using namespace basic::options;
  tag_=option[ OptionKeys::out::output_tag ]();
}

PocketGrid::PocketGrid( std::vector< core::conformation::ResidueOP > const & central_rsds ) {
  setup_default_options();
  initialize(central_rsds, size_x_, size_y_, size_z_, spacing_, markpsp_, marksps_);
  using namespace basic::options;
  tag_=option[ OptionKeys::out::output_tag ]();
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

void PocketGrid::initialize( core::Real const & xc, core::Real const & yc, core::Real const & zc, core::Real x, core::Real y, core::Real z, core::Real const & stepSize, bool psp, bool sps){
  markpsp_=psp;
  marksps_=sps;
  stepSize_=stepSize;
  xdim_=2*(core::Size)ceil(x/stepSize_)+1;
  ydim_=2*(core::Size)ceil(y/stepSize_)+1;
  zdim_=2*(core::Size)ceil(z/stepSize_)+1;
  init();
  recenter(xc,yc,zc);
  pdbno_=0;
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
  init();
  recenter(xc,yc,zc);
  pdbno_=0;
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
  init();
  recenter(xc,yc,zc);
  pdbno_=0;
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
  init();
  recenter(center);
  pdbno_=0;
}

void PocketGrid::initialize (core::conformation::Residue const & central_rsd, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize, bool psp, bool sps){
  markpsp_=psp;
  marksps_=sps;
  stepSize_=stepSize;
  xdim_=2*(core::Size)ceil(x/stepSize_)+1;
  ydim_=2*(core::Size)ceil(y/stepSize_)+1;
  zdim_=2*(core::Size)ceil(z/stepSize_)+1;
  init();
  recenter(central_rsd);
  pdbno_=0;
}

void PocketGrid::initialize (std::vector< core::conformation::ResidueOP > const & central_rsds, core::Real const & x, core::Real const & y, core::Real const & z, core::Real const & stepSize, bool psp, bool sps){
  markpsp_=psp;
  marksps_=sps;
  stepSize_=stepSize;
  xdim_=2*(core::Size)ceil(x/stepSize_)+1;
  ydim_=2*(core::Size)ceil(y/stepSize_)+1;
  zdim_=2*(core::Size)ceil(z/stepSize_)+1;
  init();
  recenter(central_rsds);
  pdbno_=0;
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
  numeric::xyzVector<core::Real> rpoint = rotatePoint(xc, yc, zc);
  xcorn_=rpoint.x()-((xdim_-1)/2*stepSize_);
  ycorn_=rpoint.y()-((ydim_-1)/2*stepSize_);
  zcorn_=rpoint.z()-((zdim_-1)/2*stepSize_);
}

void PocketGrid::recenter( core::conformation::Residue const & central_rsd ){

  assert( central_rsd.is_protein() );
  core::Vector center;
  center(1)=0.0;
  center(2)=0.0;
  center(3)=0.0;
  core::Size count=0;
  for(Size i = 1, i_end = central_rsd.nheavyatoms(); i <= i_end; ++i) {
    if (central_rsd.atom(i).type()<18||central_rsd.atom(i).type()>21){
      center(1)+=central_rsd.atom(i).xyz()(1);
      center(2)+=central_rsd.atom(i).xyz()(2);
      center(3)+=central_rsd.atom(i).xyz()(3);
      count++;
    }
  }
  if (count) {
    center(1)/=count;
    center(2)/=count;
    center(3)/=count;
  }
  recenter(center(1), center(2), center(3));

}

void PocketGrid::recenter( std::vector< core::conformation::ResidueOP > const & central_rsds ){
  core::Vector center;
  center(1)=0.0;
  center(2)=0.0;
  center(3)=0.0;
  core::Size count=0;
  int sz = central_rsds.size();
  assert (sz>0);
  for (int rnum = 0; rnum < sz; ++rnum){
    core::conformation::ResidueOP central_rsd = central_rsds[rnum];

    assert( central_rsd->is_protein() );
    for(Size i = 1, i_end = central_rsd->nheavyatoms(); i <= i_end; ++i) {
      if (central_rsd->atom(i).type()<18||central_rsd->atom(i).type()>21){
        center(1)+=central_rsd->atom(i).xyz()(1);
        center(2)+=central_rsd->atom(i).xyz()(2);
        center(3)+=central_rsd->atom(i).xyz()(3);
        count++;
      }
    }
  }
  if (count) {
    center(1)/=count;
    center(2)/=count;
    center(3)/=count;
  }
  recenter(center(1), center(2), center(3));

}
void PocketGrid::recenter(core::Vector const & center){
  recenter(center(1), center(2), center(3));
}

std::vector< core::conformation::ResidueOP > PocketGrid::getRelaxResidues( core::pose::Pose const & input_pose, std::string const & resids ) {
  std::vector< core::conformation::ResidueOP > residues;

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
      if (fpos != resid.size()-1 ) {
        chain = resid[ fpos+1 ];
      }
    } else {
      central_relax_pdb_number = ObjexxFCL::int_of( resid );
    }

    for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j ) {
      if ( input_pose.pdb_info()->number(j) == central_relax_pdb_number ) {
        if (chain != ' '){
          if ( input_pose.pdb_info()->chain(j) == chain ) {
            residues.push_back(input_pose.conformation().residue(j).clone());
            continue;
          }
        }else{
          residues.push_back(input_pose.conformation().residue(j).clone());
          continue;
        }
      }
    }
  }
  return residues;
}

void PocketGrid::findPockets(core::Size thr, core::Real max){
  thr=(core::Size)std::floor((core::Real)(thr)/stepSize_+.5);
  core::Size thr2=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(3.)+.5);
  core::Size max1=(core::Size)std::floor(max/stepSize_+.5);
  core::Size max2=(core::Size)std::floor((max/stepSize_)/sqrt(3.)+.5);
  newSearch(thr, thr2, max1, max2);
}

void PocketGrid::findPSP(core::Size thr, core::Real max){
  thr=(core::Size)std::floor((core::Real)(thr)/stepSize_+.5);
  core::Size thr2=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(2.)+.5);
  core::Size max1=(core::Size)std::floor(max/stepSize_+.5);
  core::Size max2=(core::Size)std::floor((max/stepSize_)/sqrt(2.)+.5);
  newSearch(thr, thr2, max1, max2, true);
}

void PocketGrid::findSPS(core::Size thr, core::Real max){
  thr=(core::Size)std::floor((core::Real)(thr)/stepSize_+.5);
	core::Size thr2=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(2.)+.5);
  core::Size max1=(core::Size)std::floor(max/stepSize_+.5);
  core::Size max2=(core::Size)std::floor((max/stepSize_)/sqrt(2.)+.5);
  newSearch(thr, thr2, max1, max2, false, true);
}

bool PocketGrid::fill(core::Size x, core::Size y,core::Size z){
  if (grid_[x][y][z]!=TP_POCKET && grid_[x][y][z]!=POCKET  ){
    return false;
  }
  if (pockets_[x][y][z]!=EMPTY) return false;

  grid_[x][y][z]=TP_POCKET;
  pockets_[x][y][z]=TP_POCKET;
  bool touching=false;
  core::Size x1,y1,z1;
  if (x==0) x1=0;
  else x1=x-1;
  if (y==0) y1=0;
  else y1=y-1;
  if (z==0) z1=0;
  else z1=z-1;
  for(core::Size x2=x1; x2<std::min(xdim_, x+2); ++x2){
    for(core::Size y2=y1; y2<std::min(ydim_, y+2); ++y2){
      for(core::Size z2=z1; z2<std::min(zdim_, z+2); ++z2){
        touching=touching || fill(x2,y2,z2);
      }}}
  return touching;
}

bool PocketGrid::touchesSolvent(core::Size x, core::Size y,core::Size z) const{
  bool touching=false;
  if ((x==0) || (y==0)||(z==0)) return true;
  if ((x==xdim_-1) || (y==ydim_-1)||(z==zdim_-1)) return true;
  if (grid_[x-1][y][z]==EMPTY)return true;
  if (grid_[x+1][y][z]==EMPTY)return true;
  if (grid_[x][y-1][z]==EMPTY)return true;
  if (grid_[x][y+1][z]==EMPTY)return true;
  if (grid_[x][y][z-1]==EMPTY)return true;
  if (grid_[x][y][z+1]==EMPTY)return true;

	return touching;
}

bool PocketGrid::touchesSS(core::Size x, core::Size y,core::Size z) const {
  if (x!=0) if (grid_[x-1][y][z]==EMPTY) return true;
  if (y!=0) if (grid_[x][y-1][z]==EMPTY) return true;
  if (z!=0) if (grid_[x][y][z-1]==EMPTY) return true;
  if ((x!=xdim_-1)) if (grid_[x+1][y][z]==EMPTY)return true;
  if ((y!=ydim_-1)) if (grid_[x][y+1][z]==EMPTY)return true;
  if ((z!=zdim_-1)) if (grid_[x][y][z+1]==EMPTY)return true;
  return false;
}

bool PocketGrid::touchesPS(core::Size x, core::Size y,core::Size z) const {
  if (x!=0) {
    if (grid_[x-1][y][z]==HSURFACE) return true;
    if (grid_[x-1][y][z]==PSURFACE) return true;
    if (grid_[x-1][y][z]==T_SURFACE) return true;
  }
  if (y!=0) {
    if (grid_[x][y-1][z]==HSURFACE) return true;
    if (grid_[x][y-1][z]==PSURFACE) return true;
    if (grid_[x][y-1][z]==T_SURFACE) return true;
  }
  if (z!=0) {
    if (grid_[x][y][z-1]==HSURFACE) return true;
    if (grid_[x][y][z-1]==PSURFACE) return true;
    if (grid_[x][y][z-1]==T_SURFACE) return true;
  }
  if ((x!=xdim_-1)) {
    if (grid_[x+1][y][z]==HSURFACE) return true;
    if (grid_[x+1][y][z]==PSURFACE) return true;
    if (grid_[x+1][y][z]==T_SURFACE) return true;
  }
  if ((y!=ydim_-1)) {
    if (grid_[x][y+1][z]==HSURFACE) return true;
    if (grid_[x][y+1][z]==PSURFACE) return true;
    if (grid_[x][y+1][z]==T_SURFACE) return true;
  }
  if ((z!=zdim_-1)) {
    if (grid_[x][y][z+1]==HSURFACE) return true;
    if (grid_[x][y][z+1]==PSURFACE) return true;
    if (grid_[x][y][z+1]==T_SURFACE) return true;
  }
  return false;
}


bool PocketGrid::touchesSurface(core::Size x, core::Size y,core::Size z, bool polar, bool either) const {
  bool touching=false;
  if (polar){
    if (x!=0) if (grid_[x-1][y][z]==PSURFACE)return true;
    if (x!=xdim_-1) if (grid_[x+1][y][z]==PSURFACE)return true;
    if (y!=0) if (grid_[x][y-1][z]==PSURFACE)return true;
    if (y!=ydim_-1) if (grid_[x][y+1][z]==PSURFACE)return true;
    if (z!=0) if (grid_[x][y][z-1]==PSURFACE)return true;
    if (z!=zdim_-1) if (grid_[x][y][z+1]==PSURFACE)return true;
  }else{
    if (x!=0) if (grid_[x-1][y][z]==HSURFACE)return true;
    if (x!=xdim_-1) if (grid_[x+1][y][z]==HSURFACE)return true;
    if (y!=0) if (grid_[x][y-1][z]==HSURFACE)return true;
    if (y!=ydim_-1) if (grid_[x][y+1][z]==HSURFACE)return true;
    if (z!=0) if (grid_[x][y][z-1]==HSURFACE)return true;
    if (z!=zdim_-1) if (grid_[x][y][z+1]==HSURFACE)return true;
  }
  if (either){
    if (x!=0) if (grid_[x-1][y][z]==PSURFACE)return true;
    if (x!=xdim_-1) if (grid_[x+1][y][z]==PSURFACE)return true;
    if (y!=0) if (grid_[x][y-1][z]==PSURFACE)return true;
    if (y!=ydim_-1) if (grid_[x][y+1][z]==PSURFACE)return true;
    if (z!=0) if (grid_[x][y][z-1]==PSURFACE)return true;
    if (z!=zdim_-1) if (grid_[x][y][z+1]==PSURFACE)return true;
    if (x!=0) if (grid_[x-1][y][z]==HSURFACE)return true;
    if (x!=xdim_-1) if (grid_[x+1][y][z]==HSURFACE)return true;
    if (y!=0) if (grid_[x][y-1][z]==HSURFACE)return true;
    if (y!=ydim_-1) if (grid_[x][y+1][z]==HSURFACE)return true;
    if (z!=0) if (grid_[x][y][z-1]==HSURFACE)return true;
    if (z!=zdim_-1) if (grid_[x][y][z+1]==HSURFACE)return true;
  }
  return touching;
}

void PocketGrid::fillTargetPockets(){
  for(core::Size x=0; x<xdim_; ++x){
    for(core::Size y=0; y<ydim_; ++y){
      for(core::Size z=0; z<zdim_; ++z){
        pockets_[x][y][z]=EMPTY;
      }}}
  for(core::Size x=0; x<xdim_; ++x){
    for(core::Size y=0; y<ydim_; ++y){
      for(core::Size z=0; z<zdim_; ++z){
        if (grid_[x][y][z]==TP_POCKET){
          pockets_[x][y][z]=TP_POCKET;
          core::Size x1,y1,z1;
          if (x==0) x1=0;
          else x1=x-1;
          if (y==0) y1=0;
          else y1=y-1;
          if (z==0) z1=0;
          else z1=z-1;
          for(core::Size x2=x1; x2<std::min(xdim_, x+2); ++x2){
            for(core::Size y2=y1; y2<std::min(ydim_, y+2); ++y2){
              for(core::Size z2=z1; z2<std::min(zdim_, z+2); ++z2){
                fill(x2,y2,z2);
              }}}
        }
      }}}
}


void PocketGrid::newSearch(core::Size thr1, core::Size thr2, core::Size max1, core::Size max2, bool psp, bool sps){
  //vector of deltas to add to x,y,z.  Index % 3 = 0 for x deltas, 1 for y, 2 for z
  std::vector<int> deltas;
  int dirs=7;
  deltas.push_back(0);deltas.push_back(0);deltas.push_back(1);
  deltas.push_back(0);deltas.push_back(1);deltas.push_back(0);
  deltas.push_back(1);deltas.push_back(0);deltas.push_back(0);
  deltas.push_back(1);deltas.push_back(1);deltas.push_back(1);
  deltas.push_back(-1);deltas.push_back(1);deltas.push_back(1);
  deltas.push_back(-1);deltas.push_back(-1);deltas.push_back(1);
  deltas.push_back(1);deltas.push_back(-1);deltas.push_back(1);

  core::Size thr;
  core::Size max;
  //iterate over the entire grid
  for (core::Size cx=0; cx<xdim_;cx++){
    for (core::Size cy=0; cy<ydim_;cy++){
      for (core::Size cz=0; cz<zdim_;cz++){
        bool tar_surf=false;

        //if a point is a surface, or pocket and doind psp searching, or solvent and doing a sps search, search in all directions from that point.
        if ((!sps&&(grid_[cx][cy][cz]==HSURFACE || grid_[cx][cy][cz]==PSURFACE || grid_[cx][cy][cz]==T_SURFACE||(grid_[cx][cy][cz]==POCKET&&psp)||(grid_[cx][cy][cz]==TP_POCKET&&psp)||(grid_[cx][cy][cz]==POCKET&&sps)))||(grid_[cx][cy][cz]==EMPTY&&sps)){
          if (grid_[cx][cy][cz]==T_SURFACE||(grid_[cx][cy][cz]==TP_POCKET&&psp)) tar_surf=true;
          //iterate over the searching directions, get the index associated with x
          for (int i=0;i<dirs*3;i+=3){
            if (i<=6) {
              thr=thr1;
              max=max1;
            }
            else {
              thr=thr2;
              max=max2;
            }
            bool t_surf=tar_surf;
            int x;
            int y;
            int z;
            int count;
            bool marked=true;
            //starting at that point, search outward in the direction being searched until max distance
            for (x=(int)cx+deltas[i],y=(int)cy+deltas[i+1],z=(int)cz+deltas[i+2], count=0; count<(int)max; x+=deltas[i],y+=deltas[i+1],z+=deltas[i+2], count++){
              if (x<0 || x>=(int)xdim_) break;
              if (y<0 || y>=(int)ydim_) break;
              if (z<0 || z>=(int)zdim_) break;

              //Just keep going for EMPTY (solvent)
              if (grid_[x][y][z]==EMPTY&&!sps) {
                marked=false;
                continue;
              }

              //if not doing psp, treat pockets as EMPTY coz they are still solvent
              else if((grid_[x][y][z]==POCKET&&!psp) || (grid_[x][y][z]==TP_POCKET&&!psp)) {
                if (sps) {
                  marked=false;
                }
                continue;
              }

              //if it's a surface, we want to go back to the starting point and fill the solvents in as pockets
              else if (sps&&(grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE)) break;
              else if ((!sps&&(grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE || (grid_[x][y][z]==POCKET&&psp) || (grid_[x][y][z]==TP_POCKET&&psp))) || (grid_[x][y][z]==EMPTY&&sps)){

                //if there is no EMPTY, then everything has been marked already.  No need to go further.
                if (marked) break;


                if (grid_[x][y][z]==T_SURFACE||(grid_[x][y][z]==TP_POCKET&&psp)) t_surf=true;
                if ((count>(int)thr && !psp)||psp||sps){
                  for (int c=count;c>0;c--){
                    if (grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]==HSURFACE || grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]] == PSURFACE || grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]==T_SURFACE) std::cout<<"MAJOR ERROR, overwriting surface with pocket\n";
                    if (sps){
                      grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]=EMPTY;
                    }
                    else if (markpsp_||!psp){
                      if (t_surf) grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]=TP_POCKET;
                      else if (grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]!=TP_POCKET)  grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]=POCKET;
                    }else{
                      if (grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]==EMPTY) grid_[x-c*deltas[i]][y-c*deltas[i+1]][z-c*deltas[i+2]]=PSP;
                    }
                  }
                }
                break;
              }
              //if it's not a surface, pocket, or solvent then can't be a pocket.  Stop searching this direction
              else break;
            }

          }

        }
      }
    }}

}


void PocketGrid::mark(core::Vector const & center, core::Real const & vdWd, core::Real const & buffer, bool polar, bool targetResi){
  mark(center(1),center(2),center(3), vdWd, buffer, polar, targetResi);
}

void PocketGrid::mark(core::Real x, core::Real y, core::Real z, core::Real const & vdWd, core::Real const & buffer, bool polar, bool targetResi){
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

  for (int xIter=minX; xIter<=maxX; ++xIter){
    for (int yIter=minY; yIter<=maxY; ++yIter){
      int centerZ=std::max(minZ, zcen);
      for (int zIter=centerZ;zIter <= maxZ; ++zIter){
        if (pow((xIter-xcen),2)+pow((yIter-ycen),2)+pow((zIter-zcen),2)>pow(radius, 2)) continue;
        if (pow((xIter-xcen),2)+pow((yIter-ycen),2)+pow((zIter-zcen),2)>pow(vdW, 2)){
          if (grid_[xIter][yIter][zIter]!=PROTEIN &&grid_[xIter][yIter][zIter]!=TARGET) {
            if (targetResi){
              grid_[xIter][yIter][zIter]=T_SURFACE;
            }else if (polar) { if (grid_[xIter][yIter][zIter]!=HSURFACE) grid_[xIter][yIter][zIter]=PSURFACE;
            }else{grid_[xIter][yIter][zIter]=HSURFACE;
            }
          }else if (targetResi) grid_[xIter][yIter][zIter]=TARGET;
        }else{
          if (targetResi) {
            grid_[xIter][yIter][zIter]=TARGET;
          }else{
            if (grid_[xIter][yIter][zIter]!=TARGET) grid_[xIter][yIter][zIter]=PROTEIN;
          }
        }
      }

      centerZ=std::min(maxZ, zcen-1);
      for (int zIter=centerZ;zIter >= minZ; --zIter){
        if (pow(xIter-xcen,2)+pow(yIter-ycen,2)+pow(zIter-zcen,2)>pow(radius, 2)) continue;
        if (pow(xIter-xcen,2)+pow(yIter-ycen,2)+pow(zIter-zcen,2)>pow(vdW, 2)){
          if (grid_[xIter][yIter][zIter]!=PROTEIN &&grid_[xIter][yIter][zIter]!=TARGET) {
            if (targetResi){
              grid_[xIter][yIter][zIter]=T_SURFACE;
            }else if (polar) { if (grid_[xIter][yIter][zIter]!=HSURFACE) grid_[xIter][yIter][zIter]=PSURFACE;
            }else{grid_[xIter][yIter][zIter]=HSURFACE;
            }
          }else if (targetResi) grid_[xIter][yIter][zIter]=TARGET;

        }else{
          if (targetResi) {
            grid_[xIter][yIter][zIter]=TARGET;
          }else{
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

void PocketGrid::print(){
  for (core::Size z=0; z<zdim_;z++){
    for (core::Size x=0;x<xdim_;x++){
      for (core::Size y=0;y<ydim_;y++){
        std::cout << pockets_[x][y][z];
      }
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
  }
}

void PocketGrid::dumpGridToFile(){
  std::stringstream filename;
  filename<<tag_<<"pocket"<<pdbno_<<".pdb";
  pdbno_++;
	dumpGridToFile( filename.str() );
}


void PocketGrid::dumpGridToFile( std::string const & output_filename ){

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_filename, std::ios::out);
  int counter=1;
  int counter2=1;

  core::Size x,y,z;
  for (x=0;x<(xdim_); x++){
    for (y=0;y<(ydim_); y++){
      for (z=0;z<(zdim_); z++){
        std::string concatenated_pdb_info;
        concatenated_pdb_info += "ATOM  ";
        std::stringstream  tmp;
        tmp<<counter2;
        if (counter2<10) concatenated_pdb_info += "    ";
        else if (counter2<100) concatenated_pdb_info += "   ";
        else if (counter2<1000) concatenated_pdb_info += "  ";
        else if (counter2<10000) concatenated_pdb_info += " ";
        concatenated_pdb_info += tmp.str()+"  ";
        if (grid_[x][y][z]==EMPTY) {
          continue;
        }
        if (grid_[x][y][z]==POCKET)  {
          continue;
        }
        if (grid_[x][y][z]==PO_SURF){
          continue;
        }
        if (grid_[x][y][z]==PO_BURIED)  {
          continue;
        }
        if (grid_[x][y][z]==HSURFACE) {
          continue;
        }
        if (grid_[x][y][z]==PSURFACE) {
          continue;
        }
        if (grid_[x][y][z]==TP_POCKET)  {
          continue;
        }
        if (grid_[x][y][z]==TP_SURF)  {
          continue;
        }
        if (grid_[x][y][z]==TP_BURIED)  {
          continue;
        }
        if (grid_[x][y][z]==PO_EDGE)  {
          continue;
        }
        if (grid_[x][y][z]==TP_EDGE)  {
          continue;
        }
        if (grid_[x][y][z]==T_SURFACE)  {
          concatenated_pdb_info += "TS  TS ";
        }
        if (grid_[x][y][z]==TARGET)  {
          continue;
        }
        if (grid_[x][y][z]==PROTEIN) {
          continue;
        }
        if (grid_[x][y][z]==PSP) {
          continue;
        }

        tmp.str(std::string());
        tmp<<"          "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x*stepSize_+xcorn_<<std::setw(8)<<y*stepSize_+ycorn_<<std::setw(8)<<z*stepSize_+zcorn_<<std::endl;
        concatenated_pdb_info += tmp.str();
        counter++;
        outPDB_stream<<concatenated_pdb_info;
      }}}

  int clustNo=1;
  bool smallPocket;
  for (std::list<PCluster>::iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit){
    if (cit->points_.size()*pow(stepSize_,3)<minPockSize_) smallPocket=true;
    else smallPocket=false;
    for (std::list<PCluster::Cxyz>::iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
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
      tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pit->x*stepSize_+xcorn_<<std::setw(8)<<pit->y*stepSize_+ycorn_<<std::setw(8)<<pit->z*stepSize_+zcorn_<<std::endl;
      concatenated_pdb_info += tmp.str();
      counter2++;
      outPDB_stream<<concatenated_pdb_info;

    }
    clustNo++;
  }
	outPDB_stream.close();
	outPDB_stream.clear();

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
  for (core::Size tx=0;tx<xdim_;tx++){
    cluster_grid[tx].resize(ydim_);
    for (core::Size ty=0;ty<ydim_;ty++){
      cluster_grid[tx][ty].resize(zdim_, clusters_.clusters_.end());
    }
  }

  // Search the grid for pockets
  core::Size x,y,z;
  for (x=0;x<(xdim_); x++){
    for (y=0;y<(ydim_); y++){
      for (z=0;z<(zdim_); z++){
        // If you find a pocket, check if its neighbors are already part of a cluster
        if (grid_[x][y][z]==TP_POCKET || grid_[x][y][z]==POCKET || grid_[x][y][z]==PO_EDGE || grid_[x][y][z]==TP_EDGE || grid_[x][y][z]==PO_BURIED || grid_[x][y][z]==TP_BURIED) {
          int found=0;
          for (int i=0;i<dirs*3;i+=3){
            int cx = (int)x+deltas[i];
            int cy = (int)y+deltas[i+1];
            int cz = (int)z+deltas[i+2];
            if (cx<0 || cy<0 || cz<0 || cx>=(int)xdim_ || cy>=(int)ydim_ || cz>=(int)zdim_) continue;

            // If a neighbor is found...
            if (cluster_grid[cx][cy][cz]!=clusters_.clusters_.end()){
              // If this is the first neighbor found, add it to that cluster
              if (!found){
                found=1;
                cluster_grid[x][y][z]=cluster_grid[cx][cy][cz];
                cluster_grid[x][y][z]->add(x,y,z);
              }else{
                // If this cluster differs from the one that you just added, merge them
                if (cluster_grid[cx][cy][cz]!=cluster_grid[x][y][z]){
                  std::list<PCluster>::iterator old = cluster_grid[cx][cy][cz];
                  for (std::list<PCluster::Cxyz>::iterator pit=old->points_.begin(); pit != old->points_.end(); ++pit){
                    cluster_grid[pit->x][pit->y][pit->z]=cluster_grid[x][y][z];
                  }
                  clusters_.join(cluster_grid[x][y][z], old);
                }
              }
            }
          }
          // If none of the neighbors are in a cluster, create a new cluster
          if (!found){
            cluster_grid[x][y][z]=clusters_.add(x,y,z, stepSize_);
          }
        }
      }
    }
  }
}



void PocketGrid::markPocketDepth(core::Real const & surf_d, core::Real const & bur_d){
  core::Size x,y,z;
  for (x=0;x<(xdim_); x++){
    for (y=0;y<(ydim_); y++){
      for (z=0;z<(zdim_); z++){
        if (grid_[x][y][z]==TP_POCKET || grid_[x][y][z]==POCKET) markDepth(x,y,z, surf_d, bur_d);
      }}}
}

void PocketGrid::markEdgeDepth(core::Real const & surf_d, core::Real const & bur_d ){
  for (std::list<PCluster>::iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit){
    for (std::list<PCluster::Cxyz>::iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
      if (grid_[pit->x][pit->y][pit->z]==TP_EDGE || grid_[pit->x][pit->y][pit->z]==PO_EDGE) {
        bool surf=markOneEdgeDepth(pit->x,pit->y,pit->z, surf_d, bur_d, cit->isTarget());
        if (surf){
          std::list<PCluster::Cxyz>::iterator tmp = pit;
          pit--;
          cit->points_.erase(tmp);
        }
      }
    }
    }
    }

    void PocketGrid::markDepth(core::Size x, core::Size y, core::Size z, core::Real const & surf_d, core::Real const & bur_d){
      core::Size ssteps=(core::Size)ceil(surf_d/stepSize_);
      core::Size st;
      core::Size en;
      if (ssteps>x) {
        st=0;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        //clusters_.add(x,y,z, stepSize_);
        return;
      }
      else st=x-ssteps;
      if (x+ssteps>=xdim_) {
        en=xdim_;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        //clusters_.add(x,y,z, stepSize_);
        return ;
      }
      else en=x+ssteps+1;
      for (core::Size i=st;i<en;i++){
        if (grid_[i][y][z]==EMPTY) {
          if (grid_[x][y][z]==POCKET){
            grid_[x][y][z]=PO_SURF;
          }else grid_[x][y][z]=TP_SURF;
          return ;
        }
      }
      if (ssteps>y) {
        st=0;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        //clusters_.add(x,y,z, stepSize_);
        return ;
      }
      else st=y-ssteps;
      if (y+ssteps>=ydim_) {
        en=ydim_;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        //clusters_.add(x,y,z, stepSize_);
        return ;
      }
      else en=y+ssteps+1;
      for (core::Size i=st;i<en;i++){
        if (grid_[x][i][z]==EMPTY) {
          if (grid_[x][y][z]==POCKET){
            grid_[x][y][z]=PO_SURF;
          }else grid_[x][y][z]=TP_SURF;
          return ;
        }
      }
      if (ssteps>z) {
        st=0;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        //clusters_.add(x,y,z, stepSize_);
        return ;
      }
      else st=z-ssteps;
      if (z+ssteps>=zdim_) {
        en=zdim_;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        //clusters_.add(x,y,z, stepSize_);
        return ;
      }
      else en=z+ssteps+1;
      for (core::Size i=st;i<en;i++){
        if (grid_[x][y][i]==EMPTY) {
          if (grid_[x][y][z]==POCKET){
            grid_[x][y][z]=PO_SURF;
          }else grid_[x][y][z]=TP_SURF;
          return ;
        }
      }

      ssteps=(core::Size)ceil(bur_d/stepSize_);
      if (ssteps>x) st=0;
      else st=x-ssteps;
      if (x+ssteps>=xdim_) en=xdim_;
      else en=x+ssteps+1;
      for (core::Size i=st;i<en;i++){
        if (grid_[i][y][z]==HSURFACE || grid_[i][y][z]==PSURFACE || grid_[i][y][z]==PSURFACE || grid_[i][y][z]==T_SURFACE) {
          if (grid_[x][y][z]==POCKET){
            grid_[x][y][z]=PO_BURIED;
          }else grid_[x][y][z]=TP_BURIED;
          //clusters_.add(x,y,z, stepSize_);
          return ;
        }
      }
      if (ssteps>y) st=0;
      else st=y-ssteps;
      if (y+ssteps>=ydim_) en=ydim_;
      else en=y+ssteps+1;
      for (core::Size i=st;i<en;i++){
        if (grid_[x][i][z]==HSURFACE || grid_[x][i][z]==PSURFACE || grid_[x][i][z]==T_SURFACE) {
          if (grid_[x][y][z]==POCKET){
            grid_[x][y][z]=PO_BURIED;
          }else grid_[x][y][z]=TP_BURIED;
          //clusters_.add(x,y,z, stepSize_);
          return ;
        }
      }
      if (ssteps>z) st=0;
      else st=z-ssteps;
      if (z+ssteps>=zdim_) en=zdim_;
      else en=z+ssteps+1;
      for (core::Size i=st;i<en;i++){
        if (grid_[x][y][i]==HSURFACE || grid_[x][y][i]==PSURFACE || grid_[x][y][i]==T_SURFACE) {
          if (grid_[x][y][z]==POCKET){
            grid_[x][y][z]=PO_BURIED;
          }else grid_[x][y][z]=TP_BURIED;
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
      if (ssteps>x)
        xst=0;
      else xst=x-ssteps;
      if (x+ssteps>=xdim_)
        xen=xdim_;
      else xen=x+ssteps+1;
      if (ssteps>y)
        yst=0;
      else yst=y-ssteps;
      if (y+ssteps>=ydim_)
        yen=ydim_;
      else yen=y+ssteps+1;
      if (ssteps>z)
        zst=0;
      else zst=z-ssteps;
      if (z+ssteps>=zdim_)
        zen=zdim_;
      else zen=z+ssteps+1;
      for (core::Size i=xst;i<xen;i++){
        if (grid_[i][y][z]==EMPTY) {
          if (grid_[x][y][z]==PO_EDGE||!isTarget){
            grid_[x][y][z]=PO_SURF;
          }else grid_[x][y][z]=TP_SURF;
          return true;
        }
      }
      for (core::Size i=yst;i<yen;i++){
        if (grid_[x][i][z]==EMPTY) {
          if (grid_[x][y][z]==PO_EDGE||!isTarget){
            grid_[x][y][z]=PO_SURF;
          }else grid_[x][y][z]=TP_SURF;
          return true;
        }
      }
      for (core::Size i=zst;i<zen;i++){
        if (grid_[x][y][i]==EMPTY) {
          if (grid_[x][y][z]==PO_EDGE||!isTarget){
            grid_[x][y][z]=PO_SURF;
          }else grid_[x][y][z]=TP_SURF;
          return true;
        }
      }

      ssteps=(core::Size)ceil(bur_d/stepSize_);
      if (ssteps>x)
        xst=0;
      else xst=x-ssteps;
      if (x+ssteps>=xdim_)
        xen=xdim_;
      else xen=x+ssteps+1;
      if (ssteps>y)
        yst=0;
      else yst=y-ssteps;
      if (y+ssteps>=ydim_)
        yen=ydim_;
      else yen=y+ssteps+1;
      if (ssteps>z)
        zst=0;
      else zst=z-ssteps;
      if (z+ssteps>=zdim_)
        zen=zdim_;
      else zen=z+ssteps+1;
      for (core::Size i=xst;i<xen;i++){
        if (grid_[i][y][z]==HSURFACE || grid_[i][y][z]==PSURFACE ) {
          if (grid_[x][y][z]==POCKET||!isTarget){
            grid_[x][y][z]=PO_BURIED;
          }else grid_[x][y][z]=TP_BURIED;
          return false;
        }
      }
      for (core::Size i=yst;i<yen;i++){
        if (grid_[x][i][z]==HSURFACE || grid_[x][i][z]==PSURFACE ) {
          if (grid_[x][y][z]==POCKET||!isTarget){
            grid_[x][y][z]=PO_BURIED;
          }else grid_[x][y][z]=TP_BURIED;
          return false;
        }
      }

      for (core::Size i=zst;i<zen;i++){
        if (grid_[x][y][i]==HSURFACE || grid_[x][y][i]==PSURFACE ) {
          if (grid_[x][y][z]==POCKET||!isTarget){
            grid_[x][y][z]=PO_BURIED;
          }else grid_[x][y][z]=TP_BURIED;
          return false;
        }
      }

      if (grid_[x][y][z]==PO_EDGE||!isTarget){
        grid_[x][y][z]=POCKET;
      }else grid_[x][y][z]=TP_POCKET;
      return false;

    }

    core::Real PocketGrid::targetPocketVolume(core::Real const & surf_sc, core::Real const & bur_sc) const {
      core::Real vol=0;
      core::Size x,y,z;
      for (x=0;x<(xdim_); x++){
        for (y=0;y<(ydim_); y++){
          for (z=0;z<(zdim_); z++){
            if (grid_[x][y][z]==TP_POCKET) vol+=1;
            if (grid_[x][y][z]==TP_SURF) vol+=surf_sc;
            if (grid_[x][y][z]==TP_BURIED) vol+=bur_sc;
          }}}
      return vol*pow(stepSize_,3);
    }

    core::Real PocketGrid::targetPocketSolventSurface() const {
      core::Real vol=0;
      core::Size x,y,z;
      for (x=0;x<(xdim_); x++){
        for (y=0;y<(ydim_); y++){
          for (z=0;z<(zdim_); z++){
            if (grid_[x][y][z]==TP_POCKET) {
              if (touchesSolvent(x,y,z)) vol+=1;
            }
          }}}
      vol*=pow(stepSize_,2);
      return vol;


    }

    core::Real PocketGrid::targetPocketProteinSurface() const {
      core::Real vol=0;
      core::Size x,y,z;
      for (x=0;x<(xdim_); x++){
        for (y=0;y<(ydim_); y++){
          for (z=0;z<(zdim_); z++){
            if (grid_[x][y][z]==TP_POCKET) {
              if (touchesSurface(x,y,z, true, true)) vol+=1;
            }
          }}}
      vol*=pow(stepSize_,2);
      return vol;


    }

    core::Real PocketGrid::targetPocketHydrophobicProteinSurface() const {
      core::Real vol=0;
      core::Size x,y,z;
      for (x=0;x<(xdim_); x++){
        for (y=0;y<(ydim_); y++){
          for (z=0;z<(zdim_); z++){
            if (grid_[x][y][z]==TP_POCKET) {
              if (touchesSurface(x,y,z, false)) vol+=1;
            }
          }}}
      vol*=pow(stepSize_,2);
      return vol;


    }

    core::Real PocketGrid::targetPocketPolarProteinSurface() const {
      core::Real vol=0;
      core::Size x,y,z;
      for (x=0;x<(xdim_); x++){
        for (y=0;y<(ydim_); y++){
          for (z=0;z<(zdim_); z++){
            if (grid_[x][y][z]==TP_POCKET) {
              if (touchesSurface(x,y,z, true)) vol+=1;
            }
          }}}
      vol*=pow(stepSize_,2);
      return vol;


    }

    core::Real PocketGrid::targetPocketHeuristicScore() const {
      core::Real vol=0;
      core::Size x,y,z;
      for (x=0;x<(xdim_); x++){
        for (y=0;y<(ydim_); y++){
          for (z=0;z<(zdim_); z++){
            if (grid_[x][y][z]==TP_POCKET) {
              if (!touchesSS(x,y,z)){
                core::Size protTouch=touchesPS(x,y,z);
                if (protTouch>0){
                  vol+=3*protTouch;
                }else{
                  vol+=1;
                }
              }
            }
          }}}
      vol*=pow(stepSize_,2);
      return vol;


    }

    void PocketGrid::findClusters(){
//      clusters_.findClusters();
      clusterPockets();
      //mark target clusters
      for (std::list<PCluster>::iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit){
        for (std::list<PCluster::Cxyz>::iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
          for (int x = -1; x <2; x++){
            for (int y = -1; y <2; y++){
              for (int z = -1; z <2; z++){
                if (x==-1) if (pit->x == 0) continue;
                if (y==-1) if (pit->y == 0) continue;
                if (z==-1) if (pit->z == 0) continue;
                if (x==1) if (pit->x == xdim_-1) continue;
                if (y==1) if (pit->y == ydim_-1) continue;
                if (z==1) if (pit->z == zdim_-1) continue;
                if (grid_[pit->x+x][pit->y+y][pit->z+z]==T_SURFACE) {
                  cit->target=true;
                  x=2;y=2;z=2;
                  pit=cit->points_.end();
                  --pit;
                }
              }}}
        }
      }

      //Change target clusters to target types on grid
      for (std::list<PCluster>::iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit){
        if (!cit->isTarget() ) {
          for (std::list<PCluster::Cxyz>::iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
            if (grid_[pit->x][pit->y][pit->z]==TP_POCKET) grid_[pit->x][pit->y][pit->z]=POCKET;
            if (grid_[pit->x][pit->y][pit->z]==TP_SURF) grid_[pit->x][pit->y][pit->z]=PO_SURF;
            if (grid_[pit->x][pit->y][pit->z]==TP_BURIED) grid_[pit->x][pit->y][pit->z]=PO_BURIED;
            if (grid_[pit->x][pit->y][pit->z]==TP_EDGE) grid_[pit->x][pit->y][pit->z]=PO_EDGE;
          }
          continue;
        }

        for (std::list<PCluster::Cxyz>::iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
          if (grid_[pit->x][pit->y][pit->z]==POCKET) grid_[pit->x][pit->y][pit->z]=TP_POCKET;
          if (grid_[pit->x][pit->y][pit->z]==PO_SURF) grid_[pit->x][pit->y][pit->z]=TP_SURF;
          if (grid_[pit->x][pit->y][pit->z]==PO_BURIED) grid_[pit->x][pit->y][pit->z]=TP_BURIED;
          if (grid_[pit->x][pit->y][pit->z]==PO_EDGE) grid_[pit->x][pit->y][pit->z]=TP_EDGE;
        }
      }
    }

    core::Real PocketGrid::largestTargetPocketVolume() {
      return clusters_.getLargestClusterSize( stepSize_, minPockSize_);
    }

    core::Real PocketGrid::netTargetPocketVolume() {
      if (maxPockSize_){
        return std::min(maxPockSize_, clusters_.getNetClusterSize( stepSize_ , minPockSize_));
      }else return clusters_.getNetClusterSize( stepSize_ , minPockSize_);
    }


    bool PocketGrid::isTooSmall() const {
      core::Size x,y,z;
      if (xdim_==0||ydim_==0||zdim_==0) return false;
      if ( grid_[0][0][0] == TP_POCKET) return true;
      if ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_POCKET) return true;
      if ( grid_[0][0][0] == TP_SURF) return true;
      if ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_SURF) return true;
      if ( grid_[0][0][0] == TP_BURIED) return true;
      if ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_BURIED) return true;
      if ( grid_[0][0][0] == TP_EDGE) return true;
      if ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_EDGE) return true;

      for (x=0, y=0, z=0; x<xdim_ || y<ydim_ || z<zdim_;x++,y++,z++){
        if (x==xdim_) x--;
        if (y==ydim_) y--;
        if (z==zdim_) z--;
        if ( grid_[0][y][z] == TP_POCKET) return true;
        if ( grid_[x][0][z] == TP_POCKET) return true;
        if ( grid_[x][y][0] == TP_POCKET) return true;
        if ( grid_[xdim_-1][y][z] == TP_POCKET) return true;
        if ( grid_[x][ydim_-1][z] == TP_POCKET) return true;
        if ( grid_[x][y][zdim_-1] == TP_POCKET) return true;
        if ( grid_[0][y][z] == TP_SURF) return true;
        if ( grid_[x][0][z] == TP_SURF) return true;
        if ( grid_[x][y][0] == TP_SURF) return true;
        if ( grid_[xdim_-1][y][z] == TP_SURF) return true;
        if ( grid_[x][ydim_-1][z] == TP_SURF) return true;
        if ( grid_[x][y][zdim_-1] == TP_SURF) return true;
        if ( grid_[0][y][z] == TP_BURIED) return true;
        if ( grid_[x][0][z] == TP_BURIED) return true;
        if ( grid_[x][y][0] == TP_BURIED) return true;
        if ( grid_[xdim_-1][y][z] == TP_BURIED) return true;
        if ( grid_[x][ydim_-1][z] == TP_BURIED) return true;
        if ( grid_[x][y][zdim_-1] == TP_BURIED) return true;
        if ( grid_[0][y][z] == TP_EDGE) return true;
        if ( grid_[x][0][z] == TP_EDGE) return true;
        if ( grid_[x][y][0] == TP_EDGE) return true;
        if ( grid_[xdim_-1][y][z] == TP_EDGE) return true;
        if ( grid_[x][ydim_-1][z] == TP_EDGE) return true;
        if ( grid_[x][y][zdim_-1] == TP_EDGE) return true;
      }


      return false;
    }

    core::Vector PocketGrid::whatIsTooSmall() const {
      core::Vector dims;
      dims(1)=0;
      dims(2)=0;
      dims(3)=0;
      core::Size x,y,z;
      if (xdim_==0||ydim_==0||zdim_==0) return dims;
      if (( grid_[0][0][0] == TP_POCKET) || ( grid_[0][0][0] == TP_SURF) || ( grid_[0][0][0] == TP_BURIED) || ( grid_[0][0][0] == TP_EDGE) || ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_POCKET) || ( grid_[xdim_-1][ydim_-1][zdim_-1]) || ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_BURIED) || ( grid_[xdim_-1][ydim_-1][zdim_-1] == TP_EDGE)) {
        dims(1)=1;
        dims(2)=1;
        dims(3)=1;
        return dims;
      }

      for (x=0, y=0, z=0; x<xdim_ || y<ydim_ || z<zdim_;x++,y++,z++){
        if (x==xdim_) x--;
        if (y==ydim_) y--;
        if (z==zdim_) z--;
        if (( grid_[0][y][z] == TP_POCKET) || ( grid_[0][y][z] == TP_BURIED) || ( grid_[0][y][z] == TP_EDGE) || ( grid_[xdim_-1][y][z] == TP_POCKET) || ( grid_[xdim_-1][y][z] == TP_SURF) || ( grid_[xdim_-1][y][z] == TP_BURIED) || ( grid_[xdim_-1][y][z] == TP_EDGE)){
          dims(1)=1;
        }
        if (( grid_[x][0][z] == TP_POCKET) || ( grid_[x][0][z] == TP_SURF) || ( grid_[x][0][z] == TP_BURIED) || ( grid_[x][0][z] == TP_EDGE) || ( grid_[x][ydim_-1][z] == TP_POCKET) || ( grid_[x][ydim_-1][z] == TP_SURF) ||  ( grid_[x][ydim_-1][z] == TP_BURIED) || ( grid_[x][ydim_-1][z] == TP_EDGE)) {
          dims(2)=1;
        }
        if (( grid_[x][y][0] == TP_POCKET) || ( grid_[x][y][0] == TP_SURF) || ( grid_[x][y][0] == TP_BURIED) || ( grid_[x][y][0] == TP_EDGE) || ( grid_[x][y][zdim_-1] == TP_POCKET) || ( grid_[x][y][zdim_-1] == TP_SURF) || ( grid_[x][y][zdim_-1] == TP_BURIED) || ( grid_[x][y][zdim_-1] == TP_EDGE)){
          dims(3)=1;
        }
      }


      return dims;
    }

    bool PocketGrid::autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::pose::Pose const & inPose ) {
      core::scoring::constraints::ConformationXYZ const xyz_func( inPose.conformation() );
      return autoexpanding_pocket_eval(central_rsd, xyz_func, inPose.total_residue());
    }

    bool PocketGrid::autoexpanding_pocket_eval( std::vector< core::conformation::ResidueOP > const & central_rsds, core::scoring::constraints::XYZ_Func const & xyz_func, Size const total_residues ) {
      return autoexpanding_pocket_eval(*central_rsds[0], xyz_func, total_residues);
}

    bool PocketGrid::autoexpanding_pocket_eval( std::vector< core::conformation::ResidueOP > const & central_rsds, core::pose::Pose const & inPose ) {
      core::scoring::constraints::ConformationXYZ const xyz_func( inPose.conformation() );
      return autoexpanding_pocket_eval(*central_rsds[0], xyz_func, inPose.total_residue());
    }

    bool PocketGrid::autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::scoring::constraints::XYZ_Func const & xyz_func, Size const total_residues ) {

      bool too_small=true;
      while (too_small){
        recenter( central_rsd );

        for ( Size j = 1, resnum = total_residues; j <= resnum; ++j ) {
          core::conformation::Residue const & rsd( xyz_func.residue(j) );
          bool target=false;
          if ( j == central_rsd.seqpos() ) target=true;
					core::Size total_atoms(0);
					using namespace basic::options;
					if (option[ OptionKeys::fingerprint::include_hydrogens ]()){
						total_atoms = rsd.natoms();
					} else {
						total_atoms = rsd.nheavyatoms();
					}
          for(Size i = 1, i_end = total_atoms; i <= i_end; ++i) {
            bool target_res=target;
            if (side_chains_only_){
              if ((central_rsd.atom(i).type()>=18)&&(central_rsd.atom(i).type()<=21)){
                target_res=false;
              }
            }
						numeric::xyzVector<core::Real> rpoint = rotatePoint(rsd.atom(i).xyz().x(),rsd.atom(i).xyz().y(),rsd.atom(i).xyz().z());
            mark(rpoint, rsd.atom_type(i).lj_radius(),probe_rad_, rsd.is_polar(), target_res);
          }
        }

        findPockets(0, maxLen_);
        if (markpsp_){
          findPSP(0,maxLen_);
        }

        if (marksps_){
          findSPS(0,ceil((xdim_+ydim_+zdim_)/stepSize_));
        }
        markPocketDepth(surf_dist_, bur_dist_);

        findClusters();

        //  Code to check edge of grid goes here
        if (isTooSmall()){
          core::Vector dims=whatIsTooSmall();
          if ((size_x_ + dims(1) > limit_x_) || (size_y_ + dims(2) > limit_y_) ||(size_z_ + dims(3) > limit_z_)){
            return false;
          }
          size_x_ = size_x_ + dims(1);
          size_y_ = size_y_ + dims(2);
          size_z_ = size_z_ + dims(3);
          initialize(central_rsd, size_x_, size_y_, size_z_, spacing_, markpsp_, marksps_);
        }else{
          too_small=false;
        }
      }

      markEdgeDepth(surf_dist_, bur_dist_);
      return true;
    }

	bool PocketGrid::DARC_autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::pose::Pose const & inPose, numeric::xyzVector<core::Real> grid_center ) {
		core::scoring::constraints::ConformationXYZ const xyz_func( inPose.conformation() );
		return DARC_autoexpanding_pocket_eval(central_rsd, xyz_func, inPose.total_residue(), grid_center);
	}

	bool PocketGrid::DARC_autoexpanding_pocket_eval( core::conformation::Residue const & central_rsd, core::scoring::constraints::XYZ_Func const & xyz_func, Size const total_residues, numeric::xyzVector<core::Real> grid_center ) {
		bool too_small=true;
		while (too_small){
			recenter( grid_center.x(), grid_center.y(), grid_center.z() );
			for ( Size j = 1, resnum = total_residues; j <= resnum; ++j ) {
          core::conformation::Residue const & rsd( xyz_func.residue(j) );
          bool target=false;
          if ( j == central_rsd.seqpos() ) target=true;
					core::Size total_atoms(0);
					using namespace basic::options;
					if (option[ OptionKeys::fingerprint::include_hydrogens ]()){
						total_atoms = rsd.natoms();
					} else {
						total_atoms = rsd.nheavyatoms();
					}
					for(Size i = 1, i_end = total_atoms; i <= i_end; ++i) {
            bool target_res=target;
            if (side_chains_only_){
              if ((central_rsd.atom(i).type()>=18)&&(central_rsd.atom(i).type()<=21)){
                target_res=false;
              }
            }
						numeric::xyzVector<core::Real> rpoint = rotatePoint(rsd.atom(i).xyz().x(),rsd.atom(i).xyz().y(),rsd.atom(i).xyz().z());
            mark(rpoint, rsd.atom_type(i).lj_radius(),probe_rad_, rsd.is_polar(), target_res);
          }
			}
			findPockets(0, maxLen_);
			if (markpsp_){
				findPSP(0,maxLen_);
			}
			if (marksps_){
				findSPS(0,ceil((xdim_+ydim_+zdim_)/stepSize_));
			}
			markPocketDepth(surf_dist_, bur_dist_);
			findClusters();
			//  Code to check edge of grid goes here
			if (isTooSmall()){
				core::Vector dims=whatIsTooSmall();
				if ((size_x_ + dims(1) > limit_x_) || (size_y_ + dims(2) > limit_y_) ||(size_z_ + dims(3) > limit_z_)){
					return false;
				}
				size_x_ = size_x_ + dims(1);
				size_y_ = size_y_ + dims(2);
				size_z_ = size_z_ + dims(3);
				initialize(grid_center.x(), grid_center.y(), grid_center.z(), size_x_, size_y_, size_z_, spacing_, markpsp_, marksps_);
			}else{
				too_small=false;
			}
		}
		markEdgeDepth(surf_dist_, bur_dist_);
		return true;
	}

	bool PocketGrid::DARC_pocket_eval( core::conformation::Residue const & central_rsd, core::pose::Pose const & inPose, numeric::xyzVector<core::Real> grid_center ) {
		core::scoring::constraints::ConformationXYZ const xyz_func( inPose.conformation() );
		return DARC_pocket_eval(central_rsd, xyz_func, inPose.total_residue(), grid_center);
	}

	bool PocketGrid::DARC_pocket_eval(
			core::conformation::Residue const & central_rsd,
			core::scoring::constraints::XYZ_Func const & xyz_func,
			Size const total_residues,
			numeric::xyzVector<core::Real> /*grid_center*/ )
	{
		//bool too_small=true;
		//while (too_small){
		//recenter( grid_center.x(), grid_center.y(), grid_center.z() );
		for ( Size j = 1, resnum = total_residues; j <= resnum; ++j ) {
          core::conformation::Residue const & rsd( xyz_func.residue(j) );
          bool target=false;
          if ( j == central_rsd.seqpos() ) target=true;
					core::Size total_atoms(0);
					using namespace basic::options;
					if (option[ OptionKeys::fingerprint::include_hydrogens ]()){
						total_atoms = rsd.natoms();
					} else {
						total_atoms = rsd.nheavyatoms();
					}
          for(Size i = 1, i_end = total_atoms; i <= i_end; ++i) {
            bool target_res=target;
            if (side_chains_only_){
              if ((central_rsd.atom(i).type()>=18)&&(central_rsd.atom(i).type()<=21)){
                target_res=false;
              }
            }
						numeric::xyzVector<core::Real> rpoint = rotatePoint(rsd.atom(i).xyz().x(),rsd.atom(i).xyz().y(),rsd.atom(i).xyz().z());
            mark(rpoint, rsd.atom_type(i).lj_radius(),probe_rad_, rsd.is_polar(), target_res);
          }
		}
		findPockets(0, maxLen_);
		if (markpsp_){
			findPSP(0,maxLen_);
		}
		if (marksps_){
			findSPS(0,ceil((xdim_+ydim_+zdim_)/stepSize_));
		}
		markPocketDepth(surf_dist_, bur_dist_);
		findClusters();
		/*  Code to check edge of grid goes here
				if (isTooSmall()){
				core::Vector dims=whatIsTooSmall();
				if ((size_x_ + dims(1) > limit_x_) || (size_y_ + dims(2) > limit_y_) ||(size_z_ + dims(3) > limit_z_)){
					return false;
					}
				size_x_ = size_x_ + dims(1);
				size_y_ = size_y_ + dims(2);
				size_z_ = size_z_ + dims(3);
				initialize(grid_center.x(), grid_center.y(), grid_center.z(), size_x_, size_y_, size_z_, spacing_, markpsp_, marksps_);
				}else{
				too_small=false;
				}
				}*/
		markEdgeDepth(surf_dist_, bur_dist_);
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
	for (std::list<PCluster>::const_iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit){
		if (cit->points_.size()*pow(stepSize_,3)<minPockSize_) continue; // this is a smallpocket
		for (std::list<PCluster::Cxyz>::const_iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
			if ( (grid_[pit->x][pit->y][pit->z]==TP_POCKET) || (grid_[pit->x][pit->y][pit->z]==TP_SURF) || (grid_[pit->x][pit->y][pit->z]==TP_EDGE) || (grid_[pit->x][pit->y][pit->z]==TP_BURIED) ) {
				CoM_x += pit->x*stepSize_+xcorn_;
				CoM_y += pit->y*stepSize_+ycorn_;
				CoM_z += pit->z*stepSize_+zcorn_;
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
	for (std::list<PCluster>::const_iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit){
		if (cit->points_.size()*pow(stepSize_,3)<minPockSize_) continue; // this is a smallpocket
		for (std::list<PCluster::Cxyz>::const_iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
			if ( (grid_[pit->x][pit->y][pit->z]==TP_POCKET) || (grid_[pit->x][pit->y][pit->z]==TP_SURF) || (grid_[pit->x][pit->y][pit->z]==TP_EDGE) || (grid_[pit->x][pit->y][pit->z]==TP_BURIED) ) {
				// add this pocket point to the Eigen matrix (subtracting out the pocket CoM)
				pocket_coors_matrix(pocket_inx,0) = pit->x*stepSize_+xcorn_ - CoM_x;
				pocket_coors_matrix(pocket_inx,1) = pit->y*stepSize_+ycorn_ - CoM_y;
				pocket_coors_matrix(pocket_inx,2) = pit->z*stepSize_+zcorn_ - CoM_z;
				++pocket_inx;
			}
		}
	}

	// do the SVD calculation
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(pocket_coors_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV );
	svd.computeV();
	Eigen::MatrixXf svd_result(3,3);
	svd_result = svd.matrixV();
	//	std::cout << "Matrix V is: " << std::endl << svd.matrixV() << std::endl;
	//	MatrixXf transformed = A * svd.matrixV();

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
	core::Size output_res_num = 1.;

	// Loop over all points in the template pocket
	core::Size template_num_points = 0;
	for (std::list<PCluster>::const_iterator cit=template_pocket.clusters_.clusters_.begin(); cit != template_pocket.clusters_.clusters_.end(); ++cit){
		if (cit->points_.size()*pow(template_pocket.stepSize_,3)<template_pocket.minPockSize_) continue; // this is a smallpocket
		for (std::list<PCluster::Cxyz>::const_iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
			if ( (template_pocket.grid_[pit->x][pit->y][pit->z]==template_pocket.TP_POCKET) || (template_pocket.grid_[pit->x][pit->y][pit->z]==template_pocket.TP_SURF) || (template_pocket.grid_[pit->x][pit->y][pit->z]==template_pocket.TP_EDGE) || (template_pocket.grid_[pit->x][pit->y][pit->z]==template_pocket.TP_BURIED) ) {

					++template_num_points;

				// For this template pocket point, get the cartesian coors and convert to comparison's grid
				core::Real const template_x = pit->x*template_pocket.stepSize_+template_pocket.xcorn_;
				core::Real const template_y = pit->y*template_pocket.stepSize_+template_pocket.ycorn_;
				core::Real const template_z = pit->z*template_pocket.stepSize_+template_pocket.zcorn_;

				core::Size const self_x_index = floor( ( ( template_x - xcorn_ ) / stepSize_ ) + 0.5 );
				core::Size const self_y_index = floor( ( ( template_y - ycorn_ ) / stepSize_ ) + 0.5 );
				core::Size const self_z_index = floor( ( ( template_z - zcorn_ ) / stepSize_ ) + 0.5 );

				//Check to see if template point is within comparison's grid range
				if ( ( self_x_index < xdim_ ) && ( self_y_index < ydim_ ) && ( self_z_index < zdim_ ) ) {
					// note: no need to check that self_index >= 0 since it's of type "Size"
					//				if ( ( self_x_index >= 0 ) && ( self_x_index < xdim_ ) && ( self_y_index >= 0 ) && ( self_y_index < ydim_ ) && ( self_z_index >= 0 ) && ( self_z_index < zdim_ ) ) {

					// jk note: later, we could actually match sub-pocket types - will this give better results??
					//See if template point is pocket for comparison
					if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_POCKET ) || ( grid_[self_x_index][self_y_index][self_z_index] == TP_SURF ) || ( grid_[self_x_index][self_y_index][self_z_index] == TP_EDGE ) || ( grid_[self_x_index][self_y_index][self_z_index] == TP_BURIED ) ) {

						core::Real curr_weight = partial_match_weight;
						if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_POCKET ) && (template_pocket.grid_[pit->x][pit->y][pit->z]==template_pocket.TP_POCKET) ) {
							curr_weight = full_match_weight;
						} else if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_SURF ) && (template_pocket.grid_[pit->x][pit->y][pit->z]==template_pocket.TP_SURF) ) {
							curr_weight = full_match_weight;
						} else if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_EDGE ) && (template_pocket.grid_[pit->x][pit->y][pit->z]==template_pocket.TP_EDGE) ) {
							curr_weight = full_match_weight;
						} else if ( ( grid_[self_x_index][self_y_index][self_z_index] == TP_BURIED ) && (template_pocket.grid_[pit->x][pit->y][pit->z]==template_pocket.TP_BURIED) ) {
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
	for (std::list<PCluster>::const_iterator cit=clusters_.clusters_.begin(); cit != clusters_.clusters_.end(); ++cit){
		if (cit->points_.size()*pow(stepSize_,3)<minPockSize_) continue; // this is a smallpocket
		for (std::list<PCluster::Cxyz>::const_iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
			if ( (grid_[pit->x][pit->y][pit->z]==TP_POCKET) || (grid_[pit->x][pit->y][pit->z]==TP_SURF) || (grid_[pit->x][pit->y][pit->z]==TP_EDGE) || (grid_[pit->x][pit->y][pit->z]==TP_BURIED) ) {
					++self_num_points;
			}
		}
	}

	// note: non-trivial to average number of points in each pocket, since we don't know how many points in the
	match_score = 2. * match_score / ( template_num_points + self_num_points );
	TR << "PocketGrid match score = " << match_score << std::endl;

	return 1. - sqrt( match_score );
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
	for (core::Size x=0;x<(ext_grd.xdim_); x++){
		for (core::Size y=0;y<(ext_grd.ydim_); y++){
			for (core::Size z=0;z<(ext_grd.zdim_); z++){
				if ( (ext_grd.grid_[x][y][z] == ext_grd.HSURFACE) || (ext_grd.grid_[x][y][z] == ext_grd.PSURFACE) || (ext_grd.grid_[x][y][z] == ext_grd.T_SURFACE) ){
					if (x != 0) {
						searchxmin=x-1;
					}else{
						searchxmin=x;
					}
					if (x != ext_grd.xdim_-1) {
						searchxmax=x+1;
					}else{
						searchxmax=x;
					}
					if (y != 0) {
						searchymin=y-1;
					}else{
						searchymin=y;
					}
					if (y != ext_grd.ydim_-1) {
						searchymax=y+1;
					}else{
						searchymax=y;
					}
					if (z != 0) {
						searchzmin=z-1;
					}else{
						searchzmin=z;
					}
					if (z != ext_grd.zdim_-1) {
						searchzmax=z+1;
					}else{
						searchzmax=z;
					}
					bool found = false;
					for (xx = searchxmin; xx<=searchxmax; ++xx){
						for (yy = searchymin; yy<=searchymax; ++yy){
							for (zz = searchzmin; zz<=searchzmax; ++zz){
								if (xx==x && yy==y && zz==z) continue;
								if ( ( ext_grd.grid_[xx][yy][zz]==ext_grd.EMPTY ) || (ext_grd.grid_[xx][yy][zz] == ext_grd.TP_POCKET) || (ext_grd.grid_[xx][yy][zz] == ext_grd.TP_EDGE) || (ext_grd.grid_[xx][yy][zz] == ext_grd.TP_BURIED) || (ext_grd.grid_[xx][yy][zz] == ext_grd.TP_SURF) || (ext_grd.grid_[xx][yy][zz] == ext_grd.POCKET) || (ext_grd.grid_[xx][yy][zz] == ext_grd.PO_BURIED) || (ext_grd.grid_[xx][yy][zz] == ext_grd.PO_SURF) || (ext_grd.grid_[xx][yy][zz] == ext_grd.PO_EDGE) ) {
									extra_coord.x() = (x*ext_grd.stepSize_+ext_grd.xcorn_);
									extra_coord.y() = (y*ext_grd.stepSize_+ext_grd.ycorn_);
									extra_coord.z() = (z*ext_grd.stepSize_+ext_grd.zcorn_);
									found = true;
									bool nearby = false;
									for (std::list< numeric::xyzVector<core::Real> >::const_iterator eg = eggshell_coord_list.begin(); eg != eggshell_coord_list.end(); ++eg) {
										if (extra_coord == *eg) {
											nearby = false;
											break;
										}
										if ( extra_coord.distance(*eg) < extra_dist ) {
											nearby = true;
										}
									}
									if (nearby){
										extra_coord_list_.push_back(extra_coord);
									}
								}
								if (found) break;
							}
							if (found) break;
						}
						if (found) break;
					}
				}
			}
		}
	}

//Print grid points into a PDB file
	using namespace basic::options;
	if (option[ OptionKeys::pocket_grid::print_grid ]()){
		write_grid_to_pdb(ext_grd, "ext_grd.pdb");
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
	for (core::Size x=0;x<(gr.xdim_); x++){
		for (core::Size y=0;y<(gr.ydim_); y++){
			for (core::Size z=0;z<(gr.zdim_); z++){
				if ( (gr.grid_[x][y][z] == gr.HSURFACE) || (gr.grid_[x][y][z] == gr.PSURFACE) || (gr.grid_[x][y][z] == gr.T_SURFACE) /*||(gr.grid_[x][y][z]==gr.TP_SURF)*/ ) {
					if (x != 0) {
						searchxmin=x-1;
					}else{
						searchxmin=x;
					}
					if (x != gr.xdim_-1) {
						searchxmax=x+1;
					}else{
						searchxmax=x;
					}
					if (y != 0) {
						searchymin=y-1;
					}else{
						searchymin=y;
					}
					if (y != gr.ydim_-1) {
						searchymax=y+1;
					}else{
						searchymax=y;
					}
					if (z != 0) {
						searchzmin=z-1;
					}else{
						searchzmin=z;
					}
					if (z != gr.zdim_-1) {
						searchzmax=z+1;
					}else{
						searchzmax=z;
					}
					bool found=false;
					for (xx = searchxmin; (xx<=searchxmax) && ! found; xx++){
						for (yy = searchymin; (yy<=searchymax) && ! found; yy++){
							for (zz = searchzmin; (zz<=searchzmax) && ! found; zz++){
								if (xx==x && yy==y && zz==z) continue;
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
	eggshell_CoM_ =	calculate_center_xyz(eggshell_coord_list_);

	//Print grid points into a PDB file
	using namespace basic::options;
	if (option[ OptionKeys::pocket_grid::print_grid ]()){
		write_grid_to_pdb(gr, "egg_grd.pdb");
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
	for (core::Size x=0;x<xdim_; x++) {
    for (core::Size y=0;y<ydim_; y++) {
      for (core::Size z=0;z<zdim_; z++) {
				if ( ( grid_[x][y][z]==EGGSHELL ) ) {
					outstream << x << " " << y << " " << z << std::endl;
					//					outstream << x << " " << y << " " << z << " " << x*stepSize_+xcorn_ << " " << y*stepSize_+ycorn_ << " " << z*stepSize_+zcorn_ << std::endl;
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
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = eggshell_coord_list_.begin(); pd != eggshell_coord_list_.end(); ++pd) {
		outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   EGG A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->z()<<std::endl;
	}
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = extra_coord_list_.begin(); pd != extra_coord_list_.end(); ++pd) {
		outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  O   EXT B   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->z()<<std::endl;
	}
	outPDB_stream.close();
	outPDB_stream.clear();

}

void EggshellGrid::write_grid_to_pdb( const PocketGrid& gr, std::string const & output_pdbname ) const {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_pdbname, std::ios::out);

	for (core::Size x=0;x<(gr.xdim_); x++){
    for (core::Size y=0;y<(gr.ydim_); y++){
      for (core::Size z=0;z<(gr.zdim_); z++){
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
		for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = pocketshell_coord_list.begin(); pd != pocketshell_coord_list.end(); ++pd) {
			pocketshell_CoM.x() += pd->x();
			pocketshell_CoM.y() += pd->y();
			pocketshell_CoM.z() += pd->z();
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
	core::Size output_res_num = 1.;

	// Loop over all points in the template eggshell
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = template_eggshell.eggshell_coord_list().begin(); pd != template_eggshell.eggshell_coord_list().end(); ++pd) {

    // For this template eggshell point, get the cartesian coors and convert to comparison's grid
		core::Real const template_x = pd->x();
		core::Real const template_y = pd->y();
		core::Real const template_z = pd->z();
		core::Size const self_x_index = floor( ( ( template_x - xcorn_ ) / stepSize_ ) + 0.5 );
		core::Size const self_y_index = floor( ( ( template_y - ycorn_ ) / stepSize_ ) + 0.5 );
		core::Size const self_z_index = floor( ( ( template_z - zcorn_ ) / stepSize_ ) + 0.5 );

		//Check to see if template point is within comparison's grid range
		if ( ( self_x_index < xdim_ ) && ( self_y_index < ydim_ ) && ( self_z_index < zdim_ ) ) {
			// note: no need to check that self_index >= 0 since it's of type "Size"
			//		if ( ( self_x_index >= 0 ) && ( self_x_index < xdim_ ) && ( self_y_index >= 0 ) && ( self_y_index < ydim_ ) && ( self_z_index >= 0 ) && ( self_z_index < zdim_ ) ) {
			//See if template point is eggshell for comparison
			if ( grid_[self_x_index][self_y_index][self_z_index] == EGGSHELL ) {
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

	match_score = 2. * match_score / ( template_eggshell.eggshell_coord_list().size() + eggshell_coord_list().size() );
	TR << "Eggshell match score = " << match_score << std::endl;

	return 1. - match_score;
}


  } // pockets
} // protocols
