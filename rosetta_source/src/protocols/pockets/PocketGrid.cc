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
#include <core/scoring/ScoreType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <basic/Tracer.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>
#include <fstream>
#include <iostream>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility Headers
#include <algorithm>
// Auto-header: duplicate removed #include <iostream>
#include <iomanip>
// Auto-header: duplicate removed #include <fstream>
#include <ostream>
// Auto-header: duplicate removed #include <string>
#include <sstream>
#include <cmath>
#include <map>

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

void PClusterSet::add(core::Size x, core::Size y, core::Size z, core::Real step){
  PCluster tmp(x,y,z, step);
  clusters_.push_back(tmp);
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
  core::Size thr2=(core::Size)std::floor(((core::Real)(thr)/stepSize_)/sqrt(3)+.5);
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
  core::Size max2=(core::Size)std::floor((max/stepSize_)/sqrt(2)+.5);
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

void PocketGrid::searchX(core::Size thr, bool psp){
  for (core::Size z=0; z<zdim_;z++){
    for (core::Size y=0;y<ydim_;y++){
      bool surf=false;
      bool t_surf=false;
      core::Size count=0;
      core::Size pcount=0;
      int start=-1;
      for (core::Size x=0;x<xdim_;x++){
        if (surf){
          if (psp){
            if (grid_[x][y][z]==EMPTY){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==POCKET || grid_[x][y][z]==TP_POCKET){
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              if (count<= thr && pcount>0){
                if (markpsp_){
                  for (int c=count;c>0; c--){
                    if (t_surf) grid_[x-c][y][z]=TP_POCKET;
                    else if (grid_[x][y][z]!=TP_POCKET)  grid_[x-c][y][z]=POCKET;
                  }
                }else{
                  for (int c=count;c>0; c--)
                    if (grid_[x-c][y][z]==EMPTY) grid_[x-c][y][z]=PSP;
                }
              }
              start=x;
              count=0;
              if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE) pcount=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }

          }else{
            if (grid_[x][y][z]==EMPTY||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              if (count> thr){
                for (int c=count;c>0; c--){
                  if (t_surf) grid_[x-c][y][z]=TP_POCKET;
                  else if (grid_[x][y][z]!=TP_POCKET)  grid_[x-c][y][z]=POCKET;
                }
              }
              start=x;
              count=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }
          }
        }else{
          if (psp){
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              pcount=0;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              surf=true;
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              count=0;
              start=x;
            }

          }else{
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              surf=true;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              count=0;
              start=x;
            }
          }
        }
      }
    }}
}

void PocketGrid::searchY(core::Size thr, bool psp){
  for (core::Size x=0; x<xdim_;x++){
    for (core::Size z=0;z<zdim_;z++){
      bool surf=false;
      bool t_surf=false;
      core::Size count=0;
      core::Size pcount=0;
      int start=-1;
      for (core::Size y=0;y<ydim_;y++){
        if (surf){
          if (psp){
            if (grid_[x][y][z]==EMPTY){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==POCKET || grid_[x][y][z]==TP_POCKET){
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;

              if (count<= thr && pcount>0){
                if (markpsp_){
                  for (int c=count;c>0; c--){
                    if (t_surf) grid_[x][y-c][z]=TP_POCKET;
                    else if (grid_[x][y-c][z]!=TP_POCKET)  grid_[x][y-c][z]=POCKET;
                  }
                }else{
                  for (int c=count;c>0; c--)
                    if (grid_[x][y-c][z]==EMPTY) grid_[x][y-c][z]=PSP;
                }
              }
              start=x;
              count=0;
              if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE) pcount=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }

          }else{
            if (grid_[x][y][z]==EMPTY||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              if (count> thr){
                for (int c=count;c>0; c--){
                  if (t_surf) grid_[x][y-c][z]=TP_POCKET;
                  else if (grid_[x][y][z]!=TP_POCKET)  grid_[x][y-c][z]=POCKET;
                }
              }
              start=y;
              count=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }
          }
        }else{
          if (psp){
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              pcount=0;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              surf=true;
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              count=0;
              start=x;
            }

          }else{
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              surf=true;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              count=0;
              start=y;
            }
          }
        }
      }
    }}

}

void PocketGrid::searchZ(core::Size thr, bool psp){
  for (core::Size x=0; x<xdim_;x++){
    for (core::Size y=0;y<ydim_;y++){
      bool surf=false;
      bool t_surf=false;
      core::Size count=0;
      core::Size pcount=0;
      int start=-1;
      for (core::Size z=0;z<zdim_;z++){
        if (surf){
          if (psp){
            if (grid_[x][y][z]==EMPTY){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==POCKET || grid_[x][y][z]==TP_POCKET){
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              if (count<= thr && pcount>0){
                if (markpsp_){
                  for (int c=count;c>0; c--){
                    if (t_surf) grid_[x][y][z-c]=TP_POCKET;
                    else if (grid_[x][y][z-c]!=TP_POCKET)  grid_[x][y][z-c]=POCKET;
                  }
                }else{
                  for (int c=count;c>0; c--)
                    if (grid_[x][y][z-c]==EMPTY) grid_[x][y][z-c]=PSP;
                  ;
                }
              }
              start=x;
              count=0;
              if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE) pcount=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }

          }else{
            if (grid_[x][y][z]==EMPTY||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              if (count> thr){
                for (int c=count;c>0; c--){
                  if (t_surf) grid_[x][y][z-c]=TP_POCKET;
                  else if (grid_[x][y][z]!=TP_POCKET)  grid_[x][y][z-c]=POCKET;
                }
              }
              start=z;
              count=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }
          }
        }else{
          if (psp){
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              pcount=0;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              surf=true;
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              count=0;
              start=x;
            }

          }else{
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              surf=true;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              count=0;
              start=z;
            }
          }
        }
      }
    }}

}

void PocketGrid::searchD1(core::Size thr, bool psp){
  int diag=(int)floor(sqrt(pow(zdim_,2)+pow(zdim_,2)+pow(zdim_,2)));
  int zstart=0;
  int xfact=1,yfact=1,zfact=1;
  int xend=0,yend=0,zend=0;
  if (xfact<0) xend=xdim_-1;
  if (yfact<0) yend=ydim_-1;
  if (zfact<0) zend=zdim_-1;
  for (core::Size ix=0;ix<xdim_;ix++){
    for (core::Size iy=0;iy<ydim_;iy++){
      bool surf=false;
      bool t_surf=false;
      core::Size count=0;
      core::Size pcount=0;
      core::Size count2=0;
      int start=-1;
      bool empty=false;
      int x=ix;
      int y=iy;
      int z=zstart;
      for (int c=0;c<diag;c++){
        x+=xfact;
        while (x<0) {x+=xdim_;count=0;}
        if (x>=(int)xdim_) {x=x%xdim_;count=0;}
        y+=yfact;
        while (y<0) {y+=ydim_;count=0;}
        if (y>=(int)ydim_) {y=y%ydim_;count=0;}
        z+=zfact;
        while (z<0) {z+=zdim_;count=0;}
        if (z>=(int)zdim_) {z=z%zdim_;count=0;}
        count2++;
        if (x==xend||y==yend||z==zend) surf=false;
        if (surf){
          if (psp){
            if (grid_[x][y][z]==EMPTY){
              count++;
            }else if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET||grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              if (count<= thr && pcount>0){
                if (markpsp_){
                  for (int c=count;c>0; c--){
                    if (t_surf) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=TP_POCKET;
                    else if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]!=TP_POCKET) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=POCKET;
                  }
                }else{
                  for (int c=count;c>0; c--)
                    if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]==EMPTY) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=PSP;
                }
              }
              start=c;
              count=0;
              if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE) pcount=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }

          }else{
            if (grid_[x][y][z]==EMPTY||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              count++;
              empty=true;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              if (count> thr){
                for (int c=count;c>0; c--){
                  if (t_surf) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=TP_POCKET;
                  else if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]!=TP_POCKET) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=POCKET;
                }

              }
              start=c;
              count=0;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              else t_surf=false;
              empty=false;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
              empty=false;
            }
          }
        }else{
          if (psp){
            if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET||grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              pcount=0;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              surf=true;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              count=0;
              start=c;
            }

          }else{
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              surf=true;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              count=0;
              start=c;
              empty=false;
            }
          }
        }
      }
    }}
}

void PocketGrid::searchD2(core::Size thr, bool psp){
  int diag=(int)floor(sqrt(pow(zdim_,2)+pow(zdim_,2)+pow(zdim_,2)));
  int zstart=0;
  int xfact=-1,yfact=1,zfact=1;
  int xend=0,yend=0,zend=0;
  if (xfact<0) xend=xdim_-1;
  if (yfact<0) yend=ydim_-1;
  if (zfact<0) zend=zdim_-1;
  for (core::Size ix=0;ix<xdim_;ix++){
    for (core::Size iy=0;iy<ydim_;iy++){
      bool surf=false;
      bool t_surf=false;
      core::Size count=0;
      core::Size pcount=0;
      core::Size count2=0;
      int start=-1;
      for (int c=0;c<diag;c++){
        int x=(ix+xfact*count2);
        while (x<0) {x+=xdim_;count=0;}
        if (x>=(int)xdim_) {x=x%xdim_;count=0;}
        int y=(iy+yfact*count2);
        while (y<0) {y+=ydim_;count=0;}
        if (y>=(int)ydim_) {y=y%ydim_;count=0;}
        int z=(zstart+zfact*count2);
        while (z<0) {z+=zdim_;count=0;}
        if (z>=(int)zdim_) {z=z%zdim_;count=0;}
        count2++;
        if (x==xend||y==yend||z==zend) surf=false;
        if (surf){
          if (psp){
            if (grid_[x][y][z]==EMPTY){
              count++;
            }else if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET||grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              if (count<= thr && pcount>0){
                if (markpsp_){
                  for (int c=count;c>0; c--){
                    if (t_surf) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=TP_POCKET;
                    else if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]!=TP_POCKET) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=POCKET;
                  }
                }else{
                  for (int c=count;c>0; c--)
                    if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]==EMPTY) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=PSP;

                }
              }
              start=c;
              count=0;
              if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE) pcount=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }

          }else{
            if (grid_[x][y][z]==EMPTY||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              if (count> thr){
                for (int c=count;c>0; c--){
                  if (t_surf) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=TP_POCKET;
                  else if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]!=TP_POCKET) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=POCKET;
                }
              }
              start=c;
              count=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }
          }
        }else{
          if (psp){
            if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET||grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              pcount=0;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              surf=true;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              count=0;
              start=c;
            }
          }else{
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              surf=true;
              count=0;
              start=c;
            }
          }
        }
      }
    }}

}

void PocketGrid::searchD3(core::Size thr, bool psp){
  int diag=(int)floor(sqrt(pow(zdim_,2)+pow(zdim_,2)+pow(zdim_,2)));
  int zstart=0;
  int xfact=-1,yfact=-1,zfact=1;
  int xend=0,yend=0,zend=0;
  if (xfact<0) xend=xdim_-1;
  if (yfact<0) yend=ydim_-1;
  if (zfact<0) zend=zdim_-1;
  for (core::Size ix=0;ix<xdim_;ix++){
    for (core::Size iy=0;iy<ydim_;iy++){
      bool surf=false;
      bool t_surf=false;
      core::Size count=0;
      core::Size pcount=0;
      core::Size count2=0;
      int start=-1;
      for (int c=0;c<diag;c++){
        int x=(ix+xfact*count2);
        while (x<0) {x+=xdim_;count=0;}
        if (x>=(int)xdim_) {x=x%xdim_;count=0;}
        int y=(iy+yfact*count2);
        while (y<0) {y+=ydim_;count=0;}
        if (y>=(int)ydim_) {y=y%ydim_;count=0;}
        int z=(zstart+zfact*count2);
        while (z<0) {z+=zdim_;count=0;}
        if (z>=(int)zdim_) {z=z%zdim_;count=0;}
        count2++;
        if (x==xend||y==yend||z==zend) surf=false;
        if (surf){
          if (psp){
            if (grid_[x][y][z]==EMPTY){
              count++;
            }else if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET||grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              if (count<= thr && pcount>0){
                if (markpsp_){
                  for (int c=count;c>0; c--){
                    if (t_surf) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=TP_POCKET;
                    else if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]!=TP_POCKET) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=POCKET;
                  }
                }else{
                  for (int c=count;c>0; c--)
                    if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]==EMPTY) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=PSP;
                }
              }
              count=0;
              pcount=0;
              start=c;
              count=0;
              if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE) pcount=0;
            } else{
              surf=false;
              t_surf=false;
            }

          }else{
            if (grid_[x][y][z]==EMPTY||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              if (count> thr){
                for (int c=count;c>0; c--){
                  if (t_surf) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=TP_POCKET;
                  else if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]!=TP_POCKET) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=POCKET;
                }
              }
              count=0;
              pcount=0;
              start=c;
              count=0;
            } else{
              surf=false;
              t_surf=false;
            }
          }
        }else{
          if (psp){
            if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET||grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              pcount=0;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              surf=true;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              count=0;
              start=c;
            }
          }else{
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              surf=true;
              count=0;
              start=c;
            }
          }
        }
      }
    }}

}

void PocketGrid::searchD4(core::Size thr, bool psp){
  int diag=(int)floor(sqrt(pow(zdim_,2)+pow(zdim_,2)+pow(zdim_,2)));
  int zstart=0;
  int xfact=1,yfact=-1,zfact=1;
  int xend=0,yend=0,zend=0;
  if (xfact<0) xend=xdim_-1;
  if (yfact<0) yend=ydim_-1;
  if (zfact<0) zend=zdim_-1;
  for (core::Size ix=0;ix<xdim_;ix++){
    for (core::Size iy=0;iy<ydim_;iy++){
      bool surf=false;
      bool t_surf=false;
      core::Size count=0;
      core::Size pcount=0;
      core::Size count2=0;
      int start=-1;
      for (int c=0;c<diag;c++){
        int x=(ix+xfact*count2);
        while (x<0) {x+=xdim_;count=0;}
        if (x>=(int)xdim_) {x=x%xdim_;count=0;}
        int y=(iy+yfact*count2);
        while (y<0) {y+=ydim_;count=0;}
        if (y>=(int)ydim_) {y=y%ydim_;count=0;}
        int z=(zstart+zfact*count2);
        while (z<0) {z+=zdim_;count=0;}
        if (z>=(int)zdim_) {z=z%zdim_;count=0;}
        count2++;
        if (x==xend||y==yend||z==zend) surf=false;
        if (surf){
          if (psp){
            if (grid_[x][y][z]==EMPTY){
              count++;
            }else if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET||grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE||grid_[x][y][z]==TP_POCKET) t_surf=true;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              if (count<= thr && pcount>0){
                if (markpsp_){
                  for (int c=count;c>0; c--){
                    if (t_surf) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=TP_POCKET;
                    else if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]!=TP_POCKET) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=POCKET;
                  }
                }else{
                  for (int c=count;c>0; c--)
                    if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]==EMPTY) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=PSP;

                }
              }
              start=c;
              count=0;
              if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE || grid_[x][y][z]==T_SURFACE) pcount=0;
            } else{
              count=0;
              pcount=0;
              surf=false;
              t_surf=false;
            }

          }else{
            if (grid_[x][y][z]==EMPTY||grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET){
              count++;
            }else if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              if (count> thr){
                for (int c=count;c>0; c--){
                  if (t_surf) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=TP_POCKET;
                  else if (grid_[x-xfact*c][y-yfact*c][z-zfact*c]!=TP_POCKET) grid_[x-xfact*c][y-yfact*c][z-zfact*c]=POCKET;
                }
              }
              count=0;
              pcount=0;
              start=c;
              count=0;
            } else{
              surf=false;
              t_surf=false;
            }
          }
        }else{
          if (psp){
            if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET||grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              pcount=0;
              if (grid_[x][y][z]==POCKET||grid_[x][y][z]==TP_POCKET) pcount++;
              surf=true;
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              count=0;
              start=c;
            }
          }else{
            if (grid_[x][y][z]==HSURFACE || grid_[x][y][z]==PSURFACE||grid_[x][y][z]==T_SURFACE){
              if (grid_[x][y][z]==T_SURFACE) t_surf=true;
              surf=true;
              count=0;
              start=c;
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
  std::filebuf fb;
  std::stringstream filename;
  filename<<tag_<<"pocket"<<pdbno_<<".pdb";
  pdbno_++;
  fb.open (filename.str().c_str(),std::ios::out);
  std::ostream os(&fb);
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
        os<<concatenated_pdb_info;
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
      os<<concatenated_pdb_info;

    }
    clustNo++;
  }
  fb.close();
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
        clusters_.add(x,y,z, stepSize_);
        return;
      }
      else st=x-ssteps;
      if (x+ssteps>=xdim_) {
        en=xdim_;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        clusters_.add(x,y,z, stepSize_);
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
        clusters_.add(x,y,z, stepSize_);
        return ;
      }
      else st=y-ssteps;
      if (y+ssteps>=ydim_) {
        en=ydim_;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        clusters_.add(x,y,z, stepSize_);
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
        clusters_.add(x,y,z, stepSize_);
        return ;
      }
      else st=z-ssteps;
      if (z+ssteps>=zdim_) {
        en=zdim_;
        if (grid_[x][y][z]==POCKET){
          grid_[x][y][z]=PO_EDGE;
        }else grid_[x][y][z]=TP_EDGE;
        clusters_.add(x,y,z, stepSize_);
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
          clusters_.add(x,y,z, stepSize_);
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
          clusters_.add(x,y,z, stepSize_);
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
          clusters_.add(x,y,z, stepSize_);
          return ;
        }
      }
      clusters_.add(x,y,z, stepSize_);
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
      clusters_.findClusters();

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
          for(Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i) {
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



EggshellGrid::EggshellGrid( const PocketGrid& gr ) :
	// initialize member data using PocketGrid
	PocketGrid( gr )
{

	init();
	eggshell_coord_list_.clear();
	extra_coord_list_.clear();

	pocket_CoM_.zero();
	core::Size divcountr = 0;

	// use PocketGrid to fill in Eggshell grid
	numeric::xyzVector<core::Real> grid_coord;
	core::Size searchxmin, searchxmax, searchymin, searchymax, searchzmin, searchzmax, xx,yy,zz;
	for (std::list<PCluster>::const_iterator cit=gr.clusters_.clusters_.begin(); cit != gr.clusters_.clusters_.end(); ++cit){
		if (cit->points_.size()*pow(gr.stepSize_,3)<gr.minPockSize_) continue; // this is a smallpocket
		for (std::list<PCluster::Cxyz>::const_iterator pit=cit->points_.begin(); pit != cit->points_.end(); ++pit){
			if ( (gr.grid_[pit->x][pit->y][pit->z]==gr.TP_POCKET) || (gr.grid_[pit->x][pit->y][pit->z]==gr.TP_SURF) || (gr.grid_[pit->x][pit->y][pit->z]==gr.TP_EDGE) || (gr.grid_[pit->x][pit->y][pit->z]==gr.TP_BURIED) ) {
				pocket_CoM_.x() += pit->x*gr.stepSize_+gr.xcorn_;
				pocket_CoM_.y() += pit->y*gr.stepSize_+gr.ycorn_;
				pocket_CoM_.z() += pit->z*gr.stepSize_+gr.zcorn_;
				divcountr++;

				if (gr.grid_[pit->x][pit->y][pit->z]==gr.TP_BURIED){
					if (pit->x != 0) {
						searchxmin=pit->x-1;
					}else{
						searchxmin=pit->x;
					}
					if (pit->x != gr.xdim_-1) {
						searchxmax=pit->x+1;
					}else{
						searchxmax=pit->x;
					}
					if (pit->y != 0) {
						searchymin=pit->y-1;
					}else{
						searchymin=pit->y;
					}
					if (pit->y != gr.ydim_-1) {
						searchymax=pit->y+1;
					}else{
						searchymax=pit->y;
					}
					if (pit->z != 0) {
						searchzmin=pit->z-1;
					}else{
						searchzmin=pit->z;
					}
					if (pit->z != gr.zdim_-1) {
						searchzmax=pit->z+1;
					}else{
						searchzmax=pit->z;
					}
					bool found=false;
					for (xx = searchxmin; ( xx<=searchxmax) && ! found; xx++){
						for (yy = searchymin; ( yy<=searchymax) && ! found; yy++){
							for (zz = searchzmin; ( zz<=searchzmax) && ! found; zz++){
								if (xx==pit->x && yy==pit->y && zz==pit->z) continue;
								if ( (gr.grid_[xx][yy][zz] == gr.HSURFACE) || (gr.grid_[xx][yy][zz] == gr.PSURFACE) || (gr.grid_[xx][yy][zz] == gr.T_SURFACE)) {
									found = true;
									grid_[pit->x][pit->y][pit->z] = EGGSHELL;
									grid_coord.x() = (pit->x*gr.stepSize_+gr.xcorn_);
									grid_coord.y() = (pit->y*gr.stepSize_+gr.ycorn_);
									grid_coord.z() = (pit->z*gr.stepSize_+gr.zcorn_);
									eggshell_coord_list_.push_back(grid_coord);
								}
							}
						}
					}
				}
			}
		}
	}
	pocket_CoM_ /= divcountr;

	//End eggshell

	// generate extra points around eggshell
	numeric::xyzVector<core::Real> extra_coord;
	for (core::Size x=0;x<(gr.xdim_); x++){
    for (core::Size y=0;y<(gr.ydim_); y++){
      for (core::Size z=0;z<(gr.zdim_); z++){
				if ( ( gr.grid_[x][y][z]==gr.EMPTY ) || ( gr.grid_[x][y][z]==gr.POCKET ) || ( gr.grid_[x][y][z]==gr.PO_SURF ) || ( gr.grid_[x][y][z]==gr.PO_BURIED ) ) {

					/*Print grid points
					std::string output_gridfile = "grid_points.pdb";
					std::filebuf f5;
					f5.open (output_gridfile.c_str(), std::ios::out);
					std::ostream os5(&f5);
					std::string f5_info;
					std::stringstream f5_tmp;
					f5_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x*gr.stepSize_+gr.xcorn_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y*gr.stepSize_+gr.ycorn_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z*gr.stepSize_+gr.zcorn_<<std::endl;
					f5_info += f5_tmp.str();
					f5_tmp.str(std::string());
					os5<<f5_info;
					f5.close();
					//End Printing grid points
					*/

					// check if this point touches PSURFACE or HSURFACE or T_SURFACE
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

					bool touches_surface=false;
					for (xx = searchxmin; xx<=searchxmax; xx++){
						for (yy = searchymin; yy<=searchymax; yy++){
							for (zz = searchzmin; zz<=searchzmax; zz++){
								if (xx==x && yy==y && zz==z) continue;
								if ( (gr.grid_[xx][yy][zz] == gr.HSURFACE) || (gr.grid_[xx][yy][zz] == gr.PSURFACE) || (gr.grid_[xx][yy][zz] == gr.T_SURFACE)) {
									touches_surface = true;
									break;
								}
							}
							if ( touches_surface ) break;
						}
						if ( touches_surface ) break;
					}

					if ( touches_surface ) {
						// check if this point is within 5 A of any eggshell point
						extra_coord.x() = (x*gr.stepSize_+gr.xcorn_);
						extra_coord.y() = (y*gr.stepSize_+gr.ycorn_);
						extra_coord.z() = (z*gr.stepSize_+gr.zcorn_);
						for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = eggshell_coord_list_.begin(); pd != eggshell_coord_list_.end(); ++pd) {
							//compute distance
							if ( pd->distance( extra_coord ) < 5.00001 ) {
								extra_coord_list_.push_back(extra_coord);
								break;
							}
						}
					}

				}
			}
		}
	}

	//Print Extra points to PDB file
	std::string output_pdbname = "eggshell.pdb";
	std::filebuf f4;
	f4.open (output_pdbname.c_str(), std::ios::out);
	std::ostream os4(&f4);
	std::string f4_info;
	std::stringstream f4_tmp;
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = eggshell_coord_list_.begin(); pd != eggshell_coord_list_.end(); ++pd) {
		f4_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   EGG A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->z()<<std::endl;
		f4_info += f4_tmp.str();
		f4_tmp.str(std::string());
	}
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = extra_coord_list_.begin(); pd != extra_coord_list_.end(); ++pd) {
		f4_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   EXT A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->z()<<std::endl;
		f4_info += f4_tmp.str();
		f4_tmp.str(std::string());
	}

	return;

}




    } // pockets
  } // protocols
