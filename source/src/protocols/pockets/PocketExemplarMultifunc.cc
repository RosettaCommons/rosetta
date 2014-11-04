// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/PocketExemplarMultifunc.hh
/// @brief  Pocket multifunction class
/// @author David Johnson

/// Unit headers
#include <protocols/pockets/PocketExemplarMultifunc.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <cmath>

#include <utility/vector1.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>

using namespace core;
using namespace core::pose;
using namespace basic::options;
using namespace std;
using namespace basic::options::OptionKeys;

namespace protocols{
namespace pockets {

PocketExemplarMultifunc::PocketExemplarMultifunc(std::string const input_pdb_name, std::string const resid, core::Real const c_rad, core::Real const rep_weight, utility::vector1<core::Real> & p_min, utility::vector1<core::Real> & p_max) {
  core::import_pose::pose_from_pdb( input_pose, input_pdb_name );
  residues = protocols::pockets::PocketGrid::getRelaxResidues(input_pose, resid);
  cRad = c_rad;
  repW = rep_weight; 
  vdWpen = option[ OptionKeys::pocket_grid::pocket_exemplar_vdw_pen ]();
  optimal = cRad;
  pg = protocols::pockets::PocketGrid( residues );  
  pg.autoexpanding_pocket_eval( residues, input_pose ) ;
  pg.dumpTargetPocketsToPDB("testpock.pdb");
  utility::vector1<core::Real> bounds = pg.getBounds();
  for (int i = 1; i < (int)p_min.size();i+=3){
	  for (int j=0;j<3;j++){
		p_min[i+j]=bounds[1+2*j];
		p_max[i+j]=bounds[2+2*j];
	  }
  }
}

//This is the objective function
core::Real
PocketExemplarMultifunc::operator ()( core::optimization::Multivec const & vars ) const {
  core::Real score=0;
  ComparisonGrid cg(pg);
  for (int i = 1; i < (int)vars.size();i+=3){
    score += cg.mark(pg, vars[i], vars[i+1], vars[i+2], cRad, repW);
  }
  score += cg.compareCoverage(pg);

  //hard spheres model
  for (int i = 1; i < (int)vars.size();i+=3){
    bool near=false;
    for (int j = 1; j < (int)vars.size();j+=3){
      if (i==j) continue;
      core::Real dist = std::sqrt(std::pow(vars[i]-vars[j],2) + std::pow(vars[i+1]-vars[j+1],2) + std::pow(vars[i+2]-vars[j+2],2));
      if ((dist < optimal)){
        score += repW * 8 * (3.14159/(12*dist))*std::pow(cRad-dist,2)*(pow(dist,2)+2*dist*cRad);
      }else if (dist<4) near=true;
    }
    if (!near) score += 300;
  }
	return score;
}

//This is not used in my project, is stubbd out because this can be used for other optimizations that use derivatives
void
PocketExemplarMultifunc::dfunc( core::optimization::Multivec const &, core::optimization::Multivec &) const
{
  std::cerr << "PocketExemplarMultifunc cannot be used with derivative based optimization methods."<<std::endl;
  exit(20);
}

/// @details Useful debugging code that can be re-enabled by changing the boolean
/// variables at the top of this file.
void
PocketExemplarMultifunc::dump( core::optimization::Multivec const &, core::optimization::Multivec const &) const {
  std::cout<< "In PocketExemplarMultifunc."<< std::endl;
}

} // namespace pockets
} // namespace protocols

