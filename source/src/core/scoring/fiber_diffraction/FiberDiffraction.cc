// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffraction.hh
/// @brief  Fiber diffraction data collector class
/// @author Wojciech Potrzebowski and Ingemar Andre

// Unit Headers
#include <core/scoring/fiber_diffraction/FiberDiffraction.hh>
#include <core/scoring/fiber_diffraction/util.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

// C++ headers
#include <fstream>


namespace core {
namespace scoring {
namespace fiber_diffraction {

using basic::T;
using basic::Tracer;
static thread_local basic::Tracer TR( "core.scoring.fiber_diffraction.FiberDiffraction" );

using namespace core;
using namespace basic::options;

/// null constructor
FiberDiffraction::FiberDiffraction() {
  init();
}

void FiberDiffraction::init() {
  isLoaded = false;
	c =0.0;
	res_cutoff_high = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_high ]();
	res_cutoff_low = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_low ]();
}

FiberDiffraction::~FiberDiffraction() {}

FiberDiffraction& getFiberDiffractionData( core::Real c, core::Real res_cutoff_high, core::Real res_cutoff_low,  bool force_reload ) {
	if(basic::resource_manager::ResourceManager::get_instance()->
		has_resource_with_description("FiberDiffractionLayerLines")){
		FiberDiffractionOP fiber_diffraction(
			basic::resource_manager::get_resource< FiberDiffraction >(
				"FiberDiffractionLayerLines"));
		return *fiber_diffraction;
	} else {
		return getFiberDiffractionData_legacy(c, res_cutoff_high, res_cutoff_low, force_reload);
	}
}

//
FiberDiffraction& getFiberDiffractionData_legacy(core::Real c, core::Real res_cutoff_high, core::Real res_cutoff_low, bool force_reload ) {
	static FiberDiffraction theFiberDiffractionData;

	if ( !theFiberDiffractionData.isFiberDataLoaded() || force_reload ) {
		if (!basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::layer_lines ].user() ) {
			 utility_exit_with_message("[ ERROR ] No fiber diffraction layer lines specified!");
		} else {
			bool fiber_data_loaded=false;
			bool bessel_initialized=false;
			std::string layer_lines;
			layer_lines = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::layer_lines ]();

			TR << "Loading fiber diffraction layer lines" << std::endl;
			fiber_data_loaded = theFiberDiffractionData.loadFiberDiffractionData( layer_lines, c, res_cutoff_high, res_cutoff_low );
			if (!fiber_data_loaded) {
				TR.Error << "[ ERROR ] Error loading layer_lines named: '" << layer_lines  << "'" << std::endl;
				exit(1);
			}

			//Initiallizing Bessel orders
			bessel_initialized = theFiberDiffractionData.setupBesselOrder();
			if (!bessel_initialized) {	
        TR.Error << "[ ERROR ] Error initializing bessel orders. '"<< std::endl;
        exit(1);
      }

		}
	}
	return theFiberDiffractionData;
}


bool FiberDiffraction::trimFiberDiffractionData(
	core::Real const &c, 
	core::Real const &res_cutoff_high, 
	core::Real const &res_cutoff_low 
) {
	isLoaded= false;
	utility::vector0 < utility::vector1< core::Real > > layer_lines_I_trim;
  utility::vector0 < utility::vector1< core::Real > > layer_lines_R_trim;

	core::Size Rmax_line;
  core::Size Rexc(0);
  core::Real Rinv;
  core::Size lmax_trim(lmax);
  core::Size Rmax_trim(0);
	core::Real res_cutoff_high2(res_cutoff_high*res_cutoff_high);
  core::Real res_cutoff_low2(res_cutoff_low*res_cutoff_low);
  if (lmax_trim>ceil(c*res_cutoff_high)) lmax_trim=ceil(c*res_cutoff_high);
  layer_lines_I_trim.resize(lmax_trim+1);
  layer_lines_R_trim.resize(lmax_trim+1);

  for (core::Size l=0; l<=lmax_trim; ++l) {
      Rmax_line = 0;
			core::Size Rsize(original_layer_lines_R[l].size());
      for (core::Size R=1; R<=Rsize; ++R) {
        Rinv = original_layer_lines_R[l][R];
        if ((Rinv*Rinv+(l/c)*(l/c)) <= res_cutoff_high2 && (Rinv*Rinv+(l/c)*(l/c)) >= res_cutoff_low2 && Rinv>0.0) {
          layer_lines_R_trim[l].push_back(Rinv);
          layer_lines_I_trim[l].push_back(original_layer_lines_I[l][R]);
          Rmax_line++;
        } else {
          Rexc++;
        }
      }
      if (Rmax_line>Rmax_trim) Rmax_trim = Rmax_line;
  }
  TR<<"Number of excluded Rs after threshold selection: "<<Rexc<<std::endl;

	lmax = lmax_trim;
  Rmax = Rmax_trim;
  layer_lines_I = layer_lines_I_trim;
  layer_lines_R = layer_lines_R_trim;
	
	isLoaded=true;
	return isLoaded;
}

bool FiberDiffraction::loadFiberDiffractionData(
	std::string layer_lines,
	core::Real const & c,
	core::Real const & res_cutoff_high,
	core::Real const & res_cutoff_low
) {
		std::ifstream input( layer_lines.c_str() );
		bool isLoaded( loadFiberDiffractionData(input,layer_lines, c, res_cutoff_high, res_cutoff_low) );
		return isLoaded;
}

bool FiberDiffraction::loadFiberDiffractionData(
	std::istream & input, 
	std::string layer_lines,
	core::Real const & c,
	core::Real const & res_cutoff_high,
	core::Real const & res_cutoff_low ) {

  core::Size lmax_, Rmax_;
  original_layer_lines_I.resize(80);
  original_layer_lines_R.resize(80);

	if (!input) {
    TR << "[ ERROR ]  Error opening layer lines " << layer_lines << std::endl;
    return false;
  }

  // load layer lines from disk
 // utility::io::izstream input(layer_lines.c_str());
  std::string line;

  Size nlines(0);
  core::Real R,I;
  core::Size l(0), l_prev (0);
  Rmax_ = 0;
  lmax_ = 0;
  core::Size nlines_max(0);
  while( getline( input, line ) ) {
    if ( line.substr(0,1) == "#" ) continue;
    if ( line.length() == 0 ) continue;
    std::istringstream line_stream( line );
    ++nlines;
    if ( nlines > 400 ) {
      utility_exit_with_message("[ ERROR ] No more than 400 data points per layer line can be loaded!");
    }

    line_stream >> I >> R >> l;
    TR.Debug << "    read " << l << " : " << R << " , " << I  << std::endl;

    if ( l > 80 ) {
      utility_exit_with_message("[ ERROR ] No more than 80 layer lines can be loaded!");
    }

    original_layer_lines_R[l].push_back(R); //Resolution data....
    original_layer_lines_I[l].push_back(I) ; //Amplitude data....

    if ( l_prev == (Size)-1 )
      l_prev = l;
    if ( l_prev != l ) {
      l_prev = l;
      nlines = 1;
    }
    if ( l > lmax_ ) lmax_ = l;
    if ( nlines > nlines_max ) {
      nlines_max=nlines;
		}
  }
  Rmax_ = nlines_max;
	
	//Writing to class values
	lmax = lmax_;
  Rmax = Rmax_;
  layer_lines_I = original_layer_lines_I;
  layer_lines_R = original_layer_lines_R;

	if ( c>0 ) {
		if ( res_cutoff_high==0.0 ) utility_exit_with_message("[ ERROR ] High-resolution cutoff has to be set");
  	if ( res_cutoff_low==0.0 ) utility_exit_with_message("[ ERROR ] Low-resolution cutoff has to be set");
		TR << "Adjusting layer lines to given resolution cutoff" << std::endl;
		isLoaded = trimFiberDiffractionData(c, res_cutoff_high, res_cutoff_low);
	} 
	else {
		isLoaded=true;
	}

	return isLoaded;

}

bool FiberDiffraction::setupBesselOrder()
{
  TR << "Calculating Bessel orders..." << std::endl;
	isLoaded=false;
  core::Size cn_symmetry(0);
	core::Real a(0.0), b(0.0);
	core::Size n_max(50);
	int m_max(30);

	nvals.resize(lmax+1);

	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ].user()) {
    a = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ]();
  }

  if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ].user()) {
    b = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ]();
  }

	if ( !(a > 0 ) ){
    utility_exit_with_message("The number of subunits per repeat, score::fiber_diffraction::a, must be set!");
  }

  if ( !(b > 0 ) ){
    utility_exit_with_message("The number of turns per repeat, score::fiber_diffraction::b, must be set!");
  }
	
  if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::max_bessel_order ].user()) {
    n_max = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::max_bessel_order ]();
  }
  if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::cn_symmetry ].user()) {
    cn_symmetry = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::cn_symmetry ]();
  }

	TR<<"Maximum Bessel order taken into consideration is: "<<n_max<<std::endl;
  for (Size l=0; l<= lmax; ++l ) {
    for (int m=-m_max; m<= m_max; ++m) {
      Real fn=(l-a*m)/b;
      int n(fn);
      core::Size nabs =  abs(n);
      if(std::fmod(fn, 1.0) != 0  || nabs > n_max ) continue;
      if (cn_symmetry) {
        if(std::fmod(fn, cn_symmetry) == 0 ) nvals[l].push_back(n);
      }
      else nvals[l].push_back(n);
    }
  }
	isLoaded=true;
  return isLoaded;
}


}//namespace
}//namespace
}//namespace
