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

#ifndef INCLUDED_core_scoring_fiber_diffraction_FiberDiffraction_hh
#define INCLUDED_core_scoring_fiber_diffraction_FiberDiffraction_hh

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/scoring/fiber_diffraction/util.hh>
#include <core/scoring/fiber_diffraction/xray_scattering.hh>
#include <core/scoring/fiber_diffraction/FiberDiffraction.fwd.hh>

//C++ headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace fiber_diffraction {

	class FiberDiffraction : public utility::pointer::ReferenceCount {
		public:
		 /// @brief constructor
    FiberDiffraction();

		///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
		virtual ~FiberDiffraction();

		/// @brief Initialize map from cmd line options
		void init();

		/// @brief Load fiber diffraction layer lines
		bool loadFiberDiffractionData( 
			std::string layer_lines,
			core::Real const & c, 
			core::Real const & res_cutoff_high, 
			core::Real const & res_cutoff_low
		);

		// @brief Load fiber diffraction layer lines
		bool loadFiberDiffractionData( 
			std::istream & input, 
			std::string layer_lines, 
			core::Real const & c, 
      core::Real const & res_cutoff_high, 
      core::Real const & res_cutoff_low	
		);
		
		/// @brief Load fiber diffraction layer lines
		bool trimFiberDiffractionData( 
			core::Real const &c,
  		core::Real const &res_cutoff_high,
  		core::Real const &res_cutoff_low 
		);

		bool setupBesselOrder();		

		inline bool isFiberDataLoaded() const { return this->isLoaded; };

		inline void getAllFiberData( 
			utility::vector0< utility::vector1< core::Real > >::iterator & layer_lines_I_it,
			utility::vector0< utility::vector1< core::Real > >::iterator & layer_lines_R_it,
			utility::vector0 < utility::vector0 < int > >::iterator  & nvals_it,
			core::Size & lmax_,
			core::Size & Rmax_
		) {
      layer_lines_I_it = layer_lines_I.begin();
			layer_lines_R_it = layer_lines_R.begin();
			nvals_it = nvals.begin();
			lmax_ = lmax;
			Rmax_ = Rmax;
    }
			
		inline void getIntensities( utility::vector0< utility::vector1< core::Real > >::iterator & layer_lines_I_it) {
			layer_lines_I_it = layer_lines_I.begin();
		}
	
		inline void getReciprocalRs( utility::vector0< utility::vector1< core::Real > >::iterator & layer_lines_R_it ) {
      layer_lines_R_it = layer_lines_R.begin();
    }

		inline void getNVals ( utility::vector0 < utility::vector0 < int > >::iterator  & nvals_it  ) {
			nvals_it = nvals.begin();
		}
		
		
		///////////
		private:
		bool isLoaded;
		core::Real c;
    core::Real res_cutoff_high;
    core::Real res_cutoff_low;
	
		//Layer lines, maximum layer line number and reciprocal R
		core::Size lmax, Rmax;
		utility::vector0 < utility::vector1< core::Real > > original_layer_lines_I;
 		utility::vector0 < utility::vector1< core::Real > > original_layer_lines_R;

		utility::vector0< utility::vector1< core::Real > > layer_lines_I;
  	utility::vector0< utility::vector1< core::Real > > layer_lines_R;

		//Bessel orders
		utility::vector0 < utility::vector0 < int > > nvals;
		
		//TODO: What to do with form factors
		//utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator form_factors;
		//Scattering centroids
		//utility::vector1< OneGaussianScattering >::iterator sig_centroid;	
	};

	/// @brief The EDM instance
	FiberDiffraction& getFiberDiffractionData(core::Real c = 0.0, core::Real res_cutoff_high = 0.0, core::Real res_cutoff_low = 0.0, bool force_reload = false);
	/// @brief The EDM instance
	FiberDiffraction& getFiberDiffractionData_legacy(core::Real c = 0.0, core::Real res_cutoff_high = 0.0, core::Real res_cutoff_low = 0.0, bool force_reload = false);
}
}
}
#endif
