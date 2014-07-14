// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/cryst/PhenixInterface.hh
/// @brief  Singleton class that manages interface with phenix refinement, incl. temporary storage
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_cryst_PhenixInterface_hh
#define INCLUDED_core_scoring_cryst_PhenixInterface_hh

#ifdef WITH_PYTHON
#include "Python.h"
#endif

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility headers
#include <utility/exit.hh>


// C++ headers
#include <string>
#include <map>
#include <complex>
namespace core {
namespace scoring {
namespace cryst {

class PhenixInterface {
public:
	///@brief constructor
	PhenixInterface();

	///@brief score a structure
	core::Real getScore (core::pose::Pose const & pose);

	///@brief score a structure with derivatives
	core::Real getScoreAndDerivs (
		core::pose::Pose const & pose,
		utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > & grads);

	///@brief fit bfactors
	void fitBfactors (core::pose::Pose & pose);

	///@brief dump r and rfree to a string for informational purposes
	std::string getInfoLine();

	///@brief dump r and rfree
	core::Real getR();
	core::Real getRfree();

	///@brief update fcalc
	void updateFcalc ();

	///@brief update mask
	void updateSolventMask ();

	///@brief update mask
	void updateSolventMask (core::pose::Pose const & pose);

	///@brief optimize fmask
	void optimizeSolventMask ();

	///@brief explicitly recompute ksol/bsol
	void optimizeSolvParams ();

	///@brief explicitly recompute ksol/bsol and fmask
	void optimizeSolvParamsAndMask ();

	///@brief set the res limits
	void setResLimits(core::Real res_high=0.0, core::Real res_low=0.0);

	///@brief set twin law
	void setTwinLaw(std::string twin_law);

	///@brief set sf calculation algorithm
	void setAlgorithm(std::string twin_law);

	///@brief set target function
	void set_map_type ( std::string map_type );

	///@brief set strategy for adp refinement
	void set_adp_strategy ( std::string adp_strat ) { adp_strategy_ = adp_strat; }

	///@brief set target function
	void set_target_function ( std::string tgt_val ) { target_function_ = tgt_val; }

	///@brief set target function
	void set_cif_files ( utility::vector1<std::string> cif_in ) { cif_files_ = cif_in; }

	///@brief use pose to rephase data; calculate a new density map; return map file name
	std::string calculateDensityMap (core::pose::Pose & pose, bool no_sidechain=false);

private:
	// helper function
	void stealBfactorsFromFile(core::pose::Pose & pose, std::string filename);

	// called once to initialize the python evaluator object (needs the pose so it can't be put in the constructor)
	void initialize_target_evaluator( core::pose::Pose const & pose, std::string eff_file="" );

#ifdef WITH_PYTHON
	// helper function: convert a pose to a python list
	PyObject* pose_to_pycoords( core::pose::Pose const & pose );

	// helper function: convert a python list to a vector1< xyzVector< Real > >
	void pylist_to_grads(
		core::pose::Pose const & pose,
		PyObject* pygrads,
		utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > & rosgrads );
#endif

#ifdef WITH_PYTHON
	// These two variables are only used if this code is compiled with the WITH_PYTHON flag set.
	// The ifdef guard around their inclusion here is to avoid unused-private-variable warnings
	// from clang when that flag is not set.
	core::Real res_low_, res_high_;
#endif
	std::string tempdir_;
	std::string mtzfile_;
	std::string phenix_home_;
	std::string adp_strategy_;
	std::string target_function_;
	std::string twin_law_;
	std::string algo_;
	std::string map_type_;
	utility::vector1<std::string> cif_files_;

#ifdef WITH_PYTHON
	PyObject *target_evaluator_;
#endif

	// reference pose used to initialize the evaluator
	core::pose::PoseOP ref_pose_;
};

/// @brief The EDM instance
PhenixInterface& getPhenixInterface();


}
}
}


#endif

