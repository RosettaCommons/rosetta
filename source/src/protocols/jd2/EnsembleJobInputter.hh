// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd2/EnsembleJobInputter.hh
/// @brief A Job Inputter for distributing a job based on a set of input structures that make up an ensemble
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_jd2_EnsembleJobInputter_hh
#define INCLUDED_protocols_jd2_EnsembleJobInputter_hh

#include <protocols/jd2/EnsembleJobInputter.fwd.hh>

//unit headers
//#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/PDBJobInputter.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobsContainer.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace jd2 {




///@brief A Job Inputter for distributing a job based on a set of input structures that make up an ensemble
///  Two modes.  Seed Ensemble and Grid Ensemble.
///
/// @details
///
/// Seed Ensemble
/// =============
///
/// Randomly choose the starting files using the list of structures given by -s and -l
///
/// -jd2:seed_ensemble_weights (RealVector)
///
///    - Will give weights to the input PDBs and trigger seed_ensemble mode and the JI.  We then use the weighted sampler to choose the input pdb for
///    each nstruct.  The weights could be the size of the cluster, the energy, that you like some structure better, etc. Must match number of inputs.
///
/// -jd2:seed_ensemble_weight_file (File)
///
///    - A file specifying weights to use for each input structure.  Enables seed_ensemble mode and the JI. Two columns.  basename with extension (or
///    relative path or full path), weight.  Any not given in file will be set to 0 by default.  Can give a line that is [ALL weight] to set all input
///    pdbs to a given weight. Used for example, to upweight a specific structure:
///
///    Example:
///
///    #name weight
///    ALL 1
///    awesome_model.pdb 3
///
///
///
/// -jd2:seed_ensemble (Bool)
///
///    - Enable seed ensemble mode, but simply randomly choose the input pdb for each nstruct.  For seed ensemble mode, the number of input pdbs can
///    be larger than nstruct.
///
///
/// Grid Ensemble
/// =============
/// Use the input files given in -s and -l and nstruct to cover a grid.
///
/// -in:jd2:grid_ensemble (Bool)
///
///   - Will enable the basic component of the JI.  Here, instead of sampling nstruct for every input pdb, we only sample nstruct no matter the number
///   of input PDBs (with nstruct split as evenly as possible over the input PDBs).
///
///
class EnsembleJobInputter : public protocols::jd2::PDBJobInputter {

public:

	EnsembleJobInputter();

	~EnsembleJobInputter() override;


	/// @brief this function determines what jobs exist from -s/-l
	void fill_jobs( JobsContainer & jobs ) override;

private:

	utility::vector1< core::Real >
	read_get_weights(std::string const & filename, utility::vector1< std::string > const & inputs ) const ;

};


} //protocols
} //jd2



#endif //INCLUDED_protocols_jd2_EnsembleJobInputter_hh





