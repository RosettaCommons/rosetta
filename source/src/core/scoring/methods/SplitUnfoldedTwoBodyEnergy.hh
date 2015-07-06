// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SplitUnfoldedTwoBodyEnergy.hh
/// @brief  Header for the split unfolded energy two body component energy method.
/// @author Riley Simmons-Edler (rse231@nyu.edu)


#ifndef INCLUDED_core_scoring_methods_SplitUnfoldedTwoBodyEnergy_hh
#define INCLUDED_core_scoring_methods_SplitUnfoldedTwoBodyEnergy_hh

#include <core/scoring/methods/SplitUnfoldedTwoBodyEnergy.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/SplitUnfoldedTwoBodyPotential.hh>

#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace methods {


  class SplitUnfoldedTwoBodyEnergy : public ContextIndependentOneBodyEnergy
  {

  public:
    typedef ContextIndependentOneBodyEnergy parent;

    //instantiate using the weights from the data file
    SplitUnfoldedTwoBodyEnergy(std::string const & label_type,std::string const & value_type, std::string const & score_func_type);
    //instantiate using given weight emap
    SplitUnfoldedTwoBodyEnergy(std::string const & label_type,std::string const & value_type, std::string const & score_func_type, const EnergyMap & emap_in);
		~SplitUnfoldedTwoBodyEnergy();

    virtual EnergyMethodOP clone() const;

    virtual void residue_energy(conformation::Residue const & rsd,pose::Pose const &,EnergyMap & emap) const;

		virtual	void indicate_required_context_graphs(utility::vector1<bool> &) const;

  private:
    std::string label_type_; //the atom type set for this two body energy, e.g. rosetta types, mm types, elemental types, etc.
		std::string value_type_; //the statistical value that is recorded for each energy for each atom type in the above set, e.g. mean, median, mode, boltzmann weighted average.
		std::string score_func_type_; //the base score function in use, specifies the internal two body weights that will be used.
    SplitUnfoldedTwoBodyPotential const & sutbp_;
    EnergyMap score_type_weights_;
    virtual core::Size version() const;

  };

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_methods_SplitUnfoldedTwoBodyEnergy_HH
