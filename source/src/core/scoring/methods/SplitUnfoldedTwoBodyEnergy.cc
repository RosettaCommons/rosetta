// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SplitUnfoldedTwoBodyEnergy.cc
/// @brief  Split unfolded energy two body energy method.
/// @author Riley Simmons-Edler (rse231@nyu.edu)

#include <core/scoring/methods/SplitUnfoldedTwoBodyEnergy.hh>
#include <core/scoring/methods/SplitUnfoldedTwoBodyEnergyCreator.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/SplitUnfoldedTwoBodyPotential.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>

#include <iostream>

namespace core {
namespace scoring {
namespace methods {

  SplitUnfoldedTwoBodyEnergy::SplitUnfoldedTwoBodyEnergy( std::string const & label_type, std::string const & value_type ):
	parent( methods::EnergyMethodCreatorOP( new SplitUnfoldedTwoBodyEnergyCreator ) ),
    label_type_( label_type ),
		value_type_( value_type ),
    sutbp_( ScoringManager::get_instance()->get_SplitUnfoldedTwoBodyPotential( label_type, value_type ) ),
    score_type_weights_( sutbp_.get_weights() )
	{
		//more blank space...
  }

  SplitUnfoldedTwoBodyEnergy::SplitUnfoldedTwoBodyEnergy( std::string const & label_type, std::string const & value_type, const EnergyMap & emap_in ):
	parent( methods::EnergyMethodCreatorOP( new SplitUnfoldedTwoBodyEnergyCreator ) ),
		label_type_( label_type ),
		value_type_( value_type ),
    sutbp_( ScoringManager::get_instance()->get_SplitUnfoldedTwoBodyPotential( label_type, value_type ) ),
    score_type_weights_( emap_in )
	{
		//and again, blank.
  }

	SplitUnfoldedTwoBodyEnergy::~SplitUnfoldedTwoBodyEnergy()
	{
		//this space also intentionally left blank.
	}

  EnergyMethodOP SplitUnfoldedTwoBodyEnergy::clone() const
  {
    return SplitUnfoldedTwoBodyEnergyOP( new SplitUnfoldedTwoBodyEnergy( label_type_, value_type_, score_type_weights_ ) );
  }

  void SplitUnfoldedTwoBodyEnergy::residue_energy(conformation::Residue const & rsd,pose::Pose const &,EnergyMap & emap) const
  {
    EnergyMap energies;

    sutbp_.get_restype_emap(rsd.type(),energies);

		//Now each component has its own score term, which are loaded into the appropriate terms in the emap. While weights from the database file containing the two body values are still loaded, they are no longer applied in favor of weights defined in the global weights file(.wts).
		emap[fa_atr_ref]+=energies[fa_atr_ref];
		emap[fa_rep_ref]+=energies[fa_rep_ref];
		emap[fa_sol_ref]+=energies[fa_sol_ref];
		emap[fa_elec_ref]+=energies[fa_elec_ref];
		emap[hbond_ref]+=energies[hbond_ref];
		emap[dslf_fa13_ref]+=energies[dslf_fa13_ref];

		//Originally the two body split unfolded energy score components were combined into one term using this line:
		emap[split_unfolded_two_body]+=energies.dot(score_type_weights_); //this is obsolete, use *_ref score terms instead.

  }

	void SplitUnfoldedTwoBodyEnergy::indicate_required_context_graphs(utility::vector1<bool> &) const
	{
		//this space intentionally left blank.
	}

  core::Size SplitUnfoldedTwoBodyEnergy::version() const
  {
    return 1;
  }


} // methods
} // scoring
} // core
