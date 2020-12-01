// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/SplitUnfoldedTwoBodyEnergy.cc
/// @brief  Split unfolded energy two body energy method.
/// @author Riley Simmons-Edler (rse231@nyu.edu)

#include <core/energy_methods/SplitUnfoldedTwoBodyEnergy.hh>
#include <core/energy_methods/SplitUnfoldedTwoBodyEnergyCreator.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/SplitUnfoldedTwoBodyPotential.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>

#include <iostream>

namespace core {
namespace energy_methods {


SplitUnfoldedTwoBodyEnergy::SplitUnfoldedTwoBodyEnergy( std::string const & label_type, std::string const & value_type, std::string const & score_func_type ):
	parent( utility::pointer::make_shared< SplitUnfoldedTwoBodyEnergyCreator >() ),
	label_type_( label_type ),
	value_type_( value_type ),
	score_func_type_( score_func_type ),
	sutbp_( core::scoring::ScoringManager::get_instance()->get_SplitUnfoldedTwoBodyPotential( label_type, value_type, score_func_type ) ),
	score_type_weights_( sutbp_.get_weights() )
{
}

SplitUnfoldedTwoBodyEnergy::SplitUnfoldedTwoBodyEnergy( std::string const & label_type, std::string const & value_type,  std::string const & score_func_type, const core::scoring::EnergyMap & emap_in ):
	parent( utility::pointer::make_shared< SplitUnfoldedTwoBodyEnergyCreator >() ),
	label_type_( label_type ),
	value_type_( value_type ),
	score_func_type_(score_func_type ),
	sutbp_( core::scoring::ScoringManager::get_instance()->get_SplitUnfoldedTwoBodyPotential( label_type, value_type, score_func_type ) ),
	score_type_weights_( emap_in )
{
}

SplitUnfoldedTwoBodyEnergy::~SplitUnfoldedTwoBodyEnergy() = default;

core::scoring::methods::EnergyMethodOP SplitUnfoldedTwoBodyEnergy::clone() const
{
	return utility::pointer::make_shared< SplitUnfoldedTwoBodyEnergy >( label_type_, value_type_, score_func_type_, score_type_weights_ );
}

void SplitUnfoldedTwoBodyEnergy::residue_energy(conformation::Residue const & rsd,pose::Pose const &, core::scoring::EnergyMap & emap) const
{
	core::scoring::EnergyMap energies;

	sutbp_.get_restype_emap(rsd.type(),energies);

	//Now each component has its own score term, which are loaded into the appropriate terms in the emap. While weights from the database file containing the two body values are still loaded, they are no longer applied in favor of weights defined in the global weights file(.wts).
	emap[ core::scoring::fa_atr_ref ]    += energies[ core::scoring::fa_atr_ref ];
	emap[ core::scoring::fa_rep_ref ]    += energies[ core::scoring::fa_rep_ref ];
	emap[ core::scoring::fa_sol_ref ]    += energies[ core::scoring::fa_sol_ref ];
	emap[ core::scoring::fa_elec_ref ]   += energies[ core::scoring::fa_elec_ref ];
	emap[ core::scoring::hbond_ref ]     += energies[ core::scoring::hbond_ref ];
	emap[ core::scoring::dslf_fa13_ref ] += energies[ core::scoring::dslf_fa13_ref ];

	//Originally the two body split unfolded energy score components were combined into one term using this line:
	emap[ core::scoring::split_unfolded_two_body ] += energies.dot( score_type_weights_ ); //this is obsolete, use *_ref score terms instead.
}

void SplitUnfoldedTwoBodyEnergy::indicate_required_context_graphs(utility::vector1<bool> &) const
{
}

core::Size SplitUnfoldedTwoBodyEnergy::version() const
{
	return 1;
}


} // scoring
} // core
