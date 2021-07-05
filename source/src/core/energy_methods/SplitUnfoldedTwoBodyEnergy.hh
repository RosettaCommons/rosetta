// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/SplitUnfoldedTwoBodyEnergy.hh
/// @brief  Header for the split unfolded energy two body component energy method.
/// @author Riley Simmons-Edler (rse231@nyu.edu)


#ifndef INCLUDED_core_energy_methods_SplitUnfoldedTwoBodyEnergy_hh
#define INCLUDED_core_energy_methods_SplitUnfoldedTwoBodyEnergy_hh

#include <core/energy_methods/SplitUnfoldedTwoBodyEnergy.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/SplitUnfoldedTwoBodyPotential.fwd.hh>


#include <utility/vector1.hh>



namespace core {
namespace energy_methods {



class SplitUnfoldedTwoBodyEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy
{

public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy parent;

	//instantiate using the weights from the data file
	SplitUnfoldedTwoBodyEnergy(std::string const & label_type,std::string const & value_type, std::string const & score_func_type);
	//instantiate using given weight emap
	SplitUnfoldedTwoBodyEnergy(std::string const & label_type,std::string const & value_type, std::string const & score_func_type, const core::scoring::EnergyMap & emap_in);
	~SplitUnfoldedTwoBodyEnergy() override;

	core::scoring::methods::EnergyMethodOP clone() const override;

	void residue_energy(conformation::Residue const & rsd,pose::Pose const &, core::scoring::EnergyMap & emap) const override;

	bool minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	void indicate_required_context_graphs(utility::vector1<bool> &) const override;

private:
	std::string label_type_; //the atom type set for this two body energy, e.g. rosetta types, mm types, elemental types, etc.
	std::string value_type_; //the statistical value that is recorded for each energy for each atom type in the above set, e.g. mean, median, mode, boltzmann weighted average.
	std::string score_func_type_; //the base score function in use, specifies the internal two body weights that will be used.
	core::scoring::SplitUnfoldedTwoBodyPotential const & sutbp_;
	core::scoring::EnergyMap score_type_weights_;
	core::Size version() const override;

};

} // scoring
} // core


#endif // INCLUDED_core_energy_methods_SplitUnfoldedTwoBodyEnergy_HH
