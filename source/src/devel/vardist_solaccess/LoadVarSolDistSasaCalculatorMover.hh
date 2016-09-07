// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/vardist_solaccess/LoadVarSolDistSasaCalculatorLegacyMover.hh
/// @brief  VarSolDRotamerDots classes header file
/// @author Andrew Leaver-Fay
/// @author Ron Jacak

#ifndef INCLUDED_devel_vardist_solaccess_LoadVarSolDistSasaCalculatorMover_HH
#define INCLUDED_devel_vardist_solaccess_LoadVarSolDistSasaCalculatorMover_HH

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverCreator.hh>

namespace devel {
namespace vardist_solaccess {

class LoadVarSolDistSasaCalculatorMoverCreator : public protocols::moves::MoverCreator
{
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
};

/// @brief Handles sphere-sphere overlap calculations
class LoadVarSolDistSasaCalculatorMover : public protocols::moves::Mover {

public:
	LoadVarSolDistSasaCalculatorMover(core::Real probe_radius=1.4, core::Real wobble=0.0);
	~LoadVarSolDistSasaCalculatorMover() override;

	protocols::moves::MoverOP clone() const override;
	std::string get_name() const override;
	void apply( core::pose::Pose & p ) override;

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	void parse_my_tag(
		TagCOP const,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) override;

private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real probe_radius_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real wobble_;

};


} // vardist_solaccess
} // core


#endif // INCLUDED_devel_vardist_sollaccess_LoadVarSolDistSasaCalculatorMover_HH
