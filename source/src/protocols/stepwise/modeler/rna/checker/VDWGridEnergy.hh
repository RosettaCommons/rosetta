// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/scoring/methods/rna/VDWGridEnergy.hh
/// @brief  VDWGridEnergy energy method
/// @author Kalli Kappel


#ifndef INCLUDED_protocols_stepwise_modeler_rna_checker_VDWGridEnergy_hh
#define INCLUDED_protocols_stepwise_modeler_rna_checker_VDWGridEnergy_hh

// Unit headers
#include <protocols/stepwise/modeler/rna/checker/VDWGridEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_CachedRepScreenInfo.fwd.hh>
#include <core/pose/rna/VDW_RepScreenInfo.fwd.hh>
#include <core/pose/rna/VDW_Grid.fwd.hh>
#include <core/id/AtomID_Map.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {


class VDWGridEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;

public:

	/// @brief ctor
	VDWGridEnergy();

	/// @brief dtor
	virtual ~VDWGridEnergy();

	/// clone
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for WholeStructureEnergies
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	finalize_total_energy(
		core::pose::Pose &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const;

	/// @brief VDWGridEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	virtual
	core::Size version() const;


	// data
private:

	core::Real const clash_penalty_;

};

} //checker
} //rna
} //modeler
} //stepwise
} //core


#endif // INCLUDED_protocols_stepwise_modeler_rna_checker_VDWGridEnergy_HH
