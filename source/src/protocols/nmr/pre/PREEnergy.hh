// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pre/PREEnergy.hh
/// @brief   class that calculates energy from NMR paramagnetic relaxation enhancement data
/// @details last Modified: 10/23/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pre_PREEnergy_HH
#define INCLUDED_protocols_nmr_pre_PREEnergy_HH

// Unit headers
#include <protocols/nmr/pre/PREEnergy.fwd.hh>

// Package headers
#include <core/scoring/nmr/pre/PREData.fwd.hh>
#include <core/scoring/nmr/pre/PREMultiSet.fwd.hh>
#include <core/scoring/nmr/pre/PRESingleSet.fwd.hh>
#include <core/scoring/nmr/pre/PRESingle.fwd.hh>
#include <core/scoring/nmr/NMRSpinlabel.fwd.hh>
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.fwd.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// boost headers
#include <boost/unordered/unordered_map.hpp>
#include <boost/functional/hash.hpp>

// C++ headers
#include <iostream>
#include <string>
#include <map>

namespace protocols {
namespace nmr {
namespace pre {

class PREEnergy : public core::scoring::methods::WholeStructureEnergy {

public: // Types

	typedef core::Real Real;
	typedef core::Size Size;
	typedef numeric::xyzVector<core::Real> Vector;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::EnergyMap EnergyMap;

public: // Methods

	/// @brief default constructor
	PREEnergy();

	/// @brief copy constructor
	PREEnergy(PREEnergy const & other);

	/// @brief destructor
	~PREEnergy();

	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

	/// @brief Return PREData from pose. Create PREData if not present and attach them to the pose.
	core::scoring::nmr::pre::PREData &
	get_pre_data_from_pose(Pose & pose) const;

	/// @brief Calculate the PRE score and write it into the pose's total energy EnergyMap
	void
	finalize_total_energy(
		Pose & pose,
		ScoreFunction const & /*sxfn*/,
		EnergyMap & totals
	) const;

	/// @brief Calculate the total PRE score from PREData retrieved from the pose
	Real
	calculate_total_score(
		Pose & pose
	) const;

	/// @brief Called at the beginning of atom tree minimization, this method
	///        allows the derived class the opportunity to initialize pertinent data
	///        that will be used during minimization.
	///        Here, the function creates and updates the atom_id_to_pre_xyz_deriv_map_ which
	///        is needed by the eval_atom_derivative() function.
	virtual
	void
	setup_for_minimizing(
		Pose & pose,
		ScoreFunction const & /*sxfn*/,
		core::kinematics::MinimizerMapBase const & /*minmap*/
	) const;

	/// @brief Evaluate the xyz derivative of the PRE for an atom in the pose.
	///        Called during the atomtree derivative calculation, atom_tree_minimize.cc,
	///        through the ScoreFunction::eval_atom_derivative intermediary.
	///        F1 and F2 should not zeroed, rather, setup_for_minimizing() accumulates its
	///        contribution from the xyz derivatives of atom id
	virtual
	void
	eval_atom_derivative(
		core::id::AtomID const & id,
		Pose const & pose,
		core::kinematics::DomainMap const & /*domain_map*/,
		ScoreFunction const & /*sfxn*/,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	/// @brief Indicate in the context-graphs-required list which
	///        context-graphs this energy method requires that the Pose
	///        maintain when doing neighbor evaluation. Context graphs are allowed.
	void
	indicate_required_context_graphs(
		utility::vector1< bool > &
	) const;

	/// @brief show additional information of the energy method
	void
	show_additional_info(
		std::ostream & TR,
		Pose & pose,
		bool verbose=false
	) const;

	/// @brief Return the version of the energy method
	virtual
	Size
	version() const;

private: // Methods

	/// @brief register options
	void
	register_options();

private: // Data

	// maps AtomID to the xyz derivative of the PRE
	mutable core::id::AtomID_Map< numeric::xyzVector< core::Real> > atom_id_to_pre_xyz_deriv_map_;
};

} // namespace pre
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_pre_PREEnergy_HH

