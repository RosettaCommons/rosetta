// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SWM_RMSD_Energy.hh
/// @brief  Hack to force SWM to generate native-like conformations
/// @author Arvind Kannan


#ifndef INCLUDED_core_scoring_methods_SWM_RMSD_Energy_hh
#define INCLUDED_core_scoring_methods_SWM_RMSD_Energy_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <map>
#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class SWM_RMSD_Energy : public methods::WholeStructureEnergy  {
public:
	typedef methods::WholeStructureEnergy  parent;

public:

	/// @brief Created in order to generate native-like conformations from SWM.
	SWM_RMSD_Energy();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;
	
	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

/////////////////////////////////
	void
	eval_atom_derivative(
	  id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2 ) const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}
	
	void
	add_to_atom_id_map_after_checks( std::map < core::id::AtomID, core::id::AtomID > & atom_id_map,
									std::string const & atom_name,
									Size const & n1, Size const & n2,
									pose::Pose const & pose1, pose::Pose const & pose2 ) const;

	bool
	mutate_position( pose::Pose & pose, Size const i, char const & new_seq ) const;
		
	void
	superimpose_at_fixed_res( pose::Pose & pose, pose::Pose const & native_pose,
							 Real & rmsd, Size & natoms_rmsd ) const;
	
	void
	superimpose_recursively( pose::Pose & pose, pose::Pose const & native_pose, Real & rmsd, Size & natoms ) const;
	
	Real
	superimpose_at_fixed_res_and_get_all_atom_rmsd( pose::Pose & pose, pose::Pose const & native_pose ) const;
	
private:

	virtual
	core::Size version() const;
	

	
	
/////////////////////////////////////////////////////////////////////////////
// data
/////////////////////////////////////////////////////////////////////////////
private:

	core::pose::Pose native_pose_;

};


}
}
}

#endif // INCLUDED_core_scoring_methods_SWM_RMSD_Energy_HH
