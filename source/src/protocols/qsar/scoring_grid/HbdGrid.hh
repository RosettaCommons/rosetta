// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/HbaGrid.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_HbdGrid_hh
#define INCLUDED_protocols_qsar_scoring_grid_HbdGrid_hh

#include <protocols/qsar/scoring_grid/HbdGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

#include <utility/vector1.hh>



namespace protocols {
namespace qsar {
namespace scoring_grid {

class HbdGrid : public SingleGrid
{
public:
	HbdGrid();
	virtual ~HbdGrid();
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> ligand_chain_ids_to_exclude);

	/// @brief return the current score of an UltraLightResidue using the current grid
	virtual core::Real score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapOP qsar_map);
	/// @brief return the current score of an atom using the current grid
	virtual core::Real atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapOP qsar_map);

	virtual core::Real score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map);
	/// @brief return the current score of an atom using the current grid
	virtual core::Real atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP qsar_map);

	/// @brief serialize the SingleGrid to a json_spirit object
	virtual utility::json_spirit::Value serialize();
	/// @brief deserialize a json_spirit object to a SingleGrid
	virtual void deserialize(utility::json_spirit::mObject data);
	void parse_my_tag(utility::tag::TagCOP tag);

private:
	numeric::interpolation::spline::InterpolatorOP lj_spline_;
};

}
}
}

#endif /* HBAGRID_FWD_HH_ */
