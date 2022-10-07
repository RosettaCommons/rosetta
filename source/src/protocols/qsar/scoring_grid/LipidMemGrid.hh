// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    src/protocols/qsar/scoring_grid/LipidMemGrid.hh
/// @brief   implementation of LipidMemGrid
/// @details Modified 12/11/17
/// @author  Brennica Marlow (brennica.marlow@vanderbilt.edu)

#ifndef INCLUDED_protocols_qsar_scoring_grid_LipidMemGrid_hh
#define INCLUDED_protocols_qsar_scoring_grid_LipidMemGrid_hh

#include <protocols/qsar/scoring_grid/LipidMemGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class LipidMemGrid : public SingleGrid
{
public:

	/// @brief default constructor
	LipidMemGrid();

	/// @brief construct from LipidMemGrid input file
	/// LipidMemGrid(std::string const & filename);

	/// @brief destructor
	~LipidMemGrid() override;

	/// @brief Make a copy of the grid, respecting the subclassing.
	GridBaseOP clone() const override;

	/// @brief populate the grid with LipidMemGrid values in vector based on a passed pose
	void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude) override;

	/// @brief populate the grid with LipidMemGrid values in vector based on a passed pose
	void refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> ligand_chain_ids_to_exclude) override;

	/// @brief populate the grids in the vector with Membrane ligand values based on a passed pose
	void refresh(core::pose::Pose const & pose, core::Vector const & center) override;

	/// @brief return the current score of an UltraLightResidue using the LipidMemGrid
	core::Real score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapCOP qsar_map) const override;

	/// @brief return the current score of an atom using the LipidMemGrid
	core::Real atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapCOP qsar_map) const override;

	/// @brief return the current score of a residue using the LipidMemGrid
	core::Real score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapCOP qsar_map) const override;

	/// @brief return the current score of an atom using the LipidMemGrid
	core::Real atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapCOP qsar_map) const override;

	/// @brief Serialize the LipidMemGrid object into a json_spirit Value
	utility::json_spirit::Value serialize() const override;

	/// @brief deserialize a json spirit Value into a LipidMemGrid object
	void deserialize(utility::json_spirit::mObject data) override;

	/// @setup a LipidMemGrid based on RosettaScripts input
	void parse_my_tag(utility::tag::TagCOP tag) override;

	/// @brief Print a brief summary about this grid to the provided output stream
	void show( std::ostream & out ) const override;

	/// @brief Return a string representing the settings which don't change based on reinitialization
	std::string hash_fingerprint() const override;



	static std::string grid_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	std::string kbpot_file_; /// input lipid energies from database
	numeric::interpolation::spline::InterpolatorCOP kbpot_spline_; /// interpolation
	std::string lip_atom_; /// lipid atom
	core::Real mem_weight_; /// weight

};

}
}
}


#endif /* LIPIDMEMGRID_CC_ */
