// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/ChargeGrid.hh
/// @author Sam DeLuca
/// @brief This is an implementation of equation 3 in Goodford, J. Med. Chem. 1985,28,849-857.  doi:10.1021/jm00145a002

#ifndef INCLUDED_protocols_qsar_scoring_grid_ChargeGrid_HH
#define INCLUDED_protocols_qsar_scoring_grid_ChargeGrid_HH

#include <protocols/qsar/scoring_grid/ChargeGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

/// @brief a very light representation of an atom that is just a charge and a cartesian space position
struct ChargeAtom
{
	ChargeAtom(core::Vector const &  in_xyz, core::Real const & in_charge, core::Size const & nc)
	{
		xyz = in_xyz;
		charge = in_charge;
		neighbor_count = nc;

	}

	ChargeAtom()
	{
		charge = 0.0;
		neighbor_count = 0;
	}

	utility::json_spirit::Value serialize();

	void deserialize(utility::json_spirit::mObject data);

	core::Vector xyz;
	core::Real charge;
	core::Size neighbor_count;
};


class ChargeGrid : public SingleGrid
{
public:
	ChargeGrid();
	ChargeGrid(core::Real charge);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & ligand_chain_id_to_exclude);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center);
	virtual void refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> ligand_chain_ids_to_exclude);
	/// @brief serialize the SingleGrid to a json_spirit object
	virtual utility::json_spirit::Value serialize();
	/// @brief deserialize a json_spirit object to a SingleGrid
	virtual void deserialize(utility::json_spirit::mObject data);
	void parse_my_tag(utility::tag::TagCOP const tag);
	/// @brief return the current score of an UltraLightResidue using the current grid
	virtual core::Real score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapOP qsar_map);
	/// @brief return the current score of an atom using the current grid
	virtual core::Real atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapOP qsar_map);
    virtual core::Real score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map);
	/// @brief return the current score of an atom using the current grid
	virtual core::Real atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP qsar_map);

private:

	/// @brief calculate nominal depth based on atom count within 4 Angstroms. See
	/// Table 2 of the 1985 Goodford Paper
	core::Real nominal_depth(core::Size const & n_atoms) const;

	core::Real charge_charge_electrostatic(core::pose::Pose const & pose, ChargeAtom const & atom_q, ChargeAtom const & atom_p) const;

	void setup_charge_atoms(core::pose::Pose const & pose);

private:

	core::Real zeta_;
	core::Real epsilon_;
	core::Real indirect_numerator_; // (zeta - epsilon) / (zeta + epsilon)
	core::Real epsilon_0_;
	std::list<ChargeAtom> charge_atom_list_;

};


}
}
}



#endif /* CHARGEGRID_HH_ */
