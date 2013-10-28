// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/protocols/ligand_docking/ligand_options/Interface.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_InterfaceBuilder_hh
#define INCLUDED_protocols_ligand_docking_InterfaceBuilder_hh

//// Unit Headers
#include <protocols/ligand_docking/InterfaceBuilder.fwd.hh>
#include <protocols/ligand_docking/ligand_options/Interface.fwd.hh>
// AUTO-REMOVED #include <protocols/ligand_docking/LigandArea.hh>

//// Package Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

//// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.hh>

//// C++ headers
#include <map>
// AUTO-REMOVED #include <set>

#include <protocols/ligand_docking/LigandArea.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

// vector of ligand_chain_ids in interface or -1 for residues to pack/minimize not part of interface or 0 for no pack/minimize
class InterfaceBuilder: public utility::pointer::ReferenceCount
{
public:
	InterfaceBuilder();
	virtual ~InterfaceBuilder();
	InterfaceBuilder(utility::vector1<LigandAreaOP> ligand_areas, core::Size extension_window = 0);
	InterfaceBuilder(InterfaceBuilder const & that);

	void parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & datamap
	);

	/* InterfaceBuilder(
			utility::vector1<std::string> const & ligand_chains,
			core::Size extension_window=0
	); */

	ligand_options::Interface build(core::pose::Pose const & pose) const;

	LigandAreas get_ligand_areas() const;

private:
	LigandAreas ligand_areas_;
	core::Size extension_window_;

	void enforce_minimum_length(
		ligand_options::Interface & interface,
		core::pose::Pose const & pose
	) const;

//	void extend_interface(const core::Size residue_id, const core::Size window);

	void find_interface_residues(
		ligand_options::Interface & interface,
		core::pose::Pose const & pose
	) const;

	///@brief First call find_ligand_residues
	void find_protein_residues(
			ligand_options::Interface & interface,
			core::Size ligand_residue_id,
			core::pose::Pose const & pose
	)const;

	void set_interface_residue(
			ligand_options::Interface & interface,
			core::Size const potential_interface_residue_id,
			core::Size const ligand_interface_residue_id,
			core::pose::Pose const & pose
	)const;

	bool is_interface_residue(
			core::conformation::Residue const & potential_interface_residue,
			core::conformation::Residue const & ligand_interface_residue,
			char const chain
	)const;

};

} //namespace ligand_docking
} //namespace protocols

#endif
