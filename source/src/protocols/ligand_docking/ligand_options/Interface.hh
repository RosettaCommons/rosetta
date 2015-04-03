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

#ifndef INCLUDED_protocols_ligand_docking_ligand_options_Interface_hh
#define INCLUDED_protocols_ligand_docking_ligand_options_Interface_hh

//// Unit Headers
#include <protocols/ligand_docking/ligand_options/Interface.fwd.hh>

//// Project Headers

//// Utility Headers
#include <core/types.hh>

// STL Headers
#include <sstream>

#include <utility/vector1_bool.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {
namespace ligand_options {


/// @brief info for each residue- is it part of the interface and if so, what ligands is it near
struct InterfaceInfo{
	enum Type{
		non_interface,
		near_interface,
		is_interface
	};
	Type type;
	core::Size chain_id; // chains that are causing this to be part of the interface
	char chain;

	InterfaceInfo(): type(non_interface){}
	InterfaceInfo(Type t): type(t){};
};

/// @brief For each residue is it in the interface, a mobile region or a non-mobile region?
class Interface: public utility::vector1<InterfaceInfo>{

public:
	Interface(core::Size num, InterfaceInfo info);

	core::Size find_first_interface_residue(core::Size chain_begin, core::Size const chain_end) const;

	core::Size find_start_of_next_interface_region(
			core::Size start_from,
			core::Size const chain_end
	)const;

	core::Size find_stop_of_this_interface_region(
			core::Size start_from,
			core::Size const chain_end
	)const;

	void enforce_minimum_length(
			core::Size const chain_begin,
			core::Size const chain_end,
			core::Size const window
	);

	utility::vector1<core::Size> get_interface_residues() const;
	utility::vector1<core::Size> get_near_interface_residues() const;
	std::string get_python_string() const;

	//static int test();

private:

	void set_interface_residue(
		core::Size const potential_interface_residue_id,
		core::Size const ligand_interface_residue_id
	);
};

std::ostream & operator<<(std::ostream& output, Interface const & interface);

} //namespace ligand_options
} //namespace ligand_docking
} //namespace protocols

#endif
