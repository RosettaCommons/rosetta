// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/pose/carbohydrates/util.cc
/// @brief   Utility function definitions for carbohydrate-containing poses.
/// @author  Labonte

// Unit headers
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/Pose.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>

// Project headers
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>

// Utility headers
//#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>

// External headers
#include <boost/lexical_cast.hpp>


// Construct tracer.
static thread_local basic::Tracer TR( "core.pose.carbohydrates.util" );


namespace core {
namespace pose {
namespace carbohydrates {

using namespace std;
using namespace core;

// Helper Functions ///////////////////////////////////////////////////////////
// TODO: Replace this with code using ResidueType::residue_connection_id_for_atom().
// Scan through a saccharide residue's connections to find the residue from which it follows or branches.
/// @return  The sequence position of the residue before this one (n-1) or the residue in the parent chain from which
/// the branch occurs or zero if N/A, i.e., if this is the lower terminus.
core::uint
find_seqpos_of_parent_residue(conformation::Residue const & residue) {
	uint seqpos = residue.seqpos();
	Size n_connections = residue.n_residue_connections();

	// Search backwards for speed, since "non-polymeric" connections come after polymeric ones.
	for (uint i = n_connections; i >= 1; --i) {
		uint connection_partner_position = residue.residue_connection_partner(i);
		if (connection_partner_position < seqpos) {  // assumes PDB places branches after main chains
			return connection_partner_position;
		}
	}
	TR.Warning << "This residue is a lower terminus! Returning 0." << endl;
	return 0;
}


// Return pointers to the two residues of the glycosidic bond.
/// @return  Pointers to the residue at <sequence_position> and its parent or else the same pointer twice if undefined.
std::pair<conformation::ResidueCAP, conformation::ResidueCAP>
get_glycosidic_bond_residues(Pose const & pose, uint const sequence_position)
{
	using namespace conformation;

	// Get the 1st residue of interest.
	ResidueCAP res_n = & pose.residue(sequence_position);

	if (res_n->is_lower_terminus()) {
		TR.Warning << "Glycosidic torsions are undefined for the first polysaccharide residue of a chain unless part of"
				" a branch." << endl;
		return make_pair(res_n, res_n);
	}

	// Get the 2nd residue of interest.
	// (res_n_minus_1 is a misnomer for the lower termini of branches.)
	ResidueCAP res_n_minus_1 = & pose.residue(find_seqpos_of_parent_residue(*res_n));

	return make_pair(res_n, res_n_minus_1);
}


// Scan through the list of atoms connected to a given "connect atom" and return the first heavy atom found.
/// @return  The atom index of the 1st heavy atom next to the given "connect atom" or an empty string if no heavy atom
/// is found
core::uint
atom_next_to_connect_atom(conformation::Residue const & residue, core::uint const connect_atom_index) {
	using namespace utility;

	// Get list of indices of all atoms connected to given connect atom.
	vector1<uint> const atom_indices = residue.bonded_neighbor(connect_atom_index);

	// Search for heavy atoms.  (A residue connection is not an atom.)
	Size const n_indices = atom_indices.size();
	for (uint i = 1; i <= n_indices; ++i) {
		if (!residue.atom_is_hydrogen(atom_indices[i])) {
			return atom_indices[i];
		}
	}
	return 0;
}


// Return the AtomIDs of the four phi torsion reference atoms.
/// @details For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.\n
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// Et cetera...\n
utility::vector1<id::AtomID>
get_reference_atoms_for_phi(Pose const & pose, uint const sequence_position)
{
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1<AtomID> ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair<ResidueCAP, ResidueCAP> const residues = get_glycosidic_bond_residues(pose, sequence_position);

	if (residues.first->seqpos() == residues.second->seqpos()) {  // This occurs when there is no parent residue.
		return ids;
	}

	// Set the atom names of the four reference atoms.
	// Reference 1 is O(cyclic) for cyclic saccharides, C? for linear saccharides.
	// Because the cyclic O is not connected by the atom tree,
	// we actually will return the AtomID for the virtual atom that superimposes with the real atom.
	AtomID ref1;
	if (residues.first->carbohydrate_info()->is_cyclic()) {
		//ref1 = AtomID(12, residues.first->seqpos());  // TEMP
		ref1 = AtomID(residues.first->carbohydrate_info()->virtual_cyclic_oxygen_index(), residues.first->seqpos());
	} else /* is linear */ {
		;  // TODO: Figure out how linear polysaccharides are handled by IUPAC.
	}
	ids.push_back(ref1);

	// Reference 2 is always the anomeric carbon.
	AtomID ref2(residues.first->carbohydrate_info()->anomeric_carbon_index(), residues.first->seqpos());
	ids.push_back(ref2);

	// Reference 3 is OX(n-1) for polysaccharides.
	AtomID ref3(residues.second->connect_atom(*residues.first), residues.second->seqpos());
	ids.push_back(ref3);

	// Reference 4 is CX(n-1) for polysaccharides.
	AtomID ref4(atom_next_to_connect_atom(*residues.second, ref3.atomno()), residues.second->seqpos());
	ids.push_back(ref4);

	TR.Debug << "Reference atoms for phi: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;

	return ids;
}


// Return the AtomIDs of the four psi torsion reference atoms.
/// @details For saccharides, psi is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1),\n
/// where X is the position of the glycosidic linkage.
utility::vector1<id::AtomID>
get_reference_atoms_for_psi(Pose const & pose, uint const sequence_position)
{
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1<AtomID> ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair<ResidueCAP, ResidueCAP> const residues = get_glycosidic_bond_residues(pose, sequence_position);

	if (residues.first->seqpos() == residues.second->seqpos()) {  // This occurs when there is no parent residue.
		return ids;
	}

	// Set the atom names of the four reference atoms.
	// Reference 1 is always the anomeric carbon.
	AtomID ref1(residues.first->carbohydrate_info()->anomeric_carbon_index(), residues.first->seqpos());
	ids.push_back(ref1);

	// Reference 2 is OX(n-1) for polysaccharides.
	AtomID ref2(residues.second->connect_atom(*residues.first), residues.second->seqpos());
	ids.push_back(ref2);

	// Reference 3 is CX(n-1) for polysaccharides.
	AtomID ref3(atom_next_to_connect_atom(*residues.second, ref2.atomno()), residues.second->seqpos());
	ids.push_back(ref3);

	// Reference 4 is CX-1(n-1) for polysaccharides.
	// TODO: Sadly, I can't think of a better way to do this....
	string const ref3_name = residues.second->atom_name(ref3.atomno());
	uint const ref3_position = atoi(&ref3_name[2]);  // 3rd column (index 2) is the atom number
	uint const ref4_position = ref3_position - 1;
	string const ref4_name = "C" + boost::lexical_cast<string>(ref4_position);
	AtomID ref4(residues.second->atom_index(ref4_name), residues.second->seqpos());
	ids.push_back(ref4);

	TR.Debug << "Reference atoms for psi: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;

	return ids;
}


// Return the AtomIDs of the four omega torsion reference atoms.
/// For carbohydrates glycosylated at an exocyclic position,
/// omega of residue n is defined as OX(n-1)-CX(n-1)-CX-1(n-1)-CX-2(n-1),
/// where X is the position of the glycosidic linkage.
utility::vector1<id::AtomID>
get_reference_atoms_for_1st_omega(Pose const & pose, uint const sequence_position)
{
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1<AtomID> ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair<ResidueCAP, ResidueCAP> const residues = get_glycosidic_bond_residues(pose, sequence_position);

	if (residues.first->seqpos() == residues.second->seqpos()) {  // This occurs when there is no parent residue.
		return ids;
	}
	if (!residues.second->carbohydrate_info()->has_exocyclic_linkage()) {
		TR.Warning << "Omega is undefined for this residue, because the glycosidic linkage is not exocyclic." << endl;
		return ids;
	}

	// Set the atom names of the four reference atoms.
	// Reference 1 is OX(n-1) for polysaccharides.
	AtomID ref1(residues.second->connect_atom(*residues.first), residues.second->seqpos());
	ids.push_back(ref1);

	// Reference 2 is CX(n-1) for polysaccharides.
	AtomID ref2(atom_next_to_connect_atom(*residues.second, ref1.atomno()), residues.second->seqpos());
	ids.push_back(ref2);

	// Reference 3 is CX-1(n-1) for polysaccharides.
	// TODO: Sadly, I can't think of a better way to do this....
	string const ref2_name = residues.second->atom_name(ref2.atomno());
	uint const ref2_position = atoi(&ref2_name[2]);  // 3rd column (index 2) is the atom number
	uint const ref3_position = ref2_position - 1;
	string const ref3_name = "C" + boost::lexical_cast<string>(ref3_position);
	AtomID ref3(residues.second->atom_index(ref3_name), residues.second->seqpos());
	ids.push_back(ref3);

	// Reference 4 is CX-2(n-1) for polysaccharides.
	uint const ref4_position = ref2_position - 2;
	string const ref4_name = "C" + boost::lexical_cast<string>(ref4_position);
	AtomID ref4(residues.second->atom_index(ref4_name), residues.second->seqpos());
	ids.push_back(ref4);

	TR.Debug << "Reference atoms for omega: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;

	return ids;
}


// Return the AtomIDs of the four reference atoms for the requested torsion.
utility::vector1<id::AtomID>
get_reference_atoms(uint const torsion_id, Pose const & pose, uint const sequence_position)
{
	using namespace id;
	using namespace utility;

	vector1<AtomID> ref_atoms;
	switch (torsion_id) {
		case phi_torsion:
			ref_atoms = get_reference_atoms_for_phi(pose, sequence_position);
			break;
		case psi_torsion:
			ref_atoms = get_reference_atoms_for_psi(pose, sequence_position);
			break;
		case omega_torsion:
			ref_atoms = get_reference_atoms_for_1st_omega(pose, sequence_position);
			break;
		default:
			utility_exit_with_message("An invalid torsion angle was requested.");
	}
	return ref_atoms;
}


// Virtual Atom Alignment /////////////////////////////////////////////////////
// Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given
// conformation.
/// @details  This method aligns virtual atom VOX, where X is the position of the cyclic oxygen, OY and HOY, where Y is
/// the position of the anomeric carbon, (provided the residue is not the reducing end, where OY and HOY would be real
/// atoms), and HOZ, where Z is the mainchain glycosidic bond location.  OY and HOY are aligned with the last two "main
/// chain" atoms of the parent residue.  This ensures that torsion angles with duplicate names, e.g., chi1 and phi for
/// internal linked aldoses, will always return the same values.  The same concept applies for HOZ, which aligns with
/// the anomeric carbon of the downstream residue.  This method should be called after any coordinate change for a sac-
/// charide residue and after loading a saccharide residue from a file or sequence for the first time.
/// @note     Do I need to worry about aligning the virtual atoms left over from modified sugar patches?  Can such
/// virtual atoms be deleted?
void
align_virtual_atoms_in_carbohydrate_residue(conformation::Conformation & conf, uint const sequence_position) {
	using namespace id;
	using namespace conformation;

	TR.Debug << " Aligning virtual atoms on residue " << sequence_position << "..." << endl;

	ResidueCAP res = & conf.residue(sequence_position);

	// Find and align VOX, if applicable.
	if (res->carbohydrate_info()->is_cyclic()) {
		TR.Debug << "  Aligning VOX..." << endl;
		uint x = res->carbohydrate_info()->cyclic_oxygen();
		uint OX = res->atom_index(res->carbohydrate_info()->cyclic_oxygen_name());
		uint VOX = res->atom_index("VO" + string(1, x + '0'));

		conf.set_xyz(AtomID(VOX, sequence_position), conf.xyz(AtomID(OX, sequence_position)));
		TR.Debug << "  VOX aligned." << endl;
	}

	// Find and align OY and HOY, if applicable.
	if (!res->is_lower_terminus()) {
		TR.Debug << "  Aligning OY and HOY..." << endl;
		uint y = res->carbohydrate_info()->anomeric_carbon();
		uint OY = res->atom_index("O" + string(1, y + '0'));
		uint HOY = res->atom_index("HO" + string(1, y + '0'));

		uint parent_res_seqpos = find_seqpos_of_parent_residue(*res);
		ResidueCAP parent_res = & conf.residue(parent_res_seqpos);
		uint OY_ref = parent_res->connect_atom(*res);
		uint HOY_ref = atom_next_to_connect_atom(*parent_res, OY_ref);

		conf.set_xyz(AtomID(HOY, sequence_position), conf.xyz(AtomID(HOY_ref, parent_res_seqpos)));
		TR.Debug << "  HOY aligned with atom " << parent_res->atom_name(HOY_ref) <<
				" of residue " << parent_res_seqpos << endl;

		TR.Debug << "   Updating torsions..." << endl;
		ResidueCAP dummy = & conf.residue(sequence_position);  // to trigger private method commented below
		//conf.update_residue_torsions(res->seqpos(), false);
		TR.Debug << "   Torsions updated." << endl;

		conf.set_xyz(AtomID(OY, sequence_position), conf.xyz(AtomID(OY_ref, parent_res_seqpos)));
		TR.Debug << "  OY aligned with atom " << parent_res->atom_name(OY_ref) <<
				" of residue " << parent_res_seqpos << endl;
	}

	// Find and align HOZ(s), if applicable.
	if (!res->is_upper_terminus()) {
		TR.Debug << "  Aligning HOZ..." << endl;
		uint z = res->carbohydrate_info()->mainchain_glycosidic_bond_acceptor();
		uint HOZ = res->atom_index("HO" + string(1, z + '0'));

		uint downstream_res_seqpos = sequence_position + 1;
		ResidueCAP downstream_res = & conf.residue(downstream_res_seqpos);
		uint HOZ_ref = downstream_res->atom_index(downstream_res->carbohydrate_info()->anomeric_carbon_name());

		conf.set_xyz(AtomID(HOZ, sequence_position), conf.xyz(AtomID(HOZ_ref, downstream_res_seqpos)));
		TR.Debug << "  HOZ aligned." << endl;
	}
	Size n_branches = res->carbohydrate_info()->n_branches();
	for (uint branch_num = 1; branch_num <= n_branches; ++branch_num) {
		uint z = res->carbohydrate_info()->branch_point(branch_num);\
		uint OZ = res->atom_index("O" + string(1, z + '0'));
		uint HOZ = res->atom_index("HO" + string(1, z + '0'));

		uint branch_connection_id = res->type().residue_connection_id_for_atom(OZ);
		uint branch_res_seqpos = res->residue_connection_partner(branch_connection_id);
		ResidueCAP branch_res = & conf.residue(branch_res_seqpos);
		uint HOZ_ref = branch_res->atom_index(branch_res->carbohydrate_info()->anomeric_carbon_name());

		conf.set_xyz(AtomID(HOZ, sequence_position), conf.xyz(AtomID(HOZ_ref, branch_res_seqpos)));
	}

	TR.Debug << " All virtual atoms aligned." << endl;
}


// Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given pose.
void
align_virtual_atoms_in_carbohydrate_residue(Pose & pose, uint const sequence_position) {
	align_virtual_atoms_in_carbohydrate_residue(pose.conformation(), sequence_position);
}


// Torsion Access /////////////////////////////////////////////////////////////
// Getters ////////////////////////////////////////////////////////////////////
// Return the requested torsion angle between a saccharide residue of the given pose and the previous residue.
/// @details This method is used in place of Residue::mainchain_torsion() since the main-chain torsions of saccharide
/// residues only make sense in the context of two residues.  Moreover, the reference atoms as defined by the IUPAC are
/// different from the ones that Rosetta uses by default for mainchain torsions for sugars.
/// @param   <torsion_id> is an integer representing the specific torsion angle requested, as defined in
/// core/id/types.hh:\n
///     phi_torsion = 1\n
///     psi_torsion = 2\n
///     omega_torsion = 3
/// @note    I would rather the torsion id were an enum, but as it was already defined, I'm leaving it as a constant
/// for now.
core::Angle get_glycosidic_torsion(uint const torsion_id, Pose const & pose, uint const sequence_position)
{
	using namespace id;
	using namespace numeric;
	using namespace utility;

	vector1<AtomID> const ref_atoms = get_reference_atoms(torsion_id, pose, sequence_position);

	if (ref_atoms.size() == 0) {
		// This occurs when there is no parent residue or when the glycosidic bond is not exocyclic (omega only).
		TR.Warning << "Returning zero." << endl;
		return 0.0;
	}

	Angle const angle_in_radians =
			pose.conformation().torsion_angle(ref_atoms[1], ref_atoms[2], ref_atoms[3], ref_atoms[4]);
	return principal_angle_degrees(conversions::degrees(angle_in_radians));
}


// Setters ////////////////////////////////////////////////////////////////////
// Set the requested torsion angle between a saccharide residue of the given pose and the previous residue.
/// @details This method is used in place of Conformation::set_torsion() since the reference atoms as defined by the
/// IUPAC are different from the ones that Rosetta uses by default for main-chain torsions for sugars.
/// @param   <torsion_id> is an integer representing the specific torsion angle requested, as defined in
/// core/id/types.hh:\n
///     phi_torsion = 1\n
///     psi_torsion = 2\n
///     omega_torsion = 3\n
/// <setting> is in degrees.
/// @note    I would rather the torsion id were an enum, but as it was already defined, I'm leaving it as a constant
/// for now.
void
set_glycosidic_torsion(uint const torsion_id, Pose & pose, uint const sequence_position, core::Angle const setting)
{
	using namespace id;
	using namespace numeric;
	using namespace utility;

	vector1<AtomID> const ref_atoms = get_reference_atoms(torsion_id, pose, sequence_position);

	if (ref_atoms.size() == 0) {
		// This occurs when there is no parent residue or when the glycosidic bond is not exocyclic (omega only).
		return;
	}

	Angle const setting_in_radians = conversions::radians(setting);
	pose.conformation().set_torsion_angle(ref_atoms[1], ref_atoms[2], ref_atoms[3], ref_atoms[4], setting_in_radians);
	align_virtual_atoms_in_carbohydrate_residue(pose, sequence_position);
}

}  // namespace carbohydrates
}  // namespace pose
}  // namespace core
