// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief   This protocol supercharges the surface of an input pdb with either positive or negatively charged residues.
/// @details There are two modes for supercharging.  The first is called AvNAPSA developed by the David Liu lab at Harvard.  In this approach, surface residues are defined by the Average # Neighbor Atoms Per Sidechain Atom (AvNAPSA value), with a cutoff of 150.  I think 100 is a safer cutoff.  Arg, Lys, Asp, Glu, Asn, Gln are the only residues allowed to mutated.  Lys is always chosen for positive, Glu is always chosen for negative, unless the native is Asn, then Asp is chosen.  Thus, the sequence is deterministic.  If one desires a particular net charge, the residues would be sorted from low to high AvNAPSA value and mutated one at a time until the target charge is achieved - this ignores the ceiling of 150 or 100.  The second approach uses the Rosetta score function to guide the surface mutagenesis.  The user must specifiy if Arg, Lys, or Asp, Glu are desired, and what the reference weights are.  Alternatively, the  user can specify a target net charge, and the reference weights of the charged residues will be incremented/decremented until the net charge is reached.
/// @author Bryan Der
/// @author Rocco Moretti (rmorettiase@gmail.com) (Moverization)

#ifndef INCLUDED_protocols_design_opt_Supercharge_hh
#define INCLUDED_protocols_design_opt_Supercharge_hh

#include <protocols/design_opt/Supercharge.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

#include <set>

namespace protocols {
namespace design_opt {

/// @brief Adds charged residues to a protein surface
/// @details There are two modes for supercharging.  The first is called AvNAPSA developed by the David Liu lab at Harvard.  In this approach, surface residues are defined by the Average # Neighbor Atoms Per Sidechain Atom (AvNAPSA value), with a cutoff of 150.  I think 100 is a safer cutoff.  Arg, Lys, Asp, Glu, Asn, Gln are the only residues allowed to mutated.  Lys is always chosen for positive, Glu is always chosen for negative, unless the native is Asn, then Asp is chosen.  Thus, the sequence is deterministic.  If one desires a particular net charge, the residues would be sorted from low to high AvNAPSA value and mutated one at a time until the target charge is achieved - this ignores the ceiling of 150 or 100.  The second approach uses the Rosetta score function to guide the surface mutagenesis.  The user must specifiy if Arg, Lys, or Asp, Glu are desired, and what the reference weights are.  Alternatively, the  user can specify a target net charge, and the reference weights of the charged residues will be incremented/decremented until the net charge is reached.
///
/// AvNAPSA-mode, target charge
/// 1. Define surface.  sort NQ and RK/DE residues by AvNAPSA value (low to high)
/// 2. Next residue in sorted list: Positive: mutate DENQ-->K, Negative: mutate RKQ-->E and N-->D
/// 3. If net charge = target net charge, output pdb
///
/// AvNAPSA-mode, no target charge
/// 1. Define surface by AvNAPSA value (<100 default)
/// 2. For each NQ and DE/RK residue in the surface: Positive: mutate DENQ-->K, Negative: mutate RKQ-->E and N-->D
/// 3. Output pdb
///
/// Rosetta-mode, target charge
/// 1. Define surface.  Neighbor by distance calculator (CB dist.), <16 neighbors default
///  or Define surface by AvNAPSA value (<100 default)
/// 2. Set design task
///    read user resfile, if provided
///    dont_mutate gly, pro, cys
///    dont_mutate h-bonded sidechains
///    dont_mutate correct charge residues
/// 3. Set reference energies for RK/DE, starting at user input values
/// 4. pack rotamers mover
/// 5. check net charge, increment/decrement reference energies (back to step 3.)
/// 6. Once a pack rotamers run results in the correct net charge, output pdb
///
/// Rosetta-mode, no target charge
/// 1. Define surface.  Neighbor by distance calculator (CB dist.), <16 neighbors default
///  or Define surface by AvNAPSA value (<100 default)
/// 2. Set design task
///    read user resfile, if provided
///    dont_mutate gly, pro, cys
///    dont_mutate h-bonded sidechains
///    dont_mutate correct charge residues
/// 3. Set reference energies for RK/DE, using the user input values
/// 4. pack rotamers mover
/// 5. Output pdb
///
///Note: either mode will read a user input resfile, but be sure to use ALLAA as the default, because NATAA default will make the surface residues undesignable.  Either mode will make a second resfile with NATAA as default.
///
/// CAUTION: Supercharge does it's own output.
class Supercharge : public protocols::moves::Mover {
public:
	using SizeSet = std::set<core::Size>;

	Supercharge() = default;
	~Supercharge() override = default;

	void AvNAPSA_positive(bool s) { AvNAPSA_positive_ = s; }
	void AvNAPSA_negative(bool s) { AvNAPSA_negative_ = s; }
	void target_net_charge_active(bool s) { target_net_charge_active_ = s; }
	void target_net_charge(int s) { target_net_charge_ = s; }

	/// @details if target_net_charge is specified, the AvNAPSA cutoff is ignored
	void surface_atom_cutoff(core::Size s) {
		surface_atom_cutoff_ = s;
		surface_atom_cutoff_set_ = true;
	}

	void surface_residue_cutoff(core::Size s) { surface_residue_cutoff_ = s; }
	void include_arg(bool s) { include_arg_ = s; }
	void include_lys(bool s) { include_lys_ = s; }
	void include_asp(bool s) { include_asp_ = s; }
	void include_glu(bool s) { include_glu_ = s; }
	void refweight_arg(core::Real s) { refweight_arg_ = s; }
	void refweight_lys(core::Real s) { refweight_lys_ = s; }
	void refweight_asp(core::Real s) { refweight_asp_ = s; }
	void refweight_glu(core::Real s) { refweight_glu_ = s; }
	void dont_mutate_glyprocys(bool s) { dont_mutate_glyprocys_ = s; };
	void dont_mutate_correct_charge(bool s) { dont_mutate_correct_charge_ = s; }
	void dont_mutate_hbonded_sidechains(bool s) { dont_mutate_hbonded_sidechains_ = s; }
	void pre_packminpack(bool s) { pre_packminpack_ = s; }
	void local_nstruct(core::Size s) { local_nstruct_ = s; }

	void compare_residue_energies_all(bool s) { compare_residue_energies_all_ = s; }
	void compare_residue_energies_mut(bool s) { compare_residue_energies_mut_ = s; }


	void
	apply( Pose & pose ) override;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////         BEGIN AVNAPSA MODE (average number of neighboring atoms per sidechain atom)   //////////////////////////////////////
	//////////          doesn't consider surface energetics, operates only on surface accessibility   //////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/// @brief average number of neighboring atoms per sidechain atom
	virtual
	void
	AvNAPSA_values( Pose const & pose );

	virtual
	void
	set_resfile_AvNAPSA( Pose const & pose );

	virtual
	void
	design_supercharge_AvNAPSA( Pose const & starting_pose, Pose & pose );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////         BEGIN ROSETTA MODE               /////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual
	void
	prepack_input_structure( Pose & pose );

	void set_surface( Pose const & pose );

	void set_resfile( Pose const & pose );

	utility::vector1< core::Real >
	set_reference_energies();

	void
	design_supercharge( Pose const & starting_pose, Pose & pose );

	void
	print_netcharge_and_mutations( Pose const & starting_pose, Pose const & pose );

	int
	get_net_charge( Pose const & pose );

	void
	energy_comparison( Pose & native, Pose & pose );

	std::string
	get_name() const override { return "supercharge"; }

private:
	// options:
	// AvNAPSA-mode
	bool AvNAPSA_positive_ = false;
	bool AvNAPSA_negative_ = false;
	bool target_net_charge_active_ = false;
	int target_net_charge_ = 0;

	// AvNAPSA-mode or Rosetta-mode
	// if target_net_charge is specified, the AvNAPSA cutoff is ignored
	core::Size surface_atom_cutoff_ = 120;
	bool surface_atom_cutoff_set_ = false;

	// Rosetta-mode (these will be ignored if AvNAPSA mode is on via AvNAPSA_positive or AvNAPSA_negative)
	core::Size surface_residue_cutoff_ = 16; // for choosing surface residues, cannot be done in AvNAPSA mode
	bool include_arg_ = false;
	bool include_lys_ = false;
	bool include_asp_ = false;
	bool include_glu_ = false;
	core::Real refweight_arg_ = -0.14916;
	core::Real refweight_lys_ = -0.287374;
	core::Real refweight_asp_ = -1.28682;
	core::Real refweight_glu_ = -1.55374;
	bool dont_mutate_glyprocys_ = true;
	bool dont_mutate_correct_charge_ = true;
	bool dont_mutate_hbonded_sidechains_ = true;
	bool pre_packminpack_ = false;
	core::Size local_nstruct_ = 1;

	// AvNAPSA-mode or Rosetta-mode
	bool compare_residue_energies_all_ = false;
	bool compare_residue_energies_mut_ = true;


	// Data:

	SizeSet surface_res_;
	std::string outputname_;
	std::string out_path_;
	utility::vector1< core::Real > AvNAPSA_values_;
	core::Size largest_mutated_AvNAPSA_;

	core::scoring::ScoreFunctionOP scorefxn_;
};

} // moves
} // protocols

#endif
