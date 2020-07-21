// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/PeptideStubMover.hh
/// @brief The PeptideStubMover prepends, appends, or inserts residues into an existing pose, or builds a new polymeric chain
/// @author Yifan Song
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org): Added support for stripping
/// N-acetylation and C-methylamidation when appending residues, preserving phi and the previous
/// omega in the first case and psi and the following omega in the second.
/// @modified Jack Maguire, jackmaguire1444@gmail.com: breaking apply() into smaller methods

#ifndef INCLUDED_protocols_cyclic_peptide_PeptideStubMover_hh
#define INCLUDED_protocols_cyclic_peptide_PeptideStubMover_hh

#include <protocols/moves/Mover.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.fwd.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>

namespace protocols {
namespace cyclic_peptide {

constexpr char DEFAULT_LABEL[] = "PEPTIDE_STUB_EXTENSION";

///@brief This type alias is meant to clarify certain return types.
using FirstResidAdded = core::Size;

enum PSM_StubMode {
	PSM_append,
	PSM_prepend,
	PSM_insert
};

class PeptideStubMover : public moves::Mover {

public:
	PeptideStubMover();
	~PeptideStubMover() override;
	PeptideStubMover( PeptideStubMover const &src );

	void init();

	void apply( Pose & ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override;


	/// @brief Reset mover data
	void reset_mover_data()
	{
		stub_rsd_names_.clear();
		stub_rsd_jumping_.clear();
		stub_rsd_connecting_atom_.clear();
		stub_anchor_rsd_.clear();
		stub_anchor_rsd_connecting_atom_.clear();
	}


	/// @brief Sets whether the pose gets reset (i.e. all residues deleted) or not.
	void set_reset_mode( bool reset_mode )
	{
		reset_ = reset_mode;
	}


	/// @brief Sets whether pdb numbering gets updated or not.
	void set_update_pdb_numbering_mode( bool mode )
	{
		update_pdb_numbering_ = mode;
	}


	/// @brief Adds a residue to the list of residues to be appended, prepended, or inserted.
	/// @details Calls add_residue() override that uses PSM_StubMode.
	void add_residue(
		std::string const &stubmode,
		std::string const &resname,
		core::Size const position,
		bool const jumpmode,
		std::string const &connecting_atom,
		core::Size const repeat,
		core::Size const anchor_rsd,
		core::select::residue_selector::ResidueSelectorCOP anchor_rsd_selector,
		std::string const &anchor_atom
	);

	/// @brief Adds a residue to the list of residues to be appended, prepended, or inserted.
	/// @details This version uses enums.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void add_residue(
		PSM_StubMode const stubmode,
		std::string const &resname,
		core::Size const position,
		bool const jumpmode,
		std::string const &connecting_atom,
		core::Size const repeat,
		core::Size const anchor_rsd,
		core::select::residue_selector::ResidueSelectorCOP anchor_rsd_selector,
		std::string const &anchor_atom
	);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief add label to residues that are created with this mover
	/// @author Jack Maguire, jackmaguire1444@gmail.com
	void
	set_residue_label( std::string const & label ){
		residue_label_ = label;
	}

	static
	std::string
	default_label(){
		return DEFAULT_LABEL;
	}

	static
	void
	assign_chain_ids( core::pose::Pose & pose );

protected:
	// These are all helper functions intended to make apply() smaller
	// The contexts of these functions and their arguments may not be clear in the .hh file but hopefully they are reasonable in the context of the .cc file. I was not the original author of this code so my goal was to keep veriable names unchanged, even if they are unclear as arguments.

	///@brief this performs the inner loop logic of apply()
	///@returns returns resid of the first new residue
	FirstResidAdded perform_single_iteration(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & standard_residues,
		core::Size istub
	);

	///@brief this function performs the early logic of append_by_bond()
	///@details call handle_upper_terminus and handle_lower_terminus
	void handle_termini_and_store_terminal_dihedrals(
		core::pose::Pose & pose,
		core::Size const anchor_rsd,
		core::Size const istub,
		core::Size const connecting_id,
		core::Size & anchor_connecting_id,
		core::Real & old_omega_minus1,
		core::Real & old_phi,
		core::Real & old_psi,
		core::Real & old_omega,
		bool & replace_upper_terminal_type,
		bool & replace_lower_terminal_type
	);

	///@brief determine how residues should connect when being appended by bond
	core::Size
	get_connecting_id_for_append_by_bond(
		core::conformation::Residue const & new_rsd,
		core::Size const istub
	);

	///@brief Handle the case where the new residues start a new chain
	///@returns returns resid of the first new residue added
	FirstResidAdded
	append_by_jump(
		core::pose::Pose & pose,
		core::Size const anchor_rsd,
		core::conformation::Residue & new_rsd,
		core::Size const istub
	);

	///@brief Handle the case where the input pose has zero residues
	///@returns returns resid of the first new residue
	FirstResidAdded
	add_residue_to_empty_pose(
		core::pose::Pose & pose,
		core::conformation::Residue & new_rsd,
		core::Size const istub
	);

	///@brief queries stub_anchor_rsd_ and anchor_rsd_selectors_ to find the residue to use as an achor
	///@returns resid of anchor residue residue
	core::Size
	get_anchor_rsd(
		core::pose::Pose const & pose,
		core::Size const istub
	);

	///@returns returns resid of the first new residue
	FirstResidAdded
	append_by_bond(
		core::pose::Pose & pose,
		core::Size const anchor_rsd,
		core::conformation::Residue & new_rsd,
		core::Size const istub
	);

	///@brief add additional residues if needed
	void
	handle_repeats_in_append_by_bond(
		core::pose::Pose & pose,
		core::Size const anchor_rsd,
		core::conformation::Residue & new_rsd,
		core::Size const istub
	);

private: //Functions

	/// @brief Remove terminal types from the upper terminus.  Store the old psi and omega values.
	/// @returns  Returns void, but if a terminal type was removed, old_psi will be set to the previous
	/// psi value and old_omega will be set to the previous omega value.  The replace_upper_terminal_type var
	/// is set to true if a terminal type was replaced and false otherwise.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	handle_upper_terminus(
		core::pose::Pose & pose,
		core::Size const anchor_rsd,
		core::Real & old_psi,
		core::Real & old_omega,
		bool & replace_upper_terminal_type
	) const;

	/// @brief Remove terminal types from the lower terminus.  Store the old phi and omega_nminus1 values.
	/// @returns  Returns void, but if a terminal type was removed, old_phi will be set to the previous
	/// phi value and old_omega_nminus1 will be set to the previous upstream omega value.  The
	/// replace_lower_terminal_type var is set to true if a terminal type was replaced and false otherwise.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	handle_lower_terminus(
		core::pose::Pose & pose,
		core::Size const anchor_rsd,
		core::Real & old_phi,
		core::Real & old_omega_nminus1,
		bool & replace_lower_terminal_type
	) const;

	/// @brief Update the omega-1 and phi (if we've replaced an N-acetylation) or the psi and omega
	/// (if we've replaced a C-methylamidation) to preserve these dihedral values.
	/// @details Builds a temporary foldtree rooted on the alpha carbon of the anchor residue.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	preserve_old_mainchain_torsions(
		core::pose::Pose & pose,
		core::Size const anchor_res,
		core::Real const old_omega_minus1,
		core::Real const old_phi,
		core::Real const old_psi,
		core::Real const old_omega,
		bool const replace_upper_terminal_type,
		bool const replace_lower_terminal_type
	) const;

private:
	bool reset_;


	/// @brief As residues are added, should the PDB numbering be updated?  Default true.
	bool update_pdb_numbering_;

	utility::vector1<PSM_StubMode> stub_mode_;
	utility::vector1<std::string> stub_rsd_names_;
	utility::vector1<bool> stub_rsd_jumping_;
	utility::vector1<std::string> stub_rsd_connecting_atom_;
	utility::vector1<core::Size> stub_rsd_repeat_;
	utility::vector1<core::Size> stub_insert_pos_;
	utility::vector1<core::Size> stub_anchor_rsd_;
	utility::vector1<core::select::residue_selector::ResidueSelectorCOP> anchor_rsd_selectors_;
	utility::vector1<std::string> stub_anchor_rsd_connecting_atom_;

	std::string residue_label_ = DEFAULT_LABEL;

	//Private functions:


	/// @brief Rebuilds all atoms that are dependent on bonds between residue_index and any other residues (including atoms on the other residues).
	virtual void rebuild_atoms( core::pose::Pose &pose, core::Size const residue_index) const;


	/// @brief Updates the PDB numbering (PDB number/chain ID) as residues are added.
	void update_pdb_numbering (
		core::pose::Pose &pose,
		utility::vector1< core::Size > & resids_that_we_added
	) const;

};

} // moves
} // protocols

#endif
