// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_cyclic_peptide_PeptideStubMover_hh
#define INCLUDED_protocols_cyclic_peptide_PeptideStubMover_hh

#include <protocols/moves/Mover.hh>

namespace protocols {
namespace cyclic_peptide {

enum StubMode {
    append,
    prepend,
    insert
};

class PeptideStubMover : public moves::Mover {
    
public:
	PeptideStubMover();
	virtual ~PeptideStubMover();

    void init();

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & );
    
private:
    bool reset_;

		///
		/// @brief As residues are added, should the PDB numbering be updated?  Default true.
		bool update_pdb_numbering_;
    
    utility::vector1<StubMode> stub_mode_;
    utility::vector1<std::string> stub_rsd_names_;
    utility::vector1<bool> stub_rsd_jumping_;
    utility::vector1<std::string> stub_rsd_connecting_atom_;
    utility::vector1<core::Size> stub_rsd_repeat_;
    utility::vector1<core::Size> stub_insert_pos_;
    utility::vector1<core::Size> stub_anchor_rsd_;
    utility::vector1<std::string> stub_anchor_rsd_connecting_atom_;

//Private functions:

	///
	/// @brief Rebuilds all atoms that are dependent on bonds between residue_index and any other residues (including atoms on the other residues).
	virtual void rebuild_atoms( core::pose::Pose &pose, core::Size const residue_index) const;

	///
	/// @brief Updates the PDB numbering (PDB number/chain ID) as residues are added.
	virtual void update_pdb_numbering ( core::pose::Pose &pose ) const;

};

} // moves
} // protocols

#endif
