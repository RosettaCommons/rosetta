// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/LigandSegment.hh
/// @brief a segment class with an attached ligand residue
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_data_storage_LigandSegment_hh
#define INCLUDED_protocols_sewing_data_storage_LigandSegment_hh


#include <protocols/sewing/data_storage/LigandSegment.fwd.hh>
#include <protocols/sewing/data_storage/LigandResidue.fwd.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>

namespace protocols {
namespace sewing {
namespace data_storage {
///@brief a segment that contains a single ligand residue with given attachment point(s)
class LigandSegment : public SmartSegment {

public:

	LigandSegment(); // "everything is empty" default constructor. This gets used for building segments from files
	LigandSegment(bool is_vital); // "everything is empty" default constructor. This gets used for building segments from files
	LigandSegment(bool is_vital, core::Size max_segment_length); // use this constructor to actually make segments in an assembly, since it preallocates everything
	LigandSegment(LigandSegment const & src);
	///@brief Make a LigandSegment from a non-LigandSegment
	LigandSegment(SmartSegment const & src, bool is_vital, std::set< core::Size > ligand_residues, core::Size ligand_to_attach);

	virtual ~LigandSegment();

	LigandSegmentOP
	clone() const;

	std::set< core::Size >
	get_ligand_residues() const;

	std::set< core::Size > &
	get_nonconst_ligand_residues();

	std::set< core::Size >
	get_ligand_contact_indices() const;

	std::set< core::Size > &
	get_nonconst_ligand_contact_indices();
	void
	add_ligand_contact( core::Size );

	void
	attach_unowned_ligand( LigandResidueCOP ligand_res );

	void
	attach_ligand( LigandResidueOP, bool );

	std::string
	get_name() const override{ return "LigandSegment"; }

	std::map< core::Size, LigandResidueOP > &
	get_owned_ligand_residues();

	std::map< core::Size, LigandResidueOP > const
	get_const_owned_ligand_residues() const;


private:
	std::set< core::Size > ligand_contacts_;
	std::set< core::Size > ligand_residues_;
	std::map< core::Size,  LigandResidueOP > owned_ligand_residues_;
};


} //protocols
} //sewing
} //data_storage


#endif //INCLUDED_protocols_sewing_LigandSegment_hh





