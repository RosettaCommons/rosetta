// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/data_storage/SmartSewingResidue.hh
/// @brief a minimal container for SEWING residues
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_data_storage_SmartSewingResidue_hh
#define INCLUDED_protocols_sewing_data_storage_SmartSewingResidue_hh

#include <protocols/sewing/data_storage/SmartSewingResidue.fwd.hh>

#include <utility/vector1.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/conformation/Atom.hh>
namespace protocols {
namespace sewing {
namespace data_storage {

///@brief a minimal container for SEWING residues
class SmartSewingResidue : public utility::pointer::ReferenceCount {

public:

	SmartSewingResidue();
	SmartSewingResidue(SmartSewingResidue const & src);

	virtual ~SmartSewingResidue();

	SmartSewingResidueOP
	clone() const;

	void
	set_atom_vector(utility::vector1<core::conformation::Atom> new_atom_vector);

	utility::vector1<core::conformation::Atom> &
	get_atom_vector();

	utility::vector1<core::conformation::Atom> const &
	get_const_atom_vector() const;

	core::conformation::Atom &
	get_atom(core::Size atom_number);
	void
	set_chis(utility::vector1<core::Real> new_chis);

	utility::vector1<core::Real>
	get_chis() const;

	void
	set_amino_acid_type(std::string new_amino_acid_type);

	std::string
	get_amino_acid_type() const;

	void
	set_type( PolymericType );

	PolymericType
	get_type();

	std::string
	get_full_type_name() const;

	void
	set_full_type_name( std::string const & type );



	void
	become(SmartSewingResidueCOP residue_to_become);
private:
	utility::vector1<core::conformation::Atom> atom_vector_; // for proteins NCalphaCO followed by side chain atoms
	utility::vector1<core::Real> chis_;
	std::string amino_acid_type_;
	std::string full_type_name_;
	PolymericType type_;
};


} //protocols
} //sewing
} //data_storage



#endif //INCLUDED_protocols_sewing_data_storage_SmartSewingResidue_hh





