// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CrosslinkerMoverCreator.hh
/// @brief This mover links three cysteine residues with a three-way cross-linker.  It adds the crosslinker,
/// sets up constraints, optionally packs and energy-mimizes it into place (packing/minimizing only the crosslinker and
/// the side-chains to which it connects), andthen optionally relaxes the whole structure.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_CrosslinkerMoverCreator_hh
#define INCLUDED_protocols_cyclic_peptide_CrosslinkerMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace cyclic_peptide {

class CrosslinkerMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP
	create_mover() const;

	virtual std::string
	keyname() const;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

};

} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols/cyclic_peptide_CrosslinkerMover_fwd_hh
