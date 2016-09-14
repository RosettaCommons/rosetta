// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/HelixFromSequenceCreator.hh
/// @brief Generates a (transmembrane) helix from a fasta file containing the sequence
/// @author Julia Koehler Leman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_HelixFromSequenceCreator_hh
#define INCLUDED_protocols_membrane_HelixFromSequenceCreator_hh

#include <protocols/moves/MoverCreator.hh>


namespace protocols {
namespace membrane {

class HelixFromSequenceCreator : public protocols::moves::MoverCreator {

public:

	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static std::string mover_name();



};

} //protocols
} //membrane

#endif //INCLUDED_protocols/membrane_HelixFromSequence_fwd_hh



