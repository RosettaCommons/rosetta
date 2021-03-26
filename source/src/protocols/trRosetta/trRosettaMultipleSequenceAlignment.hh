// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaMultipleSequenceAlignment.hh
/// @brief A class to store multiple sequence alignment data for input into trRosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_trRosetta_trRosettaMultipleSequenceAlignment_hh
#define INCLUDED_protocols_trRosetta_trRosettaMultipleSequenceAlignment_hh

#ifdef USE_TENSORFLOW

#include <protocols/trRosetta/trRosettaMultipleSequenceAlignment.fwd.hh>

// Basic headers
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace trRosetta {

/// @brief A class to store multiple sequence alignment data for input into trRosetta.
/// @details Must be initialized from the contents of a multiple sequence alignment file.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class trRosettaMultipleSequenceAlignment : public utility::VirtualBase {

public:

	/// @brief Default constructor: explicitly deleted.
	trRosettaMultipleSequenceAlignment() = delete;

	/// @brief Options constructor: initialize this object from the contents
	/// of a multiple sequence alignment.
	trRosettaMultipleSequenceAlignment( std::string const & file_contents );

	/// @brief Destructor.
	~trRosettaMultipleSequenceAlignment() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	trRosettaMultipleSequenceAlignmentOP clone() const;

public: //Public member functions:

	/// @brief Return the contents as input to a trRosetta model.
	utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< int32_t > > construct_input_tensors() const;

private: //Private member functions:

	/// @brief Initialize this object from the contents of a multiple sequence alignment file.
	void initialize_from_file_contents( std::string const & file_contents );

	/// @brief Remove the lowercase characters (i.e. delete them and ligate the flanking regions) from a string.
	/// @details Replaces the input.
	std::string & remove_lcase( std::string & s ) const;

private: //Private member variables:

	/// @brief The inputs into the Tensorflow model.
	/// @details Three tensors: a 1-tensor for number of columns (seq length), a 1-tensor for
	/// number of rows (sequences), and a rowsxcolumns tensor with 1-hot encodings of sequences.
	utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< int32_t > > input_tensors_; 

};

} //trRosetta
} //protocols

#endif //USE_TENSORFLOW

#endif //INCLUDED_protocols_trRosetta_trRosettaMultipleSequenceAlignment_hh
