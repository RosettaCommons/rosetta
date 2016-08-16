// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_kinematic_closure_perturbers_LogFilePerturber_HH
#define INCLUDED_protocols_kinematic_closure_perturbers_LogFilePerturber_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/LogFilePerturber.fwd.hh>

// C++ headers
#include <ostream>

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

/// @brief Read torsion angles, bond angles, and bond lengths from a log file.
/// @details This perturber is meant to be used as a debugging aid.
class LogFilePerturber : public Perturber {

public:

	/// @brief Construct using given log file.
	LogFilePerturber(std::string path);

	/// @copydoc Perturber::get_name
	std::string get_name() const { return "LogFilePerturber"; }

	/// @copydoc Perturber::get_subset
	void perturb_subset(
		Pose const & pose,
		IndexList const & residues,
		ClosureProblemOP problem);

	/// @brief Output the given torsion angles in a format that can later be read
	/// by LogFilePerturber.
	static void log_torsions(std::ostream &  out, ParameterList torsions);

	/// @brief Output the given bond angles in a format that can later be read by
	/// LogFilePerturber.
	static void log_angles(std::ostream & out, ParameterList angles);

	/// @brief Output the given bond lengths in a format that can later be read
	/// by LogFilePerturber.
	static void log_lengths(std::ostream & out, ParameterList lengths);

private:

	ParameterMatrix torsion_angles_;
	ParameterMatrix bond_angles_;
	ParameterMatrix bond_lengths_;
	Size iteration_;

};

}
}
}

#endif

