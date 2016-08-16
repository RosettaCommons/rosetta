// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/SplitUnfoldedTwoBodyPotential.hh
/// @brief  Reads in and stores the two body energies for each residue type for the split unfolded energy.
/// @author Riley Simmons-Edler (rse231@nyu.edu)


#ifndef INCLUDED_core_scoring_SplitUnfoldedTwoBodyPotential_hh
#define INCLUDED_core_scoring_SplitUnfoldedTwoBodyPotential_hh

#include <core/scoring/SplitUnfoldedTwoBodyPotential.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ResidueType.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>
#include <map>


namespace core {
namespace scoring {


class SplitUnfoldedTwoBodyPotential : public utility::pointer::ReferenceCount
{

public:

	SplitUnfoldedTwoBodyPotential(std::string filename);

	SplitUnfoldedTwoBodyPotential(std::string filename,std::string atom_type_label_name);

	//places the emap of two body energies for the given residue name into the given emap
	void get_restype_emap(const chemical::ResidueType & restype,EnergyMap & emap) const;

	//returns the weights read in from the atom type energy file(or wherever this comes from later on)
	EnergyMap get_weights() const;

private:

	//read in atom type energies from the given file.
	void read_database_file(std::string filename);

	//calculates are returns an emap for the given residue type.
	EnergyMap calculate_residue_emap(const chemical::ResidueType & restype) const;

	std::map<std::string,core::scoring::EnergyMap> residue_two_body_energies_;

	std::map<std::string,core::scoring::EnergyMap> atom_two_body_energies_;

	EnergyMap residue_score_term_weights_;

	//holds the name of the label set in use by the loaded database file, e.g. "mm", "rosetta", etc.
	std::string atom_type_label_set_used_;

};

} //scoring
} //core

#endif // INCLUDED_core_scoring_SplitUnfoldedTwoBodyPotential_HH
