// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/qsar/scoring_grid/ScoreNormalization.cc
/// @brief  Header file for Score Normalization methods
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>
#include <core/conformation/Residue.hh>
#include <utility/excn/Exceptions.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

ScoreNormalizationOP get_score_normalization_function(std::string norm_tag)
{
	if ( norm_tag == "HeavyAtomNormalization" ) return ScoreNormalizationOP( new HeavyAtomNormalization );
	if ( norm_tag == "AllAtomNormalization" ) return ScoreNormalizationOP( new AllAtomNormalization );
	if ( norm_tag == "ChiAngleNormalization" ) return ScoreNormalizationOP( new ChiAngleNormalization );
	if ( norm_tag == "MolecularWeightNormalization" ) return ScoreNormalizationOP( new MolecularWeightNormalization );

	throw utility::excn::EXCN_RosettaScriptsOption(norm_tag+" is not a valid Score Normalization method");
	return NULL;
}

core::Real HeavyAtomNormalization::operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues)
{
	core::Size n_atoms = 0;
	for ( core::Size i = 1; i <= residues.size(); ++i ) {
		n_atoms += residues[i]->nheavyatoms();
	}

	return input_score/static_cast<core::Real>(n_atoms);
}

core::Real HeavyAtomNormalization::operator()(core::Real const & input_score, core::conformation::Residue const & residue)
{
	return input_score/static_cast<core::Real>(residue.nheavyatoms());
}

core::Real AllAtomNormalization::operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues)
{
	core::Size n_atoms = 0;
	for ( core::Size i = 1; i <= residues.size(); ++i ) {
		n_atoms += residues[i]->natoms();
	}

	return input_score/static_cast<core::Real>(n_atoms);
}

core::Real AllAtomNormalization::operator()(core::Real const & input_score, core::conformation::Residue const & residue)
{
	return input_score/static_cast<core::Real>(residue.natoms());
}

core::Real ChiAngleNormalization::operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues)
{
	core::Size n_chi = 0;
	for ( core::Size i = 1; i <= residues.size(); ++i ) {
		n_chi += residues[i]->nchi();
	}

	if ( n_chi > 0 ) {
		return input_score/static_cast<core::Real>(n_chi);
	} else {
		return input_score;
	}
}

core::Real ChiAngleNormalization::operator()(core::Real const & input_score, core::conformation::Residue const & residue)
{
	if ( residue.nchi() > 0 ) {
		return input_score/static_cast<core::Real>(residue.nchi());
	} else {
		return input_score;
	}
}

core::Real MolecularWeightNormalization::operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues)
{
	core::Real total_mass = 0.0;
	for ( core::Size i = 1; i <= residues.size(); ++i ) {
		total_mass += residues[i]->type().mass();
	}
	return input_score/total_mass;
}

core::Real MolecularWeightNormalization::operator()(core::Real const & input_score, core::conformation::Residue const & residue)
{
	return input_score/residue.type().mass();
}


}
}
}
