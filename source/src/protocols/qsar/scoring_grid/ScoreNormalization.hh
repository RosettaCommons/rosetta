// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/qsar/scoring_grid/ScoreNormalization.hh
/// @brief  Header file for Score Normalization methods
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#ifndef INCLUDED_protocols_qsar_scoring_grid_ScoreNormalization_HH
#define INCLUDED_protocols_qsar_scoring_grid_ScoreNormalization_HH

#include <protocols/qsar/scoring_grid/ScoreNormalization.fwd.hh>
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <utility/VirtualBase.hh>

#include <string>

namespace protocols {
namespace qsar {
namespace scoring_grid {

ScoreNormalizationOP get_score_normalization_function(std::string const & norm_tag);

class ScoreNormalization : public utility::VirtualBase {
public:
	ScoreNormalization() {};
	~ScoreNormalization() override {};
	virtual std::string get_name() const = 0;
	virtual core::Real operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues) const = 0;
	virtual core::Real operator()(core::Real const & input_score, core::conformation::Residue const & residue) const = 0;
private:
	ScoreNormalization(ScoreNormalization const & src) = delete;
};

class HeavyAtomNormalization : public ScoreNormalization {
public:
	HeavyAtomNormalization() {};
	~HeavyAtomNormalization() override {};

	std::string get_name() const override
	{
		return "HeavyAtomNormalization";
	}

	core::Real operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues) const override;
	core::Real operator()(core::Real const & input_score, core::conformation::Residue const & residue) const override;
private:
	HeavyAtomNormalization(HeavyAtomNormalization const & src) = delete;
};

class AllAtomNormalization : public ScoreNormalization {
public:
	AllAtomNormalization() {};
	~AllAtomNormalization() override {};

	std::string get_name() const override
	{
		return "AllAtomNormalization";
	}

	core::Real operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues) const override;
	core::Real operator()(core::Real const & input_score, core::conformation::Residue const & residue) const override;
private:
	AllAtomNormalization(AllAtomNormalization const & src) = delete;
};

class ChiAngleNormalization : public ScoreNormalization {
public:
	ChiAngleNormalization() {};
	~ChiAngleNormalization() override {};

	std::string get_name() const override
	{
		return "ChiAngleNormalization";
	}

	core::Real operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues) const override;
	core::Real operator()(core::Real const & input_score, core::conformation::Residue const & residue) const override;
private:
	ChiAngleNormalization(ChiAngleNormalization const & src) = delete;
};

class MolecularWeightNormalization : public ScoreNormalization {
public:
	MolecularWeightNormalization() {};
	~MolecularWeightNormalization() override {};

	std::string get_name() const override
	{
		return "MolecularWeightNormalization";
	}

	core::Real operator()(core::Real const & input_score, core::conformation::ResidueCOPs residues) const override;
	core::Real operator()(core::Real const & input_score, core::conformation::Residue const & residue) const override;
private:
	MolecularWeightNormalization(MolecularWeightNormalization const & src) = delete;
};

}
}
}

#endif /* SCORENORMALIZATION_HH_ */
