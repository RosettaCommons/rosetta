// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/EnergyLandscapeEvaluator.hh
/// @brief Evaluates a set of score/rms points
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_toolbox_EnergyLandscapeEvaluator_hh
#define INCLUDED_protocols_toolbox_EnergyLandscapeEvaluator_hh

#include <protocols/toolbox/EnergyLandscapeEvaluator.fwd.hh>

// Protocol headers

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace toolbox {

class ScoreRmsPoint : public utility::pointer::ReferenceCount {
public:
	ScoreRmsPoint( core::Real const score_val, core::Real rms_val );
	core::Real score() const;
	core::Real rms() const;

	friend std::ostream & operator<<( std::ostream & os, ScoreRmsPoint const & score );
private:
	ScoreRmsPoint();
	core::Real score_;
	core::Real rms_;
};

class ScoreRmsPoints : public utility::vector1< ScoreRmsPoint > {
public:
	ScoreRmsPoints( ScoreRmsPoint const & bg_point ); // move-constructed
	ScoreRmsPoint const & bg() const;

private:
	ScoreRmsPoints();
	ScoreRmsPoint bg_;
};

///@brief Evaluates a set of score/rms points
class EnergyLandscapeEvaluator : public utility::pointer::ReferenceCount {
public:
	EnergyLandscapeEvaluator();
	~EnergyLandscapeEvaluator() override;

	virtual EnergyLandscapeEvaluatorOP
	clone() const = 0;

	virtual core::Real
	compute( ScoreRmsPoints const & points ) const = 0;

private:
};

///@brief Evaluates a set of score/rms points using the RotamerBoltzmann original method
class RotamerBoltzmannWeightEvaluator : public EnergyLandscapeEvaluator {
public:
	RotamerBoltzmannWeightEvaluator( core::Real const temperature, bool const include_bg_point_in_sum );
	~RotamerBoltzmannWeightEvaluator() override;

	EnergyLandscapeEvaluatorOP
	clone() const override;

	core::Real
	compute( ScoreRmsPoints const & points ) const override;

private:
	core::Real temperature_;
	bool include_bg_;
};


///@brief Evaluates a set of score/rms points using Vikram's pnear method
///@details PNear = SUM( e^(rms^2/lambda^2)*e^(-(score-bg_score)/temperature) ) / SUM( e^( -(score-bg_score)/temperature ) )
class MulliganPNearEvaluator : public EnergyLandscapeEvaluator {
public:
	MulliganPNearEvaluator( core::Real const temperature, core::Real const lambda );
	~MulliganPNearEvaluator() override;

	EnergyLandscapeEvaluatorOP
	clone() const override;

	core::Real
	compute( ScoreRmsPoints const & points ) const override;

private:
	core::Real temperature_;
	core::Real lambda_sq_;
};

/*
/// @brief computes "modified ddG" of score/rms points using the RotamerBoltzmannWeight method
class ModifiedDdgEvaluator : public EnergyLandscapeEvaluator {
public:
ModifiedDdgEvaluator( ScoreRmsPoints const & ddg_values );
virtual ~ModifiedDdgEvaluator();

virtual EnergyLandscapeEvaluatorOP
clone() const;

virtual core::Real
compute( ScoreRmsPoints const & points ) const;

private:
ScoreRmsPoints const ddg_values_;
};
*/

} //protocols
} //toolbox

#endif //INCLUDED_protocols_toolbox_EnergyLandscapeEvaluator_hh
