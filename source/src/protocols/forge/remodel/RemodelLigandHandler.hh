// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/forge/remodel/RemodelLigandHandler.hh
/// @brief
/// @author Possu Huang ( possu@uw.edu )
///

#ifndef INCLUDED_protocols_forge_remodel_RemodelLigandHandler_hh
#define INCLUDED_protocols_forge_remodel_RemodelLigandHandler_hh

//project headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

namespace protocols {
namespace forge {
namespace remodel {

class RemodelLigandHandler: public protocols::moves::Mover {

private: // typedefs

	typedef protocols::moves::Mover Super;

public: // typedefs

	typedef core::Real Real;
	typedef core::Size Size;

	typedef core::kinematics::MoveMap MoveMap;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef protocols::moves::MoverOP MoverOP;
	typedef core::scoring::constraints::ConstraintSetOP ConstraintSetOP;
	typedef core::scoring::constraints::ConstraintSet ConstraintSet;


public: //constructor/destructor

	RemodelLigandHandler();

	virtual
	~RemodelLigandHandler();

public: // virtual constructors

	virtual
	MoverOP clone() const;

	virtual
	MoverOP fresh_instance() const;

public: // options

public:

	virtual void apply(Pose& pose);
	virtual std::string get_name() const;
	void minimize( Pose & pose);


private: // data

	ScoreFunctionOP cst_sfx_;
	ScoreFunctionOP fullatom_sfx_;


};

} // remodel
} // forge
} // protocols

#endif /* INCLUDED_protocols_forge_remodel_RemodelLigandHandler_HH */

