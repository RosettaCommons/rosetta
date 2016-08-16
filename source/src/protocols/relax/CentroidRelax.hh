// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/relax/CentroidRelax.hh
/// @brief Centroid-based relax, using Frank Dimaio's updated centroid statistics
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_relax_CENTROIDRELAX_HH
#define INCLUDED_protocols_relax_CENTROIDRELAX_HH

#include <protocols/relax/CentroidRelax.fwd.hh>

//Core Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/kinematics/MoveMap.hh>

//Protocol Headers
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/moves/Mover.hh>

//Utility Headers
#include <utility/vector1.hh>

namespace protocols {
namespace relax {
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core::scoring;
using namespace protocols::simple_moves;
using namespace core::kinematics;
using namespace core::pack::task;
using core::Size;
using std::string;
using core::Real;
using core::pose::Pose;
using utility::vector1;

/// @brief Relax a pose using Frank Dimaio's smooth centroid statistics.
///Currently under optimization.
///
/// @details Minimize a centroid representation of a pose. Ramp VDW/Rama or both.
///
///May tweak structure by up to ~2.5 A without constraints.
///Use custom constraints or coordinate constraints through relax options for best results.
///Using starting coordinate constraints, structure is tweaked by ~.3/.4 A
///Use increased VDW radii option for bb to improve bb-geometry when not using constraints
///
///
class CentroidRelax : public RelaxProtocolBase {
public:
	CentroidRelax();
	CentroidRelax(MoveMapOP mm);
	CentroidRelax(MoveMapOP mm, ScoreFunctionOP cen_scorefxn_in);

	//Specific functions for Mover...
	virtual ~CentroidRelax();

	virtual string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	void set_defaults();

	void set_rounds(Size rounds);

	/// @brief use larger VDW radii for atoms - bb for now - default True (Courtesy of Frank Dimaio)
	void set_use_increased_vdw_radii(bool use);

	/// @brief Sets to use Rama2b instead of Rama - default True
	void set_use_rama2b(bool use);

	/// @brief Ramp Rama according to centroid relax parameters
	void set_ramp_rama(bool use);

	/// @brief Ramp VDW according to centroid relax parameters
	void set_ramp_vdw(bool use);

	/// @brief Sets main scorefunction used for centroid minimization.
	void set_score_function(ScoreFunctionOP cen_score);

	/// @brief Sets fullatom scorefunction - only used for scoring the full atom pose before and after protocol.
	void set_fa_score_function(ScoreFunctionOP fa_score);

	/// @brief Sets the minimizer type.
	void set_min_type(string min);

	/// @brief Sets to use the cartesian minimizer.
	void set_cartesian(bool cart);

	/// @brief If a fullatom pose is passed, should we repack sidechains according to movemap?
	void do_final_repack(bool repack_sc);

	/// @brief Applies the protocol, See notes
	///
	/// @details Setting ramp_rama and ramp_vdw to false switches to the BASIC protocol
	/// which is rounds of the centroid minmover
	virtual void apply( Pose & pose );

private:

	/// @brief Load the default parameters from the default file.
	void read_default_parameters();

	/// @brief used internally to setup extra stuff for movemap.
	void setup_class_movemap_and_constraints(Pose & pose);

	/// @brief increase VDW radii for backbone to help geometry and decrease rmsd
	void setup_increased_vdw_radii();

	bool use_increased_vdw_radii_;
	bool ramp_rama_;
	bool ramp_vdw_;
	Size rounds_;
	bool cartesian_mode_;
	bool repack_sc_;
	ScoreFunctionOP cen_scorefxn_;
	ScoreFunctionOP fa_scorefxn_;
	MoveMapOP movemap_;
	ScoreType rama_type_;

	/// @brief Container for ramp settings
	struct parameters{
		vector1< Real > vdw_params;
		vector1< Real > rama_params;
		vector1< Real > min_params;
		vector1< Real > cst_params;
	};

	parameters def_parameters;

};
}
}

#endif //#ifndef INCLUDED_protocols/relax_CENTROIDRELAX_HH
