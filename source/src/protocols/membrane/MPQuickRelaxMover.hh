// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief      Does a quick relax of a membrane protein
/// @details Uses SmallMover and ShearMover with adjustable maximum dihedral
///    angle changes, then repacking and a single round of minimization
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_MPQuickRelaxMover_hh
#define INCLUDED_protocols_membrane_MPQuickRelaxMover_hh

// Unit Headers
#include <protocols/membrane/MPQuickRelaxMover.fwd.hh>
#include <protocols/membrane/MPQuickRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

class MPQuickRelaxMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Defaults: dih angle_max = 1, nmoves = "nres", movemap = bb and all chi
	MPQuickRelaxMover();

	/// @brief Custom Constructor
	/// @details nmoves is a string because it can either be "nres" or an integer
	MPQuickRelaxMover( core::Real angle_max, std::string nmoves );

	/// @brief Custom Constructor
	MPQuickRelaxMover(
		core::Real angle_max,
		std::string nmoves,
		core::kinematics::MoveMapOP movemap
	);

	/// @brief Copy Constructor
	MPQuickRelaxMover( MPQuickRelaxMover const & src );

	/// @brief Assignment Operator
	MPQuickRelaxMover & operator = ( MPQuickRelaxMover const & src );

	/// @brief Destructor
	virtual ~MPQuickRelaxMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (MPQuickRelaxMover)
	virtual std::string get_name() const;

	/// @brief Flip the downstream partner in the membrane
	virtual void apply( core::pose::Pose & pose );

	/// @brief Run AddMembraneMover before?
	/// @details If you want to keep your anchor point for MEM, then pick no
	void add_membrane_again( bool yesno );

	/// @brief Run MembranePositionFromTopology again?
	/// @details Will change the starting membrane position
	void membrane_from_topology( bool yesno );

	/// @brief Optimize membrane position before relax?
	void optimize_membrane( bool yesno );

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	void register_options();

	/// @brief Set default values
	void set_defaults();

	/// @brief Initialize from commandline
	void init_from_cmd();

private: // data

	/// @brief Native
	core::pose::PoseOP native_;

	/// @brief Maximum allowed dihedral angle change for Small and ShearMover
	core::Real angle_max_;

	/// @brief Number of moves Small and ShearMover can make
	/// @details moves_ is a string and can take 'nres' as well as a number
	///    nmoves_ is the actual number that is taken after conversion
	std::string moves_;
	core::Size nmoves_;

	/// @brief Movemap for Small and ShearMover
	core::kinematics::MoveMapOP movemap_;

	/// @brief Scorefxn
	core::scoring::ScoreFunctionOP sfxn_;

	/// @brief constraint filename
	std::string cst_file_;

	/// @brief constraint weight
	core::Real cst_weight_;

	/// @brief Run AddMembraneMover again?
	/// @details This is stupid: we need this as a workaround for attaching the
	///   MEM at the correct anchor point, neither sliding the jump nor
	///   removing the MEM and re-attaching it works. :(
	bool addmem_;

	/// @brief Run MembranePositionFromTopology?
	bool mem_from_topo_;


	/// @brief Optimize membrane before?
	bool opt_mem_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MPQuickRelaxMover_hh
