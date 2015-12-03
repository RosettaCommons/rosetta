// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/relax/RangeRelaxMover.hh
/// @brief      Relaxes the protein by relaxing in ranges
/// @details Relaxes a protein by iteratively relaxing ranges of the protein;
///    No ramping required. Much faster than FastRelax and good for
///    large to very large proteins (tested up to 5250 residues);
///    For the membrane version, use MPRangeRelax which runs this
///    Mover in the underneath
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_relax_RangeRelaxMover_hh
#define INCLUDED_protocols_relax_RangeRelaxMover_hh

// Unit Headers
#include <protocols/relax/RangeRelaxMover.fwd.hh>
//#include <protocols/relax/RangeRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.fwd.hh>

// Project Headers

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
namespace relax {

using namespace core;
using namespace core::pose;
using namespace protocols::moves;
using namespace core::scoring;
using namespace core::kinematics;
using namespace core::optimization;

class RangeRelaxMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Defaults: dih angle_max = 0.1, nmoves = "nres", movemap = bb and all chi
	RangeRelaxMover();

	/// @brief Custom Constructor
	/// @details Takes the center residue to start relaxing from
	RangeRelaxMover( core::Size center_resnumber );

	/// @brief Copy Constructor
	RangeRelaxMover( RangeRelaxMover const & src );

	/// @brief Assignment Operator
	RangeRelaxMover & operator = ( RangeRelaxMover const & src );

	/// @brief Destructor
	virtual ~RangeRelaxMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	//  /// @brief Create a Clone of this mover
	// virtual protocols::moves::MoverOP clone() const;
	//
	//  /// @brief Create a Fresh Instance of this Mover
	// virtual protocols::moves::MoverOP fresh_instance() const;
	//
	//  /// @brief Pase Rosetta Scripts Options for this Mover
	// void parse_my_tag(
	//       utility::tag::TagCOP tag,
	//       basic::datacache::DataMap &,
	//       protocols::filters::Filters_map const &,
	//       protocols::moves::Movers_map const &,
	//       core::pose::Pose const &
	//       );

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (RangeRelaxMover)
	virtual std::string get_name() const;

	/// @brief Run RangeRelax
	virtual void apply( Pose & pose );

	/// @brief Run AddMembraneMover before?
	/// @details If you want to keep your anchor point for MEM, then pick no
	///   This is only needed for the native, not the pose itself
	void add_membrane_again( bool yesno );

	/// @brief Optimize membrane
	void optimize_membrane( bool yesno );

	/// @brief Idealize pose after run?
	/// @details Might lead to better decoy but takes longer
	void idealize( bool yesno );

	/// @brief Set maximum dihedral angle perturbation
	void set_angle_max( core::Real angle_max );

	/// @brief Set number of moves, can be "nres" or a number
	void set_nmoves( std::string nmoves );

	/// @brief Use scorefunction
	void set_scorefunction( ScoreFunctionOP sfxn );

	/// @brief Set native
	void set_native( PoseOP pose );

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

	/// @brief Finalize setup
	core::kinematics::FoldTree finalize_setup( Pose & pose );

	/// @brief Create constraints to reference pose
	void constrain_to_reference( Pose & pose, Pose & ref_pose );

	/// @brief Repack in a sequence window
	void repack_sequence_window( Pose & pose, MonteCarloOP mc );

	/// @brief Initialize from commandline
	utility::vector1< bool > get_window_repack_residues( Pose & pose, core::SSize center1, core::SSize center2, core::SSize halfrange );

	/// @brief Repack in spherical range
	void repack_spherical_range( Pose & pose, MonteCarloOP mc );

	/// @brief Initialize from commandline
	utility::vector1< bool > get_spherical_repack_residues( Pose & pose, core::Real inner_radius, core::Real outer_radius );

	/// @brief Repack all residues again
	void repack_all( Pose & pose, MonteCarloOP mc );

	/// @brief Idealize pose
	void idealize_pose( Pose & pose, MinimizerOptions minopts, AtomTreeMinimizer atm );

	/// @brief Print score to cout
	void print_score( Pose & pose );

private: // data

	/// @brief Native
	PoseOP native_;

	/// @brief Center residue number
	core::Size center_resnumber_;

	/// @brief Maximum allowed dihedral angle change for Small and ShearMover
	Real angle_max_;

	/// @brief Number of moves Small and ShearMover can make
	/// @details moves_ is a string and can take 'nres' as well as a number
	///    nmoves_ is the actual number that is taken after conversion
	std::string moves_;
	Size nmoves_;

	/// @brief Movemap for Small and ShearMover
	MoveMapOP movemap_;

	/// @brief Scorefxn for scoring
	ScoreFunctionOP sfxn_;

	/// @brief Scorefunction without constraints
	ScoreFunctionOP sfxn_nocst_;

	/// @brief constraint filename
	std::string cst_file_;

	/// @brief constraint weight
	core::Real cst_weight_;

	/// @brief constrain to starting coordinates
	bool cst_to_start_;

	/// @brief constrain to native
	bool cst_to_native_;

	/// @brief Add Membrane to native?
	bool addmem_;

	/// @brief Optimize membrane for protein?
	bool optmem_;

	/// @brief Additional round of repack which packs all residues simultaneously
	bool repack_again_;

	/// @brief Cycles for repacking and minimization
	core::Size cycles_;

	/// @brief Cycles within the Minimizer
	core::Size min_cycles_;

	/// @brief Relax in centric wave pattern
	bool spherical_wave_;

	/// @brief Idealize decoy after run?
	bool idealize_;

};
} // relax
} // protocols

#endif // INCLUDED_protocols_relax_RangeRelaxMover_hh
