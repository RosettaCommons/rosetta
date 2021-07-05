// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added support for trRosetta constraints.


#ifndef INCLUDED_protocols_abinitio_Protocol_hh
#define INCLUDED_protocols_abinitio_Protocol_hh


// Unit Headers
#include <protocols/abinitio/Protocol.fwd.hh>

// Package Headers
#include <protocols/moves/Mover.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/abinitio/KinematicControl.fwd.hh>
#ifdef USE_TENSORFLOW
#include <protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.fwd.hh>
#endif
#include <core/scoring/ScoreFunction.fwd.hh>

// Project Headers
#include <protocols/checkpoint/CheckPointer.hh>
#include <utility/exit.hh>


namespace protocols {
namespace abinitio {

class Protocol : public moves::Mover {
public:
	typedef Mover BaseClass;
	// typedef utility::vector1 < core::pose::Pose > StructureStore;

	Protocol();
	~Protocol() override = default;

	virtual
	void init( core::pose::Pose const& ) {};
	/*{
	bInitialized_ = true;
	};*/

	void
	set_evaluation( evaluation::MetaPoseEvaluatorOP ev );

	void
	add_evaluation( evaluation::PoseEvaluatorOP ev );

	void
	evaluate_pose( core::pose::Pose &pose, std::string tag, core::io::silent::SilentStruct &pss) const;

	evaluation::MetaPoseEvaluatorOP evaluator() {
		if ( !evaluator_ ) {
			evaluator_ = utility::pointer::make_shared< evaluation::MetaPoseEvaluator >();
		}
		return evaluator_;
	}

	virtual void output_debug_structure( core::pose::Pose & pose, std::string prefix );

	void set_kinematics( abinitio::KinematicControlOP kc ) {
		kinematic_control_ = kc;
	}

	abinitio::KinematicControl const& kinematics() {
		runtime_assert( kinematic_control_ != nullptr );
		return *kinematic_control_;
	}

	virtual bool start_from_centroid() const {
		return true;
	}

	virtual bool return_centroid() const {
		return return_centroid_;
	}

	virtual void return_centroid( bool setting ) {
		return_centroid_ = setting;
	}

	void apply( core::pose::Pose& ) override;
	std::string get_name() const override;

	// virtual StructureStore const& structure_store() const {
	//return structure_store_;
	// }

	// virtual StructureStore& structure_store() {
	//  return structure_store_;
	// }

	void set_fullatom_scorefxn( core::scoring::ScoreFunctionOP sfxn ) {
		scorefxn_fa_ = sfxn;
	}

	void set_centroid_scorefxn( core::scoring::ScoreFunctionOP sfxn ) {
		scorefxn_centroid_ = sfxn;
	}

	void set_silentout_file_name( std::string str ) {
		silentout_file_name_ = str;
	}

	/// @brief Set whether we're using trRosetta constraints.
	void set_use_trRosetta_constraints( bool const setting );

	core::scoring::ScoreFunctionOP fullatom_scorefxn() {
		return scorefxn_fa_;
	}

	core::scoring::ScoreFunctionOP centroid_scorefxn() {
		return scorefxn_centroid_;
	}

	static void register_options();

	/// @brief Report whether this protocol supports trRosetta constraints.  The base
	/// class returns false; derived classes can override this and return true.
	virtual bool supports_trRosetta_constraints() const;

private:
	checkpoint::CheckPointer checkpoints_;

public:
	virtual checkpoint::CheckPointer &get_checkpoints() { return checkpoints_; }

	/// @brief Get whether we're using trRosetta constraints.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	inline bool use_trRosetta_constraints() const { return use_trRosetta_constraints_; }

#ifdef USE_TENSORFLOW
	/// @brief Provide an trRosetta constraint generator as input.  Owning pointer is stored
	/// directly, not cloned!
	/// @details Must set use_trRosetta_constraints to true before calling this!
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void set_trRosetta_cst_generator(
		protocols::trRosetta_protocols::constraint_generators::trRosettaConstraintGeneratorOP cst_gen_in
	);
#endif

protected:

#ifdef USE_TENSORFLOW
	/// @brief Allow derived classes to access the constraint generator:
	inline
	protocols::trRosetta_protocols::constraint_generators::trRosettaConstraintGeneratorCOP
	trRosetta_cst_generator() const {
		return trRosetta_cst_generator_;
	}
#endif


private:
	evaluation::MetaPoseEvaluatorOP evaluator_;
	abinitio::KinematicControlOP kinematic_control_;
	//@brief initialized_ is true if init() has been called, false otherwise.
	// bool bInitialized_;

	// StructureStore structure_store_;

	core::scoring::ScoreFunctionOP scorefxn_fa_; //might be NULL
	core::scoring::ScoreFunctionOP scorefxn_centroid_;

	bool return_centroid_;

	//fast fix until the new JobDist/Mover is introduced:
	std::string silentout_file_name_;

	/// @brief Decide whether to use trRosetta constraints.
	/// @details Requires extras=tensorflow or extras=tensorflow_gpu.
	bool use_trRosetta_constraints_ = false;

#ifdef USE_TENSORFLOW
	/// @brief trRosetta constraint generator.  Only used if use_trRosetta_constraints_ is true.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	protocols::trRosetta_protocols::constraint_generators::trRosettaConstraintGeneratorOP trRosetta_cst_generator_;
#endif

};

}
}

#endif
