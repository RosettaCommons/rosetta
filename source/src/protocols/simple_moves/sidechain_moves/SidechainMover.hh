// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMover.hh
/// @brief definition of SidechainMover class and functions
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_SidechainMover_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_SidechainMover_hh

// Unit Headers
#include <protocols/simple_moves/sidechain_moves/SidechainMover.fwd.hh>

// Protocols Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/filters/Filter.fwd.hh>

// Core Headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/id/DOF_ID_Range.fwd.hh>

#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>

// Numeric Headers
// AUTO-REMOVED #include <numeric/random/random.hh>

// Utility headers
#include <utility/vector1.hh>

#include <numeric/random/random.fwd.hh>


namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

/// @brief class for non-discrete side chain sampling using Dunbrack rotamer probabilities/distributions
class SidechainMover : public protocols::canonical_sampling::ThermodynamicMover {

public:
//
	/// @brief default constructor
	SidechainMover();

	/// @brief constructor with user provided rotamer library
	SidechainMover(
		core::pack::dunbrack::RotamerLibrary const & rotamer_library
	);

	SidechainMover(
		SidechainMover const & mover
	);

	~SidechainMover();

	virtual
	protocols::moves::MoverOP
	clone() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	/// @brief initialize the packer task if necessary
	void
	init_task(
		core::pose::Pose const & pose
	);

	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	core::conformation::ResidueOP
	make_move( core::conformation::ResidueOP res );

	bool task_initialized();
	/// @brief apply a sidechain move to a Pose object
	void
	apply(
		core::pose::Pose & pose
	);

	virtual std::string get_name() const;

	core::Real
	proposal_density(
		core::conformation::Residue const & proposed_residue,
		core::Size const proposed_resnum,
		core::chemical::ResidueType const & initial_residue_type,
		utility::vector1<core::Real> const & initial_chi_angles
	) const;

	/// @brief test the backrub move
	void
	test_move(
		core::pose::Pose &
	);

	/// @brief idealize sidechains that might be altered
	void
	idealize_sidechains(
		core::pose::Pose & pose
	);

	/// @brief get the rotamer library
	core::pack::dunbrack::RotamerLibrary const &
	rotamer_library() const;

	/// @brief get the task factory
	core::pack::task::TaskFactoryCOP
	task_factory() const;

	/// @brief set the task factory
	void
	set_task_factory(
		core::pack::task::TaskFactoryCOP task_factory
	);

	/// @brief get the packer task
	core::pack::task::PackerTaskCOP
	task() const;

	/// @brief set the task
	void
	set_task(
		core::pack::task::PackerTaskCOP task
	);

	/// @brief get the probability of uniformly sampling chi angles
	core::Real
	prob_uniform() const;

	/// @brief set the probability of uniformly sampling chi angles
	void
	set_prob_uniform(
		core::Real prob_uniform
	);

	/// @brief get whether detailed balance is preserved (i.e. proposal density ratio calculated)
	bool
	preserve_detailed_balance() const;

	/// @brief set whether detailed balance is preserved (i.e. proposal density ratio calculated)
	void
	set_preserve_detailed_balance(
		bool preserve_detailed_balance
	);

	///@brief perform direct chi manipulations rather than using replace_residue to effect rotamer changes; useful if things are kinematically dependent on a sidechain.
	bool
	change_chi_without_replacing_residue() const;

	///@brief perform direct chi manipulations rather than using replace_residue to effect rotamer changes; useful if things are kinematically dependent on a sidechain.
	void
	set_change_chi_without_replacing_residue(
		bool const change_chi_without_replacing_residue
	);

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges(
		core::pose::Pose & pose
	);

	/// @brief get the DOF_IDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::DOF_ID_Range>
	dof_id_ranges(
		core::pose::Pose & pose
	);

	/// @brief get the probability of sampling within the same rotamer
	core::Real
	prob_withinrot() const;

	/// @brief set the probability of sampling within the same rotamer
	void
	set_prob_withinrot(
		core::Real prob_withinrot
	);

	core::Real
	prob_random_pert_current() const;

	void
	set_prob_random_pert_current(
		core::Real prob_pert
	);

	/// @brief get the residues that can be changed by this mover
	utility::vector1<core::Size> const &
	packed_residues() const;

	/// @brief get a vector indicating whether each residue can be packed
	utility::vector1<bool> const &
	residue_packed() const;

	/// @brief get the next residue to be changed, 0 means a random will be chosen
	core::Size
	next_resnum() const;

	/// @brief set the next residue to be changed, 0 means a random will be chosen
	void
	next_resnum(
		core::Size resnum
	);

	/// @brief get the number of chi angles sampled in the last move
	core::Size
	last_nchi() const;

	/// @brief get whether the last move mutated the residue
	bool
	last_mutation() const;

	/// @brief get whether the last move used uniform chi sampling
	bool
	last_uniform() const;

	/// @brief get whether the last move sampled within the same rotamer
	bool
	last_withinrot() const;

	/// @brief get the ratio of proposal densities for the last move
	virtual
	core::Real
	last_proposal_density_ratio();


	/// @brief update string describing the move type
	void
	update_type();

	/// @brief set temperature for bias sampling at dunbrack distribution
  void
  set_sampling_temperature( core::Real temp ){
    sampling_temperature_ = temp;
  }

  core::Real
  sampling_temperature(){
    return sampling_temperature_;
  }

private:

	void
	make_rotwell_jump(
		utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > const & rotamer_sample_data
	);

	void
	preturb_rot_and_dunbrack_eval( core::conformation::ResidueOP input_residue );

	void
	perturb_rot_within_well(
		utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > const & rotamer_sample_data,
		utility::vector1<core::Real> const & previous_chi_angles
	);

	bool
	dunbrack_accept(
		numeric::random::RandomGenerator & Rand,
		core::conformation::Residue & res,
		utility::vector1<core::Real> const & previous_chi_angles,
		utility::vector1<core::Real> const & new_chi_angles
	);

	void
	perturb_chi(
		numeric::random::RandomGenerator & Rand,
		core::Real max_deviation,
		utility::vector1<core::Real> & current_chi,
		utility::vector1<core::Real> & new_chi
	);



private:

	core::pack::dunbrack::RotamerLibrary const & rotamer_library_;
	core::pack::task::TaskFactoryCOP task_factory_;
	core::pack::task::PackerTaskCOP task_;
	core::pose::PoseOP pose_;
	utility::vector1<core::Size> packed_residues_;
	utility::vector1<bool> residue_packed_;
	core::Real prob_uniform_;
	core::Real prob_withinrot_;
	core::Real prob_random_pert_to_current_;
	bool preserve_detailed_balance_;
	bool accept_according_to_dunbrack_;
	bool sample_rotwells_unif_;
	bool change_chi_without_replacing_residue_;
	core::Size next_resnum_;
	utility::vector1<core::Real> last_chi_angles_;
	core::Size last_nchi_;
	bool last_mutation_;
	bool last_uniform_;
	bool last_withinrot_;
	bool last_pertrot_;
	core::Real last_proposal_density_ratio_;
	bool task_initialized_;
	core::pack::dunbrack::RotamerLibraryScratchSpaceOP scratch_;

protected:
	core::Real temperature0_;
	core::Real sampling_temperature_;
}; //SidechainMover


} // sidechain_moves
} // simple_moves
} // protocols

#endif
