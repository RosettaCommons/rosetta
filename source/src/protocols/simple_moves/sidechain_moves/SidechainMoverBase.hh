// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-white0space:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMoverBase.hh
/// @brief definition of SidechainMoverBase class and functions
/// @author Oliver Lange (oliver.lange@tum.de) adapted from Colin A. Smith's code


#ifndef INCLUDED_protocols_simple_moves_sidechain_moves_SidechainMoverBase_hh
#define INCLUDED_protocols_simple_moves_sidechain_moves_SidechainMoverBase_hh

// Unit Headers
#include <protocols/simple_moves/sidechain_moves/SidechainMoverBase.fwd.hh>

// Protocols Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>

// Core Headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/id/DOF_ID_Range.fwd.hh>

#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

/// @brief class for non-discrete side chain sampling using Dunbrack rotamer probabilities/distributions
class SidechainMoverBase : public protocols::canonical_sampling::ThermodynamicMover {

public:
//
	typedef utility::vector1< core::Real > ChiVector;

	/// @brief default constructor
	SidechainMoverBase();

	/// @brief constructor with user provided rotamer library
	SidechainMoverBase(
		core::pack::dunbrack::RotamerLibrary const & rotamer_library
	);

	SidechainMoverBase(
		SidechainMoverBase const & mover
	);

	~SidechainMoverBase();

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
		core::pose::Pose &pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	virtual core::conformation::ResidueOP
	make_move( core::conformation::ResidueOP res );

	virtual	void
	make_chi_move(
	  core::conformation::Residue const& residue,
		ChiVector const& old_chi,
		ChiVector&  new_chi
	) = 0;

	virtual core::Size suggest_residue_number( core::pose::Pose const& ) const;

	/// @brief apply a sidechain move to a Pose object
	virtual	void apply( core::pose::Pose& pose );

	virtual std::string get_name() const;

	virtual core::Real
	compute_proposal_density(
		core::conformation::Residue const & new_residue,
		core::Size const resnum,
		core::chemical::ResidueType const & old_res_type,
		ChiVector const & old_chi
	) const = 0;

	/// @brief idealize sidechains that might be altered
	void idealize_sidechains( core::pose::Pose& pose );

	/// @brief get the rotamer library
	core::pack::dunbrack::RotamerLibrary const& rotamer_library() const;

	/// @brief get the task factory
	core::pack::task::TaskFactoryCOP task_factory() const;

	/// @brief set the task factory
	void
	set_task_factory( core::pack::task::TaskFactoryCOP task_factory	);

	/// @brief get the packer task
	core::pack::task::PackerTaskCOP	task() const;

	/// @brief set the task
	void set_task(	core::pack::task::PackerTaskCOP task );

	/// @brief get whether detailed balance is preserved (i.e. proposal density ratio calculated)
	bool
	preserve_detailed_balance() const;

	/// @brief set whether detailed balance is preserved (i.e. proposal density ratio calculated)
	void
	set_preserve_detailed_balance( bool setting );

	///@brief perform direct chi manipulations rather than using replace_residue to effect rotamer changes; useful if things are kinematically dependent on a sidechain.
	bool change_chi_without_replacing_residue() const;

	///@brief perform direct chi manipulations rather than using replace_residue to effect rotamer changes; useful if things are kinematically dependent on a sidechain.
	void
	set_change_chi_without_replacing_residue(	bool setting );

	///@brief return true if your last move has mutated residue --- make sure the residue is replaced entirely
	virtual bool have_mutated_residue() const { return false; }

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

	/// @brief get the residues that can be changed by this mover
	utility::vector1<core::Size> const &
	packed_residues() const;

	/// @brief get a vector indicating whether each residue can be packed
	utility::vector1<bool> const &
	residue_packed() const;

	/// @brief get the ratio of proposal densities for the last move
	virtual
	core::Real
	last_proposal_density_ratio();

protected:


	//	void set_last_proposal_density_ratio( core::Real setting ) {
	//		last_proposal_density_ratio_ = setting;
	//	}

private:

	void set_defaults();
	void init_from_options();

	core::pack::dunbrack::RotamerLibrary const & rotamer_library_;
	core::pack::task::TaskFactoryCOP task_factory_;
	core::pack::task::PackerTaskCOP task_;
	utility::vector1<core::Size> packed_residues_;
	utility::vector1<bool> residue_packed_;

	bool preserve_detailed_balance_;
	bool change_chi_without_replacing_residue_;

	utility::vector1<core::Real> last_chi_angles_;
	core::Size last_nchi_;

	core::Real last_proposal_density_ratio_;

}; //SidechainMoverBase


} // sidechain_moves
} // simple_moves
} // protocols

#endif
