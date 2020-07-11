// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/PackerTask_.hh
/// @brief  Implementation class for task class to describe packer's behavior header
/// Almost all of rosetta needs to use packer tasks, but very little of rosetta needs
/// to see how it behaves internally.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Steven Lewis
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added support for multithreading.

#ifndef INCLUDED_core_pack_task_PackerTask__hh
#define INCLUDED_core_pack_task_PackerTask__hh

// Unit Headers
#include <core/pack/task/PackerTask_.fwd.hh>
#include <core/pack/task/ResidueLevelTask_.hh>

// Package Headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/rna/RNA_ResidueLevelTask.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// STL Headers
#include <iosfwd>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {

enum PackerTaskSymmetryStatus {
	NO_SYMMETRIZATION_REQUEST,
	REQUEST_SYMMETRIZE_BY_UNION,
	REQUEST_SYMMETRIZE_BY_INTERSECTION,
	ALREADY_SYMMETRIZED
};


/////////////////////////////////////////////////////////////////////////////
/// @brief the PackerTask controls what rotameric (including sequence) changes the packer is allowed to make
class PackerTask_ : public PackerTask
{
public:
	/// @brief constructor; the PackerTask will always need a pose!
	/// @note Can optionally specify the number of interaction graph threads to request.  The default (0)
	/// means to use all available threads.  In the single-threaded build, must always be 0 or 1.
	PackerTask_( core::Size const ig_threads_to_request = 0 );

	/// @brief Initialization constructor that needs a pose.  A DefaultPackerPalette
	/// is automatically created using the ResidueTypeSet of the pose.
	/// @note Can optionally specify the number of interaction graph threads to request.  The default (0)
	/// means to use all available threads.  In the single-threaded build, must always be 0 or 1.
	PackerTask_(
		pose::Pose const & pose,
		core::Size const ig_threads_to_request = 0
	);


	/// @brief Initialization constructor that needs a pose and a PackerPalette.
	/// @details The input PackerPalette is cloned.
	/// @note Can optionally specify the number of interaction graph threads to request.  The default (0)
	/// means to use all available threads.  In the single-threaded build, must always be 0 or 1.
	PackerTask_(
		pose::Pose const & pose,
		core::pack::palette::PackerPaletteCOP const & packer_palette,
		core::Size const ig_threads_to_request = 0
	);

	/// @brief dtor
	~PackerTask_() override;

	/// @brief copy method
	PackerTaskOP clone() const override;

	/// @brief Check that the number of threads to request is reasonable.
	/// @details Does nothing in multithreaded build.  In single-threaded build, checks
	/// that number of threads to request is 0 (use all available, which is 1) or 1.
	void check_threads() const;

	/// @brief replace a given residue task with a brand new one NOTE: This should be the only way to break commutativity!!!!
	void clean_residue_task( conformation::Residue const & original_residue, Size const seqpos, core::pose::Pose const & pose) override;

	/// @brief number of residues in the input pose, for convienience (PackerTask does not handle variable length)
	Size total_residue() const override;

	/// @brief get function: can this position have a rotamer change?
	bool pack_residue( int resid ) const override;
	/// @brief alias for above
	bool being_packed( Size resid ) const override { return pack_residue( resid ); }

	/// @brief get function: how many positions can have rotamer changes?
	Size num_to_be_packed() const override;

	/// @brief get function: can this position have a sequence change?
	bool design_residue( int resid ) const override;
	/// @brief alias for above
	bool being_designed( Size resid ) const override { return design_residue( resid ); }
	/// @brief get function: can any positions have a sequence change?
	bool design_any() const override;

	/// @brief set function: bump_check is activated for pack-rotamers' screening of rotamers that
	/// collide with the background.  Bump-check is not used during rotamer trials, since it is
	/// nearly as expensive as rotamer trials itself.  Energy methods may opt in to the bump check
	/// process.  The "standard" behavior is for bump-check to include only the fa_atr and fa_rep
	/// terms.
	void set_bump_check( bool setting ) override;
	/// @brief get function: has bump_check been requested?
	bool bump_check() const override;
	/// @brief Decrease the max_rotbump_energy threshold above which rotamers are rejected.
	void and_max_rotbump_energy( Real setting ) override;
	/// @brief get function: what is the energy threshold above which rotamers should be rejected
	Real max_rotbump_energy() const override;

	/// @brief for all positions, turn on include_current if false, do nothing if already true
	void or_include_current( bool setting ) override;
	/// @brief for one position, turn on include_current if false, do nothing if already true
	void or_include_current( bool setting, Size resid ) override;
	/// @brief get function: what is include_current for this residue?
	bool include_current( Size resid ) const override;

	void add_behavior( std::string const & behavior ) override;
	void add_behavior( std::string const & behavior, Size resid ) override;
	bool has_behavior( std::string const & behavior, Size resid ) const override;
	bool has_behavior( Size resid ) const override;

	/// @brief return the targeted type (may be null pointer)
	virtual chemical::ResidueTypeCOP target_type( Size resid ) const;

	/// @brief for all positions, disable adducts if false, do nothing if true
	void or_adducts( bool setting ) override;
	/// @brief for one position, disable adducts if false, do nothing if true
	void or_adducts( bool setting, Size resid ) override;
	/// @brief include adducts at this residue
	bool adducts( Size resid ) const override;

	/// @brief if setting == true, turns on optimize_H_mode for all residues
	void or_optimize_h_mode( bool setting ) override;

	/// @brief if setting == true, preserves c-beta during rotamer building for all residues
	void or_preserve_c_beta( bool setting ) override;

	/// @brief if setting == true, turns on optimize_H_mode and flip_HNQ for all residues
	void or_flip_HNQ( bool setting ) override;

	/// @brief if setting == true, fix his tautomer state for defined residues during repacking or optimizeH mode
	void or_fix_his_tautomer( utility::vector1<int> const & positions, bool setting ) override;

	/// @brief if setting == true, turns on linear-memory interaction graph usage
	void or_linmem_ig( bool setting ) override;
	/// @brief if setting == false, turns off linear-memory interaction graph usage
	void and_linmem_ig( bool setting ) override;
	/// @brief returns the linear-memory interaction graph flag
	bool linmem_ig() const override;

	/// @brief Set the linear memory interaction graph's recent history size.  Default is 10.  It can be set
	/// once to any value larger than 10, but any subsequent setting will only alter the size if it decreases
	/// it.
	void decrease_linmem_ig_history_size( Size setting ) override;

	/// @brief Return the linear memory interaction graph's recent history size.
	Size linmem_ig_history_size() const override;


	void or_precompute_ig( bool setting ) override;
	//void and_precompute_ig( bool setting ) override;
	bool precompute_ig() const override;


	/// @brief  if setting == true, turns on lazy interaction graph usage
	/// NOTE: the linear memory interaction graph takes precedence over the LazyIG when
	/// the InteractionGraphFactory examines the PackerTask.
	void or_lazy_ig( bool setting ) override;
	/// @brief returns the lazy interaction interaction graph flag
	bool lazy_ig() const override;


	/// @brief Activate the DoubleLazyInteractionGraph, which is particularly useful
	/// for multistate design, when memory and time are both limiting.  Overriden by LinMemIG.
	void or_double_lazy_ig( bool setting ) override;
	/// @brief Returns the double-lazy interaction graph flag
	bool double_lazy_ig() const override;
	/// @brief Set the memory limit, in bytes, for the storage that the DLIG should be allowed
	/// to spend on representing rotamer pair energies.  The DLIG will start deallocating
	/// rotamer pair energy blocks if it exceeds this limit.  This limit is occasionally breached by
	/// as much memory as is required to store a single amino-acid-pair submatrix block; the limit
	/// does not apply to the entirety of the packer or even the entirety of this interaction graph,
	/// merely to the amount of space for rotamer pair energies.  Remember, rotamers are expensive, too!
	/// The default value of "0" signifies an unrestricted memory limit -- 0 may not be set through
	/// this function.  The value may increase once from 0 and then may only decrease from there.
	void decrease_double_lazy_ig_memlimit( Size nbytes_for_rpes ) override;
	/// @brief the memory limit, in bytes, for the double-lazy interaction graph.  A value of 0
	/// signifies an unrestricted limit.
	Size double_lazy_ig_memlimit() const override;


	/// @brief if setting == true, turns on MultiCoolAnnealer -- so long as rotamer couplings are not
	///also turned on.
	void or_multi_cool_annealer( bool setting ) override;
	/// @brief use MultiCoolAnnealer?
	bool multi_cool_annealer() const override;

	/// @brief Increases the history size for the MultiCoolAnnealer if setting is larger than the
	///existing setting.
	void increase_multi_cool_annealer_history_size( Size setting ) override;
	/// @brief returns the requested size for the MultiCoolAnnealer
	Size multi_cool_annealer_history_size() const override;

	bool smart_annealer() const override { return smart_annealer_; }
	std::string smart_annealer_model() const override { return smart_annealer_model_; }
	core::Real smart_annealer_cutoff() const override { return smart_annealer_cutoff_; }
	bool smart_annealer_pick_again() const override { return smart_annealer_pick_again_; }
	bool smart_annealer_disable_during_quench() const override {
		return smart_annealer_disable_during_quench_;
	}

	void set_smart_annealer( bool setting ) override { smart_annealer_ = setting; }
	void set_smart_annealer_model( std::string const & setting ) override { smart_annealer_model_ = setting; }
	void set_smart_annealer_cutoff( core::Real setting ) override { smart_annealer_cutoff_ = setting; }
	void set_smart_annealer_pick_again( bool setting ) override { smart_annealer_pick_again_ = setting; }
	void set_smart_annealer_disable_during_quench( bool setting ) override { smart_annealer_disable_during_quench_ = setting; }

	void show( std::ostream & out ) const override;
	void show() const override;
	void show_residue_task( std::ostream & out, Size resid ) const override;
	void show_residue_task( Size resid ) const override;
	void show_all_residue_tasks( std::ostream & out ) const override;
	void show_all_residue_tasks() const override;

	/// @brief Has this PackerTask been initialized with a PackerPalette?
	/// @details PackerTasks must be initialized with PackerPalettes before being modified with TaskOperations.  The TaskFactory
	/// will initialize the PackerTask with a DefaultPackerPalette if no custom PackerPalette is provided.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_initialized() const override;

	/// @brief read command line options (but not resfile) to set the state of the PackerTask, NOT IN CONSTRUCTOR
	PackerTask &
	initialize_from_command_line() override;

	/// @brief read only the command line options for extra rotamer building;
	PackerTask &
	initialize_extra_rotamer_flags_from_command_line() override;

	PackerTask &
	initialize_from_options( utility::options::OptionCollection const & options ) override;

	PackerTask &
	initialize_extra_rotamer_flags_from_options( utility::options::OptionCollection const & options ) override;

	/// @brief turn off packing for residues passed false; can't turn on packing
	PackerTask &
	restrict_to_residues( utility::vector1< bool > const & residues_allowed_to_be_packed ) override;

	/// @brief turn off designing (sequence changing) all residues
	PackerTask &
	restrict_to_repacking() override;

	/// @brief const accessor for underlying ResidueLevelTask object
	ResidueLevelTask const &
	residue_task( Size resid ) const override;

	/// @brief nonconst access to underlying ResidueLevelTask object
	ResidueLevelTask &
	nonconst_residue_task( Size resid ) override;

	utility::vector1< bool >
	repacking_residues() const override;

	utility::vector1< bool >
	designing_residues() const override;

	bool
	keep_sequence_symmetry() const override {
		return keep_sequence_symmetry_;
	}

	void
	keep_sequence_symmetry( bool const setting ) override {
		keep_sequence_symmetry_ = setting;
	}

	/// @brief is there at RotamerCouplings object to worry about? (for DNA GC AT pairing, etc)
	bool
	rotamer_couplings_exist() const override;

	/// @brief const accessor for the RotamerCouplings object
	RotamerCouplingsCOP
	rotamer_couplings() const override;

	/// @brief setter for the RotamerCouplings object
	void
	rotamer_couplings( RotamerCouplingsCOP setting ) override;

	/// @brief is there at RotamerLinks object to worry about? (for repeat linking
	//of equivalent residues, etc)
	bool
	rotamer_links_exist() const override;

	/// @brief const accessor for the RotamerLinks object
	RotamerLinksCOP
	rotamer_links() const override;

	/// @brief setter for the RotamerLinks object
	void
	rotamer_links( RotamerLinksCOP setting ) override;


	/// @brief accesor for residue residue weight map
	IGEdgeReweightContainerCOP
	IGEdgeReweights() const override;

	IGEdgeReweightContainerOP
	set_IGEdgeReweights() override;

	void
	append_rotamer_operation(
		rotamer_set::RotamerOperationOP rotop
	) override;

	void
	append_rotamerset_operation(
		rotamer_set::RotamerSetOperationOP rotsetop
	) override;

	rotamer_set::RotSetsOperationListIterator
	rotamer_sets_operation_begin() const override;

	rotamer_set::RotSetsOperationListIterator
	rotamer_sets_operation_end() const override;

	void
	append_rotamersets_operation(
		rotamer_set::RotamerSetsOperationOP rotsetsop
	) override;

	//@brief return a resfile string that can recreate the packer task
	std::string
	task_string( pose::Pose const & pose ) const override;


	//some get-set functions for annealer options
	void
	low_temp( Real const & low_temp ) override;

	void
	high_temp( Real const & high_temp ) override;

	void
	disallow_quench( bool const & disallow_quench ) override;

	Real
	low_temp() const override;

	Real
	high_temp() const override;

	bool
	disallow_quench() const override;

	void remap_residue_level_tasks(
		core::id::SequenceMappingCOP seqmap,
		core::pose::Pose const & pose
	) override;

	//////////////////////// dangerous update functions ///////////////////////

	virtual
	void
	update_residue_union(
		Size resid,
		ResidueLevelTask const & res_task_in
	);

	virtual
	void
	update_residue_intersection(
		Size resid,
		ResidueLevelTask const & res_task_in
	);

	virtual
	void
	update_residue_commutative(
		Size resid,
		ResidueLevelTask const & res_task_in
	);

	void
	update_commutative(
		PackerTask const & tark_in
	) override;

	void
	request_symmetrize_by_intersection() override;

	void
	request_symmetrize_by_union() override;

	bool
	symmetrize_by_union() const override;

	bool
	symmetrize_by_intersection() const override;

	/////////////////////////////////////
	/// For multithreading            ///
	/////////////////////////////////////

	/// @brief How many threads should the packer request for interaction graph precomputation?
	/// @details Must be implemented by derived class.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	core::Size ig_threads_to_request() const override;

	/// @brief Limit the interaction graph setup threads.
	/// @details If the current ig_threads_to_request_ is zero, setting replaces it.  If setting is zero, this
	/// does nothing.  Otherwise, setting replaces ig_threads_to_request_ if and only if it is less than
	/// ig_threads_to_request_.  This preserves commutativity.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitue.org)
	void limit_ig_setup_threads( core::Size const setting ) override;

	/////////////////////////////////////
	/// For use only inside the packer //
	/////////////////////////////////////

	/// @brief turn off packing at all positions
	/// This is only to be used to control packing from within the packing
	/// functions rotamer_trials and rtmin.  Otherwise, this function should
	/// not be used.  This function is not commutative and therefor using it
	/// outside of rotamer_trials and rtmin would lead to disasterous behaviors
	void temporarily_fix_everything() override;

	/// @brief reset packer mutability arbitrarily for a given residue
	/// This is only to be used to control packing from within the packing
	/// functions rotamer_trials and rtmin.  Otherwise, this function should
	/// not be used. This function is not commutative and therefor using it
	/// outside of rotamer_trials and rtmin would lead to disasterous behaviors
	void temporarily_set_pack_residue( int resid, bool setting ) override;

private: // private methods
	void
	update_n_to_be_packed() const;

	virtual
	PackerTask &
	operator=(PackerTask const &);


private:

	/// @brief Has this PackerTask been initialized with a PackerPalette?
	/// @details PackerTasks must be initialized with PackerPalettes before being modified with TaskOperations.  The TaskFactory
	/// will initialize the PackerTask with a DefaultPackerPalette if no custom PackerPalette is provided.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_initialized_ = false;

	Size nres_;

	/// @details superficial on/off switch for repacking residues
	/// produces no impact on contents of residue_tasks_ vector.
	/// This is only to be used to control packing from within the packing
	/// functions rotamer_trials and rtmin.  Otherwise, this vector should
	/// not be used.
	utility::vector1< bool > pack_residue_;

	/// @details vector of ResidueLevelTasks
	utility::vector1< ResidueLevelTask_ > residue_tasks_;

	mutable Size n_to_be_packed_;
	mutable bool n_to_be_packed_up_to_date_ = true;

	bool precompute_ig_ = false;

	/// @details linmem_ig overrides PDInteractionGraph
	bool linmem_ig_ = false;
	bool linmem_ig_history_size_at_default_ = true;
	Size linmem_ig_history_size_ = 10;
	/// @details linmem_ig overrides LazyIG
	bool lazy_ig_ = false;

	/// @details linmem_ig overrides DoubleLazyIG
	bool double_lazy_ig_ = false;
	Size dlig_mem_limit_ = 0;

	bool multi_cool_annealer_ = false;
	Size mca_history_size_ = 1;

	bool smart_annealer_ = false;
	std::string smart_annealer_model_ = "generation2";
	core::Real smart_annealer_cutoff_ = 0.25;
	bool smart_annealer_pick_again_ = true;
	bool smart_annealer_disable_during_quench_ = true;

	bool keep_sequence_symmetry_ = false;

	/// @brief keep track: are we optimizing H at all positions?
	/// don't use linmem_ig_ if so.
	bool optimize_H_ = false;

	bool bump_check_ = true;
	Real max_rotbump_energy_ = 5.0;

	// pbhack temporary -- usually null pointer
	RotamerCouplingsCOP rotamer_couplings_;

	// psh hack, just like the pbhack above
	RotamerLinksCOP rotamer_links_;

	// rhiju -- some options that need to be sent to the annealer
	Real low_temp_ = -1.0;
	Real high_temp_ = -1.0;
	bool disallow_quench_ = false;

	/// @details holds specific residue residue weights to be used in packing
	IGEdgeReweightContainerOP IG_edge_reweights_;

	rotamer_set::RotSetsOperationList rotsetsops_;

	// sheffler
	PackerTaskSymmetryStatus symmetry_status_ = NO_SYMMETRIZATION_REQUEST;

	/// @brief The PackerPalette, which tells the PackerTask what the default set of residue types is in the absence of TaskOperations.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::pack::palette::PackerPaletteCOP packer_palette_;

	/// @brief The number of threads to request for interaction graph precomputation.
	/// @details The default, 0, means to request all available threads (multithreaded build) or
	/// 1 thread (non-threaded build).  In the non-threaded build, only 1 thread can be requested.
	core::Size ig_threads_to_request_ = 0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

//NOTE: parse_resfile is now an independent function in ResfileReader.hh, not a member function of the PackerTask hierarchy

} //namespace task
} //namespace pack
} //namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_PackerTask_ )
#endif // SERIALIZATION


#endif
