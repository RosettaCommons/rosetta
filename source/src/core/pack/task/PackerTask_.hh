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
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// STL Headers
#include <iosfwd>
#include <iostream>


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
	PackerTask_();
	PackerTask_( pose::Pose const & pose );

	/// @brief dtor
	virtual ~PackerTask_();

	/// @brief copy method
	PackerTaskOP clone() const;

	/// @brief replace a given residue task with a brand new one NOTE: This should be the only way to break commutativity!!!!
	virtual void clean_residue_task( conformation::Residue const & original_residue, Size const seqpos, core::pose::Pose const & pose);

	/// @brief number of residues in the input pose, for convienience (PackerTask does not handle variable length)
	virtual Size total_residue() const;

	/// @brief get function: can this position have a rotamer change?
	virtual bool pack_residue( int resid ) const;
	/// @brief alias for above
	virtual bool being_packed( Size resid ) const { return pack_residue( resid ); }

	/// @brief get function: how many positions can have rotamer changes?
	virtual Size num_to_be_packed() const;

	/// @brief get function: can this position have a sequence change?
	virtual bool design_residue( int resid ) const;
	/// @brief alias for above
	virtual bool being_designed( Size resid ) const { return design_residue( resid ); }
	/// @brief get function: can any positions have a sequence change?
	virtual bool design_any() const;

	/// @brief set function: bump_check is activated for pack-rotamers' screening of rotamers that
	/// collide with the background.  Bump-check is not used during rotamer trials, since it is
	/// nearly as expensive as rotamer trials itself.  Energy methods may opt in to the bump check
	/// process.  The "standard" behavior is for bump-check to include only the fa_atr and fa_rep
	/// terms.
	virtual void set_bump_check( bool setting );
	/// @brief get function: has bump_check been requested?
	virtual bool bump_check() const;
	/// @brief Decrease the max_rotbump_energy threshold above which rotamers are rejected.
	virtual void and_max_rotbump_energy( Real setting );
	/// @brief get function: what is the energy threshold above which rotamers should be rejected
	virtual Real max_rotbump_energy() const;

	/// @brief for all positions, turn on include_current if false, do nothing if already true
	virtual void or_include_current( bool setting );
	/// @brief for one position, turn on include_current if false, do nothing if already true
	virtual void or_include_current( bool setting, Size resid );
	/// @brief get function: what is include_current for this residue?
	virtual bool include_current( Size resid ) const;

	virtual void add_behavior( std::string const & behavior );
	virtual void add_behavior( std::string const & behavior, Size resid );
	virtual bool has_behavior( std::string const & behavior, Size resid ) const;
	virtual bool has_behavior( Size resid ) const;

	/// @brief return the targeted type (may be null pointer)
	virtual chemical::ResidueTypeCOP target_type( Size resid ) const;

	/// @brief for all positions, disable adducts if false, do nothing if true
	virtual void or_adducts( bool setting );
	/// @brief for one position, disable adducts if false, do nothing if true
	virtual void or_adducts( bool setting, Size resid );
	/// @brief include adducts at this residue
	virtual bool adducts( Size resid ) const;

	/// @brief if setting == true, turns on optimize_H_mode for all residues
	virtual void or_optimize_h_mode( bool setting );

	/// @brief if setting == true, preserves c-beta during rotamer building for all residues
	virtual void or_preserve_c_beta( bool setting );

	/// @brief if setting == true, turns on optimize_H_mode and flip_HNQ for all residues
	virtual void or_flip_HNQ( bool setting );

	/// @brief if setting == true, fix his tautomer state for defined residues during repacking or optimizeH mode
	virtual void or_fix_his_tautomer( utility::vector1<int> const & positions, bool setting );

	/// @brief if setting == true, turns on linear-memory interaction graph usage
	virtual void or_linmem_ig( bool setting );
	/// @brief if setting == false, turns off linear-memory interaction graph usage
	virtual void and_linmem_ig( bool setting );
	/// @brief returns the linear-memory interaction graph flag
	virtual bool linmem_ig() const;

	/// @brief Set the linear memory interaction graph's recent history size.  Default is 10.  It can be set
	/// once to any value larger than 10, but any subsequent setting will only alter the size if it decreases
	/// it.
	virtual void decrease_linmem_ig_history_size( Size setting );

	/// @brief Return the linear memory interaction graph's recent history size.
	virtual Size linmem_ig_history_size() const;

	/// @brief  if setting == true, turns on lazy interaction graph usage
	/// NOTE: the linear memory interaction graph takes precedence over the LazyIG when
	/// the InteractionGraphFactory examines the PackerTask.
	virtual void or_lazy_ig( bool setting );
	/// @brief returns the lazy interaction interaction graph flag
	virtual bool lazy_ig() const;


	/// @brief Activate the DoubleLazyInteractionGraph, which is particularly useful
	/// for multistate design, when memory and time are both limiting.  Overriden by LinMemIG.
	virtual void or_double_lazy_ig( bool setting );
	/// @brief Returns the double-lazy interaction graph flag
	virtual bool double_lazy_ig() const;
	/// @brief Set the memory limit, in bytes, for the storage that the DLIG should be allowed
	/// to spend on representing rotamer pair energies.  The DLIG will start deallocating
	/// rotamer pair energy blocks if it exceeds this limit.  This limit is occasionally breached by
	/// as much memory as is required to store a single amino-acid-pair submatrix block; the limit
	/// does not apply to the entirety of the packer or even the entirety of this interaction graph,
	/// merely to the amount of space for rotamer pair energies.  Remember, rotamers are expensive, too!
	/// The default value of "0" signifies an unrestricted memory limit -- 0 may not be set through
	/// this function.  The value may increase once from 0 and then may only decrease from there.
	virtual void decrease_double_lazy_ig_memlimit( Size nbytes_for_rpes );
	/// @brief the memory limit, in bytes, for the double-lazy interaction graph.  A value of 0
	/// signifies an unrestricted limit.
	virtual Size double_lazy_ig_memlimit() const;


	/// @brief if setting == true, turns on MultiCoolAnnealer -- so long as rotamer couplings are not
	///also turned on.
	virtual void or_multi_cool_annealer( bool setting );
	/// @brief use MultiCoolAnnealer?
	virtual bool multi_cool_annealer() const;

	/// @brief Increases the history size for the MultiCoolAnnealer if setting is larger than the
	///existing setting.
	virtual void increase_multi_cool_annealer_history_size( Size setting );
	/// @brief returns the requested size for the MultiCoolAnnealer
	virtual Size multi_cool_annealer_history_size() const;

	virtual void show( std::ostream & out ) const;
	virtual void show() const;
	virtual void show_residue_task( std::ostream & out, Size resid ) const;
	virtual void show_residue_task( Size resid ) const;
	virtual void show_all_residue_tasks( std::ostream & out ) const;
	virtual void show_all_residue_tasks() const;

	/// @brief read command line options (but not resfile) to set the state of the PackerTask, NOT IN CONSTRUCTOR
	virtual
	PackerTask &
	initialize_from_command_line();

	/// @brief read only the command line options for extra rotamer building;
	virtual PackerTask &
	initialize_extra_rotamer_flags_from_command_line();

	virtual
	PackerTask &
	initialize_from_options( utility::options::OptionCollection const & options );

	virtual PackerTask &
	initialize_extra_rotamer_flags_from_options( utility::options::OptionCollection const & options );

	/// @brief turn off packing for residues passed false; can't turn on packing
	virtual
	PackerTask &
	restrict_to_residues( utility::vector1< bool > const & residues_allowed_to_be_packed );

	/// @brief turn off designing (sequence changing) all residues
	virtual
	PackerTask &
	restrict_to_repacking();

	/// @brief const accessor for underlying ResidueLevelTask object
	virtual
	ResidueLevelTask const &
	residue_task( Size resid ) const;

	/// @brief nonconst access to underlying ResidueLevelTask object
	virtual
	ResidueLevelTask &
	nonconst_residue_task( Size resid );

	virtual
	utility::vector1< bool >
	repacking_residues() const;

	virtual
	utility::vector1< bool >
	designing_residues() const;

	/// @brief is there at RotamerCouplings object to worry about? (for DNA GC AT pairing, etc)
	virtual
	bool
	rotamer_couplings_exist() const;

	/// @brief const accessor for the RotamerCouplings object
	virtual
	RotamerCouplingsCOP
	rotamer_couplings() const;

	/// @brief setter for the RotamerCouplings object
	virtual
	void
	rotamer_couplings( RotamerCouplingsCOP setting );

	/// @brief is there at RotamerLinks object to worry about? (for repeat linking
	//of equivalent residues, etc)
	virtual
	bool
	rotamer_links_exist() const;

	/// @brief const accessor for the RotamerLinks object
	virtual
	RotamerLinksCOP
	rotamer_links() const;

	/// @brief setter for the RotamerLinks object
	virtual
	void
	rotamer_links( RotamerLinksCOP setting );


	/// @brief accesor for residue residue weight map
	virtual
	IGEdgeReweightContainerCOP
	IGEdgeReweights() const;

	virtual
	IGEdgeReweightContainerOP
	set_IGEdgeReweights();

	virtual
	void
	append_rotamer_operation(
		rotamer_set::RotamerOperationOP rotop
	);

	virtual
	void
	append_rotamerset_operation(
		rotamer_set::RotamerSetOperationOP rotsetop
	);

	//@brief return a resfile string that can recreate the packer task
	virtual
	std::string
	task_string( pose::Pose const & pose ) const;


	//some get-set functions for annealer options
	void
	low_temp( Real const & low_temp );

	void
	high_temp( Real const & high_temp );

	void
	disallow_quench( bool const & disallow_quench );

	Real
	low_temp() const;

	Real
	high_temp() const;

	bool
	disallow_quench() const;

	virtual
	void remap_residue_level_tasks(
		core::id::SequenceMappingCOP seqmap,
		core::pose::Pose const & pose
	);

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

	virtual
	void
	update_commutative(
		PackerTask const & tark_in
	);

	virtual
	void
	request_symmetrize_by_intersection();

	virtual
	void
	request_symmetrize_by_union();

	virtual
	bool
	symmetrize_by_union() const;

	virtual
	bool
	symmetrize_by_intersection() const;

	/////////////////////////////////////
	/// For use only inside the packer //
	/////////////////////////////////////

	/// @brief turn off packing at all positions
	/// This is only to be used to control packing from within the packing
	/// functions rotamer_trials and rtmin.  Otherwise, this function should
	/// not be used.  This function is not commutative and therefor using it
	/// outside of rotamer_trials and rtmin would lead to disasterous behaviors
	virtual void temporarily_fix_everything();

	/// @brief reset packer mutability arbitrarily for a given residue
	/// This is only to be used to control packing from within the packing
	/// functions rotamer_trials and rtmin.  Otherwise, this function should
	/// not be used. This function is not commutative and therefor using it
	/// outside of rotamer_trials and rtmin would lead to disasterous behaviors
	virtual void temporarily_set_pack_residue( int resid, bool setting );

private: // private methods
	void
	update_n_to_be_packed() const;

	virtual
	PackerTask &
	operator=(PackerTask const &);


private:

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
	mutable bool n_to_be_packed_up_to_date_;

	/// @details linmem_ig overrides PDInteractionGraph
	bool linmem_ig_;
	bool linmem_ig_history_size_at_default_;
	Size linmem_ig_history_size_;
	/// @details linmem_ig overrides LazyIG
	bool lazy_ig_;

	/// @details linmem_ig overrides DoubleLazyIG
	bool double_lazy_ig_;
	Size dlig_mem_limit_;

	bool multi_cool_annealer_;
	Size mca_history_size_;

	/// @brief keep track: are we optimizing H at all positions?
	/// don't use linmem_ig_ if so.
	bool optimize_H_;

	bool bump_check_;
	Real max_rotbump_energy_;

	// pbhack temporary -- usually null pointer
	RotamerCouplingsCOP rotamer_couplings_;

	// psh hack, just like the pbhack above
	RotamerLinksCOP rotamer_links_;

	// rhiju -- some options that need to be sent to the annealer
	Real low_temp_;
	Real high_temp_;
	bool disallow_quench_;

	/// @details holds specific residue residue weights to be used in packing
	IGEdgeReweightContainerOP IG_edge_reweights_;

	// sheffler
	PackerTaskSymmetryStatus symmetry_status_;

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
