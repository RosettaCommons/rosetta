// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

// STL Headers

#include <iosfwd>
#include <iostream>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {

/// @brief Residue-level task class; contained within PackerTask
class ResidueLevelTask_ : public ResidueLevelTask
{
public:
	/// @brief constructor; requires a Residue object
	ResidueLevelTask_(
		conformation::Residue const & original_residue
	);

	ResidueLevelTask_();

	/// @brief dtor
	virtual ~ResidueLevelTask_();

	/// @brief returns the extra chi sampling level
	virtual
	ExtraRotSample
	extrachi_sample_level(
		bool buried,
		int chi,
		chemical::ResidueType const & concrete_residue
	) const;

	/// @brief initialize options from command line flags
	virtual void initialize_from_command_line();

	/// @brief Initialize only the extra rotamer building flags from the command line;
	/// invoked by intialize_from_command_line();
	virtual void initialize_extra_rotamer_flags_from_command_line();

	/// @brief include the pre-existing rotamer while packing
	virtual void or_include_current( bool include_current );
	/// @brief is the pre-existing rotamer specifically allowed while packing?
	virtual bool include_current() const;

	// access to any flagged protocol-level behaviors for this residue
	virtual void add_behavior( std::string const & behavior );
	virtual bool has_behavior( std::string const & behavior ) const;
	virtual bool has_behavior() const;

	virtual void target_type( chemical::ResidueTypeCOP type );
	virtual void target_type( chemical::AA aa );
	virtual void target_type( std::string name );

	/// @brief include adducts at this residue
	virtual void or_adducts( bool setting );
	virtual bool adducts() const;

	/// @brief activate ex1 when passed true; do nothing otherwise
	virtual void or_ex1( bool ex1 );
	/// @brief activate ex2 when passed true; do nothing otherwise
	virtual void or_ex2( bool ex2 );
	/// @brief activate ex3 when passed true; do nothing otherwise
	virtual void or_ex3( bool ex3 );
	/// @brief activate ex4 when passed true; do nothing otherwise
	virtual void or_ex4( bool ex4 );

	/// @brief increase ex1 sample level; do nothing if not an increase
	virtual void or_ex1_sample_level( ExtraRotSample ex1_sample_level );
	/// @brief increase ex2 sample level; do nothing if not an increase
	virtual void or_ex2_sample_level( ExtraRotSample ex2_sample_level );
	/// @brief increase ex3 sample level; do nothing if not an increase
	virtual void or_ex3_sample_level( ExtraRotSample ex3_sample_level );
	/// @brief increase ex4 sample level; do nothing if not an increase
	virtual void or_ex4_sample_level( ExtraRotSample ex4_sample_level );

	/// @brief activate ex1 for aromatics when passed true; do nothing otherwise
	virtual void or_ex1aro( bool ex1aro );
	/// @brief activate ex2 for aromatics when passed true; do nothing otherwise
	virtual void or_ex2aro( bool ex2aro_only );
	/// @brief activate ex1 for exposed aromatics when passed true; do nothing otherwise
	virtual void or_ex1aro_exposed( bool ex1aro_exposed );
	/// @brief activate ex2 for exposed aromatics when passed true; do nothing otherwise
	virtual void or_ex2aro_exposed( bool ex2aro_exposed );

	/// @brief increase ex1aro sample level; do nothing if not an increase
	virtual void or_ex1aro_sample_level( ExtraRotSample ex1aro_sample_level );
	/// @brief increase ex2aro sample level; do nothing if not an increase
	virtual void or_ex2aro_sample_level( ExtraRotSample ex2aro_only_sample_level );
	/// @brief increase ex1aro_exposed sample level; do nothing if not an increase
	virtual void or_ex1aro_exposed_sample_level( ExtraRotSample ex1aro_exposed_sample_level );
	/// @brief increase ex2aro_exposed sample level; do nothing if not an increase
	virtual void or_ex2aro_exposed_sample_level( ExtraRotSample ex2aro_exposed_sample_level );

	virtual void or_operate_on_ex1( bool operate );
	virtual void or_operate_on_ex2( bool operate );
	virtual void or_operate_on_ex3( bool operate );
	virtual void or_operate_on_ex4( bool operate );

	virtual void or_exdna_sample_level( ExtraRotSample exdna_sample_level );

	virtual void or_optimize_h( bool setting );
	virtual bool optimize_h() const;
	virtual void or_preserve_c_beta( bool setting );
	virtual bool preserve_c_beta() const;
	virtual void or_flip_HNQ( bool setting );
	virtual bool flip_HNQ() const;
	virtual void or_fix_his_tautomer( bool setting );
	virtual bool fix_his_tautomer() const;

	virtual void or_include_virtual_side_chain( bool setting );
	virtual bool include_virtual_side_chain() const;

	/// @brief sample proton chi.
	virtual void sample_proton_chi( bool setting );

	/// @brief sample proton chi.
	virtual bool sample_proton_chi() const;

	virtual bool ex1() const;
	virtual bool ex2() const;
	virtual bool ex3() const;
	virtual bool ex4() const;

	virtual ExtraRotSample ex1_sample_level() const;
	virtual ExtraRotSample ex2_sample_level() const;
	virtual ExtraRotSample ex3_sample_level() const;
	virtual ExtraRotSample ex4_sample_level() const;

	virtual bool ex1aro() const;
	virtual bool ex2aro() const;
	virtual bool ex1aro_exposed() const;
	virtual bool ex2aro_exposed() const;

	virtual ExtraRotSample ex1aro_sample_level() const;
	virtual ExtraRotSample ex2aro_sample_level() const;
	virtual ExtraRotSample ex1aro_exposed_sample_level() const;
	virtual ExtraRotSample ex2aro_exposed_sample_level() const;

	virtual ExtraRotSample exdna_sample_level() const;

	virtual bool operate_on_ex1() const;
	virtual bool operate_on_ex2() const;
	virtual bool operate_on_ex3() const;
	virtual bool operate_on_ex4() const;

	/// @brief lower extrachi_cutoff to given value; do nothing if not a decrease
	virtual void and_extrachi_cutoff( Size num_neighbors_to_be_called_buried );

	/// @brief get function for extrachi_cutoff
	virtual Size extrachi_cutoff() const;

	/// @brief remove all ResidueTypes from the list of allowed residue types, preventing repacking
	virtual void prevent_repacking();

	/// @brief disables designing to residues not in the passed list
	virtual void restrict_absent_canonical_aas( utility::vector1< bool > const &);

	/// @brief disables designing to residues not in the passed list--and specifies the resfile command that made this list
	virtual void restrict_absent_canonical_aas( utility::vector1< bool > const & allowed_aas,
																							std::string const & mode );

	//@brief Same behavior restrict_absent_canonical_aas except that it always allows the native aa at a position even if it is not included in the allowed residues
	virtual void restrict_nonnative_canonical_aas( utility::vector1< bool > const & allowed_aas);


	/// @brief disables designing to nucleic acid residues not in the passed list
	virtual void restrict_absent_nas( utility::vector1< chemical::AA > const & keep_nas );

	/// @brief only let this residue repack -- prevent redesign
	virtual void restrict_to_repacking();

	/// @brief
	virtual bool is_original_type( chemical::ResidueTypeCOP type ) const;

	/// @brief
	virtual chemical::ResidueTypeSet const & get_original_residue_set() const;

	/// @brief
	virtual chemical::AA const & get_original_residue() const;

	/// @brief explicitly allow a NCAA
	virtual void allow_noncanonical_aa(
		std::string const & interchangeability_group,
		chemical::ResidueTypeSet const & residue_set
	);

	/// @brief explicitly allow a NCAA; assumes same ResidueTypeSet as original_residue_type_
	virtual void allow_noncanonical_aa( std::string const & aaname );

	/// @brief explicitly allow a NCAA; assumes same ResidueTypeSet as original_residue_type_
	virtual void allow_noncanonical_aa( chemical::AA aa );

	/// @brief expliciitly disallow all NCAAs
	virtual void disallow_noncanonical_aas();

	/// @brief explicitly allow a canonical AA
	virtual	void allow_aa( chemical::AA const & aa );

	virtual ResidueTypeCOPList const & allowed_residue_types() const;
	/// @brief returns iterator to beginning of allowed residue types list (traversal only)
	virtual ResidueTypeCOPListConstIter allowed_residue_types_begin() const;
	/// @brief returns iterator to end of allowed residue types list (traversal only)
	virtual ResidueTypeCOPListConstIter allowed_residue_types_end() const;

	virtual chemical::ResidueTypeCOP target_type() const;

	virtual void print_allowed_types( std::ostream & os ) const;

	/// @brief is this residue up for design (variable sequence)?
	virtual bool being_designed() const;
	/// @brief is this residue modififable at all by the packer?
	virtual bool being_packed() const;

	virtual
	rotamer_set::RotamerOperations const &
	rotamer_operations() const;

	virtual
	void
	append_rotamer_operation(
		rotamer_set::RotamerOperationOP rotop
	);

	virtual rotamer_set::RotSetOperationListIterator
	rotamer_set_operation_begin() const;

	virtual rotamer_set::RotSetOperationListIterator
	rotamer_set_operation_end() const;

	virtual
	void
	append_rotamerset_operation(
		rotamer_set::RotamerSetOperationOP rotsetop
	);


	/// @brief create a string the resfile format of all the commands applied to this residue level task
	virtual	std::string	command_string( ) const;


	///////////////////////////// dangerous update functions

	virtual
	void
	update_union(
		ResidueLevelTask const & res_task_in
	);

	virtual
	void
	update_intersection(
		ResidueLevelTask const & res_task_in
	);

	virtual
	void
	update_commutative(
		ResidueLevelTask const & res_task_in
	);


private: // private methods
	/// @brief private: bookkeeping for ex1
	void refresh_ex1_sample_levels();
	/// @brief private: bookkeeping for ex2
	void refresh_ex2_sample_levels();
	/// @brief private: bookkeeping for ex3
	void refresh_ex3_sample_levels();
	/// @brief private: bookkeeping for ex4
	void refresh_ex4_sample_levels();

	/// @brief private: bookkeeping for whether sequence can change
	void determine_if_designing();
	/// @brief private: bookkeeping for whether rotamer can change
	void determine_if_repacking();


	/// @brief private: return the EX command for the packer task
	std::string
	get_ex_flags( Size chiid,
								Size const exaro_sample_level,
								Size const ex_sample_level) const;

	/// @brief private: return the task mode that can be used to recreate
	///the task.  If the residue level task was made with POLAR it
	///should return the string "POLAR".
	std::string	task_mode() const;

	rna::RNA_ResidueLevelTask const & rna_task() const;

	rna::RNA_ResidueLevelTask & nonconst_rna_task();

	void
	do_restrict_absent_canonical_aas( utility::vector1<bool> const & allowed_aas );

private:
	/// @details is the pre-existing rotamer included for the packer to choose?
	bool include_current_;
	/// here we store any flagged protocol-level behaviors for this residue
	utility::vector1< std::string > behaviors_;
	/// @brief include adducts at this residue
	bool adducts_;
	/// @details std::list of ResidueTypeCOP objects - these are only residue types allowed at position
	ResidueTypeCOPList allowed_residue_types_;
	chemical::ResidueTypeCOP original_residue_type_; //record this on construction
	/// @details a member of the allowed types that can optionally respresent a target state
	chemical::ResidueTypeCOP target_residue_type_;

	/// @details can the sequence change?  true implies repacking_ is true as well
	bool designing_;
	/// @details can the residue change its conformation?
	bool repacking_;

	/// @details are we keeping the coordinates of the heavy atoms fixed while
	///sampling extra hydrogen placements?
	bool optimize_H_mode_;

	/// @details are c-beta positions preserved during rotamer building
	bool preserve_c_beta_;

	/// @details are we also considering heavy atom rearrangement for the crystallographically
	///symmetric seeming residues: histadine, asparagine and glutamine? Implies optimize_H_mode_
	bool flip_HNQ_;

	/// @details has this histidine tautomer been fixed?  This value is kept for bookkeeping; the tautomer's fixation is effected by removing the other tautomer from the ResidueTypeCOPList.
	bool fix_his_tautomer_;

	bool include_virtual_side_chain_;

	/// @details if this is true, this residue will be treated as part of the background.
	///a disabling takes precedence over any ResidueType additions.
	bool disabled_;
	/// @details if this is true, this residue will only be allowed to repack;
	///ResidueTypes that do not match the original_residue_type cannot be added.
	///Design disabling takes precedence over SOMETHING, ANDREW TELL ME WHAT IT IS
	bool design_disabled_;

	bool sample_proton_chi_;

	bool ex1_;
	bool ex2_;
	bool ex3_;
	bool ex4_;
	bool ex1aro_;
	bool ex2aro_;
	bool ex1aro_exposed_;
	bool ex2aro_exposed_;
	ExtraRotSample ex1_sample_level_;
	ExtraRotSample ex2_sample_level_;
	ExtraRotSample ex3_sample_level_;
	ExtraRotSample ex4_sample_level_;
	ExtraRotSample ex1aro_sample_level_;
	ExtraRotSample ex2aro_sample_level_;
	ExtraRotSample ex1aro_exposed_sample_level_;
	ExtraRotSample ex2aro_exposed_sample_level_;
	ExtraRotSample exdna_sample_level_;
	Size extrachi_cutoff_;

	bool operate_on_ex1_;
	bool operate_on_ex2_;
	bool operate_on_ex3_;
	bool operate_on_ex4_;

	rotamer_set::RotamerOperations rotamer_operations_;
	rotamer_set::RotSetOperationList rotsetops_;

	std::vector<std::string> mode_tokens_;

	rna::RNA_ResidueLevelTaskOP rna_task_;

};

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
	virtual void clean_residue_task( conformation::Residue const & original_residue, Size const seqpos);

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

};

//NOTE: parse_resfile is now an independent function in ResfileReader.hh, not a member function of the PackerTask hierarchy

} //namespace task
} //namespace pack
} //namespace core

#endif
