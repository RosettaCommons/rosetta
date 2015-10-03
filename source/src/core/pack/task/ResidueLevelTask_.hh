// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResidueLevelTask_.hh
/// @brief  Implementation class for task class to describe packer's behavior header
/// Almost all of rosetta needs to use packer tasks, but very little of rosetta needs
/// to see how it behaves internally.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Steven Lewis

#ifndef INCLUDED_core_pack_task_ResidueLevelTask__hh
#define INCLUDED_core_pack_task_ResidueLevelTask__hh

// Unit Headers
#include <core/pack/task/ResidueLevelTask.hh>

// Package Headers
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


/// @brief Residue-level task class
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
	virtual void allow_aa( chemical::AA const & aa );

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
	virtual std::string command_string( ) const;


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
	std::string task_mode() const;

	rna::RNA_ResidueLevelTask const & rna_task() const;

	rna::RNA_ResidueLevelTask & nonconst_rna_task();

	void
	do_restrict_absent_canonical_aas( utility::vector1<bool> const & allowed_aas );

	/// @brief legacy code block to handle packing at residues with custom variants.
	void
	setup_allowed_protein_residue_types_LEGACY( chemical::ResidueType const & match_residue_type );

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

} //namespace task
} //namespace pack
} //namespace core

#endif
