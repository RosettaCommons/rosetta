// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LoopHashSampler.hh
/// @brief
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loophash_LoopHashSampler_hh
#define INCLUDED_protocols_loophash_LoopHashSampler_hh

#include <protocols/loophash/LoopHashSampler.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LocalInserter.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <string>
#include <vector>

#include <utility/vector1.hh>


namespace protocols {
namespace loophash {


/// @brief Create candidate structures where some residues have been sampled by
/// loophash.
class LoopHashSampler : public utility::pointer::ReferenceCount  {
public:

	LoopHashSampler(
		LoopHashLibraryOP library,
		LocalInserterOP inserter
	);

	~LoopHashSampler() override;

	/// @brief Load default values from the command line.
	void set_defaults();

	/// @brief Create a set of structures for the given range of residues and
	/// other parameters stored in this class.
	/// @param[in] start_pose The pose to sample.
	/// @param[out] lib_structs The resulting structures.
	/// @details The algorithm, in pseudocode:
	/// - For each residue between \p start_res and \p stop_res:
	///   - For each loop length in the database
	///     - Calculate the rigid-body transformation between the ends of a loop
	///       starting at the current residue and extending the current length.
	///     - For each radius from 0 to \p max_radius:
	///       - Do a radial hashmap search for that radius.
	///       - For each hit:
	///         - Discard if the RMSD to the starting pose is either too low or
	///           too high
	///         - Discard if there are Ramachandran outliers.
	///         - Otherwise, keep!
	///       - For each remaining hit:
	///         - Insert the hit into a copy of the starting pose.
	///         - Create a silent file representation of that pose and some
	///           information describing to how it was generated, and add it to
	///           \p lib_structs.
	void build_structures(
		const core::pose::Pose& start_pose,
		std::vector< core::io::silent::SilentStructOP > &lib_structs
	);

	/// @brief Not implemented! Create a set of structures with closed gaps.
	void close_gaps(
		const core::pose::Pose& start_pose,
		std::vector< core::pose::Pose > &lib_structs,
		core::Size loop_size
	);

	/// @brief Set the first residue to sample.
	void set_start_res( core::Size  value ) {  start_res_  = value; }
	/// @brief Set the last residue to sample.
	void set_stop_res ( core::Size  value ) {  stop_res_   = value; }
	void set_min_bbrms( core::Real  value ) {  min_bbrms_  = value; }
	void set_max_bbrms( core::Real  value ) {  max_bbrms_  = value; }
	void set_min_rms  ( core::Real  value ) {  min_rms_    = value; }
	void set_max_rms  ( core::Real  value ) {  max_rms_    = value; }
	void set_max_radius  ( core::Size value ) {  max_radius_    = value; }
	void set_max_struct  ( core::Size  value ) {  max_struct_    = value; }
	void set_max_struct_per_radius  ( core::Size  value ) {  max_struct_per_radius_    = value; }
	void set_max_nstruct  ( core::Size  value ) {  max_nstruct_    = value; }
	void set_nonideal  ( bool value ) {  nonideal_  = value; }
	void set_filter_by_phipsi ( bool value) {  filter_by_phipsi_ = value; }
	// This is meant for model creation, not mpi refinement!

	core::Size get_start_res() { return  start_res_ ; }
	core::Size get_stop_res () { return  stop_res_  ; }
	core::Real get_min_bbrms() { return  min_bbrms_ ; }
	core::Real get_max_bbrms() { return  max_bbrms_ ; }
	core::Real get_min_rms  () { return  min_rms_   ; }
	core::Real get_max_rms  () { return  max_rms_   ; }
	core::Size get_max_nstruct() { return  max_nstruct_; }
	bool       get_filter_by_phipsi() { return  filter_by_phipsi_; }

	/// @brief Pre-filter structures with a scorefunction.
	/// @details This is done using a chainbroken pose (before constraint
	/// minimization!) and is useful for experimentally derived scorefunctions
	/// (eg density).
	void use_prefiltering( core::scoring::ScoreFunctionOP score_filt, core::Size nstruct ) {
		score_filt_ = score_filt;
		nprefilter_ = nstruct;
	}

private:

	/// @brief Pointer to the library used for insertion.
	LoopHashLibraryOP library_;

	/// @brief Pointer to the insertion functor which provides the peptide
	/// insertion facility.
	LocalInserterOP inserter_;

	/// @brief Parameters for insertion positions.
	core::Size start_res_;
	core::Size stop_res_ ;
	core::Real min_bbrms_;
	core::Real max_bbrms_;
	core::Real min_rms_  ;
	core::Real max_rms_  ;
	core::Size max_struct_;
	core::Size max_struct_per_radius_;
	core::Size max_radius_;
	core::Size max_nstruct_;
	bool nonideal_;
	bool filter_by_phipsi_;

	/// @brief Pre-filtering options.
	core::Size nprefilter_;
	core::scoring::ScoreFunctionOP score_filt_;
};


}
}

#endif
