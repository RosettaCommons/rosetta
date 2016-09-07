// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/BackboneGridSampler.hh
/// @brief  Headers for BackboneGridSampler.cc.  Samples conformations of a repeating chain of a residue type
/// by grid-sampling mainchain torsion values and setting all residues in a range to have the same mainchain
/// torsion values.
/// @details Note that this mover throws away the input pose and generates new geometry.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_BackboneGridSampler_hh
#define INCLUDED_protocols_helical_bundle_BackboneGridSampler_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/helical_bundle/BackboneGridSampler.fwd.hh>
#include <protocols/helical_bundle/BackboneGridSamplerHelper.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.fwd.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

//JD2:
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/grid/CartGrid.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace helical_bundle {

class BackboneGridSampler : public protocols::moves::Mover
{
public: //Typedefs

	typedef protocols::cyclic_peptide::PeptideStubMover PeptideStubMover;
	typedef protocols::cyclic_peptide::PeptideStubMoverOP PeptideStubMoverOP;
	typedef protocols::cyclic_peptide::PeptideStubMoverCOP PeptideStubMoverCOP;

	typedef protocols::helical_bundle::BackboneGridSamplerHelper BackboneGridSamplerHelper;
	typedef protocols::helical_bundle::BackboneGridSamplerHelperOP BackboneGridSamplerHelperOP;
	typedef protocols::helical_bundle::BackboneGridSamplerHelperCOP BackboneGridSamplerHelperCOP;

public:
	BackboneGridSampler();
	BackboneGridSampler( BackboneGridSampler const &src );
	~BackboneGridSampler() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	///
	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;

	std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	) override;

public:
	////////////////////////////////////////////////////////////////////////////////
	//          PUBLIC FUNCTIONS                                                  //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Set the maximum number of samples for the mover.
	/// @details If the number of gridpoints based on user options exceeds this number, an error is thrown
	/// and the mover aborts.  This is to prevent unreasonably large calculations from being attempted.
	void set_max_samples( core::Size const val ) { max_samples_=val; return; }

	/// @brief Get the maximum number of samples for the mover.
	/// @details If the number of gridpoints based on user options exceeds this number, an error is thrown
	/// and the mover aborts.  This is to prevent unreasonably large calculations from being attempted.
	core::Size max_samples() const { return max_samples_; }

	/// @brief Sets whether the selection should be for the lowest score value (true) or highest (false).
	///
	void set_selection_low( bool const val ) { select_low_=val; return; }

	/// @brief Returns whether the selection should be for the lowest score value (true) or highest (false).
	///
	bool selection_low() { return select_low_; }

	/// @brief Sets the mover that will be applied to all helical bundles generated prior to energy evaluation.
	/// @details Note: if this is used, there is no guarantee that the resulting geometry will still lie within the
	/// parameter space.  (That is, this mover could move the backbone.)
	void set_prescoring_mover ( protocols::moves::MoverOP mover );

	/// @brief Sets the filter that will be applied to all helical bundles generated prior to energy evaluation.
	/// @details See the pre_scoring_filter_ private member variable for details.
	void set_prescoring_filter ( protocols::filters::FilterOP filter );

	/// @brief Returns "true" if and only if a preselection mover has been assigned.
	///
	bool prescoring_mover_exists() const { return pre_scoring_mover_exists_; }

	/// @brief Returns "true" if and only if a preselection filter has been assigned.
	///
	bool prescoring_filter_exists() const { return pre_scoring_filter_exists_; }

	/// @brief Set whether the mover dumps pdbs or not.
	///
	void set_pdb_output( bool const val ) { dump_pdbs_=val; return; }

	/// @brief Returns whether the mover dumps pdbs or not.
	///
	bool pdb_output() const { return dump_pdbs_; }

	/// @brief Sets the filename prefix for PDB output.
	/// @details PDB files are of the form <prefix>_#####.pdb.
	void set_pdb_prefix( std::string const &prefix ) { pdb_prefix_=prefix; return; }

	/// @brief Access the filename prefix for PDB output.
	/// @details PDB files are of the form <prefix>_#####.pdb.
	std::string pdb_prefix() { return pdb_prefix_; }

	/// @brief Sets the scorefunction for this mover.
	/// @details This must be done before calling the apply() function.
	void set_sfxn( core::scoring::ScoreFunctionOP sfxn_in ) {
		runtime_assert_string_msg( sfxn_in, "In BackboneGridSampler::set_sfxn(): A null scorefunction pointer was provided.  This should not happen.  Consult a developer or a mortician." );
		sfxn_=sfxn_in;
		sfxn_set_=true;
		return;
	}

	/// @brief Returns whether the scorefunction has been set.
	///
	bool sfxn_set() const { return sfxn_set_; }

	/// @brief Set the nstruct mode.
	/// @details If true, each job samples one set of mainchain torsion values.  If false, every job samples
	/// every set of mainchain torsion values.  False by default.
	void set_nstruct_mode( bool const &val ) { nstruct_mode_=val; return; }

	/// @brief Get the nstruct mode.
	/// @details If true, each job samples one set of mainchain torsion values.  If false, every job samples
	/// every set of mainchain torsion values.  False by default.
	bool nstruct_mode( ) const { return nstruct_mode_; }

	/// @brief Set the nstruct repeats.
	/// @details This is set to 1 by default, which means that each nstruct number correspnds to a different
	/// set of mainchain torsion values.  If set greater than 1, then multiple consecutive nstruct numbers will
	/// correspond to the same mainchain torsion values.  This allows combinatorially combining this mover's sampling
	/// with another, similar mover's sampling.
	void set_nstruct_repeats( core::Size const val ) { nstruct_mode_repeats_=val; return; }

	/// @brief Get the nstruct repeats.
	/// @details This is set to 1 by default, which means that each nstruct number correspnds to a different
	/// set of mainchain torsion values.  If set greater than 1, then multiple consecutive nstruct numbers will
	/// correspond to the same mainchain torsion values.  This allows combinatorially combining this mover's sampling
	/// with another, similar mover's sampling.
	core::Size nstruct_repeats( ) const {
		if ( nstruct_mode_repeats_ < 1 ) return 1;
		return nstruct_mode_repeats_;
	}

	/// @brief Add a mainchain torsion to sample, the range of values that will be sampled, and the number of samples.
	/// @details The residue_index is the index in the repeating unit (1st residue, 2nd residue, etc.).  The torsion_index
	/// is the mainchain torsion index in this residue.  Sampled values will go from start_of_range to end_of_range, with the
	/// total number of samples given by the samples parameter.  If the range is -180 to 180, the samples are adjusted so that
	/// 180 is not sampled twice.
	void add_torsion_to_sample(
		core::Size const residue_index,
		core::Size const torsion_index,
		core::Real const &start_of_range,
		core::Real const &end_of_range,
		core::Size const samples
	);

	/// @brief Add a mainchain torsion to fix, and that torsion's value.
	/// @details The residue_index is the index of the residue in the repeating unit (1st, 2nd, 3rd, etc.).
	/// The torsion_index is the mainchain torsion index in this residue, and the torsion_value is the
	/// value at which to fix this residue.
	void add_torsion_to_fix(
		core::Size const residue_index,
		core::Size const torsion_index,
		core::Real const &torsion_value
	);

	/// @brief Set the number of residues in the pose to build.
	/// @details This is actually the number of REPEATS.  If there are two residues per repeat, you'd have twice this number of resiudes.
	void set_nres( core::Size const nres_in ) {
		if ( nres_in < 1 ) {
			utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::set_nres(): The number of residues must be greater than zero." );
		}
		nres_ = nres_in;
		return;
	}

	/// @brief Get the number of residues in the pose to build.
	/// @details This is actually the number of REPEATS.  If there are two residues per repeat, you'd have twice this number of resiudes.
	core::Size nres() const { return nres_; }

	/// @brief Set the residue type.
	/// @details Sets the residue type for the nth residue in the list of residues defining the repeating unit.
	void set_resname( core::Size const res_index, std::string const &resname_in ) {
		if ( res_index < 1 || res_index > resname_.size() ) utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::set_resname(): The residue index is out of range.  It must be greater than zero and less than or equal to the number of residues in the repeating unit." );
		resname_[res_index]=resname_in; return;
	}

	/// @brief Get the residue type.
	/// @details Returns the residue type (full name string) for the nth residue in the list of residues defining the repeating unit.
	std::string resname( core::Size const res_index ) const {
		if ( res_index < 1 || res_index > resname_.size() ) utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::resname(): The residue index is out of range.  It must be greater than zero and less than or equal to the number of residues in the repeating unit." );
		return resname_[res_index];
	}

	/// @brief Set whether then ends should be capped (N-acetyl, C-methylamide).
	///
	void set_cap_ends(bool cap) { cap_ends_=cap; return; }

	/// @brief Return whether the ends should be capped (N-acetyl, C-methylamide).
	///
	bool cap_ends() const { return cap_ends_; }

	/// @brief Set up the PeptideStubMover object that will be used to build geometry.
	///
	void set_up_peptide_stub_mover();

	/// @brief Has the PeptideStubMover been initialized?
	///
	bool peptide_stub_mover_initialized() const { return peptide_stub_mover_initialized_; }

	/// @brief Reset the torsions_to_sample_, torsions_to_fix_, and resname_ vectors, and
	/// initialize them based on the value of residues_per_repeat_.
	void reset_and_initialize_data();

	/// @brief Set the number of residues per repeat.
	/// @details This resets the torsions_to_sample_, torsions_to_fix_, and resname_ vectors, and
	/// must therefore be called BEFORE setting up the torsions to sample, torsions to fix, or
	/// residue types.
	void set_residues_per_repeat( core::Size const val );

	/// @brief Get the number of residues per repeat.
	///
	core::Size residues_per_repeat() const { return residues_per_repeat_; }

private:
	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE DATA                                                      //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Should the parallel sampling be done based on the job (nstruct number)?
	/// @details If true, each job samples one set of mainchain torsion values.  If false, every job samples
	/// every set of mainchain torsion values.  False by default.
	bool nstruct_mode_;

	/// @brief If nstruct_mode_ is true, how many times should each set of mainchain torsion values be repeated?
	/// @details This is set to 1 by default, which means that each nstruct number correspnds to a different
	/// set of mainchain torsion values.  If set greater than 1, then multiple consecutive nstruct numbers will
	/// correspond to the same mainchain torsion values.  This allows combinatorially combining this mover's sampling
	/// with another, similar mover's sampling.
	core::Size nstruct_mode_repeats_;

	/// @brief The selection type.
	/// @default If false, the pose with the highest score value is selected.  If true,
	/// the pose with the lowest score value is selected.  True by default.
	bool select_low_;

	/// @brief The maximum number of gridpoints allowed.
	/// @details If the number of gridpoints based on user options exceeds this number, an error is thrown
	/// and the mover aborts.  This is to prevent unreasonably large calculations from being attempted.
	/// Default value is ten thousand (10,000).
	core::Size max_samples_;

	/// @brief Owning pointer for an (optional) pre-selection mover applied to all helical bundles before energy evaluation.
	///
	protocols::moves::MoverOP pre_scoring_mover_;

	/// @brief Bool determining whether there exists a pre-selection mover that wlil be applied.
	///
	bool pre_scoring_mover_exists_;

	/// @brief Owning pointer for an (optional) pre-selection filter applied to all helical bundles after the pre-selection mover but before
	/// picking the lowest-energy solution.  If PDBs are dumped, only those passing filters are dumped.
	protocols::filters::FilterOP pre_scoring_filter_;

	/// @brief Bool determining whether a pre-selection filter has been set.
	///
	bool pre_scoring_filter_exists_;

	/// @brief Dump a PDB file for each bundle generated?  False by default.
	///
	bool dump_pdbs_;

	/// @brief PDB filename prefix.  Filename will be of the form <prefix>_#####.pdb.
	/// @details  Defaults to "bgs_out".
	std::string pdb_prefix_;

	/// @brief Has the scorefunction been set?
	/// @details False by default.
	bool sfxn_set_;

	/// @brief The scorefunction that this mover will use to pick the lowest-energy bundle.
	/// @details Must be set prior to calling apply() function.
	core::scoring::ScoreFunctionOP sfxn_;

	/// @brief The number of residues in the repeating unit.
	/// @details  The repeating unit in the helix need not be a single residue, particularly if the helix is composed of residues of an alternating type.
	core::Size residues_per_repeat_;

	/// @brief The list of mainchain torsions to sample, and the range of values.
	/// @details The outer vector corresponds to the list of residues in the repeating unit.
	utility::vector1 < /*Index of residue in repeating unit*/ utility::vector1 < std::pair <core::Size /*mainchain torsion index*/, std::pair < std::pair < core::Real /*start of range to sample*/, core::Real /*end of range to sample*/ >, core::Size /*samples*/ > > > > torsions_to_sample_;

	/// @brief The list of mainchain torsions to fix, and the values at which to fix them.
	/// @details  The outer vecotr corresponds to the list of residues in the repeating unit.
	utility::vector1 < /*Index of residue in repeating unit*/ utility::vector1 < std::pair <core::Size /*The index to fix*/, core::Real /*The value to fix*/ > > > torsions_to_fix_;

	/// @brief The number of residues in the pose to build.
	/// @details This is actually the number of REPEATS.  If there are two residues per repeat, you'd have twice this number of resiudes.
	core::Size nres_;

	/// @brief The residue type (full name) for the pose to build.
	/// @details This is a vector of strings, with one entry for each residue in the repeating unit.
	utility::vector1 <std::string> resname_;

	/// @brief Should we add terminal caps (N-acetyl, C-methylamidation) to the polymer?
	/// @details If false, the ends are just the regular terminal types (NH3+, COO-).
	/// False by default.
	bool cap_ends_;

	/// @brief A PeptideStubMover used to build the new pose.
	///
	PeptideStubMoverOP peptide_stub_mover_;

	/// @brief Has the PeptideStubMover been initialized?
	///
	bool peptide_stub_mover_initialized_;

private:
	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Is a value in a list?
	///
	bool is_in_list( core::Size const val, utility::vector1 < core::Size> const &list ) const;

	/// @brief Calculate the number of grid points that will be sampled, based on the options set by the user.
	///
	core::Size calculate_total_samples() const;

}; //BackboneGridSampler class

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_BackboneGridSampler_hh
