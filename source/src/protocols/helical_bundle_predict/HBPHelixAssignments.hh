// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBPHelixAssignments.hh
/// @brief A class for storing the helix assignments for a pose.  This can represent those proposed from an input file, or those at the current state of a trajectory.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_protocols_helical_bundle_predict_HBPHelixAssignments_hh
#define INCLUDED_protocols_helical_bundle_predict_HBPHelixAssignments_hh

#include <protocols/helical_bundle_predict/HBPHelixAssignments.fwd.hh>

// Protocols headers
#include <protocols/helical_bundle/BundleParametrizationCalculator.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// STL headers
#include <map>

namespace protocols {
namespace helical_bundle_predict {

typedef core::Size HBPStartPosition;
typedef core::Size HBPEndPosition;

/////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Helical parameters stored by the HBPHelix class.
/// @details Allows the HBPHelix class to *avoid* storing these parameters if it needn't.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HBPHelixParameters : public utility::VirtualBase {

	friend class ::HBPHelixAssignmentsTests; //Needed to allow unit tests to access private data.

public:

	/// @brief Default constructor.
	HBPHelixParameters();

	/// @brief Copy constructor.
	HBPHelixParameters(HBPHelixParameters const & src);

	/// @brief Destructor.
	~HBPHelixParameters() override;

	/// @brief Create a copy of this object and return an owning pointer to the copy.
	HBPHelixParametersOP clone() const;

public: //Setters

	/// @brief Set minimum r0 value.
	void set_r0_min( core::Real const &setting );

	/// @brief Set maximum r0 value.
	void set_r0_max( core::Real const &setting );

	/// @brief Set minimum omega0 value.
	/// @details Settings should be in radians.
	void set_omega0_min( core::Real const &setting );

	/// @brief Set maximum omega0 value.
	/// @details Settings should be in radians.
	void set_omega0_max( core::Real const &setting );

	/// @brief Set minimum delta_omega1 value.
	/// @details Settings should be in radians.
	void set_delta_omega1_min( core::Real const &setting );

	/// @brief Set maximum delta_omega1 value.
	/// @details Settings should be in radians.
	void set_delta_omega1_max( core::Real const &setting );

	/// @brief Set the bundle parameterization calculator.
	/// @details Input is cloned.
	void set_calculator( protocols::helical_bundle::BundleParametrizationCalculatorCOP const & calculator_in );

	/// @brief Create a new bundle parameterization calculator, and initialize it from a Crick params file.
	/// @details WARNING!  Triggers read from disk!
	void create_calculator_from_file( std::string const & filename );

	/// @brief Access the bundle calculator (const-access).
	inline protocols::helical_bundle::BundleParametrizationCalculatorCOP bundle_calculator() const { return bundle_calculator_; }

	/// @brief Access the bundle calculator (nonconst-access).
	inline protocols::helical_bundle::BundleParametrizationCalculatorOP bundle_calculator() { return bundle_calculator_; }

public: //Getters

	/// @brief Get the r0 range.
	inline std::pair< core::Real, core::Real > const & r0_range() const { return r0_range_; }

	/// @brief Get the omega0 range.
	/// @details Settings should be in radians.
	inline std::pair< core::Real, core::Real > const & omega0_range() const { return omega0_range_; }

	/// @brief Get the delta_omega1 range.
	/// @details Settings should be in radians.
	inline std::pair< core::Real, core::Real > const & delta_omega1_range() const { return delta_omega1_range_; }

private: //Data

	/// @brief The range of r0 values that are allowed.
	std::pair< core::Real, core::Real > r0_range_;

	/// @brief The range of omega0 values that are allowed.
	/// @details Settings should be in radians.
	std::pair< core::Real, core::Real > omega0_range_;

	/// @brief The range of delta_omega1 values that are allowed.
	/// @details Settings should be in radians.
	std::pair< core::Real, core::Real > delta_omega1_range_;

	/// @brief Helper object to calculate Crick parameterization.
	protocols::helical_bundle::BundleParametrizationCalculatorOP bundle_calculator_;

};

/////////////////////////////////////////////////////////////////////////////////////////////

/// @brief An individual helix.  This class stores start and end positions,
/// plus helical parameters.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HBPHelix : public utility::VirtualBase {

	friend class ::HBPHelixAssignmentsTests; //Needed to allow unit tests to access private data.

public:

	/// @brief Default constructor.
	HBPHelix();

	/// @brief Copy constructor.
	/// @details Deep-copies paramters.
	HBPHelix(HBPHelix const & src);

	/// @brief Destructor.
	~HBPHelix() override;

	/// @brief Create a copy of this object and return an owning pointer to the copy.
	HBPHelixOP clone() const;

public: //Functions

	/// @brief Get the parameters (nonconst access).
	HBPHelixParametersOP parameters();

	/// @brief Get the parameters (read-only).
	/// @details Returns nullptr if no object exists.
	inline HBPHelixParametersCOP parameters() const { return parameters_; }

	/// @brief Is a sequence position within a helix?
	/// @details Returns false if either of startpos or endpos is zero.
	bool is_in_helix( core::Size const seqpos ) const;

	/// @brief Get the number of residues per turn of this helix type.
	core::Size residues_per_turn() const;

	/// @brief Generate a map of (TorsionID->torsion value) for all mainchain torsions in the current helix.
	/// @returns True for FAILURE (Crick parameters don't generate a helix), false for success.
	bool get_torsions_for_helix( std::map< core::id::TorsionID, core::Real > & torsion_map_out, core::pose::Pose const &pose_for_reference ) const;

	/// @brief Given the start and end position of a helix, determine from the residue composition what the repeating unit offset should be.
	/// @details For example, if the stretch had pattern alpha-alpha-beta-alpha, and the Crick params file expected alpha-alpha-alpha-beta, the offset
	/// would be 3.
	core::Size determine_repeating_unit_offset( core::pose::Pose const &pose, core::Size const start_position, core::Size const end_position ) const;

public: //Setters

	/// @brief Set the start position.
	void set_start_position( core::Size const position_in );

	/// @brief Set the end position.
	void set_end_position( core::Size const position_in );

	/// @brief Set the probability of nucleating a helix.
	void set_nucleation_prob( core::Real const & setting );

	/// @brief Set the probability of extending an existing helix.
	void set_extension_prob( core::Real const & setting );

	/// @brief Set the probability of shrinking an existing helix.
	void set_retraction_prob( core::Real const & setting );

	/// @brief Set the filename for the Crick params for this helix.
	/// @details Used to determine whether two helices are of a matching type.
	void set_crick_params_filename( std::string const & name_in );

public: //Getters

	/// @brief Get the start position.
	inline core::Size start_position() const { return start_position_; }

	/// @brief Get the end position.
	inline core::Size end_position() const { return end_position_; }

	/// @brief Get the probability of nucleating a helix.
	inline core::Real const & nucleation_prob() const { return nucleation_prob_; }

	/// @brief Get the probability of extending an existing helix.
	inline core::Real const & extension_prob() const { return extension_prob_; }

	/// @brief Get the probability of shrinking an existing helix.
	inline core::Real const & retraction_prob() const { return retraction_prob_; }

	/// @brief Get the filename for the Crick params for this helix.
	/// @details Used to determine whether two helices are of a matching type.
	inline std::string const & crick_params_filename() const { return crick_params_filename_; }

private: //Data

	/// @brief The start of the helix.
	HBPStartPosition start_position_;

	/// @brief The end of the helix.
	HBPEndPosition end_position_;

	/// @brief The probability of nucleating.
	core::Real nucleation_prob_;

	/// @brief The probability of extending a nucleated helix.
	core::Real extension_prob_;

	/// @brief The probability of retracting a nucleated helix.
	core::Real retraction_prob_;

	/// @brief The helical parameters.  Will be nullptr if parameters not stored.
	HBPHelixParametersOP parameters_;

	/// @brief The name of the Crick params file.  Used to determine whether two helices are of
	/// a matching type.
	std::string crick_params_filename_;

};

/////////////////////////////////////////////////////////////////////////////////////////////

/// @brief A class for storing the helix assignments for a pose.  This can represent those
/// proposed from an input file, or those at the current state of a trajectory.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HBPHelixAssignments : public utility::VirtualBase {

	friend class ::HBPHelixAssignmentsTests; //Needed to allow unit tests to access private data.

public:

	/// @brief Default constructor.
	HBPHelixAssignments();

	/// @brief Copy constructor.
	/// @details Note explicit deep-copy.
	HBPHelixAssignments(HBPHelixAssignments const & src);

	/// @brief Destructor.
	~HBPHelixAssignments() override;

	/// @brief Create a copy of this object and return an owning pointer to the copy.
	HBPHelixAssignmentsOP clone() const;

	/// @brief Assignment operator.
	/// @details Deep-clones helices.
	HBPHelixAssignments & operator=(HBPHelixAssignments const &src);

public: //Functions

	/// @brief Reset this object.
	void clear();

	/// @brief Set this object up from the contents of a file.
	/// @details Warning!  Reads Crick params files from disk!
	void initialize_from_file_contents( std::string const & file_contents );

	/// @brief Given a sequence index, determine whether it is in a helix.
	bool is_in_helix( core::Size const sequence_index ) const;

	/// @brief Roll a die and decide whether to nucleate a helix at the present position.
	/// @details If we do it, then we add one turn of helix as a new helix (updating the helices_ list),
	/// and we copy helix parameters from reference_helix_assignments.  The probabilities for nucleation are taken
	/// from reference_helix_assignments.
	/// @note Automatically merges helices if the set of residues in question overlaps with another helix.  In this case,
	/// the parameters of the other helix are chosen as the parameters to keep.
	/// @returns True if we nucleated here, false otherwise.
	bool try_nucleating_helix( core::Size const seqpos, core::pose::Pose const &pose, HBPHelixAssignments const & reference_helix_assignments );

	/// @brief Roll a die and decide whether to elongate a helix at the present position.
	/// @details If we do it, then we add one residue to the helix (updating the helices_ list).
	/// The probabilities for elongation are taken from the adjacent helix.
	/// @note Does nothing if the present position is not adjacent to a helix.  Automatically merges
	/// helices if the present position bridges two helices, and we decide to elongate.  In this case,
	/// the parameters of one helix are randomly chosen as the parameters to keep.
	/// @returns True if we elongated here, false otherwise.
	bool try_elongating_helix( core::Size const seqpos, core::pose::Pose const &pose, HBPHelixAssignments const & reference_helix_assignments );

	/// @brief Roll a die and decide whether to contract a helix at the present position.
	/// @details If we do it, then we subtract one residue from the helix (updating the helices_ list).
	/// The probabilities for retraction are taken from the adjacent helix.  If the helix shrinks to zero size,
	/// then we remove it from the helices_ list.
	/// @note Does nothing if the present position is not at one of the ends of a helix.
	/// @returns True if we retracted here, false otherwise.
	bool try_retracting_helix( core::Size const seqpos, core::pose::Pose const &pose, HBPHelixAssignments const & reference_helix_assignments );

	/// @brief Try perturbing the parameters for all helices.
	/// @details Updates individual helices if parameters are perturbed.
	bool try_perturbing_global_helix_parameters();

	/// @brief Try perturbing the parameters of individual helices, that are not set to copy globals.
	bool try_perturbing_local_helix_parameters();

	/// @brief Are we supposed to keep sampled helices within the user-defined ranges?
	inline bool confine_to_user_defined_helices() const { return confine_to_user_defined_helices_; }

	/// @brief Generate a map of torsion ID -> torsion value for each helix, and return a vector of these maps.
	/// @returns True for FAILURE, false for success.
	bool generate_torsion_values_for_helices( utility::vector1< std::map < core::id::TorsionID, core::Real > > & torsion_values_out, core::pose::Pose const & pose_for_reference ) const;

	/// @brief Given a sequence position, get the helix index that contains it.
	/// @details Returns zero if no helix contains this sequence position.
	core::Size get_containing_helix_index( core::Size const seqpos ) const;

private: //Functions

	/// @brief Get a nonconst owning pointer to the helix containing a given sequence position.
	/// @details Returns nullptr if no helix contains that sequence position.
	HBPHelixOP helix_from_seqpos( core::Size const position );

	/// @brief Get a const owning pointer to the helix containing a given sequence position.
	/// @details Returns nullptr if no helix contains that sequence position.
	HBPHelixCOP helix_from_seqpos( core::Size const position ) const;

	/// @brief Loop through all helices and find a helix that this sequence position is adjacent to.
	/// @details.  If the position is adjacent to no helix, returns nullptr.  If the position is adjacent
	/// to two helices (i.e. it bridges two helices separated by one residue), we pick one with 50/50 odds.
	HBPHelixOP find_adjacent_helix( core::Size seqpos );

	/// @brief Are two helices of the same type?
	bool helix_types_match( HBPHelix const &first_helix, HBPHelix const &second_helix ) const;

	/// @brief Check whether the addition or elongation of current_helix causes a new contiguous stretch of helix, and
	/// if so, replace helix records in the helices_ list with a single contiguous helix.
	/// @details When helices are merged, the helix contributing parameters is chosen randomly. If current_helix is not nullptr,
	/// then we do not allow that helix to contribute the parameters -- i.e. they are chosen randomly from one of the other helices,
	/// or we just use the other helix if there's only one.
	void check_for_helix_merges( HBPHelixCOP const current_helix = nullptr );

	/// @brief Merge two helices together.
	/// @details If current_helix is not nullptr, then we randomly pick one of the two helices to provide
	/// the parameters.  If current_helix is nullptr and it matches one of the two helices, then the other
	/// provides the parameters.
	void merge_helices( core::Size const first_helix_index, core::Size const second_helix_index, HBPHelixCOP const current_helix = nullptr );

	/// @brief Given two overlapping helices of different types, prevent them from overlapping.
	/// @details The cutoff point is randomly chosen within the overlap region.
	void resolve_mismatching_helix_overlap( core::Size const helix1_index, core::Size const helix2_index );

	/// @brief Given file contents, pull out the first BEGIN_GLOBALS...END_GLOBALS block, and throw an
	/// error message if any other occurrences of this exist in the file.
	/// @details Omits the BEGIN_GLOBALS and END_GLOBALS lines.
	void pull_out_global_block( std::string const &file_contents, std::string & global_block ) const;

	/// @brief Read the BEGIN_GLOBALS...END_GLOBALS block in the file contents and set up global options.
	/// @details Warning!  Reads Crick params files from disk!
	void initialize_global_options_from_file_contents( std::string const &file_contents );

	/// @brief Given file contents, pull out all of the stuff in between the BEGIN_HELIX and END_HELIX lines.
	/// @details Omits the BEGIN_HELIX and END_HELIX lines.
	void pull_out_helix_blocks( std::string const &file_contents, utility::vector1<std::string> &helix_blocks ) const;

	/// @brief Read BEGIN_HELIX...END_HELIX blocks in the file contents and set up local helix options.
	/// @details Warning!  Triggers read from disk if Crick params file is set.
	void initialize_helices_from_file_contents( std::string const &file_contents );

	/// @brief Add a helix and return an owning pointer to it.
	/// @details Parameters are initialized from globals.
	HBPHelixOP add_helix();

	/// @brief Add a helix and return an owning pointer to it, initializing its start, end, and parameters.
	/// @details Parameters are initialized from reference_helix.
	HBPHelixOP add_helix( core::Size const startpos, core::Size const endpos, HBPHelix const &reference_helix );

	/// @brief Remove a helix.
	void remove_helix( HBPHelixOP helix_to_remove );

private: //Data

	/// @brief Stored helices.
	utility::vector1< HBPHelixOP > helices_;

	/// @brief A calculator object for the Crick equations.
	protocols::helical_bundle::BundleParametrizationCalculatorOP global_bundle_calculator_;

	/// @brief Is r0 shared among all helices?
	/// @details False by default.
	bool common_r0_;

	/// @brief Is omega0 shared among all helices?
	/// @details False by default.
	bool common_omega0_;

	/// @brief Is delta_omega1 shared among all helices?
	/// @details False by default.
	bool common_delta_omega1_;

	/// @brief Prevent expansion of helices beyond user-defined helices?
	/// @details True by default.
	bool confine_to_user_defined_helices_;

	/// @brief Global min setting for r0.  Can be overridden in individual helices.
	core::Real global_r0_min_;

	/// @brief Global max setting for r0.  Can be overridden in individual helices.
	core::Real global_r0_max_;

	/// @brief Global min setting for omega0.  Can be overridden in individual helices.
	core::Real global_omega0_min_;

	/// @brief Global max setting for omega0.  Can be overridden in individual helices.
	core::Real global_omega0_max_;

	/// @brief Global min setting for delta_omega1.  Can be overridden in individual helices.
	core::Real global_delta_omega1_min_;

	/// @brief Global max setting for delta_omega1.  Can be overridden in individual helices.
	core::Real global_delta_omega1_max_;

	/// @brief The global Crick parameters file name.
	std::string global_crick_params_file_;

	/// @brief Probability of nucleating a helix.
	/// @details Defaults to 1%.
	core::Real nucleation_prob_;

	/// @brief Probability of extending a helix.
	/// @brief Defaults to 5%
	core::Real extension_prob_;

	/// @brief Probability of retracting a helix.
	/// @brief Defaults to 3%.
	core::Real retraction_prob_;

	/// @brief The probability of doing a global Crick parameter perturbation move.
	core::Real global_perturbation_prob_;

	/// @brief The magnitude of parameter perturbations.
	/// @details This is expressed as a fraction of the difference between the minimum and maximum parameter values.  The actual
	/// magnitude of perturbation is a random Gaussian times this value.
	core::Real global_fractional_perturbation_magnitude_;

};

/////////////////////////////////////////////////////////////////////////////////////////////

} //protocols
} //helical_bundle_predict



#endif //INCLUDED_protocols_helical_bundle_predict_HBPHelixAssignments_hh





