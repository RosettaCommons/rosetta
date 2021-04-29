// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file trRosetta_protocols/movers/trRosettaProtocolMover.hh
/// @brief The full trRosetta structure prediction protocol from Yang et al, converted to C++
/// and implemented as a mover.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_trRosetta_protocols_movers_trRosettaProtocolMover_HH
#define INCLUDED_protocols_trRosetta_protocols_movers_trRosettaProtocolMover_HH

// Unit headers
#include <protocols/trRosetta_protocols/movers/trRosettaProtocolMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.fwd.hh>
#include <protocols/minimization_packing/MinMover.fwd.hh>

// Core headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/options/OptionCollection.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

namespace protocols {
namespace trRosetta_protocols {
namespace movers {

/// @brief The modes for randomizing backbone dihedrals.
/// @details If a developer adds an entry to this list, she or he should be sure to
/// update trRosettaProtocolMover::get_randomization_mode_string_from_enum().  Also,
/// the basic/options/options_rosetta.py "trRosetta" OptionGroup should have the new
/// options added to the "backbone_randomization_mode" option.
enum class trRosettaProtocolBackboneRandomizationMode {
	classic = 1, //Keep this first; this is the method used by Yang et al.
	ramachandran, //Randomization biased by the Ramachandran preferences of each amino acid.
	bins, //Randomization based on the frequency of seeing residue i in bin X and residue i+1 in bin Y.  Keep second-to-last.
	NUM_ENTRIES = bins //Keep last, and keep bins second-to-last.
};

/// @brief The modes for minimizing the backbone.
/// @details If a developer adds an entry to this list, she or he should be sure to
/// update trRosettaProtocolMover::get_minimization_mode_string_from_enum().  Also,
/// the basic/options/options_rosetta.py "trRosetta" OptionGroup should have the new
/// options added to the "backbone_minimization_mode" option.
enum class trRosettaProtocolBackboneMinimizationMode {
	classic0 = 1, //Minimize using short-range constraints, then using medium, then using long.
	classic1, //Minimize using short- and medium-range constraints, then long.
	classic2, //Minimize using short-, medium-, and long-range constraints. KEEP THIS SECOND-TO-LAST.
	NUM_ENTRIES = classic2 //KEEP THIS LAST.
};

/// @brief The full trRosetta structure prediction protocol from Yang et al, converted to
/// C++ and implemented as a mover.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class trRosettaProtocolMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor.
	/// @details Initializes from the global options system.
	trRosettaProtocolMover();

	/// @brief Options constructor.
	/// @details Initializes from a local options collection.
	trRosettaProtocolMover( utility::options::OptionCollection const & options );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~trRosettaProtocolMover() override;

	/// @brief Indicate commandline flags that are relevant to this mover.
	static void register_options( bool const only_constraint_generator_options = false );

	/// @brief Given a backbone randomization mode enum, get a string corresponding to the enum.
	/// @details This function must be updated if entries are added to the trRosettaProtocolBackboneRandomizationMode enum class!
	static
	std::string get_randomization_mode_string_from_enum( trRosettaProtocolBackboneRandomizationMode const mode_enum );

	/// @brief Given a backbone randomization mode name, get the enum corresponding to the name.
	/// @details Throws if an invalid name is provided.
	static
	trRosettaProtocolBackboneRandomizationMode get_randomization_mode_enum_from_string( std::string const & mode_name );

	/// @brief Given a backbone minimization mode enum, get a string corresponding to the enum.
	/// @details This function must be updated if entries are added to the trRosettaProtocolBackboneMinimizationMode enum class!
	static
	std::string get_minimization_mode_string_from_enum( trRosettaProtocolBackboneMinimizationMode const mode_enum );

	/// @brief Given a backbone minimization mode name, get the enum corresponding to the name.
	/// @details Throws if an invalid name is provided.
	static
	trRosettaProtocolBackboneMinimizationMode get_minimization_mode_enum_from_string( std::string const & mode_name );

public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	//trRosettaProtocolMover & operator=( trRosettaProtocolMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Setters

#ifdef USE_TENSORFLOW

	/// @brief Set the amino acid sequence of the protein being predicted.
	/// @details Must be one-letter codes; no NCAAs allowed.  Overrides any setting for
	/// fasta_file_ (i.e. sets fasta_file_ to "").
	void set_sequence( std::string const & seq_in );

	/// @brief Checks that the current sequence is valid, and throws if it is not.
	void validate_sequence() const;

	/// @brief Set the FASTA file containing the sequence of the protein to predict.
	/// @details Overrides any setting for sequence_ (i.e. sets sequence_ to "").
	void set_fasta_file( std::string const & filename );

	/// @brief Set the multiple sequence alignment filename.
	/// @details Resets the ConstraintGenerator, if already loaded.
	/// @note Read from disk of MSA file is deferred to apply time.
	void set_msa_file( std::string const & filename );

	/// @brief Set whether we are using the constraint generator to set
	/// distance constraints.
	void set_use_distance_constraints( bool const setting );

	/// @brief Set whether we are using the constraint generator to set
	/// omega dihedral constraints.
	/// @details This is NOT the backbone omega dihedral.  It is the dihedral angle
	/// between CA1-CB1-CB2-CA2.
	void set_use_omega_constraints( bool const setting );

	/// @brief Set whether we are using the constraint generator to set
	/// theta dihedral constraints.
	/// @details This is the dihedral between N1-CA1-CB1-CB2.
	void set_use_theta_constraints( bool const setting );

	/// @brief Set whether we are using the constraint generator to set
	/// phi angle constraints.
	/// @details This is NOT the backbone phi dihedral.  It is the angle between
	/// CA1-CB1-CB2.
	void set_use_phi_constraints( bool const setting );

	/// @brief Set the backbone randomization mode.
	void set_backbone_randomization_mode( trRosettaProtocolBackboneRandomizationMode const mode_in );

	/// @brief Set the backbone randomization mode, using the mode name string.
	void set_backbone_randomization_mode( std::string const & mode_in );

	/// @brief Set the backbone minimization mode.
	void set_backbone_minimization_mode( trRosettaProtocolBackboneMinimizationMode const mode_in );

	/// @brief Set the backbone minimization mode, using the mode name string.
	void set_backbone_minimization_mode( std::string const & mode_in );

	/// @brief Set the probability of sampling a cis peptide bond in 'ramachandran' mode,
	/// at non-pre-proline and pre-proline positions.
	void set_ramachandran_mode_cis_probabilities( core::Real const non_prepro_prob, core::Real const prepro_prob );

	/// @brief Set the weights file to use for minimization stage 0, 1, 2, or 3.
	/// @details Triggers read of weights from disk!
	void set_scorefunction_for_minimization_stage( std::string const & weights_filename, core::Size const stage );

	/// @brief Set the scoring function used for fullatom refinement.
	/// @details Triggers read of weights from disk!
	/// @note Turns on atom_pair_constraint, dihedral_constraint, and angle_constraint terms if they
	/// are not on in the scorefunction that is loaded.  Weights filename is deliberately passed by
	/// copy.
	void set_scorefunction_for_fullatom_refinement( std::string weights_filename );

	/// @brief Set the probability cutoffs for distance, omega, theta, and phi.
	void
	set_prob_cutoffs(
		core::Real const dist_cutoff,
		core::Real const omega_cutoff,
		core::Real const theta_cutoff,
		core::Real const phi_cutoff
	);

	/// @brief Set the constraint weights for distance, omega, theta, and phi.
	void
	set_constraint_weights(
		core::Real const dist_weight,
		core::Real const omega_weight,
		core::Real const theta_weight,
		core::Real const phi_weight
	);

	/// @brief Sets whether glycine residues should be mutated to alanine during the centroid phase.
	/// @details True by default to match the original PyRosetta protocol.
	void set_mutate_gly_to_ala( bool const setting );

	/// @brief Set whether we do fullatom refinement (with FastRelax) at the end.  Default true.
	void set_fullatom_refinement( bool const setting );

#endif //USE_TENSORFLOW

public: //Getters

#ifdef USE_TENSORFLOW

	/// @brief Get the amino acid sequence of the protein being predicted.
	inline std::string const & sequence() const { return sequence_; }

	/// @brief Get the FASTA filename (containing the sequence of the protein to predict).
	/// @details Overrides any setting for sequence_.
	inline std::string const & fasta_file() const { return fasta_file_; }

	/// @brief Get the multiple sequence alignment filename.
	inline std::string const & msa_file() const { return msa_file_; }

	/// @brief Get whether we are using the constraint generator to set
	/// distance constraints.
	inline bool use_distance_constraints() const { return use_distance_constraints_; }

	/// @brief Get whether we are using the constraint generator to set
	/// omega dihedral constraints.
	/// @details This is NOT the backbone omega dihedral.  It is the dihedral angle
	/// between CA1-CB1-CB2-CA2.
	inline bool use_omega_constraints() const { return use_omega_constraints_; }

	/// @brief Get whether we are using the constraint generator to set
	/// theta dihedral constraints.
	/// @details This is the dihedral between N1-CA1-CB1-CB2.
	inline bool use_theta_constraints() const { return use_theta_constraints_; }

	/// @brief Get whether we are using the constraint generator to set
	/// phi angle constraints.
	/// @details This is NOT the backbone phi dihedral.  It is the angle between
	/// CA1-CB1-CB2.
	inline bool use_phi_constraints() const { return use_phi_constraints_; }

	/// @brief Get the current backbone randomization mode.
	inline
	trRosettaProtocolBackboneRandomizationMode
	get_backbone_randomization_mode() const {
		return backbone_randomization_mode_;
	}

	/// @brief Get the current backbone randomization mode's name, as a string.
	std::string get_backbone_randomization_mode_name() const;

	/// @brief Get the current backbone minimization mode.
	inline
	trRosettaProtocolBackboneMinimizationMode
	get_backbone_minimization_mode() const {
		return backbone_minimization_mode_;
	}

	/// @brief Get the current backbone minimization mode's name, as a string.
	std::string get_backbone_minimization_mode_name() const;

	/// @brief Get the probability of sampling a cis peptide bond in 'ramachandran' mode.
	/// @details This is for resides that are NOT followed by proline.
	inline core::Real ramachandran_mode_cis_probability_non_prepro() const { return ramachandran_mode_cis_probability_non_prepro_; }

	/// @brief Get the probability of sampling a cis peptide bond in 'ramachandran' mode.
	/// @details This is for resides that ARE followed by proline.
	inline core::Real ramachandran_mode_cis_probability_prepro() const { return ramachandran_mode_cis_probability_prepro_; }

	/// @brief Get the probability cutoff for distance constraints.
	inline core::Real dist_prob_cutoff() const { return dist_prob_cutoff_; }

	/// @brief Get the probability cutoff for omega dihedral constraints.
	inline core::Real omega_prob_cutoff() const { return omega_prob_cutoff_; }

	/// @brief Get the probability cutoff for theta dihedral constraints.
	inline core::Real theta_prob_cutoff() const { return theta_prob_cutoff_; }

	/// @brief Get the probability cutoff for phi angle constraints.
	inline core::Real phi_prob_cutoff() const { return phi_prob_cutoff_; }

	/// @brief Get the distance constraint weight.
	inline core::Real distance_constraint_weight() const { return distance_constraint_weight_; }

	/// @brief Get the omega dihedral constraint weight.
	inline core::Real omega_constraint_weight() const { return omega_constraint_weight_; }

	/// @brief Get the theta dihedral constraint weight.
	inline core::Real theta_constraint_weight() const { return theta_constraint_weight_; }

	/// @brief Get the phi angle constraint weight.
	inline core::Real phi_constraint_weight() const { return phi_constraint_weight_; }

	/// @brief Should glycine residues be mutated to alanine during the centroid phase?
	/// @details True by default to match the original PyRosetta protocol.
	inline bool mutate_gly_to_ala() const { return mutate_gly_to_ala_; }

	/// @brief Do we do full-atom refinement at the end (with FastRelax)?
	inline bool fullatom_refinement() const { return fullatom_refinement_; }

#endif //USE_TENSORFLOW

public: //Function overrides needed for the citation manager:

	/// @brief This mover is unpublished.  It returns Vikram K. Mulligan as its author.
	/// It also adds Yang et al. (2020) as a citiation for trRosetta itself.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private: // methods

#ifdef USE_TENSORFLOW

	/// @brief Initialize this mover from an options collection.
	/// @details Intended to be called once and once only, from the
	/// constructor for this mover.
	/// @note Can trigger read from disk!
	void init_from_options( utility::options::OptionCollection const & options );

	/// @brief Create the trRosettaConstraintGenerator, apply it to the pose, and return a vector
	/// of constraint objects.
	/// @details If probability_offset is provided, the probability thresholds are all shifted by
	/// this amount.  Used during fullatom refinement.
	utility::vector1< core::scoring::constraints::ConstraintCOP >
	generate_trRosetta_constraints(
		core::pose::Pose const & pose,
		core::Real const probability_offset = 0.0
	) const;

	/// @brief Remove previously-added trRosetta constraints from the pose.
	void
	remove_trRosetta_constraints(
		utility::vector1< core::scoring::constraints::ConstraintCOP > const & constraints,
		core::pose::Pose & pose
	) const;

	/// @brief Add constraints from a list to the pose, if the constraints are between residues that are
	/// separated by at least min_seqsep but less than max_seqsep residues.
	/// @details This does not clear constraints from the pose.
	/// @note It is "less than" max_seqsep, not "less than or equal to".
	void
	add_constraints_to_pose(
		core::pose::Pose & pose,
		utility::vector1< core::scoring::constraints::ConstraintCOP > const & trRosetta_constraints,
		core::Size const min_seqsep,
		core::Size const max_seqsep,
		bool const skip_glycine_positions = false
	) const;

	/// @brief Given a constraint, determine if it is an AtomPairConstraint, an
	/// AngleConstraint, or a DihedralConstraint, pull out the pair of residues
	/// that are constrained, and return the pair.
	/// @details Throws if type is unrecognized or if more than two residues are
	/// constrained.  Values of res1 and res2 are overwritten by this operation.
	void
	get_residues_from_constraint(
		core::Size & res1,
		core::Size & res2,
		core::scoring::constraints::ConstraintCOP const & cst
	) const;

	/// @brief Given a constraint, determine if it is an AtomPairConstraint, pull
	/// out the pair of residues that are constrained, and return the pair.
	/// @details If successful, values of res1 and res2 are overwritten by this
	/// operation, and the function returns "true".  Otherwise, values are not
	/// altered, and the function returns "false".
	bool
	get_residues_from_atom_pair_constraint(
		core::Size & res1,
		core::Size & res2,
		core::scoring::constraints::ConstraintCOP const & cst
	) const;

	/// @brief Given a constraint, determine if it is an AngleConstraint, pull
	/// out the pair of residues that are constrained, and return the pair.
	/// @details If successful, values of res1 and res2 are overwritten by this
	/// operation, and the function returns "true".  Otherwise, values are not
	/// altered, and the function returns "false".
	/// @note Throws if the middle residue doesn't match either the first or
	/// last.
	bool
	get_residues_from_angle_constraint(
		core::Size & res1,
		core::Size & res2,
		core::scoring::constraints::ConstraintCOP const & cst
	) const;

	/// @brief Given a constraint, determine if it is a DihedralConstraint, pull
	/// out the pair of residues that are constrained, and return the pair.
	/// @details If successful, values of res1 and res2 are overwritten by this
	/// operation, and the function returns "true".  Otherwise, values are not
	/// altered, and the function returns "false".
	/// @note Throws if the middle residues doesn't match either the first or
	/// last.
	bool
	get_residues_from_dihedral_constraint(
		core::Size & res1,
		core::Size & res2,
		core::scoring::constraints::ConstraintCOP const & cst
	) const;

	/// @brief Perform minimization using short-range constraints, then medium-range, then long-range.
	/// @details Note that constraints are cumulative (short-range only in round 1, short- and medium-range
	/// in round 2, and short-, medium-, and long-range in round 3).
	void
	perform_classic0_minimization_protocol(
		core::pose::Pose & pose,
		utility::vector1< core::scoring::constraints::ConstraintCOP > const & trRosetta_constraints,
		protocols::moves::MoverOP repeat_minmover0,
		protocols::moves::MoverOP minmover1,
		protocols::moves::MoverOP minmover3
	) const;

	/// @brief Perform minimization using short-range and medium-range constraints, then long-range.
	/// @details Note that constraints are cumulative (short-range and medium-range
	/// in round 1, and short-, medium-, and long-range in round 2).
	void
	perform_classic1_minimization_protocol(
		core::pose::Pose & pose,
		utility::vector1< core::scoring::constraints::ConstraintCOP > const & trRosetta_constraints,
		protocols::moves::MoverOP repeat_minmover0,
		protocols::moves::MoverOP minmover1,
		protocols::moves::MoverOP minmover3
	) const;


	/// @brief Perform minimization using short-, medium-, and long-range constraints all in one go.
	void
	perform_classic2_minimization_protocol(
		core::pose::Pose & pose,
		utility::vector1< core::scoring::constraints::ConstraintCOP > const & trRosetta_constraints,
		protocols::moves::MoverOP repeat_minmover0,
		protocols::moves::MoverOP minmover1,
		protocols::moves::MoverOP minmover3
	) const;

	/// @brief Perform the minimization steps common to the above three protocols.
	void
	perform_classic_minsteps(
		core::pose::Pose & pose,
		protocols::moves::MoverOP repeat_minmover0,
		protocols::moves::MoverOP minmover1,
		protocols::moves::MoverOP minmover3
	) const;

	/// @brief Set up a MinMover for use by this protocol.
	protocols::minimization_packing::MinMoverOP
	configure_minmover(
		core::kinematics::MoveMapOP movemap,
		core::scoring::ScoreFunctionCOP sfxn,
		core::Size const maxiters,
		bool const do_cartesian
	) const;

	/// @brief Convert a pose to fullatom from centroid.
	void
	convert_to_fullatom(
		core::pose::Pose & pose
	) const;

	/// @brief Do all-atom FastRelax refinement.
	void
	do_fullatom_refinement(
		core::pose::Pose & pose
	) const;

	/// @brief Store the constraints score in the pose.
	void
	store_constraints_score(
		core::pose::Pose & pose,
		std::string const & metric_suffix
	) const;

	/// @brief Store the RMSD to native in the pose.
	/// @details Does nothing if native_pose == nullptr.
	void
	store_rmsd(
		core::pose::Pose & pose,
		core::pose::PoseCOP const & native_pose,
		std::string const & metric_suffix
	) const;

	/// @brief Make a centroid pose from the sequence stored in sequence_.
	/// @details Does not check that the sequence is valid first!
	core::pose::Pose make_centroid_pose_from_sequence() const;

	/// @brief Set or randomize the backbone phi, psi, and omega angles based on the backbone_randomization_mode_ setting.
	/// @details Calls one of randomize_backbone_dihedrals_classic(), randomize_backbone_dihedrals_rama_prepro(), or
	/// randomize_backbone_dihedrals_by_bins().
	void randomize_backbone_dihedrals( core::pose::Pose & pose ) const;

	/// @brief Set or randomize the backbone phi, psi, and omega angles using the method described in Yang et al.
	/// @details This sets each residue's dihedrals to one of the following phi/psi combinations, chosen randomly:
	/// -140  153 180 0.135 B
	///  -72  145 180 0.155 B
	/// -122  117 180 0.073 B
	///  -82  -14 180 0.122 A
	///  -61  -41 180 0.497 A
	///   57   39 180 0.018 L
	void randomize_backbone_dihedrals_classic( core::pose::Pose & pose ) const;

	/// @brief Set or randomize the backbone phi, psi, and omega angles based on the Ramachandran preferences of
	/// each amino acid type (using the rama_prepro scoreterm).
	void randomize_backbone_dihedrals_rama_prepro( core::pose::Pose & pose ) const;

	/// @brief Set or randomize the backbone phi, psi, and omega angles based on the probability of observing
	/// residue type i in backbone bin X and residue type i+1 in backbone bin Y.
	void randomize_backbone_dihedrals_by_bins( core::pose::Pose & pose ) const;

	/// @brief Mutate all the glycines to alanine.
	/// @returns A StoredResidueSubsetSelector that preserves the glycine positions, for later mutating
	/// them all back to alanines.
	core::select::residue_selector::ResidueSelectorCOP do_mutate_gly_to_ala( core::pose::Pose & pose ) const;

	/// @brief Mutate the alanines that were previously glycines back to glycine.
	void do_mutate_ala_to_gly( core::select::residue_selector::ResidueSelectorCOP gly_selector, core::pose::Pose & pose ) const;

	/// @brief Ensure that the termini have the N-terminal and C-terminal protein types.
	void add_termini( core::pose::Pose & pose ) const;

	/// @brief Set energy method options for a scorefunction (to turn on hb_cen_soft).
	void set_emethod_options( core::scoring::ScoreFunction & sfxn ) const;

	/// @brief Perform up to 5 rounds of minimization, or until the score is less than 10.0.
	void
	remove_clashes(
		core::pose::Pose & pose,
		protocols::moves::Mover & minmover,
		core::scoring::ScoreFunction const & sfxn
	) const;

#endif //USE_TENSORFLOW

private: // data

#ifdef USE_TENSORFLOW

	/// @brief The sequence of the protien being predicted.
	/// @details One-letter codes; only the 20 canonical amino acids allowed.
	mutable std::string sequence_;

	/// @brief A file containing the one-letter sequence of the protein to predict.
	std::string fasta_file_;

	/// @brief Filename for the multiple sequence alignment file.
	std::string msa_file_;

	/// @brief A constraint generator for adding trRosetta constraints to the pose.
	mutable protocols::trRosetta_protocols::constraint_generators::trRosettaConstraintGeneratorOP constraint_generator_;

	/// @brief Are we using the constraint generator to set distance constraints?
	bool use_distance_constraints_ = true;

	/// @brief Are we using the constraint generator to set omega dihedral constraints?
	/// @details This is NOT the backbone omega dihedral.  It is the dihedral angle
	/// between CA1-CB1-CB2-CA2.
	bool use_omega_constraints_ = true;

	/// @brief Are we using the constraint generator to set theta dihedral constraints?
	/// @details This is the dihedral between N1-CA1-CB1-CB2.
	bool use_theta_constraints_ = true;

	/// @brief Are we using the constraint generator to set phi angle constraints?
	/// @details This is NOT the backbone phi dihedral.  It is the angle between
	/// CA1-CB1-CB2.
	bool use_phi_constraints_ = true;

	/// @brief The mode for randomizing the backbone.  Defaults to 'classic'.
	trRosettaProtocolBackboneRandomizationMode backbone_randomization_mode_ = trRosettaProtocolBackboneRandomizationMode::classic;

	/// @brief The mode for minimizing the backbone.  Defaults to 'classic2'.
	trRosettaProtocolBackboneMinimizationMode backbone_minimization_mode_ = trRosettaProtocolBackboneMinimizationMode::classic2;

	/// @brief The probability of sampling a cis peptide bond at a non-preproline position
	/// when 'ramachandran' backbone randomization mode is used.
	core::Real ramachandran_mode_cis_probability_non_prepro_ = 0.0005;

	/// @brief The probability of sampling a cis peptide bond at a preproline position
	/// when 'ramachandran' backbone randomization mode is used.
	core::Real ramachandran_mode_cis_probability_prepro_ = 0.05;

	/// @brief The ScoreFunction used for initial minimization (stage 0).
	core::scoring::ScoreFunctionOP sfxn0_;

	/// @brief The ScoreFunction used for stage 1 minimization.
	core::scoring::ScoreFunctionOP sfxn1_;

	/// @brief The ScoreFunction used for stage 2 (Van der Waals) minimization.
	core::scoring::ScoreFunctionOP sfxn2_;

	/// @brief The ScoreFunction used for stage 3 (Cartesian) minimization.
	core::scoring::ScoreFunctionOP sfxn3_;

	/// @brief The ScoreFunction used for fullatom refinement steps.
	core::scoring::ScoreFunctionOP sfxn_fullatom_;

	/// @brief The cutoff for considering a constraint to be medium-range.
	/// @details If sequence separation (in residues) is LESS than this, it is a short-range constraint.  If
	/// GREATER THAN or EQUAL to this, it is medium-range (or long-range).
	core::Size medium_range_seqsep_cutoff_ = 12;

	/// @brief The cutoff for considering a constraint to be long-range.
	/// @details If sequence separation (in residues) is LESS than this, it is a medium-range constraint (or
	/// short-range).  If GREATER THAN or EQUAL to this, it is long-range.
	core::Size long_range_seqsep_cutoff_ = 24;

	/// @brief Probability cutoff for distance constraints.  Default 0.05.
	core::Real dist_prob_cutoff_ = 0.05;

	/// @brief Probability cutoff for omega dihedral constraints.  Default 0.55.
	core::Real omega_prob_cutoff_ = 0.55;

	/// @brief Probability cutoff for theta dihedral constraints.  Default 0.55.
	core::Real theta_prob_cutoff_ = 0.55;

	/// @brief Probability cutoff for phi angle constraints.  Default 0.65.
	core::Real phi_prob_cutoff_ = 0.65;

	/// @brief Should glycine residues be mutated to alanine during the centroid phase?
	/// @details True by default to match the original PyRosetta protocol.
	bool mutate_gly_to_ala_ = true;

	/// @brief Do we do full-atom refinement with FastRelax at the end?  Default true.
	bool fullatom_refinement_ = true;

	/// @brief Weighting factor on the distance constraints.  Default 1.0.
	core::Real distance_constraint_weight_ = 1.0;

	/// @brief Weighting factor on the phi angle constraints.  Default 1.0.
	core::Real phi_constraint_weight_ = 1.0;

	/// @brief Weighting factor on the omega dihedral angle constraints.  Default 1.0.
	core::Real omega_constraint_weight_ = 1.0;

	/// @brief Weighting factor on the theta dihedral angle constraints.  Default 1.0.
	core::Real theta_constraint_weight_ = 1.0;

#endif //USE_TENSORFLOW

};

std::ostream &
operator<<( std::ostream & os, trRosettaProtocolMover const & mover );

} //movers
} //trRosetta_protocols
} //protocols

#endif //protocols_trRosetta_protocols_movers_trRosettaProtocolMover_HH
