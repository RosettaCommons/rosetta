// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.hh
/// @brief A module that runs a trRosetta neural network on an input multiple sequence alignment and
/// uses the output to apply distance and/or angle constraints to a pose for subsequent structure
/// prediction or refinement.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_trRosetta_protocols_constraint_generators_trRosettaConstraintGenerator_hh
#define INCLUDED_protocols_trRosetta_protocols_constraint_generators_trRosettaConstraintGenerator_hh

// Unit headers
#include <protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Protocols headers
#ifdef USE_TENSORFLOW
#include <protocols/trRosetta/trRosettaOutputsBase.hh>
#include <protocols/trRosetta/trRosettaOutputs_v1.fwd.hh>
#include <core/id/AtomID.hh>
#endif //USE_TENSORFLOW

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

// Basic headers

// Utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace trRosetta_protocols {
namespace constraint_generators {

///@brief A module that runs a trRosetta neural network on an input multiple sequence alignment and uses the output to apply distance and/or angle constraints to a pose for subsequent structure prediction or refinement.
class trRosettaConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:

	/// @brief Default constructor.
	trRosettaConstraintGenerator();

#ifdef USE_TENSORFLOW
	/// @brief Options constructor.
	trRosettaConstraintGenerator( std::string const & msa_file );
#endif //USE_TENSORFLOW

	/// @brief OptionsCollection constructor.
	trRosettaConstraintGenerator( utility::options::OptionCollection const & opts );

	~trRosettaConstraintGenerator() override;

	static std::string
	class_name();

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const override;

	core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Provide citations to the passed CitationCollectionList.
	/// This allows the constraint generator to provide citations for itself
	/// and for any modules that it invokes.
	/// @details Calls the static version of this function, to allow citation info
	/// to be obtained without an instance.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

	/// @brief Provide citations for hte trRosettaConstraintGenerator, without needing an instance.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	static void provide_citation_info_static(basic::citation_manager::CitationCollectionList & citations);

public: //Setters

#ifdef USE_TENSORFLOW
	/// @brief Set the name of the multiple sequence alignment file.
	/// @details Calls reset().
	void set_msa_file( std::string const & filename_in );

	/// @brief Reset (delete) the trRosetta neural network output.
	/// @details Forces a re-run of the trRosetta neural network when
	/// the apply() function is called.
	void reset();

	/// @brief Set whether this will generate distance constraints.
	/// @details Distance constraints constrain the distance between pairs of amino acids.
	/// These are symmetric, and are only generated once per amino acid pair (since
	/// dist(a,b) == dist(b,a)).
	void set_generate_dist_constraints(bool const setting);

	/// @brief Set whether this will generate omega dihedral constraints.
	/// @details Omega constraints constrain the dihedral between CA1-CB1-CB2-CA2 in pairs
	/// of amino acids.  These are symmetric, and are only generated once per amino acid
	/// pair (since omega(a,b) == omega(b,a)).
	/// @note This is NOT omega the backbone dihedral torsion!
	void set_generate_omega_constraints(bool const setting);

	/// @brief Set whether this will generate theta dihedral constraints.
	/// @details Theta constraints constrain the dihedral between N1-CA2-CB1-CB2 in pairs
	/// of amino acids.  These are asymmetric (i.e. theta(a,b)!=theta(b,a)), so there are
	/// two per amino acid pair (unless a == b, which is skipped).
	void set_generate_theta_constraints(bool const setting);

	/// @brief Set whether this will generate phi angle constraints.
	/// @details Phi constraints constrain the angle between CA1-CB1-CB2 in pairs
	/// of amino acids.  These are asymmetric (i.e. phi(a,b)!=phi(b,a)), so there are
	/// two per amino acid pair (unless a == b, which is skipped).
	/// @note This is NOT phi the backbone dihedral torsion!
	void set_generate_phi_constraints(bool const setting);

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
#endif //USE_TENSORFLOW

public: //Getters

#ifdef USE_TENSORFLOW

	/// @brief Get the name of the multiple sequence alignment file.
	std::string const & msa_file() const { return msa_file_; }

	/// @brief Get whether this will generate distance constraints.
	/// @details Distance constraints constrain the distance between pairs of amino acids.
	/// These are symmetric, and are only generated once per amino acid pair (since
	/// dist(a,b) == dist(b,a)).
	inline bool generate_dist_constraints() const { return generate_dist_constraints_; }

	/// @brief Get whether this will generate omega dihedral constraints.
	/// @details Omega constraints constrain the dihedral between CA1-CB1-CB2-CA2 in pairs
	/// of amino acids.  These are symmetric, and are only generated once per amino acid
	/// pair (since omega(a,b) == omega(b,a)).
	/// @note This is NOT omega the backbone dihedral torsion!
	inline bool generate_omega_constraints() const { return generate_omega_constraints_; }

	/// @brief Get whether this will generate theta dihedral constraints.
	/// @details Theta constraints constrain the dihedral between N1-CA2-CB1-CB2 in pairs
	/// of amino acids.  These are asymmetric (i.e. theta(a,b)!=theta(b,a)), so there are
	/// two per amino acid pair (unless a == b, which is skipped).
	inline bool generate_theta_constraints() const { return generate_theta_constraints_; }

	/// @brief Get whether this will generate phi angle constraints.
	/// @details Phi constraints constrain the angle between CA1-CB1-CB2 in pairs
	/// of amino acids.  These are asymmetric (i.e. phi(a,b)!=phi(b,a)), so there are
	/// two per amino acid pair (unless a == b, which is skipped).
	/// @note This is NOT phi the backbone dihedral torsion!
	inline bool generate_phi_constraints() const { return generate_phi_constraints_; }

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

#endif //USE_TENSORFLOW

protected:

	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

private: //Functions

#ifdef USE_TENSORFLOW
	/// @brief Initialize this constraint generator from an options collection.
	/// @details Called once by constructor.
	void init_from_commandline( utility::options::OptionCollection const & opts );

	/// @brief Load the multiple sequence alignment file and run the neural network.
	/// @details Assumes that trRosetta_outputs_ is nullptr; throws if msa_file_ is empty.
	/// @note TRIGGERS READ FROM DISK.
	void load_msa_and_run_neural_network() const;

	/// @brief Given a pose, generate the constraints.
	/// @details Assumes that trRosetta_outputs_ is not a nullptr.  Throws if it is, in debug mode.
	core::scoring::constraints::ConstraintCOPs
	generate_constraints(
		core::pose::Pose const & pose
	) const;

	/// @brief Given a pose, generate the constraints using the version 1 model.
	/// @details Assumes that trRosetta_outputs_ is not a nullptr.  Throws if it is, in debug mode.
	core::scoring::constraints::ConstraintCOPs
	generate_constraints_v1(
		core::pose::Pose const & pose
	) const;

	/// @brief Generate the distance constraints using the version 1 model.
	/// @details New constraints are appended to outputvec.  Outputvec is not cleared
	/// by this operation.
	void
	generate_dist_constraints_v1(
		protocols::trRosetta::trRosettaOutputs_v1COP trRosetta_outputs,
		utility::vector1 <core::scoring::constraints::ConstraintCOP> & outputvec,
		utility::vector1< core::Real > const & dist_bins_background,
		utility::vector1< core::Real > const & dist_bins_vect,
		core::id::AtomID const & cb_atom_i,
		core::id::AtomID const & cb_atom_j,
		core::Size const n_dist_bins,
		core::Size const n_dist_bins_model,
		core::Size const ires,
		core::Size const jres,
		core::Real const prob_cutoff
	) const;

	/// @brief Generate the omega or theta dihedral constraints using the version 1 model.
	/// @details New constraints are appended to outputvec.  Outputvec is not cleared
	/// by this operation.
	void
	generate_omega_or_theta_constraints_v1(
		bool const generate_theta_constraints,
		protocols::trRosetta::trRosettaOutputs_v1COP trRosetta_outputs,
		utility::vector1 <core::scoring::constraints::ConstraintCOP> & outputvec,
		utility::vector1< core::Real > const & bins_vect,
		core::id::AtomID const & at1,
		core::id::AtomID const & at2,
		core::id::AtomID const & at3,
		core::id::AtomID const & at4,
		core::Size const n_bins,
		core::Size const n_bins_model,
		core::Size const ires,
		core::Size const jres,
		core::Real const prob_cutoff
	) const;

	/// @brief Generate the phi angle constraints using the version 1 model.
	/// @details New constraints are appended to outputvec.  Outputvec is not cleared
	/// by this operation.
	void
	generate_phi_constraints_v1(
		protocols::trRosetta::trRosettaOutputs_v1COP trRosetta_outputs,
		utility::vector1 <core::scoring::constraints::ConstraintCOP> & outputvec,
		utility::vector1< core::Real > const & phi_bins_vect,
		core::id::AtomID const & at1,
		core::id::AtomID const & at2,
		core::id::AtomID const & at3,
		core::Size const n_phi_bins,
		core::Size const n_phi_bins_model,
		core::Size const ires,
		core::Size const jres,
		core::Real const prob_cutoff
	) const;

#endif //USE_TENSORFLOW

private: //Variables

#ifdef USE_TENSORFLOW

	/// @brief The filename of a multiple sequence alignment file.
	std::string msa_file_;

	/// @brief Will this generate distance constraints?
	/// @details Distance constraints constrain the distance between pairs of amino acids.
	/// These are symmetric, and are only generated once per amino acid pair (since
	/// dist(a,b) == dist(b,a)).
	bool generate_dist_constraints_ = true;

	/// @brief Will this generate omega dihedral constraints?
	/// @details Omega constraints constrain the dihedral between CA1-CB1-CB2-CA2 in pairs
	/// of amino acids.  These are symmetric, and are only generated once per amino acid
	/// pair (since omega(a,b) == omega(b,a)).
	/// @note This is NOT omega the backbone dihedral torsion!
	bool generate_omega_constraints_ = true;

	/// @brief Will this generate theta dihedral constraints?
	/// @details Theta constraints constrain the dihedral between N1-CA2-CB1-CB2 in pairs
	/// of amino acids.  These are asymmetric (i.e. theta(a,b)!=theta(b,a)), so there are
	/// two per amino acid pair (unless a == b, which is skipped).
	bool generate_theta_constraints_ = true;

	/// @brief Will this generate phi angle constraints?
	/// @details Phi constraints constrain the angle between CA1-CB1-CB2 in pairs
	/// of amino acids.  These are asymmetric (i.e. phi(a,b)!=phi(b,a)), so there are
	/// two per amino acid pair (unless a == b, which is skipped).
	/// @note This is NOT phi the backbone dihedral torsion!
	bool generate_phi_constraints_ = true;

	/// @brief The version of the trRosetta model to use.  TODO make this an option when there's
	/// more than one version available.
	core::Size trRosetta_model_version_ = 1;

	/// @brief The outputs from the trRosetta neural network.
	mutable protocols::trRosetta::trRosettaOutputsBaseCOP trRosetta_outputs_;

	/// @brief Probability cutoff for distance constraints.  Default 0.05.
	core::Real dist_prob_cutoff_ = 0.05;

	/// @brief Probability cutoff for omega dihedral constraints.  Default 0.55.
	core::Real omega_prob_cutoff_ = 0.55;

	/// @brief Probability cutoff for theta dihedral constraints.  Default 0.55.
	core::Real theta_prob_cutoff_ = 0.55;

	/// @brief Probability cutoff for phi angle constraints.  Default 0.65.
	core::Real phi_prob_cutoff_ = 0.65;

	/// @brief Constraint weight for distance constraints.  Default 1.0.
	core::Real distance_constraint_weight_ = 1.0;

	/// @brief Constraint weight for omega dihedral constraints.  Default 1.0.
	core::Real omega_constraint_weight_ = 1.0;

	/// @brief Constraint weight for theta dihedral constraints.  Default 1.0.
	core::Real theta_constraint_weight_ = 1.0;

	/// @brief Constraint weight for phi anlge constraints.  Default 1.0.
	core::Real phi_constraint_weight_ = 1.0;

#endif //USE_TENSORFLOW

};

} //constraint_generators
} //trRosetta_protocols
} //protocols

#endif //INCLUDED_protocols_trRosetta_protocols_constraint_generators_trRosettaConstraintGenerator_hh
