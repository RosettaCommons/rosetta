// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/NMerSVMEnergy.hh
/// @brief  SVM energy method declaration
/// @author Chris King (dr.chris.king@gmail.com)


#ifndef INCLUDED_core_scoring_methods_NMerSVMEnergy_hh
#define INCLUDED_core_scoring_methods_NMerSVMEnergy_hh

// Unit headers
#include <core/scoring/methods/NMerSVMEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/NMerPSSMEnergy.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/libsvm/Svm_rosetta.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <map>

// Utility headers


namespace core {
namespace scoring {
namespace methods {

class NMerSVMEnergy : public ContextIndependentOneBodyEnergy
{
public:
	typedef ContextIndependentOneBodyEnergy parent;

public:


	NMerSVMEnergy();


	NMerSVMEnergy( utility::vector1< std::string > const & svm_fnames );


	virtual ~NMerSVMEnergy();

	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }


	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;

	/// @brief context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	core::Real
	get_residue_energy_by_svm(
		core::pose::Pose const &,
		core::Size const &,
		core::Real &,
		utility::vector1< core::Real > &
	) const;
	void read_nmer_svm_list( std::string const );
	void read_nmer_svm( std::string const );
	core::Size n_svms() const;
	void read_aa_encoding_matrix( std::string const );
	void nmer_length( core::Size const );
	core::Size nmer_length() const;
	void term_length( core::Size const );
	core::Size term_length() const;
	void use_pssm_features( bool const );
	void gate_svm_scores( bool const );
	void nmer_svm_scorecut( core::Real const );
	utility::vector1< utility::libsvm::Svm_node_rosettaOP > get_svm_nodes( utility::vector1< core::Real > const ) const;
	utility::vector1< core::Real > encode_aa_string( std::string const ) const;
	utility::vector1< core::Real > encode_wtd_avg_aa_string( std::string const, utility::vector1< core::Real > const ) const;
	utility::vector1< core::Real > encode_nmer( std::string const, core::Size const, core::Size const ) const;
	void add_encoded_termini( std::string const, core::Size const, utility::vector1< core::Real > & ) const;
	void add_pssm_features( std::string const, core::Size const, utility::vector1< core::Real > & ) const;

private:
	utility::vector1< utility::libsvm::Svm_rosettaOP > all_nmer_svms_;
	std::map< char, utility::vector1< core::Real > > aa_encoder_;
	core::Size nmer_length_;
	core::Size nmer_cterm_;
	core::Size term_length_;
	bool gate_svm_scores_;
	bool use_pssm_features_;
	core::Real nmer_svm_scorecut_;
	core::scoring::methods::NMerPSSMEnergy nmer_pssm_;

	void read_nmer_svms_from_options();
	void initialize_from_options();
	virtual
	core::Size version() const;
};

} // methods
} // scoring
} // core


#endif
