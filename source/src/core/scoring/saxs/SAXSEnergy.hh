// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/SAXSEnergy.hh
/// @brief  "Energy" based on a similarity of theoretical SAXS spectrum computed for a pose and the experimental data
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_core_scoring_saxs_SAXSEnergy_hh
#define INCLUDED_core_scoring_saxs_SAXSEnergy_hh

// Package headers
#include <core/scoring/saxs/DistanceHistogram.hh>
#include <core/scoring/saxs/FormFactor.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorCEN.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorFA.hh>
#include <core/scoring/saxs/SAXSEnergyCreator.hh>


#include <core/scoring/methods/WholeStructureEnergy.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <string>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/saxs/FormFactorManager.fwd.hh>
#include <utility/vector1.hh>
#include <map>


namespace core {
namespace scoring {
namespace saxs {


class SAXSEnergy : public methods::WholeStructureEnergy  {
public:

    static std::string fa_cfg_file_;
    static std::string cen_cfg_file_;

    SAXSEnergy(std::string &,core::chemical::ResidueTypeSetCAP, ScoreType, methods::EnergyMethodCreatorOP);

    SAXSEnergy(const std::string &,const utility::vector1<Real> &,const utility::vector1<Real> &,
	    ScoreType, methods::EnergyMethodCreatorOP);

    virtual ~SAXSEnergy() {}

    virtual methods::EnergyMethodOP clone() const {

	if(saxs_score_variant_ == saxs_fa_score)
	    return methods::EnergyMethodOP( new SAXSEnergy(the_config_file_,q_, reference_intensities_,saxs_score_variant_, methods::EnergyMethodCreatorOP( new SAXSEnergyCreatorFA ) ) );
	else
	    return methods::EnergyMethodOP( new SAXSEnergy(the_config_file_,q_, reference_intensities_,saxs_score_variant_, methods::EnergyMethodCreatorOP( new SAXSEnergyCreatorCEN ) ) );
    }

    virtual void finalize_total_energy(pose::Pose & pose,ScoreFunction const &,EnergyMap & totals) const;

    virtual void indicate_required_context_graphs(utility::vector1< bool > & /*context_graphs_required*/
    ) const {}

    methods::EnergyMethodOP create_energy_method(methods::EnergyMethodOptions const &) const {

	if(saxs_score_variant_ == saxs_fa_score)
	    return methods::EnergyMethodOP( new SAXSEnergy(the_config_file_,q_, reference_intensities_,saxs_score_variant_, methods::EnergyMethodCreatorOP( new SAXSEnergyCreatorFA ) ) );
	else
	    return methods::EnergyMethodOP( new SAXSEnergy(the_config_file_,q_, reference_intensities_,saxs_score_variant_, methods::EnergyMethodCreatorOP( new SAXSEnergyCreatorCEN ) ) );
    }

    utility::vector1<Real>& get_reference_intensities() { return reference_intensities_; }
    utility::vector1<Real>& get_pose_intensities() { return pose_intensities_; }
    utility::vector1<Real>& get_q() { return q_; }
    Size count_scoring_atoms() { return ff_ops_.size(); }

    Real total_energy(const pose::Pose & pose) const;

    Real compute_zero_intensity() const;

protected:

    ScoreType saxs_score_variant_;
    std::string the_config_file_;


    void init_ff(const std::string &);

    /// @brief Sets a vector of q arguments to be used by this class in all SAXS-related calculations
    void set_up_q(const utility::vector1<Real> &);

    /// @brief Sets a vector of q arguments using the command-line flags to be used by this class in all SAXS-related calculations
    void set_up_q();

    /// @brief Reads two vectors from a file: q and I(q) (meant to be experimental values)
    void read_spectrum(std::string &,utility::vector1<Real> &,utility::vector1<Real> &) const;

    /// @brief Takes two vectors: q and I(q) and returns  a SAXS spectrum as a function of Q as defined in this class
    /// @details This is basically done by fiting  a spline and re-evaluating I(q) at different q values
    void fit_intensities(const utility::vector1<Real> &,const utility::vector1<Real> &,utility::vector1<Real> &) const;

    /// @brief A shortcut that calls read_spectrum() and then fit_intensities()
    void read_intensities(std::string &,utility::vector1<Real> &) const;

    /// @brief computes I(q) from a pose
    void compute_intensities(const core::pose::Pose &,utility::vector1<Real> &) const;

    /// @brief ce-calculates form factors 
    /// @details This is necessary when :
    ///    - this is the first time
    ///    - an amino acid sequence within the scored pose has been changed
    ///    - teh q-value set has been affected
    void rehash_form_factors(const core::pose::Pose & pose) const;

private:

    mutable utility::vector1<Size> atom_ff_types_;		// its size is nAtoms
    mutable utility::vector1<FormFactorOP> ff_ops_;		// its size is nUniqFFs
    mutable std::map<FormFactorOP,Size> ff_map_;		// its size is nUniqFFs
    mutable utility::vector1< utility::vector1< DistanceHistogramOP > > dhist_;	// its size is nUniqFFs x nUniqFFs
    mutable utility::vector1<Size> r_ids_;			// its size is nAtoms
    mutable utility::vector1<Size> a_ids_;			// its size is nAtoms

    bool if_hydrogens_;
    mutable Real zero_;
    Real score_bias_;
    utility::vector1<Real> q_;
    mutable utility::vector1<Real> pose_intensities_;
    mutable utility::vector1<Real> reference_intensities_;
    FormFactorManager* ff_manager_;

    Real compute_chi(utility::vector1<Real> const &, utility::vector1<Real> const &) const;
    Real compute_chi_with_fit(utility::vector1<Real> const &, utility::vector1<Real> const &) const;
    Real compute_L1(utility::vector1<Real> const &, utility::vector1<Real> const &) const;
    void compute_distance_histogram(const core::pose::Pose & pose) const;

    virtual core::Size version() const;
};


}
}
}

#endif
