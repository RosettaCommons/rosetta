// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/protocols/scoring/methods/saxs/PDDFEnergy.hh
/// @brief  "Energy" based on a similarity of theoretical PDDF (pairwise distance distribution function)
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_protocols_scoring_methods_saxs_PDDFEnergy_hh
#define INCLUDED_protocols_scoring_methods_saxs_PDDFEnergy_hh

// Package headers
#include <protocols/scoring/methods/saxs/PDDFEnergy.fwd.hh>

#include <core/scoring/methods/WholeStructureEnergy.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/saxs/FormFactorManager.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace saxs {


class PDDFEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:

    PDDFEnergy();

    PDDFEnergy(utility::vector1<core::Real> const &,utility::vector1<core::Real> const &);

    virtual ~PDDFEnergy() {}

    virtual core::scoring::methods::EnergyMethodOP clone() const { return core::scoring::methods::EnergyMethodOP( new PDDFEnergy() ); }

    virtual void finalize_total_energy(core::pose::Pose & pose,core::scoring::ScoreFunction const &,core::scoring::EnergyMap & totals) const;

    virtual void indicate_required_context_graphs(utility::vector1< bool > & /*context_graphs_required*/
    ) const {}

	core::scoring::methods::EnergyMethodOP create_energy_method(core::scoring::methods::EnergyMethodOptions const &) const {
	return core::scoring::methods::EnergyMethodOP( new PDDFEnergy() );
    }

    utility::vector1<core::Real>&  get_pddf() {
	return pose_pddf_;
    }

    utility::vector1<core::Real>&  get_dist_bins() {
	return d_;
    }

    utility::vector1<core::Real> & compute_pddf(const core::pose::Pose &) const;
    utility::vector1<core::Real> & compute_pddf_without_ff(const core::pose::Pose &) const;
    core::Real compute_chi(utility::vector1<core::Real> const &, utility::vector1<core::Real> const &) const;
    core::Real compute_L1(utility::vector1<core::Real> const &, utility::vector1<core::Real> const &) const;
    void create_pddf(core::pose::Pose &,core::Real,core::Real,core::Real);

    core::Real evaluate_pddf_energy(const core::pose::Pose & pose) const;

private:
    mutable utility::vector1< utility::vector1<core::Real> > factors_;
    mutable utility::vector1<core::Size> r_ids_;
    mutable utility::vector1<core::Size> a_ids_;
    mutable utility::vector1< utility::vector1<core::Real> > dmatrix_;
    mutable utility::vector1<bool> is_glob_;

    core::Real norm_;
    bool if_fit_area_;
    utility::vector1<core::Real> d_;
    mutable utility::vector1<core::Real> pose_pddf_;
    utility::vector1<core::Real> reference_pddf_;
    core::scoring::saxs::FormFactorManager* ff_manager_;
    core::Real bin_size_;
    core::Size min_bin_;
    core::Size max_bin_;
    bool if_hydrogens_;

    void read_pddf(std::string);
virtual
core::Size version() const;
};


}
}
}
}

#endif
