// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/data/RNA_DMS_Potential.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_rna_data_RNA_DMS_Potential_HH
#define INCLUDED_core_scoring_rna_data_RNA_DMS_Potential_HH

#include <core/scoring/rna/data/RNA_DMS_Potential.fwd.hh>
#include <core/scoring/rna/data/RNA_DataInfo.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <numeric/MathNTensor.hh>
#include <numeric/MathTensor.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.hh>
#include <numeric/interpolation/spline/PolycubicSpline.tmpl.hh>
#include <numeric/interpolation/InterpolatedPotential.hh>
#include <numeric/interpolation/interpolation.hh>

#include <numeric/numeric.functions.hh>
#include <set>

#include <utility/fixedsizearray1.hh>


namespace core {
namespace scoring {
namespace rna {
namespace data {

class RNA_DMS_Potential: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_DMS_Potential();

	//destructor
	~RNA_DMS_Potential();

public:

	void
	initialize( core::pose::Pose const & pose );

	core::Real
	evaluate( core::pose::Pose const & pose,
		RNA_Reactivity const & rna_reactivity );

	Real
	get_binding_energy( Size const i,
		core::Vector const & probe_xyz,
		core::scoring::ScoreFunction const & scorefxn );

	Real
	get_N1_lonepair_donor_angle(  core::conformation::Residue const & acc_rsd,
		core::conformation::Residue const & don_rsd,
		Size const don_h_atm ) const;

	bool
	check_hbonded( pose::Pose const & pose,
		Size const & i /*residue number*/,
		std::string const & atom_name,
		bool is_acceptor ) const;

	bool
	check_chbonded( pose::Pose const & pose,
		Size const & i /*residue number*/,
		std::string const & atom_name ) const;

	core::Vector
	get_probe_xyz( core::conformation::Residue const & rsd, Distance const probe_dist  ) const;

	core::scoring::ScoreFunctionOP
	get_probe_scorefxn( bool const soft_rep, bool const just_atr_rep ) const;

	void
	get_occupancy_densities( utility::vector1< Real > & occupancy_densities,
		pose::Pose const & pose,
		Size const i /*for exclusion*/,
		core::Vector const & probe_xyz,
		utility::vector1< Distance > const & shells ) const;

	Real
	get_occupancy_density( pose::Pose const & pose,
		Size const i /*for exclusion*/,
		core::Vector const & probe_xyz,
		std::pair< Distance, Distance > const & shells ) const;

	utility::vector1< Real >
	get_logL_values( pose::Pose const & pose, Size const i /*, utility::vector1< Real > const & DMS_values */  );

	std::set< Real > const &
	DMS_values() const { return DMS_values_; }

private:

	void
	initialize_DMS_potential();

	numeric::MathNTensor< Real, 3 > // this is silly -- should use a grid object
	read_DMS_stats_file( std::string const & potential_file );

	void
	figure_out_potential();

	bool
	get_features( pose::Pose const & pose,
		Size const i,
		bool & ade_n1_bonded,
		Real & binding_energy,
		Real & occupancy_density );

	void
	add_probe_to_pose( pose::Pose & pose );

	void
	update_virtual_base_if_necessary( pose::Pose & pose, Size const i );

private:

	numeric::interpolation::InterpolatedPotential< 4 > interpolated_potential_;

	bool const separate_scores_;
	core::Distance occ_dist_, methyl_probe_dist_, oxygen_probe_dist_;
	std::pair< Distance, Distance > occ_shells_;

	// Values for the axes don't need to be mathmatrices or anything
	std::set< Real > is_bonded_values_;
	std::set< Real > occ_values_, binding_energy_values_, DMS_values_;

	// useful stats to hold on to.
	numeric::MathNTensor< Real, 4 > DMS_stats_;
	// 'model' = (is_bonded, occ, binding_energy).
	numeric::MathTensor< Real > p_model_;
	numeric::MathVector< Real > p_DMS_;

	// Briefly, DMS_potential = -kT log( p / ( p_DMS  p_model ) );
	// this is ridiculous -- there should be a universal grid object in Rosetta -- but will do the job for now:
	numeric::MathNTensor< Real, 4 > DMS_potential_;
	numeric::MathMatrix< Real > DMS_potential_is_bonded_;
	numeric::MathMatrix< Real > DMS_potential_occ_;
	numeric::MathMatrix< Real > DMS_potential_binding_energy_;

	// these are a bit silly -- currently using a very slow computation of binding energy
	// of a 'mock' DMS molecule, based on virtualizing residues, scoring poses, etc.
	core::pose::PoseOP working_pose_, working_pose_with_probe_;
	hbonds::HBondSetOP hbond_set_;
	core::scoring::ScoreFunctionOP probe_scorefxn_;
};

} //data
} //rna
} //scoring
} //core

#endif
