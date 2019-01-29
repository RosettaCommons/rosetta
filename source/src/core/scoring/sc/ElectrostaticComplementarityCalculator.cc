// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/scoring/sc/ElectrostaticComplementarityCalculator.cc
/// @brief    Headers for the Electrostatic Complementarity Calculator
/// @details
/// @author   Brian Coventry (bcov@uw.edu)

#include <core/scoring/sc/ElectrostaticComplementarityCalculator.hh>


#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>


#include <numeric/statistics/functions.hh>
#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/FArray1D.hh>


static basic::Tracer TR( "core.scoring.sc.ElectrostaticComplementarityCalculator" );


namespace core {
namespace scoring {
namespace sc {


ElectrostaticComplementarityCalculator::ElectrostaticComplementarityCalculator()
{
	Reset();
}

ElectrostaticComplementarityCalculator::~ElectrostaticComplementarityCalculator() = default;


int
ElectrostaticComplementarityCalculator::Init( core::pose::Pose const & pose ) {

	molecule_atoms_.resize(2);

	pose::initialize_atomid_map( molecule_atoms_[0], pose, false );
	pose::initialize_atomid_map( molecule_atoms_[1], pose, false );

	ignore_radius_ = ElectrostaticComplementarityDefaults::IGNORE_RADIUS;
	interface_trim_radius_ = ElectrostaticComplementarityDefaults::INTERFACE_TRIM_RADIUS;
	partially_solvated_ = ElectrostaticComplementarityDefaults::PARTIALLY_SOLVATED;
	ignore_radius_resolution_ = ElectrostaticComplementarityDefaults::IGNORE_RADIUS_RESOLUTION;

	return sc_calc_.Init();
}

void
ElectrostaticComplementarityCalculator::Reset() {

	molecule_atoms_.clear();

	// Unfortunately, -1 is a valid output
	results_.ec_0_p = -2;
	results_.ec_0_s = -2;
	results_.ec_1_p = -2;
	results_.ec_1_s = -2;
	results_.polarity_mag_avg0 = 0;
	results_.polarity_avg0 = 0;
	results_.polarity_mag_avg1 = 0;
	results_.polarity_avg1 = 0;


	sc_calc_.Reset();
}

core::Size
ElectrostaticComplementarityCalculator::AddResidue( core::pose::Pose const & pose, int molecule, Size seqpos ) {

	if ( molecule_atoms_.size() == 0 ) return 0;

	conformation::Residue const & res = pose.residue( seqpos );

	for ( Size ia = 1; ia <= res.natoms(); ia++ ) {

		if ( res.is_virtual( ia ) ) continue;

		molecule_atoms_[molecule]( seqpos, ia ) = true;
	}

	return sc_calc_.AddResidue( molecule, res );
}

int
ElectrostaticComplementarityCalculator::Calc( core::pose::Pose const & pose, core::Size jump_id )
{
	if ( molecule_atoms_.size() == 0 ) return 0;

	if ( jump_id > pose.num_jump() ) {
		TR.Error << "Jump ID out of bounds (pose has " << pose.num_jump() << " jumps)" << std::endl;
		return 0;
	}

	// Partition pose by jump_id
	ObjexxFCL::FArray1D_bool is_upstream ( pose.size(), true );

	if ( jump_id > 0 ) {
		pose.fold_tree().partition_by_jump( jump_id, is_upstream );
	}

	for ( Size i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue const & residue = pose.residue(i);
		if ( residue.type().name() == "VRT" ) {
			continue;
		}
		if ( !AddResidue( pose, is_upstream(i) ? 0 : 1, i ) ) {
			return 0;
		}
	}

	return Calc( pose );
}


int ElectrostaticComplementarityCalculator::Calc( core::pose::Pose const & pose )
{
	if ( molecule_atoms_.size() == 0 ) return 0;

	sc_calc_.settings.band = interface_trim_radius_;
	// Do this first because it's fast compared to PB
	int sc_success = sc_calc_.Calc();

	if ( ! sc_success ) {
		TR.Error << "Failed because internal ShapeComplementarityCalculator failed" << std::endl;
		return 0;
	}

	id::AtomID_Map<bool> atoms_within_radius;
	if ( ignore_radius_ > 0 ) {

		TR << "Ignoring atoms farther than " << ignore_radius_ << " from interface." << std::endl;
		std::vector<const DOT*> all_dots = sc_calc_.GetTrimmedDots(0);
		all_dots.insert( all_dots.end(), sc_calc_.GetTrimmedDots(1).begin(), sc_calc_.GetTrimmedDots(1).end() );
		std::vector<const DOT*> sparse_dots = trim_dots_to_resl( all_dots );

		TR << "Trimmed " << all_dots.size() << " dots to " << sparse_dots.size() << " dots." << std::endl;

		get_atoms_within_radius( pose, atoms_within_radius, sparse_dots, molecule_atoms_[0], molecule_atoms_[1] );
	} else {

		pose::initialize_atomid_map( atoms_within_radius, pose, true );
	}

	// The check to see if the file exists doesn't work here because we don't know the extension
	// Just use "." even if there's a temp folder and hope we don't get lucky
	std::string temp_tag = utility::file_basename( utility::file::create_temp_filename( ".", "PB_" ) );

	PoissonBoltzmannPotential pb_of_0;
	{
		id::AtomID_Map<bool> present_atoms = get_present_atoms( pose, molecule_atoms_[0], molecule_atoms_[1], atoms_within_radius );
		pb_of_0.solve_pb( pose, temp_tag, molecule_atoms_[0], present_atoms, true );
	}

	PoissonBoltzmannPotential pb_of_1;
	{
		id::AtomID_Map<bool> present_atoms = get_present_atoms( pose, molecule_atoms_[1], molecule_atoms_[0], atoms_within_radius );
		pb_of_1.solve_pb( pose, temp_tag, molecule_atoms_[1], present_atoms, true );
	}

	// Calculate the actual electrostatic complementarities

	utility::vector1<Real> charge_of_0_on_0;
	utility::vector1<Real> charge_of_1_on_0;
	prepare_correlation_lists( sc_calc_.GetTrimmedDots(0), pb_of_0, pb_of_1, charge_of_0_on_0, charge_of_1_on_0 );

	results_.ec_0_p = -numeric::statistics::corrcoef( charge_of_0_on_0, charge_of_1_on_0 );
	results_.ec_0_s = -numeric::statistics::spearman_r( charge_of_0_on_0, charge_of_1_on_0 );

	utility::vector1<Real> charge_of_0_on_1;
	utility::vector1<Real> charge_of_1_on_1;
	prepare_correlation_lists( sc_calc_.GetTrimmedDots(1), pb_of_0, pb_of_1, charge_of_0_on_1, charge_of_1_on_1 );

	results_.ec_1_p = -numeric::statistics::corrcoef( charge_of_0_on_1, charge_of_1_on_1 );
	results_.ec_1_s = -numeric::statistics::spearman_r( charge_of_0_on_0, charge_of_1_on_0 );

	// Do some averaging that bcov made up.

	{
		utility::vector1<Real> all_charge_of_0 = charge_of_0_on_0;
		all_charge_of_0.insert( all_charge_of_0.end(), charge_of_0_on_1.begin(), charge_of_0_on_1.end() );

		utility::vector1<Real> all_charge_of_1 = charge_of_1_on_0;
		all_charge_of_1.insert( all_charge_of_1.end(), charge_of_1_on_1.begin(), charge_of_1_on_1.end() );


		results_.polarity_mag_avg0 = 0;
		results_.polarity_avg0 = 0;
		results_.polarity_mag_avg1 = 0;
		results_.polarity_avg1 = 0;

		for ( Size i = 1; i <= all_charge_of_0.size(); i++ ) {
			results_.polarity_mag_avg0 += std::abs(all_charge_of_0[i]);
			results_.polarity_avg0 += all_charge_of_0[i];
			results_.polarity_mag_avg1 += std::abs(all_charge_of_1.at(i));
			results_.polarity_avg1 += all_charge_of_1.at(i);
		}

		results_.polarity_mag_avg0 /= all_charge_of_0.size();
		results_.polarity_avg0 /= all_charge_of_0.size();
		results_.polarity_mag_avg1 /= all_charge_of_0.size();
		results_.polarity_avg1 /= all_charge_of_0.size();
	}

	return 1;
}


void
ElectrostaticComplementarityCalculator::prepare_correlation_lists(
	std::vector<const DOT*> const & dots,
	PoissonBoltzmannPotential const & pb_a,
	PoissonBoltzmannPotential const & pb_b,
	utility::vector1<Real> & charges_from_a,
	utility::vector1<Real> & charges_from_b ) const {

	for ( auto dot_p : dots ) {

		charges_from_a.push_back( pb_a.get_potential( dot_p->coor ) );
		charges_from_b.push_back( pb_b.get_potential( dot_p->coor ) );
	}
}


id::AtomID_Map<bool>
ElectrostaticComplementarityCalculator::get_present_atoms(
	core::pose::Pose const & pose,
	id::AtomID_Map<bool> const & charged_side,
	id::AtomID_Map<bool> const & uncharged_side,
	id::AtomID_Map<bool> const & atoms_within_radius
) const {

	id::AtomID_Map<bool> present_atoms;
	pose::initialize_atomid_map( present_atoms, pose, false );

	for ( Size ires = 1; ires <= uncharged_side.n_residue(); ires++ ) {
		for ( Size iatom = 1; iatom <= uncharged_side.n_atom(ires); iatom++ ) {

			if ( ! atoms_within_radius( ires, iatom ) ) continue;

			if ( charged_side( ires, iatom ) ) {
				present_atoms( ires, iatom ) = true;
			}

			if ( partially_solvated_ && uncharged_side( ires, iatom ) ) {
				present_atoms( ires, iatom ) = true;
			}
		}
	}

	return present_atoms;
}



bool
ElectrostaticComplementarityCalculator::is_within_radius(
	numeric::xyzVector< Real > xyz,
	std::vector<const DOT*> const & dots
) const {

	if ( dots.size() == 0 ) return true;

	Real radius = ignore_radius_ + ignore_radius_resolution_;
	Real radius2 = radius*radius;

	for ( auto dotp : dots ) {
		if ( dotp->coor.distance_squared( xyz ) < radius2 ) return true;
	}

	return false;
}

/// @brief Figure out which atoms are within the radius and part of one of the molecules
void
ElectrostaticComplementarityCalculator::get_atoms_within_radius(
	core::pose::Pose const & pose,
	id::AtomID_Map<bool> & atoms_within_radius,
	std::vector<const DOT*> const & dots,
	id::AtomID_Map<bool> const & molecule1,
	id::AtomID_Map<bool> const & molecule2
) const {

	pose::initialize_atomid_map( atoms_within_radius, pose, false );

	for ( Size ires = 1; ires <= atoms_within_radius.n_residue(); ires++ ) {
		for ( Size iatom = 1; iatom <= atoms_within_radius.n_atom(ires); iatom++ ) {

			if ( ! ( molecule1( ires, iatom ) || molecule2( ires, iatom ) ) ) continue;

			atoms_within_radius( ires, iatom ) = is_within_radius( pose.residue(ires).xyz( iatom ), dots );
		}
	}
}

/// @brief Do something similar to k-means to reduce the number of dots
std::vector<const DOT*>
ElectrostaticComplementarityCalculator::trim_dots_to_resl( std::vector<const DOT*> const & dots ) const {

	std::vector<const DOT*> good_dots;

	Real resl = ignore_radius_resolution_;
	Real resl2 = resl*resl;

	for ( auto dotp : dots ) {
		bool is_unique = true;

		for ( auto good_dotp : good_dots ) {
			if ( good_dotp->coor.distance_squared( dotp->coor ) < resl2 ) {
				is_unique = false;
				break;
			}
		}

		if ( is_unique ) {
			good_dots.push_back( dotp );
		}
	}

	return good_dots;
}




} // namespace sc
} // namespace scoring
} // namespace core




