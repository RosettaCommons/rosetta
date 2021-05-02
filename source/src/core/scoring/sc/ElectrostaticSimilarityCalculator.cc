// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/scoring/sc/ElectrostaticSimilarityCalculator.cc
/// @brief    Headers for the Electrostatic Similarity Calculator
/// @details  The code closely follows Brian Coventry's implementation
///    of electrostatic complementarity with a added features
///    that in turn is based on the method:
///    McCoy, A. J., Epa, V. C., & Colman, P. M. (1997).
///    Electrostatic complementarity at protein/protein interfaces.
///    Journal of molecular biology, 268(2), 570-584.
/// @author   Andreas Scheck (andreas.scheck@epfl.ch)

// Project Headers
#include <core/scoring/sc/ElectrostaticSimilarityCalculator.hh>
#include <core/select/residue_selector/ResidueVector.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/init_id_map.hh>
#include <core/conformation/Residue.hh>

#include <numeric/statistics/functions.hh>
#include <ObjexxFCL/FArray1D.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <utility/file/file_sys_util.hh>


static basic::Tracer TR( "core.scoring.sc.ElectrostaticSimilarityCalculator" );


namespace core {
namespace scoring {
namespace sc {


ElectrostaticSimilarityCalculator::ElectrostaticSimilarityCalculator()
{
	Reset();
}

ElectrostaticSimilarityCalculator::~ElectrostaticSimilarityCalculator() = default;


int
ElectrostaticSimilarityCalculator::Init() {

	partially_solvated_ = ElectrostaticSimilarityDefaults::PARTIALLY_SOLVATED;
	reference_pose_= nullptr;
	selector1_ =  nullptr;
	selector2_ = nullptr;

	return sc_calc_.Init();
}

void
ElectrostaticSimilarityCalculator::Reset() {

	molecule_atoms_.clear();

	// -1 is a valid output
	results_.es_0_p = -2;
	results_.es_0_s = -2;
	results_.es_1_p = -2;
	results_.es_1_s = -2;
	results_.polarity_mag_avg0 = 0;
	results_.polarity_avg0 = 0;
	results_.polarity_mag_avg1 = 0;
	results_.polarity_avg1 = 0;

	sc_calc_.Reset();
}

/// @brief Add a residue to the molecular surface that is being evaluated
/// @details
core::Size
ElectrostaticSimilarityCalculator::AddResidue( core::pose::Pose const & pose, int molecule, Size seqpos ) {

	if ( molecule_atoms_.size() == 0 ) return 0;

	conformation::Residue const & res = pose.residue( seqpos );

	for ( Size ia = 1; ia <= res.natoms(); ia++ ) {
		if ( res.is_virtual( ia ) ) continue;
		molecule_atoms_[molecule]( seqpos, ia ) = true;
	}

	return sc_calc_.AddResidue( molecule, res );
}

/// @brief Run the ES calculation for previously defined molecules (via AddResidue)
/// @details
/// This is a function calcualtes the electrostatic similarity for the given pose in respect
/// to a provided reference pose. Residue selectors for the pose and reference pose are used
/// to specify the evaluated surface. The residues are added with AddResdiue().
int ElectrostaticSimilarityCalculator::Calc( core::pose::Pose const & pose ) {
	using core::select::residue_selector::ResidueVector;

	ShapeSimilarityCalculator ssc;
	ShapeSimilarityCalculator ssc_full;
	ShapeSimilarityCalculator ssc_ref;
	ShapeSimilarityCalculator ssc_ref_full;

	ResidueVector const residues_pose( selector1_->apply( pose ) );
	ResidueVector const residues_native( selector2_->apply( *reference_pose_ ) );

	// this assumes all residues are from the same chain
	Size pose_chain = pose.residue(residues_pose[1]).chain();
	Size native_chain = reference_pose_->residue(residues_native[1]).chain();

	TR << "chain pose: " << pose_chain << " and chain native: " << native_chain << std::endl;



	// Dump information about residues

	TR << "Using residues for molecule surface (rosetta numbering):" << std::endl;
	TR << "  Scaffold: ";
	for ( auto r = residues_pose.begin(); r != residues_pose.end(); ++r ) {
		TR << (r == residues_pose.begin() ? "" : ", ") << *r;
	}
	TR << std::endl;
	TR << "  Native: ";
	for ( auto r = residues_native.begin(); r != residues_native.end(); ++r ) {
		TR << (r == residues_native.begin() ? "" : ", ") << *r;
	}
	TR << std::endl;


	for ( core::Size r : residues_pose ) {
		ssc.AddResidue( 0, pose.residue( r ) );
	}
	core::pose::PoseOP pose_subset = pose.split_by_chain(pose_chain);
	for ( auto r = pose_subset->begin(); r != pose_subset->end(); ++r ) {
		ssc_full.AddResidue( 0, *r);
	}

	for ( core::Size r : residues_native ) {
		ssc_ref.AddResidue( 0, reference_pose_->residue( r ) );
	}
	core::pose::PoseOP native_subset = reference_pose_->split_by_chain(native_chain);
	for ( auto r = native_subset->begin(); r != native_subset->end(); ++r ) {
		ssc_ref_full.AddResidue( 0, *r);
	}


	if ( !ssc.Calc() ) {
		return -1;
	}
	if ( !ssc_full.Calc() ) {
		return -1;
	}
	if ( !ssc_ref.Calc() ) {
		return -1;
	}
	if ( !ssc_ref_full.Calc() ) {
		return -1;
	}


	// To get all dots from the selected surface, we first compute the overall surface of the protein
	// followed by the surface of the selected residues, for the target and reference respectively.
	// The two dot clouds of the full surface and selected residues's surface are compared and only
	// dots that intersect are kept. This is to avoid dots that are an artefact of the surface
	// generation of selected residues and are actually buried in the protein core.

	std::vector<DOT const *> intersecting_dots[1];
	std::vector<DOT const *> intersecting_dots_native[1];

	for ( auto r1 = ssc.GetDots(0).begin(); r1 != ssc.GetDots(0).end(); ++r1 ) {
		for ( auto r2 = ssc_full.GetDots(0).begin(); r2 != ssc_full.GetDots(0).end(); ++r2 ) {
			core::Real dist = r1->coor.distance(r2->coor);
			if ( dist == 0 ) {
				DOT const &dot1 = *r1;
				intersecting_dots[0].push_back(&dot1);
				break;
			}
		}
	}

	for ( auto r1 = ssc_ref.GetDots(0).begin(); r1 != ssc_ref.GetDots(0).end(); ++r1 ) {
		for ( auto r2 = ssc_ref_full.GetDots(0).begin(); r2 != ssc_ref_full.GetDots(0).end(); ++r2 ) {
			core::Real dist = r1->coor.distance(r2->coor);
			if ( dist == 0 ) {
				DOT const &dot2 = *r1;
				intersecting_dots_native[0].push_back(&dot2);
				break;
			}
		}
	}

	sort( intersecting_dots[0].begin(), intersecting_dots[0].end() );
	sort( intersecting_dots_native[0].begin(), intersecting_dots_native[0].end() );
	intersecting_dots[0].erase( std::unique( intersecting_dots[0].begin(), intersecting_dots[0].end() ), intersecting_dots[0].end() );
	intersecting_dots_native[0].erase( std::unique( intersecting_dots_native[0].begin(), intersecting_dots_native[0].end() ), intersecting_dots_native[0].end() );


	id::AtomID_Map<bool> atoms_within_radius;
	id::AtomID_Map<bool> atoms_within_radius_native;
	pose::initialize_atomid_map( atoms_within_radius, pose, true );
	pose::initialize_atomid_map( atoms_within_radius_native, *reference_pose_, true );


	/////////////////////////////////////////////////////////////////
	/////////// Calculate electrostatics and correlations ///////////
	/////////////////////////////////////////////////////////////////

	// The check to see if the file exists doesn't work here because we don't know the extension
	// Just use "." even if there's a temp folder and hope we don't get lucky
	std::string temp_tag = utility::file_basename( utility::file::create_temp_filename( ".", "PB_" ) );


	molecule_atoms_.resize(4);

	pose::initialize_atomid_map( molecule_atoms_[0], pose, false );
	pose::initialize_atomid_map( molecule_atoms_[1], pose, false );
	pose::initialize_atomid_map( molecule_atoms_[2], *reference_pose_, false );
	pose::initialize_atomid_map( molecule_atoms_[3], *reference_pose_, false );

	for ( Size i = 1; i <= pose.size(); ++i ) {
		conformation::Residue const & res = pose.residue( i );
		for ( Size ia = 1; ia <= res.natoms(); ia++ ) {
			if ( res.is_virtual( ia ) ) continue;
			if ( res.chain() == pose_chain ) {
				molecule_atoms_[0]( i, ia ) = true;
			} else {
				molecule_atoms_[1]( i, ia ) = true;
			}
		}
	}

	for ( Size i = 1; i <= reference_pose_->size(); ++i ) {
		conformation::Residue const & res = reference_pose_->residue( i );
		for ( Size ia = 1; ia <= res.natoms(); ia++ ) {
			if ( res.is_virtual( ia ) ) continue;
			if ( res.chain() == native_chain ) {
				molecule_atoms_[2]( i, ia ) = true;
			} else {
				molecule_atoms_[3]( i, ia ) = true;

			}
		}
	}

	PoissonBoltzmannPotential pb_of_0;
	{
		id::AtomID_Map<bool> present_atoms = get_present_atoms( pose, molecule_atoms_[0], molecule_atoms_[1], atoms_within_radius );
		pb_of_0.solve_pb( pose, temp_tag, molecule_atoms_[0], present_atoms, true );
	}

	PoissonBoltzmannPotential pb_of_1;
	{
		id::AtomID_Map<bool> present_atoms_native = get_present_atoms( *reference_pose_, molecule_atoms_[2], molecule_atoms_[3], atoms_within_radius_native );
		pb_of_1.solve_pb( *reference_pose_, temp_tag, molecule_atoms_[2], present_atoms_native, true );
	}


	// Calculate the actual electrostatic similarities

	utility::vector1<Real> charge_of_0_on_0;
	utility::vector1<Real> charge_of_1_on_0;
	prepare_correlation_lists( intersecting_dots[0], pb_of_0, pb_of_1, charge_of_0_on_0, charge_of_1_on_0 );

	results_.es_0_p = -numeric::statistics::corrcoef( charge_of_0_on_0, charge_of_1_on_0 );
	results_.es_0_s = -numeric::statistics::spearman_r( charge_of_0_on_0, charge_of_1_on_0 );


	utility::vector1<Real> charge_of_0_on_1;
	utility::vector1<Real> charge_of_1_on_1;
	prepare_correlation_lists( intersecting_dots_native[0], pb_of_0, pb_of_1, charge_of_0_on_1, charge_of_1_on_1 );

	results_.es_1_p = -numeric::statistics::corrcoef( charge_of_0_on_1, charge_of_1_on_1 );
	results_.es_1_s = -numeric::statistics::spearman_r( charge_of_0_on_1, charge_of_1_on_1 );


	TR << "ES: " << ( ( ( (results_.es_0_p + results_.es_1_p) / 2 ) + ( (results_.es_0_s + results_.es_1_s) / 2 ) ) / 2) << std::endl;

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

/// @brief Prepare lists containing the point clouds for correlation analysis
void
ElectrostaticSimilarityCalculator::prepare_correlation_lists(
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

/// @brief Get atoms of the selected surface patch
id::AtomID_Map<bool>
ElectrostaticSimilarityCalculator::get_present_atoms(
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


} // namespace sc
} // namespace scoring
} // namespace core
