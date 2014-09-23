// INCLUDED_core_scoring_ScoreFunction_HH
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_Mg_Energy.cc
/// @brief  Statistically derived Mg(2+) binding potential for RNA.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_Mg_Energy.hh>
#include <core/scoring/rna/RNA_Mg_EnergyCreator.hh>
#include <core/scoring/rna/RNA_Mg_KnowledgeBasedPotential.hh>

// Package headers
#include <core/chemical/AtomType.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rescore.OptionKeys.gen.hh>

#include <numeric/conversions.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

using ObjexxFCL::format::I;
using namespace core::chemical::rna;

static thread_local basic::Tracer tr( "core.scoring.rna.RNA_Mg_Energy" );

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_Mg_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_Mg_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RNA_Mg_Energy;
}

ScoreTypes
RNA_Mg_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_mg );
	sts.push_back( rna_mg_indirect );
	return sts;
}


RNA_Mg_Energy::RNA_Mg_Energy() :
	parent( methods::EnergyMethodCreatorOP( new RNA_Mg_EnergyCreator ) ),
	rna_mg_knowledge_based_potential_( new RNA_Mg_KnowledgeBasedPotential ),
	hbond_options_( new core::scoring::hbonds::HBondOptions ), // useful helper functions
	verbose_( basic::options::option[ basic::options::OptionKeys::rescore::verbose ]() )
{}


/// clone
methods::EnergyMethodOP
RNA_Mg_Energy::clone() const
{
	return new RNA_Mg_Energy;
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Mg_Energy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	rna_mg_knowledge_based_potential_->setup_info_for_mg_calculation( pose );
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Mg_Energy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	rna_mg_knowledge_based_potential_->setup_info_for_mg_calculation( pose );
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Mg_Energy::setup_for_packing( pose::Pose & pose,	utility::vector1< bool > const &, utility::vector1< bool > const &  ) const
{
	pose.update_residue_neighbors();
	rna_mg_knowledge_based_potential_->setup_info_for_mg_calculation( pose );
}



//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Mg_Energy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

	rna::RNA_ScoringInfo const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	utility::vector1< bool > const & is_magnesium = rna_scoring_info.is_magnesium();

	if ( rsd1.is_RNA() && is_magnesium[ rsd2.seqpos() ] ) {
		residue_pair_energy_one_way( rsd1, rsd2, pose, emap );
	}

	if ( rsd2.is_RNA() && is_magnesium[ rsd1.seqpos() ] ) {
		residue_pair_energy_one_way( rsd2, rsd1, pose, emap );
	}

	return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_Mg_Energy::get_cos_theta( core::conformation::Residue const & rsd1,
	   Size const i, Vector const & i_xyz, Vector const & j_xyz ) const
{

  Vector dummy, xyz_base;

  if ( rsd1.heavyatom_is_an_acceptor( i ) ){
    chemical::Hybridization acc_hybrid( rsd1.atom_type( i ).hybridization() );
    make_hbBasetoAcc_unitvector(
				*hbond_options_,
				acc_hybrid,
				rsd1.atom( i ).xyz(),
				rsd1.xyz( rsd1.atom_base( i ) ),
				rsd1.xyz( rsd1.abase2( i ) ),
			      xyz_base, dummy );
  } else if ( rsd1.atom_type( i ).name() == "Hpol" ) {
    xyz_base = rsd1.xyz( rsd1.atom_base( i ) );
  } else {
    return -999; // bogus value
  }

  Vector const a = j_xyz - i_xyz;
  Vector const b = xyz_base - i_xyz;

  Real cos_theta = dot( a, b )/( a.length() * b.length() );

  //  if ( a.length() < 3.0 ) std::cout << rsd1.name1() << " " << rsd1.seqpos() << " " << rsd1.atom_name(i) << " " << cos_theta << std::endl;

  return cos_theta;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Mg_Energy::residue_pair_energy_one_way(
					   conformation::Residue const & rsd1, // The RNA residue
					   conformation::Residue const & rsd2, // The Mg(2+)
					   pose::Pose const & pose,
					   EnergyMap & emap
					   ) const{

  rna::RNA_ScoringInfo const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
  utility::vector1< bool > const & is_magnesium = rna_scoring_info.is_magnesium();

  runtime_assert( rsd1.is_RNA() );
  runtime_assert( is_magnesium[ rsd2.seqpos() ] );

  Real score( 0.0 ), score_indirect( 0.0 );

  // get magnesium position
  Size const j = 1;  //First atom of Mg2+ residue is assumed to be Mg2+ atom.
  runtime_assert( rsd2.atom_name( j ) ==  "MG  " );
  if ( rsd2.is_virtual( j ) ) return;

  Vector const & j_xyz = rsd2.xyz( j );

  // Loop over potential ligand positions.
  Size const pos1 = rsd1.seqpos();

  // These should be filled in during a pre-scoring step (setup_info_for_mg_calculation ).
  utility::vector1< utility::vector1< Size > > const
    atom_numbers_for_mg_calculation( rna_scoring_info.atom_numbers_for_mg_calculation() );
  utility::vector1< Size > const & atom_numbers1   ( atom_numbers_for_mg_calculation[ pos1 ]  );


  // apply angle cut.
  bool const apply_angle_potential_( true );

  //Distance const direct_interaction_cutoff_( 4.0 ), indirect_interaction_cutoff_( 8.0 );
  Distance const direct_interaction_cutoff_( 999.0 ), indirect_interaction_cutoff_( 999.0 );

  // direct interactions, including special non-pair-wise handling of phosphate oxygens
  utility::vector1< Real > phosphate_scores;
  bool is_phosphate_oxygen( false );
  for ( Size m = 1; m <= atom_numbers1.size(); ++m ) {

    Size const i = atom_numbers1[ m ];
    if ( rsd1.is_virtual( i ) ) continue;
    Vector const & i_xyz( rsd1.xyz( i ) );

    Distance d = ( i_xyz - j_xyz ).length();

    GaussianParameter const & mg_potential_gaussian_parameter = rna_mg_knowledge_based_potential_->get_mg_potential_gaussian_parameter( rsd1, i, is_phosphate_oxygen );

    Real cos_theta( -999.0 );
    bool get_angle_form_factor = false;
    if ( apply_angle_potential_ && ( rsd1.heavyatom_is_an_acceptor( i ) || rsd1.atom_type( i ).name() == "Hpol" ) ) get_angle_form_factor = true;

    if ( d < direct_interaction_cutoff_  && mg_potential_gaussian_parameter.center > 0.0 ){

      Real binding_score = get_gaussian_potential_score( mg_potential_gaussian_parameter, i_xyz, j_xyz );

      if ( get_angle_form_factor ){
	cos_theta = get_cos_theta( rsd1, i, i_xyz, j_xyz );
	GaussianParameter const & mg_potential_costheta_gaussian_parameter = rna_mg_knowledge_based_potential_->get_mg_potential_costheta_gaussian_parameter( rsd1, i );
	Real const angle_potential = get_gaussian_score( mg_potential_costheta_gaussian_parameter, cos_theta );
	binding_score *= angle_potential;
      }

      if ( verbose_ && std::abs( binding_score ) >  0.1 ) tr <<  "Mg " << rsd2.seqpos() << "   direct to ligand " << pos1 << ' ' << rsd1.atom_name( i )  << "   cos_angle " << cos_theta << "  score: " <<  binding_score  << std::endl;

      if ( is_phosphate_oxygen ) {
	phosphate_scores.push_back( binding_score ); // will get added in later
      } else {
	score += binding_score;
      }
    }


    // indirect interactions, simple pair-wise
    GaussianParameter const & mg_potential_indirect_gaussian_parameter = rna_mg_knowledge_based_potential_->get_mg_potential_indirect_gaussian_parameter( rsd1, i );

    if ( d < indirect_interaction_cutoff_ && mg_potential_indirect_gaussian_parameter.center > 0.0 ){

      if ( mg_potential_indirect_gaussian_parameter.center > 0.0 ){
	Real binding_score_indirect = get_gaussian_potential_score( mg_potential_indirect_gaussian_parameter, i_xyz, j_xyz );

	if ( get_angle_form_factor ){

	  if ( cos_theta < -1.0 ) cos_theta = get_cos_theta( rsd1, i, i_xyz, j_xyz );

	  GaussianParameter const & mg_potential_costheta_indirect_gaussian_parameter = rna_mg_knowledge_based_potential_->get_mg_potential_costheta_indirect_gaussian_parameter( rsd1, i );
	  Real const angle_potential_indirect = get_gaussian_score( mg_potential_costheta_indirect_gaussian_parameter, cos_theta );
	  binding_score_indirect *= angle_potential_indirect;

	}
	score_indirect += binding_score_indirect;

	if ( verbose_ && std::abs( binding_score_indirect ) >  0.1 ) tr <<  "Mg " << rsd2.seqpos() << " indirect to ligand " << pos1 << ' ' << rsd1.atom_name( i ) << "  score: " <<  binding_score_indirect << std::endl;

      }
    }


  }

  // there appear to be very few (if any) metal ions that bind to both the OP2 and OP1 of a single phosphate -- so let's use the best of those scores.
  if ( phosphate_scores.size() == 1 ) phosphate_scores.push_back( 0.0 );

  if ( phosphate_scores.size() > 1 ){
    runtime_assert( phosphate_scores.size() == 2 ); // OP2 and OP1
    score += std::min( phosphate_scores[1], phosphate_scores[2] );
  }

  emap[ rna_mg ] += score;

  emap[ rna_mg_indirect ] += score_indirect;


}


/////////////////////////////////
Real
RNA_Mg_Energy::get_gaussian_potential_score(
					    GaussianParameter const & mg_potential_gaussian_parameter,
					    Vector const & pos1,
					    Vector const & pos2 ) const
{ // later expand to do derivative calculation

  Distance const d = ( pos1 - pos2 ).length();

  return get_gaussian_score( mg_potential_gaussian_parameter, d );
}

///////////////////////////////////////////////////
Real
RNA_Mg_Energy::get_gaussian_score(
					    GaussianParameter const & mg_potential_gaussian_parameter,
					    Real const d ) const
{ // later expand to do derivative calculation

  Real const a     = mg_potential_gaussian_parameter.amplitude;
  Real const d0    = mg_potential_gaussian_parameter.center;
  Real const sigma = mg_potential_gaussian_parameter.width;

  Real const score = a * exp( -0.5 * std::pow( ( d - d0 )/sigma, 2 ) );
  return score;
}

/////////////////////////////////
void
RNA_Mg_Energy::eval_atom_derivative(
	id::AtomID const & /*atom_id*/,
	pose::Pose const & /*pose*/,
	kinematics::DomainMap const & /*domain_map*/,
	ScoreFunction const &,
	EnergyMap const & /*weights*/,
	Vector & /*F1*/,
	Vector & /*F2*/
 	) const
{
	// NEED TO FILL THIS IN LATER!!
}



/// @brief RNA_Mg_Energy distance cutoff
Distance
RNA_Mg_Energy::atomic_interaction_cutoff() const
{
	return 6.0;
}

/// @brief RNA_Mg_Energy
void
RNA_Mg_Energy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}

core::Size
RNA_Mg_Energy::version() const
{
	return 1; // Initial versioning
}


} //rna
} //scoring
} //core
