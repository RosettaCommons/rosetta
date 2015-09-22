// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/rna_rotamer_building_functions.cc
/// @brief  RNA nucleotide rotamer set class implementation
/// @author Rhiju Das (rhiju@stanford.edu), Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/rotamer_set/rotamer_building_functions.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/rna/RNA_ResidueLevelTask.hh>
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.hh>

// Project Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/rna/RNA_SamplerUtil.hh>

#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/Stub.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/database/open.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray3D.hh>

// C++ headers
#include <string>
#include <iostream>
#include <fstream>


namespace core {
namespace pack {
namespace rotamer_set {

using namespace core::conformation;

static THREAD_LOCAL basic::Tracer TR( "core.pack.rotamer_set.rna_rotamer_building_functions", basic::t_info );

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RNA-specific methods, pulled out of rna_rotamer_building_functions.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
fill_chi_rotamers_with_center_and_stddev( conformation::ResidueOP const & rot,
	utility::vector1< conformation::ResidueOP > & rotamers,
	utility::vector1< Real > const & chi_steps,
	Real const & center, Real const & width )
{
	using namespace conformation;
	for ( Size n = 1 ; n <= chi_steps.size(); n++ ) {
		ResidueOP new_rot = rot->clone();
		new_rot->set_chi( 1, center + chi_steps[n] * width );
		rotamers.push_back( new_rot );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_rna_chi_rotamers( conformation::ResidueOP const & rot,
	utility::vector1< conformation::ResidueOP > & rotamers,
	pack::task::ExtraRotSample const & level,
	utility::vector1< chemical::rna::GaussianParameter > const & gaussian_parameter_set )
{
	using namespace conformation;
	using namespace chemical::rna;
	using namespace pack::task;

	utility::vector1< Real > chi_steps;
	chi_steps.push_back( 0.0 );

	//This copies code, unfortunately, from RotamerSet_.cc
	switch ( level ) {
	case NO_EXTRA_CHI_SAMPLES : // 0
		//  return;
		break;
	case EX_ONE_STDDEV : // 1
		chi_steps.push_back(1);
		chi_steps.push_back(-1);
		break;
	case EX_ONE_HALF_STEP_STDDEV : // 2
		chi_steps.push_back(0.5);
		chi_steps.push_back(-0.5);
		break;
	case EX_TWO_FULL_STEP_STDDEVS : // 3
		chi_steps.push_back(1);
		chi_steps.push_back(2);
		chi_steps.push_back(-1);
		chi_steps.push_back(-2);
		break;
	case EX_TWO_HALF_STEP_STDDEVS : // 4
		chi_steps.push_back(0.5);
		chi_steps.push_back(1);
		chi_steps.push_back(-0.5);
		chi_steps.push_back(-1);
		break;
	case EX_FOUR_HALF_STEP_STDDEVS : // 5
		chi_steps.push_back(0.5);
		chi_steps.push_back(1);
		chi_steps.push_back(1.5);
		chi_steps.push_back(2.0);
		chi_steps.push_back(-0.5);
		chi_steps.push_back(-1);
		chi_steps.push_back(-1.5);
		chi_steps.push_back(-2);
		break;
	case EX_THREE_THIRD_STEP_STDDEVS : // 6
		chi_steps.push_back(0.33);
		chi_steps.push_back(0.67);
		chi_steps.push_back(1);
		chi_steps.push_back(-0.33);
		chi_steps.push_back(-0.67);
		chi_steps.push_back(-1);
		break;
	case EX_SIX_QUARTER_STEP_STDDEVS : // 7
		chi_steps.push_back(0.25);
		chi_steps.push_back(0.5);
		chi_steps.push_back(0.75);
		chi_steps.push_back(1);
		chi_steps.push_back(1.25);
		chi_steps.push_back(1.5);
		chi_steps.push_back(-0.25);
		chi_steps.push_back(-0.5);
		chi_steps.push_back(-0.75);
		chi_steps.push_back(-1);
		chi_steps.push_back(-1.25);
		chi_steps.push_back(-1.5);
		break;
	case ExtraRotSampleCardinality :
	default :
		std::cerr << "Error in RotamerSet_::set_extrachi_samples, invalid ExtraChiSample type" << std::endl;
		utility_exit();
		break;
	}

	GaussianParameter const & gaussian_parameter( gaussian_parameter_set[1] );
	fill_chi_rotamers_with_center_and_stddev( rot, rotamers,
		chi_steps,
		gaussian_parameter.center, gaussian_parameter.width );

	//Add in syn configurations for purines.
	if ( rot->aa() == chemical::na_rad || rot->aa() == chemical::na_rgu ) {
		GaussianParameter const & gaussian_parameter2( gaussian_parameter_set[2] );
		fill_chi_rotamers_with_center_and_stddev( rot, rotamers,
			chi_steps,
			gaussian_parameter2.center, gaussian_parameter2.width );
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
build_rna_chi_rotamers(
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	pack::task::ExtraRotSample const & level,
	bool const sample_rna_chi, bool const & include_current,
	utility::vector1< conformation::ResidueOP > & rotamers
)
{
	using numeric::principal_angle_degrees;

	using namespace pose;
	using namespace conformation;
	using namespace id;
	using namespace pack::task;

	// Real const max_rotation( 5.0 ); // degrees
	// bool const tether_bond_distance( true );
	// bool const tether_bond_angle( false );

	Residue const & existing_residue( pose.residue( resid ) );
	debug_assert( existing_residue.is_RNA() && concrete_residue->is_RNA() );

	// a residue of the correct type aligned to the existing backbone
	ResidueOP rot = ResidueFactory::create_residue( *concrete_residue, existing_residue, pose.conformation(), false /*true*/ /*preserve c_beta*/ );

	//To define base, need a torsion that depends on sugar pucker...
	static chemical::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;
	static Real const delta_cutoff = rna_fitted_torsion_info.delta_cutoff();

	//Different chi angles depending on sugar pucker.
	Real const delta( principal_angle_degrees( existing_residue.mainchain_torsion( chemical::rna::DELTA ) ) );

	// Note: do not necessarily include current chi, unless specified by user!
	if ( include_current && sample_rna_chi ) {
		//Can't make a direct copy, but use ideal bonds/angles for base, and copy chi
		ResidueOP new_rot = rot->clone();
		rotamers.push_back( new_rot );
	}

	if ( !sample_rna_chi && level == 0 ) return;

	//Different chi angles depending on sugar pucker.
	if ( delta <= delta_cutoff ) {
		add_rna_chi_rotamers( rot, rotamers, level, rna_fitted_torsion_info.gaussian_parameter_set_chi_north() );
	} else {
		add_rna_chi_rotamers( rot, rotamers, level, rna_fitted_torsion_info.gaussian_parameter_set_chi_south() );
	}

}


/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Expand proton chi? For RNA, this covers the 2'-OH torsion -- critical!.
//  This is a pretty generic operation  Should be its own method?
// Might be worth re-integrating with main-stream proton chi code --
// would just need to be smart about 2'-OH virtualization below.
//
void
build_proton_chi_rotamers(
	Size const /* resid */,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	pack::task::ResidueLevelTask const & residue_task,
	utility::vector1< conformation::ResidueOP > & rotamers )
{
	using namespace conformation;
	Size const n_proton_chi = concrete_residue->n_proton_chi();
	bool const & include_virtual_side_chain = residue_task.include_virtual_side_chain(); /* this means 2'-OH for RNA*/

	if ( n_proton_chi > 0 ) {
		// This seems a little silly -- suck out the chi's, then put them back into rotamers.
		utility::vector1< pack::dunbrack::ChiSetOP > proton_chi_chisets;
		proton_chi_chisets.push_back( dunbrack::ChiSetOP( new pack::dunbrack::ChiSet( concrete_residue->nchi() ) ) );

		for ( Size ii = 1; ii <= n_proton_chi; ++ii ) {
			pack::dunbrack::expand_proton_chi( residue_task.extrachi_sample_level( true /*nneighb test*/,
				concrete_residue->proton_chi_2_chi( ii ), *concrete_residue ),
				concrete_residue, ii, proton_chi_chisets);
		}

		Size const number_of_starting_rotamers = rotamers.size();
		for ( Size n = 1 ; n <= number_of_starting_rotamers; ++n ) {

			Size n_min( 0 ); Real rna_torsion_min( 0.0 );
			for ( Size ii = 1; ii <= proton_chi_chisets.size(); ++ii ) {

				rotamers.push_back( rotamers[n]->clone() );
				Size const new_rotamer_number = rotamers.size();

				for ( Size jj = 1; jj <= n_proton_chi; ++jj ) {
					Size const jj_protchi( concrete_residue->proton_chi_2_chi( jj ) );
					ResidueOP rot = rotamers[ new_rotamer_number ];
					rot->set_chi( jj_protchi,
						proton_chi_chisets[ ii ]->chi[ jj_protchi ] );

					if ( include_virtual_side_chain ) {
						runtime_assert( !rot->has_variant_type( chemical::VIRTUAL_O2PRIME_HYDROGEN ) );
						static core::scoring::rna::RNA_EnergyMethodOptions const rna_options;
						static core::scoring::rna::RNA_TorsionPotential rna_torsion_potential( rna_options ); // I think this is OK.
						Real const rna_torsion = rna_torsion_potential.eval_intrares_energy( *rot, pose );
						if ( n_min == 0 || rna_torsion < rna_torsion_min ) {
							n_min = new_rotamer_number; rna_torsion_min = rna_torsion;
						}
					}
				}

			}
			if ( include_virtual_side_chain ) {
				ResidueOP rot = core::pose::add_variant_type_to_residue( *rotamers[ n_min ], chemical::VIRTUAL_O2PRIME_HYDROGEN, pose);
				rotamers.push_back( rot );
			}
		}
	}

} //proton chi


/////////////////////////////////////////////////////////////////////////////
void
build_five_prime_phosphate_rotamers( utility::vector1< conformation::ResidueOP > & rotamers,
	pose::Pose const & pose, bool const allow_phosphate_virtualization ){
	using namespace core::conformation;

	utility::vector1< Real > const full_torsions  = chemical::rna::get_full_torsions();

	Size const number_of_starting_rotamers = rotamers.size();
	for ( Size n = 1 ; n <= number_of_starting_rotamers; ++n ) {

		runtime_assert( rotamers[n]->has_variant_type( chemical::FIVE_PRIME_PACKABLE_PHOSPHATE ) );
		// this nested iteration is totally ghetto -- would be better to unify with fang's protocols/rotamer_sample/rna/RNA_SuiteRotamer code.
		for ( Size i = 1; i <= full_torsions.size(); i++ ) {
			for ( Size j = 1; j <= full_torsions.size(); j++ ) {
				for ( Size k = 1; k <= full_torsions.size(); k++ ) {
					ResidueOP  new_rsd = rotamers[n]->clone();
					Size const chi_number_i = rotamers[n]->type().RNA_type().chi_number_pseudoalpha();
					Size const chi_number_j = rotamers[n]->type().RNA_type().chi_number_pseudobeta();
					Size const chi_number_k = rotamers[n]->type().RNA_type().chi_number_pseudogamma();
					runtime_assert( chi_number_i > 0 );
					runtime_assert( chi_number_j > 0 );
					runtime_assert( chi_number_k > 0 );
					new_rsd->set_chi( chi_number_i, full_torsions[ i ] );
					new_rsd->set_chi( chi_number_j, full_torsions[ j ] );
					new_rsd->set_chi( chi_number_k, full_torsions[ k ] );
					rotamers.push_back( new_rsd );
				} //k
			} //j
		} //i

		if ( allow_phosphate_virtualization ) {
			//also add a rotamer that removes this variant!
			ResidueOP new_rsd = remove_variant_type_from_residue( *rotamers[n], chemical::FIVE_PRIME_PACKABLE_PHOSPHATE, pose );
			rotamers.push_back( new_rsd );
		}

	}

}


/////////////////////////////////////////////////////////////////////////////
void
build_three_prime_phosphate_rotamers( utility::vector1< conformation::ResidueOP > & rotamers,
	pose::Pose const & pose,
	bool const allow_phosphate_virtualization ){

	using namespace core::conformation;
	using namespace core::chemical::rna;
	utility::vector1< Real > const full_torsions  = get_full_torsions( 10.0 /*bin size*/ );

	Size const number_of_starting_rotamers = rotamers.size();
	for ( Size n = 1 ; n <= number_of_starting_rotamers; ++n ) {
		utility::vector1< Real > const epsilon_torsions  = get_epsilon_torsions( rotamers[n]->mainchain_torsion( DELTA ),
			true /*extra epsilon*/, 10.0 /*bin size*/ );
		runtime_assert( rotamers[n]->has_variant_type( chemical::THREE_PRIME_PACKABLE_PHOSPHATE ) );
		// this nested iteration is totally ghetto -- would be better to unify with fang's protocols/rotamer_sample/rna/RNA_SuiteRotamer code.
		for ( Size i = 1; i <= epsilon_torsions.size(); i++ ) {
			for ( Size j = 1; j <= full_torsions.size(); j++ ) {
				ResidueOP  new_rsd = rotamers[n]->clone();
				Size const chi_number_i = rotamers[n]->type().RNA_type().chi_number_pseudoepsilon();
				Size const chi_number_j = rotamers[n]->type().RNA_type().chi_number_pseudozeta();
				runtime_assert( chi_number_i > 0 );
				runtime_assert( chi_number_j > 0 );
				new_rsd->set_chi( chi_number_i, epsilon_torsions[ i ] );
				new_rsd->set_chi( chi_number_j, full_torsions[ j ] );
				rotamers.push_back( new_rsd );
			} //j
		} //i

		if ( allow_phosphate_virtualization ) {
			//also add a rotamer that removes this variant!
			ResidueOP new_rsd = remove_variant_type_from_residue( *rotamers[n], chemical::THREE_PRIME_PACKABLE_PHOSPHATE, pose );
			rotamers.push_back( new_rsd );
		}

	}

}

/////////////////////////////////////////////////////////////////////////////
void
build_rna_rotamers(
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	pack::task::PackerTask const & task,
	utility::vector1< conformation::ResidueOP > & new_rotamers,
	Size & id_for_current_rotamer
)
{
	using namespace conformation;

	pack::task::ExtraRotSample const level
		( task.residue_task( resid ).extrachi_sample_level( true /*buried*/, 1 /*chi*/, *concrete_residue ) );

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	// Basic set up...
	build_rna_chi_rotamers( resid, pose, concrete_residue, level,
		task.residue_task( resid ).rna_task().sample_rna_chi(),
		task.include_current( resid ),
		new_rotamers );

	//////////////////////////////////////////////////////////////////////////
	// include current... this doesn't use ideal bond lengths/angles.
	Residue const & existing_residue( pose.residue( resid ));
	if ( task.include_current( resid ) && existing_residue.aa() == concrete_residue->aa() ) {
		ResidueOP rot =  remove_variant_type_from_residue( existing_residue, chemical::VIRTUAL_O2PRIME_HYDROGEN, pose );
		new_rotamers.push_back( rot );
		id_for_current_rotamer = new_rotamers.size();
	}

	///////////////////////////////////////////////////////////
	core::pack::task::rna::RNA_ResidueLevelTask const & rna_task = task.residue_task( resid ).rna_task();
	if ( rna_task.sample_five_prime_phosphate() )  {
		TR << "Going to build rotamers for  5' phosphate at " << resid << "  before: " << new_rotamers.size() << std::endl;
		build_five_prime_phosphate_rotamers( new_rotamers, pose, rna_task.allow_phosphate_virtualization() );
		TR << "  after building rotamers: " << new_rotamers.size() << std::endl;
	}
	if ( rna_task.sample_three_prime_phosphate() ) {
		TR << "Going to build rotamers for  3' phosphate at " << resid << "  before: " << new_rotamers.size() << std::endl;
		build_three_prime_phosphate_rotamers( new_rotamers, pose, rna_task.allow_phosphate_virtualization() );
		TR << "  after building rotamers: " << new_rotamers.size() << std::endl;
		TR << " starting pseudo-epsilon: " << new_rotamers[1]->chi( new_rotamers[1]->type().RNA_type().chi_number_pseudoepsilon() ) << std::endl;
		TR << " starting    pseudo-zeta: " << new_rotamers[1]->chi( new_rotamers[1]->type().RNA_type().chi_number_pseudozeta() ) << std::endl;
	}

	///////////////////////////////////////////////////////////
	if ( task.residue_task( resid ).sample_proton_chi() ) {
		build_proton_chi_rotamers( resid, pose, concrete_residue,
			task.residue_task( resid ),
			new_rotamers );
	}

}


} // rotamer_set
} // pack
} // core
