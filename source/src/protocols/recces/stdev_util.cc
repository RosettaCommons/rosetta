// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/stdev_util.cc
/// @brief encodes 'reasonable' assumptions about sampling stdevs.
/// @detailed stdevs change with temperature.
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/recces/stdev_util.hh>
#include <protocols/recces/sampler/MC_Any.hh>
#include <protocols/recces/sampler/MC_Comb.hh>
#include <protocols/recces/sampler/MC_Loop.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>
#include <protocols/recces/sampler/rna/MC_RNA_OneJump.hh>
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>
#include <protocols/recces/options/RECCES_Options.hh>
#include <protocols/recces/params/RECCES_Parameters.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.recces.stdev_util" );

using namespace core;
using namespace utility;

namespace protocols {
namespace recces {

///////////////////////////////////////////////////////////////////////////////////
// TODO: Unify with set_gaussian_stdevs_thermal_sampler (see note below)
void
set_gaussian_stdevs_legacy_turner( protocols::recces::sampler::MC_CombOP sampler,
	core::Real const & temperature,
	core::pose::Pose const & pose,
	protocols::recces::params::RECCES_Parameters const & params )
{
	using namespace core::id;
	using namespace core::chemical::rna;
	using namespace core::pose::rna;
	using namespace protocols::recces::sampler;

	Size const n_rsd( pose.total_residue() );
	Real const bp_stdev(       gaussian_stdev( n_rsd, temperature, true ) );
	Real const dangling_stdev( gaussian_stdev( n_rsd, temperature, false ) );
	Real stdev( 0.0 );
	utility::vector1< Real > const & bp_res( params.bp_res() );
	utility::vector1< Real > const & dangling_res( params.dangling_res() );

	runtime_assert( n_rsd == bp_res.size() + dangling_res.size() );

	////////////////////////////////////////////////
	// update gaussian stdev for all chi's.
	////////////////////////////////////////////////
	for ( Size i = 1; i <= n_rsd; ++ i ) {
		if ( bp_res.has_value(i) )            stdev = bp_stdev;
		else  stdev = dangling_stdev;
		MC_SamplerOP torsion_sampler = sampler->find( TorsionID( i, TorsionType::CHI, 1 ) );
		runtime_assert( torsion_sampler != 0 ); // these all move in RECCES
		std::dynamic_pointer_cast< MC_OneTorsion >(torsion_sampler)->set_gaussian_stdev( stdev );
	}

	//////////////////////////////////////////////////////////////////
	// no need to update sugars -- they do not have gaussian ranges.
	//////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////
	// update gaussian stdev for all suites
	////////////////////////////////////////////////
	for ( Size i = 1; i < n_rsd; ++ i ) { // watch out: may need to get last residue if we cyclize
		if ( pose.fold_tree().is_cutpoint( i ) ) continue; // watch out: later generalize to cutpoint_closed
		if ( bp_res.has_value(i) && bp_res.has_value( i + 1) )  stdev = bp_stdev;
		else stdev = dangling_stdev;
		vector1< TorsionID > suite_torsion_ids = get_suite_torsion_ids( i );
		for ( auto bb_torsion_id : suite_torsion_ids ) {
			MC_SamplerOP torsion_sampler = sampler->find( bb_torsion_id );
			runtime_assert( torsion_sampler != 0 ); // these all move in RECCES
			std::dynamic_pointer_cast< MC_OneTorsion >(torsion_sampler)->set_gaussian_stdev( stdev );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Simple heuristic for gaussian stdev
core::Real gaussian_stdev( core::Real const n_rsd, core::Real const temp, bool const is_bp ) {
	// Negative temp is infinite
	if ( temp < 0 ) return -1;
	if ( is_bp ) return 5 * temp / n_rsd;
	return 6 * std::pow( temp / n_rsd, 0.75 );
}

///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// replace this eventually using sampler::find()
void set_gaussian_stdevs_thermal_sampler(
	protocols::recces::sampler::MC_CombOP internal_bb_sampler,
	protocols::recces::sampler::MC_CombOP chi_sampler,
	protocols::recces::sampler::rna::MC_RNA_MultiSuiteOP standard_bb_sampler,
	Real const & temp,
	Size const & total_rsd,
	Size const & sampled_rsd,
	utility::vector1<bool> const & is_free )
{
	using namespace protocols::recces::sampler;
	using namespace protocols::recces::sampler::rna;
	Real internal_bb_stdev( 0.1 * pow( temp, 0.25 ) + 0.1);
	Real free_chi_stdev( 55 * pow( temp, 0.5 ) + 50 );
	Real chi_stdev( 5 * pow( temp , 0.5) + 15 );
	Real standard_bb_stdev( 8 * pow( temp, 0.5 ) / (2 * total_rsd + sampled_rsd) );
	if ( temp < 0 ) {
		internal_bb_stdev = 0.5 ;
		free_chi_stdev = -1 ;
		chi_stdev = -1 ;
		standard_bb_stdev = -1 ;
	}
	if ( internal_bb_sampler != 0 ) {
		for ( Size i = 1; i <= internal_bb_sampler->num_rotamers(); ++i ) {
			runtime_assert(  (*internal_bb_sampler)[i]->type() == toolbox::MC_RNA_KIC );
			std::dynamic_pointer_cast< MC_RNA_KIC_Sampler >( (*internal_bb_sampler)[i] )->set_gaussian_stdev( internal_bb_stdev );
		}
	}
	if ( chi_sampler != 0 ) {
		for ( Size i = 1; i <= chi_sampler->num_rotamers(); ++i ) {
			runtime_assert(  (*chi_sampler)[i]->type() == toolbox::MC_ONE_TORSION );
			if ( is_free[i] ) {
				std::dynamic_pointer_cast< MC_OneTorsion >( (*chi_sampler)[i] )->set_gaussian_stdev( free_chi_stdev );
			} else {
				std::dynamic_pointer_cast< MC_OneTorsion >( (*chi_sampler)[i] )->set_gaussian_stdev( chi_stdev );
			}
		}
	}
	if ( standard_bb_sampler != 0 ) standard_bb_sampler->set_gaussian_stdev( standard_bb_stdev );
}

//////////////////////////////////////////////////////////////////////////////
// TODO: replace this eventually using sampler::find()
// TODO: unify with set_gaussian_stdevs_recces_turner, based on whether a residue
//        is specified as HELIX or DANGLE or FREE or SAMPLE
//        [all inside pose full_model or in Options]
void
set_gaussian_stdevs_thermal_sampler(
	protocols::recces::sampler::MC_SamplerOP sampler,
	core::Real const & temperature,
	core::pose::Pose const & pose,
	protocols::recces::options::RECCES_Options const & options )
{

	using namespace protocols::toolbox;
	using namespace protocols::recces::sampler;
	using namespace protocols::recces::sampler::rna;

	// assume that this is thermal_sampler type. Later generalize, based on sampler.find( TorsionID ) and sampler.find( TorsionIDs ).
	runtime_assert( sampler->type() == MC_LOOP );
	MC_Loop & loop_sampler( *std::dynamic_pointer_cast< MC_Loop >( sampler ) );

	MC_CombOP chi_sampler, internal_bb_sampler;
	MC_RNA_MultiSuiteOP standard_bb_sampler;
	if ( loop_sampler.num_rotamers() >= 1 ) {
		runtime_assert( loop_sampler[ 1 ]->type() == MC_ANY );
		chi_sampler = MC_CombOP( std::dynamic_pointer_cast< MC_Comb >( loop_sampler[ 1 ] ) );
	}
	if ( loop_sampler.num_rotamers() >= 2 ) {
		runtime_assert( loop_sampler[ 2 ]->type() == MC_ANY );
		internal_bb_sampler = MC_CombOP( std::dynamic_pointer_cast< MC_Comb >( loop_sampler[ 2 ] ) );
	}
	if ( loop_sampler.num_rotamers() >= 10 ) {
		runtime_assert( loop_sampler[ 10 ]->type() == MC_RNA_MULTI_SUITE );
		standard_bb_sampler = MC_RNA_MultiSuiteOP( std::dynamic_pointer_cast< MC_RNA_MultiSuite >( loop_sampler[ 10 ] ) );
	}

	utility::vector1< bool > is_free;
	for ( auto sample_res : options.sample_residues() ) is_free.push_back( options.free_residues().has_value( sample_res ) );

	set_gaussian_stdevs_thermal_sampler( internal_bb_sampler, chi_sampler, standard_bb_sampler, temperature, pose.size(), options.sample_residues().size(), is_free );

	// TODO: at higher temperatures, use larger translation/rotation magnitude for *jump_sampler* if present (as 5th of the 10 cycles in each loop, turned on with -sample_jump)

}


} //recces
} //protocols
