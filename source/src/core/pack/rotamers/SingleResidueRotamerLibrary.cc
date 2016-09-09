// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleResidueRotamerLibrary.cc
/// @brief  SingleResidueRotamerLibrary class
/// @author Andrew Leaver-Fay
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit Headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/BumpSelector.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/dunbrack/ChiSet.hh>

#include <utility/graph/Graph.hh>

#include <utility/LexicographicalIterator.hh>
#include <numeric/random/reservoir_sample.hh>

#include <basic/Tracer.hh>

namespace core {
namespace pack {
namespace rotamers {

static THREAD_LOCAL basic::Tracer TR( "core.pack.rotamers.SingleResidueRotamerLibrary" );

SingleResidueRotamerLibrary::~SingleResidueRotamerLibrary()
{}

/* OL: I thought copying the whole chi_set_vector is unnecessary and made a new version of this function */
/*
void
expand_proton_chi_oldversion(
pack::task::ExtraRotSample ex_samp_level,
chemical::ResidueTypeCOP concrete_residue,
Size proton_chi,
utility::vector1< dunbrack::ChiSetOP > & chi_set_vector
)
{
using namespace pack::task;
using namespace pack::dunbrack;
// count the number of extra hydroxyl samples -- n
// copy the chi_set_vector n times into a temporary
// set the hydroxyl_chi value for these n samples in the temporary
// assign the temporary to the input chi_set_vector

utility::vector1< Real > const & samples = concrete_residue->proton_chi_samples( proton_chi );
// i.e., -60, 60, 180

// i.e., 10, 20  to model -40 -50 -60 -70 -80 etc.
utility::vector1< Real > const & extra_samples = concrete_residue->proton_chi_extra_samples( proton_chi );

Size chi_id = concrete_residue->proton_chi_2_chi( proton_chi );

bool const include_extra( ex_samp_level != NO_EXTRA_CHI_SAMPLES );

Size nsamples = samples.size() * ( 1 + ( include_extra ? extra_samples.size() * 2 : 0 ) );
utility::vector1< ChiSetOP > newchi_vect( nsamples * chi_set_vector.size() );

// copy old chi_set_vector nsample times
for ( Size ii = 1; ii <= nsamples; ++ii ) {
Size offset = (ii-1) * chi_set_vector.size();
for ( Size jj = 1; jj <= chi_set_vector.size(); ++jj ) {
newchi_vect[ jj + offset ] = ChiSetOP( new pack::dunbrack::ChiSet( *(chi_set_vector[ jj ]) ) );
}
}

// add extra chi samples
Size count( 1 );
for ( Size ii = 1; ii <= samples.size(); ++ii ) {
Real ii_sample = samples[ ii ]; //chi-angle
for ( Size jj = 1; jj <= chi_set_vector.size(); ++jj ) {
newchi_vect[ count ]->chi[ chi_id ] = ii_sample;
++count;
if ( ! include_extra ) continue;

for ( Size kk = 1; kk <= extra_samples.size(); ++kk ) {
newchi_vect[ count ]->chi[ chi_id ] = ii_sample + extra_samples[ kk ];
++count;
newchi_vect[ count ]->chi[ chi_id ] = ii_sample - extra_samples[ kk ];
++count;
}
}
}

debug_assert( count - 1 == nsamples * chi_set_vector.size() );
chi_set_vector = newchi_vect;
}

/// olli -- I think this is slightly simpler to read than the original version
/// apl -- needs to find a new residence
void
SingleResidueRotamerLibrary::expand_proton_chi(
pack::task::ExtraRotSample ex_samp_level,
chemical::ResidueTypeCOP concrete_residue,
Size proton_chi,
utility::vector1< dunbrack::ChiSetOP > & chi_set_vector
)
{
using namespace pack::task;
using namespace pack::dunbrack;
// count the number of extra hydroxyl samples -- n
// copy the chi_set_vector n times into a temporary
// set the hydroxyl_chi value for these n samples in the temporary
// assign the temporary to the input chi_set_vector

utility::vector1< Real > const & samples = concrete_residue->proton_chi_samples( proton_chi );

// i.e., -60, 60, 180
debug_assert( samples.size() > 0 ); // or less harsh and just a return ?

// i.e., 10, 20  to model -40 -50 -60 -70 -80 etc.
utility::vector1< Real > const & extra_samples = concrete_residue->proton_chi_extra_samples( proton_chi );

Size chi_id = concrete_residue->proton_chi_2_chi( proton_chi );

bool const include_extra( ex_samp_level != NO_EXTRA_CHI_SAMPLES );

Size nsamples = samples.size() * ( 1 + ( include_extra ? extra_samples.size() * 2 : 0 ) );
chi_set_vector.reserve( nsamples * chi_set_vector.size() ); //preallocate the necessary memory

// add extra chi samples
ChiSetOP new_chi_vec; ChiSetOP base_chi_vec;
Size nr_of_old_elem = chi_set_vector.size();
for ( Size jj = 1; jj <= nr_of_old_elem; ++jj ) {
for ( Size ii = 1; ii <= samples.size(); ++ii ) {
Real ii_sample = samples[ ii ]; //chi-angle
if ( ii == 1 ) {  //change first chi in place:
base_chi_vec = new_chi_vec = chi_set_vector[ jj ];
} else { // make copies for all others
new_chi_vec = ChiSetOP( new pack::dunbrack::ChiSet( *base_chi_vec ) );
chi_set_vector.push_back( new_chi_vec );
}
new_chi_vec->chi[ chi_id ] = ii_sample;

if ( ! include_extra )  continue;

for ( Size kk = 1; kk <= extra_samples.size(); ++kk ) {
chi_set_vector.push_back( new_chi_vec = ChiSetOP( new pack::dunbrack::ChiSet( *base_chi_vec  ) ) );
new_chi_vec->chi[ chi_id ] = ii_sample  + extra_samples[ kk ];
chi_set_vector.push_back( new_chi_vec = ChiSetOP( new pack::dunbrack::ChiSet( *base_chi_vec  ) ) );
new_chi_vec->chi[ chi_id ] = ii_sample  - extra_samples[ kk ];
} // for extra_samples
} // for sample.size()
} // for jj (chi_set_vector)

debug_assert( chi_set_vector.size()  == nsamples * nr_of_old_elem );
}

utility::vector1< dunbrack::ChiSetOP >
SingleResidueRotamerLibrary::expand_proton_chis_old(
chemical::ResidueTypeCOP concrete_residue,
pack::task::ResidueLevelTask const & rlt,
bool buried
) {
utility::vector1< pack::dunbrack::ChiSetOP > proton_chi_chisets;
proton_chi_chisets.push_back(
dunbrack::ChiSetOP( new pack::dunbrack::ChiSet( concrete_residue->nchi() ) ) );
for ( Size ii = 1; ii <= concrete_residue->n_proton_chi(); ++ii ) {
SingleResidueRotamerLibrary::expand_proton_chi(
rlt.extrachi_sample_level(
buried,
concrete_residue->proton_chi_2_chi( ii ),
*concrete_residue ),
concrete_residue,
ii, proton_chi_chisets);
}
return proton_chi_chisets;
}
*/

utility::vector1< utility::vector1< core::Real > >
SingleResidueRotamerLibrary::compute_proton_chi_samplings(
	chemical::ResidueType const & concrete_residue,
	pack::task::ResidueLevelTask const & rlt,
	bool buried
) const {

	core::Size n_proton_chi( concrete_residue.n_proton_chi() );

	utility::vector1< utility::vector1< core::Real > > sampling( n_proton_chi ); // A proton chi-indexed vector of vectors of samples to use.

	for ( Size proton_chi = 1; proton_chi <= n_proton_chi; ++proton_chi ) {

		utility::vector1< core::Real > & chi_samples( sampling[ proton_chi ] );

		// i.e., -60, 60, 180
		utility::vector1< Real > const & base_samples = concrete_residue.proton_chi_samples( proton_chi );

		bool include_extra( rlt.extrachi_sample_level(
			buried,
			concrete_residue.proton_chi_2_chi( proton_chi ),
			concrete_residue )
			!= pack::task::NO_EXTRA_CHI_SAMPLES );

		// i.e., 10, 20  to model -40 -50 -60 -70 -80 etc.
		utility::vector1< Real > const & extra_samples = concrete_residue.proton_chi_extra_samples( proton_chi );

		for ( Size jj = 1; jj <= base_samples.size(); ++jj ) {
			chi_samples.push_back( base_samples[jj] );

			if ( ! include_extra ) { continue; }

			for ( Size kk = 1; kk <= extra_samples.size(); ++kk ) {
				chi_samples.push_back( base_samples[jj] + extra_samples[ kk ] );
				chi_samples.push_back( base_samples[jj] - extra_samples[ kk ] );
			}
		}
	}

	return sampling;
}

utility::vector1< dunbrack::ChiSetOP >
SingleResidueRotamerLibrary::expand_proton_chis(
	utility::vector1< utility::vector1< core::Real > > const & sampling, // A proton chi-indexed vector of vectors of samples to use.
	chemical::ResidueType const & concrete_residue,
	core::Size max_rotamers // The maximum number of proton-chi-expanded rotamers to produce
) const {
	debug_assert( sampling.size() == concrete_residue.n_proton_chi() );

	// For pathological cases, even building a lightweight representation of all the combinations kills you in memory.
	// Avoid building as much as possible.

	utility::vector1< core::Size > dim_sizes;
	for ( core::Size ii(1); ii <= sampling.size(); ++ii ) {
		dim_sizes.push_back( sampling[ii].size() );
	}

	utility::LexicographicalIterator lexi_it( dim_sizes );

	TR.Trace << "Untrimmed factor for proton chi expansions with " << concrete_residue.name() << " is " << lexi_it.num_states_total() << std::endl;

	if ( lexi_it.num_states_total() > max_rotamers ) {
		TR.Warning << "Subsampling proton chi expansion for " << concrete_residue.name() << ", as the estimated factor of "
			<< lexi_it.num_states_total() << " exceeds " << max_rotamers << " proton rotamers per base residue!" << std::endl;
	}

	// The ReservoirSampler will accumulate at most the max_rotamers states in memory.
	// Also, it won't try to (re)store each and every state, which minimizes copying
	// But the order will not necessarily be preserved for items which exceed max_rotamers
	numeric::random::ReservoirSampler< utility::vector1< core::Size > > subsampler( max_rotamers );

	for ( ; !lexi_it.at_end(); ++lexi_it ) {
		subsampler.add_value( *lexi_it );
	}

	utility::vector1< utility::vector1< core::Size > > const & indicies( subsampler.values() );

	using namespace core::pack::dunbrack;
	utility::vector1< ChiSetOP > chi_set_vector;
	chi_set_vector.reserve( indicies.size() );

	for ( core::Size ii(1); ii <= indicies.size(); ++ii ) {
		utility::vector1< core::Size > const & indexvect( indicies[ii] );
		debug_assert( indexvect.size() == concrete_residue.n_proton_chi() );
		ChiSetOP new_chi_vec( new ChiSet( concrete_residue.nchi() ) );
		for ( core::Size proton_chi(1); proton_chi <= indexvect.size(); ++proton_chi ) {
			utility::vector1< core::Real > const & chi_samples( sampling[proton_chi] );
			Size chi_id = concrete_residue.proton_chi_2_chi( proton_chi );
			new_chi_vec->chi[ chi_id ] = chi_samples[ indexvect[proton_chi] ];
		}
		chi_set_vector.push_back( new_chi_vec );
	}

	return chi_set_vector;
}


void
SingleResidueRotamerLibrary::bump_filter(
	RotamerVector & rotamers,
	core::Size resid,
	scoring::ScoreFunction const & scorefxn,
	pose::Pose const & pose,
	task::PackerTask const & task,
	utility::graph::GraphCOP packer_neighbor_graph
) const
{
	using namespace core::pack::rotamer_set; // Temporary for BumpSelectors

	if ( rotamers.size() == 0 || ! task.bump_check() ) { return; } //Nothing to do
	BumpSelector bump_selector( task.max_rotbump_energy() );

	RotamerVector passing;

	TR.Debug << "Filtering rotamers on bump energy for residue " << resid << " " << rotamers[1]->name() << std::endl;
	for ( Size ii = 1; ii <= rotamers.size(); ++ii ) {
		core::conformation::ResidueOP rot = rotamers[ ii ];
		TR.Debug << " Suggested rotamer " << ii << ": ";
		for ( Size jj = 1; jj <= rot->nchi(); ++jj ) {
			TR.Debug << rot->chi()[jj] << ' ';
		}

		core::PackerEnergy bumpenergy = bump_check( rot, resid, scorefxn, pose, task, packer_neighbor_graph );
		TR.Debug << " Bump energy: " << bumpenergy;
		BumpSelectorDecision decision =  bump_selector.iterate_bump_selector( bumpenergy );
		switch ( decision ) {
		case KEEP_ROTAMER :
			TR.Debug << " ... added" << std::endl;
			passing.push_back( rot );
			break;
		case DELETE_PREVIOUS_ROTAMER :
			TR.Debug << " ... replace previous" << std::endl;
			if ( passing.size() == 0 ) {
				utility_exit_with_message("Internal consistency error: cannot replace non-existant previous residue.");
			}
			passing[ passing.size() ] = rot;
			break;
		case DELETE_ROTAMER : // Do nothing.
			TR.Debug << " ... deleted" << std::endl;
			break;
		}
	}
	// Swap the in/out value with the contents of the temporary vector to pass back the value
	TR.Debug << " N rotamers before: " << rotamers.size() << " after: " << passing.size() << std::endl;
	rotamers.swap(passing);
}

/// @details Bump check does not include long range energies,
/// though, maybe this should change.
core::PackerEnergy
SingleResidueRotamerLibrary::bump_check(
	core::conformation::ResidueCOP rotamer,
	core::Size resid,
	scoring::ScoreFunction const & sf,
	pose::Pose const & pose,
	task::PackerTask const & task,
	utility::graph::GraphCOP packer_neighbor_graph
) const
{
	using namespace scoring;
	using namespace conformation;

	EnergyMap emap;

	for ( utility::graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( resid )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( resid )->const_edge_list_end();
			ir != ire; ++ir ) {
		int const neighbor_id( (*ir)->get_other_ind( resid ) );
		Residue const & neighbor( pose.residue( neighbor_id ) );

		if ( ! task.pack_residue( neighbor_id ) ) {
			sf.bump_check_full( *rotamer, neighbor, pose, emap);
		} else {
			sf.bump_check_backbone( *rotamer, neighbor, pose, emap);
		}
	}
	return static_cast< core::PackerEnergy > (sf.weights().dot( emap ));
}

core::Size
SingleResidueRotamerLibrary::current_rotamer(
	RotamerVector & rotamers,
	core::Size resid,
	task::PackerTask const & task,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue
) const
{
	if ( task.include_current( resid ) && existing_residue.name() == concrete_residue->name() ) {
		TR.Debug << " Including current rotamer: " <<
			existing_residue.name() << ' ' << existing_residue.seqpos() << std::endl;
		conformation::ResidueOP rot = existing_residue.create_rotamer();
		rotamers.push_back( rot );
		return rotamers.size();
	} else {
		return 0;
	}
}

void
SingleResidueRotamerLibrary::emergency_rotamer(
	RotamerVector & rotamers,
	core::Size /*resid*/,
	pose::Pose const & pose,
	task::PackerTask const & /*task*/,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue
) const
{
	if ( rotamers.size() == 0 && concrete_residue->nchi() == 0 ) {
		TR.Debug << " Using emergency rotamer" << std::endl;
		conformation::ResidueOP rot = conformation::ResidueFactory::create_residue( *concrete_residue, existing_residue, pose.conformation() );
		rotamers.push_back( rot );
	}
}

RotamerVector
SingleResidueRotamerLibrary::virtual_sidechain(
	RotamerVector const & rotamers,
	core::Size resid,
	pose::Pose const & pose,
	task::PackerTask const & task,
	chemical::ResidueTypeCOP /*concrete_residue*/,
	conformation::Residue const & existing_residue
) const {
	RotamerVector retval;
	if ( task.residue_task( resid ).include_virtual_side_chain() ) {
		if ( existing_residue.nchi() > 0 &&
				existing_residue.aa() != chemical::aa_pro &&
				existing_residue.n_non_polymeric_residue_connections() == 0 ) {
			dunbrack::RotamerLibraryScratchSpace scratch;
			Size n_min( 0 );
			Real fa_dun_min( 0.0 );
			for ( Size n = 1; n <= rotamers.size(); n++ ) {
				Real const fa_dun = this->rotamer_energy( *rotamers[ n ], scratch );
				if ( n_min == 0 || fa_dun < fa_dun_min ) {
					n_min = n; fa_dun_min = fa_dun;
				}
				// hack for now. PackerTask must be instantiated with a pose without virtual sidechains.
				runtime_assert( !rotamers[n]->has_variant_type( chemical::VIRTUAL_SIDE_CHAIN ) );
			}
			conformation::ResidueOP rot = core::pose::add_variant_type_to_residue( *rotamers[ n_min ], chemical::VIRTUAL_SIDE_CHAIN, pose);
			retval.push_back( rot );
		}
	}
	return retval;
}

/// @brief Equality test for equivalence.
/// Two SingleResidueRotamerLibraries test equal if and only if they represent the exact same behavior
bool
SingleResidueRotamerLibrary::operator ==( SingleResidueRotamerLibrary const & ) const {
	// If you're comparing arbitrary SingleResidueRotamerLibrary, chances are they aren't equal.
	// (Override your subclass if this doesn't work for you.)
	TR.Warning << "[ WARNING ] Program is trying to compare two arbitrary SingleResidueRotamerLibraries - this is probably a bug." << std::endl;
	return false;
}

} // namespace rotamers
} // namespace pack
} // namespace core
