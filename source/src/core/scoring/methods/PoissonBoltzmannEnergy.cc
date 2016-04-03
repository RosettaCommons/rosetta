// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/PoissonBoltzmannEnergy.cc
/// @brief  Ramachandran energy method class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/methods/PoissonBoltzmannEnergy.hh>
#include <core/scoring/methods/PoissonBoltzmannEnergyCreator.hh>

// Package Headers
#include <core/scoring/PoissonBoltzmannPotential.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/OneToAllEnergyContainer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/pb_potential.OptionKeys.gen.hh>

// numeric
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.PoissonBoltzmannEnergy" );

namespace core {
namespace scoring {
namespace methods {

//*****************************************************************
//
// PBLifetimeCache: carries cache for the entire runtime cycle.
//
//*****************************************************************
PBLifetimeCache::PBLifetimeCache()
{}
PBLifetimeCache::~PBLifetimeCache()
{}
basic::datacache::CacheableDataOP
PBLifetimeCache::clone() const
{
	return basic::datacache::CacheableDataOP( new PBLifetimeCache(*this) );
}
void
PBLifetimeCache::set_charged_residues_map( const std::map<std::string, bool> & charged_residues_map )
{
	charged_residues_map_ = charged_residues_map;
}
void
PBLifetimeCache::set_energy_state( const std::string& energy_state )
{
	energy_state_ = energy_state;
}
void
PBLifetimeCache::set_conformational_data( const std::string& energy_state,
	const core::pose::Pose & pose,
	PoissonBoltzmannPotentialOP pbp )
{
	pose_by_state_[energy_state] = core::pose::PoseOP( new core::pose::Pose(pose) );
	pb_by_state_[energy_state] = pbp;
}
std::map<std::string, bool> &
PBLifetimeCache::get_charged_residues_map()
{
	return charged_residues_map_;
}
const std::string &
PBLifetimeCache::get_energy_state() const
{
	return energy_state_;
}

core::pose::PoseCOP
PBLifetimeCache::get_pose( const std::string& energy_state )
{
	//TR << "Looking for cached pose for state \"" << energy_state << "\" in PBLifetimeCache " << "{ ";
	//std::map<std::string, pose::PoseCOP>::const_iterator itr = poses_by_state_.begin();
	//for(; itr != poses_by_state_.end(); ++itr ){
	//  TR << "(key=" << itr->first << ": val=" << itr->second->pdb_info()->name() << "), ";
	//}
	//TR << " }" << std::endl;
	return pose_by_state_[energy_state];
}
PoissonBoltzmannPotentialOP
PBLifetimeCache::get_pbp( const std::string& energy_state )
{
	return pb_by_state_[energy_state];
}

bool
PBLifetimeCache::has_cache( const std::string& energy_state ) const
{
	return pb_by_state_.count(energy_state) && pose_by_state_.count( energy_state);
}
//*****************************************************************
//
// PoissonBoltzmannEnergyCreator
//
//*****************************************************************
/// @details This must return a fresh instance of the PoissonBoltzmannEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
PoissonBoltzmannEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new PoissonBoltzmannEnergy );
}

ScoreTypes
PoissonBoltzmannEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( PB_elec );
	return sts;
}

//*****************************************************************
//
// PoissonBoltzmannEnergy
//
//*****************************************************************
/// ctor
PoissonBoltzmannEnergy::PoissonBoltzmannEnergy() :
	parent( methods::EnergyMethodCreatorOP( new PoissonBoltzmannEnergyCreator ) ),
	fixed_residue_(1),
	epsilon_(2.0)
	//,
	//poisson_boltzmann_potential_( new scoring::PoissonBoltzmannPotential )

{
	if (  basic::options::option[basic::options::OptionKeys::pb_potential::epsilon].user() ) {
		epsilon_ = basic::options::option[basic::options::OptionKeys::pb_potential::epsilon];
	}
}

/// clone
EnergyMethodOP
PoissonBoltzmannEnergy::clone() const
{
	return EnergyMethodOP( new PoissonBoltzmannEnergy( *this ) );
}
methods::LongRangeEnergyType
PoissonBoltzmannEnergy::long_range_type() const { return methods::PB_elec_lr; }

void
PoissonBoltzmannEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const & scorefxn
) const {
	using namespace methods;

	// If the cached object is not in the pose data-cache, it suggests that the setup protocol has not
	// been called in prior.  Warn and get out!
	if ( ! pose.data().has( pose::datacache::CacheableDataType::PB_LIFETIME_CACHE ) ) {
		TR << "PB_LIFETIME_CACHE object is not initialized.  Did you call SetupPoissonBoltzmannPotential mover?  Terminaing the program..." << std::endl;
		TR.flush();
		// Register the empty cache holder if not done so yet.
		PoissonBoltzmannEnergy::PBLifetimeCacheOP new_cache( new PoissonBoltzmannEnergy::PBLifetimeCache() );
		pose.data().set( pose::datacache::CacheableDataType::PB_LIFETIME_CACHE, new_cache );

		//runtime_assert(false);
	}
	PBLifetimeCacheOP cached_data =
		static_cast< PBLifetimeCacheOP > (pose.data().get_ptr< PBLifetimeCache >(pose::datacache::CacheableDataType::PB_LIFETIME_CACHE ));

	// Do we have the map of charged residues by name?
	charged_residues_ = cached_data->get_charged_residues_map();
	if ( charged_residues_.size() == 0 ) {
		utility::vector1<int> charged_chains
			= basic::options::option[basic::options::OptionKeys::pb_potential::charged_chains];
		TR << "Charged residues: [ ";
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			core::conformation::Residue const & rsd( pose.residue(i) );
			bool residue_charged = false;
			if ( std::find(charged_chains.begin(), charged_chains.end(), (int)rsd.chain()) != charged_chains.end() ) {
				residue_charged = true;
				TR << rsd.type().name() << ",";
			}
			charged_residues_[rsd.type().name()] = residue_charged;
		}
		cached_data->set_charged_residues_map( charged_residues_ );
		TR << "]" << std::endl;
	}

	core::scoring::methods::EnergyMethodOptions emoptions = scorefxn.energy_method_options();
	std::string energy_state = cached_data->get_energy_state();

	if ( energy_state == "" ) {
		energy_state = "stateless";
	}
	TR << "Energy state: \"" << energy_state << "\" for scorefxn: " << scorefxn.get_name() << std::endl;

	const Size atom_index = 2;  // alpha carbon


	// Solve PB
	//
	// use cached data if the bound comformation hasn't changed.
	if ( cached_data->has_cache( energy_state ) ) {
		core::pose::PoseCOP prev_pose = cached_data->get_pose( energy_state );

		// switch to the state's pb
		poisson_boltzmann_potential_ = cached_data->get_pbp(energy_state);
		debug_assert(poisson_boltzmann_potential_ != 0);

		TR << "Found cached pose for state: " << energy_state << std::endl;

		// re-evaluate only for bound-state.
		if ( energy_state == emoptions.pb_bound_tag() ) {

			if ( !protein_position_equal_within( pose, *prev_pose, atom_index, epsilon_ ) ) {
				TR << "Atoms (" << atom_index << ") in charged chains moved more than " ;
				TR << epsilon_ << "A" << std::endl;
				poisson_boltzmann_potential_->solve_pb(pose, energy_state, charged_residues_);
			}
		}
	} else {
		TR << "No cached pose for state: " << energy_state << std::endl;
		poisson_boltzmann_potential_ = scoring::PoissonBoltzmannPotentialOP( new PoissonBoltzmannPotential() );
		poisson_boltzmann_potential_->solve_pb(pose, energy_state, charged_residues_);
	}


	// Update the cache
	cached_data->set_conformational_data( energy_state, pose, poisson_boltzmann_potential_ );

	// create LR energy container
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		OneToAllEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::OneToAllEnergyContainer > ( lrc ) );
		// make sure size or root did not change
		if ( dec->size() != pose.total_residue() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		TR << "Creating new one-to-all energy container (" << pose.total_residue() << ")" << std::endl;
		LREnergyContainerOP new_dec( new OneToAllEnergyContainer( fixed_residue_, pose.total_residue(), PB_elec ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}


/////////////////////////////////////////////////////////////////////////////
// methods
/////////////////////////////////////////////////////////////////////////////

bool PoissonBoltzmannEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size res1,
	Size res2
) const {
	return ( res1 == fixed_residue_ || res2 == fixed_residue_ );
}

void
PoissonBoltzmannEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const {
	return;
}

bool
PoissonBoltzmannEnergy::residue_in_chains(
	conformation::Residue const & rsd,
	utility::vector1 <Size> chains
) const {
	for ( Size ichain=1; ichain<=chains.size(); ++ichain ) {
		if ( rsd.chain() == chains[ichain] ) return true;
	}
	return false;
}

Real
PoissonBoltzmannEnergy::revamp_weight_by_burial(
	conformation::Residue const & rsd,
	pose::Pose const & pose
) const {
	utility::vector1 <Size> chains = basic::options::option[basic::options::OptionKeys::pb_potential::revamp_near_chain]();
	Real weight = 1.;
	Real neighbor_count = 0.;
	Real threshold = 4.;
	for ( Size j_res = 1; j_res <= pose.total_residue(); ++j_res ) {
		if ( pose.residue(j_res).is_virtual_residue() ) continue;
		if ( !residue_in_chains(pose.residue(j_res), chains) ) continue;

		bool found_neighbor = false;
		for ( Size i_atom = rsd.last_backbone_atom() + 1; i_atom <= rsd.nheavyatoms(); ++i_atom ) {
			for ( Size j_atom = 1; j_atom <= rsd.nheavyatoms(); ++j_atom ) {
				if ( pose.residue(j_res).is_virtual(j_atom) ) continue;
				Real distance = rsd.xyz(i_atom).distance(pose.residue(j_res).xyz(j_atom));
				if ( distance < threshold ) {
					neighbor_count += 1.;
					found_neighbor = true;
					break;
				}
			}
			if ( found_neighbor ) break;
		}
	}
	if ( neighbor_count > 0 ) weight = 1./neighbor_count;
	//using namespace ObjexxFCL::format;
	//TR << "PB_weight:" << I(4,rsd.seqpos()) << F(8, 3, weight) << std::endl;
	return weight;
}

void
PoissonBoltzmannEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	//check fixed_residue_
	conformation::Residue const &rsd (rsd1.seqpos() == fixed_residue_? rsd2 : rsd1 );

	// If it is part of charged chain, skip.
	if ( const_cast<std::map<std::string, bool>&>(charged_residues_)[rsd.type().name()] ) {
		return;
	}

	Real PB_score_residue, PB_score_backbone, PB_score_sidechain;
	Real PB_burial_weight(1.0);
	if ( basic::options::option[basic::options::OptionKeys::pb_potential::revamp_near_chain].user() ) {
		PB_burial_weight = revamp_weight_by_burial(rsd, pose);
	}
	poisson_boltzmann_potential_->eval_PB_energy_residue(
		rsd,
		PB_score_residue,
		PB_score_backbone,
		PB_score_sidechain,
		PB_burial_weight );

	if ( basic::options::option[basic::options::OptionKeys::pb_potential::sidechain_only]() ) {
		emap[ PB_elec ] += PB_score_sidechain;
	} else {
		emap[ PB_elec ] += PB_score_residue;
	}
}

/// @brief Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
PoissonBoltzmannEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}


/// Compare if two poses are close in fold within tolerance.
///
/// To be specific, it returns True if for all a in A and b in B,
/// sqrt( (a.x-b.x)^2 + (a.y-b.y)^2 + (a.z-b.z)^2 ) <= tol,
/// for A and B are sets of all atoms in protein 1 and protein 2, respectively.
///
/// @param pose1 A protein's pose
/// @param pose2 Another protein's pose
/// @param atom_num The atom number
/// @param tol   Tolerable distance in Angstrom, >= 0.A
///
bool
PoissonBoltzmannEnergy::protein_position_equal_within(
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	Size atom_num,
	Real tol) const {

	// error: reference cannot be bound to dereferenced null pointer in well-defined C++ code; comparison may be assumed to always evaluate to false [-Werror,-Wtautological-undefined-compare]
	//if( &pose1 == NULL || &pose2 == NULL ) {
	// return false;
	//}
	//if( &pose1 == NULL && &pose2 == NULL ) {
	// return true;
	//}

	if ( pose1.total_residue() != pose2.total_residue() ) {
		return false;
	}

	for ( Size i=1; i<=pose1.total_residue(); ++i ) {

		// Skip nonprotein pose 1.
		if ( !pose1.residue(i).is_protein() ) continue;

		// If a protein residue from pose 1 matches seqpos of a nonprotein
		// residue from pose 2, must be false.
		if ( !pose2.residue(i).is_protein() ) return false;

		// Nonmatching by definition.
		if ( pose1.residue_type(i).natoms() != pose2.residue(i).natoms()
				|| pose1.residue_type(i).aa() != pose2.residue_type(i).aa() ) {
			return false;
		}

		const core::conformation::Atom atom1 = pose1.residue(i).atom(atom_num);
		const core::conformation::Atom atom2 = pose2.residue(i).atom(atom_num);
		const numeric::xyzVector<Real> xyz1 = atom1.xyz();
		const numeric::xyzVector<Real> xyz2 = atom2.xyz();

		const Real delta = xyz1.distance(xyz2);
		if ( delta > tol ) {
			// exceed tolerance
			return false;
		}
	}
	return true;
}

core::Size
PoissonBoltzmannEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

