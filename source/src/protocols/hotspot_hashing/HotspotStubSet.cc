// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pose/hotspot_hashing/HotspotStubSet.cc
/// @brief  HotspotStubSet class
/// @author Jacob Corn (jecorn@u.washington.edu), John Karanicolas, Sarel Fleishman

// Unit headers
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/id/AtomID_Map.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/chemical/AA.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/TenANeighborGraph.hh>

#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>


#include <protocols/cluster/APCluster.hh>


#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/constraints/BackboneStubLinearConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>


#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <core/types.hh>

#include <core/pose/Remarks.hh>
#include <core/conformation/Atom.hh>
#include <numeric/random/random.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <sstream>
#include <map>
#include <set>
#include <utility/assert.hh>
#include <algorithm>
#include <cmath>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <boost/foreach.hpp>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <utility/vector0.hh>
#include <basic/options/option.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>


using basic::T;
using basic::Error;
using basic::Warning;

using ObjexxFCL::lead_zero_string_of;


namespace protocols {
namespace hotspot_hashing {
typedef platform::Size Size;
static thread_local basic::Tracer TR( "protocols.hotspot_hashing" );

HotspotStubSet::HotspotStubSet() :
	ReferenceCount(),
	sc_only_(true),
	target_resnum_(0),
	target_distance_(15.0),
	score_threshold_(-1.0),
	filter_( protocols::filters::FilterCOP( protocols::filters::FilterOP( new protocols::filters::TrueFilter ) ) ),
	hotspot_length_( 1 )
{ set_chain_to_design(); }

HotspotStubSet::HotspotStubSet( HotspotStubSet const & init ) :
	ReferenceCount( init ),
	stub_set_( init.stub_set_ ),
	stub_set_vec_( init.stub_set_vec_ ),
	sc_only_( init.sc_only_ ),
	target_resnum_( init.target_resnum_ ),
	target_distance_( init.target_distance_ ),
	chain_to_design_( init.chain_to_design_ ),
	score_threshold_( init.score_threshold_ ),
	filter_( init.filter_ ),
	hotspot_length_( init.hotspot_length_ )
{}

HotspotStubSet::~HotspotStubSet() {}

void HotspotStubSet::clear() {
	stub_set_.clear();
	handshake_stub_sets();
}

Size HotspotStubSet::hotspot_length( ) const { return hotspot_length_; }
void HotspotStubSet::hotspot_length( core::Size const length ) {
	runtime_assert( length > 0 );
	hotspot_length_ = length;
}
void HotspotStubSet::filter( protocols::filters::FilterCOP filter ) { filter_ = filter; }
bool HotspotStubSet::sc_only() const { return sc_only_; }
void HotspotStubSet::sc_only( bool const sc_switch ) { sc_only_ = sc_switch; }

void HotspotStubSet::add_stub_set( HotspotStubSet const & stubset ){
	BOOST_FOREACH ( Hs_data const hs_data, stubset ) add_stub_( hs_data.second.second );
}

void HotspotStubSet::score_threshold( core::Real const threshold ) { score_threshold_ = threshold; }

/// @brief returns a new stub_set with stub scores recalculated by colony energy (Xiang, Soto, and Honig PNAS 2002)
/// @details E = -ln (sum exp( -E(j) - rmsd(ij)^3/6L ) )
HotspotStubSetOP HotspotStubSet::colonyE( ) {
	HotspotStubSetOP new_set( new HotspotStubSet );
	core::pose::PoseOP nonconstpose; // for making our new stubs.
	if ( pose_ ) nonconstpose = core::pose::PoseOP( new core::pose::Pose( *pose_ ) );

	utility::vector1< std::string > amino_acids;
	amino_acids.push_back( "ALA" );
	amino_acids.push_back( "ARG" );
	amino_acids.push_back( "ASN" );
	amino_acids.push_back( "ASP" );
	amino_acids.push_back( "GLU" );
	amino_acids.push_back( "GLN" );
	amino_acids.push_back( "HIS" );
	amino_acids.push_back( "ILE" );
	amino_acids.push_back( "LEU" );
	amino_acids.push_back( "LYS" );
	amino_acids.push_back( "MET" );
	amino_acids.push_back( "PHE" );
	amino_acids.push_back( "SER" );
	amino_acids.push_back( "THR" );
	amino_acids.push_back( "TRP" );
	amino_acids.push_back( "TYR" );
	amino_acids.push_back( "VAL" );

	TR << "Calculating colony energy..." << std::endl;
	core::Size const nres(1); // number of residues in our "loop" (single residue stub) this is nonsense in our case, but equates to a constant so shouldn't matter
	for ( utility::vector1< std::string >::const_iterator it = amino_acids.begin(); it != amino_acids.end(); ++it ) {
		std::string const resname = *it;
		// Loop over all stubs with this restype
		Hotspots res_stub_set( retrieve( resname ) );
		utility::vector1< HotspotStubCOP > stub_vec;
		for ( std::multimap<core::Real, HotspotStubOP >::const_iterator i = res_stub_set.begin(); i != res_stub_set.end(); ++i ) {
			stub_vec.push_back( i->second );
		}
		if ( stub_vec.size() == 0 ) continue; // in case we don't have any stubs of that type

		for ( std::multimap<core::Real, HotspotStubOP >::const_iterator i = res_stub_set.begin(); i != res_stub_set.end(); ++i ) {
			core::Real stubi_E(0);
			HotspotStubCOP stubi = i->second;
			for ( std::multimap<core::Real, HotspotStubOP >::const_iterator j = i; j != res_stub_set.end(); ++j ) {
				HotspotStubCOP stubj = j->second;
				core::Real const rms = residue_sc_rmsd_no_super( stubi->residue(), stubj->residue() );
				stubi_E += exp( 0-stubi->bonus_value() - pow(rms,3)/6*nres );
			}

			stubi_E = 0-log( stubi_E );
			HotspotStubCOP new_stub( HotspotStubOP( new HotspotStub( stubi->residue(), stubi_E, nonconstpose, chain_to_design_, filter_ ) ) );
			new_set->add_stub_( new_stub );
		}
	}
	return new_set;
}

/// @brief clusters all residues within this stubset on a per-restype basis and returns a new clustered stubset
HotspotStubSetOP HotspotStubSet::cluster( ) {

	bool const fxnal_group = basic::options::option[ basic::options::OptionKeys::hotspot::fxnal_group ]();

	HotspotStubSetOP new_set( new HotspotStubSet );
	utility::vector1< std::string > amino_acids;
	amino_acids.push_back( "ALA" );
	amino_acids.push_back( "ARG" );
	amino_acids.push_back( "ASN" );
	amino_acids.push_back( "ASP" );
	amino_acids.push_back( "GLU" );
	amino_acids.push_back( "GLN" );
	amino_acids.push_back( "HIS" );
	amino_acids.push_back( "ILE" );
	amino_acids.push_back( "LEU" );
	amino_acids.push_back( "LYS" );
	amino_acids.push_back( "MET" );
	amino_acids.push_back( "PHE" );
	amino_acids.push_back( "SER" );
	amino_acids.push_back( "THR" );
	amino_acids.push_back( "TRP" );
	amino_acids.push_back( "TYR" );
	amino_acids.push_back( "VAL" );

	for ( utility::vector1< std::string >::const_iterator it = amino_acids.begin(); it != amino_acids.end(); ++it ) {
		std::string const resname = *it;
		// Loop over all stubs with this restype
		utility::vector1< HotspotStubCOP > stub_vec;
		{
			Hotspots res_stub_set( retrieve( resname ) );
			for ( std::multimap<core::Real, HotspotStubOP >::const_iterator i = res_stub_set.begin(); i != res_stub_set.end(); ++i ) {
				stub_vec.push_back( i->second );
			}
		}
		if ( stub_vec.size() == 0 ) continue; // in case we don't have any stubs of that type

		utility::vector1< core::Real > flat_rms; // used to determine minimum rmsd for APCluster self-similarity
		flat_rms.reserve( (core::Size)pow(stub_vec.size(), 2.0f) );

		protocols::cluster::APCluster apcluster( stub_vec.size() );

		core::Size stubidx_i( 1 );
		// using i,k per APCluster notation
		for ( utility::vector1< HotspotStubCOP >::const_iterator i = stub_vec.begin(); i != stub_vec.end(); ++i ) {
			HotspotStubCOP stub1 = *i;
			core::Size stubidx_k( stubidx_i );
			for ( utility::vector1< HotspotStubCOP >::const_iterator k = i; k != stub_vec.end(); ++k ) {
				HotspotStubCOP stub2 = *k;
				core::Real const rms = residue_sc_rmsd_no_super( stub1->residue(), stub2->residue(), fxnal_group ); // do clustering only based on functional group
				TR.Debug << "rmsd " << stubidx_i << " " << stubidx_k << " " << rms << std::endl;

				flat_rms.push_back( rms );

				apcluster.set_sim( stubidx_i, stubidx_k, 0-rms ); // use negative rmsd for clustering (larger values are closer)
				apcluster.set_sim( stubidx_k, stubidx_i, 0-rms ); // go both ways for symmetry
				++stubidx_k;
			}
			++stubidx_i;
		}

		// find median rmsd
		sort( flat_rms.begin(), flat_rms.end() );
		core::Real const median_rmsd = *(flat_rms.begin() + flat_rms.size() / 2);
		//core::Real const min_rmsd = *(std::min_element( flat_rms.begin(), flat_rms.end() ) );
		for ( core::Size i = 1; i <= stub_vec.size(); ++i ) {
			apcluster.set_sim( i, i, 0 - median_rmsd ); // use median_rmsd for self-similarity (see APCluster.hh). Better to use stub score to weight for better stubs?
			TR.Debug << "Self-sim " << i << " " << 0-median_rmsd  << std::endl;
		}
		/*
		for( core::Size i = 1; i <= all_rms.size(); ++i ) {
		if( all_rms[i][1]==all_rms[i][2] ) { // if self-comparison (see APCluster.hh for special meaning)
		//TR << all_rms[i][1] << " " << all_rms[i][2] << " " << 0 - median_rmsd << std::endl; // debug output

		}
		else {
		//TR << all_rms[i][1] << " " << all_rms[i][2] << " " << 0 - all_rms[i][3] << std::endl; // debug output
		apcluster.set_sim( (core::Size)all_rms[i][1], (core::Size) all_rms[i][2], 0 - all_rms[i][3] ); // use negative rmsd for clustering (larger values are closer)
		apcluster.set_sim( (core::Size)all_rms[i][2], (core::Size)all_rms[i][1], 0 - all_rms[i][3] ); // go both ways for symmetry
		}
		}
		*/
		// these settings are relatively low. clustering may not converge completely, but should be OK
		core::Size maxits( 1000 );
		core::Size convits( 10 );
		core::Real lambda( 0.75 );
		TR << "Clustering " << resname << " stubs..." << std::endl;
		bool cluster_success = apcluster.cluster( maxits, convits, lambda );
		if ( !cluster_success ) TR << "Hotspot APClustering did not fully converge... This happens all the time, but just FYI." << std::endl;

		//the following  was used to test convergence of clustering. in the end, just using the max values as defaults worked best
		/*
		while( !cluster_success ) {
		TR << "Hotspot APClustering did not fully converge... This happens all the time, but just FYI." << std::endl;
		maxits += 500;
		convits += 100;
		lambda += 0.1;
		TR << "maxits=" << maxits << " convits=" << convits << " lambda=" << lambda << std::endl;
		cluster_success = apcluster.cluster( maxits, convits, lambda );
		if( maxits >= 4000 || convits >= 400 || lambda >= 0.95 ) {
		TR << "Hotspot APClustering never fully converged. Using current cluster centers." << std::endl;
		break;
		}
		}
		*/

		TR << "Finding low energy cluster members..." << std::endl;
		utility::vector1< core::Size > exemplars; // all exembplars
		utility::vector1< core::Size > lowE_members; // lowest energy cluster member for each exemplar
		apcluster.get_all_exemplars( exemplars );
		// find the lowest energy cluster member for each exemplar
		for ( utility::vector1<core::Size>::const_iterator it = exemplars.begin(); it != exemplars.end(); ++it ) {
			core::Size lowE_index(0);
			utility::vector1< core::Size > cluster_members;
			apcluster.get_cluster_for(*it, cluster_members);
			core::Real best_score(0);
			for ( utility::vector1<core::Size>::const_iterator it = cluster_members.begin(); it != cluster_members.end(); ++it ) {
				HotspotStubCOP stub = stub_vec[*it];
				core::Real const score = stub->bonus_value();
				if ( score < best_score ) {
					lowE_index = *it;
					best_score = score;
					TR.Debug << lowE_index << " " << stub->bonus_value() << " " << score << " " << best_score << std::endl;
				}
			}
			assert( lowE_index != 0 );
			lowE_members.push_back( lowE_index );
		}

		//TR << exemplars << std::endl; // debug output
		for ( utility::vector1<core::Size>::const_iterator it = lowE_members.begin(); it != lowE_members.end(); ++it ) {
			// TR << "Adding stub number " << *it << std::endl; // debug output
			new_set->add_stub_( (stub_vec[*it]) );
		}
	}
	return new_set;
}

// @brief retrieve stubs with a given residue name3
HotspotStubSet::Hotspots HotspotStubSet::retrieve( std::string const & resname ) {
	std::map< std::string, std::multimap< core::Real, HotspotStubOP > >::iterator ss_iter;
	ss_iter = stub_set_.find( resname );
	if ( ss_iter == stub_set_.end() ) {
		Hotspots empty_set;
		stub_set_.insert( std::make_pair( resname, empty_set ) );
		ss_iter = stub_set_.find( resname );
	}
	TR.Debug << "Retrieving " << ss_iter->second.size() << " " << resname << " stubs." << std::endl;
	return (ss_iter->second);
}

// @brief retrieve a subset of stubs with a given residue name3 and score cutoff. If scorecut is negative, treat it as a score. If scorecut is positive, treat it as a fraction.
HotspotStubSetOP HotspotStubSet::subset( std::string const & residue_name3, core::Real const scorecut ) {

	HotspotStubSetOP new_set( new HotspotStubSet );
	new_set->sc_only_ = sc_only_;

	Hotspots all_stubs = retrieve(residue_name3);
	// make a dummy stub with the same residue as the first stub and a user-supplied score
	//HotspotStub scorecut_stub( all_stubs.begin()->residue(), scorecut );

	if ( ( scorecut > 0 ) && ( scorecut <=1 ) ) {
		// add 0.5 to get rounding up
		Size n_return = static_cast<Size>(all_stubs.size() * scorecut + 0.5);
		if ( n_return < 1 ) n_return = 1;
		TR << "Finding the top " << n_return << " stubs." << std::endl;
		Size i = 1;
		for ( Hotspots::const_iterator stub_iter = all_stubs.begin(); stub_iter != all_stubs.end() ; ++stub_iter ) {
			if ( i <= n_return ) {
				//core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_BB", pose.total_residue() );
				new_set->add_stub_( stub_iter->second );
				++i;
			} else break;
		}
	} else if ( scorecut <=0 ) {
		for ( Hotspots::const_iterator stub_iter = all_stubs.begin(); stub_iter != all_stubs.upper_bound( scorecut ) ; ++stub_iter ) {
			if ( stub_iter->second->bonus_value() <= scorecut ) {
				new_set->add_stub_( stub_iter->second );
			}
		}
	}
	TR << "Found " << new_set->size() << " stubs better than " << scorecut << std::endl;
	return new_set;
}


// @brief retrieve all stubs that pass a score cutoff. If scorecut is negative, treat it as a score. If scorecut is positive, treat it as a fraction.
HotspotStubSetOP HotspotStubSet::subset( core::Real const scorecut ) const {

	HotspotStubSetOP new_set( new HotspotStubSet );
	new_set->sc_only_ = sc_only_;
	// make a dummy stub with the same residue as the first stub and a user-supplied score
	//HotspotStub scorecut_stub( stub_set_.begin()->second.begin()->residue(), scorecut );

	for ( Hs_map::const_iterator ss_iter= stub_set_.begin(); ss_iter != stub_set_.end(); ++ss_iter ) {

		if ( ( scorecut > 0) && ( scorecut <= 1 ) ) {
			// add 0.5 to get rounding up
			Size n_return = static_cast<Size>(ss_iter->second.size() * scorecut + 0.5);
			TR << "Finding the top " << n_return << " stubs." << std::endl;
			if ( n_return < 1 ) n_return = 1;
			Size i = 1;
			for ( Hotspots::const_iterator stub_iter = ss_iter->second.begin(); stub_iter != ss_iter->second.end(); ++stub_iter ) {
				if ( i <= n_return ) {
					new_set->add_stub_( stub_iter->second );
					++i;
				}
			}
		} else if ( scorecut <= 0 ) {
			for ( Hotspots::const_iterator stub_iter = ss_iter->second.begin(); stub_iter != ss_iter->second.upper_bound( scorecut ) ; ++stub_iter ) {
				//    if( stub_iter->bonus_value() <= scorecut ) {
				new_set->add_stub_( stub_iter->second );
				//    }
			}
		}
	}
	TR << "Found " << new_set->size() << " stubs better than " << scorecut << std::endl;
	return new_set;
}


HotspotStubCOP
HotspotStubSet::get_best_energy_stub() const {
	core::Real min_energy( 100000.0 );
	HotspotStubOP ret( NULL );
	for ( Hs_map::const_iterator hs_map_it=stub_set_.begin(); hs_map_it!=stub_set_.end(); ++hs_map_it ) {
		//typedef std::multimap< core::Real, HotspotStubOP > Hs_multimap;
		for ( Hotspots::const_iterator hs_it=hs_map_it->second.begin(); hs_it!=hs_map_it->second.end(); ++hs_it ) {
			if ( min_energy > hs_it->first ) {
				min_energy = hs_it->first;
				ret = HotspotStubOP( new HotspotStub( *hs_it->second ) );
			}
		}
	}
	runtime_assert( ret != 0 );
	return( ret );
}

/// @details find stub nearest to residue based on CA-CA distance
HotspotStubCOP
HotspotStubSet::get_nearest_stub( core::conformation::ResidueCOP residue ) const {
	using namespace core::conformation;

	HotspotStubCOP nearest_stub( NULL );
	core::Real nearest_distance( 100000.0 );
	numeric::xyzVector< core::Real > const residue_CA( residue->xyz( "CA" ) );
	for ( Hs_map::const_iterator map_it=stub_set_.begin(); map_it!=stub_set_.end(); ++map_it ) {
		//typedef std::multimap< core::Real, HotspotStubOP > Hs_multimap;
		for ( Hotspots::const_iterator stub_it=map_it->second.begin(); stub_it!=map_it->second.end(); ++stub_it ) {
			numeric::xyzVector< core::Real > const stub_CA( stub_it->second->residue()->xyz( "CA" ) );

			core::Real const distance( stub_CA.distance( residue_CA ) );
			if ( distance <= nearest_distance ) {
				nearest_distance = distance;
				nearest_stub = stub_it->second;
			}
		}
	}
	runtime_assert( nearest_stub != 0 );
	return( nearest_stub );
}

std::set< std::pair< std::string, core::Real > >
HotspotStubSet::find_neighboring_stubs( HotspotStubCOP stub ) const {
	using namespace core::conformation;

	Residue const stub_rsd( *stub->residue() );

	std::set< std::pair< std::string, core::Real > > stub_subset;

	core::Real const dist_threshold( 3.0 );

	for ( Hs_map::const_iterator hs_map_it=stub_set_.begin(); hs_map_it!=stub_set_.end(); ++hs_map_it ) {
		//typedef std::multimap< core::Real, HotspotStubOP > Hs_multimap;
		for ( Hotspots::const_iterator hs_it=hs_map_it->second.begin(); hs_it!=hs_map_it->second.end(); ++hs_it ) {
			Residue const potential_neighbor( *(hs_it->second->residue()) );
			for ( Atoms::const_iterator stub_atom_it=stub_rsd.atom_begin(); stub_atom_it!=stub_rsd.atom_end(); ++stub_atom_it ) {
				for ( Atoms::const_iterator pot_neigh_it=potential_neighbor.atom_begin(); pot_neigh_it!=potential_neighbor.atom_end(); ++pot_neigh_it ) {
					core::Real const dist( stub_atom_it->xyz().distance( pot_neigh_it->xyz() ) );
					if ( dist <= dist_threshold ) {
						stub_subset.insert( std::make_pair( hs_map_it->first, hs_it->first ) );
						break;
					}
				}
			}
		}
	}
	return( stub_subset );
}

HotspotStubSet::Hotspots::const_iterator
HotspotStubSet::get_stub( std::string const residue_name3, core::Real const score ) const
{
	Hs_map::const_iterator hs_it( stub_set_.find( residue_name3 ) );
	return( hs_it->second.find( score ) );
}

/// @details removes the first occurence of stub in the stubset
bool
HotspotStubSet::remove_stub( HotspotStubCOP stub ){
	for ( Hs_map::iterator datum( stub_set_.begin() ); datum!=stub_set_.end(); ++datum ) {
		Hotspots::iterator it( datum->second.begin() );
		for ( ; it!= datum->second.end(); ++it ) {
			if ( it->second == stub ) break;
		}
		if ( it != datum->second.end() ) {
			datum->second.erase( it );
			handshake_stub_sets();
			return( true );
		}
	}
	TR.Warning<<"WARNING: Stub "<<stub<<" not found in remove stub"<<std::endl;
	return( false );
}

/// @details removes a set of stubs from the stub_set_
void
HotspotStubSet::remove_stubs_from_set( std::set< std::pair< std::string, core::Real > > const stubs ) {
	for ( std::set< std::pair< std::string, core::Real > >::const_iterator stub_it=stubs.begin(); stub_it!=stubs.end(); ++stub_it ) {
		Hs_map::iterator hs_it( stub_set_.find( stub_it->first ) );
		hs_it->second.erase( get_stub( stub_it->first, stub_it->second )->first );
	}
	handshake_stub_sets();
}

core::scoring::constraints::ConstraintCOPs HotspotStubSet::constraints() const {
	return constraints_;
}

void HotspotStubSet::add_stub_( HotspotStubCOP stub ) {
	std::string const resname = stub->residue()->name3();
	Hs_map::iterator ss_iter;
	ss_iter = stub_set_.find( resname );
	if ( ss_iter == stub_set_.end() ) {
		Hotspots empty_set;
		stub_set_.insert( std::make_pair( resname, empty_set ) );
		ss_iter = stub_set_.find( resname );
	}
	(ss_iter->second).insert( std::make_pair(stub->bonus_value(), HotspotStubOP( new HotspotStub( *stub ) ) ) );
	handshake_stub_sets();
}

void HotspotStubSet::read_data( std::string const filename ) {
	// keep PDB header to preserve REMARKs
	basic::options::option[ basic::options::OptionKeys::run::preserve_header ].value(true);

	// read all the poses in filename
	utility::vector1< core::pose::Pose > poses;
	core::import_pose::pose_from_pdb( poses, filename );

	core::pose::PoseOP nonconstpose( NULL );//SJF pose will be attached later
	for ( utility::vector1<core::pose::Pose>::iterator it = poses.begin(); it!= poses.end(); ++it ) {
		// get REMARKS associated with the pose
		core::pose::PDBInfoCOP pdbinfo = it->pdb_info();
		core::pose::Remarks const & remarks ( pdbinfo->remarks() );
		core::Real score = 0;
		for ( std::vector< core::pose::RemarkInfo >::const_iterator remark_it = remarks.begin(); remark_it != remarks.end(); ++remark_it ) {
			// special remark code for theoretical scores
			if ( remark_it->num == 221 ) {
				score = std::atof( remark_it->value.c_str() );
				if ( score >= -0.0001 ) {
					TR<<"****WARNING WARNING**** stub score has zero or higher weight. Reading nonetheless."<<std::endl;
				}
				break;
			}
		}

		// only keep the pose/residue if we found a valid score
		//if ( score < 0 )
		//{
		// make a stub from the last residue on the pose (which should only have one residue) and add it to the set
		//using namespace core::chemical;
		core::conformation::ResidueCOP residue( core::conformation::ResidueOP( new core::conformation::Residue( it->residue( it->total_residue() ) ) ) );
		HotspotStubCOP stub( HotspotStubOP( new HotspotStub( residue, score, nonconstpose, chain_to_design_, filter_ ) ) );
		add_stub_( stub );
		//}
	}
	TR << "Read " << size() << " stubs from " << filename << std::endl;
	return;
}

void HotspotStubSet::remove_random_stubs_from_set( int const num_to_remove ){
	if ( num_to_remove <= 0 ) return;
	if ( num_to_remove >= ( int ) size() ) {
		TR<<"ERROR: Trying to remove "<< num_to_remove<<" stubs from a set containing "<<size()<<" stubs."<<std::endl;
		utility_exit();
	}

	std::vector< core::Size > to_remove;
	for ( core::Size i = 1; i <= size(); ++i ) {
		to_remove.push_back( i );
	}

	//std::random__shuffle( to_remove.begin(), to_remove.end() );
	numeric::random::random_permutation( to_remove.begin(), to_remove.end(), numeric::random::rg() );

	std::vector< HotspotStubOP > stubs_to_remove;
	for ( int i = 1; i <= num_to_remove; ++i ) {
		stubs_to_remove.push_back( stub_set_vec_[ i ].second.second );
	}

	BOOST_FOREACH ( HotspotStubOP hs, stubs_to_remove ) {
		remove_stub( hs );
	}

	TR<<"Removed stubs from set. Current number "<<size()<<std::endl;
}

void HotspotStubSet::autofill( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn, Size const n_stubs )
{
	utility::vector1< std::string > amino_acids;
	amino_acids.push_back( "ALA" );
	amino_acids.push_back( "ARG" );
	amino_acids.push_back( "ASN" );
	amino_acids.push_back( "ASP" );
	amino_acids.push_back( "GLU" );
	amino_acids.push_back( "GLN" );
	amino_acids.push_back( "HIS" );
	amino_acids.push_back( "ILE" );
	amino_acids.push_back( "LEU" );
	amino_acids.push_back( "LYS" );
	amino_acids.push_back( "MET" );
	amino_acids.push_back( "PHE" );
	amino_acids.push_back( "SER" );
	amino_acids.push_back( "THR" );
	amino_acids.push_back( "TRP" );
	amino_acids.push_back( "TYR" );
	amino_acids.push_back( "VAL" );
	for ( utility::vector1<std::string>::iterator it = amino_acids.begin(); it != amino_acids.end(); ++it ) {
		fill( pose, scorefxn, *it, n_stubs );
	}
	return;
}

void HotspotStubSet::autofill( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn, core::Size const target, core::Real const distance, Size const n_stubs )
{
	if ( ( target <= pose.total_residue() ) &&  ( distance > 0 ) ) {
		target_resnum_ = target;
		target_distance_ = distance;

		autofill( pose, scorefxn, n_stubs );
	} else {
		utility_exit_with_message( "Unable to find target residue, or distance less than zero!\n");
	}

	return;

}

// convenience fillers
void HotspotStubSet::fill( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn_in, core::Size const target, core::Real const distance, std::string const residue_name3, Size const n_stubs )
{
	if ( ( target <= pose.total_residue() ) &&  ( distance > 0 ) ) {
		target_resnum_ = target;
		target_distance_ = distance;

		fill( pose, scorefxn_in, residue_name3, n_stubs );
	} else {
		utility_exit_with_message( "Unable to find target residue, or distance less than zero!\n");
	}

	return;
}

// MAIN FILLING MACHINERY HERE
void HotspotStubSet::fill( core::pose::Pose const & reference_pose, core::scoring::ScoreFunctionCOP scorefxn_in, std::string const residue_name3, Size const n_stubs )
{

	// set up scorefxn's. DISABLE environment-dependent H-bonds
	core::scoring::ScoreFunctionOP centroid_scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" ) );

	core::scoring::ScoreFunctionOP scorefxn = scorefxn_in->clone();

	/*
	core::scoring::ScoreFunctionOP noenvhbond_scorefxn( core::scoring::get_score_function_legacy( "score13" ) );
	core::scoring::methods::EnergyMethodOptions options( noenvhbond_scorefxn->energy_method_options() );
	options.hbond_options()->use_hb_env_dep( false );
	noenvhbond_scorefxn->set_energy_method_options( options );
	noenvhbond_scorefxn->set_weight( core::scoring::fa_dun, 0.1 );
	noenvhbond_scorefxn->set_weight( core::scoring::envsmooth, 0 );
	core::scoring::ScoreFunctionOP scorefxn( noenvhbond_scorefxn );
	*/

	core::pose::Pose pose = reference_pose;

	for ( Size i = 1; i <= n_stubs; ++i ) {
		core::Real score = 0;
		// some counting for output
		std::stringstream index;
		index.str("");
		index << i;

		// Keep docking until we get a good score
		runtime_assert( score_threshold_ != 0 );
		while ( score > score_threshold_ ) {

			pose = reference_pose;

			// make our hotspot
			create_hotspot_after_pose( pose, residue_name3 );

			TR.Debug << "old fold tree: " << pose.fold_tree() << std::endl;

			setup_hotspot_foldtree_( pose );

			TR.Debug << "Hotspot fold tree: " << pose.fold_tree() << std::endl;

			// placed_seqpos = position of the hotspot residue of peptide
			Size const placed_seqpos = pose.total_residue() - stub_offset(); //(hotspot_length()/2.0 == 0.0) ? pose.total_residue() - (hotspot_length()/2) + 1 : pose.total_residue() - (hotspot_length()/2);
			Size const jump_number = pose.num_jump();
			core::kinematics::FoldTree old_tree = pose.fold_tree(); // need the old tree to reset after docking
			core::pose::PoseOP nonconstpose( new core::pose::Pose( pose ) ); // dummy for making stubs and as a pseudo-native for docking

			// set up docking mover
			// includes an initial perturbation, so we don't need to keep making a new residue each time.
			protocols::docking::DockingProtocolOP innerdock( new protocols::docking::DockingProtocol( jump_number, false /*low_res_protocol_only*/, false /*local_refine*/, false /*set_autofoldtree*/, centroid_scorefxn /*lowres*/, scorefxn /*highres*/ ) );
			innerdock->set_reporting( false ); // avoid rescoring steps
			innerdock->set_no_filters( true ); // don't use docking filters
			/// APL -- reverting a portion of 38676 -- innerdock->set_no_filters( true ); // don't use filters
			protocols::moves::MoverOP dock = innerdock;
			dock->set_input_pose( nonconstpose );
			dock->set_native_pose( nonconstpose );

			//(*scorefxn)(pose);
			if ( target_resnum_ && target_distance_ ) { // randomize residue and translate target_distance
				protocols::rigid::RigidBodyPerturbMover pert( jump_number, 0, target_distance_ );
				protocols::rigid::RigidBodyRandomizeMover rand2( pose, jump_number, protocols::rigid::partner_downstream );
				protocols::docking::FaDockingSlideIntoContact slide_together( jump_number );
				pert.apply( pose );
				rand2.apply( pose );
				slide_together.apply( pose );
			} else { // manually do -randomize1/2 (without commandline flags)
				protocols::rigid::RigidBodyRandomizeMover rand1( pose, jump_number, protocols::rigid::partner_upstream );
				protocols::rigid::RigidBodyRandomizeMover rand2( pose, jump_number, protocols::rigid::partner_downstream );
				protocols::docking::FaDockingSlideIntoContact slide_together( jump_number );
				rand1.apply( pose );
				rand2.apply( pose );
				rand1.apply( pose ); // double-randomize, since I've been seeing weird nonrandom behavior
				rand2.apply( pose );
				slide_together.apply( pose );
			}

			// check initial distance to target
			if ( target_resnum_ && target_distance_ ) {
				// convenience pointers
				core::conformation::ResidueCOP const res_target( pose.residue( target_resnum_ ).get_self_ptr() );
				core::conformation::ResidueCOP const stub( pose.residue( placed_seqpos ).get_self_ptr() );
				// distance check
				core::Real distance( res_target->xyz( res_target->nbr_atom() ).distance( stub->xyz( stub->nbr_atom() )) );
				TR << "InitDist: " << distance << "A from res " << target_resnum_;
				if ( !( distance <= target_distance_+5 ) ) {
					TR << ".  Reject." << std::endl;
					continue;
				} else TR << ". ";
			}

			// reset fold tree from docking modifications. necessary for compatibility with repeated dock->apply calls.
			dock->apply( pose );

			// check final distance to target
			if ( target_resnum_ && target_distance_ ) {
				// convenience pointers
				core::conformation::ResidueCOP const res_target( pose.residue( target_resnum_ ).get_self_ptr() );
				core::conformation::ResidueCOP const stub( pose.residue( placed_seqpos ).get_self_ptr() );
				// distance check
				core::Real distance( res_target->xyz( res_target->nbr_atom() ).distance( stub->xyz( stub->nbr_atom() )) );
				TR << "FinalDist: " << distance << "A from res " << target_resnum_;
				if ( !( distance <= target_distance_ ) ) {
					TR << ".  Reject." << std::endl;
					continue;
				} else TR << ". ";
			}

			// check angle towards CoM
			if ( basic::options::option[ basic::options::OptionKeys::hotspot::angle ].user() ) {
				core::conformation::ResidueCOP const stub( pose.residue( placed_seqpos ).get_self_ptr() );
				core::Real const threshold = basic::options::option[ basic::options::OptionKeys::hotspot::angle ].value();
				runtime_assert( threshold > 0 );
				// angle_res defaults to 0, which causes stub_tgt_angle() to calc and use target center of mass
				core::Size const angle_res = basic::options::option[ basic::options::OptionKeys::hotspot::angle_res ].value();

				core::Real const angle = stub_tgt_angle( pose, stub, angle_res );
				TR << " Stub-target angle: " << angle << "o";
				if ( angle > threshold ) {
					TR << ".  Reject." << std::endl;
					continue;
				} else TR << ". ";
			}

			score = get_residue_score_( pose, scorefxn, placed_seqpos );
			TR << "Stub score: " << score;
			// only keep stub if scores better than/eq a threshold (defaults to -1.0)
			if ( score <= score_threshold_ ) {
				core::conformation::ResidueCOP found_residue( core::conformation::ResidueOP( new core::conformation::Residue( pose.residue( placed_seqpos ) ) ) );
				HotspotStubCOP stub( HotspotStubOP( new HotspotStub( found_residue, score, nonconstpose, chain_to_design_, filter_ ) ) );
				add_stub_( stub );
				TR << ". Accept." << std::endl;
				//pose.dump_pdb( "placed_stub.pdb"); // debug pdb output
				break; // break out of while( score > -1.0)
			} else {
				TR << ". Reject." << std::endl;
			}

		} // end while ( score > threshold )
		TR.flush();
	} // end for n_stubs

	return;
}

HotspotStubSetOP
HotspotStubSet::rescore( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn ) {

	HotspotStubSetOP new_set( new HotspotStubSet );
	new_set->sc_only( sc_only_ );

	using namespace core;
	pose::PoseOP nonconstpose( new core::pose::Pose( pose ) );
	TR << "Rescoring hotspots...\n";
	TR << "Original Rescored\n";
	for ( Hs_map::const_iterator it = stub_set_.begin(); it != stub_set_.end(); ++it ) {
		for ( Hotspots::const_iterator stub_it = it->second.begin(); stub_it != it->second.end(); ++stub_it ) {
			pose::Pose working_pose = pose;
			conformation::ResidueCOP residue = stub_it->second->residue();

			// append the residue
			working_pose.append_residue_by_jump( *residue, working_pose.total_residue(), "", "", true );
			Size const placed_seqpos = working_pose.total_residue();

			// option to use make backbone invisible
			if ( new_set->sc_only() ) {
				core::pose::add_variant_type_to_pose_residue( working_pose, core::chemical::SHOVE_BB, pose.total_residue() );
			}
			core::Real const score = get_residue_score_( working_pose, scorefxn, placed_seqpos );
			TR << stub_it->first << " " << score << "\n";
			HotspotStubCOP stub( HotspotStubOP( new HotspotStub( residue, score, nonconstpose, chain_to_design_, filter_ ) ) );
			new_set->add_stub_( stub );
		}
	}
	TR.flush();
	return new_set;
}

void HotspotStubSet::write_all( std::string const & filename ) const
{
	utility::io::ozstream outstream;
	outstream.open( filename.c_str(), std::ios::app );
	// convenience number. would be better to read the last atom number prior to appending the new residue.
	Size i = 0;
	std::string tag( "" );
	for ( Hs_map::const_iterator it = stub_set_.begin(); it != stub_set_.end(); ++it ) {
		for ( Hotspots::const_iterator stub_it = it->second.begin(); stub_it != it->second.end(); ++stub_it ) {
			tag = "S_" + stub_it->second->residue()->name3() + "_" + lead_zero_string_of( i, 9 );
			write_stub( outstream, stub_it->second, tag );
			++i;
		}
	}
	outstream.close();
	TR << "Wrote " << i << " stubs to " << filename << std::endl;
	return;
}

void HotspotStubSet::write_stub( utility::io::ozstream & outstream, HotspotStubCOP stub, std::string const & tag ) const
{
	// convenience number. would be better to read the last atom number prior to appending the new residue.
	Size atom_number = 10000;

	outstream << "MODEL " << tag << "\n";
	outstream << "REMARK 221 " << stub->bonus_value() << " \n";
	core::io::pdb::dump_pdb_residue( *stub->residue(), atom_number, outstream );
	outstream << "TER\n";
	outstream << "ENDMDL\n";
	return;
}

/// @brief set up scaffold_status_ for each stub included in this set
void HotspotStubSet::pair_with_scaffold( core::pose::Pose const & pose, core::Size const partner, protocols::filters::FilterCOP filter )
{
	pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) );
	chain_to_design_ = partner;
	filter_ = filter;
	core::Size const chain_beg( pose_->conformation().chain_begin( chain_to_design_ ) );
	core::Size const chain_end( pose_->conformation().chain_end( chain_to_design_ ) );
	std::vector< StubStatus > temp_status( chain_end-chain_beg+1, unchecked );

	core::pose::PoseOP nonconstpose( new core::pose::Pose( *pose_ ) );
	for ( Hs_map::iterator set_it=stub_set_.begin(); set_it!=stub_set_.end(); ++set_it ) {
		for ( Hotspots::iterator stub_it = set_it->second.begin(); stub_it != set_it->second.end(); ++stub_it ) {
			// makes sure contained stubs know their StubSet parent
			//   stub_it->second.set_stub_parent_( *this );
			// readies stubs for setting of scaffold_status vector
			stub_it->second->pair_with_scaffold( nonconstpose, filter_, chain_to_design_ );
		}
	}
	TR << "Associated stubs with scaffold chain " << partner << std::endl;
}

/// @brief how many total stubs are in the set (all residues)?
core::Size HotspotStubSet::size() const
{
	core::Size n_stubs(0);
	for ( Hs_map::const_iterator it = stub_set_.begin(); it != stub_set_.end(); ++it ) {
		n_stubs += it->second.size();
	}
	return n_stubs;
}

/// @brief how many stubs are in the set by residue?
core::Size HotspotStubSet::size( std::string const resname )
{
	core::Size n_stubs(0);
	if ( resname != "ALL" ) {
		Hotspots const stubs = retrieve( resname );
		n_stubs = stubs.size();
		return n_stubs;
	} else {
		n_stubs = size();
		return n_stubs;
	}
}

/// @details remove all stubs from stub_set_vec_ and repopulate it with stub_set_
void
HotspotStubSet::handshake_stub_sets(){
	stub_set_vec_.clear();
	Hs_map::iterator stubset_it;
	for ( stubset_it=stub_set_.begin(); stubset_it!=stub_set_.end(); ++stubset_it ) {
		std::multimap<core::Real,HotspotStubOP >::iterator hs_it;
		hs_it = stubset_it->second.begin();
		for ( hs_it = stubset_it->second.begin(); hs_it!=stubset_it->second.end(); ++hs_it ) {
			stub_set_vec_.push_back( std::make_pair( stubset_it->first, std::make_pair( hs_it->first, hs_it->second ) ));
		}
	}
}

HotspotStubOP
HotspotStubSet::random_stub()
{
	core::Size const stubset_size( size() );
	core::Size const random_element( ( core::Size )( numeric::random::rg().uniform() * stubset_size ) + 1 );

	runtime_assert( random_element <= stubset_size );
	core::Size access_point( random_element );
	Hs_map::iterator stubset_it;
	for ( stubset_it=stub_set_.begin(); stubset_it!=stub_set_.end(); ++stubset_it ) {
		core::Size const stub_number( size( stubset_it->first ) );
		if ( stub_number < access_point ) {
			access_point -= stub_number;
		} else {
			break;
		}
	}
	std::multimap<core::Real,HotspotStubOP >::iterator hs_it;
	runtime_assert( access_point <= stubset_it->second.size() );
	runtime_assert( access_point > 0 );
	hs_it = stubset_it->second.begin();
	for ( core::Size i=1; i < access_point; ++i, ++hs_it ) {}
	HotspotStubOP returnop( new HotspotStub( *hs_it->second ) );
	return( returnop );
}

HotspotStubOP
HotspotStubSet::random_stub( std::string const resname )
{
	std::multimap<core::Real,HotspotStubOP > stubs = retrieve( resname );
	core::Size const subset_size( stubs.size() );
	core::Size const random_element( ( core::Size ) ( numeric::random::rg().uniform() * subset_size ) + 1 );

	runtime_assert( random_element <= subset_size );
	std::multimap<core::Real,HotspotStubOP>::iterator hs_it;
	hs_it = stubs.begin();
	for ( core::Size i=1; i < random_element; ++i, ++hs_it ) {}
	HotspotStubOP returnop( new HotspotStub( *hs_it->second ) );
	return( returnop );
}

void HotspotStubSet::create_hotspot_after_pose(core::pose::Pose & pose, std::string const & resname )
{
	core::conformation::ResidueOP residue;
	core::chemical::ResidueTypeSet const & residue_set ( pose.residue(1).residue_type_set() );

	// if we're targeting a residue, make a copy of it to get a nearby location (we'll switch identity later)
	if ( target_resnum_ && target_distance_ ) {
		residue = core::conformation::ResidueOP( new core::conformation::Residue( pose.residue(target_resnum_) ) );
		//residue->phi( -150 );
		//residue->psi( 150 );
		utility::vector1< core::Real > sheet_phipsi;
		sheet_phipsi.push_back( -150.0 );
		sheet_phipsi.push_back( 150.0 );
		sheet_phipsi.push_back( 180.0 );
		utility::vector1< core::Real > helix_phipsi; // phi/psi taken from most probable helix, not geometrically pure one
		helix_phipsi.push_back( -64.0 );
		helix_phipsi.push_back( -41.0 );
		helix_phipsi.push_back( 180.0 );

		residue->mainchain_torsions( helix_phipsi );
	} else {
		// otherwise, make a new residue from scratch at 0,0,0
		core::chemical::ResidueType const & restype( residue_set.name_map( resname ) );
		residue = core::conformation::ResidueFactory::create_residue( restype );
	}
	core::conformation::ResidueOP ala;
	core::chemical::ResidueType const & alatype( residue_set.name_map( "ALA" ) );
	ala = core::conformation::ResidueFactory::create_residue( alatype );

	// this scorefxn doesn't really matter, since we're just using it to repack the new residue and make sure the pose is scored.
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );

	for ( core::Size i = 1; i <= hotspot_length(); ++i ) {
		// add multimer and set extended chain parameters, eg - ala-Hotspot-ala-ala
		// place hotspot in middle of peptide
		if ( i == 1 ) {
			if ( hotspot_length() == 1 ) pose.append_residue_by_jump( *residue, pose.total_residue(), "", "", true );
			else pose.append_residue_by_jump( *ala, pose.total_residue(), "", "", true );
		} else if ( i == (hotspot_length()/2) + 1 ) { // middle of the peptide
			pose.append_residue_by_bond( *residue, true /*ideal geometry*/ );
		} else {
			pose.append_residue_by_bond( *ala, true );
		}
	}
	for ( core::Size i = pose.total_residue() - hotspot_length() + 1; i <= pose.total_residue(); ++i ) {
		TR.Debug << "i phi/psi/omega " << i << " " << pose.phi(i) << " " << pose.psi(i) << " " << pose.omega(i) << std::endl;
		pose.set_phi( i, -64.0 ); // phi/psi for helix
		pose.set_psi( i, -41.0 );
		pose.set_omega( i, 180 );
		TR.Debug << "after i phi/psi/omega " << i << " " << pose.phi(i) << " " << pose.psi(i) << " " << pose.omega(i) << std::endl;
	}
	// option to use make backbone invisible (defaults true), overridden by threemer
	if ( sc_only_ ) {
		for ( core::Size i = pose.total_residue() - hotspot_length() + 1; i <= pose.total_residue(); ++i ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::SHOVE_BB, i );
		}
	}

	// move hotspots slightly off current location, in case we started precisely on another residue (target_res case)
	protocols::rigid::RigidBodyTransMover trans_mover( pose, pose.num_jump() );
	//trans_mover.trans_axis( trans_mover.trans_axis() );
	trans_mover.step_size(1e-5);
	trans_mover.apply( pose );


	//(*scorefxn)(pose);
	core::pack::task::PackerTaskOP packer_task = core::pack::task::TaskFactory::create_packer_task( pose );
	for ( Size ii=1; ii <= pose.total_residue() - hotspot_length(); ++ii ) {
		packer_task->nonconst_residue_task( ii ).prevent_repacking();
	}
	//make sure our appended aa has the correct identity and has realistic coordinates
	utility::vector1<bool> restrict_to_aa(20,false);
	restrict_to_aa[ core::chemical::aa_from_name(resname) ] = true;

	for ( core::Size i = pose.total_residue()-hotspot_length()+1 ; i <= pose.total_residue(); ++i ) {
		packer_task->nonconst_residue_task( i ).restrict_absent_canonical_aas( restrict_to_aa );
	}

	//core::pack::pack_rotamers( pose, *scorefxn, packer_task );
	protocols::simple_moves::PackRotamersMover packrotamersmover( scorefxn, packer_task, 1 /*nloops*/ );
	packrotamersmover.apply( pose );
	return;
}

/// @brief set up foldtree for residue docking. If target_res is defined, takes that as connection point.
/// Assumes hotspot to dock is last res in the pose! Assumes monomeric fold_tree!
void HotspotStubSet::setup_hotspot_foldtree_( core::pose::Pose & pose ) {
	core::Size cutpoint = pose.total_residue()-hotspot_length();
	runtime_assert( cutpoint );

	core::Size jump_pos1(0);
	core::Size jump_pos2(0);

	// calculate jump points from partners
	TR.Debug << "cutpoint: " << cutpoint << std::endl;
	if ( target_resnum_ && target_distance_ ) {
		jump_pos1 = target_resnum_;
	} else jump_pos1 = core::pose::residue_center_of_mass( pose, 1, cutpoint );
	TR.Debug << "jump1: " << jump_pos1 << std::endl;
	// hotspot created stub_offset from the end of the pose
	jump_pos2 = pose.total_residue() - stub_offset();
	TR.Debug << "jump2: " << jump_pos2 << std::endl;

	runtime_assert( jump_pos1 );
	runtime_assert( jump_pos2 );

	core::kinematics::FoldTree f( pose.fold_tree() );
	f.slide_jump( pose.num_jump(), jump_pos1, jump_pos2 ); // make our last jump

	f.reorder( 1 );
	f.check_fold_tree();
	runtime_assert( f.check_fold_tree() );
	pose.fold_tree( f );
}

/// @brief utility function to find distance of a hotspot from the end of the current pose. Should be in middle or more C-terminal if hotspot_length is even
core::Size HotspotStubSet::stub_offset() {
	if ( hotspot_length() == 0 ) return 0;
	core::Size const offset = (hotspot_length()%2 == 0.0) ? (hotspot_length()/2) - 1 : (hotspot_length()/2);
	return offset;
}

/// @brief utility function to calculate interaction energy for a single residue
core::Real HotspotStubSet::get_residue_score_ ( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn, Size const seqpos) {
	protocols::simple_filters::DdgFilter const ddg( 100/*ddg_threshold*/, scorefxn, 1 /*jump__number*/, 3 /*repeats*/ );
	/*(*scorefxn)(pose);
	core::Real weighted_fa_atr = pose.energies().residue_total_energies( placed_seqpos )[core::scoring::fa_atr] * scorefxn->weights()[core::scoring::fa_atr];
	core::Real weighted_fa_rep = pose.energies().residue_total_energies( placed_seqpos )[core::scoring::fa_rep] * scorefxn->weights()[core::scoring::fa_rep];
	core::Real weighted_hbond_bb_sc = pose.energies().residue_total_energies( placed_seqpos )[core::scoring::hbond_bb_sc] * scorefxn->weights()[core::scoring::hbond_bb_sc];
	core::Real weighted_hbond_sc = pose.energies().residue_total_energies( placed_seqpos )[core::scoring::hbond_sc] * scorefxn->weights()[core::scoring::hbond_sc];
	core::Real weighted_fa_sol = pose.energies().residue_total_energies( placed_seqpos )[core::scoring::fa_sol] * scorefxn->weights()[core::scoring::fa_sol];
	core::Real weighted_contact_score = weighted_fa_atr + weighted_fa_rep + weighted_hbond_bb_sc + weighted_hbond_sc + weighted_fa_sol;
	TR.Debug << pose.energies().residue_total_energies( placed_seqpos ).weighted_string_of( scorefxn->weights() ) << std::endl;
	return weighted_contact_score;
	*/

	core::pose::Pose copy_pose( pose );

	// extract hotspot res from peptide, remove peptide from pose, re-append hotspot res, then score
	core::Size const hotspot_start = (hotspot_length() == 1) ? pose.total_residue() : pose.total_residue() - hotspot_length() + 1;
	core::Size const hotspot_stop = pose.total_residue();
	core::conformation::ResidueOP hotspot( new core::conformation::Residue( copy_pose.residue( seqpos ) ) );
	copy_pose.conformation().delete_residue_range_slow( hotspot_start, hotspot_stop );
	copy_pose.append_residue_by_jump( *hotspot, copy_pose.total_residue(), "","", true );
	(*scorefxn)(copy_pose); // necessary ro get neighbors for the Interface class, which is needed by DdgFilter
	core::Real const score = ddg.compute( copy_pose );

	return score;
}

/// @brief utility function to add hotspot bbcst's to a pose
void
HotspotStubSet::add_hotspot_constraints_to_pose(
	core::pose::Pose & pose,
	core::Size const partner,
	HotspotStubSetOP hotspot_stub_set,
	core::Real const & CB_force_constant,
	core::Real const & worst_allowed_stub_bonus,
	bool const apply_self_energies,
	core::Real const & bump_cutoff,
	bool const apply_ambiguous_constraints // = false
) {

	// Assign a fixed residue (for the stub constraints)
	core::Size fixed_res(1);
	if ( partner == 1 ) fixed_res = pose.total_residue();
	core::id::AtomID fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );

	core::pack::task::PackerTaskOP packer_task = prepare_hashing_packer_task_( pose, partner );

	add_hotspot_constraints_to_pose(
		pose,
		fixed_atom_id,
		packer_task,
		hotspot_stub_set,
		CB_force_constant,
		worst_allowed_stub_bonus,
		apply_self_energies,
		bump_cutoff,
		apply_ambiguous_constraints );
}


/// @brief utility function to add hotspot bbcst's to a pose
void
HotspotStubSet::add_hotspot_constraints_to_wholepose(
	core::pose::Pose & pose,
	core::Size const partner,
	HotspotStubSetOP hotspot_stub_set,
	core::Real const & CB_force_constant,
	core::Real const & worst_allowed_stub_bonus,
	bool const apply_self_energies,
	core::Real const & bump_cutoff,
	bool const apply_ambiguous_constraints // = false
) {

	// Assign a fixed residue (for the stub constraints)
	core::Size fixed_res(1);
	if ( partner == 1 ) fixed_res = pose.total_residue();
	core::id::AtomID fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );

	core::pack::task::PackerTaskOP packer_task = prepare_hashing_packer_task_( pose, partner );

	add_hotspot_constraints_to_wholepose(
		pose,
		fixed_atom_id,
		packer_task,
		hotspot_stub_set,
		CB_force_constant,
		worst_allowed_stub_bonus,
		apply_self_energies,
		bump_cutoff,
		apply_ambiguous_constraints );
}

void
HotspotStubSet::add_hotspot_constraints_to_pose(
	core::pose::Pose & pose,
	core::id::AtomID const & fixed_atom,
	core::pack::task::PackerTaskCOP const packer_task,
	HotspotStubSetOP hotspot_stub_set,
	core::Real const & CB_force_constant,
	core::Real const & worst_allowed_stub_bonus,
	bool const apply_self_energies,
	core::Real const & bump_cutoff,
	bool const apply_ambiguous_constraints // = false
) {

	// Take in a scaffold pose (with PackerTask, for its DesignMap), and a set of stubs.
	// Each repacked residue will get one "AmbiguousConstraint".
	// This AmbiguousConstraint will contain a series of BackboneStubConstraints (one for each valid stub)

	runtime_assert( CB_force_constant > -1E-6 ); // these can't be negative
	runtime_assert( worst_allowed_stub_bonus < 1E-6 ); // these can't be positive
	runtime_assert( bump_cutoff > -1E-6 ); // these can't be negative

	// setup for evaluating the stub in the context of the scaffold
	core::pose::Pose unbound_pose = pose;
	generate_unbound_pose_(unbound_pose);

	// scorefxn for the bump check
	core::scoring::ScoreFunctionOP bump_scorefxn( new core::scoring::ScoreFunction );
	bump_scorefxn->reset();
	bump_scorefxn->set_weight( core::scoring::fa_rep, 1.0 );

	// set scorefunction.
	core::scoring::ScoreFunctionOP noenvhbond_scorefxn;
	//if( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
	// noenvhbond_scorefxn = core::scoring::get_score_function();
	//}
	//else {
	// just use sc12 with no env hbonding. Legacy reasons for this scorefxn, but it's only used for building neighbor graphs, packing ala poses, and bump checking.
	noenvhbond_scorefxn = core::scoring::get_score_function_legacy( "score13" );
	noenvhbond_scorefxn->set_weight( core::scoring::fa_dun, 0.1 );
	noenvhbond_scorefxn->set_weight( core::scoring::envsmooth, 0 );
	//}
	core::scoring::methods::EnergyMethodOptions options( noenvhbond_scorefxn->energy_method_options() );
	options.hbond_options().use_hb_env_dep( basic::options::option[ basic::options::OptionKeys::hotspot::envhb]() );
	noenvhbond_scorefxn->set_energy_method_options( options );

	core::scoring::ScoreFunctionOP full_scorefxn( noenvhbond_scorefxn );

	// core::scoring::ScoreFunctionOP full_scorefxn( core::scoring::get_score_function() );

	// score the pose, to update the tenA_neighbor_graph and setup Hbond stuff
	(*full_scorefxn)(unbound_pose);
	core::scoring::TenANeighborGraph const unbound_neighbor_graph = unbound_pose.energies().tenA_neighbor_graph();

	// make an alanine pose to use for checking self energies (since we'll almost always be evaluating stubs in the absence of sidechains)
	core::pose::Pose ala_pose = unbound_pose;
	utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, false );
	allowed_aas[ core::chemical::aa_ala ] = true;
	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().or_include_current( true );
	for ( core::Size i=1; i<=ala_pose.total_residue(); ++i ) {
		if ( packer_task->pack_residue(i) ) {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
		} else {
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		core::pack::task::parse_resfile(pose, *task);
	}
	core::pack::pack_rotamers( ala_pose, *full_scorefxn, task);
	(*full_scorefxn)( ala_pose ); // to ensure that 10Aneighborgraph_state==GOOD
	core::scoring::TenANeighborGraph const ala_neighbor_graph = ala_pose.energies().tenA_neighbor_graph();


	// *****associate the stub set with the unbound pose
	// *****hotspot_stub_set->pair_with_scaffold( pose, partner );

	protocols::filters::FilterCOP true_filter( protocols::filters::FilterOP( new protocols::filters::TrueFilter ) );
	for ( core::Size resnum=1; resnum <= pose.total_residue(); ++resnum ) {
		if ( packer_task->pack_residue(resnum) ) {
			hotspot_stub_set->pair_with_scaffold( pose, pose.chain( resnum ), true_filter );
			break;
		}
	}

	TR << "Making hotspot constraints..." << std::endl;
	//TR << pose.total_residue() << " residues" << std::endl;
	//Size scaffold_seqpos(0);  // unused ~Labonte
	for ( core::Size resnum=1; resnum <= pose.total_residue(); ++resnum ) {

		// Check that this position is allowed to be used for stub constraints
		if ( ! packer_task->pack_residue(resnum) ) continue;

		// sets the index used by the hotspot for its associated scaffold
		//scaffold_seqpos = resnum - pose.conformation().chain_begin( pose.chain( resnum ) );  // unused ~Labonte

		// Start the vector which will become a single AmbiguousConstraint, if apply_ambiguous_constraints is true
		utility::vector1< core::scoring::constraints::ConstraintCOP > ambig_csts;
		// Loop over all allowed AAs at this position
		std::list< core::chemical::ResidueTypeCOP > allowed_aas = packer_task->residue_task( resnum ).allowed_residue_types();
		for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator restype = allowed_aas.begin();
				restype != allowed_aas.end(); ++restype ) {

			// Loop over all stubs with this restype
			Hotspots res_stub_set( hotspot_stub_set->retrieve( (*restype )->name3() ) );
			for ( std::multimap<core::Real,HotspotStubOP >::iterator hs_stub = res_stub_set.begin();
					hs_stub != res_stub_set.end(); ++hs_stub ) {

				// prevent Gly/Pro constraints
				if ( (hs_stub->second->residue()->aa() == core::chemical::aa_gly) || (hs_stub->second->residue()->aa() == core::chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] ) ) {
					TR << "ERROR - Gly/Pro stubs cannot be used for constraints." << std::endl;
					continue;
				}

				// prevent Gly/Pro constraints
				if ( (pose.residue(resnum).aa() == core::chemical::aa_gly) || (pose.residue(resnum).aa() == core::chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline]) ) {
					TR.Debug << "ERROR - Position " << resnum << " is currently Gly/Pro and cannot be used for stub constraints." << std::endl;
					continue;
				}

				core::Real stub_bonus_value = hs_stub->second->bonus_value();
				//TR << "stub: " << hs_stub->second->residue()->name3() << " pose" << pose.residue(resnum).name3() << " " << resnum << " StubEnergy=" << stub_bonus_value << std::endl;
				// Evaluate how this stub fits on the scaffold
				// orient the stub onto the pose
				core::conformation::Residue placed_stub_res = *(hs_stub->second->residue());
				placed_stub_res.orient_onto_residue( unbound_pose.residue(resnum) );
				// set phi,psi of the placed stub res to match the pose
				placed_stub_res.mainchain_torsions()[1] = unbound_pose.phi(resnum);
				placed_stub_res.mainchain_torsions()[2] = unbound_pose.psi(resnum);

				core::Real bump_energy = evaluate_stub_bumps_( placed_stub_res, unbound_pose, resnum, unbound_neighbor_graph, bump_scorefxn, bump_cutoff );
				if ( bump_energy > bump_cutoff ) {
					// force cutoff to fail
					stub_bonus_value = worst_allowed_stub_bonus + 0.1;
					hs_stub->second->set_scaffold_status( resnum, protocols::hotspot_hashing::reject );
					//TR << " FailBumpEnergy=" << bump_energy;
				} else if ( apply_self_energies ) {
					//TR << " SuccBumpEnergy=" << bump_energy;
					// note: if the pose can move, these energies won't be totally accurate (but they'll still be good estimates if the movement is small)
					stub_bonus_value += evaluate_stub_self_energy_( placed_stub_res, unbound_pose, resnum, unbound_neighbor_graph, full_scorefxn );
					TR.Debug << resnum << " InitBonus " << hs_stub->second->bonus_value() << "    SelfEBonus " << stub_bonus_value << std::endl;
				}
				if ( stub_bonus_value < worst_allowed_stub_bonus ) {
					hs_stub->second->set_scaffold_status( resnum, protocols::hotspot_hashing::accept );
					//TR << " SuccSelfEnergy=" << stub_bonus_value << std::endl;
					// ****** accept the pairing -- do we really want this? better to just reject, since bb fit doesn't necessarily mean good pair
					// ****** hs_stub->scaffold_status( resnum, accept );

					// Build a BackboneStubConstraint from this stub
					if ( apply_ambiguous_constraints ) {
						// Push it onto ambig_csts for this residue
						ambig_csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::BackboneStubConstraint( pose, resnum, fixed_atom, *(hs_stub->second->residue()), stub_bonus_value, CB_force_constant ) ) );
					} else {
						// Apply it directly
						constraints_.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::BackboneStubConstraint( pose, resnum, fixed_atom, *(hs_stub->second->residue()), stub_bonus_value, CB_force_constant ) ) );
					}
				} else hs_stub->second->set_scaffold_status( resnum, protocols::hotspot_hashing::reject );
				//else TR << " FailSelfEnergy=" << stub_bonus_value << std::endl;
				// ****** reject the pairing
				// ******else hs_stub->scaffold_status( resnum, reject );
			}
		}

		// Finally, add the constraint corresponding to this resnum to the main set
		if ( ( apply_ambiguous_constraints ) && ( ambig_csts.size() > 0 ) ) {
			constraints_.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AmbiguousConstraint(ambig_csts) ) );
		}
	}

	constraints_ = pose.add_constraints( constraints_ );
	TR << "Applied " << constraints_.size() << " hotspot constraints to the pose." << std::endl;
	return;
}

void
HotspotStubSet::add_hotspot_constraints_to_wholepose(
	core::pose::Pose & pose,
	core::id::AtomID const & fixed_atom,
	core::pack::task::PackerTaskCOP const packer_task,
	HotspotStubSetOP hotspot_stub_set,
	core::Real const & CB_force_constant,
	core::Real const & worst_allowed_stub_bonus,
	bool const apply_self_energies,
	core::Real const & bump_cutoff,
	bool const apply_ambiguous_constraints // = false
) {

	// Take in a scaffold pose (with PackerTask, for its DesignMap), and a set of stubs.
	// Each repacked residue will get one "AmbiguousConstraint".
	// This AmbiguousConstraint will contain a series of BackboneStubLinearConstraints (one for each valid stub)

	runtime_assert( CB_force_constant > -1E-6 ); // these can't be negative
	runtime_assert( worst_allowed_stub_bonus < 1E-6 ); // these can't be positive
	runtime_assert( bump_cutoff > -1E-6 ); // these can't be negative

	// setup for evaluating the stub in the context of the scaffold
	core::pose::Pose unbound_pose = pose;
	generate_unbound_pose_(unbound_pose);

	// scorefxn for the bump check
	core::scoring::ScoreFunctionOP bump_scorefxn( new core::scoring::ScoreFunction );
	bump_scorefxn->reset();
	bump_scorefxn->set_weight( core::scoring::fa_rep, 1.0 );

	// set scorefunction.
	core::scoring::ScoreFunctionOP noenvhbond_scorefxn;
	//if( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
	// noenvhbond_scorefxn = core::scoring::get_score_function();
	//}
	//else {
	// just use sc12 with no env hbonding. Legacy reasons for this scorefxn, but it's only used for building neighbor graphs, packing ala poses, and bump checking.
	noenvhbond_scorefxn = core::scoring::get_score_function_legacy( "score13" );
	noenvhbond_scorefxn->set_weight( core::scoring::fa_dun, 0.1 );
	noenvhbond_scorefxn->set_weight( core::scoring::envsmooth, 0 );
	//}
	core::scoring::methods::EnergyMethodOptions options( noenvhbond_scorefxn->energy_method_options() );
	options.hbond_options().use_hb_env_dep( basic::options::option[ basic::options::OptionKeys::hotspot::envhb]() );
	noenvhbond_scorefxn->set_energy_method_options( options );

	core::scoring::ScoreFunctionOP full_scorefxn( noenvhbond_scorefxn );

	// core::scoring::ScoreFunctionOP full_scorefxn( core::scoring::get_score_function() );

	// score the pose, to update the tenA_neighbor_graph and setup Hbond stuff
	(*full_scorefxn)(unbound_pose);
	core::scoring::TenANeighborGraph const unbound_neighbor_graph = unbound_pose.energies().tenA_neighbor_graph();

	// make an alanine pose to use for checking self energies (since we'll almost always be evaluating stubs in the absence of sidechains)
	core::pose::Pose ala_pose = unbound_pose;
	utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, false );
	allowed_aas[ core::chemical::aa_ala ] = true;
	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().or_include_current( true );
	for ( core::Size i=1; i<=ala_pose.total_residue(); ++i ) {
		if ( packer_task->pack_residue(i) ) {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
		} else {
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		core::pack::task::parse_resfile(pose, *task);
	}
	core::pack::pack_rotamers( ala_pose, *full_scorefxn, task);
	(*full_scorefxn)( ala_pose ); // to ensure that 10Aneighborgraph_state==GOOD
	core::scoring::TenANeighborGraph const ala_neighbor_graph = ala_pose.energies().tenA_neighbor_graph();


	// *****associate the stub set with the unbound pose
	// *****hotspot_stub_set->pair_with_scaffold( pose, partner );

	protocols::filters::FilterCOP true_filter( protocols::filters::FilterOP( new protocols::filters::TrueFilter ) );
	for ( core::Size resnum=1; resnum <= pose.total_residue(); ++resnum ) {
		if ( packer_task->pack_residue(resnum) ) {
			hotspot_stub_set->pair_with_scaffold( pose, pose.chain( resnum ), true_filter );
			break;
		}
	}

	TR << "Making hotspot constraints..." << std::endl;
	//TR<< pose.total_residue() << " residues" << std::endl;
	utility::vector1< core::scoring::constraints::ConstraintCOP > ambig_csts;
	//Size scaffold_seqpos(0);
	for ( core::Size resnum=1; resnum <= pose.total_residue(); ++resnum ) {

		// Check that this position is allowed to be used for stub constraints
		if ( ! packer_task->pack_residue(resnum) ) continue;

		// sets the index used by the hotspot for its associated scaffold
		//scaffold_seqpos = resnum - pose.conformation().chain_begin( pose.chain( resnum ) );

		// Start the vector which will become a single AmbiguousConstraint, if apply_ambiguous_constraints is true
		// Loop over all allowed AAs at this position
		std::list< core::chemical::ResidueTypeCOP > allowed_aas = packer_task->residue_task( resnum ).allowed_residue_types();
		for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator restype = allowed_aas.begin();
				restype != allowed_aas.end(); ++restype ) {

			// Loop over all stubs with this restype
			Hotspots res_stub_set( hotspot_stub_set->retrieve( (*restype )->name3() ) );
			for ( std::multimap<core::Real,HotspotStubOP >::iterator hs_stub = res_stub_set.begin();
					hs_stub != res_stub_set.end(); ++hs_stub ) {

				// prevent Gly/Pro constraints
				if ( (hs_stub->second->residue()->aa() == core::chemical::aa_gly) || (hs_stub->second->residue()->aa() == core::chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] ) ) {
					TR << "ERROR - Gly/Pro stubs cannot be used for constraints." << std::endl;
					continue;
				}

				// prevent Gly/Pro constraints
				if ( (pose.residue(resnum).aa() == core::chemical::aa_gly) || (pose.residue(resnum).aa() == core::chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline]) ) {
					TR.Debug << "ERROR - Position " << resnum << " is currently Gly/Pro and cannot be used for stub constraints." << std::endl;
					continue;
				}

				core::Real stub_bonus_value = hs_stub->second->bonus_value();
				//TR << "stub: " << hs_stub->second->residue()->name3() << " pose" << pose.residue(resnum).name3() << " " << resnum << " StubEnergy=" << stub_bonus_value;
				// Evaluate how this stub fits on the scaffold
				// orient the stub onto the pose
				core::conformation::Residue placed_stub_res = *(hs_stub->second->residue());
				placed_stub_res.orient_onto_residue( unbound_pose.residue(resnum) );
				// set phi,psi of the placed stub res to match the pose
				placed_stub_res.mainchain_torsions()[1] = unbound_pose.phi(resnum);
				placed_stub_res.mainchain_torsions()[2] = unbound_pose.psi(resnum);

				core::Real bump_energy = evaluate_stub_bumps_( placed_stub_res, unbound_pose, resnum, unbound_neighbor_graph, bump_scorefxn, bump_cutoff );
				if ( bump_energy > bump_cutoff ) {
					// force cutoff to fail
					stub_bonus_value = worst_allowed_stub_bonus + 0.1;
					hs_stub->second->set_scaffold_status( resnum, protocols::hotspot_hashing::reject );
					//TR << " FailBumpEnergy=" << bump_energy;
				} else if ( apply_self_energies ) {
					//TR << " SuccBumpEnergy=" << bump_energy;
					// note: if the pose can move, these energies won't be totally accurate (but they'll still be good estimates if the movement is small)
					stub_bonus_value += evaluate_stub_self_energy_( placed_stub_res, unbound_pose, resnum, unbound_neighbor_graph, full_scorefxn );
					TR.Debug << resnum << " InitBonus " << hs_stub->second->bonus_value() << "    SelfEBonus " << stub_bonus_value << std::endl;
				}
				if ( stub_bonus_value < worst_allowed_stub_bonus ) {
					hs_stub->second->set_scaffold_status( resnum, protocols::hotspot_hashing::accept );
					//TR << " SuccSelfEnergy=" << stub_bonus_value << std::endl;
					// ****** accept the pairing -- do we really want this? better to just reject, since bb fit doesn't necessarily mean good pair
					// ****** hs_stub->scaffold_status( resnum, accept );

					// Build a BackboneStubLinearConstraint from this stub
					if ( apply_ambiguous_constraints ) {
						// Push it onto ambig_csts for this residue
						ambig_csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::BackboneStubLinearConstraint( pose, resnum, fixed_atom, *(hs_stub->second->residue()), stub_bonus_value, CB_force_constant ) ) );
					} else {
						// Apply it directly
						constraints_.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::BackboneStubLinearConstraint( pose, resnum, fixed_atom, *(hs_stub->second->residue()), stub_bonus_value, CB_force_constant ) ) );
					}
				} else hs_stub->second->set_scaffold_status( resnum, protocols::hotspot_hashing::reject );
				//else TR << " FailSelfEnergy=" << stub_bonus_value << std::endl;
				// ****** reject the pairing
				// ******else hs_stub->scaffold_status( resnum, reject );
			} //end all stubs
			//TR << "Add ambiguous cst of size for all stubs for a residue: " << ambig_csts.size() << std::endl;
		}//end all residues
	}

	// Finally, add the constraint corresponding of all resnum to the main set
	if ( ( apply_ambiguous_constraints ) && ( ambig_csts.size() > 0 ) ) {
		//TR << "Add ambiguous cst of size: " << ambig_csts.size() << std::endl;
		//TR << "Add ambiguous cst of size for all stubs and all residues: " << ambig_csts.size() << std::endl;
		constraints_.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AmbiguousConstraint(ambig_csts) ) );
	}

	constraints_ = pose.add_constraints( constraints_ );
	TR << "Applied " << constraints_.size() << " ambiguous hotspot constraints to whole pose." << std::endl;
	return;
}
bool
HotspotStubSet::remove_all_hotspot_constraints( core::pose::Pose & pose ) const
{
	bool return_val( pose.remove_constraints( constraints_ ) );
	constraints_.empty();
	return( return_val ); // returns the state of constraints removal
}

void HotspotStubSet::set_chain_to_design( core::Size const chain_to_design ) {
	chain_to_design_ = chain_to_design;
}

core::Real
HotspotStubSet::evaluate_stub_bumps_(
	core::conformation::Residue const & placed_stub_residue,
	core::pose::Pose const & unbound_pose,
	core::Size const resnum,
	core::scoring::TenANeighborGraph const & unbound_neighbor_graph,
	core::scoring::ScoreFunctionCOP const & bump_scorefxn,
	core::Real const & max_bump_energy
) {

	core::Real bump_energy = 0.;
	core::scoring::EnergyMap emap;
	for ( core::graph::Graph::EdgeListConstIter
			iter = unbound_neighbor_graph.get_node( resnum )->const_edge_list_begin(),
			iter_end = unbound_neighbor_graph.get_node( resnum )->const_edge_list_end();
			iter != iter_end; ++iter ) {
		Size const neighbor_resnum( (*iter)->get_other_ind( resnum ) );
		emap.clear();
		bump_scorefxn->bump_check_backbone( placed_stub_residue, unbound_pose.residue(neighbor_resnum), unbound_pose, emap );
		bump_energy += emap[ core::scoring::fa_rep ];
		if ( bump_energy > max_bump_energy ) break;
	}

	return bump_energy;
}


core::Real
HotspotStubSet::evaluate_stub_self_energy_(
	core::conformation::Residue const & placed_stub_residue,
	core::pose::Pose const & unbound_pose,
	core::Size const resnum,
	core::scoring::TenANeighborGraph const & unbound_neighbor_graph,
	core::scoring::ScoreFunctionCOP const & full_scorefxn
) {

	core::Real self_energy = 0.;
	TR.Debug << resnum << " InitStubSelfE " << self_energy << " ";
	core::scoring::EnergyMap weights = full_scorefxn->weights();

	// probability of the stub AA given the phi/psi of the (starting) pose
	core::scoring::EnergyMap emap;
	emap.clear();
	full_scorefxn->eval_ci_1b( placed_stub_residue, unbound_pose, emap );
	self_energy += weights.dot( emap );
	TR.Debug << self_energy << " ";

	for ( core::graph::Graph::EdgeListConstIter
			iter = unbound_neighbor_graph.get_node( resnum )->const_edge_list_begin(),
			iter_end = unbound_neighbor_graph.get_node( resnum )->const_edge_list_end();
			iter != iter_end; ++iter ) {
		Size const neighbor_resnum( (*iter)->get_other_ind( resnum ) );

		// how well the sidechain of the stub fits into the pose
		core::scoring::EnergyMap emap_2b;
		emap_2b.clear();
		// backbone of neighbor, sidechain of stub
		full_scorefxn->eval_ci_2b_bb_sc( unbound_pose.residue(neighbor_resnum), placed_stub_residue, unbound_pose, emap_2b );
		full_scorefxn->eval_cd_2b_bb_sc( unbound_pose.residue(neighbor_resnum), placed_stub_residue, unbound_pose, emap_2b );
		// sidechain of neighbor, sidechain of stub
		full_scorefxn->eval_ci_2b_sc_sc( unbound_pose.residue(neighbor_resnum), placed_stub_residue, unbound_pose, emap_2b );
		full_scorefxn->eval_cd_2b_sc_sc( unbound_pose.residue(neighbor_resnum), placed_stub_residue, unbound_pose, emap_2b );
		self_energy += weights.dot( emap_2b );
		TR.Debug << self_energy << " ";
	}
	TR.Debug << std::endl;
	return self_energy;
}

core::pack::task::PackerTaskOP
HotspotStubSet::prepare_hashing_packer_task_(
	core::pose::Pose const & pose,
	core::Size const chain_to_redesign // default = 2
) {

	core::pack::task::PackerTaskOP const hotspot_hash_packer_taskOP =
		core::pack::task::TaskFactory::create_packer_task( pose );

	//in case there is a resfile, information in this resfile can set the packer task, but can be overridden below.
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		core::pack::task::parse_resfile(pose, *hotspot_hash_packer_taskOP);
	}

	// Setup the unbound pose, from which to compute residue sasa
	// (move the chains apart first, in case we're starting from a good binding mode)
	core::pose::Pose unbound_pose = pose;
	generate_unbound_pose_(unbound_pose);

	// Compute residue sasa
	core::Real const sasa_thres(2.0);
	core::Real const probe_radius(1.4);
	core::id::AtomID_Map< core::Real > atom_sasa;
	utility::vector1< core::Real > residue_sasa;
	core::scoring::calc_per_atom_sasa( unbound_pose, atom_sasa, residue_sasa, probe_radius);

	// Specify which chain is the designed one
	utility::vector1<bool> allow_stub_constraint( pose.total_residue(), false );
	for ( int ii = pose.conformation().chain_begin( chain_to_redesign ),
			lastres = pose.conformation().chain_end( chain_to_redesign ); ii <= lastres; ++ii ) {
		// don't allow residues that are currently Gly or Pro to land on stubs
		if ( pose.residue(ii).aa() != core::chemical::aa_gly && (pose.residue(ii).aa() != core::chemical::aa_pro || basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] ) ) {
			if ( residue_sasa.at(ii) > sasa_thres ) {
				allow_stub_constraint.at(ii) = true;
			}
		}
	}
	hotspot_hash_packer_taskOP->restrict_to_residues( allow_stub_constraint ); //overriding the resfile

	return hotspot_hash_packer_taskOP;
}

void HotspotStubSet::generate_unbound_pose_( core::pose::Pose & pose ) {
	core::Real const unbound_dist = 40.;
	Size const rb_jump = 1; // use the first jump as the one between partners
	protocols::rigid::RigidBodyTransMover trans_mover( pose, rb_jump );
	trans_mover.trans_axis( trans_mover.trans_axis() );
	trans_mover.step_size(unbound_dist);
	trans_mover.apply( pose );
	return;
}


/// @details Searches for all ambiguous constraints that have a backbone stub active constraint.
/// and removes them from the pose. Returns a vector of the removed constraints.
core::scoring::constraints::ConstraintCOPs
remove_hotspot_constraints_from_pose( core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;

	ConstraintCOPs original_csts = pose.constraint_set()->get_all_constraints() ;
	ConstraintCOPs bb_csts;
	for ( ConstraintCOPs::const_iterator it = original_csts.begin(), end = original_csts.end(); it != end; ++it ) {
		ConstraintCOP cst( *it );
		if ( cst->type() == "AmbiguousConstraint" ) {
			AmbiguousConstraintCOP ambiguous_cst = AmbiguousConstraintCOP( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::AmbiguousConstraint const > ( cst ) ); //downcast to derived ambiguous constraint
			if ( ambiguous_cst ) { // safety check for downcasting
				if ( ambiguous_cst->active_constraint()->type() == "BackboneStubLinear" || ambiguous_cst->active_constraint()->type() == "BackboneStub" ) {
					bb_csts.push_back( ambiguous_cst ); // add the entire ambiguous cst, since it contained a bbcst
				}
			}
		} else if ( cst->type() == "BackboneStubLinear" || cst->type() == "BackboneStub" ) {
			//tricky; some hotspot stub constraints are added directly with no ambiguity
			bb_csts.push_back( cst );
		}
	}
	pose.remove_constraints( bb_csts ); // remove all the ambigcsts that contain a bbcst
	return( bb_csts );
}

/// @details Iterates over all non-hydrogen sidechain atoms of two residues and returns their rmsd without superposition.
core::Real
residue_sc_rmsd_no_super( core::conformation::ResidueCOP res1, core::conformation::ResidueCOP res2, bool const fxnal_group_only /*false*/ )
{
	return core::scoring::residue_sc_rmsd_no_super( res1, res2, fxnal_group_only );
}

/// @brief utility function to make sure stub's Cbeta is not pointed away from the target.
/// Returns true if stub is pointed generally towards target's center of mass. Returns false otherwise.
/// algorithm: get vector between CA and target center of mass, then CA-CB vector. Check alignment between the two.
/// doesn't take into account hotspots at the end of the pose so might be a little off in CoM. But probably pretty OK
core::Real stub_tgt_angle( core::pose::Pose const & pose, core::conformation::ResidueCOP stub, core::Size const target_res ) {
	runtime_assert( target_res <= pose.total_residue() );
	core::Size angle_res = target_res;
	if ( angle_res == 0 ) {
		angle_res = core::pose::residue_center_of_mass( pose, 1, int(pose.total_residue() ) );
	}

	core::Size const CA_stub_index = stub->atom_index("CA");
	core::Size const CB_stub_index = stub->atom_index("CB");
	core::Size const CA_tgt_index = pose.residue( angle_res ).atom_index("CA");

	core::Vector const CA_tgt = stub->xyz( CA_stub_index ) - pose.residue( angle_res ).xyz( CA_tgt_index );
	core::Vector const CA_CB = stub->xyz( CA_stub_index ) - stub->xyz( CB_stub_index );
	TR.Debug << "CA-target " << CA_tgt.length() << "\n";
	TR.Debug << "CA-CB " << CA_CB.length() << "\n"; // should always be 1.533 A
	core::Real const angle = numeric::conversions::degrees(angle_of( CA_tgt, CA_CB ));
	TR.Debug << "target-CA-CB angle " << angle;

	return angle;
}

/*
core::scoring::constraints::AngleConstraintOP make_angle_cst( core::pose::Pose const & pose, core::conformation::ResidueCOP stub ) {
core::Size const CoM = core::pose::residue_center_of_mass( pose, 1, int(pose.total_residue()) );
core::Size const CA_stub_index = stub->atom_index("CA");
core::Size const CB_stub_index = stub->atom_index("CB");
core::Size const CA_CoM_index = pose.residue( CoM ).atom_index("CA");

// note: PeriodicFunc has functional form y = ( k * cos(n * (x - x0) ) ) + C.
FuncOP cos_func = new PeriodicFunc(0., 1., 1., 0.);
anglecst = new AngleConstraint( cos_func );
pose.add_constraint( anglecst );
}
*/


} // hotspot_hashing
} // protocols
