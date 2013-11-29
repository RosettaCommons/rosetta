// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/AddCavitiesMover.cc
///
/// @brief
/// @author

#include <protocols/simple_moves/AddCavitiesMover.hh>

#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/compute_sasa.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/chemical/residue_io.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheablePoseRawPtr.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

	using namespace core;
	using namespace core::scoring::packstat;
	using core::id::AtomID;


	AddCavitiesMover::AddCavitiesMover(
		core::Size max_cav ,
		core::Real min_size,
		core::Size min_nb  ,
		core::Real min_sep
	) : max_cav_(max_cav), min_size_(min_size), min_nb_(min_nb), min_sep_(min_sep) {}

	core::id::AtomID
	AddCavitiesMover::get_closest_heavy_atom( Pose & pose, numeric::xyzVector<Real> xyz ) {
		id::AtomID closest;
		Real closest_dist = 9e9;
		for( size_t ir = 1; ir <= pose.total_residue(); ++ir ) {
			for( size_t ia = 1; ia <= pose.residue(ir).natoms(); ++ia ) {
				id::AtomID id( ia, ir );
				if( 22 <= pose.residue(id.rsd()).atom(id.atomno()).type() ) continue;
				if( pose.xyz(id).distance_squared( xyz ) < closest_dist ) {
					closest_dist = pose.xyz(id).distance_squared( xyz );
					closest = id;
				}
			}
		}
		return closest;
	}

	core::conformation::ResidueOP
	AddCavitiesMover::get_suck_res() {
		using namespace core;
		using namespace chemical;
		using namespace conformation;
		ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		return ResidueFactory::create_residue( residue_set->name_map("SUCK") );
	}

	CavBalls
	AddCavitiesMover::get_cavities( Pose const & pose, Real nbdis, int nbcount, Real minsep ) {
		Spheres spheres = core::scoring::packstat::pose_to_pack_data(pose).spheres;

		SasaOptions opts;
		opts.prune_max_iters = 5;
		opts.prune_max_delta = 0;
		opts.num_cav_ball_layers = 2;
		opts.frac_cav_ball_required_exposed = 0.0;
		opts.area_cav_ball_required_exposed = 0.0;
		opts.surrounding_sasa_smoothing_window = 1;
		for( PackstatReal pr = 3.0; pr >= 0.9; pr -= 0.1 )
			opts.probe_radii.push_back(pr);
		opts.prune_cavity_burial_probe_radii.push_back( 1.6 );

		CavBalls cball = compute_sasa( spheres, opts )->cavballs;
		CavBalls cbph  = prune_hidden_cavity_balls( cball, opts );
		CavBalls cbpr  = prune_cavity_balls( spheres, cbph, opts );
		compute_cav_ball_neighbor_count( spheres, cbpr, nbdis );
		CavBalls cbbur;
		for( CavBallIter i = cbpr.begin(); i != cbpr.end(); ++i ) if( i->anb > nbcount ) cbbur.push_back( *i );
		CavBalls selcb = select_cav_balls(cbbur,minsep);

		CavBalls cb;

		return selcb;
	}

	void
	AddCavitiesMover::clear_suckers( Pose & pose ) {
		//using namespace core::pose::datacache::CacheableDataType;
		using namespace basic::datacache;
		using namespace core::pose::datacache;
		if( !pose.data().has(core::pose::datacache::CacheableDataType::POSE_BEFORE_CAVITIES_ADDED) ) {
			return;
		}
		CacheableDataOP cd = pose.data().get_ptr( core::pose::datacache::CacheableDataType::POSE_BEFORE_CAVITIES_ADDED );
		core::pose::PoseOP cache_pose = dynamic_cast< core::pose::datacache::CacheablePoseRawPtr*>(cd())->pose();
		pose.data().set( core::pose::datacache::CacheableDataType::POSE_BEFORE_CAVITIES_ADDED, NULL );
		Pose orig_pose = *cache_pose;
		orig_pose.copy_segment( orig_pose.total_residue(), pose, 1, 1 );
		pose = orig_pose;
	}

	void
	AddCavitiesMover::add_suckers( Pose & pose ) {
		using namespace core::scoring;
		using namespace constraints;
		//using namespace core::pose::datacache::CacheableDataType;
		using namespace basic::datacache;
		using namespace core::pose::datacache;

		clear_suckers(pose);

		using namespace basic;
		pose.data().set( core::pose::datacache::CacheableDataType::POSE_BEFORE_CAVITIES_ADDED, new CacheablePoseRawPtr( new Pose(pose)) );

		CavBalls cbs = get_cavities(pose, 10.0, min_nb_, 3.0 );
		int Ncb = max_cav_;

		// add VRT res for coord constraints
		pose.append_residue_by_jump
				( *conformation::ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ),
					pose.total_residue()/2 );
		int virt_resno = pose.total_residue();
		FuncOP func( new HarmonicFunc( 0.0, 1.0 ) );

		// std::cerr << "add sucker atoms" << std::endl;
		int count = 0;
		for( int i = 1; i <= std::min(Ncb,(int)cbs.size()); ++i ) {
			if( cbs[i].radius() < min_size_ ) {
				Ncb = i-1;
				break;
			}
			// std::cerr << "adding cb" << cbs[i].str() << std::endl;
			conformation::ResidueOP sucker = get_suck_res();
			sucker->set_xyz( 1, cbs[i].xyz() );
			AtomID closest = get_closest_heavy_atom( pose, cbs[i].xyz() );
			// std::cerr << "closest " << closest << " score before " << (*sf)(pose) <<	std::endl;
			pose.append_residue_by_jump( *sucker, closest.rsd() );
			int suck_resno = pose.total_residue();

			pose.add_constraint( new CoordinateConstraint( AtomID(1,suck_resno),
																												 AtomID(1,virt_resno),
																												 sucker->xyz(1),
																												 func ) );
		  ++count;
		}
		//std::cerr << "added " << count << " suckers" << std::endl;

	}
	// strip off the suckers and virt res

	void
	AddCavitiesMover::apply( core::pose::Pose & pose ) {
		clear_suckers( pose );
		if( 0 < max_cav_ ) {
			add_suckers( pose );
		}
	}

	std::string
	AddCavitiesMover::get_name() const {
		return "AddCavitiesMover";
	}



} // end namespace simple_moves
} // end namespace protocols
