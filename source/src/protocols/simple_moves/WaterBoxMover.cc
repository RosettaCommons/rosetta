// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file WaterBoxMover.cc

// Unit headers
#include <protocols/simple_moves/WaterBoxMover.hh>
#include <protocols/simple_moves/WaterBoxMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator_Private.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/chains_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/annealer/FixbbPwatSimAnnealer.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/xml_util.hh>

#include <core/scoring/lkball/LK_BallInfo.hh>

#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/packer_neighbors.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/io/izstream.hh>
#include <utility/graph/Graph.hh>

#include <utility/thread/threadsafe_creation.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>


namespace protocols {
namespace simple_moves {

static basic::Tracer TR("protocols.moves.WaterBoxMover");


WaterRotsDB WaterBoxMover::water_rots_db_;

#ifdef MULTI_THREADED
utility::thread::ReadWriteMutex WaterBoxMover::db_mutex_;
#endif

// adds a single water to the water hash
void
AtomHash::add_point( core::Size residx, core::Vector xyz ) {
	AtomHashNode w_i;
	w_i.idx = nodes_.size();
	w_i.residx = residx;
	w_i.xyz = xyz;

	int x,y,z;
	x = (int) std::floor( (xyz[0]+0.5)/bindis_ );
	x = x%1024; if ( x<0 ) x+=1024;
	y = (int) std::floor( (xyz[1]+0.5)/bindis_ );
	y = y%1024; if ( y<0 ) y+=1024;
	z = (int) std::floor( (xyz[2]+0.5)/bindis_ );
	z = z%1024; if ( z<0 ) z+=1024;
	int hashid = (x<<20)+(y<<10)+z;

	// ignore redundant points
	core::Real dcut = 0.04; // we're pretty strict here, d<0.2 is redundant
	auto it_cell = nodes_.equal_range( hashid );
	for ( auto it=it_cell.first; it!=it_cell.second; it++ ) {
		core::Real distAct2 = (it->second.xyz - xyz).length_squared();
		if ( distAct2 < dcut ) return;
	}

	nodes_.insert(std::make_pair(hashid, w_i));
}

// get the number of atom neighbors
core::Size
AtomHash::get_neighborcount( core::Vector xyz, core::Real dlim ) {
	int x,y,z;
	x = (int) std::floor( (xyz[0]+0.5)/bindis_ );
	x = x%1024; if ( x<0 ) x+=1024;
	y = (int) std::floor( (xyz[1]+0.5)/bindis_ );
	y = y%1024; if ( y<0 ) y+=1024;
	z = (int) std::floor( (xyz[2]+0.5)/bindis_ );
	z = z%1024; if ( z<0 ) z+=1024;

	core::Size retval = 0;
	core::Real maxDis2 = dlim*dlim;
	int boundcheck = std::ceil( dlim/bindis_ );
	for ( int ix=x-boundcheck; ix<=x+boundcheck; ++ix ) {
		int xx = ix%1024; if ( xx<0 ) xx+=1024;
		for ( int iy=y-boundcheck; iy<=y+boundcheck; ++iy ) {
			int yy = iy%1024; if ( yy<0 ) yy+=1024;
			for ( int iz=z-boundcheck; iz<=z+boundcheck; ++iz ) {
				int zz = iz%1024; if ( zz<0 ) zz+=1024;

				int hashid = (xx<<20)+(yy<<10)+zz;
				auto it_cell = nodes_.equal_range( hashid );
				for ( auto it=it_cell.first; it!=it_cell.second; it++ ) {
					core::Real distAct2 = (it->second.xyz - xyz).length_squared();
					if ( distAct2 <= maxDis2 ) {
						retval++;
					}
				}
			}
		}
	}
	return (retval);
}


// removes all points from the water hash and returns clusters
// 3 different distance criteria checks:
//    1) identify i,j pairs with different resids and d<cut1 (1.5A)
//    2) throw out pairs w.i cut2 of each other (0.35A)
//    3) group into "rotamer clouds" w/i cut3 of each other (3A)
// this function carries out 1
void
AtomHash::trim_to_heterogeneous_clusters(core::Real dlim) {
	int x,y,z;

	utility::vector1<bool> marked( nodes_.size(), false );
	utility::vector1<core::Vector> newpoints;

	for ( auto it1=nodes_.begin(); it1 != nodes_.end(); ++it1 ) {
		x = (it1->first>>20) & 1023;
		y = (it1->first>>10) & 1023;
		z = (it1->first) & 1023;

		core::Real maxDis2 = dlim*dlim;
		int boundcheck = std::ceil( dlim/bindis_ );
		for ( int ix=x-boundcheck; ix<=x+boundcheck; ++ix ) {
			int xx = ix%1024; if ( xx<0 ) xx+=1024;
			for ( int iy=y-boundcheck; iy<=y+boundcheck; ++iy ) {
				int yy = iy%1024; if ( yy<0 ) yy+=1024;
				for ( int iz=z-boundcheck; iz<=z+boundcheck; ++iz ) {
					int zz = iz%1024; if ( zz<0 ) zz+=1024;

					int hashid = (xx<<20)+(yy<<10)+zz;
					auto it_cell = nodes_.equal_range( hashid );
					for ( auto it2=it_cell.first; it2!=it_cell.second; it2++ ) {
						if ( it2->second.idx <= it1->second.idx ) continue; // one-sided check
						if ( it2->second.residx == it1->second.residx ) continue; // diff residues check
						core::Real distAct2 = (it1->second.xyz - it2->second.xyz).length_squared();
						if ( distAct2 > maxDis2 ) continue;
						core::Vector newhit = 0.5 * (it1->second.xyz + it2->second.xyz);
						newpoints.push_back( newhit );
					}
				}
			}
		}
	}

	nodes_.clear();
	for ( auto it=newpoints.begin(); it!=newpoints.end(); ++it ) {
		add_point( 0, *it );
	}
}

// 3 different distance criteria checks:
//    1) identify i,j pairs with different resids and d<cut1 (1.5A)
//    2) throw out pairs w.i cut2 of each other (0.35A)
//    3) group into "rotamer clouds" w/i cut3 of each other (3A)
// this function carries out 2&3
void
AtomHash::reset_and_get_clusters(
	utility::vector1< utility::vector1<core::Vector> > &clusters,
	core::Real dredundant,
	core::Real dcluster
) {
	clusters.clear();
	int x,y,z;
	core::Real dredundant2 = dredundant*dredundant;
	core::Real dcluster2 = dcluster*dcluster;

	while ( nodes_.size()>0 ) {
		auto it1 = nodes_.begin();

		utility::vector1< core::Vector > clust_i;
		clust_i.push_back( it1->second.xyz );

		x = (it1->first>>20) & 1023;
		y = (it1->first>>10) & 1023;
		z = (it1->first) & 1023;

		nodes_.erase(it1);

		// now find all neighboring points
		int boundcheck = std::ceil( std::max(dredundant,dcluster)/bindis_ );
		for ( int ix=x-boundcheck; ix<=x+boundcheck; ++ix ) {
			int xx = ix%1024; if ( xx<0 ) xx+=1024;
			for ( int iy=y-boundcheck; iy<=y+boundcheck; ++iy ) {
				int yy = iy%1024; if ( yy<0 ) yy+=1024;
				for ( int iz=z-boundcheck; iz<=z+boundcheck; ++iz ) {
					int zz = iz%1024; if ( zz<0 ) zz+=1024;
					int hashid = (xx<<20)+(yy<<10)+zz;
					auto it_cell = nodes_.equal_range( hashid );
					for ( auto it2=it_cell.first; it2!=it_cell.second; ) {
						core::Real distAct2 = (clust_i[1] - it2->second.xyz).length_squared();
						if ( distAct2 > dcluster2 ) {
							it2++;
							continue;
						}

						// now make sure this point is not redundant
						bool isnotredund = true;
						for ( auto it_cl=clust_i.begin(); it_cl!=clust_i.end() && isnotredund; ++it_cl ) {
							distAct2 = (it2->second.xyz - *it_cl).length_squared();
							isnotredund = (distAct2>dredundant2);
						}

						if ( isnotredund ) {
							clust_i.push_back( it2->second.xyz );
						}
						it2 = nodes_.erase(it2);
					}
				}
			}
		}
		clusters.push_back( clust_i );
	}
}


void
WaterRotsDB::initialize() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::io::izstream data;
	std::string filename = option[corrections::water::wat_rot_params]();
	basic::database::open( data, filename );

	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream l( line );
		std::string aatm,abase1,abase2, tag, aatag;
		numeric::xyzVector< core::Real> params;

		l >> tag;
		if ( tag == "PROTEIN" ) {
			l >> aatag >> aatm >> abase1 >> abase2 >> params[0] >> params[1] >> params[2];
			if ( l.fail() ) {
				utility_exit_with_message("parse error: "+line);
			}
			protein_rots_.push_back( WaterRot(params,aatag,aatm,abase1,abase2) );
		} else if ( tag == "LIGAND" ) {
			l >> aatag >> params[0] >> params[1] >> params[2];
			if ( l.fail() ) {
				utility_exit_with_message("parse error: "+line);
			}
			ligand_rots_.push_back( WaterRot(params,aatag,"","","") );
		}
	}
	TR << "Read " << protein_rots_.size() << " protein water rotamers and " << ligand_rots_.size() << " ligand water rotamers." << std::endl;
	init_ = true;
}


std::string
WaterBoxMoverCreator::keyname() const {
	return WaterBoxMoverCreator::mover_name();
}

protocols::moves::MoverOP
WaterBoxMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new WaterBoxMover );
}

std::string
WaterBoxMover::mover_name() {
	return WaterBoxMoverCreator::mover_name();
}

std::string
WaterBoxMoverCreator::mover_name() {
	return "WaterBoxMover";
}

void WaterBoxMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	WaterBoxMover::provide_xml_schema( xsd );
}

void WaterBoxMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	auto ct_gen = define_water_box_mover_schema();
	ct_gen->write_complex_type_to_schema(xsd);
}

void
WaterBoxMover::init() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	mode_ = "replace";

	lkb_overlap_dist_ = option[ OptionKeys::corrections::water::lkb_overlap_distance ]();  // default 1.25
	lkb_clust_rad_ = option[ OptionKeys::corrections::water::lkb_cluster_radius ]();       // default 0.35
	lkb_rotset_radius_ = option[ OptionKeys::corrections::water::lkb_rotset_radius ]();    // default 3.0

	clust_radius_ = option[ OptionKeys::corrections::water::cluster_radius ]();  // default 1.0
	clust_cutoff_ = option[ OptionKeys::corrections::water::cluster_cutoff ]();  // default 0.4

	dwell_cutoff_ = option[ OptionKeys::corrections::water::dwell_cutoff ]();  // default 0.12
	watlim_scale_ = option[ OptionKeys::corrections::water::watlim_scale ]();  // default 0.0

	gen_fixed_ = false;  // when enabled, generate water for fixed regions

	if ( option[ in::file::native ].user() ) {
		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::ImportPoseOptions options;
		options.set_ignore_waters(false);
		core::import_pose::pose_from_file( *native_, option[ in::file::native ](), options, false, core::import_pose::PDB_file);
	}
}

std::string
WaterBoxMover::get_name() const {
	return WaterBoxMoverCreator::mover_name();
}

protocols::moves::MoverOP
WaterBoxMover::clone() const {
	return protocols::moves::MoverOP( new WaterBoxMover( *this ) );
}

protocols::moves::MoverOP
WaterBoxMover::fresh_instance() const {
	return protocols::moves::MoverOP( new WaterBoxMover );
}

// helper function to delete (all or virtualized) waters from a pose
void
WaterBoxMover::delete_waters( core::pose::Pose &pose, bool virt_only ) {
	core::Size nres = pose.total_residue();

	for ( core::Size i=nres; i>=1; --i ) {
		core::conformation::Residue const &res = pose.residue(i);
		if ( res.name() == "HOH_V" || res.name() == "PWAT_V" ) {
			pose.conformation().delete_residue_slow( i );
		} else if ( !virt_only && (res.name() == "HOH" || res.name() == "PWAT") ) {
			pose.conformation().delete_residue_slow( i );
		}
	}
	core::pose::renumber_pdbinfo_based_on_conf_chains(pose);
}

/// setup the initial packing step
void
WaterBoxMover::setup_pack( core::pose::Pose & pose ) {
	using namespace core::pack::interaction_graph;

	pose.update_residue_neighbors();
	runtime_assert ( sf_ );

	// process task factory & validate task
	if ( task_factory_ != 0 ) {
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	} else if ( task_ == 0 ) {
		TR << "No packer task defined!  Generating waters on fixed pose only!" << std::endl;
		core::pack::task::PackerTaskOP task_new = core::pack::task::TaskFactory::create_packer_task( pose );
		for ( Size i(1); i <= pose.total_residue(); ++i ) {
			task_new->nonconst_residue_task(i).prevent_repacking();
		}
		task_ = task_new;
		gen_fixed_ = true; // override setting for this (otherwise mover will do nothing)
	} else {
		bool task_valid = true;
		if ( task_->total_residue() != pose.total_residue() ) task_valid=false;
		for ( Size i(1); i <= pose.total_residue(); ++i ) {
			core::chemical::ResidueTypeCOP r = pose.residue_type(i).get_self_ptr();
			if ( ! task_->residue_task(i).is_original_type( r ) ) task_valid=false;
		}
		runtime_assert( task_valid );
	}

	// build initial rotset (no water)
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		rotamer_sets_ = core::pack::rotamer_set::RotamerSetsOP(new core::pack::rotamer_set::symmetry::SymmetricRotamerSets());
	} else {
		rotamer_sets_ = core::pack::rotamer_set::RotamerSetsOP(new core::pack::rotamer_set::RotamerSets());
	}
	sf_->setup_for_packing( pose, task_->repacking_residues(), task_->designing_residues() );
	utility::graph::GraphOP packer_neighbor_graph = core::pack::create_packer_graph( pose, *sf_, task_ );
	rotamer_sets_->set_task( task_ );
	rotamer_sets_->initialize_pose_for_rotsets_creation(pose);
	rotamer_sets_->build_rotamers( pose, *sf_, packer_neighbor_graph );
	//rotamer_sets_->prepare_sets_for_packing( pose, *sf_ );

	// backbone, 1-sided clouds
	PWatRotamerCloud watercloud_bb, watercloud_sc;
	build_backbone_rotamer_clouds( pose, watercloud_bb );

	// rotamer overlap clouds
	build_lkboverlap_rotamer_clouds( pose, watercloud_sc );

	// finally we need to update:
	//   - rotamersets
	//   - packertask
	//   - pose
	// with our new water rotamer clouds
	watercloud_bb.insert( watercloud_bb.end(), watercloud_sc.begin(), watercloud_sc.end() );
	attach_rotamer_clouds_to_pose_and_rotset( pose, watercloud_bb );

	// finally everything is resolved, regenerate derived data and pack
	rotamer_sets_->prepare_sets_for_packing( pose, *sf_ );
	sf_->setup_for_packing( pose, task_->repacking_residues(), task_->designing_residues() );
	packer_neighbor_graph = core::pack::create_packer_graph( pose, *sf_, task_ );

	ig_ = InteractionGraphFactory::create_and_initialize_annealing_graph(
		*task_, *rotamer_sets_, pose, *sf_, packer_neighbor_graph );
}


void
WaterBoxMover::build_backbone_rotamer_clouds(
	core::pose::Pose const & pose,
	PWatRotamerCloud & retval
) {
	using core::conformation::Residue;
	using core::conformation::ResidueOP;
	using numeric::conversions::radians;

	retval.clear();

	bool dbloaded = false;

	{ // Scope for read lock
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard readlock(db_mutex_);
#endif
		dbloaded = water_rots_db_.is_initialized();
	}

	if ( !dbloaded ) {
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writelock(db_mutex_);
#endif
		if ( !water_rots_db_.is_initialized() ) {  // check again
			water_rots_db_.initialize();
		}
	}

#ifdef MULTI_THREADED
	// now get read lock for remainder of this function
	utility::thread::ReadLockGuard readlock(db_mutex_);
#endif

	// now build waters for all residues
	// TODO: do we apply the taskop here????
	core::Size nwaterres=0, nwaterrot=0;
	for ( int i=1; i<=(int)pose.total_residue(); ++i ) {
		core::conformation::Residue const &res_i = pose.residue(i);

		// split protein/ligand logic
		if ( res_i.is_protein() ) {
			if ( !gen_fixed_ && !task_->design_residue( i ) && !task_->pack_residue( i ) ) continue;

			// store water clouds indexed by base atom type
			std::map< core::Size, utility::vector1< core::Vector > > waters_i;

			for ( int k=1; k<=(int)water_rots_db_.n_protein_rots(); ++k ) {
				WaterRot const & wat_gen_k = water_rots_db_.protein_rot(k);

				// rsd typecheck
				if ( wat_gen_k.aa_ != "*" && core::chemical::aa_from_name(wat_gen_k.aa_) != res_i.aa() ) continue;

				// atom typecheck
				if ( !res_i.has( wat_gen_k.aatm_ ) || !res_i.has( wat_gen_k.abase1_ ) || !res_i.has( wat_gen_k.abase2_ ) ) continue;

				core::Size aatm_idx = res_i.atom_index( wat_gen_k.aatm_ );
				core::Size abase1_idx = res_i.atom_index( wat_gen_k.abase1_ );
				core::Size abase2_idx = res_i.atom_index( wat_gen_k.abase2_ );

				core::kinematics::Stub stub( res_i.xyz(aatm_idx), res_i.xyz(abase1_idx), res_i.xyz(abase2_idx) );
				if ( waters_i.find( aatm_idx ) == waters_i.end() ) {
					waters_i[aatm_idx] = utility::vector1< core::Vector >();
				}

				waters_i[aatm_idx].push_back( stub.spherical(
					radians( wat_gen_k.coords_[2] ),
					radians( 180.0 - wat_gen_k.coords_[1] ),
					wat_gen_k.coords_[0]
					) );
			}

			for ( auto mapit=waters_i.begin(); mapit != waters_i.end(); ++mapit ) {
				if ( mapit->second.size() > 0 ) { // this check is unnecessary I think...
					retval.push_back( mapit->second );
					nwaterres++;
					nwaterrot+=mapit->second.size();
				}
			}
		} else if ( res_i.is_ligand() || res_i.is_NA() ) {
			// for ligands individually loop over donor/acceptor positions
			for ( core::chemical::AtomIndices::const_iterator
					hnum  = res_i.Hpos_polar().begin(),
					hnume = res_i.Hpos_polar().end(); hnum != hnume; ++hnum ) {
				core::Size hatm(*hnum);
				core::Size datm( res_i.atom_base( hatm ) );
				core::Size datm_base( res_i.atom_base( datm ) );

				// sanity check for HOH
				if ( datm_base == hatm ) {
					datm_base = 0;
					core::chemical::AtomIndices const & datm_nbrs( res_i.type().nbrs( datm ) );
					for ( Size ii=1; ii<= datm_nbrs.size(); ++ii ) {
						if ( datm_nbrs[ii] == hatm ) continue;
						else datm_base = datm_nbrs[ii];
					}
					runtime_assert( datm_base );
				}

				utility::vector1< core::Vector > waters_i;
				for ( int k=1; k<=(int)water_rots_db_.n_ligand_rots(); ++k ) {
					WaterRot const & wat_gen_k = water_rots_db_.ligand_rot(k);
					if ( wat_gen_k.aa_ != "DON" ) continue;

					core::kinematics::Stub stub( res_i.xyz(hatm), res_i.xyz(datm), res_i.xyz(datm_base) );
					waters_i.push_back( stub.spherical(
						radians( wat_gen_k.coords_[2] ),
						radians( 180.0 - wat_gen_k.coords_[1] ),
						wat_gen_k.coords_[0]
						) );
				}

				if ( waters_i.size() > 0 ) { // this one is necessary though
					retval.push_back( waters_i );
					nwaterres++;
					nwaterrot+=waters_i.size();
				}
			}

			for ( core::chemical::AtomIndices::const_iterator
					anum  = res_i.accpt_pos().begin(),
					anume = res_i.accpt_pos().end(); anum != anume; ++anum ) {
				core::Size const aatm( *anum );
				core::Size const abase1( res_i.atom_base( aatm ) ), abase2( res_i.abase2( aatm ) );
				core::chemical::Hybridization hyb = res_i.atom_type( *anum ).hybridization();
				std::string tgt_name;
				if ( hyb == core::chemical::SP2_HYBRID ) tgt_name = "SP2_ACC";
				else if ( hyb == core::chemical::SP3_HYBRID ) tgt_name = "SP3_ACC";
				else if ( hyb == core::chemical::RING_HYBRID ) tgt_name = "RING_ACC";
				else utility_exit_with_message("Unknown hybridization!");

				utility::vector1< core::Vector > waters_i;
				for ( int k=1; k<=(int)water_rots_db_.n_ligand_rots(); ++k ) {
					WaterRot const & wat_gen_k = water_rots_db_.ligand_rot(k);
					if ( wat_gen_k.aa_ != tgt_name ) continue;

					if ( hyb == core::chemical::RING_HYBRID ) {
						numeric::xyzVector< core::Real > virtual_abase1 = 0.5 * ( res_i.xyz(abase1) + res_i.xyz(abase2) );
						core::kinematics::Stub stub( res_i.xyz(aatm), virtual_abase1, res_i.xyz(abase2) );
						waters_i.push_back( stub.spherical(
							radians( wat_gen_k.coords_[2] ),
							radians( 180.0 - wat_gen_k.coords_[1] ),
							wat_gen_k.coords_[0]
							) );
					} else {
						core::kinematics::Stub stub( res_i.xyz(aatm), res_i.xyz(abase1), res_i.xyz(abase2) );
						waters_i.push_back( stub.spherical(
							radians( wat_gen_k.coords_[2] ),
							radians( 180.0 - wat_gen_k.coords_[1] ),
							wat_gen_k.coords_[0]
							) );
					}
				}

				if ( waters_i.size() > 0 ) { // this one is necessary though
					retval.push_back( waters_i );
					nwaterres++;
					nwaterrot+=waters_i.size();
				}
			}
		}
	}

	TR << "build_backbone_rotamer_clouds generated " << nwaterrot << " rotamers in " << nwaterres << " clusters" << std::endl;
}


void
WaterBoxMover::build_lkboverlap_rotamer_clouds(
	core::pose::Pose const & pose,
	PWatRotamerCloud & retval
) {
	retval.clear();

	// (0) build protein occupancy map for fast neighborchecks
	//core::Real neighbrad = 10.0;
	//core::Size neighbcountcut = 16;
	//AtomHash allprotein( neighbrad+0.5 );
	// for ( Size i=1; i <= pose.total_residue(); ++i ) {
	//  if (!pose.residue(i).is_water()) {
	//   allprotein.add_point( i, pose.residue(i).xyz( pose.residue(i).nbr_atom() ) );
	//  }
	//}

	// (1) get lkball positions from packable side chains
	AtomHash allwaters( lkb_rotset_radius_+1.0 );

	for ( Size i=1; i <= task_->total_residue(); ++i ) {

		if ( task_->design_residue( i ) || task_->pack_residue( i ) ) {
			using namespace core::pack::rotamer_set;
			RotamerSetOP res_rotset = rotamer_sets_->rotamer_set_for_residue( i );

			for ( Size ii=1; ii <= res_rotset->num_rotamers(); ++ii ) {
				core::conformation::Residue res_rot = *res_rotset->rotamer( ii );
				core::scoring::lkball::LKB_ResidueInfoOP lkb_resinfo( new core::scoring::lkball::LKB_ResidueInfo( res_rot ) );

				utility::vector1< core::Size > rsdin_waters = lkb_resinfo->n_attached_waters();
				utility::vector1< core::scoring::lkball::WaterCoords > rsdiwaters = lkb_resinfo->waters();

				for ( Size j=1; j <= rsdin_waters.size(); ++j ) {
					for ( Size k=1; k <= rsdin_waters[j]; ++k ) {
						// "exposure" check
						//if ( allprotein.get_neighborcount( rsdiwaters[j][k], neighbrad ) > neighbcountcut) {
						allwaters.add_point( i, rsdiwaters[j][k] );
						//}
					}
				}
			}
		} else if ( gen_fixed_ ) {  // residue not being packed/designed
			// get lkball positions for fixed residues
			core::scoring::lkball::LKB_ResidueInfoOP lkb_resinfo( new core::scoring::lkball::LKB_ResidueInfo( pose.residue(i) ) );

			utility::vector1< core::Size > rsdin_waters = lkb_resinfo->n_attached_waters();
			utility::vector1< utility::fixedsizearray1< core::Vector, 4 > > rsdiwaters = lkb_resinfo->waters();

			for ( Size j=1; j <= rsdiwaters.size(); ++j ) {
				for ( Size k=1; k <= rsdin_waters[j]; ++k ) {
					// "exposure" check
					//if ( allprotein.get_neighborcount( rsdiwaters[j][k], neighbrad ) > neighbcountcut) {
					allwaters.add_point( i, rsdiwaters[j][k] );
					//}
				}
			}
		}
	}  // loop over nres
	TR << "build_lkboverlap_rotamer_clouds initially generated " << allwaters.npoints() << " waters" << std::endl;

	// (2) identify overlaps
	allwaters.trim_to_heterogeneous_clusters(lkb_overlap_dist_);
	TR << "build_lkboverlap_rotamer_clouds identified " << allwaters.npoints() << " overlap sites" << std::endl;

	// (3) cluster into "rotamers"
	allwaters.reset_and_get_clusters( retval, lkb_clust_rad_, lkb_rotset_radius_ );
	core::Size npts=0;
	for ( core::Size i=1; i<=retval.size(); ++i ) npts+=retval[i].size();
	TR << "build_lkboverlap_rotamer_clouds generated " << npts << " rotamers in " << retval.size() << " clusters" << std::endl;
}

// updates: rotamersets, task, and pose with additional clouds
//   applies bump check
void
WaterBoxMover::attach_rotamer_clouds_to_pose_and_rotset(
	core::pose::Pose &pose,
	PWatRotamerCloud &watercloud
) {
	using namespace core::scoring;
	using namespace core::conformation;

	core::Size nposorig = pose.total_residue();
	core::Size nposnew = watercloud.size();

	// O preload some data
	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
	core::conformation::ResidueOP pwat = core::conformation::ResidueFactory::create_residue( rsd_set->name_map("PWAT") );  // point water
	core::conformation::ResidueOP vwat = core::conformation::ResidueFactory::create_residue( rsd_set->name_map("PWAT_V") ); // virtualized
	core::Vector const OH1( pwat->xyz("V1") - pwat->xyz("O") );
	core::Vector const OH2( pwat->xyz("V2") - pwat->xyz("O") );

	// 1 update the pose
	for ( core::Size i=1; i<=nposnew; ++i ) {
		core::conformation::ResidueOP vrt_wat = core::conformation::ResidueOP( new core::conformation::Residue( *vwat ) );
		vrt_wat->set_xyz("O",  watercloud[i][1]);
		vrt_wat->set_xyz("V1", watercloud[i][1]+OH1);
		vrt_wat->set_xyz("V2", watercloud[i][1]+OH2);
		pose.append_residue_by_jump( *vrt_wat, nposorig ); //?
	}
	//pose.dump_pdb("out1.pdb");

	// 2 update the packertask
	task_ = update_packer_task( pose, task_ );

	// 3 add residues to the rotset
	rotamer_sets_->set_task( task_ );

	//  3b setup some stuff for bump check
	utility::graph::GraphOP packer_neighbor_graph = core::pack::create_packer_graph( pose, *sf_, task_ );

	// NOTE NOTE NOTE:
	//   this does not try to regenetate rotamers, as set task only calls resize and not clear()
	//   I'm not sure if this is intentional however, so if that behavior changes,
	//   protein rotamers will need to be regenerated here

	core::Size nrot=0, nfailbump=0;
	for ( core::Size i=1; i<=nposnew; ++i ) {
		//utility::vector1< core::conformation::ResidueOP > rotset_i;
		core::Size resid = nposorig+i;
		core::pack::rotamer_set::RotamerSetOP rotset_i( core::pack::rotamer_set::RotamerSetFactory::create_rotamer_set( pose ) );
		rotset_i->set_resid( resid );

		for ( core::Size j=1; j<=watercloud[i].size(); ++j ) {
			core::conformation::ResidueOP pt_wat = core::conformation::ResidueOP( new core::conformation::Residue( *pwat ) );
			pt_wat->set_xyz("O",  watercloud[i][j]);
			pt_wat->set_xyz("V1", watercloud[i][j]+OH1);
			pt_wat->set_xyz("V2", watercloud[i][j]+OH2);

			// bumpcheck
			EnergyMap emap;
			for ( utility::graph::Graph::EdgeListConstIter
					ir  = packer_neighbor_graph->get_node( resid )->const_edge_list_begin(),
					ire = packer_neighbor_graph->get_node( resid )->const_edge_list_end();
					ir != ire; ++ir ) {
				int const neighbor_id( (*ir)->get_other_ind( resid ) );
				if ( neighbor_id >=(int)nposorig ) continue;
				Residue const & neighbor( pose.residue( neighbor_id ) );

				if ( ! task_->pack_residue( neighbor_id ) ) {
					sf_->bump_check_full( *pt_wat, neighbor, pose, emap);
				} else {
					sf_->bump_check_backbone( *pt_wat, neighbor, pose, emap);
				}
			}
			core::Real bumpE = (sf_->weights().dot( emap ));

			if ( bumpE < 0 ) {
				rotset_i->add_rotamer( *pt_wat );
				nrot++;
			} else {
				nfailbump++;
			}
		}

		// add current (the virtualized residue) as the last rotamer (this ensures we always have >=1 rotamer as well)
		core::conformation::ResidueOP vrt_wat = core::conformation::ResidueOP(
			new core::conformation::Residue( pose.residue( resid ) ) );
		rotset_i->add_rotamer( *vrt_wat );
		nrot++;

		// and update
		core::Size moltenid = rotamer_sets_->resid_2_moltenres( resid );
		rotamer_sets_->set_explicit_rotamers( moltenid , rotset_i );
	}

	TR << "Built " << nrot << " water rotamers at " << nposnew << " positions (" << nfailbump << " fail bumpcheck)" << std::endl;

	// we're done adding rotamers, update derived data
	rotamer_sets_->update_offset_data();
}

/// run packing trajectory
void
WaterBoxMover::run_pack(
	core::pose::Pose & pose,
	utility::vector0< int > rot_to_pack,
	utility::vector1< core::pack::annealer::PointDwell > & all_rot
) {
	ObjexxFCL::FArray1D_int bestrotamer_at_seqpos( pose.total_residue() );
	core::PackerEnergy bestenergy( 0.0 );
	bool start_with_current = false;
	ObjexxFCL::FArray1D_int current_rot_index( pose.total_residue(), 0 );
	ObjexxFCL::FArray1D< core::PackerEnergy > rot_freq( ig_->get_num_total_states(), 0.0 );

	core::pack::annealer::FixbbPwatSimAnnealer annealer(
		rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig_,
		rotamer_sets_, current_rot_index, false, rot_freq );

	annealer.set_min_dwell( dwell_cutoff_ );

	time_t const pwatanneal_start = time(NULL);
	annealer.run();
	TR << "time spent packing pwat = " << time(NULL) - pwatanneal_start << " seconds -- runtime" << std::endl;

	all_rot = annealer.get_dwell_times();
}

/// update task to include all waters
core::pack::task::PackerTaskOP
WaterBoxMover::update_packer_task(
	Pose const & pose,
	core::pack::task::PackerTaskCOP & packer_task
) {
	utility::vector1< bool > residues_allowed_to_be_packed;
	Size nres = pose.total_residue();
	for ( Size i=1; i <= nres; ++i ) {
		std::string resname = pose.residue(i).name();
		if ( resname == "PWAT" || resname == "PWAT_V" || resname == "HOH" || resname == "HOH_V" ) {
			residues_allowed_to_be_packed.push_back(true);
		} else if ( packer_task->pack_residue(i) ) {
			residues_allowed_to_be_packed.push_back(true);
		} else {
			residues_allowed_to_be_packed.push_back(false);
		}
	}

	core::pack::task::PackerTaskOP updated_task( new core::pack::task::PackerTask_( pose ) );
	updated_task->restrict_to_residues(residues_allowed_to_be_packed);
	updated_task->restrict_to_repacking();
	updated_task->initialize_from_command_line();

	// copy extra_chi information from original task to new task
	for ( Size i=1; i <= packer_task->total_residue(); ++i ) {
		if ( i > nres ) continue; // protection in case of smaller pose than what entered the WaterDdg mover
		updated_task->nonconst_residue_task( i ).or_ex1( packer_task->residue_task( i ).ex1() );
		updated_task->nonconst_residue_task( i ).or_ex2( packer_task->residue_task( i ).ex2() );
		updated_task->nonconst_residue_task( i ).or_ex3( packer_task->residue_task( i ).ex3() );
		updated_task->nonconst_residue_task( i ).or_ex4( packer_task->residue_task( i ).ex4() );
		updated_task->nonconst_residue_task( i ).or_ex1aro( packer_task->residue_task( i ).ex1aro() );
		updated_task->nonconst_residue_task( i ).or_ex2aro( packer_task->residue_task( i ).ex2aro() );
		updated_task->nonconst_residue_task( i ).or_ex1aro_exposed( packer_task->residue_task( i ).ex1aro_exposed() );
		updated_task->nonconst_residue_task( i ).or_ex2aro_exposed( packer_task->residue_task( i ).ex2aro_exposed() );
		updated_task->nonconst_residue_task( i ).and_extrachi_cutoff( packer_task->residue_task( i ).extrachi_cutoff() );
		updated_task->nonconst_residue_task( i ).or_include_current( packer_task->residue_task( i ).include_current() );
	}

	return updated_task;
}


/// cluster water clouds
/// inefficient but probably okay (there should not be many waters
///      passing filters else rotatible waters would take excessively long)
void
WaterBoxMover::cluster_rotset( utility::vector1< core::pack::annealer::PointDwell > &rotset ) {
	if ( rotset.size() == 0 ) return;

	TR << "Clustering " << rotset.size() << " PWAT rotamers " << std::endl;

	utility::vector1< core::pack::annealer::PointDwell > cluster_centers;
	utility::vector1< core::pack::annealer::PointDwell > new_cluster_centers;
	cluster_centers.push_back(rotset[1]);
	core::pack::annealer::PointDwell center_i = rotset[1];
	center_i.xyz *= center_i.dwell;
	new_cluster_centers.push_back(center_i);

	// pass 1: rough cluster assignment
	//   - first member added to a cluster is starting point
	//   - "soft" membership: points may belong to multiple clusters
	for ( Size i = 2; i <= rotset.size(); ++i ) {
		bool unique = true;
		for ( Size j = 1; j <= cluster_centers.size(); ++j ) {
			core::Vector const &clust_cent = cluster_centers[j].xyz;
			if ( (clust_cent - rotset[i].xyz).length_squared() <= clust_radius_*clust_radius_ ) {
				new_cluster_centers[j].dwell += rotset[i].dwell;
				new_cluster_centers[j].xyz += rotset[i].dwell * rotset[i].xyz;
				unique = false;
			}
		}
		if ( unique ) {
			cluster_centers.push_back(rotset[i]);
			center_i = rotset[i];
			center_i.xyz *= center_i.dwell;
			new_cluster_centers.push_back(center_i);
		}
	}
	for ( Size j = 1; j <= cluster_centers.size(); ++j ) {
		new_cluster_centers[j].xyz /= new_cluster_centers[j].dwell;
	}

	// pass 2+: iterate (EM)
	//   - use centroids from previous round
	//   - "hard" membership: points may only belong to a single cluster
	bool done=false;
	while ( !done ) {
		done = true;
		cluster_centers = new_cluster_centers;

		for ( Size j = 1; j <= cluster_centers.size(); ++j ) {
			new_cluster_centers[j].dwell = 0;
			new_cluster_centers[j].xyz = 0;
		}

		for ( Size i = 1; i <= rotset.size(); ++i ) {
			core::Real mindist = 1e6;
			core::Size minclust = 1;
			for ( Size j = 1; j <= cluster_centers.size(); ++j ) {
				core::Real dist2 = (cluster_centers[j].xyz - rotset[i].xyz).length_squared();
				if ( dist2 < mindist ) {
					mindist = dist2;
					minclust = j;
				}
			}

			if ( mindist <= clust_radius_*clust_radius_ ) {
				new_cluster_centers[minclust].dwell += rotset[i].dwell;
				new_cluster_centers[minclust].xyz += rotset[i].dwell * rotset[i].xyz;
			} else {
				cluster_centers.push_back(rotset[i]);
				center_i = rotset[i];
				center_i.xyz *= center_i.dwell;
				new_cluster_centers.push_back(center_i);
				done = false;
			}
		}

		// check if changed
		for ( Size j = 1; j <= cluster_centers.size(); ++j ) {
			if ( std::fabs( new_cluster_centers[j].dwell - cluster_centers[j].dwell ) > 1e-6 ) {
				done = false;
			}
		}

		// delete unoccupied
		for ( auto it = new_cluster_centers.begin(); it != new_cluster_centers.end(); ) {
			if ( it->dwell == 0 ) {
				it = new_cluster_centers.erase( it );
			} else {
				it++;
			}
		}

		for ( Size j = 1; j <= new_cluster_centers.size(); ++j ) {
			new_cluster_centers[j].xyz /= new_cluster_centers[j].dwell;
		}
	}

	TR << "Total number of clusters = " << cluster_centers.size() << std::endl;

	// filter on clusters that don't meet cumulative dwell time cutoff
	rotset.clear();
	for ( Size i = 1; i <= cluster_centers.size(); ++i ) {
		if ( cluster_centers[i].dwell >= clust_cutoff_ ) {
			rotset.push_back(cluster_centers[i]);
		}
	}
	TR << "Total number of clusters passing dwell-time filter = " << rotset.size() << std::endl;
}


/// find closest protein residue
core::Size
WaterBoxMover::find_closest( core::pose::Pose const & pose, core::Vector Ocoord ) {
	core::Real mindist = 1e6;
	core::Size attach_to = 0;
	for ( core::Size i(1); i <= pose.total_residue(); ++i ) {
		core::conformation::Residue const &res = pose.residue(i);
		if ( res.is_water() ) continue;
		core::Real dist = ( Ocoord-res.atom( res.nbr_atom() ).xyz() ).length_squared();
		if ( dist < mindist ) {
			mindist = dist;
			attach_to = i;
		}
	}
	TR << "Attaching new water to residue " << attach_to << " with distance = " << mindist << std::endl;
	return attach_to;
}

void
WaterBoxMover::get_water_recovery( core::pose::Pose const & pose, bool incl_vrt ) {
	if ( !native_ ) return;

	utility::vector1< core::Vector > native_hoh, predicted_hoh;
	for ( core::Size i=1; i<pose.total_residue(); ++i ) {
		bool is_wat = ( pose.residue(i).name() == "PWAT" || pose.residue(i).name() == "HOH" );
		bool is_vwat = ( pose.residue(i).name() == "PWAT_V" || pose.residue(i).name() == "HOH_V" );
		if ( is_wat || (incl_vrt&&is_vwat) ) {
			predicted_hoh.push_back( pose.residue(i).xyz(1) );
		}
	}

	for ( core::Size i=1; i<native_->total_residue(); ++i ) {
		if ( native_->residue(i).name() == "PWAT" || native_->residue(i).name() == "HOH" ) {
			native_hoh.push_back( native_->residue(i).xyz(1) );
		}
	}

	core::Real ncorr_05 = 0, ncorr_10 = 0;
	for ( core::Size i=1; i<native_hoh.size(); ++i ) {
		core::Real mindist_i = 100.0;
		for ( core::Size j=1; j<predicted_hoh.size(); ++j ) {
			core::Real dist_ij = (predicted_hoh[j] - native_hoh[i]).length_squared();
			mindist_i = std::min( mindist_i, dist_ij );
		}
		if ( mindist_i<0.5*0.5 ) ncorr_05+=1.0;
		if ( mindist_i<1.0*1.0 ) ncorr_10+=1.0;
	}

	core::Size npred = predicted_hoh.size();
	core::Size nnative = native_hoh.size();
	TR << "Number of predicted/native waters = " << npred << "/" << nnative << std::endl;
	TR << "Prec/Rec (0.5A) = " << (npred==0?0:ncorr_05/npred) << "/" << (nnative==0?0:ncorr_05/nnative) << std::endl;
	TR << "Prec/Rec (1.0A) = " << (npred==0?0:ncorr_10/npred) << "/" << (nnative==0?0:ncorr_10/nnative) << std::endl;
}

// main apply function
void
WaterBoxMover::apply( Pose & pose ) {
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pack::task;

	if ( sf_->get_weight( core::scoring::pointwater ) == 0 ) {
		TR << "Setting pointwater weight to 1.0." << std::endl;
		sf_->set_weight( core::scoring::pointwater, 1.0 );
	}

	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );

	// clean up existing waters
	if ( mode_ == "remove" ) {
		delete_waters(pose, false);
		return;
	}

	if ( mode_ == "eval" ) {
		get_water_recovery( pose );
		return;
	}

	// mode is either replace or append
	if ( mode_ == "replace" ) {
		delete_waters(pose, false);
	}

	core::Size orig_nres = pose.total_residue();

	// setup packing task
	Pose working_pose = pose;
	this->setup_pack( working_pose );

	utility::vector1< core::pack::annealer::PointDwell > all_rot;
	utility::vector0< int > rot_to_pack;
	this->run_pack( working_pose, rot_to_pack, all_rot );

	// truncate rotamer list to those with dwell times within the defined cutoff
	utility::vector1< core::pack::annealer::PointDwell > kept_rots;
	for ( Size x = 1; x <= all_rot.size(); ++x ) {
		if ( all_rot[x].dwell >= dwell_cutoff_ ) {
			kept_rots.push_back(all_rot[x]);
		}
	}

	std::sort(kept_rots.begin(), kept_rots.end(),
		[](core::pack::annealer::PointDwell const &a, core::pack::annealer::PointDwell const &b) {return a.dwell > b.dwell;});

	TR << "Rotamers with dwell time > " << dwell_cutoff_ << " = " << kept_rots.size() << std::endl;
	for ( Size x = 1; x <= kept_rots.size(); ++x ) {
		TR << "   " << x << " " << kept_rots[x].xyz.to_string() << " " << kept_rots[x].dwell << std::endl;
	}

	// cluster and resort
	cluster_rotset(kept_rots);
	std::sort(kept_rots.begin(), kept_rots.end(),
		[](core::pack::annealer::PointDwell const &a, core::pack::annealer::PointDwell const &b) {return a.dwell > b.dwell;});

	// kind of hacky .. should query orig task, but this option is not very useful...
	core::Size nkeep = kept_rots.size();
	if ( watlim_scale_ > 0 ) {
		nkeep = (core::Size) std::ceil( watlim_scale_ * orig_nres );
		if ( nkeep < kept_rots.size() ) {
			TR << "Limiting generated waters to " << nkeep << std::endl;
		}
	}

	// all that work to go back to original pose
	//   but now we add centroids as rotatable waters
	for ( Size x = 1; x <= nkeep; ++x ) {
		TR << "  centroid #" << x << ": " << kept_rots[x].xyz.to_string() << " with cumulative dwell time of " << kept_rots[x].dwell << std::endl;

		ResidueOP vrt_wat = ResidueFactory::create_residue( rsd_set->name_map("HOH_V") );
		core::Vector const OH1( vrt_wat->xyz("H1") - vrt_wat->xyz("O") );
		core::Vector const OH2( vrt_wat->xyz("H2") - vrt_wat->xyz("O") );
		ResidueOP new_res = ResidueOP( new Residue( *vrt_wat ) );
		new_res->set_xyz("O", kept_rots[x].xyz);
		new_res->set_xyz("H1", kept_rots[x].xyz+OH1);
		new_res->set_xyz("H2", kept_rots[x].xyz+OH2);
		core::Size attach_to = find_closest( pose, kept_rots[x].xyz );
		pose.append_residue_by_jump( *new_res, attach_to );

		pose.pdb_info()->set_resinfo( pose.total_residue(), pose.pdb_info()->chain(attach_to) , pose.total_residue(), ' ');
	}

	get_water_recovery( pose );

	core::pose::renumber_pdbinfo_based_on_conf_chains(pose); // rep: update pdb_info to account for added waters
	(*sf_)(pose); // some downstream stuff needs this....
}


void
WaterBoxMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &  ) {

	using namespace core::conformation;
	using namespace core::pack::task;

	sf_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );

	if ( tag->hasOption( core::pack::task::TASK_OPERATIONS_TAG ) ) {
		task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, datamap );
	}

	if ( tag->hasOption( "mode" ) ) {
		mode_ = tag->getOption<std::string>("mode");
		if ( mode_ != "remove" && mode_ != "append" && mode_ != "replace" && mode_ != "eval" ) {
			utility_exit_with_message("Unsupported mode!  Must be one of remove, append, replace, or eval");
		}
	}

	if ( tag->hasOption( "lkb_overlap_dist" ) ) {
		lkb_overlap_dist_ = tag->getOption<core::Real>("lkb_overlap_dist");
	}
	if ( tag->hasOption( "lkb_clust_rad" ) ) {
		lkb_clust_rad_ = tag->getOption<core::Real>("lkb_clust_rad");
	}
	if ( tag->hasOption( "lkb_rotset_radius" ) ) {
		lkb_rotset_radius_ = tag->getOption<core::Real>("lkb_rotset_radius");
	}

	if ( tag->hasOption( "dwell_cutoff" ) ) {
		dwell_cutoff_ = tag->getOption<core::Real>("dwell_cutoff");
	}


	if ( tag->hasOption( "clust_radius" ) ) {
		clust_radius_ = tag->getOption<core::Real>("clust_radius");
	}
	if ( tag->hasOption( "clust_cutoff" ) ) {
		clust_cutoff_ = tag->getOption<core::Real>("clust_cutoff");
	}
	if ( tag->hasOption( "gen_fixed" ) ) {
		gen_fixed_ = tag->getOption<bool>("gen_fixed");
	}

	if ( tag->hasOption( "native" ) ) {
		native_ = core::pose::PoseOP( new core::pose::Pose );
		std::string nativepdb = tag->getOption<std::string>("native");
		core::import_pose::ImportPoseOptions options;
		options.set_ignore_waters(false);
		core::import_pose::pose_from_file( *native_, nativepdb, options, false, core::import_pose::PDB_file);
	}
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
WaterBoxMover::define_water_box_mover_schema() {
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute( "mode", xs_string, "The mode to run in (one of: remove, append, replace, or eval)");
	attlist + XMLSchemaAttribute( "lkb_overlap_dist", xsct_real, "The overlap distance of sidechain lkb groups necessary to generate a water");
	attlist + XMLSchemaAttribute( "lkb_clust_rad", xsct_real, "The redundancy cutoff for sidechain lkb groups");
	attlist + XMLSchemaAttribute( "lkb_rotset_radius", xsct_real, "The rotameric cutoff for lkb groups (that is the distance cutoff for which two point waters are considered rotamers of a single water)");
	attlist + XMLSchemaAttribute( "dwell_cutoff", xsct_real, "The dwell cutoff for water packing (before clustering)");
	attlist + XMLSchemaAttribute( "clust_radius", xsct_real, "The radius of clustering of packed waters");
	attlist + XMLSchemaAttribute( "clust_cutoff", xsct_real, "The dwell cutoff for water packing (after clustering)");
	attlist + XMLSchemaAttribute( "gen_fixed", xsct_rosetta_bool, "If set, consider generating waters over fixed regions");
	attlist + XMLSchemaAttribute( "native", xs_string, "Native pose (overrides in:file:native)");

	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );

	ct_gen->complex_type_naming_func(&moves::complex_type_name_for_mover)
		.element_name( mover_name() )
		.description(
		"This mover solvates a pose based on statistics of where waters are found in hi-res crystal "
		"structures about polar side chains and backbone groups.  The input task operation is used to"
		"generate potential hydration sites; subsequent packing steps allow waters to be added at these"
		"positions.")
		.add_attributes( attlist )
		.add_optional_name_attribute();
	return ct_gen;
}

} // moves
} // protocols
