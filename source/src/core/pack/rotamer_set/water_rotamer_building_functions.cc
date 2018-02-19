// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/water_rotamer_building_functions.hh
/// @brief  a few water rotamer building functions
/// @author Frank DiMaio & Ryan Pavlovicz

// Unit Headers
#include <core/pack/rotamer_set/water_rotamer_building_functions.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/types.hh>

#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/Stub.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/graph/Graph.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/database/open.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/constants.hh>
#include <numeric/UniformRotationSampler.hh>
#include <numeric/AxisRotationSampler.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// C++ headers
#include <string>
#include <iostream>
#include <fstream>


namespace core {
namespace pack {
namespace rotamer_set {

static basic::Tracer TR("core.pack.rotamer_set.water_rotamer_building_functions");

// ugly
inline
std::string
strip_whitespace( std::string const & name ) {
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;

	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to dothis?
	return trimmed_name;
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

/// duplicate code -- consolidate?
/// @details Bump check does not include long range energies,
/// though, maybe this should change.
core::PackerEnergy
bump_check(
	core::conformation::ResidueCOP rotamer,
	core::Size resid,
	scoring::ScoreFunction const & sf,
	pose::Pose const & pose,
	task::PackerTask const & task,
	utility::graph::GraphCOP packer_neighbor_graph
)
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

//fpd
void
build_rotated_water_rotamers(
	Size const seqpos,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers,
	bool incl_vrt
) {
	using conformation::Residue;
	using conformation::ResidueOP;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	Residue const & currres( pose.residue(seqpos) );
	core::Vector const O( currres.xyz("O") );
	core::Vector const H1( currres.xyz("H1") );
	core::Vector const H2( currres.xyz("H2") );

	ResidueOP rsd_canonic = conformation::ResidueFactory::create_residue( residue_set->name_map("HOH") );
	core::Vector OH1c( rsd_canonic->xyz("H1") - rsd_canonic->xyz("O") );
	core::Vector OH2c( rsd_canonic->xyz("H2") - rsd_canonic->xyz("O") );

	// uniformly sample rotations
	//
	core::Real rotstep=option[corrections::water::wat_rot_sampling]();
	bool sample_axis=option[corrections::water::wat_axis_sampling]();
	core::Real axis_rotstep=option[corrections::water::axis_rot_sampling]();
	bool wat_bump=option[corrections::water::water_bump_check]();

	core::Real const CO_cutoff = 3.5; // pretty lenient?, look into shortening

	// collected all backbone O atoms within 3.5 A of the water position
	utility::vector1< core::Vector > COs;
	if ( sample_axis ) {
		Size nres = pose.total_residue();
		for ( Size i=1; i <= nres; ++i ) {
			Residue const &res = pose.residue(i);
			if ( res.is_protein() ) {
				Real dist = O.distance( res.xyz("O") );
				if ( dist <= CO_cutoff ) {
					COs.push_back( res.xyz("O") );
				}
			}
		}
	}

	// if sample_axis is set and the water O atom is within 3.5 A
	// of a C=O oxygen, then the O-H1 vector of the water is aligned with the
	// C=O oxygen, and rotamers a built at 'axis_rot_sampling' intervals
	// about the O-O axis
	if ( COs.size() > 0 ) {
		for ( Size i=1; i <= COs.size(); ++i ) {
			Vector axis = ( O - COs[i] );

			Vector init_OH1 = ( currres.xyz("H1") - currres.xyz("O") );
			Vector init_OH2 = ( currres.xyz("H2") - currres.xyz("O") );
			// align O-H1 water vector to O-O axis
			numeric::xyzMatrix<Real> align_matrix = ( numeric::xyzTransform<Real>::align( axis, O-H1 ) ).R ;
			OH1c = align_matrix*init_OH1;
			OH2c = align_matrix*init_OH2;

			// apply rotation matrices about O-O axis to water
			numeric::AxisRotationSampler ars( axis, axis_rotstep );

			numeric::xyzMatrix<core::Real> R;
			for ( Size ctr=1; ctr<=ars.nrots(); ++ctr ) {
				ars.get(ctr, R);

				// rotated OH vectors
				core::Vector newOH1 = R * OH1c;
				core::Vector newOH2 = R * OH2c;

				// make new residue
				ResidueOP new_res = ResidueOP( new Residue( *rsd_canonic ) );

				// add <x,y,z> to 'O'; compute H's from this
				new_res->set_xyz("O", O);
				new_res->set_xyz("H1", O+newOH1);
				new_res->set_xyz("H2", O+newOH2);

				if ( wat_bump ) {
					core::PackerEnergy bumpenergy = bump_check( new_res, seqpos, scorefxn, pose, task, packer_neighbor_graph );

					Real max_bump_energy = task.max_rotbump_energy();
					if ( bumpenergy < max_bump_energy ) {
						new_rotamers.push_back( new_res );
					}
				} else {
					new_rotamers.push_back( new_res );
				}
			}
		}
	} else {
		// if axis_rot=false or no backbone O atoms found within 3.5 A of the water position
		numeric::UniformRotationSampler urs( rotstep );

		// optionally sample shifts
		//int box[9][3] = {{-1,-1,-1},{-1,-1,1},{-1,1,1},{1,1,1},{1,-1,-1},{1,1,-1},{1,-1,1},{-1,1,-1},{0,0,0}};
		//int box[7][3] = {{0,-1,0},{-1,0,0},{0,0,-1},{0,1,0},{1,0,0},{0,0,1},{0,0,0}};
		//int box[5][3] = {{1,1,1},{-1,-1,1},{-1,1,-1},{1,-1,-1},{0,0,0}};
		//int box[1][3] = {{0,0,0}};

		// rotate so H's are symm about [1,0,0]
		core::Real cos_rot_angle = std::sqrt( 0.5* ( 1+ (OH1c.dot(OH2c)/(OH1c.length()*OH2c.length())) ) );
		core::Real sin_rot_angle = std::sqrt( 1-cos_rot_angle*cos_rot_angle );
		numeric::xyzMatrix<Real> rotH = numeric::xyzMatrix<Real>::rows(
			cos_rot_angle,sin_rot_angle,0,   -sin_rot_angle,cos_rot_angle,0,   0,0,1);
		OH1c = rotH*OH1c;
		OH2c = rotH*OH2c;

		// system is now invariant in rotation of 180 about X.  Eliminate 1/2 of rotations
		numeric::xyzMatrix<Real> flipZ = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1);
		urs.remove_redundant( flipZ );

		numeric::xyzMatrix<core::Real> R;
		for ( core::Size ctr=1; ctr<=urs.nrots() ; ++ctr ) {
			urs.get(ctr, R);

			// rotated OH vectors
			core::Vector newOH1 = R * OH1c;
			core::Vector newOH2 = R * OH2c;

			// make new residue
			ResidueOP new_res = ResidueOP( new Residue( *rsd_canonic ) );

			// add <x,y,z> to 'O'; compute H's from this
			new_res->set_xyz("O", O);
			new_res->set_xyz("H1", O+newOH1);
			new_res->set_xyz("H2", O+newOH2);

			if ( wat_bump ) {
				core::PackerEnergy bumpenergy = bump_check( new_res, seqpos, scorefxn, pose, task, packer_neighbor_graph );

				Real max_bump_energy = task.max_rotbump_energy();
				if ( bumpenergy < max_bump_energy ) {
					new_rotamers.push_back( new_res );
				}
			} else {
				new_rotamers.push_back( new_res );
			}
		}
	}

	// allow water to "virtualize"
	if ( incl_vrt ) {
		ResidueOP vrt_wat = conformation::ResidueFactory::create_residue( residue_set->name_map("HOH_V") );
		vrt_wat->set_xyz("O", O);
		vrt_wat->set_xyz("H1", O+OH1c);
		vrt_wat->set_xyz("H2", O+OH2c);
		new_rotamers.push_back( vrt_wat );
	}

}


// rpav
void
build_backbone_point_water_rotamers(
	Size const seqpos,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	utility::graph::GraphCOP , //packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers,
	bool incl_vrt
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using conformation::Residue;
	using conformation::ResidueOP;
	using numeric::conversions::radians;

	static WaterRotsDB water_rots;
	if ( !water_rots.is_initialized() ) {
		water_rots.initialize();
	}

	bool exclude_exposed = option[corrections::water::exclude_exposed]();

	// all atoms to clash check against
	core::Real clash = 1.0;
	core::Size nclash = 0;
	core::Size nexposed = 0;
	core::Size nres = pose.total_residue();
	utility::vector1< core::Vector > clashcheck;
	utility::vector1< core::Vector > neighbor_atoms; // CB for all except GLY which uses CA
	for ( int i=1; i<=(int)nres; ++i ) {
		core::conformation::Residue const &res = pose.residue(i);
		if ( res.is_protein() ) {
			// get all atoms for buriedness check
			if ( res.name3() != "GLY" ) {
				neighbor_atoms.push_back( res.atom("CB").xyz() );
			} else {
				neighbor_atoms.push_back( res.atom("CA").xyz() );
			}
			for ( int j=1; j<=(int)res.last_backbone_atom(); ++j ) {
				clashcheck.push_back( res.atom(j).xyz() );
			}

			// if gly is not a target type, add CB to the clash list
			bool gly_allowed = false;
			core::pack::task::ResidueLevelTask const &task_i = task.residue_task(i);
			for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator type( task_i.allowed_residue_types_begin() );
					type != task_i.allowed_residue_types_end() && !gly_allowed; ++type ) {
				gly_allowed = ((**type).aa() == core::chemical::aa_gly);
			}
			if ( !gly_allowed && res.aa() != core::chemical::aa_gly ) {
				clashcheck.push_back( res.atom("CB").xyz() );
			}
		} else {
			for ( int j=1; j<=(int)res.nheavyatoms(); ++j ) {
				clashcheck.push_back( res.atom(j).xyz() );
			}
		}
	}

	chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

	// the water "center"
	Residue const & currres( pose.residue(seqpos) );
	core::Vector const O( currres.xyz("O") );

	// find centers of all rotamers to build
	utility::vector1< core::Vector > waters;

	// determine which residue we are building off of and get base atoms for that residue
	// this could definitiely be made more efficient ...
	for ( int i=1; i<=(int)nres; ++i ) {
		core::conformation::Residue const &res_i = pose.residue(i);

		// split protein/ligand logic
		if ( res_i.is_protein() ) {

			// figure out what restypes are allowed at this pos
			utility::vector1< bool > is_aa_allowed ((core::Size) core::chemical::num_canonical_aas, false);
			/* rpav -- commented out for now -- is not used and was causing problems when solvating non-canonicals
			if ( !task.pack_residue(i) && !task.design_residue(i) )  {
			is_aa_allowed[ (core::Size) core::chemical::aa_from_name( res_i.name3() ) ];  // not packing, just solvate current rot
			} else {
			core::pack::task::ResidueLevelTask const &task_i = task.residue_task(i);
			for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator type( task_i.allowed_residue_types_begin() );
			type != task_i.allowed_residue_types_end(); ++type ) {
			is_aa_allowed[(**type).aa()] = true;  // all allowed residue types
			}
			}
			*/
			core::chemical::AtomIndices const & bbatoms = res_i.all_bb_atoms();
			// need to iterate over CB atoms in addtion to BB to allow for SC solvation
			core::chemical::AtomIndices bbatoms_2;
			for ( int j=1; j<=(int)bbatoms.size(); ++j ) {
				bbatoms_2.push_back(bbatoms[j]);
			}
			if ( res_i.name3() != "GLY" ) {
				bbatoms_2.push_back(res_i.atom_index("CB"));
			}

			for ( int j=1; j<=(int)bbatoms_2.size(); ++j ) {
				if (  ( res_i.xyz( bbatoms_2[j] )-O ).length() > 1e-6 ) continue;
				// we have a water/atom correspondence.... now find all matching rotamers
				for ( int k=1; k<=(int)water_rots.n_protein_rots(); ++k ) {
					if ( water_rots.protein_rot(k).aa_ != "*" && core::chemical::aa_from_name(water_rots.protein_rot(k).aa_) != res_i.aa() ) continue;
					if ( strip_whitespace(res_i.atom_name(bbatoms_2[j])) != water_rots.protein_rot(k).aatm_ ) continue;

					// add clash check before adding rotamers

					// add the rotamer!
					core::kinematics::Stub stub( res_i.xyz(water_rots.protein_rot(k).aatm_), res_i.xyz(water_rots.protein_rot(k).abase1_), res_i.xyz(water_rots.protein_rot(k).abase2_) );
					waters.push_back( stub.spherical( radians( water_rots.protein_rot(k).coords_[2] ), radians( 180.0 - water_rots.protein_rot(k).coords_[1] ), water_rots.protein_rot(k).coords_[0] ) );
				}
			}

		} else {
			// add PWAT rotamers for ligands
			// donors
			for ( chemical::AtomIndices::const_iterator
					hnum  = res_i.Hpos_polar().begin(),
					hnume = res_i.Hpos_polar().end(); hnum != hnume; ++hnum ) {
				if (  ( res_i.xyz( *hnum )-O ).length() > 1e-6 ) continue;

				core::Size hatm(*hnum);
				core::Size datm( res_i.atom_base( hatm ) );
				core::Size datm_base( res_i.atom_base( datm ) );
				if ( datm_base == hatm ) { // could happen with water
					datm_base = 0;
					chemical::AtomIndices const & datm_nbrs( res_i.type().nbrs( datm ) );
					for ( Size ii=1; ii<= datm_nbrs.size(); ++ii ) {
						if ( datm_nbrs[ii] == hatm ) continue;
						else datm_base = datm_nbrs[ii];
					}
					runtime_assert( datm_base );
				}

				// we have a water/atom correspondence.... now find all matching rotamers
				for ( int k=1; k<=(int)water_rots.n_ligand_rots(); ++k ) {
					if ( water_rots.ligand_rot(k).aa_ != "DON" ) continue;

					// add the rotamer!
					core::kinematics::Stub stub( res_i.xyz(hatm), res_i.xyz(datm), res_i.xyz(datm_base) );
					core::Real hd_len((res_i.xyz(hatm)-res_i.xyz(datm)).length()); // subtract bond length from pwat placement distance
					waters.push_back( stub.spherical( radians( water_rots.ligand_rot(k).coords_[2] ), radians( 180.0 - water_rots.ligand_rot(k).coords_[1] ), water_rots.ligand_rot(k).coords_[0]-hd_len ) );
				}
			}

			// acceptors
			for ( chemical::AtomIndices::const_iterator
					anum  = res_i.accpt_pos().begin(),
					anume = res_i.accpt_pos().end(); anum != anume; ++anum ) {
				if (  ( res_i.xyz( *anum )-O ).length() > 1e-6 ) continue;

				core::Size const aatm( *anum );
				core::Size const abase1( res_i.atom_base( aatm ) ), abase2( res_i.abase2( aatm ) );

				// get hybridization
				core::chemical::Hybridization hyb = res_i.atom_type( *anum ).hybridization();
				std::string tgt_name;
				if ( hyb == core::chemical::SP2_HYBRID ) tgt_name = "SP2_ACC";
				if ( hyb == core::chemical::SP3_HYBRID ) tgt_name = "SP3_ACC";
				if ( hyb == core::chemical::RING_HYBRID ) tgt_name = "RING_ACC";


				// we have a water/atom correspondence.... now find all matching rotamers
				for ( int k=1; k<=(int)water_rots.n_ligand_rots(); ++k ) {
					if ( water_rots.ligand_rot(k).aa_ != tgt_name ) continue;

					// add the rotamer!
					if ( hyb == core::chemical::RING_HYBRID ) {
						numeric::xyzVector< core::Real > virtual_abase1 = 0.5 * ( res_i.xyz(abase1) + res_i.xyz(abase2) );
						core::kinematics::Stub stub( res_i.xyz(aatm), virtual_abase1, res_i.xyz(abase2) );
						waters.push_back( stub.spherical( radians( water_rots.ligand_rot(k).coords_[2] ), radians( 180.0 - water_rots.ligand_rot(k).coords_[1] ), water_rots.ligand_rot(k).coords_[0] ) );
					} else {
						core::kinematics::Stub stub( res_i.xyz(aatm), res_i.xyz(abase1), res_i.xyz(abase2) );
						waters.push_back( stub.spherical( radians( water_rots.ligand_rot(k).coords_[2] ), radians( 180.0 - water_rots.ligand_rot(k).coords_[1] ), water_rots.ligand_rot(k).coords_[0] ) );
					}
				}
			}
		}
	}

	// now build the rotamers
	ResidueOP rsd_canonic = conformation::ResidueFactory::create_residue( residue_set->name_map("PWAT") );
	core::Vector OH1c( rsd_canonic->xyz("V1") - rsd_canonic->xyz("O") );
	core::Vector OH2c( rsd_canonic->xyz("V2") - rsd_canonic->xyz("O") );
	for ( Size i=1; i<=waters.size(); ++i ) {
		ResidueOP new_res = ResidueOP( new Residue( *rsd_canonic ) );
		new_res->set_xyz("O", waters[i]);
		new_res->set_xyz("V1", waters[i]+OH1c);
		new_res->set_xyz("V2", waters[i]+OH2c);

		// clash check with backbone before adding rotamer
		bool isclash = false;
		for ( int ii=1; ii<=(int)clashcheck.size() && !isclash; ++ii ) {
			core::Real dist = (new_res->xyz("O")-clashcheck[ii]).length_squared();
			isclash = (dist < clash*clash);
		}
		if ( isclash ) {
			nclash++;
			continue;
		}

		if ( exclude_exposed ) {
			bool isexposed = true;
			Size nneighbors = 0;
			for ( Size jj=1; jj <= neighbor_atoms.size() && isexposed; ++jj ) {
				core::Real dist = new_res->xyz("O").distance(neighbor_atoms[jj]);
				if ( dist <= 10.0 ) ++nneighbors;
				if ( nneighbors > 16 ) isexposed = false;
			}
			if ( isexposed ) {
				++nexposed;
				continue;
			}
		}

		new_rotamers.push_back( new_res );
	}

	// allow water to "virtualize"
	if ( incl_vrt ) {
		ResidueOP vrt_wat = conformation::ResidueFactory::create_residue( residue_set->name_map("PWAT_V") );
		vrt_wat->set_xyz("O", O);
		vrt_wat->set_xyz("V1", O+OH1c);
		vrt_wat->set_xyz("V2", O+OH2c);
		new_rotamers.push_back( vrt_wat );
	}
	if ( ( exclude_exposed ) && ( nexposed > 0 ) ) {
		TR << "Removed " << nexposed << " point water rotamers due to surface exposure for residue " << seqpos << std::endl;
	}
}




} // rotamer_set
} // pack
} // core
