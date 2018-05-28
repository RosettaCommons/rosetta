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
	using namespace scoring;
	using namespace conformation;
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

	// collected all backbone O atoms within 3.5 A of this water position
	// fd: perhaps we should so the same around ligands?
	utility::vector1< core::Vector > COs;
	if ( sample_axis ) {
		EnergyMap emap;

		for ( utility::graph::Graph::EdgeListConstIter
				ir  = packer_neighbor_graph->get_node( seqpos )->const_edge_list_begin(),
				ire = packer_neighbor_graph->get_node( seqpos )->const_edge_list_end();
				ir != ire; ++ir ) {
			int const neighbor_id( (*ir)->get_other_ind( seqpos ) );
			Residue const & neighbor( pose.residue( neighbor_id ) );
			Vector neighborO =  neighbor.xyz("O");
			Real dist = O.distance( neighborO );
			if ( dist <= CO_cutoff ) {
				COs.push_back( neighborO ); // save all incase there are > 1
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

			Vector init_OH1 = ( H1-O );
			Vector init_OH2 = ( H2-O );

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

} // rotamer_set
} // pack
} // core
