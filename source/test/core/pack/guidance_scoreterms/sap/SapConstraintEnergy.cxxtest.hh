// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/ApproximateBuriedUnsatPenalty.cxxtest.hh
/// @brief  test suite for ApproximateBuriedUnsatPenalty
/// @author Brian Coventry (bcov@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/pack/guidance_scoreterms/sap/SapConstraintEnergy.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapDatabase.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh>
#include <core/pack/guidance_scoreterms/sap/SapScoreMetric.hh>
#include <core/pack/guidance_scoreterms/sap/PerResidueSapScoreMetric.hh>
#include <core/pack/guidance_scoreterms/sap/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/init_id_map.hh>
#include <core/import_pose/import_pose.hh>

#include <util/pose_funcs.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/AddSapConstraintMover.hh>
#include <protocols/simple_moves/AddSapMathConstraintMover.hh>

#include <core/types.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>

#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>
#include <utility/io/izstream.hh>

#include <boost/format.hpp>

#include <basic/Tracer.hh> // AUTO IWYU For Tracer
#include <numeric/xyzMatrix.hh> // AUTO IWYU For xyzMatrix
#include <core/pack/task/PackerTask.hh> // AUTO IWYU For PackerTask

using namespace core;
using namespace core::pack::guidance_scoreterms::sap;


static basic::Tracer TR("core.pack.guidance_scoreterms.sap.SapConstraintEnergyTests.cxxtest");

#define SYMM 3


struct TestBlockParam {
	Real max_sasa_score;
	Real no_block;
	Real full_block;

	TestBlockParam( Real max, Real no, Real full )
	:
		max_sasa_score( max ),
		no_block( no ),
		full_block( full )
	{}

};

class SapConstraintEnergyTests : public CxxTest::TestSuite
{

	bool old_symm_ = false;

public:
	void setUp()
	{
		core_init();
	}

	utility::pointer::shared_ptr< std::unordered_map< std::string, TestBlockParam > >
	load_block_data( ) {

		utility::pointer::shared_ptr< std::unordered_map< std::string, TestBlockParam > > atomtype_to_block_param
			= utility::pointer::make_shared< std::unordered_map< std::string, TestBlockParam > >();

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/sap_sasa_calib.dat" );

		std::string tmp;
		stream >> tmp;
		runtime_assert( tmp == "atom_type" );
		stream >> tmp;
		runtime_assert( tmp == "full_block" );
		stream >> tmp;
		runtime_assert( tmp == "max_sasa" );
		stream >> tmp;
		runtime_assert( tmp == "no_block" );

		std::string atom_type;
		Real full_block;
		Real max_sasa;
		Real no_block;

		while ( !stream.eof() ) {

			stream >> atom_type;
			stream >> full_block;
			stream >> max_sasa;
			stream >> no_block;

			atomtype_to_block_param->emplace( std::make_pair( atom_type, TestBlockParam( max_sasa, no_block, full_block ) ));

		}
		return atomtype_to_block_param;
	}


	utility::vector1<utility::vector1<Real>>
	manually_calculate_fast_blocks(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets
	) {
		utility::vector1<utility::vector1<Real>> blocks( rotamer_sets.nrotamers() + pose.size() );    // rotamer set then pose

		utility::vector1<bool> designing_residues = rotamer_sets.task()->designing_residues();

		for ( Size seqpos1 = 1; seqpos1 <= pose.size(); seqpos1++ ) {
			// for res1 we use the pose res
			core::conformation::Residue const & res1 = pose.residue(seqpos1);
			if ( res1.is_virtual_residue() ) continue;

			for ( Size seqpos2 = 1; seqpos2 <= pose.size(); seqpos2++ ) {

				if ( seqpos1 == seqpos2 ) continue;
				if ( pose.residue( seqpos2).is_virtual_residue() ) continue;

				// for res2, we iterate through the rotamer set
				Size begin = 0;
				Size number = 0;

				if ( rotamer_sets.has_rotamer_set_for_residue( seqpos2 ) ) {
					begin = rotamer_sets.nrotamer_offset_for_moltenres( rotamer_sets.resid_2_moltenres( seqpos2 ) ) + 1;
					number = rotamer_sets.nrotamers_for_moltenres( rotamer_sets.resid_2_moltenres( seqpos2 ) );
				}

				for ( Size ilocal = 0; ilocal <= number; ilocal++ ) {
					Size irot = ilocal + begin - 1;
					if ( ilocal == 0 ) irot = seqpos2 + rotamer_sets.nrotamers();

					core::conformation::Residue const & res2 = ilocal == 0 ? pose.residue( seqpos2 ) : *rotamer_sets.rotamer( irot );

					if ( blocks[irot].size() == 0 ) blocks[irot].resize( res2.natoms() );


					for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
						Real lj1 = res1.atom_type(at1).lj_radius();

						// This should mirror the code in helper
						if ( designing_residues[seqpos1] && ! ( (int)at1 - (int)res1.first_sidechain_atom() <= 0 ) ) continue;

						for ( Size at2 = 1; at2 <= res2.natoms(); at2++ ) {
							Real lj2 = res2.atom_type(at2).lj_radius();

							Real full_block = lj1 + lj2;
							Real full_unblock = lj1 + lj2 + 3;

							Real dist = res1.xyz(at1).distance(res2.xyz(at2));

							Real frac = ( full_unblock - dist ) / ( full_unblock - full_block );
							frac = std::min<Real>( frac, 1 );
							frac = std::max<Real>( frac, 0 );

							blocks[irot][at2] += 5 * lj2*lj2 * frac;
						}
					}
				}
			}
		}
		return blocks;
	}

	utility::vector1<utility::vector1<float>>
	manually_calculate_fast_sasa_scores(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & rotamer_sets,
		utility::vector1<utility::vector1<Real>> const & blocks
	) {
		utility::vector1<utility::vector1<float>> sasas( rotamer_sets.nrotamers() + pose.size() );    // rotamer set then pose

		SapDatabase * db = SapDatabase::get_instance();

		utility::pointer::shared_ptr< std::unordered_map< std::string, TestBlockParam > > atomtype_to_block_param
			= load_block_data( );

		for ( Size irot = 1; irot <= blocks.size(); irot++ ) {

			core::conformation::Residue const & res1 = irot <= rotamer_sets.nrotamers() ? *rotamer_sets.rotamer( irot )
				: pose.residue( irot - rotamer_sets.nrotamers() );
			if ( res1.is_virtual_residue() ) continue;

			sasas[irot].resize( res1.natoms() );

			float max_sasa = 0;
			for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
				if ( res1.atom_is_backbone(at1) ) continue;

				std::string atom_type = res1.name3() + "_" + utility::trim( res1.atom_name(at1) );
				if ( atomtype_to_block_param->count( atom_type ) == 0 ) {
					TR << "Missing atom type: " << atom_type << std::endl;
					continue;
				}
				TestBlockParam const & param = atomtype_to_block_param->at( atom_type );

				max_sasa += param.max_sasa_score;
			}

			float hydrophobic = db->hydrophobic_weight(res1.name1());

			for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
				if ( res1.atom_is_backbone(at1) ) continue;

				std::string atom_type = res1.name3() + "_" + utility::trim( res1.atom_name(at1) );
				if ( atomtype_to_block_param->count( atom_type ) == 0 ) {
					continue;
				}
				TestBlockParam const & param = atomtype_to_block_param->at( atom_type );

				Real frac = ( param.full_block - blocks[irot][at1] ) / ( param.full_block - param.no_block );
				frac = std::min<Real>( frac, 1 );
				frac = std::max<Real>( frac, 0 );

				sasas[ irot ][ at1 ] = frac * param.max_sasa_score / max_sasa * hydrophobic;
			}
		}
		return sasas;
	}



	// There aren't that many atoms in trpcage...
	core::id::AtomID_Map<Real>
	manually_calculate_blocks( core::pose::Pose const & pose ) {
		core::id::AtomID_Map<Real> atoms;
		core::pose::initialize_atomid_map( atoms, pose, Real(0) );


		for ( Size seqpos1 = 1; seqpos1 <= pose.size(); seqpos1++ ) {
			core::conformation::Residue const & res1 = pose.residue(seqpos1);
			if ( res1.is_virtual_residue() ) continue;

			for ( Size seqpos2 = 1; seqpos2 <= pose.size(); seqpos2++ ) {
				core::conformation::Residue const & res2 = pose.residue(seqpos2);
				if ( res2.is_virtual_residue() ) continue;

				if ( seqpos1 == seqpos2 ) continue;

				for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
					Real lj1 = res1.atom_type(at1).lj_radius();

					for ( Size at2 = 1; at2 <= res2.natoms(); at2++ ) {
						Real lj2 = res2.atom_type(at2).lj_radius();

						Real full_block = lj1 + lj2;
						Real full_unblock = lj1 + lj2 + 3;

						Real dist = res1.xyz(at1).distance(res2.xyz(at2));

						Real frac = ( full_unblock - dist ) / ( full_unblock - full_block );
						frac = std::min<Real>( frac, 1 );
						frac = std::max<Real>( frac, 0 );


						atoms(seqpos1, at1) += 5 * lj1*lj1 * frac / 2; // divide by 2 because we are double counting
						atoms(seqpos2, at2) += 5 * lj2*lj2 * frac / 2;

					}
				}
			}
		}
		return atoms;
	}

	core::id::AtomID_Map<Real>
	manually_calculate_sasa_scores(
		core::pose::Pose const & pose,
		core::id::AtomID_Map<Real> const & blocks
	) {
		core::id::AtomID_Map<Real> atoms;
		core::pose::initialize_atomid_map( atoms, pose, Real(0) );

		SapDatabase * db = SapDatabase::get_instance();

		utility::pointer::shared_ptr< std::unordered_map< std::string, TestBlockParam > > atomtype_to_block_param
			= load_block_data( );

		for ( Size seqpos1 = 1; seqpos1 <= pose.size(); seqpos1++ ) {
			core::conformation::Residue const & res1 = pose.residue(seqpos1);
			if ( res1.is_virtual_residue( ) ) continue;

			float max_sasa = 0;
			for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
				if ( res1.atom_is_backbone(at1) ) continue;

				std::string atom_type = res1.name3() + "_" + utility::trim( res1.atom_name(at1) );
				if ( atomtype_to_block_param->count( atom_type ) == 0 ) {
					TR << "Missing atom type: " << atom_type << std::endl;
					continue;
				}
				TestBlockParam const & param = atomtype_to_block_param->at( atom_type );

				max_sasa += param.max_sasa_score;
			}

			float hydrophobic = db->hydrophobic_weight(res1.name1());

			for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
				if ( res1.atom_is_backbone(at1) ) continue;

				std::string atom_type = res1.name3() + "_" + utility::trim( res1.atom_name(at1) );
				if ( atomtype_to_block_param->count( atom_type ) == 0 ) {
					continue;
				}
				TestBlockParam const & param = atomtype_to_block_param->at( atom_type );

				Real frac = ( param.full_block - blocks(seqpos1, at1) ) / ( param.full_block - param.no_block );
				frac = std::min<Real>( frac, 1 );
				frac = std::max<Real>( frac, 0 );


				atoms( seqpos1, at1 ) = frac * param.max_sasa_score / max_sasa * hydrophobic;
			}
		}
		return atoms;
	}

	core::id::AtomID_Map<Real>
	manually_calculate_sap_scores(
		core::pose::Pose const & pose,
		core::id::AtomID_Map<Real> const & sasa_scores,
		Real & sap_score
	) {
		core::id::AtomID_Map<Real> atoms;
		core::pose::initialize_atomid_map( atoms, pose, Real(0) );


		for ( Size seqpos1 = 1; seqpos1 <= pose.size(); seqpos1++ ) {
			core::conformation::Residue const & res1 = pose.residue(seqpos1);
			if ( res1.is_virtual_residue() ) continue;

			for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
				if ( res1.atom_is_backbone(at1) ) continue;

				for ( Size seqpos2 = 1; seqpos2 <= pose.size(); seqpos2++ ) {
					core::conformation::Residue const & res2 = pose.residue(seqpos2);
					if ( res2.is_virtual_residue() ) continue;

					for ( Size at2 = 1; at2 <= res2.natoms(); at2++ ) {
						if ( res2.atom_is_backbone(at2) ) continue;

						Real dist2 = res1.xyz(at1).distance_squared(res2.xyz(at2));

						if ( dist2 > 5*5 ) continue;

						atoms(seqpos1, at1) += sasa_scores(seqpos2, at2);
					}
				}
				if ( atoms(seqpos1, at1) > 0 ) {
					sap_score += atoms(seqpos1, at1);
				}
			}
		}
		return atoms;
	}


	void assert_blocks(
		utility::vector1< core::conformation::ResidueCOP > const & res_vector,
		SapConstraintHelperOP const & helper
	) {

		for ( Size seqpos1 = 1; seqpos1 <= res_vector.size(); seqpos1++ ) {
			core::conformation::Residue const & res1 = *res_vector[seqpos1];
			if ( res1.is_virtual_residue() ) continue;

			for ( Size seqpos2 = 1; seqpos2 <= res_vector.size(); seqpos2++ ) {
				core::conformation::Residue const & res2 = *res_vector[seqpos2];
				if ( res2.is_virtual_residue() ) continue;

				if ( seqpos1 >= seqpos2 ) continue;

				utility::vector1<Real> blocks_on_1( res1.natoms(), 0 );
				utility::vector1<Real> blocks_on_2( res2.natoms(), 0 );

				bool any_blocks = false;

				for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
					Real lj1 = res1.atom_type(at1).lj_radius();

					for ( Size at2 = 1; at2 <= res2.natoms(); at2++ ) {
						Real lj2 = res2.atom_type(at2).lj_radius();

						Real full_block = lj1 + lj2;
						Real full_unblock = lj1 + lj2 + 3;

						Real dist = res1.xyz(at1).distance(res2.xyz(at2));

						Real frac = ( full_unblock - dist ) / ( full_unblock - full_block );
						frac = std::min<Real>( frac, 1 );
						frac = std::max<Real>( frac, 0 );

						blocks_on_1[at1] += 5 * lj1*lj1 * frac;
						blocks_on_2[at2] += 5 * lj2*lj2 * frac;

						any_blocks |= frac > 0 && ! (res1.atom_is_backbone(at1) && res2.atom_is_backbone(at2));

					}
				}

				Size first_sidechain1 = res1.first_sidechain_atom();
				Size natoms1 = res1.natoms();
				Size first_sidechain2 = res2.first_sidechain_atom();
				Size natoms2 = res2.natoms();

				std::pair< core::conformation::Residue const *, core::conformation::Residue const * > key(
					&*res_vector[seqpos1], &*res_vector[seqpos2] );

				if ( ! any_blocks ) {
					TS_ASSERT( helper->interacting_block_offset_.count( key ) == 0 );
				} else {
					TS_ASSERT( helper->interacting_block_offset_.count( key ) > 0 );
					if ( helper->interacting_block_offset_.count( key ) == 0 ) {
						std::cout << "MISS: " << seqpos1 << " " << seqpos2 << std::endl;
					}
				}
				if ( helper->interacting_block_offset_.count( key ) == 0 ) continue;

				Size offset = helper->interacting_block_offset_.at( key );

				for ( Size iat = 0; (int)iat <= (int)natoms1 - (int)first_sidechain1; iat++ ) {
					Size iatom = iat + first_sidechain1;
					uint8_t expected = (uint8_t)std::min<Size>( 255, blocks_on_1[iatom] / SAP_BLOCK_STORE_SCALE );

					TS_ASSERT_DELTA( expected, helper->interacting_block_[offset++], 2 );

				}

				for ( Size iat = 0; (int)iat <= (int)natoms2 - (int)first_sidechain2; iat++ ) {
					Size iatom = iat + first_sidechain2;
					uint8_t expected = (uint8_t)std::min<Size>( 255, blocks_on_2[iatom] / SAP_BLOCK_STORE_SCALE );

					TS_ASSERT_DELTA( expected, helper->interacting_block_[offset++], 2 );

				}

			}
		}
	}

	struct WeirdIter {
		Size offset;
		utility::vector1<uint8_t> const & atom_within_5;
		uint8_t const * value_data;

		Size cur_pos;

		WeirdIter(
			Size _offset,
			utility::vector1<uint8_t> const & _atom_within_5,
			uint8_t const _value_data[ATOM_WITHIN_5_ELEMS],
			Size pos
		) :
			offset( _offset ),
			atom_within_5( _atom_within_5 ),
			value_data( _value_data ),
			cur_pos( pos )
		{}

		uint8_t
		operator *() {
			if ( cur_pos < ATOM_WITHIN_5_ELEMS ) return value_data[ cur_pos ];
			else return atom_within_5[ offset + cur_pos - ATOM_WITHIN_5_ELEMS ];
		}

		WeirdIter &
		operator++(int) {
			cur_pos++;
			return *this;
		}
	};

	void
	assert_nearest_5(
		utility::vector1< core::conformation::ResidueCOP > const & res_vector,
		SapConstraintHelperOP const & helper
	) {


		for ( Size seqpos1 = 1; seqpos1 <= res_vector.size(); seqpos1++ ) {
			core::conformation::Residue const & res1 = *res_vector[seqpos1];
			if ( res1.is_virtual_residue() ) continue;

			for ( Size seqpos2 = 1; seqpos2 <= res_vector.size(); seqpos2++ ) {
				core::conformation::Residue const & res2 = *res_vector[seqpos2];
				if ( res2.is_virtual_residue() ) continue;

				utility::vector1<utility::vector1<bool>> is_close( res1.natoms(), utility::vector1<bool>( res2.natoms(), false ) );
				utility::vector1<utility::vector1<bool>> is_close_helper( res1.natoms(), utility::vector1<bool>( res2.natoms(), false ) );

				bool any_close = false;

				for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
					if ( res1.atom_is_backbone(at1) ) continue;


					for ( Size at2 = 1; at2 <= res2.natoms(); at2++ ) {
						if ( res2.atom_is_backbone(at2) ) continue;

						Real dist2 = res1.xyz(at1).distance_squared(res2.xyz(at2));

						if ( dist2 > 5*5 ) continue;

						is_close[at1][at2] = true;
						any_close = true;
					}
				}

				std::pair< core::conformation::Residue const *, core::conformation::Residue const * > key(
					&*res_vector[seqpos1], &*res_vector[seqpos2] );

				if ( any_close ) {
					TS_ASSERT( helper->atom_within_5_map_.count( key ) != 0 );
				} else {
					TS_ASSERT( helper->atom_within_5_map_.count( key ) == 0 );
				}

				if ( helper->atom_within_5_map_.count( key ) == 0 ) {
					continue;
				}

				Size first_sidechain1 = res1.first_sidechain_atom();
				Size natoms1 = res1.natoms();
				Size first_sidechain2 = res2.first_sidechain_atom();
				// Size natoms2 = res2.natoms();


				typename SapConstraintHelper::atom_within_5_value const & value = helper->atom_within_5_map_.at( key );

				WeirdIter iter( value.first, helper->atom_within_5_, value.second, natoms1-first_sidechain1+1 );


				for ( Size iat = 0; (int)iat <= (int)natoms1 - (int)first_sidechain1; iat++ ) {
					Size iatom = iat + first_sidechain1;
					Size elems = value.second[iat];
					for ( Size ielem = 0; ielem < elems; ielem++ ) {
						Size iatom2 = *iter + first_sidechain2;
						iter++;
						is_close_helper[iatom][iatom2] = true;
					}
				}

				for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
					if ( res1.atom_is_backbone(at1) ) continue;
					for ( Size at2 = 1; at2 <= res2.natoms(); at2++ ) {
						if ( res2.atom_is_backbone(at2) ) continue;

						TS_ASSERT_EQUALS((int)is_close[at1][at2], (int)is_close_helper[at1][at2]);
					}
				}
			}
		}
	}


	void
	assert_sasa_scores(
		utility::vector1< core::conformation::ResidueCOP > const & res_vector,
		SapConstraintHelperOP const & helper
	) {

		SapDatabase * db = SapDatabase::get_instance();

		utility::pointer::shared_ptr< std::unordered_map< std::string, TestBlockParam > > atomtype_to_block_param
			= load_block_data( );

		for ( Size seqpos1 = 1; seqpos1 <= res_vector.size(); seqpos1++ ) {
			core::conformation::Residue const & res1 = *res_vector[seqpos1];
			if ( res1.is_virtual_residue() ) continue;

			float max_sasa = 0;
			for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
				if ( res1.atom_is_backbone(at1) ) continue;

				std::string atom_type = res1.name3() + "_" + utility::trim( res1.atom_name(at1) );
				if ( atomtype_to_block_param->count( atom_type ) == 0 ) {
					TR << "Missing atom type: " << atom_type << std::endl;
					continue;
				}
				TestBlockParam const & param = atomtype_to_block_param->at( atom_type );

				max_sasa += param.max_sasa_score;
			}

			float hydrophobic = db->hydrophobic_weight(res1.name1());

			for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
				if ( res1.atom_is_backbone(at1) ) continue;
				Size iat = at1 - res1.first_sidechain_atom();

				std::string atom_type = res1.name3() + "_" + utility::trim( res1.atom_name(at1) );
				if ( atomtype_to_block_param->count( atom_type ) == 0 ) {
					continue;
				}
				TestBlockParam const & param = atomtype_to_block_param->at( atom_type );

				Real frac = Real( int16_t(param.full_block / SAP_BLOCK_STORE_SCALE ) - helper->sasa_blocks_[seqpos1][ iat ] )
					/ Real( int16_t( param.full_block / SAP_BLOCK_STORE_SCALE ) - int16_t( param.no_block / SAP_BLOCK_STORE_SCALE ) );
				frac = std::min<Real>( frac, 1 );
				frac = std::max<Real>( frac, 0 );


				Real sasa_score = frac * param.max_sasa_score / max_sasa * hydrophobic;


				TS_ASSERT_DELTA( sasa_score, helper->atom_sasa_score_[seqpos1][iat], 0.001 );

				// if ( std::abs( sasa_score - helper->atom_sasa_score_[seqpos1][iat]) > 0.001 ) {

				//     std::cout << seqpos1 << "-" << at1 << " " << param.max_sasa_score / max_sasa * hydrophobic << " "
				//               << +int16_t( param.full_block / SAP_BLOCK_STORE_SCALE ) << " "
				//               << +int16_t( param.no_block / SAP_BLOCK_STORE_SCALE ) << " " << sasa_score << " "
				//               << +helper->sasa_blocks_[seqpos1][ iat ] << std::endl;
				// }

			}
		}
	}


	void
	assert_saps(
		core::pose::Pose const & pose,
		utility::vector1< core::conformation::ResidueCOP > const & res_vector,
		SapConstraintHelperOP const & helper,
		bool fast = false
	) {

		core::id::AtomID_Map<Real> atoms;
		core::pose::initialize_atomid_map( atoms, pose, Real(0) );

		core::id::AtomID_Map<bool> was_used;
		core::pose::initialize_atomid_map( was_used, pose, false );

		Real sap_score = 0;

		Real manual_sap_sum = 0;

		for ( Size seqpos1 = 1; seqpos1 <= res_vector.size(); seqpos1++ ) {
			core::conformation::Residue const & res1 = *res_vector[seqpos1];
			if ( res1.is_virtual_residue() ) continue;

			for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
				if ( res1.atom_is_backbone(at1) ) continue;

				was_used( seqpos1, at1 ) = true;

				for ( Size seqpos2 = 1; seqpos2 <= res_vector.size(); seqpos2++ ) {
					core::conformation::Residue const & res2 = *res_vector[seqpos2];
					if ( res2.is_virtual_residue() ) continue;

					for ( Size at2 = 1; at2 <= res2.natoms(); at2++ ) {
						if ( res2.atom_is_backbone(at2) ) continue;

						Real dist2 = res1.xyz(at1).distance_squared(res2.xyz(at2));

						if ( dist2 > 5*5 ) continue;

						if ( fast ) {
							float * sasas = helper->atom_sasa_score_fast_[seqpos2];
							assert( sasas );
							atoms(seqpos1, at1) += sasas[at2-res2.first_sidechain_atom()];
						} else {
							atoms(seqpos1, at1) += helper->atom_sasa_score_[seqpos2][at2-res2.first_sidechain_atom()];
						}
					}
				}
				if ( atoms(seqpos1, at1) > 0 ) {
					sap_score += atoms(seqpos1, at1);
				}


				TS_ASSERT_DELTA( atoms(seqpos1, at1), helper->atom_sap_[seqpos1][at1-res1.first_sidechain_atom()], 0.001 );

				if ( helper->atom_sap_[seqpos1][at1-res1.first_sidechain_atom()] > 0 ) {
					manual_sap_sum += helper->atom_sap_[seqpos1][at1-res1.first_sidechain_atom()];
				}

			}
		}
		TS_ASSERT_DELTA( sap_score, helper->current_score_, 0.1 );
		TS_ASSERT_DELTA( sap_score, manual_sap_sum, 0.1 );
		TS_ASSERT_DELTA( manual_sap_sum, helper->current_score_, 0.1 );


		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
			core::conformation::ResidueCOP const & rotamer = res_vector[ seqpos ];
			for ( Size at = 0; at < helper->atom_sap_[seqpos].size(); at++ ) {
				Size iatom = at + rotamer->first_sidechain_atom();
				if ( iatom > rotamer->natoms() || !was_used( seqpos, iatom ) ) {
					TS_ASSERT( helper->atom_sap_[seqpos][at] == 0 );
				}
			}
		}



	}

	// Center the pose with center of mass on origin then apply c3 symmetry for maximum ASU-ASU interaction
	void terrible_C3_symmetry( core::pose::Pose & symm_pose, bool actually_do_it ) {
		if ( ! actually_do_it ) return;

		SapDatabase * db = SapDatabase::get_instance();
		db->symm_debug_ = true;

		numeric::xyzVector< Real > com = core::pose::center_of_mass( symm_pose, 1, symm_pose.size() );
		symm_pose.apply_transform_Rx_plus_v( numeric::xyzMatrix<Real>::identity(), -com );


		std::string symdef_file = basic::database::full_name( "symmetry/cyclic/C3_Z.sym" );
		protocols::symmetry::SetupForSymmetryMover setup( symdef_file, basic::options::option );
		setup.apply( symm_pose );
	}


	// void test_helper(){
	//        TR << "test_helper" << std::endl;
	//        for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
	//            bool using_symm = bool(int_sym);
	//            TR << (using_symm ? "Symm pass" : "first pass") << std::endl;


	//      core::pose::Pose pose = create_trpcage_ideal_pose();

	//            terrible_C3_symmetry( pose, using_symm );

	//      core::select::residue_selector::TrueResidueSelectorOP sel =
	//       utility::pointer::make_shared<core::select::residue_selector::TrueResidueSelector>();


	//      SapConstraintOptionsOP options = utility::pointer::make_shared<SapConstraintOptions>();


	//      utility::vector1< core::conformation::ResidueCOP > res_vector;
	//      SapConstraintHelperOP helper = common_setup( options, pose, sel, sel, sel, res_vector );
	//      Real helper_sap_score = helper->calculate_energy( res_vector, 0 );

	//      Real sap_score = 0;


	//            assert_blocks( res_vector, helper );
	//            assert_nearest_5( res_vector, helper );

	//            // utility::vector1<Size> positions;
	//            // for ( Size i = 1; i <= pose.size(); i++ ) {
	//            //     positions.push_back(i);
	//            //     helper->shadow_mismatch_[i] = true;
	//            //     for ( Size j = 0; j < helper->max_rotamer_atoms_; j++ ) {
	//            //         helper->dirty_sasa_[i][j] = true;
	//            //     }

	//            // }

	//            // helper->recalculate_sasa( );
	//            // helper->recalculate_saps( positions );


	//            assert_sasa_scores( res_vector, helper );
	//            assert_saps( pose, res_vector, helper );


	//      core::id::AtomID_Map<Real> blocks = manually_calculate_blocks( pose );
	//      core::id::AtomID_Map<Real> sasa_scores = manually_calculate_sasa_scores( pose, blocks );
	//      core::id::AtomID_Map<Real> sap_scores = manually_calculate_sap_scores( pose, sasa_scores, sap_score );


	//      TS_ASSERT_DELTA( helper_sap_score, sap_score, 1 );

	//      // for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

	//      //  core::conformation::Residue const & res = pose.residue(seqpos);

	//      //        Size first_sidechain = res.first_sidechain_atom();
	//      //        Size natoms = res.natoms();
	//      //  for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {

	//      //   Size iatom = iat + first_sidechain;

	//      //   std::cout << "T " << seqpos << " " << iatom << " " << blocks(seqpos, iatom) / SAP_BLOCK_STORE_SCALE << " " << sasa_scores(seqpos, iatom)
	//      //       << " " << sap_scores(seqpos, iatom) << std::endl;

	//      //   std::cout << "H " << seqpos << " " << iatom << " " << helper->sasa_blocks_[seqpos][iat] << " " << helper->atom_sasa_score_[seqpos][iat]
	//      //       << " " << helper->atom_sap_[seqpos][iat] << std::endl;

	//      //  }
	//      // }


	//      for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

	//       core::conformation::Residue const & res = pose.residue(seqpos);
	//                if ( res.is_virtual_residue() ) continue;

	//             Size first_sidechain = res.first_sidechain_atom();
	//             Size natoms = res.natoms();
	//       for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
	//        Size iatom = iat + first_sidechain;

	//        uint16_t expected_block = std::min<uint16_t>( 255, blocks(seqpos, iatom) / SAP_BLOCK_STORE_SCALE );
	//        uint16_t clipped_helper = std::min<uint16_t>( 255, helper->sasa_blocks_[seqpos][iat] );

	//        TS_ASSERT_DELTA( expected_block, clipped_helper, 3 );


	//        TS_ASSERT_DELTA( sasa_scores(seqpos, iatom), helper->atom_sasa_score_[seqpos][iat], 0.03 );

	//        TS_ASSERT_DELTA( sap_scores(seqpos, iatom), helper->atom_sap_[seqpos][iat], 0.05 );

	//       }
	//      }

	//        }

	// }

	void prepare_packing_pose_and_rotsets(
		core::pose::Pose & pose,
		utility::vector1< core::conformation::ResidueCOP > & res_vector,
		pack::rotamer_set::RotamerSetsOP & rotsets,
		SapConstraintOptionsOP & options,
		select::residue_selector::ResidueSelectorCOP & sel,
		Real & pack_size,
		bool symmetry,
		bool doubles = false,
		select::residue_selector::ResidueSelectorCOP design_sel = utility::pointer::make_shared<core::select::residue_selector::TrueResidueSelector>()
	) {

		// No cysteine because it crashes when we do symmetry
		std::string letters = "ADEFGHIKLMNPQRSTVWY";

		utility::vector1<core::pose::Pose> poses( letters.size() );
		for ( Size i = 1; i <= poses.size(); i++ ) {
			std::string sequence = std::string( 6, letters[i-1] );
			core::pose::make_pose_from_sequence( poses[i], sequence, "fa_standard" );

			terrible_C3_symmetry( poses[i], symmetry );
		}


		// core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "AAAAAA", "fa_standard" );
		pack_size = pose.size();

		terrible_C3_symmetry( pose, symmetry );


		core::pack::task::operation::OperateOnResidueSubsetOP restrict_task = utility::pointer::make_shared<core::pack::task::operation::OperateOnResidueSubset>(
			utility::pointer::make_shared<core::pack::task::operation::PreventRepackingRLT>(), design_sel, true );


		pack::task::TaskFactoryOP tf = utility::pointer::make_shared< pack::task::TaskFactory >();
		tf->push_back( restrict_task );

		// utility::vector1< core::conformation::ResidueCOP > res_vector;
		pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );
		rotsets = utility::pointer::make_shared< pack::rotamer_set::RotamerSets >();
		rotsets->set_task( task );


		core::select::residue_selector::ResidueSubset design_sub = design_sel->apply( pose );
		for ( core::Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
			res_vector.push_back( pose.residue(seqpos).get_self_ptr() );

			if ( ! design_sub[ seqpos ] ) continue;

			pack::rotamer_set::RotamerSetOP rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();
			rotset->set_resid( seqpos );
			if ( pose.residue(seqpos).is_virtual_residue() ) {
				rotset->add_rotamer( pose.residue(seqpos) );
			} else {

				for ( core::pose::Pose const & letter_pose : poses ) {
					rotset->add_rotamer( letter_pose.residue(seqpos) );

					if ( doubles ) {
						core::conformation::Residue backwards = *letter_pose.residue(seqpos).clone();
						if ( backwards.nchi() > 0 ) {
							backwards.set_chi( 1, backwards.chi(1) - 180 );
						}
						rotset->add_rotamer( backwards );
					}
				}
			}
			// std::cout << seqpos << " " << rotsets->resid_2_moltenres( seqpos ) << " " << rotsets->nmoltenres() << std::endl;
			rotsets->set_explicit_rotamers( rotsets->resid_2_moltenres( seqpos ), rotset );

		}

		rotsets->update_offset_data();

		sel = utility::pointer::make_shared<core::select::residue_selector::TrueResidueSelector>();
		options = utility::pointer::make_shared<SapConstraintOptions>();
		options->score_selector( sel );
		options->sap_calculate_selector( sel );


	}




	// void test_packing() {
	//     TR << "test_packing" << std::endl;
	//     for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
	//         bool using_symm = bool(int_sym);
	//         TR << (using_symm ? "Symm pass" : "first pass") << std::endl;



	//         core::pose::Pose pose;
	//         utility::vector1< core::conformation::ResidueCOP > res_vector;
	//         pack::rotamer_set::RotamerSetsOP initial_rotsets;
	//         SapConstraintOptionsOP options;
	//         select::residue_selector::ResidueSelectorCOP sel;
	//         Real pack_size = 0;

	//         prepare_packing_pose_and_rotsets( pose, res_vector, initial_rotsets, options, sel, pack_size, using_symm );


	//         SapConstraintHelperOP helper = utility::pointer::make_shared<SapConstraintHelper>( options );
	//         pack::rotamer_set::RotamerSets rotsets = helper->init_with_pose( pose, *initial_rotsets );

	//         Real init_score = helper->calculate_energy( res_vector, 0 );

	//         utility::vector1<Size> positions;
	//         for ( Size i = 1; i <= pose.size(); i++ ) {
	//             if ( pose.residue(i).is_virtual_residue() ) continue;
	//             positions.push_back(i);
	//             helper->shadow_mismatch_[i] = true;
	//             for ( Size j = 0; j < helper->max_rotamer_atoms_; j++ ) {
	//                 helper->dirty_sasa_[i][j] = true;
	//             }

	//         }

	//         helper->recalculate_sasa( );
	//         helper->commit_considered_substitution();

	//         assert_blocks( res_vector, helper );
	//         assert_nearest_5( res_vector, helper );
	//         assert_sasa_scores( res_vector, helper );
	//         assert_saps( pose, res_vector, helper );

	//         TS_ASSERT_DELTA( init_score, calculate_slow_approx_sap( pose, sel, sel, sel ), 0.001 );


	//         // random but repeatable
	//         utility::vector1<Size> order = { 9, 7, 14, 16, 6, 10, 13, 5, 17, 15, 8, 19, 2, 18, 12, 4, 1, 3, 11,
	//                                         13, 10, 15, 6, 18, 3, 8, 17, 1, 9, 4, 5, 16, 12, 19, 11, 14, 7, 2,
	//                                         9, 7, 4, 12, 6, 1, 3, 11, 19, 10, 13, 5, 2, 8, 17, 14, 18, 16, 15,
	//                                         11, 18, 6, 19, 14, 8, 10, 12, 3, 16, 5, 15, 7, 4, 1, 17, 2, 9, 13,
	//                                         5, 15, 9, 11, 14, 8, 17, 13, 18, 4, 3, 7, 1, 12, 2, 19, 6, 16, 10,
	//                                         2, 5, 7, 17, 14, 8, 6, 1, 9, 15, 3, 4, 19, 12, 11, 13, 18, 16, 10
	//         };



	//         for ( Size repeat = 1; repeat <= 3; repeat++ ) {

	//             for ( Size letter1 : order ) {
	//                 for ( core::Size seqpos = 1; seqpos <= pack_size; seqpos++ ) {

	//                     // std::cout << " " << repeat << " " << letter1 << " " << seqpos << std::endl;
	//                     core::conformation::ResidueCOP rotamer = rotsets.rotamer_set_for_residue( seqpos )->rotamer( letter1 );
	//                     utility::vector1< core::conformation::ResidueCOP > old_res_vector = res_vector;

	//                     // Do a repeatable random to get commit status
	//                     double temp = std::sqrt( repeat * letter1 * seqpos );
	//                     double deci = temp - (long)temp;
	//                     bool accept = deci > 0.5;

	//                     pose.replace_residue( seqpos, *rotamer, false );
	//                     res_vector[seqpos] = rotamer;
	//                     if ( pack_size != pose.size() ) {
	//                         for ( Size i = 1; i < SYMM; i++ ) {
	//                             Size sympos = seqpos + i * pack_size;
	//                             res_vector[ sympos ] = rotsets.rotamer_set_for_residue( sympos )->rotamer( letter1 );
	//                         }
	//                     }

	//                     Real score = helper->calculate_energy( res_vector, seqpos );

	//                     // It's really slow to calculate approx_sap, so we can't do it every time
	//                     double temp2 = std::sqrt( repeat * letter1 * seqpos * seqpos );
	//                     double deci2 = temp2 - (long)temp2;
	//                     bool accept2 = deci2 < 0.05;

	//                     if ( accept2 ) {
	//                         assert_blocks( res_vector, helper );
	//                         assert_nearest_5( res_vector, helper );
	//                         assert_sasa_scores( res_vector, helper );
	//                         assert_saps( pose, res_vector, helper );

	//                         TS_ASSERT_DELTA( score, calculate_slow_approx_sap( pose, sel, sel, sel ), 0.001 );
	//                     }

	//                     if ( ! accept ) {
	//                         pose.replace_residue( seqpos, *old_res_vector[seqpos], false );
	//                         res_vector = old_res_vector;
	//                     } else {
	//                         helper->commit_considered_substitution();
	//                     }
	//                 }
	//             }
	//         }
	//     }
	// }



	void test_partial_packing() {
		TR << "test_partial_packing" << std::endl;
		for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
			bool using_symm = bool(int_sym);
			TR << (using_symm ? "Symm pass" : "first pass") << std::endl;


			core::pose::Pose pose;
			utility::vector1< core::conformation::ResidueCOP > res_vector;
			pack::rotamer_set::RotamerSetsOP initial_rotsets;
			SapConstraintOptionsOP options;
			select::residue_selector::ResidueSelectorCOP sel;
			select::residue_selector::ResidueSelectorCOP design_sel =
				utility::pointer::make_shared<select::residue_selector::ResidueIndexSelector>("1,3,5");
			Real pack_size = 0;

			prepare_packing_pose_and_rotsets( pose, res_vector, initial_rotsets, options, sel, pack_size, using_symm, false, design_sel );

			core::select::residue_selector::ResidueSubset design_sub = design_sel->apply( pose );


			SapConstraintHelperOP helper = utility::pointer::make_shared<SapConstraintHelper>( options );
			pack::rotamer_set::RotamerSets rotsets = helper->init_with_pose( pose, *initial_rotsets );

			Real init_score = helper->calculate_energy( res_vector, 0 );

			utility::vector1<Size> positions;
			for ( Size i = 1; i <= pose.size(); i++ ) {
				if ( pose.residue(i).is_virtual_residue() ) continue;
				positions.push_back(i);
				helper->shadow_mismatch_[i] = true;
				for ( Size j = 0; j < helper->max_rotamer_atoms_; j++ ) {
					helper->dirty_sasa_[i][j] = true;
				}

			}
			helper->recalculate_sasa( );
			helper->commit_considered_substitution();

			assert_blocks( res_vector, helper );
			assert_nearest_5( res_vector, helper );
			assert_sasa_scores( res_vector, helper );
			assert_saps( pose, res_vector, helper );

			TS_ASSERT_DELTA( init_score, calculate_slow_approx_sap( pose, sel, sel, sel ), 0.001 );


			// random but repeatable
			utility::vector1<Size> order = { 9, 7, 14, 16, 6, 10, 13, 5, 17, 15, 8, 19, 2, 18, 12, 4, 1, 3, 11,
				// 13, 10, 15, 6, 18, 3, 8, 17, 1, 9, 4, 5, 16, 12, 19, 11, 14, 7, 2,
				// 9, 7, 4, 12, 6, 1, 3, 11, 19, 10, 13, 5, 2, 8, 17, 14, 18, 16, 15,
				// 11, 18, 6, 19, 14, 8, 10, 12, 3, 16, 5, 15, 7, 4, 1, 17, 2, 9, 13,
				// 5, 15, 9, 11, 14, 8, 17, 13, 18, 4, 3, 7, 1, 12, 2, 19, 6, 16, 10,
				2, 5, 7, 17, 14, 8, 6, 1, 9, 15, 3, 4, 19, 12, 11, 13, 18, 16, 10
				};



			for ( Size repeat = 1; repeat <= 1; repeat++ ) {

				for ( Size letter1 : order ) {
					for ( core::Size seqpos = 1; seqpos <= pack_size; seqpos++ ) {
						if ( ! design_sub[ seqpos ] ) continue;

						// std::cout << " " << repeat << " " << letter1 << " " << seqpos << std::endl;
						core::conformation::ResidueCOP rotamer = rotsets.rotamer_set_for_residue( seqpos )->rotamer( letter1 );
						utility::vector1< core::conformation::ResidueCOP > old_res_vector = res_vector;

						// Do a repeatable random to get commit status
						double temp = std::sqrt( repeat * letter1 * seqpos );
						double deci = temp - (long)temp;
						bool accept = deci > 0.5;

						pose.replace_residue( seqpos, *rotamer, false );
						res_vector[seqpos] = rotamer;
						if ( pack_size != pose.size() ) {
							for ( Size i = 1; i < SYMM; i++ ) {
								Size sympos = seqpos + i * pack_size;
								res_vector[ sympos ] = rotsets.rotamer_set_for_residue( sympos )->rotamer( letter1 );
							}
						}

						Real score = helper->calculate_energy( res_vector, seqpos );

						// It's really slow to calculate approx_sap, so we can't do it every time
						double temp2 = std::sqrt( repeat * letter1 * seqpos * seqpos );
						double deci2 = temp2 - (long)temp2;
						bool accept2 = deci2 > 0.05;

						if ( accept2 ) {
							assert_blocks( res_vector, helper );
							assert_nearest_5( res_vector, helper );
							assert_sasa_scores( res_vector, helper );
							assert_saps( pose, res_vector, helper );

							TS_ASSERT_DELTA( score, calculate_slow_approx_sap( pose, sel, sel, sel ), 0.001 );
						}

						if ( ! accept ) {
							pose.replace_residue( seqpos, *old_res_vector[seqpos], false );
							res_vector = old_res_vector;
						} else {
							helper->commit_considered_substitution();
						}
					}
				}
			}
		}
	}


	scoring::ScoreFunctionOP
	get_scorefunction( bool symmetric ) {
		if ( symmetric && old_symm_ ) {
			return utility::pointer::make_shared<core::scoring::symmetry::SymmetricScoreFunction>();
		} else {
			return utility::pointer::make_shared<scoring::ScoreFunction>();
		}
	}


	void test_mover() {
		TR << "test_mover" << std::endl;
		for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
			bool using_symm = bool(int_sym);
			TR << (using_symm ? "Symm pass" : "first pass") << std::endl;


			core::pose::Pose pose = create_trpcage_ideal_pose();
			terrible_C3_symmetry( pose, using_symm );

			core::select::residue_selector::TrueResidueSelectorOP sel =
				utility::pointer::make_shared<core::select::residue_selector::TrueResidueSelector>();


			protocols::simple_moves::AddSapConstraintMover mover;
			mover.set_score_selector( sel );
			mover.set_sap_calculate_selector( sel );

			mover.set_sap_goal( 10 );
			mover.set_penalty_per_sap( 2 );

			mover.apply( pose );


			scoring::ScoreFunctionOP sfxn = get_scorefunction( using_symm );
			sfxn->set_weight( core::scoring::sap_constraint, 1.0 );

			Real score = (*sfxn)(pose);

			Real sap = calculate_sap( pose, sel, sel, sel );
			Real expected_score = (sap - 10) * 2;

			TS_ASSERT_DELTA( score, expected_score, 0.001 );

			pose.remove_constraints();

			mover.set_packing_correction( -2.3 );
			expected_score = (sap-2.3 - 10)*2;
			mover.apply( pose );
			score = (*sfxn)(pose);
			TS_ASSERT_DELTA( score, expected_score, 0.001 );
			pose.remove_constraints();

			SapScoreMetric my_metric( sel, sel );
			TS_ASSERT_DELTA( sap, my_metric.calculate( pose ), 0.001 );

			PerResidueSapScoreMetric my_per_res_metric( sel, sel );
			std::map< Size, Real> per_res_sap = my_per_res_metric.calculate( pose );
			Real my_sum = 0;
			for ( auto const & pair : per_res_sap ) my_sum += pair.second;
			TS_ASSERT_DELTA( sap, my_sum, 0.001 );


			// Make sure the different speeds all end up calculating sap the same
			{
				pose.remove_constraints();
				protocols::simple_moves::AddSapConstraintMover mover3;
				mover3.set_score_selector( sel );
				mover3.set_sap_calculate_selector( sel );

				mover3.set_sap_goal( 0 );
				mover3.set_penalty_per_sap( 1 );
				mover3.set_speed("slow");

				mover3.apply( pose );

				TS_ASSERT_DELTA( sap, sfxn->score( pose ), 0.001 );
			}
			{
				pose.remove_constraints();
				protocols::simple_moves::AddSapConstraintMover mover3;
				mover3.set_score_selector( sel );
				mover3.set_sap_calculate_selector( sel );

				mover3.set_sap_goal( 0 );
				mover3.set_penalty_per_sap( 1 );
				mover3.set_speed("fast");

				mover3.apply( pose );

				TS_ASSERT_DELTA( sap, sfxn->score( pose ), 0.001 );
			}
			{
				pose.remove_constraints();
				protocols::simple_moves::AddSapConstraintMover mover3;
				mover3.set_score_selector( sel );
				mover3.set_sap_calculate_selector( sel );

				mover3.set_sap_goal( 0 );
				mover3.set_penalty_per_sap( 1 );
				mover3.set_speed("lightning");

				mover3.apply( pose );

				TS_ASSERT_DELTA( sap, sfxn->score( pose ), 0.001 );
				pose.remove_constraints();
			}



			protocols::simple_moves::AddSapConstraintMover mover2;
			mover2.set_score_selector( sel );
			mover2.set_sap_calculate_selector( sel );

			mover2.set_sap_goal( 1000 );
			mover2.set_sap_lb_goal( 100 );
			mover2.set_penalty_per_sap( 2 );

			mover2.apply( pose );


			score = (*sfxn)(pose);

			sap = calculate_sap( pose, sel, sel, sel );
			expected_score = (100 - sap) * 2;

			TS_ASSERT_DELTA( score, expected_score, 0.001 );





			utility::vector1<bool> subset( pose.size(), true );
			for ( Size i = 1; i<= pose.size(); i++ ) {
				if ( pose.residue(i).is_virtual_residue() ) subset[i] = false;
			}

			// This just makes sure it doesn't crash on packing, not a great test
			protocols::toolbox::pose_manipulation::repack_these_residues( subset, pose, sfxn, false, "AV" );

		}

	}


	void test_math_mover() {
		TR << "test_math_mover" << std::endl;
		for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
			bool using_symm = bool(int_sym);
			TR << (using_symm ? "Symm pass" : "first pass") << std::endl;


			core::pose::Pose pose = create_trpcage_ideal_pose();
			terrible_C3_symmetry( pose, using_symm );

			core::select::residue_selector::TrueResidueSelectorOP true_sel =
				utility::pointer::make_shared<core::select::residue_selector::TrueResidueSelector>();

			core::select::residue_selector::ResidueIndexSelectorOP sel1 =
				utility::pointer::make_shared<core::select::residue_selector::ResidueIndexSelector>("1,2,3");
			core::select::residue_selector::ResidueIndexSelectorOP sel2 =
				utility::pointer::make_shared<core::select::residue_selector::ResidueIndexSelector>("4,5,6");


			scoring::ScoreFunctionOP sfxn = get_scorefunction( using_symm );
			sfxn->set_weight( core::scoring::sap_constraint, 1.0 );

			Real sap1 = calculate_sap( pose, sel1, true_sel, true_sel );
			Real sap2 = calculate_sap( pose, sel2, true_sel, true_sel );

			std::cout << "Test: " << sap1 << " " << sap2 << std::endl;
			Real score = 0;
			Real expected_score = 0;

			{

				protocols::simple_moves::AddSapConstraintMover mover1;
				mover1.set_score_selector( sel1 );
				mover1.set_sap_calculate_selector( true_sel );
				mover1.set_name("sap1");
				mover1.set_penalty_per_sap( 0 );

				protocols::simple_moves::AddSapConstraintMover mover2;
				mover2.set_score_selector( sel2 );
				mover2.set_sap_calculate_selector( true_sel );
				mover2.set_name("sap2");
				mover2.set_penalty_per_sap( 0 );

				protocols::simple_moves::AddSapMathConstraintMover math_mover;
				math_mover.add_constraint( 1.4, "sap1" );
				math_mover.add_constraint( -1.7, "sap2" );
				math_mover.set_upper_bound( 0 );
				math_mover.set_lower_bound( 0 );
				math_mover.set_penalty_per_unit( 1 );

				mover1.apply( pose );
				mover2.apply( pose );
				math_mover.apply( pose );


				score = sfxn->score( pose );
				expected_score = std::abs(sap1 * 1.4 - sap2 * 1.7);

				TS_ASSERT_DELTA( score, expected_score, 0.001 );

				{
					pose.remove_constraints();

					math_mover.set_upper_bound( utility::get_undefined_real() );
					math_mover.set_lower_bound( 1000 );

					mover1.apply( pose );
					mover2.apply( pose );
					math_mover.apply( pose );

					score = sfxn->score( pose );
					expected_score = 1000 - (sap1 * 1.4 - sap2 * 1.7);

					TS_ASSERT_DELTA( score, expected_score, 0.001 );
				}

				{
					pose.remove_constraints();

					math_mover.set_upper_bound( utility::get_undefined_real() );
					math_mover.set_lower_bound( -1000 );

					mover1.apply( pose );
					mover2.apply( pose );
					math_mover.apply( pose );

					score = sfxn->score( pose );
					expected_score = 0;

					TS_ASSERT_DELTA( score, expected_score, 0.001 );
				}
				{
					pose.remove_constraints();

					math_mover.set_upper_bound( 1000 );
					math_mover.set_lower_bound( utility::get_undefined_real() );

					mover1.apply( pose );
					mover2.apply( pose );
					math_mover.apply( pose );

					score = sfxn->score( pose );
					expected_score = 0;

					TS_ASSERT_DELTA( score, expected_score, 0.001 );
				}

				{
					pose.remove_constraints();

					math_mover.set_upper_bound( -1000 );
					math_mover.set_lower_bound( utility::get_undefined_real() );

					mover1.apply( pose );
					mover2.apply( pose );
					math_mover.apply( pose );

					score = sfxn->score( pose );
					expected_score = (sap1 * 1.4 - sap2 * 1.7) + 1000;

					TS_ASSERT_DELTA( score, expected_score, 0.001 );
				}
			}



			// This just makes sure it doesn't crash on packing, not a great test
			protocols::toolbox::pose_manipulation::repack_these_residues( true_sel->apply(pose), pose, sfxn, false, "AV" );

		}

	}


	void test_selectors() {
		TR << "test_selectors" << std::endl;
		for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
			bool using_symm = bool(int_sym);
			TR << (using_symm ? "Symm pass" : "first pass") << std::endl;


			core::pose::Pose pose = create_trpcage_ideal_pose();
			terrible_C3_symmetry( pose, using_symm );

			core::select::residue_selector::TrueResidueSelectorOP true_sel =
				utility::pointer::make_shared<core::select::residue_selector::TrueResidueSelector>();

			core::select::residue_selector::ResidueNameSelectorOP hydrophobic_sel =
				utility::pointer::make_shared<core::select::residue_selector::ResidueNameSelector>(
				"ALA,CYS,PHE,ILE,LEU,MET,THR,PRO,VAL,TRP,TYR", true );

			core::select::residue_selector::ResidueSubset hydrophobic_sub = hydrophobic_sel->apply( pose );

			Real total_sap = calculate_sap( pose, true_sel, true_sel, true_sel );
			Real hydrophobic_sap = calculate_sap( pose, true_sel, hydrophobic_sel, hydrophobic_sel );

			// Mostly just make sure they aren't equal
			TS_ASSERT( hydrophobic_sap > total_sap );

			Real total_sum = 0;
			Real hydrophobic_sum = 0;

			for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
				if ( pose.residue(seqpos).is_virtual_residue() ) continue;

				core::select::residue_selector::ResidueIndexSelectorOP sel =
					utility::pointer::make_shared<core::select::residue_selector::ResidueIndexSelector>(
					boost::str(boost::format("%i")%seqpos) );

				Real this_sap = calculate_sap( pose, sel, true_sel, true_sel );
				Real this_hydrophobic_sap = calculate_sap( pose, sel, hydrophobic_sel, hydrophobic_sel );

				if ( ! hydrophobic_sub[ seqpos ] ) {
					TS_ASSERT( this_hydrophobic_sap == 0 );
				}

				total_sum += this_sap;
				hydrophobic_sum += this_hydrophobic_sap;
			}

			TS_ASSERT_DELTA( total_sap, total_sum, 0.001 );
			TS_ASSERT_DELTA( hydrophobic_sap, hydrophobic_sum, 0.001 );


			// We're simulating an interface here by splitting 1ubq in half
			if ( ! using_symm ) {
				pose.clear();
				core::import_pose::pose_from_file( pose, "core/pack/dunbrack/1UBQ_repack.pdb" );

				core::select::residue_selector::ResidueIndexSelectorOP first_half =
					utility::pointer::make_shared<core::select::residue_selector::ResidueIndexSelector>( "1-38" );

				Real complex_sap = calculate_sap( pose, first_half, first_half, true_sel );
				Real split_sap = calculate_sap( pose, first_half, first_half, first_half );

				TS_ASSERT( split_sap > complex_sap );
			}





		}

	}



	void assert_fast_sasa_scores(
		SapConstraintHelperOP const & helper,
		core::pose::Pose const & pose,
		pack::rotamer_set::RotamerSets const & rotsets,
		utility::vector1<utility::vector1<float>> & sasa_scores
	) {

		for ( Size irot = 1; irot <= rotsets.nrotamers() + pose.size(); irot++ ) {

			core::conformation::ResidueCOP rotamer = irot <= rotsets.nrotamers() ? rotsets.rotamer( irot )
				: pose.residue( irot - rotsets.nrotamers() ).get_self_ptr();

			if ( rotamer->is_virtual_residue() ) continue;

			TS_ASSERT( helper->rotamer_to_sasa_data_.count( &*rotamer ) == 1 );
			if ( helper->rotamer_to_sasa_data_.count( &*rotamer ) == 0 ) continue;

			float * helper_scores = helper->rotamer_to_sasa_data_.at( &*rotamer );

			Size first_sidechain = rotamer->first_sidechain_atom();
			for ( Size iat = 0; (int)iat <= (int)rotamer->natoms() - (int)first_sidechain; iat++ ) {
				Size iatom = first_sidechain + iat;
				if ( rotamer->atom_is_backbone( iatom ) ) continue;
				TS_ASSERT_DELTA( sasa_scores[irot][iatom], helper_scores[iat], 0.001 );
				sasa_scores[irot][iatom] = helper_scores[iat]; // clean out small numerical differences
			}

		}
	}

	void
	assert_correct_sasa_ptrs(
		utility::vector1< core::conformation::ResidueCOP > const & res_vector,
		SapConstraintHelperOP const & helper
	) {
		for ( Size seqpos = 1; seqpos <= res_vector.size(); seqpos++ ) {
			if ( res_vector[ seqpos ]->is_virtual_residue() ) continue;

			float * ptr = helper->rotamer_to_sasa_data_.at( &*res_vector[ seqpos ] );

			assert( ptr );

			TS_ASSERT( ptr == helper->atom_sasa_score_fast_[seqpos] );
		}
	}

	Real
	calculate_fast_sap(
		core::pose::Pose const & pose,
		utility::vector1<core::conformation::ResidueCOP> const & all_rotamers,
		utility::vector1<Size> irot_at_position,
		utility::vector1<utility::vector1<float>> & sasa_scores, // const but whatever
		SapConstraintHelperOP const & check_helper,
		SapConstraintHelperOP const & scratch_helper


	) {



		SapConstraintHelperOP helper = scratch_helper;

		utility::vector1<Size> to_recalc;
		utility::vector1<core::conformation::ResidueCOP> my_resvect;

		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

			Size irot = irot_at_position[ seqpos ];
			core::conformation::ResidueCOP const & rotamer = all_rotamers[irot]; //irot <= rotsets->nrotamers() ? rotsets->rotamer( irot )
			//: pose.residue( irot - rotsets->nrotamers() ).get_self_ptr();

			my_resvect.push_back( rotamer );
			assert( rotamer->natoms() == pose.residue(seqpos).natoms() );
		}

		helper->calculate_energy( my_resvect, 0 );

		// Now we insert the new sasa scores
		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
			if ( my_resvect[ seqpos ]->is_virtual_residue() ) continue;

			Size irot = irot_at_position[ seqpos ];
			core::conformation::ResidueCOP const & rotamer = my_resvect[ seqpos ];
			Size first_sidechain = rotamer->first_sidechain_atom();
			Size natoms = rotamer->natoms();

			float * sasa_ptr = (&sasa_scores[irot].front()) + first_sidechain - 1;
			helper->atom_sasa_score_fast_[ seqpos ] = sasa_ptr;

			to_recalc.push_back( seqpos );

			std::fill( helper->atom_sap_[seqpos].begin(), helper->atom_sap_[seqpos].end(), 0);

			for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
				Size iatom = iat + first_sidechain;
				if ( rotamer->atom_is_backbone( iatom ) ) continue;
				TS_ASSERT_DELTA( sasa_ptr[iat], check_helper->atom_sasa_score_fast_[seqpos][iat], 0.0001 );
			}
		}

		assert( helper->fast_ );

		helper->current_score_ = 0;
		helper->recalculate_saps( to_recalc );

		// std::cout << "before" << std::endl;
		assert_saps( pose, my_resvect, helper, true );
		// std::cout << "after" << std::endl;

		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

			core::conformation::ResidueCOP const & rotamer = my_resvect[ seqpos ];
			if ( rotamer->is_virtual_residue() ) continue;
			Size first_sidechain = rotamer->first_sidechain_atom();
			Size natoms = rotamer->natoms();

			for ( Size iat = 0; (int)iat <= (int)natoms - (int)first_sidechain; iat++ ) {
				Size iatom = iat + first_sidechain;
				if ( rotamer->atom_is_backbone( iatom ) ) continue;
				TS_ASSERT_DELTA( helper->atom_sap_[seqpos][iat], check_helper->atom_sap_[seqpos][iat], 0.001 );
			}
		}

		TS_ASSERT_DELTA( helper->current_score_, check_helper->current_score_, 0.01 );


		return helper->current_score_;
	}




	// void test_fast_packing() {
	//     TR << "test_fast_packing" << std::endl;
	//     for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
	//         bool using_symm = bool(int_sym);
	//         TR << (using_symm ? "Symm pass" : "first pass") << std::endl;


	//         core::pose::Pose pose;
	//         utility::vector1< core::conformation::ResidueCOP > res_vector;
	//         pack::rotamer_set::RotamerSetsOP initial_rotsets;
	//         SapConstraintOptionsOP options;
	//         select::residue_selector::ResidueSelectorCOP sel;
	//         Real pack_size = 0;

	//         prepare_packing_pose_and_rotsets( pose, res_vector, initial_rotsets, options, sel, pack_size, using_symm );
	//         options->fast( true );

	//         SapConstraintHelperOP scratch_helper = utility::pointer::make_shared<SapConstraintHelper>( options );
	//         scratch_helper->init_with_pose( pose, *initial_rotsets );

	//         SapConstraintHelperOP helper = utility::pointer::make_shared<SapConstraintHelper>( options );
	//         pack::rotamer_set::RotamerSets rotsets = helper->init_with_pose( pose, *initial_rotsets );

	//         utility::vector1<core::conformation::ResidueCOP> all_rotamers;
	//         for ( Size i = 1; i<= rotsets.nrotamers(); i++ ) {
	//             all_rotamers.push_back( rotsets.rotamer( i ) );
	//         }
	//         for ( Size i = 1; i <= pose.size(); i++ ) {
	//             all_rotamers.push_back( pose.residue(i).get_self_ptr() );
	//         }


	//         utility::vector1<utility::vector1<Real>> blocks = manually_calculate_fast_blocks( pose, rotsets );
	//         utility::vector1<utility::vector1<float>> sasa_scores = manually_calculate_fast_sasa_scores( pose, rotsets, blocks );

	//         assert_fast_sasa_scores( helper, pose, rotsets, sasa_scores );

	//         Real init_score = helper->calculate_energy( res_vector, 0 );


	//         assert_nearest_5( res_vector, helper );
	//         assert_correct_sasa_ptrs( res_vector, helper );
	//         assert_saps( pose, res_vector, helper, true );

	//         utility::vector1<Size> to_recalc;

	//         utility::vector1<Size> irot_at_position;
	//         for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
	//             irot_at_position.push_back( rotsets.nrotamers() + seqpos );
	//             to_recalc.push_back( seqpos );
	//         }

	//         // for ( Size i = 1; i <= pose.size(); i++ ) {
	//         //     core::conformation::Residue const & pose_res = pose.residue(i);
	//         //     core::conformation::Residue const & vector_res = *res_vector[i];
	//         //     Size irot = irot_at_position[ i ];
	//         //     core::conformation::Residue const & irot_res = irot <= rotsets->nrotamers() ? *rotsets->rotamer( irot )
	//         //                                                                         : pose.residue(irot - rotsets->nrotamers() );

	//         //     std::cout << pose_res.natoms() << " " << vector_res.natoms() << " " << irot_res.natoms() << std::endl;
	//         // }

	//         TS_ASSERT_DELTA( init_score, calculate_fast_sap( pose, all_rotamers, irot_at_position, sasa_scores, helper, scratch_helper ), 0.001 );


	//         // std::cout << "After" << std::endl;
	//         // random but repeatable
	//         utility::vector1<Size> order = { 9, 7, 14, 16, 6, 10, 13, 5, 17, 15, 8, 19, 2, 18, 12, 4, 1, 3, 11,
	//                                         13, 10, 15, 6, 18, 3, 8, 17, 1, 9, 4, 5, 16, 12, 19, 11, 14, 7, 2,
	//                                         9, 7, 4, 12, 6, 1, 3, 11, 19, 10, 13, 5, 2, 8, 17, 14, 18, 16, 15,
	//                                         11, 18, 6, 19, 14, 8, 10, 12, 3, 16, 5, 15, 7, 4, 1, 17, 2, 9, 13,
	//                                         5, 15, 9, 11, 14, 8, 17, 13, 18, 4, 3, 7, 1, 12, 2, 19, 6, 16, 10,
	//                                         2, 5, 7, 17, 14, 8, 6, 1, 9, 15, 3, 4, 19, 12, 11, 13, 18, 16, 10
	//         };



	//         for ( Size repeat = 1; repeat <= 3; repeat++ ) {

	//             for ( Size letter1 : order ) {
	//                 for ( core::Size seqpos = 1; seqpos <= pack_size; seqpos++ ) {




	//                     // std::cout << " " << repeat << " " << letter1 << " " << seqpos << std::endl;
	//                     core::conformation::ResidueCOP rotamer = rotsets.rotamer_set_for_residue( seqpos )->rotamer( letter1 );
	//                     utility::vector1< core::conformation::ResidueCOP > old_res_vector = res_vector;
	//                     utility::vector1<Size> old_irot_at_position = irot_at_position;

	//                     irot_at_position[ seqpos ] = rotsets.nrotamer_offset_for_moltenres( rotsets.resid_2_moltenres( seqpos ) ) + letter1;
	//                     assert( &*rotamer == &*rotsets.rotamer( irot_at_position[ seqpos ] ) );

	//                     // Do a repeatable random to get commit status
	//                     double temp = std::sqrt( repeat * letter1 * seqpos );
	//                     double deci = temp - (long)temp;
	//                     bool accept = deci > 0.5;
	//                     accept = true;

	//                     pose.replace_residue( seqpos, *rotamer, false );
	//                     res_vector[seqpos] = rotamer;

	//                     if ( pack_size != pose.size() ) {
	//                         for ( Size i = 1; i < SYMM; i++ ) {
	//                             Size sympos = seqpos + i * pack_size;
	//                             res_vector[ sympos ] = rotsets.rotamer_set_for_residue( sympos )->rotamer( letter1 );

	//                             irot_at_position[ sympos ] = rotsets.nrotamer_offset_for_moltenres( rotsets.resid_2_moltenres( sympos ) ) + letter1;
	//                             assert( &*res_vector[ sympos ] == &*rotsets.rotamer( irot_at_position[ sympos ] ) );
	//                         }
	//                     }



	//                     Real score = helper->calculate_energy( res_vector, seqpos );

	//                     // for ( Size i = 1; i<= pose.size(); i++ ) {
	//                     //     for ( Size j = 0; j< helper->dirty_sasa_[i].size(); j++ ) {
	//                     //         assert( helper->dirty_sasa_[i][j] );
	//                     //         helper->atom_sap_[i][j] = 0;
	//                     //     }
	//                     // }
	//                     // helper->current_score_ = 0;
	//                     // helper->recalculate_saps( to_recalc );

	//                     assert_nearest_5( res_vector, helper );
	//                     assert_correct_sasa_ptrs( res_vector, helper );
	//                     // std::cout << "new before" << std::endl;
	//                     assert_saps( pose, res_vector, helper, true );
	//                     // std::cout << "new after" << std::endl;

	//                     score = helper->current_score_;


	//                     // for ( Size i = 1; i <= pose.size(); i++ ) {
	//                     //     core::conformation::Residue const & pose_res = pose.residue(i);
	//                     //     core::conformation::Residue const & vector_res = *res_vector[i];
	//                     //     Size irot = irot_at_position[ i ];
	//                     //     core::conformation::Residue const & irot_res = irot <= rotsets->nrotamers() ? *rotsets->rotamer( irot )
	//                     //                                                                         : pose.residue(irot - rotsets->nrotamers() );

	//                     //     std::cout << pose_res.natoms() << " " << vector_res.natoms() << " " << irot_res.natoms() << std::endl;
	//                     // }

	//                     TS_ASSERT_DELTA( score, calculate_fast_sap( pose, all_rotamers, irot_at_position, sasa_scores, helper, scratch_helper ), 0.01 );

	//                     if ( ! accept ) {
	//                         pose.replace_residue( seqpos, *res_vector[seqpos], false );
	//                         res_vector = old_res_vector;
	//                         irot_at_position = old_irot_at_position;
	//                     } else {
	//                         helper->commit_considered_substitution();
	//                     }
	//                 }
	//             }
	//         }
	//     }
	// }


	void test_fast_partial_packing() {
		TR << "test_fast_partial_packing" << std::endl;
		for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
			bool using_symm = bool(int_sym);
			TR << (using_symm ? "Symm pass" : "first pass") << std::endl;

			core::pose::Pose pose;
			utility::vector1< core::conformation::ResidueCOP > res_vector;
			pack::rotamer_set::RotamerSetsOP initial_rotsets;
			SapConstraintOptionsOP options;
			select::residue_selector::ResidueSelectorCOP sel;
			Real pack_size = 0;

			select::residue_selector::ResidueSelectorCOP design_sel =
				utility::pointer::make_shared<select::residue_selector::ResidueIndexSelector>("1,3,5");
			// utility::pointer::make_shared<select::residue_selector::ResidueIndexSelector>("1-6");

			prepare_packing_pose_and_rotsets( pose, res_vector, initial_rotsets, options, sel, pack_size, using_symm, false, design_sel );

			core::select::residue_selector::ResidueSubset design_sub = design_sel->apply( pose );
			options->fast( true );

			SapConstraintHelperOP scratch_helper = utility::pointer::make_shared<SapConstraintHelper>( options );
			scratch_helper->init_with_pose( pose, *initial_rotsets );

			SapConstraintHelperOP helper = utility::pointer::make_shared<SapConstraintHelper>( options );
			pack::rotamer_set::RotamerSets rotsets = helper->init_with_pose( pose, *initial_rotsets );

			utility::vector1<core::conformation::ResidueCOP> all_rotamers;
			for ( Size i = 1; i<= rotsets.nrotamers(); i++ ) {
				all_rotamers.push_back( rotsets.rotamer( i ) );
			}
			for ( Size i = 1; i <= pose.size(); i++ ) {
				all_rotamers.push_back( pose.residue(i).get_self_ptr() );
			}


			utility::vector1<utility::vector1<Real>> blocks = manually_calculate_fast_blocks( pose, rotsets );
			utility::vector1<utility::vector1<float>> sasa_scores = manually_calculate_fast_sasa_scores( pose, rotsets, blocks );


			// std::cout << "T1" << std::endl;
			assert_fast_sasa_scores( helper, pose, rotsets, sasa_scores );
			// std::cout << "End T1" << std::endl;

			Real init_score = helper->calculate_energy( res_vector, 0 );

			assert_nearest_5( res_vector, helper );
			assert_correct_sasa_ptrs( res_vector, helper );
			assert_saps( pose, res_vector, helper, true );

			utility::vector1<Size> to_recalc;

			utility::vector1<Size> irot_at_position;
			for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
				irot_at_position.push_back( rotsets.nrotamers() + seqpos );
				to_recalc.push_back( seqpos );
			}

			// for ( Size i = 1; i <= pose.size(); i++ ) {
			//     core::conformation::Residue const & pose_res = pose.residue(i);
			//     core::conformation::Residue const & vector_res = *res_vector[i];
			//     Size irot = irot_at_position[ i ];
			//     core::conformation::Residue const & irot_res = irot <= rotsets->nrotamers() ? *rotsets->rotamer( irot )
			//                                                                         : pose.residue(irot - rotsets->nrotamers() );

			//     std::cout << pose_res.natoms() << " " << vector_res.natoms() << " " << irot_res.natoms() << std::endl;
			// }

			TS_ASSERT_DELTA( init_score, calculate_fast_sap( pose, all_rotamers, irot_at_position, sasa_scores, helper, scratch_helper ), 0.001 );


			// std::cout << "After" << std::endl;
			// random but repeatable
			utility::vector1<Size> order = { 9, 7, 14, 16, 6, 10, 13, 5, 17, 15, 8, 19, 2, 18, 12, 4, 1, 3, 11,
				// 13, 10, 15, 6, 18, 3, 8, 17, 1, 9, 4, 5, 16, 12, 19, 11, 14, 7, 2,
				// 9, 7, 4, 12, 6, 1, 3, 11, 19, 10, 13, 5, 2, 8, 17, 14, 18, 16, 15,
				// 11, 18, 6, 19, 14, 8, 10, 12, 3, 16, 5, 15, 7, 4, 1, 17, 2, 9, 13,
				// 5, 15, 9, 11, 14, 8, 17, 13, 18, 4, 3, 7, 1, 12, 2, 19, 6, 16, 10,
				2, 5, 7, 17, 14, 8, 6, 1, 9, 15, 3, 4, 19, 12, 11, 13, 18, 16, 10
				};



			for ( Size repeat = 1; repeat <= 1; repeat++ ) {

				for ( Size letter1 : order ) {
					for ( core::Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
						if ( ! design_sub[seqpos] ) continue;




						// std::cout << " " << repeat << " " << letter1 << " " << seqpos << std::endl;
						core::conformation::ResidueCOP rotamer = rotsets.rotamer_set_for_residue( seqpos )->rotamer( letter1 );
						utility::vector1< core::conformation::ResidueCOP > old_res_vector = res_vector;
						utility::vector1<Size> old_irot_at_position = irot_at_position;

						irot_at_position[ seqpos ] = rotsets.nrotamer_offset_for_moltenres( rotsets.resid_2_moltenres( seqpos ) ) + letter1;
						assert( &*rotamer == &*rotsets.rotamer( irot_at_position[ seqpos ] ) );

						// Do a repeatable random to get commit status
						//double temp = std::sqrt( repeat * letter1 * seqpos );
						//double deci = temp - (long)temp;
						bool accept = true; //deci > 0.5;

						pose.replace_residue( seqpos, *rotamer, false );
						res_vector[seqpos] = rotamer;

						if ( pack_size != pose.size() ) {
							for ( Size i = 1; i < SYMM; i++ ) {
								Size sympos = seqpos + i * pack_size;
								res_vector[ sympos ] = rotsets.rotamer_set_for_residue( sympos )->rotamer( letter1 );

								irot_at_position[ sympos ] = rotsets.nrotamer_offset_for_moltenres( rotsets.resid_2_moltenres( sympos ) ) + letter1;
								assert( &*res_vector[ sympos ] == &*rotsets.rotamer( irot_at_position[ sympos ] ) );
							}
						}


						helper->calculate_energy( res_vector, seqpos );

						// for ( Size i = 1; i<= pose.size(); i++ ) {
						//     for ( Size j = 0; j< helper->dirty_sasa_[i].size(); j++ ) {
						//         assert( helper->dirty_sasa_[i][j] );
						//         helper->atom_sap_[i][j] = 0;
						//     }
						// }
						// helper->current_score_ = 0;
						// helper->recalculate_saps( to_recalc );

						assert_nearest_5( res_vector, helper );
						assert_correct_sasa_ptrs( res_vector, helper );
						// std::cout << "new before" << std::endl;
						assert_saps( pose, res_vector, helper, true );
						// std::cout << "new after" << std::endl;

						Real score = helper->current_score_;


						// for ( Size i = 1; i <= pose.size(); i++ ) {
						//     core::conformation::Residue const & pose_res = pose.residue(i);
						//     core::conformation::Residue const & vector_res = *res_vector[i];
						//     Size irot = irot_at_position[ i ];
						//     core::conformation::Residue const & irot_res = irot <= rotsets->nrotamers() ? *rotsets->rotamer( irot )
						//                                                                         : pose.residue(irot - rotsets->nrotamers() );

						//     std::cout << pose_res.natoms() << " " << vector_res.natoms() << " " << irot_res.natoms() << std::endl;
						// }

						TS_ASSERT_DELTA( score, calculate_fast_sap( pose, all_rotamers, irot_at_position, sasa_scores, helper, scratch_helper ), 0.01 );

						if ( ! accept ) {
							pose.replace_residue( seqpos, *res_vector[seqpos], false );
							res_vector = old_res_vector;
							irot_at_position = old_irot_at_position;
						} else {
							helper->commit_considered_substitution();
						}
					}
				}
			}
		}

	}

	void assert_lightning_sasa_scores(
		SapConstraintHelperOP const & helper,
		pack::rotamer_set::RotamerSets const & rotsets,
		utility::vector1<utility::vector1<float>> & sasa_scores
	) {

		for ( Size irot = 1; irot <= rotsets.nrotamers(); irot++ ) {

			core::conformation::ResidueCOP rotamer = rotsets.rotamer( irot );
			if ( rotamer->is_virtual_residue() ) continue;

			TS_ASSERT( helper->rotamer_to_sasa_data_.count( &*rotamer ) == 1 );
			if ( helper->rotamer_to_sasa_data_.count( &*rotamer ) == 0 ) continue;

			float * helper_scores = helper->rotamer_to_sasa_data_.at( &*rotamer );

			Size first_sidechain = rotamer->first_sidechain_atom();
			for ( Size iat = 0; (int)iat <= (int)rotamer->natoms() - (int)first_sidechain; iat++ ) {
				Size iatom = first_sidechain + iat;
				// std::cout << irot << " " << iatom << " " << rotamer->atom_name(iatom) << std::endl;
				if ( rotamer->atom_is_backbone( iatom ) ) continue;
				TS_ASSERT_DELTA( sasa_scores[irot][iatom], helper_scores[iat], 0.001 );
				// sasa_scores[irot][iatom] = helper_scores[iat]; // clean out small numerical differences
			}

		}
	}

	Real
	lightning_sap_from_sasa_scores_and_assert_1b_2b(
		core::pose::Pose const & pose,
		utility::vector1<utility::vector1<float>> const & sasa_scores,
		utility::vector1<Size> const & irot_at_position,
		std::unordered_map< Size, core::conformation::ResidueCOP > const & seq100_aa_to_rotamer,
		SapConstraintHelperOP helper
	) {
		Real sap_score = 0;

		for ( Size seqpos1 = 1; seqpos1 <= pose.size(); seqpos1++ ) {
			core::conformation::Residue const & pre_res1 = pose.residue(seqpos1);
			if ( pre_res1.is_virtual_residue() ) continue;
			Size lookup = seqpos1 * 100 + Size(pre_res1.aa());
			core::conformation::Residue const & res1 = seq100_aa_to_rotamer.count( lookup ) ? *seq100_aa_to_rotamer.at( lookup ) : pre_res1;

			Real seqpos1_sum = 0;

			for ( Size seqpos2 = 1; seqpos2 <= pose.size(); seqpos2++ ) {
				core::conformation::Residue const & pre_res2 = pose.residue(seqpos2);
				if ( pre_res2.is_virtual_residue() ) continue;
				Size lookup = seqpos2 * 100 + Size(pre_res2.aa());
				core::conformation::Residue const & res2 = seq100_aa_to_rotamer.count( lookup ) ? *seq100_aa_to_rotamer.at( lookup ) : pre_res2;

				// std::cout << "test from " << seqpos2 << " to " << seqpos1 << std::endl;

				Real two_onto_1 = 0;

				// This is actually an assert to test another function in this test
				TS_ASSERT( sasa_scores[irot_at_position[seqpos2]].size() == res2.natoms() );

				for ( Size at2 = 1; at2 <= res2.natoms(); at2++ ) {
					if ( res2.atom_is_backbone(at2) ) continue;

					for ( Size at1 = 1; at1 <= res1.natoms(); at1++ ) {
						if ( res1.atom_is_backbone(at1) ) continue;

						Real dist2 = res1.xyz(at1).distance_squared(res2.xyz(at2));

						if ( dist2 > 5*5 ) continue;

						two_onto_1 += sasa_scores[irot_at_position[seqpos2]][at2];

						// std::cout << "cast " << at2 - res2.first_sidechain_atom() << " " << sasa_scores[seqpos2][at2] << std::endl;
					}

					// check 1b and 2b here

				}
				if ( seqpos1 == seqpos2 ) {
					Real stored = helper->lightning_1b_lookup( res1.aa(), seqpos1 );
					TS_ASSERT_DELTA( two_onto_1, stored, 0.1 );
				} else {
					Real stored = 0;
					if ( seqpos1 < seqpos2 ) {
						stored = helper->lightning_2b_lookup( res1.aa(), seqpos1, res2.aa(), seqpos2 ).first;
						// std::cout << "A " << res1.name1() << " " << seqpos1 << " " << res2.name1() << " " << seqpos2 << " -- " << stored << ", " << two_onto_1 << std::endl;
					} else {
						stored = helper->lightning_2b_lookup( res2.aa(), seqpos2, res1.aa(), seqpos1 ).second;
						// std::cout << "B " << res2.name1() << " " << seqpos2 << " " << res1.name1() << " " << seqpos1 << " -- " << stored << ", " << two_onto_1 << std::endl;
					}
					TS_ASSERT_DELTA( two_onto_1, stored, 0.1 );
				}

				seqpos1_sum += two_onto_1;
			}
			if ( seqpos1_sum > 0 ) {
				sap_score += seqpos1_sum;
			}
		}
		return sap_score;
	}




	void test_helper_lightning(){
		TR << "test_helper_lightning" << std::endl;
		for ( Size int_sym = 0; int_sym <= 2; int_sym++ ) {
			bool using_symm = bool(int_sym);
			bool using_map = int_sym == 2;
			TR << (using_symm ? "Symm pass" : "first pass") << std::endl;
			SapDatabase::get_instance()->symm_debug_force_map_ = using_map;

			core::pose::Pose pose = create_trpcage_ideal_pose();
			terrible_C3_symmetry( pose, using_symm );


			core::select::residue_selector::TrueResidueSelectorOP sel =
				utility::pointer::make_shared<core::select::residue_selector::TrueResidueSelector>();


			SapConstraintOptionsOP options = utility::pointer::make_shared<SapConstraintOptions>();
			options->lightning( true );


			utility::vector1< core::conformation::ResidueCOP > res_vector;
			SapConstraintHelperOP helper = common_setup( options, pose, sel, sel, sel, res_vector );
			Real helper_sap_score = helper->calculate_energy( res_vector, 0 );

			// need to check
			// 1. sap_score
			// 2. 1b table
			// 3. 2b table


			// Take a diversion and assert some really hard to debug stuff
			utility::vector1< core::conformation::ResidueCOP > res_vector2;
			pack::rotamer_set::RotamerSetsOP initial_rotsets = rotamer_sets_from_pose( pose, res_vector2 );

			SapConstraintHelperOP helper2 = utility::pointer::make_shared<SapConstraintHelper>( options );
			helper2->apply_residue_selectors( pose );
			pack::rotamer_set::RotamerSets symm_rotsets = helper2->setup_for_symmetry( pose, *initial_rotsets );
			pack::rotamer_set::RotamerSets & rotsets = using_symm ? symm_rotsets : *initial_rotsets;
			helper2->resize_arrays( pose );
			core::pack::rotamer_set::RotamerSets lean_rotamer_sets = helper2->setup_lightning( pose, rotsets );

			// First lets make sure the lean_rotamer_sets makes sense

			for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
				if ( pose.residue(seqpos).is_virtual_residue() ) continue;
				pack::rotamer_set::RotamerSetOP seq_rotamer_set = lean_rotamer_sets.rotamer_set_for_residue( seqpos );
				TS_ASSERT( seq_rotamer_set->num_rotamers() == 1 );

				TS_ASSERT( seq_rotamer_set->rotamer(1)->name() == pose.residue(seqpos).name() );

				Size last_atom = pose.residue(seqpos).natoms();
				TS_ASSERT_DELTA( seq_rotamer_set->rotamer(1)->xyz(last_atom).distance(pose.residue(seqpos).xyz(last_atom)), 0, 0.01 );
			}

			// Fast assumes the pose is poly-cys. So we have to continue using fast methods
			utility::vector1<utility::vector1<Real>> blocks = manually_calculate_fast_blocks( pose, lean_rotamer_sets );
			utility::vector1<utility::vector1<float>> sasa_scores = manually_calculate_fast_sasa_scores( pose, lean_rotamer_sets, blocks );

			// Next, we need to look at the sasa_scores

			assert_lightning_sasa_scores( helper2, lean_rotamer_sets, sasa_scores );


			// Sap score is easy, just do the clipping at the residue level

			utility::vector1<Size> irot_at_position;
			std::unordered_map< Size, core::conformation::ResidueCOP > seq100_aa_to_rotamer;
			for ( Size i = 1; i <= pose.size(); i++ ) {
				irot_at_position.push_back( i + lean_rotamer_sets.nrotamers() );
				seq100_aa_to_rotamer[ i * 100 + Size(pose.residue(i).aa()) ] = pose.residue(i).get_self_ptr();
			}

			Real sap_score = lightning_sap_from_sasa_scores_and_assert_1b_2b( pose, sasa_scores, irot_at_position, seq100_aa_to_rotamer, helper );

			TS_ASSERT_DELTA( helper_sap_score, sap_score, 0.1 );

			SapDatabase::get_instance()->symm_debug_force_map_ = false;
		}

	}


	Real
	calculate_lightning_sap(
		core::pose::Pose const & pose,
		utility::vector1<core::conformation::ResidueCOP> const & all_rotamers,
		core::pack::rotamer_set::RotamerSets const & lean_rotamer_sets,
		utility::vector1<Size> irot_at_position,
		std::unordered_map< Size, core::conformation::ResidueCOP > const & seq100_aa_to_rotamer,
		utility::vector1<utility::vector1<float>> const & sasa_scores, // const but whatever
		// SapConstraintHelperOP const & check_helper,
		SapConstraintHelperOP const & scratch_helper

	) {

		utility::vector1<Size> to_recalc;
		utility::vector1<core::conformation::ResidueCOP> my_resvect;

		pack::rotamer_set::RotamerSetOP fake_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();

		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

			Size irot = irot_at_position[ seqpos ];
			core::conformation::ResidueCOP const & rotamer = all_rotamers[irot]; //irot <= rotsets->nrotamers() ? rotsets->rotamer( irot )
			//: pose.residue( irot - rotsets->nrotamers() ).get_self_ptr();

			pack::rotamer_set::RotamerSetCOP lean_rotset = lean_rotamer_sets.has_rotamer_set_for_residue( seqpos ) ?
				lean_rotamer_sets.rotamer_set_for_residue( seqpos ) : fake_rotset;
			bool found = false;
			Size start_rot = lean_rotamer_sets.has_rotamer_set_for_residue( seqpos ) ? 1 : 0;
			for ( Size i = start_rot; i <= lean_rotset->num_rotamers(); i++ ) {
				core::conformation::ResidueCOP const & rotset_rot = i > 0 ? lean_rotset->rotamer(i) : pose.residue( seqpos ).get_self_ptr();
				if ( rotset_rot->aa() == rotamer->aa() ) {
					TS_ASSERT( ! found );
					found = true;
					my_resvect.push_back( rotset_rot );
				}
			}
			TS_ASSERT( found );

			assert( rotamer->natoms() == pose.residue(seqpos).natoms() );
		}

		Real other_method = lightning_sap_from_sasa_scores_and_assert_1b_2b( pose, sasa_scores, irot_at_position, seq100_aa_to_rotamer, scratch_helper );

		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			return other_method;    // for symmetry they do have to be perfect clones
		}

		// This works because they don't have to be perfect clones
		scratch_helper->calculate_energy( my_resvect, 0 );

		TS_ASSERT_DELTA( scratch_helper->current_score_, other_method, 0.2 );


		return scratch_helper->current_score_;
	}



	utility::vector1<utility::vector1<float>>
	expand_sasa_scores(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const & lean_rotamer_sets,
		core::pack::rotamer_set::RotamerSets const & rotsets,
		utility::vector1<utility::vector1<float>> const & sasa_scores,
		std::unordered_map< Size, core::conformation::ResidueCOP > & seq100_aa_to_rotamer
	) {
		utility::vector1<utility::vector1<float>> sasa_scores_for_full;

		std::unordered_map< Size, Size > aa_2_offset;
		Size cur_seq = 0;
		Size cur_rotamer = 1;
		Size cur_lean_rotamer = 1;
		while ( cur_rotamer <= rotsets.nrotamers() && cur_lean_rotamer <= lean_rotamer_sets.nrotamers() ) {

			aa_2_offset.clear();
			cur_seq ++;

			// Build a table for this seqpos
			while ( cur_lean_rotamer <= lean_rotamer_sets.nrotamers() && lean_rotamer_sets.rotamer( cur_lean_rotamer )->seqpos() == cur_seq ) {
				Size aa = lean_rotamer_sets.rotamer( cur_lean_rotamer )->aa();
				TS_ASSERT( aa_2_offset.count( aa ) == 0 );
				aa_2_offset[ aa ] = cur_lean_rotamer;
				seq100_aa_to_rotamer[ cur_seq * 100 + aa ] = lean_rotamer_sets.rotamer( cur_lean_rotamer );
				cur_lean_rotamer++;
			}

			while ( cur_rotamer <= rotsets.nrotamers() && rotsets.rotamer( cur_rotamer )->seqpos() == cur_seq ) {
				Size offset = aa_2_offset.at( Size(rotsets.rotamer( cur_rotamer )->aa()) );
				sasa_scores_for_full.push_back( sasa_scores[offset] );
				cur_rotamer++;
			}
		}

		TS_ASSERT( cur_rotamer > rotsets.nrotamers() );
		TS_ASSERT( cur_lean_rotamer > lean_rotamer_sets.nrotamers() );

		for ( Size i = 1; i <= pose.size(); i++ ) {
			sasa_scores_for_full.push_back( sasa_scores[ lean_rotamer_sets.nrotamers() + i ] );
		}

		return sasa_scores_for_full;
	}




	// void test_lightning_packing() {
	//     TR << "test_lightning_packing" << std::endl;
	//     for ( Size int_sym = 0; int_sym <= 1; int_sym++ ) {
	//         bool using_symm = bool(int_sym);
	//         TR << (using_symm ? "Symm pass" : "first pass") << std::endl;

	//         core::pose::Pose pose;
	//         utility::vector1< core::conformation::ResidueCOP > res_vector;
	//         pack::rotamer_set::RotamerSetsOP initial_rotsets;
	//         SapConstraintOptionsOP options;
	//         select::residue_selector::ResidueSelectorCOP sel;
	//         Real pack_size = 0;

	//         prepare_packing_pose_and_rotsets( pose, res_vector, initial_rotsets, options, sel, pack_size, using_symm, true );
	//         options->lightning( true );

	//         SapConstraintHelperOP scratch_helper = utility::pointer::make_shared<SapConstraintHelper>( options );
	//         scratch_helper->init_with_pose( pose, *initial_rotsets );

	//         SapConstraintHelperOP helper = utility::pointer::make_shared<SapConstraintHelper>( options );
	//         pack::rotamer_set::RotamerSets rotsets = helper->init_with_pose( pose, *initial_rotsets );

	//         utility::vector1<core::conformation::ResidueCOP> all_rotamers;
	//         for ( Size i = 1; i<= rotsets.nrotamers(); i++ ) {
	//             all_rotamers.push_back( rotsets.rotamer( i ) );
	//         }
	//         for ( Size i = 1; i <= pose.size(); i++ ) {
	//             all_rotamers.push_back( pose.residue(i).get_self_ptr() );
	//         }

	//         // Get the lean rotamer sets
	//         SapConstraintHelperOP helper2 = utility::pointer::make_shared<SapConstraintHelper>( options );
	//         helper2->apply_residue_selectors( pose );
	//         pack::rotamer_set::RotamerSets symm_rotsets2 = helper2->setup_for_symmetry( pose, rotsets );
	//         pack::rotamer_set::RotamerSets & rotsets2 = using_symm ? symm_rotsets2 : rotsets;
	//         helper2->resize_arrays( pose );
	//         core::pack::rotamer_set::RotamerSets lean_rotamer_sets = helper2->setup_lightning( pose, rotsets2 );

	//         // blocks have to be based on the lean rotamer sets
	//         utility::vector1<utility::vector1<Real>> blocks = manually_calculate_fast_blocks( pose, lean_rotamer_sets );
	//         utility::vector1<utility::vector1<float>> sasa_scores = manually_calculate_fast_sasa_scores( pose, lean_rotamer_sets, blocks );

	//         assert_lightning_sasa_scores( helper2, lean_rotamer_sets, sasa_scores );

	//         std::unordered_map< Size, core::conformation::ResidueCOP > seq100_aa_to_rotamer;
	//         utility::vector1<utility::vector1<float>> sasa_scores_for_full = expand_sasa_scores( pose, lean_rotamer_sets, rotsets, sasa_scores, seq100_aa_to_rotamer );

	//         // start actually testing

	//         Real init_score = helper->calculate_energy( res_vector, 0 );


	//         utility::vector1<Size> to_recalc;

	//         utility::vector1<Size> irot_at_position;
	//         for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
	//             irot_at_position.push_back( rotsets.nrotamers() + seqpos );
	//             to_recalc.push_back( seqpos );
	//         }



	//         TS_ASSERT_DELTA( init_score, calculate_lightning_sap( pose, all_rotamers, lean_rotamer_sets, irot_at_position, seq100_aa_to_rotamer, sasa_scores_for_full, scratch_helper ), 0.001 );


	//         // std::cout << "After" << std::endl;
	//         // random but repeatable
	//         utility::vector1<Size> order = { 9, 7, 14, 16, 6, 10, 13, 5, 17, 15, 8, 19, 2, 18, 12, 4, 1, 3, 11,
	//                                         13, 10, 15, 6, 18, 3, 8, 17, 1, 9, 4, 5, 16, 12, 19, 11, 14, 7, 2,
	//                                         9, 7, 4, 12, 6, 1, 3, 11, 19, 10, 13, 5, 2, 8, 17, 14, 18, 16, 15,
	//                                         11, 18, 6, 19, 14, 8, 10, 12, 3, 16, 5, 15, 7, 4, 1, 17, 2, 9, 13,
	//                                         5, 15, 9, 11, 14, 8, 17, 13, 18, 4, 3, 7, 1, 12, 2, 19, 6, 16, 10,
	//                                         2, 5, 7, 17, 14, 8, 6, 1, 9, 15, 3, 4, 19, 12, 11, 13, 18, 16, 10
	//         };


	//         for ( Size repeat = 1; repeat <= 3; repeat++ ) {

	//             for ( Size letter1 : order ) {
	//                 for ( core::Size seqpos = 1; seqpos <= pack_size; seqpos++ ) {


	//                     Size use_backwards = (repeat + letter1 + seqpos) % 2;

	//                     // std::cout << " " << repeat << " " << letter1 << " " << seqpos << std::endl;
	//                     core::conformation::ResidueCOP rotamer = rotsets.rotamer_set_for_residue( seqpos )->rotamer( (letter1-1) * 2 + use_backwards + 1 );
	//                     utility::vector1< core::conformation::ResidueCOP > old_res_vector = res_vector;
	//                     utility::vector1<Size> old_irot_at_position = irot_at_position;

	//                     irot_at_position[ seqpos ] = rotsets.nrotamer_offset_for_moltenres( rotsets.resid_2_moltenres( seqpos ) ) + (letter1-1) * 2 + use_backwards + 1;
	//                     assert( &*rotamer == &*rotsets.rotamer( irot_at_position[ seqpos ] ) );

	//                     // Do a repeatable random to get commit status
	//                     double temp = std::sqrt( repeat * letter1 * seqpos );
	//                     double deci = temp - (long)temp;
	//                     bool accept = deci > 0.5;
	//                     accept = true;

	//                     pose.replace_residue( seqpos, *rotamer, false );
	//                     res_vector[seqpos] = rotamer;

	//                     if ( pack_size != pose.size() ) {
	//                         for ( Size i = 1; i < SYMM; i++ ) {
	//                             Size sympos = seqpos + i * pack_size;
	//                             res_vector[ sympos ] = rotsets.rotamer_set_for_residue( sympos )->rotamer( (letter1-1) * 2 + use_backwards + 1 );

	//                             irot_at_position[ sympos ] = rotsets.nrotamer_offset_for_moltenres( rotsets.resid_2_moltenres( sympos ) ) + (letter1-1) * 2 + use_backwards + 1;
	//                             assert( &*res_vector[ sympos ] == &*rotsets.rotamer( irot_at_position[ sympos ] ) );
	//                         }
	//                     }

	//                     Real score = helper->calculate_energy( res_vector, seqpos );

	//                     // for ( Size i = 1; i<= pose.size(); i++ ) {
	//                     //     for ( Size j = 0; j< helper->dirty_sasa_[i].size(); j++ ) {
	//                     //         assert( helper->dirty_sasa_[i][j] );
	//                     //         helper->atom_sap_[i][j] = 0;
	//                     //     }
	//                     // }
	//                     // helper->current_score_ = 0;
	//                     // helper->recalculate_saps( to_recalc );

	//                     // std::cout << "new after" << std::endl;

	//                     score = helper->current_score_;


	//                     // for ( Size i = 1; i <= pose.size(); i++ ) {
	//                     //     core::conformation::Residue const & pose_res = pose.residue(i);
	//                     //     core::conformation::Residue const & vector_res = *res_vector[i];
	//                     //     Size irot = irot_at_position[ i ];
	//                     //     core::conformation::Residue const & irot_res = irot <= rotsets->nrotamers() ? *rotsets->rotamer( irot )
	//                     //                                                                         : pose.residue(irot - rotsets->nrotamers() );

	//                     //     std::cout << pose_res.natoms() << " " << vector_res.natoms() << " " << irot_res.natoms() << std::endl;
	//                     // }

	//                     // Different thresholds because we can't perfectly test the symmetry version
	//                     Real threshold = using_symm ? 0.5 : 0.01;

	//                     TS_ASSERT_DELTA( score, calculate_lightning_sap( pose, all_rotamers, lean_rotamer_sets, irot_at_position, seq100_aa_to_rotamer, sasa_scores_for_full, scratch_helper ), threshold );

	//                     if ( ! accept ) {
	//                         pose.replace_residue( seqpos, *old_res_vector[seqpos], false );
	//                         res_vector = old_res_vector;
	//                         irot_at_position = old_irot_at_position;
	//                     } else {
	//                         helper->commit_considered_substitution();
	//                     }

	//                 }
	//             }
	//         }
	//     }

	// }

	void test_lightning_packing_partial() {
		TR << "test_lightning_packing" << std::endl;
		for ( Size int_sym = 0; int_sym <= 2; int_sym++ ) {
			bool using_symm = bool(int_sym);
			bool using_map = int_sym == 2;
			TR << (using_symm ? "Symm pass" : "first pass") << std::endl;
			SapDatabase::get_instance()->symm_debug_force_map_ = using_map;

			core::pose::Pose pose;
			utility::vector1< core::conformation::ResidueCOP > res_vector;
			pack::rotamer_set::RotamerSetsOP initial_rotsets;
			SapConstraintOptionsOP options;
			select::residue_selector::ResidueSelectorCOP sel;
			select::residue_selector::ResidueSelectorCOP design_sel =
				utility::pointer::make_shared<select::residue_selector::ResidueIndexSelector>("1,3,5");
			// utility::pointer::make_shared<select::residue_selector::ResidueIndexSelector>("1-6");
			Real pack_size = 0;

			prepare_packing_pose_and_rotsets( pose, res_vector, initial_rotsets, options, sel, pack_size, using_symm, true, design_sel );

			core::select::residue_selector::ResidueSubset design_sub = design_sel->apply( pose );

			options->lightning( true );

			SapConstraintHelperOP scratch_helper = utility::pointer::make_shared<SapConstraintHelper>( options );
			scratch_helper->init_with_pose( pose, *initial_rotsets );

			SapConstraintHelperOP helper = utility::pointer::make_shared<SapConstraintHelper>( options );
			pack::rotamer_set::RotamerSets rotsets = helper->init_with_pose( pose, *initial_rotsets );

			utility::vector1<core::conformation::ResidueCOP> all_rotamers;
			for ( Size i = 1; i<= rotsets.nrotamers(); i++ ) {
				all_rotamers.push_back( rotsets.rotamer( i ) );
			}
			for ( Size i = 1; i <= pose.size(); i++ ) {
				all_rotamers.push_back( pose.residue(i).get_self_ptr() );
			}


			// Get the lean rotamer sets
			SapConstraintHelperOP helper2 = utility::pointer::make_shared<SapConstraintHelper>( options );
			helper2->apply_residue_selectors( pose );
			pack::rotamer_set::RotamerSets symm_rotsets2 = helper2->setup_for_symmetry( pose, rotsets );
			pack::rotamer_set::RotamerSets & rotsets2 = using_symm ? symm_rotsets2 : rotsets;
			helper2->resize_arrays( pose );
			core::pack::rotamer_set::RotamerSets lean_rotamer_sets = helper2->setup_lightning( pose, rotsets2 );

			// blocks have to be based on the lean rotamer sets
			utility::vector1<utility::vector1<Real>> blocks = manually_calculate_fast_blocks( pose, lean_rotamer_sets );
			utility::vector1<utility::vector1<float>> sasa_scores = manually_calculate_fast_sasa_scores( pose, lean_rotamer_sets, blocks );

			assert_lightning_sasa_scores( helper2, lean_rotamer_sets, sasa_scores );

			std::unordered_map< Size, core::conformation::ResidueCOP > seq100_aa_to_rotamer;
			utility::vector1<utility::vector1<float>> sasa_scores_for_full = expand_sasa_scores( pose, lean_rotamer_sets, rotsets, sasa_scores, seq100_aa_to_rotamer );


			// start actually testing


			Real init_score = helper->calculate_energy( res_vector, 0 );


			utility::vector1<Size> to_recalc;

			utility::vector1<Size> irot_at_position;
			for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
				irot_at_position.push_back( rotsets.nrotamers() + seqpos );
				to_recalc.push_back( seqpos );
			}



			TS_ASSERT_DELTA( init_score, calculate_lightning_sap( pose, all_rotamers, lean_rotamer_sets, irot_at_position, seq100_aa_to_rotamer, sasa_scores_for_full, scratch_helper ), 0.001 );


			// std::cout << "After" << std::endl;
			// random but repeatable
			utility::vector1<Size> order = { 9, 7, 14, 16, 6, 10, 13, 5, 17, 15, 8, 19, 2, 18, 12, 4, 1, 3, 11,
				// 13, 10, 15, 6, 18, 3, 8, 17, 1, 9, 4, 5, 16, 12, 19, 11, 14, 7, 2,
				// 9, 7, 4, 12, 6, 1, 3, 11, 19, 10, 13, 5, 2, 8, 17, 14, 18, 16, 15,
				// 11, 18, 6, 19, 14, 8, 10, 12, 3, 16, 5, 15, 7, 4, 1, 17, 2, 9, 13,
				// 5, 15, 9, 11, 14, 8, 17, 13, 18, 4, 3, 7, 1, 12, 2, 19, 6, 16, 10,
				2, 5, 7, 17, 14, 8, 6, 1, 9, 15, 3, 4, 19, 12, 11, 13, 18, 16, 10
				};



			for ( Size repeat = 1; repeat <= 1; repeat++ ) {

				for ( Size letter1 : order ) {
					for ( core::Size seqpos = 1; seqpos <= pack_size; seqpos++ ) {
						if ( ! design_sub[ seqpos ] ) continue;


						Size use_backwards = (repeat + letter1 + seqpos) % 2;

						// std::cout << " " << repeat << " " << letter1 << " " << seqpos << std::endl;
						core::conformation::ResidueCOP rotamer = rotsets.rotamer_set_for_residue( seqpos )->rotamer( (letter1-1) * 2 + use_backwards + 1 );
						utility::vector1< core::conformation::ResidueCOP > old_res_vector = res_vector;
						utility::vector1<Size> old_irot_at_position = irot_at_position;

						irot_at_position[ seqpos ] = rotsets.nrotamer_offset_for_moltenres( rotsets.resid_2_moltenres( seqpos ) ) + (letter1-1) * 2 + use_backwards + 1;
						assert( &*rotamer == &*rotsets.rotamer( irot_at_position[ seqpos ] ) );

						// Do a repeatable random to get commit status
						//double temp = std::sqrt( repeat * letter1 * seqpos );
						//double deci = temp - (long)temp;
						bool accept = true; //deci > 0.5;

						pose.replace_residue( seqpos, *rotamer, false );
						res_vector[seqpos] = rotamer;

						if ( pack_size != pose.size() ) {
							for ( Size i = 1; i < SYMM; i++ ) {
								Size sympos = seqpos + i * pack_size;
								res_vector[ sympos ] = rotsets.rotamer_set_for_residue( sympos )->rotamer( (letter1-1) * 2 + use_backwards + 1 );

								irot_at_position[ sympos ] = rotsets.nrotamer_offset_for_moltenres( rotsets.resid_2_moltenres( sympos ) ) + (letter1-1) * 2 + use_backwards + 1;
								assert( &*res_vector[ sympos ] == &*rotsets.rotamer( irot_at_position[ sympos ] ) );
							}
						}

						helper->calculate_energy( res_vector, seqpos );

						// for ( Size i = 1; i<= pose.size(); i++ ) {
						//     for ( Size j = 0; j< helper->dirty_sasa_[i].size(); j++ ) {
						//         assert( helper->dirty_sasa_[i][j] );
						//         helper->atom_sap_[i][j] = 0;
						//     }
						// }
						// helper->current_score_ = 0;
						// helper->recalculate_saps( to_recalc );

						// std::cout << "new after" << std::endl;

						Real score = helper->current_score_;


						// for ( Size i = 1; i <= pose.size(); i++ ) {
						//     core::conformation::Residue const & pose_res = pose.residue(i);
						//     core::conformation::Residue const & vector_res = *res_vector[i];
						//     Size irot = irot_at_position[ i ];
						//     core::conformation::Residue const & irot_res = irot <= rotsets->nrotamers() ? *rotsets->rotamer( irot )
						//                                                                         : pose.residue(irot - rotsets->nrotamers() );

						//     std::cout << pose_res.natoms() << " " << vector_res.natoms() << " " << irot_res.natoms() << std::endl;
						// }

						// Different thresholds because we can't perfectly test the symmetry version
						Real threshold = using_symm ? 0.5 : 0.01;

						TS_ASSERT_DELTA( score, calculate_lightning_sap( pose, all_rotamers, lean_rotamer_sets, irot_at_position, seq100_aa_to_rotamer, sasa_scores_for_full, scratch_helper ), threshold );

						if ( ! accept ) {
							pose.replace_residue( seqpos, *old_res_vector[seqpos], false );
							res_vector = old_res_vector;
							irot_at_position = old_irot_at_position;
						} else {
							helper->commit_considered_substitution();
						}
					}
				}
			}
			SapDatabase::get_instance()->symm_debug_force_map_ = false;
		}

	}



};
