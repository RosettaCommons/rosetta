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

#include <platform/types.hh>
#include <core/scoring/lkball/LK_DomeEnergy.hh>
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomeHelper.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <util/pose_funcs.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/residue_datacache.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/pack/packer_neighbors.hh>
#include <numeric/numeric.functions.hh>

#include <numeric/xyzVector.io.hh>
#include <numeric/xyzMatrix.io.hh>


#include <core/types.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>

#include <test/UTracer.hh>
#include <test/core/init_util.hh>

#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

#include <boost/format.hpp>

using namespace core;
using namespace core::scoring;
using namespace core::scoring::lkball;
using namespace core::pack;
using namespace core::pack::guidance_scoreterms::lk_dome;
using namespace core::conformation::residue_datacache;


static basic::Tracer TR("core.scoring.lkball.LK_DomeEnergyTests.cxxtest");

#define SYMM 3



class LK_DomeEnergyTests : public CxxTest::TestSuite
{


public:
	void setUp()
	{
		core_init();
	}

	void test_bit_functions() {


		utility::vector1<uint32_t> spare_data(10);

		bit_set( 0, &spare_data[1] );
		bit_set( 32*1 + 1, &spare_data[1] );
		bit_set( 32*2 + 2, &spare_data[1] );
		bit_set( 32*3 + 3, &spare_data[1] );
		bit_set( 32*4 + 4, &spare_data[1] );
		bit_set( 32*5 + 5, &spare_data[1] );
		bit_set( 32*6 + 6, &spare_data[1] );
		bit_set( 32*7 + 7, &spare_data[1] );
		bit_set( 32*8 + 30, &spare_data[1] );
		bit_set( 32*9 + 31, &spare_data[1] );

		TS_ASSERT_EQUALS( spare_data[1], 0x1 );
		TS_ASSERT_EQUALS( spare_data[2], 0x2 );
		TS_ASSERT_EQUALS( spare_data[3], 0x4 );
		TS_ASSERT_EQUALS( spare_data[4], 0x8 );
		TS_ASSERT_EQUALS( spare_data[5], 0x10 );
		TS_ASSERT_EQUALS( spare_data[6], 0x20 );
		TS_ASSERT_EQUALS( spare_data[7], 0x40 );
		TS_ASSERT_EQUALS( spare_data[8], 0x80 );
		TS_ASSERT_EQUALS( spare_data[9], 0x40000000 );
		TS_ASSERT_EQUALS( spare_data[10], 0x80000000 );


		utility::vector1<uint32_t> data { 0xF0F0F0F0, 0xAAAAAAAA };

		TS_ASSERT_EQUALS( bit_test( 31, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 30, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 29, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 28, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 27, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 26, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 25, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 24, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 1, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 0, &data[1]), false );

		TS_ASSERT_EQUALS( bit_test( 32 + 31, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 32 + 30, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 32 + 29, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 32 + 28, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 32 + 27, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 32 + 26, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 32 + 25, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 32 + 24, &data[1]), false );
		TS_ASSERT_EQUALS( bit_test( 32 + 1, &data[1]), true );
		TS_ASSERT_EQUALS( bit_test( 32 + 0, &data[1]), false );

		TS_ASSERT_EQUALS( popcount_until_bit( 0, &data[1]), 0 );
		TS_ASSERT_EQUALS( popcount_until_bit( 4, &data[1]), 0 );
		TS_ASSERT_EQUALS( popcount_until_bit( 8, &data[1]), 4 );
		TS_ASSERT_EQUALS( popcount_until_bit( 31, &data[1]), 15 );
		TS_ASSERT_EQUALS( popcount_until_bit( 32, &data[1]), 16 );
		TS_ASSERT_EQUALS( popcount_until_bit( 38, &data[1]), 19 );
		TS_ASSERT_EQUALS( popcount_until_bit( 63, &data[1]), 31 );
	}

	utility::vector1<utility::vector0<LKD_ResidueInfo>>
	get_lkd_infos(
		core::pose::Pose const & pose,
		pack::rotamer_set::RotamerSets const & rotamer_sets,
		LK_DomeEnergyCOP const & lk_dome
	) {

		pack::rotamer_set::RotamerSetOP fake_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();

		utility::vector1<utility::vector0<LKD_ResidueInfo>> infos(pose.size());

		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

			pack::rotamer_set::RotamerSetCOP const & rotset = rotamer_sets.has_rotamer_set_for_residue( seqpos ) ?
				rotamer_sets.rotamer_set_for_residue( seqpos )
				: fake_rotset;

			for ( Size irot = 0; irot <= rotset->num_rotamers(); irot++ ) {

				conformation::Residue const & res = irot > 0 ? *rotset->rotamer(irot) : pose.residue( seqpos );

				infos[seqpos].emplace_back( res, &*lk_dome );
			}
		}
		return infos;
	}

	std::map<std::pair<Size, MyWaterHolder>, Size>
	get_water_map( LK_DomeHelper const & helper ) {
		utility::vector1<utility::vector0<MyWaterHolder>> const & waters = *helper.debug_waters_;

		std::map<std::pair<Size, MyWaterHolder>, Size> water_map;

		for ( Size seqpos = 1; seqpos <= waters.size(); seqpos++ ) {
			for ( Size water_id = 0; water_id < waters[seqpos].size(); water_id++ ) {

				std::pair<Size, MyWaterHolder> key( seqpos, waters[seqpos][water_id] );

				water_map[key] = water_id;
			}
		}
		return water_map;
	}

	Size
	iter_pwp_till_seqpos(
		LK_DomeHelper const & helper,
		Size pwp_cursor,
		Size pwp_elements,
		Size other_seqpos
	) {

		for ( Size ielem = 1; ielem <= pwp_elements; ielem++ ) {
			Size seqpos = helper.pos_water_pair_to_score_[pwp_cursor + PWP_OTHER_SEQPOS];
			if ( seqpos == other_seqpos ) return pwp_cursor;

			Size num_bitmasks = helper.seqpos_num_rotamers_32_[ seqpos ];
			pwp_cursor += PWP_BITFIELDS + num_bitmasks;
		}

		pwp_cursor = 0;
		return pwp_cursor;
	}

	std::pair<int8_t, int16_t>
	fetch_score_data(
		LK_DomeHelper const & helper,
		Size water_seqpos,
		Size water_id,
		Size other_seqpos,
		Size other_rotamer,
		bool both_ways,
		Size & offset
	) {

		std::pair<int8_t, int16_t> result( 0, 0 );

		Size score_data_offset = 0;
		{
			Size major_index = helper.seqpos_to_water_to_pwp_index_list_[water_seqpos];
			Size minor_index = water_id * 2;

			Size pwp_index = helper.water_to_pwp_index_list_[ major_index + minor_index ];
			Size pwp_elements = helper.water_to_pwp_index_list_[ major_index + minor_index + 1 ];

			Size pwp_cursor = iter_pwp_till_seqpos( helper, pwp_index, pwp_elements, other_seqpos );

			if ( pwp_cursor != 0 ) {
				if ( bit_test( other_rotamer, &helper.pos_water_pair_to_score_[pwp_cursor + PWP_BITFIELDS] ) ) {

					Size starting_offset = helper.pos_water_pair_to_score_[pwp_cursor + PWP_OFFSET_INTO_SCORE_DATA];
					Size popcount = popcount_until_bit( other_rotamer, &helper.pos_water_pair_to_score_[pwp_cursor + PWP_BITFIELDS]);
					score_data_offset = starting_offset + popcount * 3;

					result.first = helper.score_data_[score_data_offset];
					result.second = *((int16_t *)&helper.score_data_[score_data_offset+1]);
				}
			}
		}
		if ( both_ways ) {
			Size major_index = helper.seqpos_to_riww_at_pos_index_list_[other_seqpos];
			Size minor_index = other_rotamer * 2;

			Size riww_cursor = helper.riww_at_pos_index_list_[major_index + minor_index];
			Size riww_elements = helper.riww_at_pos_index_list_[major_index + minor_index + 1];

			bool found_interaction = false;
			for ( Size ielement = 1; ielement <= riww_elements; ielement ++ ) {
				Size seqpos = helper.rotamer_interacts_with_waters_at_pos_[ riww_cursor ];
				Size num_waters_32 = ((helper.seqpos_num_unique_waters_[seqpos]-1) >> 5) + 1;
				if ( seqpos == water_seqpos ) {
					found_interaction = bit_test( water_id, &helper.rotamer_interacts_with_waters_at_pos_[ riww_cursor + 1 ] );
					break;
				}
				riww_cursor += 1 + num_waters_32;
			}

			TS_ASSERT_EQUALS( found_interaction, score_data_offset != 0 );
		}

		offset = score_data_offset;
		return result;
	}



	// super duper slow because we're not taking shortcuts
	void assert_score_data(
		core::pose::Pose const & pose,
		ScoreFunction const & sfxn,
		LK_DomeHelper const & helper,
		pack::rotamer_set::RotamerSets const & rotamer_sets,
		utility::vector1<utility::vector0<LKD_ResidueInfo>> const & lkd_infos,
		std::map<std::pair<Size, MyWaterHolder>, Size> const & water_map
	) {


		Real occlude_inflater = helper.lk_dome()->occlusion_max_ / 127;

		pack::rotamer_set::RotamerSetOP fake_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();

		for ( Size water_seqpos = 1; water_seqpos <= pose.size(); water_seqpos++ ) {

			pack::rotamer_set::RotamerSetCOP const & water_rotset = rotamer_sets.has_rotamer_set_for_residue( water_seqpos ) ?
				rotamer_sets.rotamer_set_for_residue( water_seqpos )
				: fake_rotset;

			for ( Size wat_rot = 0; wat_rot <= water_rotset->num_rotamers(); wat_rot++ ) {

				conformation::Residue const & wat_res = wat_rot > 0 ? *water_rotset->rotamer(wat_rot) : pose.residue( water_seqpos );
				LKD_ResidueInfo const & wat_info = lkd_infos[water_seqpos][wat_rot];


				for ( Size wat_atom = 1; wat_atom <= wat_res.nheavyatoms(); wat_atom++ ) {
					Vector base_atom_xyz = wat_res.xyz( wat_atom );
					Size wat_offset = wat_info.water_offset_for_atom()[wat_atom];

					for ( Size wat_iwat = 1; wat_iwat <= wat_info.n_attached_waters()[wat_atom]; wat_iwat++ ) {

						Vector wat_wat_xyz = wat_info.waters()[wat_iwat + wat_offset];
						Real sol_value = wat_info.water_sol_values()[wat_atom];

						MyWaterHolder holder( wat_res.atom_is_backbone(wat_atom), sol_value, wat_wat_xyz, base_atom_xyz );
						Size water_id = water_map.at( std::pair<Size, MyWaterHolder>( water_seqpos, holder ) );

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						for ( Size other_seqpos = 1; other_seqpos <= pose.size(); other_seqpos++ ) {
							if ( other_seqpos == water_seqpos ) continue;

							pack::rotamer_set::RotamerSetCOP const & other_rotset = rotamer_sets.has_rotamer_set_for_residue( other_seqpos ) ?
								rotamer_sets.rotamer_set_for_residue( other_seqpos )
								: fake_rotset;

							for ( Size other_rot = 0; other_rot <= other_rotset->num_rotamers(); other_rot++ ) {

								conformation::Residue const & other_res = other_rot > 0 ? *other_rotset->rotamer(other_rot) : pose.residue( other_seqpos );
								LKD_ResidueInfo const & other_info = lkd_infos[other_seqpos][other_rot];

								Real interact_sum = 0;
								Real occlude_sum = 0;

								for ( Size other_atom = 1; other_atom <= other_res.nheavyatoms(); other_atom++ ) {

									Vector const & other_atom_xyz = other_res.xyz(other_atom);
									Size other_atom_type_index = other_res.atom_type_index(other_atom);

									occlude_sum += helper.lk_dome()->lk_ball()->get_lk_fractional_contribution_for_single_water(
										other_atom_xyz, other_atom_type_index, wat_wat_xyz );

									Real lk_dome_iso_frac = 0;
									Real lk_dome_frac = 0;
									helper.lk_dome()->single_water_atom_fractions( base_atom_xyz, wat_wat_xyz, other_atom_xyz, other_atom_type_index,
										lk_dome_iso_frac, lk_dome_frac );

									// Real saved = interact_sum;
									interact_sum += -1 * lk_dome_iso_frac * sol_value * sfxn.weights()[core::scoring::lk_dome_iso];
									interact_sum += -1 * lk_dome_frac     * sol_value * sfxn.weights()[core::scoring::lk_dome];



									// if ( wat_rot == 0 && other_rot == 0 ) {

									//  std::cout << "D " << water_seqpos << " " << other_seqpos << " " << wat_atom << " " << other_atom << " " << wat_iwat << " " << interact_sum - saved << std::endl;
									// }
									Size other_offset = other_info.water_offset_for_atom()[other_atom];

									for ( Size other_iwat = 1; other_iwat <= other_info.n_attached_waters()[other_atom]; other_iwat++ ) {

										Vector other_wat_xyz = other_info.waters()[other_iwat + other_offset];

										Real bridge_frac = helper.lk_dome()->single_water_water_fraction_1way( base_atom_xyz, wat_wat_xyz, other_wat_xyz );

										interact_sum += bridge_frac * sol_value * sfxn.weights()[core::scoring::lk_dome_bridge];
										interact_sum += bridge_frac * sfxn.weights()[core::scoring::lk_dome_bridge_uncpl];


										Real ball_bridge2_frac = helper.lk_dome()->single_water_water_bridge2_fraction_1way( wat_wat_xyz, other_wat_xyz );

										interact_sum += ball_bridge2_frac * sol_value * sfxn.weights()[core::scoring::lk_ball_bridge2];
										interact_sum += ball_bridge2_frac * sfxn.weights()[core::scoring::lk_ball_bridge_uncpl2];
									}
								}

								occlude_sum = std::min( occlude_sum, helper.lk_dome()->occlusion_max_ );

								Size offset = 0;
								std::pair<int8_t, int16_t> score_data_values = fetch_score_data( helper, water_seqpos, water_id, other_seqpos,
									other_rot, true, offset);


								Real inflated_occlude = score_data_values.first * occlude_inflater;
								Real inflated_interact = score_data_values.second * helper.partial_score_scaler_;
								(void)inflated_occlude;
								(void)inflated_interact;

								// std::cout << interact_sum << " " << inflated_interact << " -- " << offset << std::endl;
								// }
								TS_ASSERT_DELTA( occlude_sum, inflated_occlude, occlude_inflater );
								TS_ASSERT_DELTA( interact_sum, inflated_interact, helper.partial_score_scaler_*2 );



							}
						}


						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					}
				}
			}
		}
	}


	rotamer_set::RotamerSetsOP
	trpcage_rotsets(
		core::pose::Pose & pose
	) {

		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
		sfxn->score(pose);
		// This is just trpcage sequence + 2 ascending letters
		// All 20 aas are covered
		utility::vector1<std::string> allowed_aas {
			"NAC",
			"LDE",
			"YFG",
			"IHI",
			"QKL",
			"WMN",
			"LPW",
			"KRS",
			"DTV",
			"GWY",
			"GAC",
			"PDE",
			"SFG",
			"SHI",
			"GKL",
			"RMN",
			"PPW",
			"PRS",
			"PTV",
			"SWY"
			};

		core::pack::task::TaskFactoryOP tf = utility::pointer::make_shared<core::pack::task::TaskFactory>();
		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
			core::pack::task::operation::RestrictAbsentCanonicalAASOP restr =
				utility::pointer::make_shared<core::pack::task::operation::RestrictAbsentCanonicalAAS>();
			restr->keep_aas( allowed_aas[seqpos] );
			restr->include_residue( seqpos );
			tf->push_back( restr );
		}
		pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );
		core::pack::interaction_graph::AnnealableGraphBaseOP ig = nullptr;
		rotamer_set::RotamerSetsOP rotsets = utility::pointer::make_shared<rotamer_set::RotamerSets>();

		pack_rotamers_setup( pose, *sfxn, task, rotsets, ig );

		return rotsets;
	}


	void assert_correct_waters(
		LK_DomeHelper const & helper,
		utility::vector1< Size > const & current_rotamers,
		core::pose::Pose const & pose,
		rotamer_set::RotamerSets const & rotamer_sets,
		std::map<std::pair<Size, MyWaterHolder>, Size> const & water_map,
		utility::vector1<utility::vector0<LKD_ResidueInfo>> const & lkd_infos
	) {

		pack::rotamer_set::RotamerSetOP fake_rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();

		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

			pack::rotamer_set::RotamerSetCOP const & rotset = rotamer_sets.has_rotamer_set_for_residue( seqpos ) ?
				rotamer_sets.rotamer_set_for_residue( seqpos )
				: fake_rotset;

			Size irot = current_rotamers[seqpos];

			conformation::Residue const & res = irot > 0 ? *rotset->rotamer(irot) : pose.residue( seqpos );

			std::map<Size, bool> water_accounted_for;
			for ( Size iwat = 0; iwat < helper.waters_per_seqpos_; iwat++ ) {
				int water_id = helper.current_waters_[ (seqpos-1)*helper.waters_per_seqpos_ + 1 + iwat ];
				if ( water_id == -1 ) continue;
				TS_ASSERT( water_accounted_for.count(water_id) == 0 );
				water_accounted_for[water_id] = false;
			}

			LKD_ResidueInfo const & info = lkd_infos[seqpos][irot];

			for ( Size iatom = 1; iatom <= res.nheavyatoms(); iatom++ ) {
				Vector base_atom_xyz = res.xyz( iatom );
				Size atom_offset = info.water_offset_for_atom()[iatom];

				for ( Size iwat = 1; iwat <= info.n_attached_waters()[iatom]; iwat++ ) {

					Vector wat_wat_xyz = info.waters()[iwat + atom_offset];
					Real sol_value = info.water_sol_values()[iatom];

					MyWaterHolder holder( res.atom_is_backbone(iatom), sol_value, wat_wat_xyz, base_atom_xyz );
					Size water_id = water_map.at( std::pair<Size, MyWaterHolder>( seqpos, holder ) );

					TS_ASSERT( water_accounted_for.count( water_id ) == 1 );
					if ( water_accounted_for.count( water_id ) == 1 ) {
						TS_ASSERT( !water_accounted_for.at( water_id ) );
						water_accounted_for[water_id] = true;
					}
				}
			}

			for ( auto const & pair : water_accounted_for ) {
				TS_ASSERT( pair.second );
			}
		}
	}

	void assert_correct_scores(
		LK_DomeHelper const & helper,
		utility::vector1< Size > const & current_rotamers
	) {

		Size _;

		for ( Size water_seqpos = 1; water_seqpos <= current_rotamers.size(); water_seqpos++ ) {

			Size water_offset = (water_seqpos-1)*helper.waters_per_seqpos_ + 1;

			for ( Size iwat = 0; iwat < helper.waters_per_seqpos_; iwat++ ) {
				int water_id = helper.current_waters_[ water_offset + iwat ];

				int32_t occlude_sum = 0;
				int32_t interact_sum = 0;

				if ( water_id != -1 ) {
					for ( Size other_seqpos = 1; other_seqpos <= current_rotamers.size(); other_seqpos++ ) {

						Size other_irot = current_rotamers[other_seqpos];

						// This is so inefficient
						std::pair<int8_t, int16_t> score_data = fetch_score_data( helper, water_seqpos, water_id, other_seqpos, other_irot, false , _ );

						occlude_sum += score_data.first;
						interact_sum += score_data.second;
					}
				}

				// std::cout << "Position " << water_seqpos << " " << helper.current_water_interact_sum_[ water_offset + iwat ] << " " <<
				// interact_sum << std::endl;
				TS_ASSERT_EQUALS( helper.current_water_occl_sum_[ water_offset + iwat ], occlude_sum );
				TS_ASSERT_EQUALS( helper.current_water_interact_sum_[ water_offset + iwat ], interact_sum );
			}
		}
	}

	void
	assert_score_sum(
		LK_DomeHelper const & helper,
		Real reported_score
	) {

		int32_t score_sum = 0;
		for ( Size seqpos = 1; seqpos <= helper.current_rotamer_at_seqpos_.size(); seqpos++ ) {

			Size water_offset = (seqpos-1)*helper.waters_per_seqpos_ + 1;

			int32_t seqpos_score_sum = 0;
			for ( Size iwat = 0; iwat < helper.waters_per_seqpos_; iwat++ ) {

				int32_t occlude = helper.current_water_occl_sum_[ water_offset + iwat ];
				int32_t interact = helper.current_water_interact_sum_[ water_offset + iwat ];

				int32_t available;
				if ( occlude >= 127 ) {
					available = 0;
				} else {
					if ( occlude > helper.occlusion_min_ ) {
						available = helper.occlusion_span_ - (occlude - helper.occlusion_min_);
					} else {
						available = helper.occlusion_span_;
					}
				}

				seqpos_score_sum += interact * available;
			}

			TS_ASSERT_EQUALS( helper.current_score_at_seqpos_[seqpos], seqpos_score_sum );

			score_sum += seqpos_score_sum;
		}

		TS_ASSERT_EQUALS( helper.total_score_, score_sum );

		TS_ASSERT_DELTA( reported_score, score_sum * helper.score_scaler_, helper.score_scaler_ );

	}

	// void test_single_water_iso_derivs() {

	//  core::pose::Pose pose;
	//  core::pose::make_pose_from_sequence( pose, "KF", "fa_standard" );

	//  pose.set_chi(1, 1, 180);
	//  pose.set_chi(2, 1, 180);
	//  pose.set_chi(3, 1, 180);
	//  pose.set_chi(4, 1, 180);


	//  ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
	//  sfxn->set_weight(core::scoring::lk_dome_iso, 1);

	//  sfxn->score(pose);

	//  LK_DomeEnergy lk_dome(sfxn->energy_method_options());

	//  conformation::Residue const & lys = pose.residue(1);
	//  conformation::Residue const & phe = pose.residue(2);

	//  // Creating a copy
	//  LKD_ResidueInfo lys_info = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( 1 ).get_raw_ptr( LK_DOME_INFO )));
	//  LKD_ResidueInfo phe_info = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( 2 ).get_raw_ptr( LK_DOME_INFO )));

	//  Vector NZ = lys.xyz("NZ");
	//  Size NZ_no = lys.atom_index("NZ");
	//  Size NZ_type_index = lys.atom_type_index(NZ_no);
	//  Vector _3HZ = lys.xyz("3HZ");

	//  Vector unit = (_3HZ - NZ).normalized();
	//  Vector suspected_water_pos = NZ + unit * 2.65;

	//  Size CZ_no = phe.atom_index("CZ");
	//  Size CZ_type_index = phe.atom_type_index(CZ_no);

	//  bool found = false;
	//  Size water_iatom = 0;
	//  Size water_iwat = 0;
	//  for ( Size iatom = 1; iatom <= lys.nheavyatoms(); iatom++ ) {
	//   for ( Size iwat = 1; iwat <= lys_info.n_attached_waters()[iatom]; iwat++ ) {
	//    Real dist = suspected_water_pos.distance( lys_info.waters()[iatom][iwat] );
	//    if ( dist < 0.1 ) {
	//     found = true;
	//     water_iatom = iatom;
	//     water_iwat = iwat;
	//    } else {
	//     lys_info.water_occlusions()[iatom][iwat] = lk_dome.occlusion_scale();
	//    }
	//   }
	//  }

	//  runtime_assert(found);
	//  Vector water_xyz = lys_info.waters()[water_iatom][water_iwat];

	//  etable::count_pair::CountPairFunctionOP cp = core::scoring::etable::count_pair::CountPairFactory::create_count_pair_function(
	//   lys, phe, etable::count_pair::CP_CROSSOVER_4 );


	//  pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * 3.65 );


	//  EnergyMap emap;
	//  evaluate_lk_dome_energy_for_atom_ranges( lk_dome, lys, lys_info, phe, phe_info,
	//   *sfxn, *cp, NZ_no, NZ_no, CZ_no, CZ_no, emap );


	//  EnergyMap emap_manual;
	//  lk_dome.accumulate_single_atom_contributions( NZ_no, NZ_type_index, lys_info.n_attached_waters()[NZ_no],
	//   lys_info.waters()[NZ_no], lys_info.water_occlusions()[NZ_no], lys,
	//   CZ_type_index, phe.xyz("CZ"), lys_info.water_sol_values()[NZ_no], emap_manual );

	//  TS_ASSERT_EQUALS( emap[core::scoring::lk_dome_iso], emap_manual[core::scoring::lk_dome_iso] );


	//  Real deriv_size = 0.001;
	//  for ( Real dist = 1; dist < 10; dist += 0.05 ) {
	//   if ( std::abs<Real>(dist - 2.65) < deriv_size ) {
	//    continue;
	//   }


	//   pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * dist );

	//   utility::vector1< DerivVectorPair > r1_at_derivs( lys.natoms() );
	//   utility::vector1< DerivVectorPair > r2_at_derivs( phe.natoms() );

	//   lk_dome.sum_deriv_contributions_for_heavyatom_pair_one_way( NZ_no, lys, lys_info, CZ_no, phe, phe_info,
	//    sfxn->weights(), 1, NZ.distance_squared(phe.xyz("CZ")), r1_at_derivs, r2_at_derivs );

	//   Vector f2w = water_xyz - phe.xyz("CZ");
	//   Real r = f2w.norm();
	//   Vector reported_f2 = r2_at_derivs[CZ_no].f2();
	//   Real deriv = - reported_f2[0] / f2w[0] * r;

	//   Real direction = (dist < 2.65 ? -1 : 1 );


	//   Real low_score;
	//   Real center_score;
	//   Real high_score;
	//   {
	//    EnergyMap em;
	//    lk_dome.accumulate_single_atom_contributions( NZ_no, NZ_type_index, lys_info.n_attached_waters()[NZ_no],
	//     lys_info.waters()[NZ_no], lys_info.water_occlusions()[NZ_no], lys,
	//     CZ_type_index, phe.xyz("CZ"), lys_info.water_sol_values()[NZ_no], em );
	//    center_score = em[core::scoring::lk_dome_iso];
	//   }
	//   {
	//    pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * (dist-deriv_size*direction) );

	//    EnergyMap em;
	//    lk_dome.accumulate_single_atom_contributions( NZ_no, NZ_type_index, lys_info.n_attached_waters()[NZ_no],
	//     lys_info.waters()[NZ_no], lys_info.water_occlusions()[NZ_no], lys,
	//     CZ_type_index, phe.xyz("CZ"), lys_info.water_sol_values()[NZ_no], em );
	//    low_score = em[core::scoring::lk_dome_iso];
	//   }
	//   {
	//    pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * (dist+deriv_size*direction) );

	//    EnergyMap em;
	//    lk_dome.accumulate_single_atom_contributions( NZ_no, NZ_type_index, lys_info.n_attached_waters()[NZ_no],
	//     lys_info.waters()[NZ_no], lys_info.water_occlusions()[NZ_no], lys,
	//     CZ_type_index, phe.xyz("CZ"), lys_info.water_sol_values()[NZ_no], em );
	//    high_score = em[core::scoring::lk_dome_iso];
	//   }

	//   Real low_deriv = (center_score - low_score) / deriv_size;
	//   Real high_deriv = (high_score - center_score) / deriv_size;

	//   TS_ASSERT_DELTA( (high_deriv + low_deriv)/2.0, deriv, std::abs(deriv) / 100 );

	//   // std::cout << "" << std::endl;
	//   // std::cout << "Distance: " << dist << std::endl;
	//   // std::cout << low_deriv << std::endl;
	//   // std::cout << deriv << std::endl;
	//   // std::cout << high_deriv << std::endl;




	//  }



	// }

	void
	lkdome_all_pairwise(
		LK_DomeEnergy const & lk_dome,
		core::pose::Pose & pose,
		core::scoring::ScoreFunction const & sfxn,
		EnergyMap & emap
	) {

		EnergyGraph & energy_graph( pose.energies().energy_graph() );

		for ( Size i=1, i_end = pose.size(); i<= i_end; ++i ) {
			conformation::Residue const & resl( pose.residue( i ) );

			// LKD_ResidueInfo & info_l = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( i ).get_raw_ptr( LK_DOME_INFO )));

			for ( utility::graph::Graph::EdgeListIter
					iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->upper_edge_list_end();
					iru != irue; ++iru ) {
				auto & edge( static_cast< EnergyEdge & > (**iru) );

				Size const j( edge.get_second_node_ind() );
				conformation::Residue const & resu( pose.residue( j ) );

				lk_dome.residue_pair_energy( resl, resu, pose, sfxn, emap );

			}

		}

	}


	// void test_single_dofs()
	// {
	//  core::pose::Pose pose = create_trpcage_ideal_pose();
	//  core::scoring::ScoreFunction sfxn;
	//  // core::scoring::methods::EnergyMethodOptions options;
	//  // options.analytic_etable_evaluation( false );
	//  // sfxn.set_energy_method_options( options );
	//  // sfxn.set_weight( fa_atr, 0.5 );
	//  // sfxn.set_weight( fa_rep, 0.25 );
	//  // sfxn.set_weight( fa_sol, 0.125 );
	//  sfxn.set_weight( lk_dome_iso, 0.625 );


	//  sfxn(pose);

	//  LK_DomeEnergy lk_dome(sfxn.energy_method_options());
	//  core::optimization::MinimizerMap base;
	//  lk_dome.setup_for_minimizing( pose, sfxn, base);

	//  EnergyMap start_emap;
	//  lk_dome.setup_for_scoring(pose, sfxn);
	//  lkdome_all_pairwise(lk_dome, pose, sfxn, start_emap);
	//  Real start_score = start_emap[core::scoring::lk_dome_iso];

	//  Real start_chi = pose.chi(2, 20);

	//  pose.set_chi(2, 20, start_chi - 20);
	//  EnergyMap low_emap;
	//  lk_dome.setup_for_scoring(pose, sfxn);
	//  lkdome_all_pairwise(lk_dome, pose, sfxn, low_emap);
	//  Real low_score = low_emap[core::scoring::lk_dome_iso];

	//  pose.set_chi(2, 20, start_chi + 20);
	//  EnergyMap high_emap;
	//  lk_dome.setup_for_scoring(pose, sfxn);
	//  lkdome_all_pairwise(lk_dome, pose, sfxn, high_emap);
	//  Real high_score = low_emap[core::scoring::lk_dome_iso];

	//  std::cout << "ASDASD" << std::endl;
	//  std::cout << low_score << std::endl;
	//  std::cout << start_score << std::endl;
	//  std::cout << high_score << std::endl;

	// }

	void test_simple_derivs() {

		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
		sfxn->set_weight(core::scoring::lk_dome, 1);

		LK_DomeEnergy lk_dome(sfxn->energy_method_options());

		Real type_index = 6;
		Real deriv_size = 0.0001;

		for ( Real dist = 0.05; dist < 10; dist += 0.05 ) {

			Real dist2 = dist*dist;

			Real d2_low = lk_dome.d2_low_[ type_index ];
			Real d2_delta = dist2 - d2_low;

			Real ddome_frac_dr = 0;
			if ( d2_delta < lk_dome.ramp_width_A2_ && d2_delta > 0  ) {
				ddome_frac_dr = lk_dome.eval_d_lk_fraction_dr_over_r( d2_delta, lk_dome.ramp_width_A2_ ) * dist;
			}

			Real diso_frac_dr =  lk_dome.dlk_frac_sigmoid( dist );


			Real iso_frac = lk_dome.lk_frac_sigmoid( dist );
			Real dome_frac = lk_dome.get_lkd_frac_dome_dist( dist2, type_index );

			Real low_iso_frac = lk_dome.lk_frac_sigmoid( dist - deriv_size );
			Real low_dome_frac = lk_dome.get_lkd_frac_dome_dist( numeric::square( dist - deriv_size ), type_index );

			Real high_iso_frac = lk_dome.lk_frac_sigmoid( dist + deriv_size );
			Real high_dome_frac = lk_dome.get_lkd_frac_dome_dist( numeric::square( dist + deriv_size ), type_index );


			Real iso_low_deriv = (iso_frac - low_iso_frac) / deriv_size;
			Real iso_high_deriv = (high_iso_frac - iso_frac) / deriv_size;

			// std::cout << "Iso: " << iso_low_deriv << std::endl;
			// std::cout << "Iso: " << diso_frac_dr << std::endl;
			// std::cout << "Iso: " << iso_high_deriv << std::endl;
			TS_ASSERT_DELTA( (iso_high_deriv + iso_low_deriv) / 2, diso_frac_dr, std::abs<Real>(diso_frac_dr / 100) );


			Real dome_low_deriv = (dome_frac - low_dome_frac) / deriv_size;
			Real dome_high_deriv = (high_dome_frac - dome_frac) / deriv_size;

			// std::cout << "Dome: " << dome_low_deriv << std::endl;
			// std::cout << "Dome: " << ddome_frac_dr << std::endl;
			// std::cout << "Dome: " << dome_high_deriv << std::endl;
			TS_ASSERT_DELTA( (dome_high_deriv + dome_low_deriv) / 2, ddome_frac_dr, std::abs<Real>(ddome_frac_dr / 100) );

		}



	}



	void test_dlk_dome_E_dR() {

		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "KF", "fa_standard" );

		pose.set_chi(1, 1, 180);
		pose.set_chi(2, 1, 180);
		pose.set_chi(3, 1, 180);
		pose.set_chi(4, 1, 180);


		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
		sfxn->set_weight(core::scoring::lk_dome, 1);

		sfxn->score(pose);

		LK_DomeEnergy lk_dome(sfxn->energy_method_options());

		conformation::Residue const & lys = pose.residue(1);
		conformation::Residue const & phe = pose.residue(2);

		// Creating a copy
		LKD_ResidueInfo lys_info = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( 1 ).get_raw_ptr( LK_DOME_INFO )));
		LKD_ResidueInfo phe_info = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( 2 ).get_raw_ptr( LK_DOME_INFO )));

		Vector NZ = lys.xyz("NZ");
		Size NZ_no = lys.atom_index("NZ");
		Size NZ_type_index = lys.atom_type_index(NZ_no);
		Vector _3HZ = lys.xyz("3HZ");

		Vector unit = (_3HZ - NZ).normalized();
		Vector suspected_water_pos = NZ + unit * lk_dome.w_dist();

		Size CZ_no = phe.atom_index("CZ");
		Size CZ_type_index = phe.atom_type_index(CZ_no);

		bool found = false;
		Size water_iatom = 0;
		Size water_iwat = 0;
		Size water_offset = 0;
		for ( Size iatom = 1; iatom <= lys.nheavyatoms(); iatom++ ) {
			Size offset = lys_info.water_offset_for_atom()[iatom];
			for ( Size iwat = 1; iwat <= lys_info.n_attached_waters()[iatom]; iwat++ ) {
				Real dist = suspected_water_pos.distance( lys_info.waters()[iwat + offset] );
				if ( dist < 0.1 ) {
					found = true;
					water_iatom = iatom;
					water_iwat = iwat;
					water_offset = offset;
				} else {
					lys_info.water_occlusions()[iatom][iwat] = lk_dome.occlusion_max();
				}
			}
		}

		runtime_assert(found);
		Vector water_xyz = lys_info.waters()[water_iwat + water_offset];



		Real sol_value = lys_info.water_sol_values()[water_iatom];
		Real avail = lk_dome.get_avail( lys_info.water_occlusions()[water_iatom][water_iwat] );





		etable::count_pair::CountPairFunctionOP cp = core::scoring::etable::count_pair::CountPairFactory::create_count_pair_function(
			lys, phe, etable::count_pair::CP_CROSSOVER_4 );


		pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * 3.65 );


		EnergyMap emap;
		evaluate_lk_dome_energy_for_atom_ranges( lk_dome, lys, lys_info, phe, phe_info,
			*sfxn, *cp, NZ_no, NZ_no, CZ_no, CZ_no, emap );


		EnergyMap emap_manual;
		lk_dome.accumulate_single_atom_contributions( NZ_no, NZ_type_index, lys_info.n_attached_waters()[NZ_no], lys_info.water_offset_for_atom()[NZ_no],
			lys_info.waters(), lys_info.water_occlusions()[NZ_no], lys,
			CZ_type_index, phe.xyz("CZ"), lys_info.water_sol_values()[NZ_no], emap_manual );

		TS_ASSERT_EQUALS( emap[core::scoring::lk_dome], emap_manual[core::scoring::lk_dome] );

		unit += Vector(0.3, 0, 0);
		unit = unit.normalized();

		Real deriv_size = 0.00001;
		for ( Real dist = 1; dist < 10; dist += 0.05 ) {
			if ( std::abs<Real>(dist - lk_dome.w_dist()) < deriv_size ) {
				continue;
			}



			pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * dist );



			Real lk_dome_iso_frac = 0;
			Real center_frac = 0;
			lk_dome.single_water_atom_fractions( NZ, water_xyz, phe.xyz("CZ"), CZ_type_index, lk_dome_iso_frac, center_frac );

			Real water_other_distance;
			Vector wadj_water;
			Real dw_dist2 = lk_dome.dome_water_dist2( NZ, water_xyz, phe.xyz("CZ"), water_other_distance, wadj_water );
			Real dw_dist = std::sqrt( dw_dist2 );
			Real _;
			Vector _2;
			Vector dome_water = lk_dome.get_dome_water( NZ, water_xyz, phe.xyz("CZ"), _, _2);

			Real d2_low = lk_dome.d2_low_[ CZ_type_index ];
			Real d2_delta = dw_dist2 - d2_low;

			// std::cout << dw_dist2 << " " << d2_low << " " << d2_delta << " " << lk_dome.ramp_width_A2_ << std::endl;

			Real dome_frac = lk_dome.get_lkd_frac_dome_dist( dw_dist2, CZ_type_index );
			Real iso_frac = lk_dome.my_lk_fraction( water_other_distance );

			Real diso_frac_dr =  lk_dome.dmy_lk_fraction( water_other_distance );

			Real ddome_frac_dr_over_r = 0;
			if ( d2_delta < lk_dome.ramp_width_A2_ && d2_delta > 0 ) {
				// std::cout << "ddome nonzero" << std::endl;
				ddome_frac_dr_over_r = lk_dome.eval_d_lk_fraction_dr_over_r( d2_delta, lk_dome.ramp_width_A2_ );
			}


			Real center_deriv_iso_dome = -avail * sol_value * diso_frac_dr * dome_frac;

			// Real direction = (dist < lk_dome.w_dist() ? -1 : 1 );
			// center_deriv_iso_dome *= direction;


			Real center_deriv_dome_dr = -avail * sol_value * iso_frac * ddome_frac_dr_over_r * dw_dist;
			DerivativeFinder finder( NZ, water_xyz, phe.xyz("CZ"),
				lk_dome.min_angle_cos_, lk_dome.min_angle_sin_,
				lk_dome.max_angle_cos_, lk_dome.max_angle_sin_, lk_dome.w_dist(), lk_dome.water_adjust() );

			Vector f2d = dome_water - phe.xyz("CZ");
			numeric::xyzMatrix<Real> ddome_dother = finder.dDome_dOther();
			Vector ddome_dotherx ( ddome_dother(1,1), ddome_dother(2,1), ddome_dother(3,1) );
			Vector ddome_dothery ( ddome_dother(1,2), ddome_dother(2,2), ddome_dother(3,2) );
			Vector ddome_dotherz ( ddome_dother(1,3), ddome_dother(2,3), ddome_dother(3,3) );
			Vector dRdother_times_r( f2d.dot( ddome_dotherx ), f2d.dot( ddome_dothery ), f2d.dot( ddome_dotherz ) );
			dRdother_times_r -= f2d;
			Vector dr_dother = dRdother_times_r / dw_dist;

			Real dr_dderiv_unit = dr_dother.dot( unit );

			Real disor_dderive_unit = unit.dot( (phe.xyz("CZ") - wadj_water ).normalized() );


			Real center_deriv = center_deriv_iso_dome * disor_dderive_unit + dr_dderiv_unit * center_deriv_dome_dr;




			Real low_score;
			Real center_score;
			Real high_score;
			{
				lk_dome.single_water_atom_fractions( NZ, water_xyz, phe.xyz("CZ"), CZ_type_index, lk_dome_iso_frac, center_score );
				center_score *= -1 * sol_value * avail;
			}
			{
				pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * (dist-deriv_size) );

				lk_dome.single_water_atom_fractions( NZ, water_xyz, phe.xyz("CZ"), CZ_type_index, lk_dome_iso_frac, low_score );
				low_score *= -1 * sol_value * avail;
			}
			{
				pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * (dist+deriv_size) );

				lk_dome.single_water_atom_fractions( NZ, water_xyz, phe.xyz("CZ"), CZ_type_index, lk_dome_iso_frac, high_score );
				high_score *= -1 * sol_value * avail;
			}

			Real low_deriv = (center_score - low_score) / deriv_size;
			Real high_deriv = (high_score - center_score) / deriv_size;

			TS_ASSERT_DELTA( (high_deriv + low_deriv)/2.0, center_deriv, std::abs(center_deriv) / 100 );

			// std::cout << "" << std::endl;
			// std::cout << "Distance: " << dist << std::endl;
			// std::cout << low_deriv << std::endl;
			// std::cout << center_deriv << std::endl;
			// std::cout << high_deriv << std::endl;

		}


	}


	Real
	get_R(
		LK_DomeEnergy const & lk_dome,
		Vector const & base,
		Vector const & water,
		Vector const & other
	) {
		Real _;
		Vector _2;
		return std::sqrt( lk_dome.dome_water_dist2( base, water, other, _, _2 ));
	}

	void test_dr_datom() {



		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "KF", "fa_standard" );

		pose.set_chi(1, 1, 180);
		pose.set_chi(2, 1, 180);
		pose.set_chi(3, 1, 180);
		pose.set_chi(4, 1, 180);


		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
		sfxn->set_weight(core::scoring::lk_dome_iso, 1);

		sfxn->score(pose);

		LK_DomeEnergy lk_dome(sfxn->energy_method_options());

		conformation::Residue const & lys = pose.residue(1);
		conformation::Residue const & phe = pose.residue(2);

		// Creating a copy
		LKD_ResidueInfo lys_info = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( 1 ).get_raw_ptr( LK_DOME_INFO )));
		LKD_ResidueInfo phe_info = *( static_cast< LKD_ResidueInfo * >( pose.residue_data( 2 ).get_raw_ptr( LK_DOME_INFO )));
		lys_info.build_waters( lys, true );

		Vector NZ = lys.xyz("NZ");
		Size NZ_no = lys.atom_index("NZ");
		// Size NZ_type_index = lys.atom_type_index(NZ_no);
		Vector _3HZ = lys.xyz("3HZ");

		Vector unit = (_3HZ - NZ).normalized();
		Vector suspected_water_pos = NZ + unit * lk_dome.w_dist();

		Size CZ_no = phe.atom_index("CZ");
		// Size CZ_type_index = phe.atom_type_index(CZ_no);

		bool found = false;
		//Size water_iatom = 0;
		//(void)water_iatom;
		Size water_iwat = 0;
		Size water_offset = 0;
		for ( Size iatom = 1; iatom <= lys.nheavyatoms(); iatom++ ) {
			Size offset = lys_info.water_offset_for_atom()[iatom];
			for ( Size iwat = 1; iwat <= lys_info.n_attached_waters()[iatom]; iwat++ ) {
				Real dist = suspected_water_pos.distance( lys_info.waters()[iwat + offset] );
				if ( dist < 0.1 ) {
					found = true;
					//water_iatom = iatom;
					water_iwat = iwat;
					water_offset = offset;
				} else {
					lys_info.water_occlusions()[iatom][iwat] = lk_dome.occlusion_max();
				}
			}
		}

		runtime_assert(found);
		Vector water_xyz = lys_info.waters()[water_iwat + water_offset];

		etable::count_pair::CountPairFunctionOP cp = core::scoring::etable::count_pair::CountPairFactory::create_count_pair_function(
			lys, phe, etable::count_pair::CP_CROSSOVER_4 );


		pose.set_xyz(core::id::AtomID(CZ_no, 2), NZ + unit * 3.65 + Vector(1, 0, 0) );

		Vector base = NZ;
		Vector water = water_xyz;
		Vector other = phe.xyz("CZ");


		DerivativeFinder finder( base, water, other, lk_dome.min_angle_cos_, lk_dome.min_angle_sin_,
			lk_dome.max_angle_cos_, lk_dome.max_angle_sin_, lk_dome.w_dist(), lk_dome.water_adjust() );

		numeric::xyzMatrix<Real> ddome_dother = finder.dDome_dOther();
		numeric::xyzMatrix<Real> ddome_dwater = finder.dDome_dWater();
		numeric::xyzMatrix<Real> ddome_dbase = finder.dDome_dBase();

		Real _;
		Vector wadj_water;
		Vector dome_water = lk_dome.get_dome_water( base, water, other, _, wadj_water);
		Vector f2d = dome_water - other;

		Real water_other_distance;
		Vector _2;
		Real dw_dist2 = lk_dome.dome_water_dist2( base, water, other, water_other_distance, _2 );
		Real dw_dist = std::sqrt(dw_dist2);


		core::pose::Pose save_pose = pose;


		{
			Vector ddome_dotherx ( ddome_dother(1,1), ddome_dother(2,1), ddome_dother(3,1) );
			Vector ddome_dothery ( ddome_dother(1,2), ddome_dother(2,2), ddome_dother(3,2) );
			Vector ddome_dotherz ( ddome_dother(1,3), ddome_dother(2,3), ddome_dother(3,3) );
			Vector dRdother_times_r( f2d.dot( ddome_dotherx ), f2d.dot( ddome_dothery ), f2d.dot( ddome_dotherz ) );
			dRdother_times_r -= f2d;

			Vector dRdother = dRdother_times_r / dw_dist;

			Real r = 0;
			Real dr = 0;
			Real expected_dr;
			Vector pert;
			pert = Vector(0.00, 00.00, 0.001);
			r = get_R( lk_dome, base, water, other + pert );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdother );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );

			pert = Vector(0.00, 00.001, 0.00);
			r = get_R( lk_dome, base, water, other + pert );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdother );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );

			pert = Vector(0.001, 00.00, 0.00);
			r = get_R( lk_dome, base, water, other + pert );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdother );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );

			pert = Vector(0.00, 00.00, -0.001);
			r = get_R( lk_dome, base, water, other + pert );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdother );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );

			pert = Vector(0.00, -00.01, 0.00);
			r = get_R( lk_dome, base, water, other + pert );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdother );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );

			pert = Vector(-0.001, 00.00, 0.00);
			r = get_R( lk_dome, base, water, other + pert );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdother );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );



		}



		WaterBuilders const & rsd1_wb( lys_info.get_water_builder( lys , NZ_no ) );

		pose = save_pose;

		{
			Size atom3 = rsd1_wb[water_iwat].atom3();
			// Vector const & r1_atom3_xyz( rsd1.xyz( atom3 ) );

			numeric::xyzMatrix< Real >const & dwater_datom3 = lys_info.atom3_derivs()[water_iwat + water_offset];
			numeric::xyzMatrix< Real > ddome_datom3 = ddome_dwater * dwater_datom3;

			if ( atom3 == NZ_no ) ddome_datom3 += ddome_dbase;

			Vector ddome_datom3x ( ddome_datom3(1,1), ddome_datom3(2,1), ddome_datom3(3,1) );
			Vector ddome_datom3y ( ddome_datom3(1,2), ddome_datom3(2,2), ddome_datom3(3,2) );
			Vector ddome_datom3z ( ddome_datom3(1,3), ddome_datom3(2,3), ddome_datom3(3,3) );
			Vector dRdatom_times_R( f2d.dot( ddome_datom3x ), f2d.dot( ddome_datom3y ), f2d.dot( ddome_datom3z ) );

			Vector dRdatom = dRdatom_times_R / dw_dist;

			Real r = 0;
			Real dr = 0;
			Real expected_dr;
			Vector pert;
			Vector new_water;

			LKD_ResidueInfo info( pose.residue(1), &lk_dome );


			pert = Vector(0.00, 00.00, 0.001);

			pose.set_xyz( core::id::AtomID(atom3, 1), save_pose.residue(1).xyz(atom3) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );

			pert = Vector(0.00, 00.001, 0.00);

			pose.set_xyz( core::id::AtomID(atom3, 1), save_pose.residue(1).xyz(atom3) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );

			pert = Vector(0.001, 00.00, 0.00);

			pose.set_xyz( core::id::AtomID(atom3, 1), save_pose.residue(1).xyz(atom3) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );




		}

		pose = save_pose;

		{
			Size atom2 = rsd1_wb[water_iwat].atom2();
			// Vector const & r1_atom2_xyz( rsd1.xyz( atom2 ) );

			numeric::xyzMatrix< Real >const & dwater_datom2 = lys_info.atom2_derivs()[water_iwat + water_offset];
			numeric::xyzMatrix< Real > ddome_datom2 = ddome_dwater * dwater_datom2;

			if ( atom2 == NZ_no ) ddome_datom2 += ddome_dbase;

			Vector ddome_datom2x ( ddome_datom2(1,1), ddome_datom2(2,1), ddome_datom2(3,1) );
			Vector ddome_datom2y ( ddome_datom2(1,2), ddome_datom2(2,2), ddome_datom2(3,2) );
			Vector ddome_datom2z ( ddome_datom2(1,3), ddome_datom2(2,3), ddome_datom2(3,3) );
			Vector dRdatom_times_R( f2d.dot( ddome_datom2x ), f2d.dot( ddome_datom2y ), f2d.dot( ddome_datom2z ) );

			Vector dRdatom = dRdatom_times_R / dw_dist;

			Real r = 0;
			Real dr = 0;
			Real expected_dr;
			Vector pert;
			Vector new_water;

			LKD_ResidueInfo info( pose.residue(1), &lk_dome );


			pert = Vector(0.00, 00.00, 0.001);

			pose.set_xyz( core::id::AtomID(atom2, 1), save_pose.residue(1).xyz(atom2) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );


			pert = Vector(0.00, 00.001, 0.00);

			pose.set_xyz( core::id::AtomID(atom2, 1), save_pose.residue(1).xyz(atom2) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );


			pert = Vector(0.001, 00.00, 0.00);

			pose.set_xyz( core::id::AtomID(atom2, 1), save_pose.residue(1).xyz(atom2) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );




		}

		pose = save_pose;


		{
			Size atom1 = rsd1_wb[water_iwat].atom1();
			// Vector const & r1_atom1_xyz( lys.xyz( atom1 ) );

			numeric::xyzMatrix< Real >const & dwater_datom1 = lys_info.atom1_derivs()[water_iwat + water_offset];
			numeric::xyzMatrix< Real > ddome_datom1 = ddome_dwater * dwater_datom1;

			if ( atom1 == NZ_no ) ddome_datom1 += ddome_dbase;

			Vector ddome_datom1x ( ddome_datom1(1,1), ddome_datom1(2,1), ddome_datom1(3,1) );
			Vector ddome_datom1y ( ddome_datom1(1,2), ddome_datom1(2,2), ddome_datom1(3,2) );
			Vector ddome_datom1z ( ddome_datom1(1,3), ddome_datom1(2,3), ddome_datom1(3,3) );
			Vector dRdatom_times_R( f2d.dot( ddome_datom1x ), f2d.dot( ddome_datom1y ), f2d.dot( ddome_datom1z ) );

			Vector dRdatom = dRdatom_times_R / dw_dist;

			Real r = 0;
			Real dr = 0;
			Real expected_dr;
			Vector pert;
			Vector new_water;

			LKD_ResidueInfo info( pose.residue(1), &lk_dome );


			pert = Vector(0.00, 00.00, 0.001);

			pose.set_xyz( core::id::AtomID(atom1, 1), save_pose.residue(1).xyz(atom1) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );


			pert = Vector(0.00, 00.001, 0.00);

			pose.set_xyz( core::id::AtomID(atom1, 1), save_pose.residue(1).xyz(atom1) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );


			pert = Vector(0.001, 00.00, 0.00);

			pose.set_xyz( core::id::AtomID(atom1, 1), save_pose.residue(1).xyz(atom1) + pert );
			info.build_waters( pose.residue(1), true );

			new_water = info.waters()[water_iwat + water_offset];
			runtime_assert( new_water.distance(water) < 0.1 );

			r = get_R( lk_dome, pose.residue(1).xyz("NZ"), new_water, other );
			dr = r - dw_dist;
			expected_dr = pert.dot( dRdatom );

			TS_ASSERT_DELTA( dr, expected_dr, pert.norm()/100 );



		}




	}




	void test_lkball_numeric_deriv_check()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		// sfxn.set_weight( fa_atr, 0.5 );
		// sfxn.set_weight( fa_rep, 0.25 );
		// sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( lk_dome, 0.625 );
		sfxn.set_weight( lk_dome_iso, 0.665 );
		sfxn.set_weight( lk_dome_bridge, 0.669 );
		sfxn.set_weight( lk_dome_bridge_uncpl, 0.73 );
		sfxn.set_weight( lk_ball_bridge2, 0.83 );
		sfxn.set_weight( lk_ball_bridge_uncpl2, 0.89 );

		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );

		adv.simple_deriv_check( true, 5e-3 );
	}

	void my_derivative_line_scan(
		LK_DomeEnergy const & lk_dome,
		Vector base,
		Vector water,
		Vector other,
		int scan_who,
		Vector end
	) {

		Vector unit;
		Size steps = 100;
		Real step_size = 0;

		switch (scan_who) {
		case 0 : {
			unit = (end - base).normalized();
			step_size = (end - base).norm() / (steps-1);
			break;
		}
		case 1 : {
			unit = (end - water).normalized();
			step_size = (end - water).norm() / (steps-1);
			break;
		}
		case 2 : {
			unit = (end - other).normalized();
			step_size = (end - other).norm() / (steps-1);
			break;
		}
		}

		Real deriv_size = 0.0001;

		for ( Size i = 0; i < steps; i++ ) {

			DerivativeFinder finder( base, water, other, lk_dome.min_angle_cos_, lk_dome.min_angle_sin_,
				lk_dome.max_angle_cos_, lk_dome.max_angle_sin_, lk_dome.w_dist(), lk_dome.water_adjust() );

			Real colinear_center;
			Vector _2;
			Vector old_dome_water = lk_dome.get_dome_water( base, water, other, colinear_center, _2 );

			Real colinear_pre;
			Real colinear_post;
			numeric::xyzMatrix<Real> ddome_datom;
			Vector dome_low;
			Vector dome_high;
			switch (scan_who) {
			case 0 : {
				dome_low = lk_dome.get_dome_water( base - deriv_size * unit , water, other, colinear_pre, _2 );
				dome_high= lk_dome.get_dome_water( base + deriv_size * unit , water, other, colinear_post, _2 );
				ddome_datom = finder.dDome_dBase();
				base += step_size;
				break;
			}
			case 1 : {
				dome_low = lk_dome.get_dome_water( base, water - deriv_size * unit, other, colinear_pre, _2 );
				dome_high= lk_dome.get_dome_water( base, water + deriv_size * unit, other, colinear_post, _2 );
				ddome_datom = finder.dDome_dWater();
				water += step_size;
				break;
			}
			case 2 : {
				dome_low = lk_dome.get_dome_water( base, water, other - deriv_size * unit, colinear_pre, _2 );
				dome_high= lk_dome.get_dome_water( base, water, other + deriv_size * unit, colinear_post, _2 );
				ddome_datom = finder.dDome_dOther();
				other += step_size;
				break;
			}
			}

			Vector dome_movement = ddome_datom * deriv_size * unit;

			Vector deriv_dome_low = old_dome_water - dome_movement;
			Vector deriv_dome_high = old_dome_water + dome_movement;

			Real distance_low = dome_low.distance( deriv_dome_low );
			Real distance_high = dome_high.distance( deriv_dome_high );

			// std::cout << " " << std::endl;

			// std::cout << "     " << dome_movement << std::endl;
			// std::cout << old_dome_water - dome_low << std::endl;
			// std::cout << old_dome_water << std::endl; //(colinear_center ? " colinear" : " fringe") <<
			// //  (colinear_pre ? " colinear" : " fringe") <<  (colinear_post ? " colinear" : " fringe") << std::endl;
			// std::cout << dome_high - old_dome_water << std::endl;

			TS_ASSERT_DELTA( distance_low, 0, deriv_size * 0.001 );
			TS_ASSERT_DELTA( distance_high, 0, deriv_size * 0.001 );


		}


	}


	void test_derivative_finder() {

		core::scoring::methods::EnergyMethodOptions opts;
		LK_DomeEnergy lk_dome( opts );


		std::cout << "ASDFASDFASDF " << lk_dome.min_angle_cos_ << " " << lk_dome.min_angle_sin_ << std::endl;

		{
			Vector base(0, 0, 0);
			Vector water(lk_dome.w_dist(), 0, 0);
			Vector other(-3, 5, 0);
			Vector end(13, 5, 0);
			my_derivative_line_scan( lk_dome, base, water, other, 2, end );
		}
		{
			Vector base(0, 0, 0);
			Vector water(lk_dome.w_dist(), 0, 0);
			Vector other(5, -3, 2);
			Vector end(5, 2, -2);
			my_derivative_line_scan( lk_dome, base, water, other, 2, end );
		}

		{
			Vector base(0, -5, -5);
			Vector end(0, 5, 5);
			Vector water(lk_dome.w_dist(), 0, 0);
			Vector other(lk_dome.w_dist(), 5, 0);
			my_derivative_line_scan( lk_dome, base, water, other, 0, end );
		}

		{
			Vector base(0, 0, 0);
			Vector water(lk_dome.w_dist(), -5, 5);
			Vector end(lk_dome.w_dist(), 5, -5);
			Vector other(lk_dome.w_dist(), 6, 0);
			my_derivative_line_scan( lk_dome, base, water, other, 1, end );
		}
		{
			Vector base(0, 0, 0);
			Vector water(lk_dome.w_dist(), -5, 5);
			Vector end(lk_dome.w_dist(), 5, -5);
			Vector other(4, 6, 0);
			my_derivative_line_scan( lk_dome, base, water, other, 1, end );
		}




	}




	void test_actual_packing() {

		// core::pose::PoseOP pose_p( new core::pose::Pose() );
		// import_pose::pose_from_file( *pose_p, "/home/bcov/sc/random/trp_cage.pdb" );
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::pose::Pose og_pose = pose;

		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
		// sfxn->score(pose);
		// This is just trpcage sequence + 2 ascending letters
		// All 20 aas are covered
		utility::vector1<std::string> allowed_aas {
			"NAC",
			"LDE",
			"YFG",
			"IHI",
			"QKL",
			"WMN",
			"LPW",
			"KRS",
			"DTV",
			"GWY",
			"GAC",
			"PDE",
			"SFG",
			"SHI",
			"GKL",
			"RMN",
			"PPW",
			"PRS",
			"PTV",
			"SWY"
			};

		core::pack::task::TaskFactoryOP tf = utility::pointer::make_shared<core::pack::task::TaskFactory>();
		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
			core::pack::task::operation::RestrictAbsentCanonicalAASOP restr =
				utility::pointer::make_shared<core::pack::task::operation::RestrictAbsentCanonicalAAS>();
			restr->keep_aas( allowed_aas[seqpos] );
			restr->include_residue( seqpos );
			tf->push_back( restr );
		}
		pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );

		// core::scoring::ScoreFunctionOP sfxn = utility::pointer::make_shared<scoring::ScoreFunction>();

		LK_DomeEnergyOP lk_dome = utility::pointer::make_shared<LK_DomeEnergy>( sfxn->energy_method_options() );


		sfxn->set_weight(core::scoring::lk_ball, 1);
		sfxn->set_weight(core::scoring::lk_ball_iso, 1);
		sfxn->set_weight(core::scoring::lk_dome, 1);
		sfxn->set_weight(core::scoring::lk_dome_iso, 1);
		sfxn->set_weight(core::scoring::lk_dome_bridge, 1);
		sfxn->set_weight(core::scoring::lk_dome_bridge_uncpl, 1);
		sfxn->set_weight(core::scoring::lk_ball_bridge2, 1);
		sfxn->set_weight(core::scoring::lk_ball_bridge_uncpl2, 1);

		pack_rotamers( pose, *sfxn, task );

		sfxn->score( pose );



	}

	void test_interaction_range() {


		std::cout << "Testing interaction" << std::endl;

		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "AAAAAAAAAAAAAAAAAAAA", "fa_standard" );


		core::scoring::ScoreFunctionOP sfxn = utility::pointer::make_shared<scoring::ScoreFunction>();
		sfxn->set_weight(scoring::fa_rep, 1);
		scoring::methods::EnergyMethodOptions opts = sfxn->energy_method_options();
		opts.etable_options().max_dis = 6;
		sfxn->set_energy_method_options( opts );

		LK_DomeEnergyOP lk_dome = utility::pointer::make_shared<LK_DomeEnergy>( sfxn->energy_method_options() );


		sfxn->set_weight(core::scoring::lk_dome, 1);
		sfxn->set_weight(core::scoring::lk_dome_iso, 1);
		sfxn->set_weight(core::scoring::lk_dome_bridge, 1);
		sfxn->set_weight(core::scoring::lk_dome_bridge_uncpl, 1);
		sfxn->set_weight(core::scoring::lk_ball_bridge2, 1);
		sfxn->set_weight(core::scoring::lk_ball_bridge_uncpl2, 1);

		sfxn->score(pose);

		Real expected_range = lk_dome->atomic_interaction_cutoff();

		Energies const & energies( pose.energies() );
		EnergyGraph const & energy_graph( energies.energy_graph() );

		for ( Size seqpos = 2; seqpos <= pose.size(); seqpos++ ) {

			Real distance = pose.residue(seqpos).nbr_atom_xyz().distance( pose.residue(1).nbr_atom_xyz() );
			Real buffer = pose.residue(seqpos).nbr_radius() + pose.residue(1).nbr_radius();

			if ( distance < buffer + expected_range ) {
				TS_ASSERT( energy_graph.find_edge( 1, seqpos ) != nullptr );
			} else {
				TS_ASSERT( energy_graph.find_edge( 1, seqpos ) == nullptr );
			}
		}

		utility::vector1<bool> residues( pose.size() );

		pack::task::TaskFactoryOP tf = utility::pointer::make_shared< pack::task::TaskFactory >();
		tf->push_back( utility::pointer::make_shared<core::pack::task::operation::PreventRepacking>() );
		pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );


		sfxn->setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );
		utility::graph::GraphOP packer_neighbor_graph = core::pack::create_packer_graph( pose, *sfxn, task );


		expected_range = 6;
		for ( Size seqpos = 2; seqpos <= pose.size(); seqpos++ ) {

			Real distance = pose.residue(seqpos).nbr_atom_xyz().distance( pose.residue(1).nbr_atom_xyz() );
			Real buffer = pose.residue(seqpos).nbr_radius() + pose.residue(1).nbr_radius();

			if ( distance < buffer + expected_range ) {
				TS_ASSERT( energy_graph.find_edge( 1, seqpos ) != nullptr );
			} else {
				TS_ASSERT( energy_graph.find_edge( 1, seqpos ) == nullptr );
			}
		}

		std::cout << "Interaction ok" << std::endl;

	}




	void test_packing() {


		// core::pose::PoseOP pose_p( new core::pose::Pose() );
		// import_pose::pose_from_file( *pose_p, "/home/bcov/sc/random/trp_cage.pdb" );
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::pose::Pose og_pose = pose;

		rotamer_set::RotamerSetsOP rotamer_sets_p = trpcage_rotsets( pose );
		rotamer_set::RotamerSets & rotamer_sets = *rotamer_sets_p;

		// Magic numbers that disable countpair
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
		sfxn->set_weight(core::scoring::lk_dome, 1.737050808);
		sfxn->set_weight(core::scoring::lk_dome_iso, 2.241067977);
		sfxn->set_weight(core::scoring::lk_dome_bridge, 2.454489743);
		sfxn->set_weight(core::scoring::lk_dome_bridge_uncpl, 2.650751311);
		sfxn->set_weight(core::scoring::lk_ball_bridge2, 1.1);
		sfxn->set_weight(core::scoring::lk_ball_bridge_uncpl2, 2.2);


		LK_DomeEnergyOP lk_dome = utility::pointer::make_shared<LK_DomeEnergy>( sfxn->energy_method_options() );
		LK_DomeHelper helper = LK_DomeHelper(
			lk_dome,
			sfxn->weights()[core::scoring::lk_dome],
			sfxn->weights()[core::scoring::lk_dome_iso],
			sfxn->weights()[core::scoring::lk_dome_bridge],
			sfxn->weights()[core::scoring::lk_dome_bridge_uncpl],
			sfxn->weights()[core::scoring::lk_ball_bridge2],
			sfxn->weights()[core::scoring::lk_ball_bridge_uncpl2]
		);
		helper.debug_store_waters_ = true;
		helper.init_with_pose( pose, rotamer_sets );

		utility::vector1<utility::vector0<LKD_ResidueInfo>> lkd_infos = get_lkd_infos( pose, rotamer_sets, lk_dome );
		std::map<std::pair<Size, MyWaterHolder>, Size> water_map = get_water_map( helper );

		utility::vector1< core::conformation::ResidueCOP > res_vector;
		utility::vector1< Size > current_rotamers;
		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
			res_vector.push_back( pose.residue(seqpos).get_self_ptr() );
			current_rotamers.push_back(0);
		}

		// positions = pose.size()
		// total waters = pose.size() * helper.waters_per_seqpos
		// error in interaction = partial_score_scaler = score_scaler * occlusion_span // units = score
		// error in avail = 1 / occlusion_span // units = fraction


		// Error in conversion from int16 and int8 * number of possible interactions
		Real score_delta = helper.score_scaler_ * helper.occlusion_span_ * (1 + 1.0/helper.occlusion_span_) * pose.size() * pose.size() * helper.waters_per_seqpos_ * 2 ;

		Real init_score = helper.calculate_energy( res_vector, current_rotamers, 0 );

		TS_ASSERT_DELTA( sfxn->score(pose), init_score, score_delta );

		assert_correct_waters( helper, current_rotamers, og_pose, rotamer_sets, water_map, lkd_infos );
		assert_correct_scores( helper, current_rotamers );
		assert_score_sum( helper, init_score );

		// std::cout << "Pack!" << std::endl;

		// random but repeatable
		utility::vector1<Size> order = { 9, 7, 14, 16, 6, 10, 13, 5, 17, 15, 8, 19, 2, 18, 12, 4, 1, 3, 11,
			13, 10, 15, 6, 18, 3, 8, 17, 1, 9, 4, 5, 16, 12, 19, 11, 14, 7, 2,
			9, 7, 4, 12, 6, 1, 3, 11, 19, 10, 13, 5, 2, 8, 17, 14, 18, 16, 15,
			11, 18, 6, 19, 14, 8, 10, 12, 3, 16, 5, 15, 7, 4, 1, 17, 2, 9, 13,
			5, 15, 9, 11, 14, 8, 17, 13, 18, 4, 3, 7, 1, 12, 2, 19, 6, 16, 10,
			2, 5, 7, 17, 14, 8, 6, 1, 9, 15, 3, 4, 19, 12, 11, 13, 18, 16, 10
			};



		for ( Size repeat = 1; repeat <= 3; repeat++ ) {

			for ( Size letter1 : order ) {
				for ( core::Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

					// repeatable random
					Size num_rots = rotamer_sets.rotamer_set_for_residue( seqpos )->num_rotamers();
					Size choice = (((seqpos << letter1) + 100 - repeat - letter1 ) % num_rots) + 1;

					// std::cout << " " << repeat << " " << letter1 << " " << seqpos << std::endl;
					core::conformation::ResidueCOP rotamer = rotamer_sets.rotamer_set_for_residue( seqpos )->rotamer( choice );
					utility::vector1< core::conformation::ResidueCOP > old_res_vector = res_vector;
					utility::vector1< Size > old_current_rotamers = current_rotamers;

					// Do a repeatable random to get commit status
					double temp = std::sqrt( repeat * letter1 * seqpos );
					double deci = temp - (long)temp;
					bool accept = deci > 0.5;

					pose.replace_residue( seqpos, *rotamer, false );
					res_vector[seqpos] = rotamer;
					current_rotamers[seqpos] = choice;

					Real score = helper.calculate_energy( res_vector, current_rotamers, seqpos );

					// Don't check every time so that we can make this semi-fast
					double temp2 = std::sqrt( repeat * letter1 * seqpos * seqpos );
					double deci2 = temp2 - (long)temp2;
					bool accept2 = deci2 < 0.05;

					// accept2 = true;
					// accept = true;
					// std::cout << "Mutated position " << seqpos << std::endl;

					if ( accept2 ) {
						assert_correct_waters( helper, current_rotamers, og_pose, rotamer_sets, water_map, lkd_infos );
						assert_correct_scores( helper, current_rotamers );
						assert_score_sum( helper, score );

						TS_ASSERT_DELTA( score, sfxn->score(pose), score_delta );
					}

					if ( ! accept ) {
						pose.replace_residue( seqpos, *old_res_vector[seqpos], false );
						res_vector = old_res_vector;
						current_rotamers = old_current_rotamers;
					} else {
						helper.commit_considered_substitution();
					}
				}
			}
		}
	}







	void test_helper_setup() {


		// core::pose::PoseOP pose_p( new core::pose::Pose() );
		// import_pose::pose_from_file( *pose_p, "/home/bcov/sc/random/trp_cage.pdb" );
		core::pose::Pose pose = create_trpcage_ideal_pose();

		rotamer_set::RotamerSetsOP rotamer_sets_p = trpcage_rotsets( pose );
		rotamer_set::RotamerSets & rotamer_sets = *rotamer_sets_p;

		// Use weird numbers to avoid accidentally being correct
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
		sfxn->set_weight(core::scoring::lk_dome, 0.113);
		sfxn->set_weight(core::scoring::lk_dome_iso, 0.134);
		sfxn->set_weight(core::scoring::lk_dome_bridge, 0.167);
		sfxn->set_weight(core::scoring::lk_dome_bridge_uncpl, 0.189);
		sfxn->set_weight(core::scoring::lk_ball_bridge2, 0.195);
		sfxn->set_weight(core::scoring::lk_ball_bridge_uncpl2, 2.01);

		sfxn->score(pose);

		LK_DomeEnergyOP lk_dome = utility::pointer::make_shared<LK_DomeEnergy>( sfxn->energy_method_options() );
		LK_DomeHelper helper = LK_DomeHelper(
			lk_dome,
			sfxn->weights()[core::scoring::lk_dome],
			sfxn->weights()[core::scoring::lk_dome_iso],
			sfxn->weights()[core::scoring::lk_dome_bridge],
			sfxn->weights()[core::scoring::lk_dome_bridge_uncpl],
			sfxn->weights()[core::scoring::lk_ball_bridge2],
			sfxn->weights()[core::scoring::lk_ball_bridge_uncpl2]
		);
		helper.debug_store_waters_ = true;
		helper.init_with_pose( pose, rotamer_sets );

		utility::vector1<utility::vector0<LKD_ResidueInfo>> lkd_infos = get_lkd_infos( pose, rotamer_sets, lk_dome );
		std::map<std::pair<Size, MyWaterHolder>, Size> water_map = get_water_map( helper );


		assert_score_data( pose, *sfxn, helper, rotamer_sets, lkd_infos, water_map );



	}



	void test_simple() {

		// core::pose::PoseOP pose_p( new core::pose::Pose() );
		// import_pose::pose_from_file( *pose_p, "/home/bcov/sc/random/trp_cage.pdb" );
		core::pose::Pose pose = create_trpcage_ideal_pose();

		{
			ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
			sfxn->set_weight(core::scoring::lk_dome, 1);

			sfxn->score( pose );

			Energies const & energies = pose.energies();
			EnergyGraph const & graph = energies.energy_graph();

			for ( Size i = 1; i <= pose.size(); i++ ) {

				for ( Size j = 1; j <= pose.size(); j++ ) {
					if ( i > j ) continue;

					EnergyEdge const * edge = graph.find_energy_edge( i, j );
					if ( ! edge ) continue;

					Real score = edge->dot( sfxn->weights() );

					if ( score != 0 ) {
						std::cout << boost::str(boost::format( "%2i%s -- %2i%s -- %6.3f")%i%pose.residue(i).name1()%j%pose.residue(j).name1()%score) << std::endl;
					}
				}

			}

			EnergyMap emap;
			LK_DomeEnergy lkd = LK_DomeEnergy( sfxn->energy_method_options() );
			lkd.dump_dome_waters_ = true;
			lkd.residue_pair_energy( pose.residue(1), pose.residue(2), pose, *sfxn, emap );

			std::ofstream f("test.pdb");
			for ( Size i = 1; i <= lkd.debug_dome_waters_.size(); i++ ) {
				numeric::xyzVector<Real> xyz = lkd.debug_dome_waters_[i];
				f << boost::str(boost::format("HETATM%5i VOXL VOX A%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s")
					%i%i%xyz.x()%xyz.y()%xyz.z()%1.0%1.0%"HB") << std::endl;
			}
			f.close();
		}

		{
			ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "none" );
			sfxn->set_weight(core::scoring::lk_dome_bridge, 1);
			// sfxn->set_weight(core::scoring::lk_dome_iso, 1);

			sfxn->score( pose );

			Energies const & energies = pose.energies();
			EnergyGraph const & graph = energies.energy_graph();

			for ( Size i = 1; i <= pose.size(); i++ ) {

				for ( Size j = 1; j <= pose.size(); j++ ) {
					if ( i > j ) continue;

					EnergyEdge const * edge = graph.find_energy_edge( i, j );
					if ( ! edge ) continue;

					Real score = edge->dot( sfxn->weights() );

					if ( score != 0 ) {
						std::cout << boost::str(boost::format( "%2i%s -- %2i%s -- %6.3f")%i%pose.residue(i).name1()%j%pose.residue(j).name1()%score) << std::endl;
					}
				}

			}



			EnergyMap emap;
			LK_DomeEnergy lkd = LK_DomeEnergy( sfxn->energy_method_options() );
			lkd.dump_dome_bridge_waters_ = true;

			for ( Size i = 1; i <= pose.size(); i++ ) {

				for ( Size j = 1; j <= pose.size(); j++ ) {
					if ( i >= j ) continue;

					EnergyEdge const * edge = graph.find_energy_edge( i, j );
					if ( ! edge ) continue;


					lkd.residue_pair_energy( pose.residue(i), pose.residue(j), pose, *sfxn, emap );
				}

			}


			std::ofstream f("bridge_waters.pdb");
			for ( Size i = 1; i <= lkd.debug_dome_waters_.size(); i++ ) {
				numeric::xyzVector<Real> xyz = lkd.debug_dome_waters_[i];
				f << boost::str(boost::format("HETATM%5i VOXL VOX A%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s")
					%i%i%xyz.x()%xyz.y()%xyz.z()%1.0%1.0%"HB") << std::endl;
			}
			f.close();
		}

	}



};
