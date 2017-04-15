// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/NPDHBondInteractionGraph.cxxtest.hh
/// @brief  test suite for the hpatch score based interaction graph
/// @author Ron Jacak

// Test framework headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>



// Core Headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>

#include <basic/options/keys/packing.OptionKeys.gen.hh>

//#include <core/pack/annealer/AnnealerFactory.hh>
//#include <core/pack/annealer/SimAnnealerBase.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/interaction_graph/DensePDInteractionGraph.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/NPDHBondInteractionGraph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>


#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondEnergy.hh>
#include <core/scoring/hbonds/NPDHBondSet.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// ObjexxFCL Header
#include <ObjexxFCL/FArray1.io.hh>
#include <ObjexxFCL/FArray1D.hh>

// Utility Headers
#include <utility/graph/Graph.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/LexicographicalIterator.hh>

// Numeric headers

// Test headers
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <assert.h>
#include <basic/options/option.hh>

// C++ headers
#include <list>
#include <map>
#include <tuple>

static basic::Tracer TR("test.core.pack.interactiongraph.NPDHBondInteractionGraph");

using namespace core;
using namespace core::conformation;
using namespace core::pose;
using namespace core::pack;
using namespace core::pack::interaction_graph;
using namespace core::pack::rotamer_set;
using namespace core::scoring;
using namespace core::scoring::hbonds;
using namespace ObjexxFCL::format;


// --------------- Test Class --------------- //

class NPDHBondInteractionGraphTests : public CxxTest::TestSuite {

	void setUp() {
		core_init();
	}

	void tearDown() {}


public:

	// --------------- Test Cases --------------- //


	/// @details
	/// Tests to make sure when doing a design on only some residues that certain residues are indeed being treated and set
	/// as background nodes. If this array returns the wrong indices, background nodes are not being set properly.
	///
	void dont_test_find_networks_with_multiple_hbs_to_donors_and_acceptors() {
		// first we need to find two rotamers that hbond to a third residue on the same atom
		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunctionOP sfxn( new ScoreFunction );
		//sfxn->set_weight( fa_atr, 1.0 );
		sfxn->set_weight( hbond, 1.0 );
		(*sfxn)( *trpcage ); // score the pose to get things started

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *trpcage );
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		utility::graph::GraphOP neighbors = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation( *trpcage );
		rotsets->build_rotamers( *trpcage, *sfxn, neighbors );

		HBondOptions hb_opts;
		//hb_opts.decompose_bb_hb_into_pair_energies( true );
		hb_opts.bb_donor_acceptor_check( false );
		HBondDatabaseCOP hbdb( HBondDatabase::get_database( hb_opts.params_database_tag() ));
		HBondEnergy hbe( hb_opts );
		HBondSet hbset( hb_opts, *trpcage, false );

		Size zero( 0 );
		std::tuple< Size, Size, Size > uninit = std::make_tuple( zero, zero, zero );
		utility::vector1< std::tuple< Size, Size, Size > > empty_vect( trpcage->total_residue(), uninit );

		std::map< std::tuple< Size, Size, Size >, utility::vector1< std::tuple< Size, Size, Size > > > hbmap;
		std::map< std::tuple< Size, Size, Size >, Size > nhbmap;

		EnergyMap emap;
		for ( Size ii = 1; ii <= trpcage->total_residue(); ++ii ) {
			RotamerSetCOP iirots( rotsets->rotamer_set_for_residue( ii ));
			for ( Size jj = ii+1; jj <= trpcage->total_residue(); ++jj ) {
				RotamerSetCOP jjrots( rotsets->rotamer_set_for_residue( jj ));
				for ( Size kk = 1; kk <= iirots->num_rotamers(); ++kk ) {
					ResidueCOP kkrot = iirots->rotamer( kk );
					for ( Size ll = 1; ll <= jjrots->num_rotamers(); ++ll ) {
						ResidueCOP llrot = jjrots->rotamer( ll );
						// look and see if we have any hydrogen bonds between kkrot and llrot
						emap[ hbond ] = 0;
						hbe.residue_pair_energy( *kkrot, *llrot, *trpcage, *sfxn, emap );
						if ( emap[ hbond ] != 0 ) {
							for ( auto kkacc : kkrot->accpt_pos() ) {
								for ( auto lldon : llrot->Hpol_index() ) {
									if ( kkrot->xyz( kkacc ).distance_squared( llrot->xyz( lldon ) ) > MAX_R2 ) continue;

									if ( hb_energy( *hbdb, hb_opts, hbset, *kkrot, kkacc, *llrot, lldon ) ) {
										auto kktup = std::make_tuple( ii, kk, kkacc );
										auto lltup = std::make_tuple( jj, ll, lldon );
										if ( hbmap.find( kktup ) == hbmap.end() ) { hbmap[ kktup ] = empty_vect; nhbmap[ kktup ] = 0; }
										if ( hbmap.find( lltup ) == hbmap.end() ) { hbmap[ lltup ] = empty_vect; nhbmap[ lltup ] = 0; }
										if ( std::get< 0 >( hbmap[ kktup ][ jj ] ) == 0 ) { hbmap[ kktup ][ jj ] = lltup; ++nhbmap[ kktup ]; }
										if ( std::get< 0 >( hbmap[ lltup ][ ii ] ) == 0 ) { hbmap[ lltup ][ ii ] = kktup; ++nhbmap[ lltup ]; }
									}
								}
							}
							for (  auto kkdon : kkrot->Hpol_index() ) {
								for ( auto llacc : llrot->accpt_pos() ) {
									if ( kkrot->xyz( kkdon ).distance_squared( llrot->xyz( llacc ) ) > MAX_R2 ) continue;
									if ( hb_energy( *hbdb, hb_opts, hbset, *llrot, llacc, *kkrot, kkdon ) ) {
										auto kktup = std::make_tuple( ii, kk, kkdon );
										auto lltup = std::make_tuple( jj, ll, llacc );
										if ( hbmap.find( kktup ) == hbmap.end() ) { hbmap[ kktup ] = empty_vect; nhbmap[ kktup ] = 0; }
										if ( hbmap.find( lltup ) == hbmap.end() ) { hbmap[ lltup ] = empty_vect; nhbmap[ lltup ] = 0; }
										if ( std::get< 0 >( hbmap[ kktup ][ jj ] ) == 0 ) { hbmap[ kktup ][ jj ] = lltup; ++nhbmap[ kktup ]; }
										if ( std::get< 0 >( hbmap[ lltup ][ ii ] ) == 0 ) { hbmap[ lltup ][ ii ] = kktup; ++nhbmap[ lltup ]; }
									}
								}
							}
						}
					}
				}
			}
		}

		for ( auto const & tup : hbmap ) {
			Size count_hbs = 0;
			for ( auto const & hbonder : tup.second ) {
				if ( std::get<0>( hbonder ) ) {
					++count_hbs;
				}
				if ( count_hbs >= 2 ) break;
			}
			if ( count_hbs < 2 ) continue;

			if ( tup.second.size() >= 2 ) {
				Size seqpos, rotid, atomno;
				std::tie( seqpos, rotid, atomno ) = tup.first;
				ResidueCOP central_res = rotsets->rotamer_set_for_residue( seqpos )->rotamer( rotid );

				std::cout << " HBonding group '" << central_res->atom_name( atomno ) << "' on residue " << central_res->name() << " rotamer #" << rotid << " at position " << seqpos << " is forming two or more hbonds" << std::endl;
				for ( auto const & other_res_info : tup.second ) {
					Size other_seqpos, other_rotid, other_atomno;
					std::tie( other_seqpos, other_rotid, other_atomno ) = other_res_info;
					if ( other_seqpos == 0 ) continue;
					ResidueCOP other_res = rotsets->rotamer_set_for_residue( other_seqpos )->rotamer( other_rotid );
					std::cout << "   to '" << other_res->atom_name( other_atomno ) << "' on residue " << other_res->name() << " rotamer #" << other_rotid << " at " << other_seqpos << std::endl;
				}


				// OK! Here are two useful networks that I will use as test cases for my unit tests
				// Network #1
				//    HBonding group ' OD1' on residue ASP rotamer #9 at position 2 is forming two or more hbonds
				//      to '1HE2' on residue GLN rotamer #161 at 5
				//      to ' HE2' on residue HIS rotamer #52 at 6
				//      to '1HZ ' on residue LYS rotamer #79 at 19

				// Network #2
				//    HBonding group '3HZ ' on residue LYS rotamer #64 at position 2 is forming two or more hbonds
				//      to ' OH ' on residue TYR rotamer #298 at 5
				//      to ' OE2' on residue GLU rotamer #29 at 6
				//      to ' OE2' on residue GLU rotamer #19 at 19


			}
		}
	}

	void output_coords_for_rotamer(
		Size resid,
		Size rotid,
		ResidueCOP rotamer
	) {
		Size precision( std::cout.precision() );
		std::cout.precision( 16 );
		std::cout << "{\n";
		std::cout << "// Next rotamer for " + rotamer->name() << "\n";
		std::cout << "std::map< std::string, core::Vector > coords;\n";
		for ( Size ii = 1; ii <= rotamer->natoms(); ++ii ) {
			std::cout << "coords[ \"" + rotamer->atom_name( ii ) + "\" ] = core::Vector( ";
			for ( Size jj = 1; jj <= 3; ++jj ) {
				std::cout << rotamer->xyz( ii )( jj ) << ( jj < 3 ? "," : "" );
			}
			std::cout << " );\n";
		}
		std::cout << "rotamers[ std::make_tuple< Size, Size, std::string >( " << resid << ", " << rotid << ", \"" + rotamer->name() + "\" ) ] = coords;\n";
		std::cout << "}\n";
		std::cout.precision( precision );
	}

	void dont_test_output_coordinates_for_networks() {
		// OK! Here are two useful networks that I will use as test cases for my unit tests
		// Network #1
		//    HBonding group ' OD1' on residue ASP rotamer #9 at position 2 is forming two or more hbonds
		//      to '1HE2' on residue GLN rotamer #161 at 5
		//      to ' HE2' on residue HIS rotamer #52 at 6
		//      to '1HZ ' on residue LYS rotamer #79 at 19

		// Network #2
		//    HBonding group '3HZ ' on residue LYS rotamer #64 at position 2 is forming two or more hbonds
		//      to ' OH ' on residue TYR rotamer #298 at 5
		//      to ' OE2' on residue GLU rotamer #29 at 6
		//      to ' OE2' on residue GLU rotamer #19 at 19

		// first we need to find two rotamers that hbond to a third residue on the same atom
		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunctionOP sfxn( new ScoreFunction );
		//sfxn->set_weight( fa_atr, 1.0 );
		sfxn->set_weight( hbond, 1.0 );
		(*sfxn)( *trpcage ); // score the pose to get things started

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *trpcage );
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		utility::graph::GraphOP neighbors = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation( *trpcage );
		rotsets->build_rotamers( *trpcage, *sfxn, neighbors );

		utility::vector1< std::pair< Size, Size > > rots_to_output;
		rots_to_output.push_back( std::make_pair< Size, Size > ( 2, 9 ) );
		rots_to_output.push_back( std::make_pair< Size, Size > ( 5, 161 ));
		rots_to_output.push_back( std::make_pair< Size, Size > ( 6, 52 ));
		rots_to_output.push_back( std::make_pair< Size, Size > ( 19, 79 ));

		rots_to_output.push_back( std::make_pair< Size, Size > ( 2, 65 ));
		rots_to_output.push_back( std::make_pair< Size, Size > ( 5, 298 ));
		rots_to_output.push_back( std::make_pair< Size, Size > ( 6, 29 ));
		rots_to_output.push_back( std::make_pair< Size, Size > ( 19, 19 ));

		std::map< Size, Size > nrots_for_res;
		for ( auto rot : rots_to_output ) {
			if ( ! nrots_for_res.count( rot.first ) ) { nrots_for_res[ rot.first ] = 0; }
			output_coords_for_rotamer( rot.first, ++nrots_for_res[ rot.first ], rotsets->rotamer_set_for_residue( rot.first )->rotamer( rot.second ) );
		}
	}


	void dont_test_find_extra_leu_rots() {
		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunctionOP sfxn( new ScoreFunction );
		//sfxn->set_weight( fa_atr, 1.0 );
		sfxn->set_weight( hbond, 1.0 );
		(*sfxn)( *trpcage ); // score the pose to get things started

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *trpcage );
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		utility::graph::GraphOP neighbors = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation( *trpcage );
		rotsets->build_rotamers( *trpcage, *sfxn, neighbors );

		std::set< Size > residues_to_find_leu_for = {2, 5, 6, 19 };
		for ( Size ii = 1; ii <= trpcage->total_residue(); ++ii ) {
			if ( residues_to_find_leu_for.count( ii ) == 0 ) continue;
			RotamerSetOP iirots = rotsets->rotamer_set_for_residue( ii );

			for ( Size jj = 1; jj <= iirots->num_rotamers(); ++jj ) {
				if ( iirots->rotamer( jj )->name() == "LEU" ) {
					output_coords_for_rotamer( ii, 3, iirots->rotamer( jj ) );
					break;
				}
			}
		}
	}

	std::map< std::tuple< Size, Size, std::string >, std::map< std::string, Vector > >
	get_rotamer_coords_for_trpcage_hb_networks()
	{
		std::map< std::tuple< Size, Size, std::string >, std::map< std::string, Vector > > rotamers;
		{
			// Next rotamer for ASP
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -6.414,4.03,-2.127 );
			coords[ " CA " ] = core::Vector( -4.993,3.959,-2.449 );
			coords[ " C  " ] = core::Vector( -4.238,3.117,-1.429 );
			coords[ " O  " ] = core::Vector( -3.492,2.207,-1.791 );
			coords[ " CB " ] = core::Vector( -4.385625144324645,5.362672106342716,-2.509794306252662 );
			coords[ " CG " ] = core::Vector( -2.928460180918675,5.359156886684338,-2.952008494081177 );
			coords[ " OD1" ] = core::Vector( -2.647382857910334,4.827880892253254,-4.000150740744638 );
			coords[ " OD2" ] = core::Vector( -2.110655139366914,5.887700680738954,-2.237497878048868 );
			coords[ " H  " ] = core::Vector( -6.865,4.932,-2.088 );
			coords[ " HA " ] = core::Vector( -4.862,3.469,-3.413 );
			coords[ "1HB " ] = core::Vector( -4.959215522798488,5.97766779272022,-3.203619413973113 );
			coords[ "2HB " ] = core::Vector( -4.45079287268357,5.830345026896177,-1.526947553476557 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 2, 1, "ASP" ) ] = coords;
		}
		{
			// Next rotamer for GLN
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -4.861,-0.541,-1.607 );
			coords[ " CA " ] = core::Vector( -4.136,-1.183,-2.697 );
			coords[ " C  " ] = core::Vector( -2.671,-1.395,-2.335 );
			coords[ " O  " ] = core::Vector( -2.035,-2.336,-2.808 );
			coords[ " CB " ] = core::Vector( -4.242985106969453,-0.3486099722784353,-3.976261726497569 );
			coords[ " CG " ] = core::Vector( -3.572172191482379,1.011644468070044,-3.890127189452449 );
			coords[ " CD " ] = core::Vector( -4.485694352687442,2.071403788082159,-3.304198542743064 );
			coords[ " OE1" ] = core::Vector( -5.351135184871185,1.775728749192312,-2.475495208453527 );
			coords[ " NE2" ] = core::Vector( -4.298496355589819,3.3147968505922,-3.731754524885219 );
			coords[ " H  " ] = core::Vector( -5.203,0.401,-1.731 );
			coords[ " HA " ] = core::Vector( -4.555,-2.171,-2.883 );
			coords[ "1HB " ] = core::Vector( -3.792177992359282,-0.8941628315319847,-4.805264153018222 );
			coords[ "2HB " ] = core::Vector( -5.293275875798802,-0.1908564002725024,-4.220958277857048 );
			coords[ "1HG " ] = core::Vector( -2.690437892732531,0.9308640340519853,-3.254423525873286 );
			coords[ "2HG " ] = core::Vector( -3.282243955028009,1.326447558177028,-4.892594461818761 );
			coords[ "1HE2" ] = core::Vector( -4.871339727584427,4.056738079285267,-3.380629760483347 );
			coords[ "2HE2" ] = core::Vector( -3.585130418373254,3.510620998629973,-4.404755230472079 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 5, 1, "GLN" ) ] = coords;
		}
		{
			// Next rotamer for HIS
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -2.142,-0.514,-1.493 );
			coords[ " CA " ] = core::Vector( -0.771,-0.642,-1.012 );
			coords[ " C  " ] = core::Vector( -0.643,-1.78,-0.008999999999999999 );
			coords[ " O  " ] = core::Vector( 0.338,-2.524,-0.021 );
			coords[ " CB " ] = core::Vector( -0.2975628480929781,0.6663713475232996,-0.3705944203312281 );
			coords[ " CG " ] = core::Vector( -0.086834047703648,1.776210576277158,-1.35311683834596 );
			coords[ " ND1" ] = core::Vector( 1.159861695392433,2.115799361855943,-1.835325587831322 );
			coords[ " CD2" ] = core::Vector( -0.962520527900795,2.623482672770471,-1.942748673602035 );
			coords[ " CE1" ] = core::Vector( 1.040898495012494,3.124877666096042,-2.680958154903031 );
			coords[ " NE2" ] = core::Vector( -0.2359038103592044,3.45115121757672,-2.762930744746846 );
			coords[ " H  " ] = core::Vector( -2.703,0.264,-1.177 );
			coords[ " HA " ] = core::Vector( -0.109,-0.883,-1.844 );
			coords[ "1HB " ] = core::Vector( -1.030869667011066,0.9970035887663173,0.3658213364114248 );
			coords[ "2HB " ] = core::Vector( 0.6406353518005781,0.4945242342074473,0.1562900682670793 );
			coords[ " HD2" ] = core::Vector( -2.042455007264204,2.645365415675913,-1.794103933563863 );
			coords[ " HE1" ] = core::Vector( 1.857730980734911,3.605325181288884,-3.219090314914457 );
			coords[ " HE2" ] = core::Vector( -0.6216588181237066,4.189004327737486,-3.334797432476333 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 6, 1, "HIS" ) ] = coords;
		}
		{
			// Next rotamer for LYS
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( 1.708,5.373,0.359 );
			coords[ " CA " ] = core::Vector( 1.211,6.549,-0.346 );
			coords[ " C  " ] = core::Vector( 0.757,7.626,0.631 );
			coords[ " O  " ] = core::Vector( -0.091,7.381,1.489 );
			coords[ " CB " ] = core::Vector( 0.06120836281959496,6.167746080554171,-1.279819769570765 );
			coords[ " CG " ] = core::Vector( 0.4590753963959231,5.242853276752158,-2.422448338948159 );
			coords[ " CD " ] = core::Vector( -0.7048387577230075,5.006906780877911,-3.373307334819931 );
			coords[ " CE " ] = core::Vector( -1.766522343284927,4.119420815022755,-2.740529010852899 );
			coords[ " NZ " ] = core::Vector( -2.873891215291822,3.815138892527663,-3.686896863549124 );
			coords[ " H  " ] = core::Vector( 1.129017179709997,4.546802000766096,0.4067260921501928 );
			coords[ " HA " ] = core::Vector( 1.983,7.033,-0.963 );
			coords[ "1HB " ] = core::Vector( -0.7240613699099732,5.673439636289669,-0.7072535332998706 );
			coords[ "2HB " ] = core::Vector( -0.3692217288497304,7.070218149169488,-1.71436022154702 );
			coords[ "1HG " ] = core::Vector( 1.286927398609709,5.685453673120805,-2.977688916788921 );
			coords[ "2HG " ] = core::Vector( 0.7867503486275691,4.285594112296414,-2.018249042573859 );
			coords[ "1HD " ] = core::Vector( -1.155358026272526,5.96341918149519,-3.641798661113676 );
			coords[ "2HD " ] = core::Vector( -0.3406290460160286,4.528915882278047,-4.282643987659799 );
			coords[ "1HE " ] = core::Vector( -1.312914819790993,3.183326450526847,-2.418252080973174 );
			coords[ "2HE " ] = core::Vector( -2.181496824049027,4.616986141521774,-1.864072114705411 );
			coords[ "1HZ " ] = core::Vector( -3.555719672891684,3.226189857737639,-3.229964290934857 );
			coords[ "2HZ " ] = core::Vector( -3.314736111497143,4.676430700001008,-3.978390813808839 );
			coords[ "3HZ " ] = core::Vector( -2.504238047090861,3.336234697392992,-4.495548479486056 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 19, 1, "LYS" ) ] = coords;
		}
		{
			// Next rotamer for LYS
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -6.414,4.03,-2.127 );
			coords[ " CA " ] = core::Vector( -4.993,3.959,-2.449 );
			coords[ " C  " ] = core::Vector( -4.238,3.117,-1.429 );
			coords[ " O  " ] = core::Vector( -3.492,2.207,-1.791 );
			coords[ " CB " ] = core::Vector( -4.389229435592382,5.362407115948574,-2.52183058783101 );
			coords[ " CG " ] = core::Vector( -2.904442460537195,5.391704694021104,-2.859283122440771 );
			coords[ " CD " ] = core::Vector( -2.384529745465096,6.819073899751222,-2.941791645295829 );
			coords[ " CE " ] = core::Vector( -0.8724987751809419,6.849537631624214,-3.109261310208411 );
			coords[ " NZ " ] = core::Vector( -0.4437645852300135,6.253600512822128,-4.403676763015196 );
			coords[ " H  " ] = core::Vector( -6.865,4.932,-2.088 );
			coords[ " HA " ] = core::Vector( -4.862,3.469,-3.413 );
			coords[ "1HB " ] = core::Vector( -4.915596016305317,5.946313833115372,-3.277350724727484 );
			coords[ "2HB " ] = core::Vector( -4.524819182597342,5.867013749197883,-1.564992634446567 );
			coords[ "1HG " ] = core::Vector( -2.346008868022482,4.854322282153491,-2.091886084472511 );
			coords[ "2HG " ] = core::Vector( -2.739179789302174,4.897964520217069,-3.816363027000588 );
			coords[ "1HD " ] = core::Vector( -2.846155840932461,7.326294692833922,-3.79008380676622 );
			coords[ "2HD " ] = core::Vector( -2.650248616009092,7.356277934512105,-2.031394771225268 );
			coords[ "1HE " ] = core::Vector( -0.5224306212944905,7.879729861368302,-3.0639720007255 );
			coords[ "2HE " ] = core::Vector( -0.4051225582932677,6.294317803131941,-2.296093321069158 );
			coords[ "1HZ " ] = core::Vector( 0.5631262026913825,6.292498562706917,-4.47576192101039 );
			coords[ "2HZ " ] = core::Vector( -0.7473419908564751,5.290820954539404,-4.448763392020134 );
			coords[ "3HZ " ] = core::Vector( -0.8555549402038504,6.772009102786612,-5.166298141051216 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 2, 2, "LYS" ) ] = coords;
		}
		{
			// Next rotamer for TYR
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -4.861,-0.541,-1.607 );
			coords[ " CA " ] = core::Vector( -4.136,-1.183,-2.697 );
			coords[ " C  " ] = core::Vector( -2.671,-1.395,-2.335 );
			coords[ " O  " ] = core::Vector( -2.035,-2.336,-2.808 );
			coords[ " CB " ] = core::Vector( -4.251473106590077,-0.353704840575134,-3.977979969948759 );
			coords[ " CG " ] = core::Vector( -3.453728654308245,0.9311276127973455,-3.947360859079676 );
			coords[ " CD1" ] = core::Vector( -2.13267048979636,0.9375601447814781,-4.370532120549265 );
			coords[ " CD2" ] = core::Vector( -4.043076675491849,2.10270138849191,-3.4962105792086 );
			coords[ " CE1" ] = core::Vector( -1.404067659054359,2.111034381899273,-4.342566810617988 );
			coords[ " CE2" ] = core::Vector( -3.314575101701091,3.276012543037695,-3.4682491557322 );
			coords[ " CZ " ] = core::Vector( -2.000380331442816,3.282411655597731,-3.889221883267187 );
			coords[ " OH " ] = core::Vector( -1.274724499885924,4.451139510359296,-3.861369685366363 );
			coords[ " H  " ] = core::Vector( -5.203,0.401,-1.731 );
			coords[ " HA " ] = core::Vector( -4.555,-2.171,-2.883 );
			coords[ "1HB " ] = core::Vector( -3.909658259878626,-0.9465307477935512,-4.8276517542439 );
			coords[ "2HB " ] = core::Vector( -5.296496814390581,-0.1009745179851631,-4.154671727997792 );
			coords[ " HD1" ] = core::Vector( -1.669382771464004,0.01658353691999026,-4.725182303927683 );
			coords[ " HD2" ] = core::Vector( -5.080817481808779,2.097430906861058,-3.162725914839952 );
			coords[ " HE1" ] = core::Vector( -0.3663590171244058,2.116087220758864,-4.674973427404346 );
			coords[ " HE2" ] = core::Vector( -3.777632186202212,4.196857445778932,-3.112796819427033 );
			coords[ " HH " ] = core::Vector( -1.826710985303299,5.161365595761237,-3.52529117571094 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 5, 2, "TYR" ) ] = coords;
		}
		{
			// Next rotamer for GLU
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -2.142,-0.514,-1.493 );
			coords[ " CA " ] = core::Vector( -0.771,-0.642,-1.012 );
			coords[ " C  " ] = core::Vector( -0.643,-1.78,-0.008999999999999999 );
			coords[ " O  " ] = core::Vector( 0.338,-2.524,-0.021 );
			coords[ " CB " ] = core::Vector( -0.3049778570357484,0.6680370281497833,-0.3728579274492325 );
			coords[ " CG " ] = core::Vector( -0.1223839970005628,1.816601245662128,-1.354831697053954 );
			coords[ " CD " ] = core::Vector( 0.283299324012068,3.099944640231869,-0.685103355985445 );
			coords[ " OE1" ] = core::Vector( 1.233769436345559,3.084670095048803,0.05962984493567058 );
			coords[ " OE2" ] = core::Vector( -0.3582882416237854,4.097149568595722,-0.9186017333509082 );
			coords[ " H  " ] = core::Vector( -2.703,0.264,-1.177 );
			coords[ " HA " ] = core::Vector( -0.109,-0.883,-1.844 );
			coords[ "1HB " ] = core::Vector( -1.0279759972999,0.9831405583451993,0.3797387479014938 );
			coords[ "2HB " ] = core::Vector( 0.6463879018175743,0.5071999585720603,0.1343554729367327 );
			coords[ "1HG " ] = core::Vector( 0.6421072330318428,1.541877329840707,-2.081129558914558 );
			coords[ "2HG " ] = core::Vector( -1.056330912238951,1.972508610433322,-1.892965446564986 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 6, 2, "GLU" ) ] = coords;
		}
		{
			// Next rotamer for GLU
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( 1.708,5.373,0.359 );
			coords[ " CA " ] = core::Vector( 1.211,6.549,-0.346 );
			coords[ " C  " ] = core::Vector( 0.757,7.626,0.631 );
			coords[ " O  " ] = core::Vector( -0.091,7.381,1.489 );
			coords[ " CB " ] = core::Vector( 0.05429619796640539,6.166648265742737,-1.272131582137489 );
			coords[ " CG " ] = core::Vector( 0.4484285435755422,5.258576512797375,-2.42836605200568 );
			coords[ " CD " ] = core::Vector( -0.7067813651876076,4.926152020704196,-3.331177937126844 );
			coords[ " OE1" ] = core::Vector( -1.811449865256254,5.294235159849608,-3.01105828081542 );
			coords[ " OE2" ] = core::Vector( -0.484224601290146,4.303140125977175,-4.342562350510866 );
			coords[ " H  " ] = core::Vector( 1.129017509555413,4.546802963544627,0.4067255980679335 );
			coords[ " HA " ] = core::Vector( 1.983,7.033,-0.963 );
			coords[ "1HB " ] = core::Vector( -0.7199021992237609,5.658705108969328,-0.6967804783832312 );
			coords[ "2HB " ] = core::Vector( -0.3892869772495784,7.069475534333179,-1.692093120084509 );
			coords[ "1HG " ] = core::Vector( 1.223377161638183,5.750767395468097,-3.015419226124192 );
			coords[ "2HG " ] = core::Vector( 0.8662664595441574,4.336668854833608,-2.026309307997114 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 19, 2, "GLU" ) ] = coords;
		}
		{
			// Next rotamer for LEU
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -6.414,4.03,-2.127 );
			coords[ " CA " ] = core::Vector( -4.993,3.959,-2.449 );
			coords[ " C  " ] = core::Vector( -4.238,3.117,-1.429 );
			coords[ " O  " ] = core::Vector( -3.492,2.207,-1.791 );
			coords[ " CB " ] = core::Vector( -4.390531327182051,5.368537470546533,-2.503227695942783 );
			coords[ " CG " ] = core::Vector( -2.885663432848604,5.438105568053777,-2.792692971066834 );
			coords[ " CD1" ] = core::Vector( -2.605719101870337,4.854929528256522,-4.171124321709845 );
			coords[ " CD2" ] = core::Vector( -2.417894788329695,6.883184515310088,-2.705008443916748 );
			coords[ " H  " ] = core::Vector( -6.865,4.932,-2.088 );
			coords[ " HA " ] = core::Vector( -4.862,3.469,-3.413 );
			coords[ "1HB " ] = core::Vector( -4.903340146051702,5.935817797194991,-3.278279685351964 );
			coords[ "2HB " ] = core::Vector( -4.566939808519612,5.859062860919908,-1.546051570125166 );
			coords[ " HG " ] = core::Vector( -2.346012827018542,4.836890249997746,-2.060477801514174 );
			coords[ "1HD1" ] = core::Vector( -1.53632223378118,4.904192875725371,-4.376743608057858 );
			coords[ "2HD1" ] = core::Vector( -2.93226633426701,3.815300483359948,-4.199556699477272 );
			coords[ "3HD1" ] = core::Vector( -3.146424110256977,5.426604226899061,-4.924486989007987 );
			coords[ "1HD2" ] = core::Vector( -1.348001382541205,6.932539745802514,-2.90999420760376 );
			coords[ "2HD2" ] = core::Vector( -2.955919595293047,7.484954362834577,-3.437790078169681 );
			coords[ "3HD2" ] = core::Vector( -2.613824892722565,7.269005067829715,-1.704471545223325 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 2, 3, "LEU" ) ] = coords;
		}
		{
			// Next rotamer for LEU
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -4.861,-0.541,-1.607 );
			coords[ " CA " ] = core::Vector( -4.136,-1.183,-2.697 );
			coords[ " C  " ] = core::Vector( -2.671,-1.395,-2.335 );
			coords[ " O  " ] = core::Vector( -2.035,-2.336,-2.808 );
			coords[ " CB " ] = core::Vector( -4.236906645350885,-0.3358621853193817,-3.971707965563555 );
			coords[ " CG " ] = core::Vector( -5.631669296083212,-0.2497155990267943,-4.604537451982502 );
			coords[ " CD1" ] = core::Vector( -5.603257811969437,0.7355320154769536,-5.765147651817539 );
			coords[ " CD2" ] = core::Vector( -6.063972618782808,-1.631846154095095,-5.070990247910963 );
			coords[ " H  " ] = core::Vector( -5.203,0.401,-1.731 );
			coords[ " HA " ] = core::Vector( -4.555,-2.171,-2.883 );
			coords[ "1HB " ] = core::Vector( -3.914750808010703,0.6783000663497081,-3.741108781356134 );
			coords[ "2HB " ] = core::Vector( -3.559121133682176,-0.7492972309666466,-4.718433825519288 );
			coords[ " HG " ] = core::Vector( -6.343262808781459,0.1236611632363818,-3.867620797503421 );
			coords[ "1HD1" ] = core::Vector( -6.594435193959708,0.7968902161054188,-6.214723450151425 );
			coords[ "2HD1" ] = core::Vector( -5.309944810480971,1.719792556795579,-5.399813946768315 );
			coords[ "3HD1" ] = core::Vector( -4.886636745795849,0.3966447579133716,-6.512368402210797 );
			coords[ "1HD2" ] = core::Vector( -7.055776338019824,-1.570746269131183,-5.52012184155142 );
			coords[ "2HD2" ] = core::Vector( -5.353813199751154,-2.005509209859701,-5.808971335821273 );
			coords[ "3HD2" ] = core::Vector( -6.092775557980173,-2.311001047272155,-4.218793549977763 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 5, 3, "LEU" ) ] = coords;
		}
		{
			// Next rotamer for LEU
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( -2.142,-0.514,-1.493 );
			coords[ " CA " ] = core::Vector( -0.771,-0.642,-1.012 );
			coords[ " C  " ] = core::Vector( -0.643,-1.78,-0.008999999999999999 );
			coords[ " O  " ] = core::Vector( 0.338,-2.524,-0.021 );
			coords[ " CB " ] = core::Vector( -0.3108881928227309,0.6707366846198047,-0.3656838642286033 );
			coords[ " CG " ] = core::Vector( -0.1413639798117423,1.858964059440474,-1.321002916702062 );
			coords[ " CD1" ] = core::Vector( 0.1565907467030179,3.118065056145086,-0.5182428209567841 );
			coords[ " CD2" ] = core::Vector( 0.9780665811568019,1.561745597879457,-2.307555322901021 );
			coords[ " H  " ] = core::Vector( -2.703,0.264,-1.177 );
			coords[ " HA " ] = core::Vector( -0.109,-0.883,-1.844 );
			coords[ "1HB " ] = core::Vector( -1.036908576489176,0.9572238161455209,0.3934636615995429 );
			coords[ "2HB " ] = core::Vector( 0.6476928231643392,0.4998359107186782,0.1240472946669605 );
			coords[ " HG " ] = core::Vector( -1.07129338432383,2.023339763605501,-1.866010824961243 );
			coords[ "1HD1" ] = core::Vector( 0.2768697665670552,3.962467548272209,-1.197089798169653 );
			coords[ "2HD1" ] = core::Vector( -0.6688730918713807,3.31775251502024,0.1651343926361204 );
			coords[ "3HD1" ] = core::Vector( 1.07437946260508,2.977393139128162,0.05149787476946177 );
			coords[ "1HD2" ] = core::Vector( 1.098042232123132,2.406070540978128,-2.987149492591749 );
			coords[ "2HD2" ] = core::Vector( 1.908717243758543,1.398619480633151,-1.763638252138312 );
			coords[ "3HD2" ] = core::Vector( 0.7307861906957946,0.6675574996142046,-2.879919056566269 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 6, 3, "LEU" ) ] = coords;
		}
		{
			// Next rotamer for LEU
			std::map< std::string, core::Vector > coords;
			coords[ " N  " ] = core::Vector( 1.708,5.373,0.359 );
			coords[ " CA " ] = core::Vector( 1.211,6.549,-0.346 );
			coords[ " C  " ] = core::Vector( 0.757,7.626,0.631 );
			coords[ " O  " ] = core::Vector( -0.091,7.381,1.489 );
			coords[ " CB " ] = core::Vector( 0.04715894397188913,6.162323935195833,-1.267227790950715 );
			coords[ " CG " ] = core::Vector( -0.5507894421323898,7.304342100498562,-2.098798282708859 );
			coords[ " CD1" ] = core::Vector( 0.4917961687576367,7.815868983511477,-3.083629030703742 );
			coords[ " CD2" ] = core::Vector( -1.791948417712211,6.807742837035008,-2.825202311996351 );
			coords[ " H  " ] = core::Vector( 1.129018658464556,4.546803402843235,0.4067259782655011 );
			coords[ " HA " ] = core::Vector( 1.983,7.033,-0.963 );
			coords[ "1HB " ] = core::Vector( 0.3915381998497804,5.394678324960459,-1.958320104037286 );
			coords[ "2HB " ] = core::Vector( -0.7529866222183514,5.741426876344327,-0.6585127628133689 );
			coords[ " HG " ] = core::Vector( -0.8215508810298116,8.130553523760048,-1.440858607024661 );
			coords[ "1HD1" ] = core::Vector( 0.06694864472282358,8.627524293506605,-3.674421776833213 );
			coords[ "2HD1" ] = core::Vector( 1.360240537876136,8.182670173456657,-2.536334493105228 );
			coords[ "3HD1" ] = core::Vector( 0.7949020732843373,7.005687576835241,-3.745793719099434 );
			coords[ "1HD2" ] = core::Vector( -2.217454567061064,7.61976688516121,-3.415700799693013 );
			coords[ "2HD2" ] = core::Vector( -1.522009812814972,5.982463987735024,-3.484455325853874 );
			coords[ "3HD2" ] = core::Vector( -2.527048044986266,6.464869962977874,-2.096925633516018 );
			rotamers[ std::make_tuple< Size, Size, std::string >( 19, 3, "LEU" ) ] = coords;
		}
		return rotamers;
	}

	void append_network_rotamers_to_rotsets( Pose const & pose, RotamerSetsOP rotsets ) {
		std::map< std::tuple< Size, Size, std::string >, std::map< std::string, Vector > > rots =
			get_rotamer_coords_for_trpcage_hb_networks();
		chemical::ResidueTypeSetCOP fa_standard( core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ));

		for ( auto const & rot_tuple : rots ) {
			Size ii_seqpos = std::get<0>( rot_tuple.first );
			if ( rotsets->resid_2_moltenres( ii_seqpos ) == 0 ) continue;
			RotamerSetOP rset = rotsets->rotamer_set_for_residue( ii_seqpos );
			ResidueOP rot( new Residue(
				fa_standard->name_map( std::get<2>( rot_tuple.first )),
				pose.residue( ii_seqpos ),
				pose.conformation(),
				false ));
			for ( auto const & atom_coord : rot_tuple.second ) {
				rot->set_xyz( atom_coord.first, atom_coord.second );
			}
			rset->add_rotamer( *rot );
		}
	}


	void test_create_overcommitted_hbond_networks_standard_npdhb_ig() {
		Size const leu_rot( 4 );
		Size const acc_net_rot( 2 );
		Size const don_net_rot( 3 );
		Size const input_sc_rot( 1 );

		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunctionOP sfxn( new ScoreFunction );
		HBondOptions hb_opts;
		hb_opts.decompose_bb_hb_into_pair_energies( true );
		hb_opts.bb_donor_acceptor_check( false );
		hb_opts.use_hb_env_dep( false );
		methods::EnergyMethodOptions enmeth_opts = sfxn->energy_method_options();
		enmeth_opts.hbond_options( hb_opts );
		sfxn->set_energy_method_options( enmeth_opts );

		sfxn->set_weight( fa_atr, 1.0 );
		sfxn->set_weight( npd_hbond, 1.25 );
		(*sfxn)( *trpcage ); // score the pose to get things started

		utility::vector1< Size > res_in_network_v1( 4 );
		res_in_network_v1[ 1 ] = 2; res_in_network_v1[ 2 ] = 5;
		res_in_network_v1[ 3 ] = 6; res_in_network_v1[ 4 ] = 19;
		std::set< Size > res_in_network;
		for ( auto res : res_in_network_v1 ) res_in_network.insert( res );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *trpcage );
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			task->nonconst_residue_task( ii ).or_include_current( true );
		}
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		utility::graph::GraphOP neighbors = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation( *trpcage );
		rotsets->build_rotamers( *trpcage, *sfxn, neighbors );

		HBondDatabaseCOP hbdb( HBondDatabase::get_database( hb_opts.params_database_tag() ));
		HBondEnergy hbe( hb_opts );
		HBondSet hbset( hb_opts, *trpcage, false );

		// Drop all of the existing rotamers
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			RotamerSetOP iirotset = rotsets->rotamer_set_for_residue( ii );
			Size ii_nrots = iirotset->num_rotamers();
			utility::vector1< bool > rots_to_delete( ii_nrots, true );
			rots_to_delete[ iirotset->id_for_current_rotamer() ] = false;
			iirotset->drop_rotamers( rots_to_delete );
		}

		// Now append the rotamers that are part of the network
		append_network_rotamers_to_rotsets( *trpcage, rotsets );

		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );

		Size count_rotamer( 0 );
		TS_ASSERT_EQUALS( rotsets->nrotamers(), 32 );
		for ( Size ii = 1; ii <= trpcage->total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( rotsets->rotamer_set_for_residue( ii )->id_for_current_rotamer(), 1 );
			if ( res_in_network.count( ii ) ) {
				TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres( ii ), 4 );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					++count_rotamer;
					TS_ASSERT_EQUALS( rotsets->rotid_on_moltenresidue( count_rotamer ), jj );
					TS_ASSERT_EQUALS( rotsets->moltenres_for_rotamer( count_rotamer ), ii );
				}
			} else {
				TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres( ii ), 1 );
				++count_rotamer;
				TS_ASSERT_EQUALS( rotsets->rotid_on_moltenresidue( count_rotamer ), 1 );
				TS_ASSERT_EQUALS( rotsets->moltenres_for_rotamer( count_rotamer ), ii );
			}
		}

		StandardNPDHBondInteractionGraphOP ig( new StandardNPDHBondInteractionGraph( task->num_to_be_packed() ));
		//DensePDInteractionGraphOP ig( new DensePDInteractionGraph( task->num_to_be_packed() ));
		ig->set_pose( *trpcage );
		ig->set_score_function( *sfxn );
		ig->set_packer_neighbor_graph( *neighbors );
		ig->set_packer_task( *task );
		ig->set_rotamer_sets( *rotsets );

		rotsets->compute_energies( *trpcage, *sfxn, neighbors, ig );

		ig->prepare_for_simulated_annealing();

		ObjexxFCL::FArray1D_int state1( trpcage->total_residue(), input_sc_rot );
		for ( Size netres : res_in_network ) {
			state1( netres ) = leu_rot;
			trpcage->replace_residue( netres, *rotsets->rotamer_set_for_residue( netres )->rotamer( leu_rot ), false );
		}
		ig->set_network_state( state1 );
		(*sfxn)(*trpcage);

		// OK -- let's test some rotamer substitutions
		// Size count_output_pdbs( 0 );
		utility::vector1< Size > twos( 4, 2 );
		utility::LexicographicalIterator lex( twos );
		while ( ! lex.at_end() ) {
			for ( Size ii = 1; ii <= 4; ++ii ) {
				Size iires = res_in_network_v1[ ii ];

				// "leave one out" -- look at all combinations of leu/network hbonds
				if ( lex[ ii ] == 2 ) continue;

				Pose copy1( *trpcage );
				ObjexxFCL::FArray1D_int start_state( trpcage->total_residue(), input_sc_rot );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					if ( jj == ii ) {
						start_state( iires ) = leu_rot;
						continue;
					}
					Size jjres = res_in_network_v1[ jj ];
					if ( lex[ jj ] == 2 ) {
						copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( acc_net_rot ), false );
						start_state( jjres ) = acc_net_rot;
					} else {
						start_state( jjres ) = leu_rot;
					}
				}
				// std::cout << "network state: ";
				// for ( Size ii = 1; ii <= start_state.size(); ++ii ) std::cout << " " << start_state( ii );
				// std::cout << std::endl;

				// copy1.dump_pdb( "trpcage_npdhb_before_" + utility::to_string( ++count_output_pdbs ) + ".pdb" );
				ig->set_network_state( start_state );
				// std::cout << "Start score: " << std::endl;
				Real start_score( (*sfxn)( copy1 ) );
				// EnergyMap start_tot = copy1.energies().total_energies();
				// for ( auto iter = copy1.energies().energy_graph().get_node( iires )->edge_list_begin();
				//    iter !=  copy1.energies().energy_graph().get_node( iires )->edge_list_end(); ++iter ) {
				//  EnergyEdge const * ee = static_cast< EnergyEdge const * > (*iter);
				//  EnergyMap emap = ee->fill_energy_map();
				//  std::cout << "Edge " << ee->get_first_node_ind() << " " << ee->get_second_node_ind() << ": ";
				//  emap.show_weighted( std::cout, sfxn->weights() );
				//  std::cout << std::endl;
				// }
				//Pose copy2( copy1 );
				copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( acc_net_rot ), false );
				// copy1.dump_pdb( "trpcage_npdhb_after_" + utility::to_string( count_output_pdbs ) + ".pdb" );
				// std::cout << "End score: " << std::endl;
				Real end_score( (*sfxn )( copy1 ));
				// for ( auto iter = copy1.energies().energy_graph().get_node( iires )->edge_list_begin();
				//    iter !=  copy1.energies().energy_graph().get_node( iires )->edge_list_end(); ++iter ) {
				//  EnergyEdge const * ee = static_cast< EnergyEdge const * > (*iter);
				//  EnergyMap emap = ee->fill_energy_map();
				//  std::cout << "Edge " << ee->get_first_node_ind() << " " << ee->get_second_node_ind() << ": ";
				//  emap.show_weighted( std::cout, sfxn->weights() );
				//  std::cout << std::endl;
				// }
				// EnergyMap end_tot = copy1.energies().total_energies();
				// std::cout << "DeltaE: ";
				// EnergyMap diff( end_tot ); diff -= start_tot;
				// diff.show_nonzero( std::cout );
				// std::cout << std::endl;

				float prev_e, deltaE;
				//ig->print();
				ig->consider_substitution( iires, acc_net_rot, deltaE, prev_e );
				//ig->print();
				TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );

			}
			++lex;
		}

		// And for the network with an over-committed donor
		lex.begin();
		while ( ! lex.at_end() ) {
			for ( Size ii = 1; ii <= 4; ++ii ) {
				Size iires = res_in_network_v1[ ii ];

				// "leave one out" -- look at all combinations of leu/network hbonds
				if ( lex[ ii ] == 2 ) continue;

				Pose copy1( *trpcage );
				ObjexxFCL::FArray1D_int start_state( trpcage->total_residue(), input_sc_rot );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					if ( jj == ii ) {
						start_state( iires ) = leu_rot;
						continue;
					}
					Size jjres = res_in_network_v1[ jj ];
					if ( lex[ jj ] == 2 ) {
						copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( don_net_rot ), false );
						start_state( jjres ) = don_net_rot;
					} else {
						start_state( jjres ) = leu_rot;
					}
				}
				ig->set_network_state( start_state );
				Real start_score( (*sfxn)( copy1 ) );
				copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( don_net_rot ), false );
				Real end_score( (*sfxn )( copy1 ));
				float prev_e, deltaE;
				ig->consider_substitution( iires, don_net_rot, deltaE, prev_e );
				TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );

			}
			++lex;
		}


	}

	void test_create_overcommitted_hbond_networks_dense_npdhb_ig() {
		Size const leu_rot( 4 );
		Size const acc_net_rot( 2 );
		Size const don_net_rot( 3 );
		Size const input_sc_rot( 1 );

		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunctionOP sfxn( new ScoreFunction );
		HBondOptions hb_opts;
		hb_opts.decompose_bb_hb_into_pair_energies( true );
		hb_opts.bb_donor_acceptor_check( false );
		hb_opts.use_hb_env_dep( false );
		methods::EnergyMethodOptions enmeth_opts = sfxn->energy_method_options();
		enmeth_opts.hbond_options( hb_opts );
		sfxn->set_energy_method_options( enmeth_opts );

		sfxn->set_weight( fa_atr, 1.0 );
		sfxn->set_weight( npd_hbond, 1.25 );
		(*sfxn)( *trpcage ); // score the pose to get things started

		utility::vector1< Size > res_in_network_v1( 4 );
		res_in_network_v1[ 1 ] = 2; res_in_network_v1[ 2 ] = 5;
		res_in_network_v1[ 3 ] = 6; res_in_network_v1[ 4 ] = 19;
		std::set< Size > res_in_network;
		for ( auto res : res_in_network_v1 ) res_in_network.insert( res );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *trpcage );
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			task->nonconst_residue_task( ii ).or_include_current( true );
		}
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		utility::graph::GraphOP neighbors = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation( *trpcage );
		rotsets->build_rotamers( *trpcage, *sfxn, neighbors );

		HBondDatabaseCOP hbdb( HBondDatabase::get_database( hb_opts.params_database_tag() ));
		HBondEnergy hbe( hb_opts );
		HBondSet hbset( hb_opts, *trpcage, false );

		// Drop all of the existing rotamers
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			RotamerSetOP iirotset = rotsets->rotamer_set_for_residue( ii );
			Size ii_nrots = iirotset->num_rotamers();
			utility::vector1< bool > rots_to_delete( ii_nrots, true );
			rots_to_delete[ iirotset->id_for_current_rotamer() ] = false;
			iirotset->drop_rotamers( rots_to_delete );
		}

		// Now append the rotamers that are part of the network
		append_network_rotamers_to_rotsets( *trpcage, rotsets );

		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );

		Size count_rotamer( 0 );
		TS_ASSERT_EQUALS( rotsets->nrotamers(), 32 );
		for ( Size ii = 1; ii <= trpcage->total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( rotsets->rotamer_set_for_residue( ii )->id_for_current_rotamer(), 1 );
			if ( res_in_network.count( ii ) ) {
				TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres( ii ), 4 );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					++count_rotamer;
					TS_ASSERT_EQUALS( rotsets->rotid_on_moltenresidue( count_rotamer ), jj );
					TS_ASSERT_EQUALS( rotsets->moltenres_for_rotamer( count_rotamer ), ii );
				}
			} else {
				TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres( ii ), 1 );
				++count_rotamer;
				TS_ASSERT_EQUALS( rotsets->rotid_on_moltenresidue( count_rotamer ), 1 );
				TS_ASSERT_EQUALS( rotsets->moltenres_for_rotamer( count_rotamer ), ii );
			}
		}

		DenseNPDHBondInteractionGraphOP ig( new DenseNPDHBondInteractionGraph( task->num_to_be_packed() ));
		ig->set_pose( *trpcage );
		ig->set_score_function( *sfxn );
		ig->set_packer_neighbor_graph( *neighbors );
		ig->set_packer_task( *task );
		ig->set_rotamer_sets( *rotsets );

		rotsets->compute_energies( *trpcage, *sfxn, neighbors, ig );

		ig->prepare_for_simulated_annealing();

		ObjexxFCL::FArray1D_int state1( trpcage->total_residue(), input_sc_rot );
		for ( Size netres : res_in_network ) {
			state1( netres ) = leu_rot;
			trpcage->replace_residue( netres, *rotsets->rotamer_set_for_residue( netres )->rotamer( leu_rot ), false );
		}
		ig->set_network_state( state1 );
		(*sfxn)(*trpcage);

		// OK -- let's test some rotamer substitutions
		// Size count_output_pdbs( 0 );
		utility::vector1< Size > twos( 4, 2 );
		utility::LexicographicalIterator lex( twos );
		while ( ! lex.at_end() ) {
			for ( Size ii = 1; ii <= 4; ++ii ) {
				Size iires = res_in_network_v1[ ii ];

				// "leave one out" -- look at all combinations of leu/network hbonds
				if ( lex[ ii ] == 2 ) continue;

				Pose copy1( *trpcage );
				ObjexxFCL::FArray1D_int start_state( trpcage->total_residue(), input_sc_rot );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					if ( jj == ii ) {
						start_state( iires ) = leu_rot;
						continue;
					}
					Size jjres = res_in_network_v1[ jj ];
					if ( lex[ jj ] == 2 ) {
						copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( acc_net_rot ), false );
						start_state( jjres ) = acc_net_rot;
					} else {
						start_state( jjres ) = leu_rot;
					}
				}
				ig->set_network_state( start_state );
				Real start_score( (*sfxn)( copy1 ) );
				copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( acc_net_rot ), false );
				Real end_score( (*sfxn )( copy1 ));
				float prev_e, deltaE;
				ig->consider_substitution( iires, acc_net_rot, deltaE, prev_e );
				TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );
			}
			++lex;
		}

		// And for the network with an over-committed donor
		lex.begin();
		while ( ! lex.at_end() ) {
			for ( Size ii = 1; ii <= 4; ++ii ) {
				Size iires = res_in_network_v1[ ii ];

				// "leave one out" -- look at all combinations of leu/network hbonds
				if ( lex[ ii ] == 2 ) continue;

				Pose copy1( *trpcage );
				ObjexxFCL::FArray1D_int start_state( trpcage->total_residue(), input_sc_rot );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					if ( jj == ii ) {
						start_state( iires ) = leu_rot;
						continue;
					}
					Size jjres = res_in_network_v1[ jj ];
					if ( lex[ jj ] == 2 ) {
						copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( don_net_rot ), false );
						start_state( jjres ) = don_net_rot;
					} else {
						start_state( jjres ) = leu_rot;
					}
				}
				ig->set_network_state( start_state );
				Real start_score( (*sfxn)( copy1 ) );
				copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( don_net_rot ), false );
				Real end_score( (*sfxn )( copy1 ));
				float prev_e, deltaE;
				ig->consider_substitution( iires, don_net_rot, deltaE, prev_e );
				TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );

			}
			++lex;
		}


	}

	void test_create_overcommitted_hbond_networks_linmem_npdhb_ig() {
		Size const leu_rot( 4 );
		Size const acc_net_rot( 2 );
		Size const don_net_rot( 3 );
		Size const input_sc_rot( 1 );

		PoseOP trpcage = create_trpcage_ideal_poseop();
		ScoreFunctionOP sfxn( new ScoreFunction );
		HBondOptions hb_opts;
		hb_opts.decompose_bb_hb_into_pair_energies( true );
		hb_opts.bb_donor_acceptor_check( false );
		hb_opts.use_hb_env_dep( false );
		methods::EnergyMethodOptions enmeth_opts = sfxn->energy_method_options();
		enmeth_opts.hbond_options( hb_opts );
		sfxn->set_energy_method_options( enmeth_opts );

		sfxn->set_weight( fa_atr, 1.0 );
		sfxn->set_weight( npd_hbond, 1.25 );
		(*sfxn)( *trpcage ); // score the pose to get things started

		utility::vector1< Size > res_in_network_v1( 4 );
		res_in_network_v1[ 1 ] = 2; res_in_network_v1[ 2 ] = 5;
		res_in_network_v1[ 3 ] = 6; res_in_network_v1[ 4 ] = 19;
		std::set< Size > res_in_network;
		for ( auto res : res_in_network_v1 ) res_in_network.insert( res );

		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *trpcage );
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			task->nonconst_residue_task( ii ).or_include_current( true );
		}
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		utility::graph::GraphOP neighbors = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation( *trpcage );
		rotsets->build_rotamers( *trpcage, *sfxn, neighbors );

		HBondDatabaseCOP hbdb( HBondDatabase::get_database( hb_opts.params_database_tag() ));
		HBondEnergy hbe( hb_opts );
		HBondSet hbset( hb_opts, *trpcage, false );

		// Drop all of the existing rotamers
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			RotamerSetOP iirotset = rotsets->rotamer_set_for_residue( ii );
			Size ii_nrots = iirotset->num_rotamers();
			utility::vector1< bool > rots_to_delete( ii_nrots, true );
			rots_to_delete[ iirotset->id_for_current_rotamer() ] = false;
			iirotset->drop_rotamers( rots_to_delete );
		}

		// Now append the rotamers that are part of the network
		append_network_rotamers_to_rotsets( *trpcage, rotsets );

		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );

		Size count_rotamer( 0 );
		TS_ASSERT_EQUALS( rotsets->nrotamers(), 32 );
		for ( Size ii = 1; ii <= trpcage->total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( rotsets->rotamer_set_for_residue( ii )->id_for_current_rotamer(), 1 );
			if ( res_in_network.count( ii ) ) {
				TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres( ii ), 4 );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					++count_rotamer;
					TS_ASSERT_EQUALS( rotsets->rotid_on_moltenresidue( count_rotamer ), jj );
					TS_ASSERT_EQUALS( rotsets->moltenres_for_rotamer( count_rotamer ), ii );
				}
			} else {
				TS_ASSERT_EQUALS( rotsets->nrotamers_for_moltenres( ii ), 1 );
				++count_rotamer;
				TS_ASSERT_EQUALS( rotsets->rotid_on_moltenresidue( count_rotamer ), 1 );
				TS_ASSERT_EQUALS( rotsets->moltenres_for_rotamer( count_rotamer ), ii );
			}
		}

		LinearMemoryNPDHBondInteractionGraphOP ig( new LinearMemoryNPDHBondInteractionGraph( task->num_to_be_packed() ));
		ig->set_pose( *trpcage );
		ig->set_score_function( *sfxn );
		ig->set_packer_neighbor_graph( *neighbors );
		ig->set_packer_task( *task );
		ig->set_rotamer_sets( *rotsets );

		rotsets->compute_energies( *trpcage, *sfxn, neighbors, ig );

		ig->prepare_for_simulated_annealing();

		ObjexxFCL::FArray1D_int state1( trpcage->total_residue(), input_sc_rot );
		for ( Size netres : res_in_network ) {
			state1( netres ) = leu_rot;
			trpcage->replace_residue( netres, *rotsets->rotamer_set_for_residue( netres )->rotamer( leu_rot ), false );
		}
		ig->set_network_state( state1 );
		(*sfxn)(*trpcage);

		// OK -- let's test some rotamer substitutions
		// Size count_output_pdbs( 0 );
		utility::vector1< Size > twos( 4, 2 );
		utility::LexicographicalIterator lex( twos );
		while ( ! lex.at_end() ) {
			for ( Size ii = 1; ii <= 4; ++ii ) {
				Size iires = res_in_network_v1[ ii ];

				// "leave one out" -- look at all combinations of leu/network hbonds
				if ( lex[ ii ] == 2 ) continue;

				Pose copy1( *trpcage );
				ObjexxFCL::FArray1D_int start_state( trpcage->total_residue(), input_sc_rot );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					if ( jj == ii ) {
						start_state( iires ) = leu_rot;
						continue;
					}
					Size jjres = res_in_network_v1[ jj ];
					if ( lex[ jj ] == 2 ) {
						copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( acc_net_rot ), false );
						start_state( jjres ) = acc_net_rot;
					} else {
						start_state( jjres ) = leu_rot;
					}
				}
				ig->set_network_state( start_state );
				Real start_score( (*sfxn)( copy1 ) );
				copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( acc_net_rot ), false );
				Real end_score( (*sfxn )( copy1 ));
				float prev_e, deltaE;
				ig->consider_substitution( iires, acc_net_rot, deltaE, prev_e );
				TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );
			}
			++lex;
		}

		// And for the network with an over-committed donor
		lex.begin();
		while ( ! lex.at_end() ) {
			for ( Size ii = 1; ii <= 4; ++ii ) {
				Size iires = res_in_network_v1[ ii ];

				// "leave one out" -- look at all combinations of leu/network hbonds
				if ( lex[ ii ] == 2 ) continue;

				Pose copy1( *trpcage );
				ObjexxFCL::FArray1D_int start_state( trpcage->total_residue(), input_sc_rot );
				for ( Size jj = 1; jj <= 4; ++jj ) {
					if ( jj == ii ) {
						start_state( iires ) = leu_rot;
						continue;
					}
					Size jjres = res_in_network_v1[ jj ];
					if ( lex[ jj ] == 2 ) {
						copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( don_net_rot ), false );
						start_state( jjres ) = don_net_rot;
					} else {
						start_state( jjres ) = leu_rot;
					}
				}
				ig->set_network_state( start_state );
				Real start_score( (*sfxn)( copy1 ) );
				copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( don_net_rot ), false );
				Real end_score( (*sfxn )( copy1 ));
				float prev_e, deltaE;
				ig->consider_substitution( iires, don_net_rot, deltaE, prev_e );
				TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );

			}
			++lex;
		}


	}

	void test_create_overcommitted_hbond_networks_w_bgnode_linmem_npdhb_ig_acc_network() {
		Size const leu_rot( 4 );
		Size const acc_net_rot( 2 );
		//Size const don_net_rot( 3 );
		//Size const input_sc_rot( 1 );

		// What we'll do is create a Pose that has one of the acc_net_rots
		// subbed in, and then we'll repack the other acc_net_rots around it.
		PoseOP orig_trpcage = create_trpcage_ideal_poseop();
		ScoreFunctionOP sfxn( new ScoreFunction );
		HBondOptions hb_opts;
		hb_opts.decompose_bb_hb_into_pair_energies( true );
		hb_opts.bb_donor_acceptor_check( false );
		hb_opts.use_hb_env_dep( false );
		methods::EnergyMethodOptions enmeth_opts = sfxn->energy_method_options();
		enmeth_opts.hbond_options( hb_opts );
		sfxn->set_energy_method_options( enmeth_opts );

		sfxn->set_weight( fa_atr, 1.0 );
		sfxn->set_weight( npd_hbond, 1.25 );
		(*sfxn)( *orig_trpcage ); // score the pose to get things started

		utility::vector1< Size > res_in_network_v1( 4 );
		utility::vector1< Size > res_in_subnetwork_v1; res_in_subnetwork_v1.reserve( 3 );
		res_in_network_v1[ 1 ] = 2; res_in_network_v1[ 2 ] = 5;
		res_in_network_v1[ 3 ] = 6; res_in_network_v1[ 4 ] = 19;
		std::set< Size > res_in_network, res_in_subnetwork;
		for ( auto res : res_in_network_v1 ) {
			res_in_network.insert( res );
			if ( res != 6 ) {
				res_in_subnetwork_v1.push_back( res );
				res_in_subnetwork.insert( res );
			}
		}


		task::PackerTaskOP orig_task = task::TaskFactory::create_packer_task( *orig_trpcage );
		for ( Size ii = 1; ii <= orig_trpcage->size(); ++ii ) {
			if ( ii == 6 ) {
				orig_task->nonconst_residue_task( ii ).or_include_current( true );
			} else {
				orig_task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}
		sfxn->setup_for_packing( *orig_trpcage, orig_task->repacking_residues(), orig_task->designing_residues() );
		utility::graph::GraphOP neighbors = create_packer_graph( *orig_trpcage, *sfxn, orig_task );

		RotamerSetsOP orig_rotsets( new RotamerSets() );
		orig_rotsets->set_task( orig_task );
		orig_rotsets->initialize_pose_for_rotsets_creation( *orig_trpcage );
		orig_rotsets->build_rotamers( *orig_trpcage, *sfxn, neighbors );

		HBondDatabaseCOP hbdb( HBondDatabase::get_database( hb_opts.params_database_tag() ));
		HBondEnergy hbe( hb_opts );
		HBondSet hbset( hb_opts, *orig_trpcage, false );

		// Drop all of the existing rotamers for residue 6
		RotamerSetOP res6_rotset = orig_rotsets->rotamer_set_for_residue( 6 );
		Size res6_nrots = res6_rotset->num_rotamers();
		utility::vector1< bool > res6_rots_to_delete( res6_nrots, true );
		res6_rots_to_delete[ res6_rotset->id_for_current_rotamer() ] = false;
		res6_rotset->drop_rotamers( res6_rots_to_delete );

		// Now append the rotamers that are part of the network
		append_network_rotamers_to_rotsets( *orig_trpcage, orig_rotsets );

		// orig_rotsets->prepare_sets_for_packing( *orig_trpcage, *sfxn );

		// Now let's create a Pose and put an acceptor residue in the middle of it
		PoseOP trpcage( new Pose( *orig_trpcage ));
		trpcage->replace_residue( 6, *orig_rotsets->rotamer_set_for_residue( 6 )->rotamer( acc_net_rot ), false );
		(*sfxn)(*trpcage);

		// Now, then. Let's create a packing job where residue 6 is held fixed along with the other,
		// non-network residues.
		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *trpcage );
		for ( Size ii = 1; ii <= orig_trpcage->size(); ++ii ) {
			if ( ii == 6 || res_in_network.count( ii ) == 0 ) {
				task->nonconst_residue_task( ii ).prevent_repacking();
			} else {
				task->nonconst_residue_task( ii ).or_include_current( true );
			}
		}
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		neighbors = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation( *trpcage );
		rotsets->build_rotamers( *trpcage, *sfxn, neighbors );

		// Drop all of the existing rotamers
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			if ( rotsets->resid_2_moltenres( ii ) == 0 ) continue;
			RotamerSetOP iirotset = rotsets->rotamer_set_for_residue( ii );
			Size ii_nrots = iirotset->num_rotamers();
			utility::vector1< bool > rots_to_delete( ii_nrots, true );
			rots_to_delete[ iirotset->id_for_current_rotamer() ] = false;
			iirotset->drop_rotamers( rots_to_delete );
		}

		// Now append the rotamers that are part of the network
		append_network_rotamers_to_rotsets( *trpcage, rotsets );

		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), 3 );
		for ( Size netres : res_in_subnetwork ) {
			TS_ASSERT_EQUALS( rotsets->rotamer_set_for_residue( netres )->num_rotamers(), 4 );
		}

		LinearMemoryNPDHBondInteractionGraphOP ig( new LinearMemoryNPDHBondInteractionGraph( task->num_to_be_packed() ));
		ig->set_pose( *trpcage );
		ig->set_score_function( *sfxn );
		ig->set_packer_neighbor_graph( *neighbors );
		ig->set_packer_task( *task );
		ig->set_rotamer_sets( *rotsets );

		rotsets->compute_energies( *trpcage, *sfxn, neighbors, ig );

		ig->prepare_for_simulated_annealing();

		ObjexxFCL::FArray1D_int state1( res_in_subnetwork.size(), leu_rot );
		for ( Size netres : res_in_subnetwork ) {
			trpcage->replace_residue( netres, *rotsets->rotamer_set_for_residue( netres )->rotamer( leu_rot ), false );
		}
		ig->set_network_state( state1 );
		(*sfxn)(*trpcage);

		// OK -- let's test some rotamer substitutions
		//Size count_output_pdbs( 0 );
		utility::vector1< Size > twos( 3, 2 );
		utility::LexicographicalIterator lex( twos );
		while ( ! lex.at_end() ) {
			for ( Size ii = 1; ii <= 3; ++ii ) {
				Size iires = res_in_subnetwork_v1[ ii ];

				// "leave one out" -- look at all combinations of leu/network hbonds
				if ( lex[ ii ] == 2 ) continue;

				Pose copy1( *trpcage );
				ObjexxFCL::FArray1D_int start_state( 3, leu_rot );
				for ( Size jj = 1; jj <= 3; ++jj ) {
					if ( ii == jj ) continue;
					Size jjres = res_in_subnetwork_v1[ jj ];
					if ( lex[ jj ] == 2 ) {
						copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( acc_net_rot ), false );
						start_state( jj ) = acc_net_rot;
					}
				}
				//++count_output_pdbs;
				//std::cout << "network state #" << count_output_pdbs << ": ";
				//for ( Size ii = 1; ii <= start_state.size(); ++ii ) std::cout << " " << start_state( ii );
				//std::cout << std::endl;
				//
				//copy1.dump_pdb( "trpcage_npdhb_before_" + utility::to_string( count_output_pdbs ) + ".pdb" );
				ig->set_network_state( start_state );
				//std::cout << "Start" << std::endl;
				Real start_score( (*sfxn)( copy1 ) );
				copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( acc_net_rot ), false );
				//copy1.dump_pdb( "trpcage_npdhb_after_" + utility::to_string( count_output_pdbs ) + ".pdb" );
				//std::cout << "End" << std::endl;
				Real end_score( (*sfxn )( copy1 ));
				float prev_e, deltaE;
				ig->consider_substitution( ii, acc_net_rot, deltaE, prev_e );
				//if ( std::abs( (end_score - start_score) - deltaE ) < 1e-5 ) {
				// std::cout << "Passed!" << std::endl;
				//}
				TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );

			}
			++lex;
		}

		//  // OK -- now we have to start over to do the same thing with the donor network;
		//  // when reside 6 gets made a donor residue, the one-body energies for the other
		//  // residues in the network will have changed.
		//  trpcage->replace_residue( 6, *orig_rotsets->rotamer_set_for_residue( 6 )->rotamer( don_net_rot ), false );
		//  (*sfxn)(*trpcage);
		//  sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		//  neighbors = create_packer_graph( *trpcage, *sfxn, task );
		//  ig = LinearMemoryNPDHBondInteractionGraphOP( new LinearMemoryNPDHBondInteractionGraph( task->num_to_be_packed() ));
		//  ig->set_pose( *trpcage );
		//  ig->set_packer_neighbor_graph( *neighbors );
		//  ig->set_score_function( *sfxn );
		//  ig->set_packer_task( *task );
		//  ig->set_rotamer_sets( *rotsets );
		//  ig->set_score_weight( 1.0 );
		//
		//  rotsets->compute_energies( *trpcage, *sfxn, neighbors, ig );
		//
		//  ig->prepare_for_simulated_annealing();
		//
		//  lex.begin();
		//  while ( ! lex.at_end() ) {
		//   for ( Size ii = 1; ii <= 3; ++ii ) {
		//    Size iires = res_in_subnetwork_v1[ ii ];
		//
		//    // "leave one out" -- look at all combinations of leu/network hbonds
		//    if ( lex[ ii ] == 2 ) continue;
		//
		//    Pose copy1( *trpcage );
		//    ObjexxFCL::FArray1D_int start_state( 3, leu_rot );
		//    for ( Size jj = 1; jj <= 3; ++jj ) {
		//     if ( ii == jj ) continue;
		//     Size jjres = res_in_subnetwork_v1[ jj ];
		//     if ( lex[ jj ] == 2 ) {
		//      copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( don_net_rot ), false );
		//      start_state( jj ) = don_net_rot;
		//     }
		//    }
		//    ig->set_network_state( start_state );
		//    Real start_score( (*sfxn)( copy1 ) );
		//    copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( don_net_rot ), false );
		//    Real end_score( (*sfxn )( copy1 ));
		//    float prev_e, deltaE;
		//    ig->consider_substitution( ii, don_net_rot, deltaE, prev_e );
		//    TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );
		//
		//   }
		//   ++lex;
		//  }


	}

	void test_create_overcommitted_hbond_networks_w_bgnode_linmem_npdhb_ig_don_network() {
		Size const leu_rot( 4 );
		//Size const acc_net_rot( 2 );
		Size const don_net_rot( 3 );
		//Size const input_sc_rot( 1 );

		// What we'll do is create a Pose that has one of the don_net_rots
		// subbed in, and then we'll repack the other don_net_rots around it.
		PoseOP orig_trpcage = create_trpcage_ideal_poseop();
		ScoreFunctionOP sfxn( new ScoreFunction );
		HBondOptions hb_opts;
		hb_opts.decompose_bb_hb_into_pair_energies( true );
		hb_opts.bb_donor_acceptor_check( false );
		hb_opts.use_hb_env_dep( false );
		methods::EnergyMethodOptions enmeth_opts = sfxn->energy_method_options();
		enmeth_opts.hbond_options( hb_opts );
		sfxn->set_energy_method_options( enmeth_opts );

		sfxn->set_weight( fa_atr, 1.0 );
		sfxn->set_weight( npd_hbond, 1.25 );
		(*sfxn)( *orig_trpcage ); // score the pose to get things started

		utility::vector1< Size > res_in_network_v1( 4 );
		utility::vector1< Size > res_in_subnetwork_v1; res_in_subnetwork_v1.reserve( 3 );
		res_in_network_v1[ 1 ] = 2; res_in_network_v1[ 2 ] = 5;
		res_in_network_v1[ 3 ] = 6; res_in_network_v1[ 4 ] = 19;
		std::set< Size > res_in_network, res_in_subnetwork;
		for ( auto res : res_in_network_v1 ) {
			res_in_network.insert( res );
			if ( res != 6 ) {
				res_in_subnetwork_v1.push_back( res );
				res_in_subnetwork.insert( res );
			}
		}


		task::PackerTaskOP orig_task = task::TaskFactory::create_packer_task( *orig_trpcage );
		for ( Size ii = 1; ii <= orig_trpcage->size(); ++ii ) {
			if ( ii == 6 ) {
				orig_task->nonconst_residue_task( ii ).or_include_current( true );
			} else {
				orig_task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}
		sfxn->setup_for_packing( *orig_trpcage, orig_task->repacking_residues(), orig_task->designing_residues() );
		utility::graph::GraphOP neighbors = create_packer_graph( *orig_trpcage, *sfxn, orig_task );

		RotamerSetsOP orig_rotsets( new RotamerSets() );
		orig_rotsets->set_task( orig_task );
		orig_rotsets->initialize_pose_for_rotsets_creation( *orig_trpcage );
		orig_rotsets->build_rotamers( *orig_trpcage, *sfxn, neighbors );

		HBondDatabaseCOP hbdb( HBondDatabase::get_database( hb_opts.params_database_tag() ));
		HBondEnergy hbe( hb_opts );
		HBondSet hbset( hb_opts, *orig_trpcage, false );

		// Drop all of the existing rotamers for residue 6
		RotamerSetOP res6_rotset = orig_rotsets->rotamer_set_for_residue( 6 );
		Size res6_nrots = res6_rotset->num_rotamers();
		utility::vector1< bool > res6_rots_to_delete( res6_nrots, true );
		res6_rots_to_delete[ res6_rotset->id_for_current_rotamer() ] = false;
		res6_rotset->drop_rotamers( res6_rots_to_delete );

		// Now append the rotamers that are part of the network
		append_network_rotamers_to_rotsets( *orig_trpcage, orig_rotsets );

		orig_rotsets->prepare_sets_for_packing( *orig_trpcage, *sfxn );

		// Now let's create a Pose and put an acceptor residue in the middle of it
		PoseOP trpcage( new Pose( *orig_trpcage ));
		trpcage->replace_residue( 6, *orig_rotsets->rotamer_set_for_residue( 6 )->rotamer( don_net_rot ), false );
		(*sfxn)(*trpcage);

		// Now, then. Let's create a packing job where residue 6 is held fixed along with the other,
		// non-network residues.
		task::PackerTaskOP task = task::TaskFactory::create_packer_task( *trpcage );
		for ( Size ii = 1; ii <= orig_trpcage->size(); ++ii ) {
			if ( ii == 6 || res_in_network.count( ii ) == 0 ) {
				task->nonconst_residue_task( ii ).prevent_repacking();
			} else {
				task->nonconst_residue_task( ii ).or_include_current( true );
			}
		}
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		neighbors = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation( *trpcage );
		rotsets->build_rotamers( *trpcage, *sfxn, neighbors );

		// Drop all of the existing rotamers
		for ( Size ii = 1; ii <= trpcage->size(); ++ii ) {
			if ( rotsets->resid_2_moltenres( ii ) == 0 ) continue;
			RotamerSetOP iirotset = rotsets->rotamer_set_for_residue( ii );
			Size ii_nrots = iirotset->num_rotamers();
			utility::vector1< bool > rots_to_delete( ii_nrots, true );
			rots_to_delete[ iirotset->id_for_current_rotamer() ] = false;
			iirotset->drop_rotamers( rots_to_delete );
		}

		// Now append the rotamers that are part of the network
		append_network_rotamers_to_rotsets( *trpcage, rotsets );

		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );

		TS_ASSERT_EQUALS( rotsets->nmoltenres(), 3 );
		for ( Size netres : res_in_subnetwork ) {
			TS_ASSERT_EQUALS( rotsets->rotamer_set_for_residue( netres )->num_rotamers(), 4 );
		}

		LinearMemoryNPDHBondInteractionGraphOP ig( new LinearMemoryNPDHBondInteractionGraph( task->num_to_be_packed() ));
		ig->set_pose( *trpcage );
		ig->set_score_function( *sfxn );
		ig->set_packer_neighbor_graph( *neighbors );
		ig->set_packer_task( *task );
		ig->set_rotamer_sets( *rotsets );

		rotsets->compute_energies( *trpcage, *sfxn, neighbors, ig );

		ig->prepare_for_simulated_annealing();

		ObjexxFCL::FArray1D_int state1( res_in_subnetwork.size(), leu_rot );
		for ( Size netres : res_in_subnetwork ) {
			trpcage->replace_residue( netres, *rotsets->rotamer_set_for_residue( netres )->rotamer( leu_rot ), false );
		}
		ig->set_network_state( state1 );
		(*sfxn)(*trpcage);

		// OK -- let's test some rotamer substitutions
		// Size count_output_pdbs( 0 );
		utility::vector1< Size > twos( 3, 2 );
		utility::LexicographicalIterator lex( twos );
		while ( ! lex.at_end() ) {
			for ( Size ii = 1; ii <= 3; ++ii ) {
				Size iires = res_in_subnetwork_v1[ ii ];

				// "leave one out" -- look at all combinations of leu/network hbonds
				if ( lex[ ii ] == 2 ) continue;

				Pose copy1( *trpcage );
				ObjexxFCL::FArray1D_int start_state( 3, leu_rot );
				for ( Size jj = 1; jj <= 3; ++jj ) {
					if ( ii == jj ) continue;
					Size jjres = res_in_subnetwork_v1[ jj ];
					if ( lex[ jj ] == 2 ) {
						copy1.replace_residue( jjres, *rotsets->rotamer_set_for_residue( jjres )->rotamer( don_net_rot ), false );
						start_state( jj ) = don_net_rot;
					}
				}
				ig->set_network_state( start_state );
				Real start_score( (*sfxn)( copy1 ) );
				copy1.replace_residue( iires, *rotsets->rotamer_set_for_residue( iires )->rotamer( don_net_rot ), false );
				Real end_score( (*sfxn )( copy1 ));
				float prev_e, deltaE;
				ig->consider_substitution( ii, don_net_rot, deltaE, prev_e );
				TS_ASSERT_DELTA( end_score - start_score, deltaE, 1e-5 );

			}
			++lex;
		}
	}

};
