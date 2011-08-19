// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/EtableEnergy.cxxtest.hh
/// @brief  Unit tests for the lennard-jones and EEF1 solvation model.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/types.hh>

// Unit headers
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>

// Package headers
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

#include <basic/Tracer.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <utility/stream_util.hh>
#include <string>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.etable.EtableEnergy.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// @name EtableEnergyTest
/// @brief: Test the functionality of the EtableEnergy class
///////////////////////////////////////////////////////////////////////////
class EtableEnergyTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_eval_residue_pair_energy()
	{
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		EtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );

		EnergyMap emap;
		ScoreFunction sfxn;
		etab_energy.residue_pair_energy( pose.residue( 4 ), pose.residue( 5 ), pose, sfxn, emap );

		//int const old_precision( TR.precision() );
		//TR.precision( 16 );
		//TR << "etable energy, trpcage 4 and 5: " << emap[ fa_atr ] << " " << emap[ fa_rep ] << " " << emap[ fa_sol ] << std::endl;
		//TR.precision( old_precision );

		//core.scoring.etable.EtableEnergy.cxxtest: etable energy, trpcage 4 and 5: -0.6079854186280232 0.0778507865338765 0.4275630856840362
		TS_ASSERT_DELTA( emap[ fa_atr ], -0.6079854186280232, 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_rep ],  0.0778507865338765, 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_sol ],  0.4275630856840362, 1e-12 );
	}

	void test_eval_residue_pair_energy_w_minimization_data()
	{
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		EtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );

		EnergyMap emap;
		ScoreFunction sfxn;
		etab_energy.residue_pair_energy( pose.residue( 4 ), pose.residue( 5 ), pose, sfxn, emap );

		optimization::MinimizerMap minmap;
		EnergyMap emap2;
		ResSingleMinimizationData r4dat, r5dat;
		ResPairMinimizationData min_data;
		etab_energy.setup_for_minimizing_for_residue( pose.residue( 4 ), pose, sfxn, minmap, r4dat );
		etab_energy.setup_for_minimizing_for_residue( pose.residue( 5 ), pose, sfxn, minmap, r5dat );
		etab_energy.setup_for_minimizing_for_residue_pair( pose.residue( 4 ), pose.residue( 5 ), pose, sfxn,  minmap, r4dat, r5dat, min_data );
		etab_energy.residue_pair_energy_ext( pose.residue( 4 ), pose.residue( 5 ), min_data, pose, sfxn, emap2 );

		TS_ASSERT_DELTA( emap[ fa_atr ], emap2[ fa_atr ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_rep ], emap2[ fa_rep ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_sol ], emap2[ fa_sol ], 1e-12 );
	}

	void test_eval_intra_residue_energy_w_minimization_data()
	{
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		EtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );

		EnergyMap emap;
		ScoreFunction sfxn;
		sfxn.set_weight( fa_intra_atr, 0.1 ); // enable intra-residue energy evaluation; skipped if all intra-residue weights are 0
		etab_energy.eval_intrares_energy( pose.residue( 4 ), pose, sfxn, emap );

		optimization::MinimizerMap minmap;
		EnergyMap emap2;
		ResSingleMinimizationData min_data;
		etab_energy.setup_for_minimizing_for_residue( pose.residue( 4 ), pose, sfxn,  minmap, min_data );
		etab_energy.eval_intrares_energy_ext( pose.residue( 4 ), min_data, pose, sfxn, emap2 );

		TS_ASSERT_DELTA( emap[ fa_intra_atr ], emap2[ fa_intra_atr ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_intra_rep ], emap2[ fa_intra_rep ], 1e-12 );
		TS_ASSERT_DELTA( emap[ fa_intra_sol ], emap2[ fa_intra_sol ], 1e-12 );

	}

	/*void dont_test_atom_deriv_old_vs_new()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		EtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );

		EnergyMap emap;
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );

		//sfxn.set_weight( fa_intra_rep, 0.02 );

		sfxn( pose );
		optimization::MinimizerMap minmap;
		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );
		minmap.setup( pose, movemap );

		MinimizationGraph g( pose.total_residue() );
		g.copy_connectivity( pose.energies().energy_graph() );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			etab_energy.setup_for_minimizing_for_residue(
				pose.residue( ii ), pose, sfxn, minmap, g.get_minimization_node( ii )->res_min_data() );
		}
		for ( Graph::EdgeListIter edge_iter = g.edge_list_begin(),
				edge_iter_end = g.edge_list_end(); edge_iter != edge_iter_end;
				++edge_iter ) {
			Size const node1 = (*edge_iter)->get_first_node_ind();
			Size const node2 = (*edge_iter)->get_second_node_ind();
			MinimizationEdge * minedge = dynamic_cast< MinimizationEdge * > ( *edge_iter );
			//std::cout << "Setting up minimization data for pair " << node1 << " " << node2 << std::endl;
			etab_energy.setup_for_minimizing_for_residue_pair(
				pose.residue( node1 ), pose.residue( node2 ), pose, sfxn,  minmap,
				g.get_minimization_node( node1 )->res_min_data(),
				g.get_minimization_node( node2 )->res_min_data(),
				minedge->res_pair_min_data() );
		}
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			MinimizationNode * iinode = g.get_minimization_node( ii );
			etab_energy.setup_for_minimizing_for_residue( pose.residue( ii ), pose, sfxn, minmap, iinode->res_min_data() );
		}

		/// Let's measure the derivatie for the CD2 atom on residue 6
		id::AtomID cd2( pose.residue(6).atom_index("CD2"), 6 );

		pose.energies().set_use_nblist( pose, minmap.domain_map(), false ); // false = autoupdate
		sfxn.setup_for_minimizing( pose, minmap );
		sfxn.setup_for_derivatives( pose );

		Vector goldF1(0.0), goldF2(0.0);
		sfxn.eval_atom_derivative( cd2, pose, minmap.domain_map(), goldF1, goldF2 );


		Vector decompF1(0.0), decompF2(0.0);
		EnergyMap weights;
		weights[ fa_atr ] = 0.5;
		weights[ fa_rep ] = 0.25;
		weights[ fa_sol ] = 0.125;
		//weights[ fa_intra_rep ] = 0.02;
		for ( Graph::EdgeListConstIter edge_iter = g.get_node( 6 )->const_edge_list_begin(),
				edge_iter_end = g.get_node( 6 )->const_edge_list_end(); edge_iter != edge_iter_end; ++edge_iter ) {
			MinimizationEdge const * minedge = dynamic_cast< MinimizationEdge const * > ( *edge_iter );
			Size const node2 = (*edge_iter)->get_other_ind( 6 );
			//std::cout << "Evaluating derivative for pair " << (*edge_iter)->get_first_node_ind() << " " << node2 << std::endl;
			etab_energy.eval_atom_derivative_for_residue_pair(
				cd2.atomno(), pose.residue( 6 ), pose.residue( node2 ),
				g.get_minimization_node( 6 )->res_min_data(), g.get_minimization_node( node2 )->res_min_data(),
				minedge->res_pair_min_data(), pose, minmap.domain_map(),
				sfxn, weights, decompF1, decompF2 );
		}
		etab_energy.eval_intrares_atom_derivative( cd2.atomno(), pose.residue( 6 ),
			g.get_minimization_node(6)->res_min_data(), pose, minmap.domain_map(),
			sfxn, weights, decompF1, decompF2 );

		TS_ASSERT_DELTA( goldF1.x(), decompF1.x(), 1e-12 );
		TS_ASSERT_DELTA( goldF1.y(), decompF1.y(), 1e-12 );
		TS_ASSERT_DELTA( goldF1.z(), decompF1.z(), 1e-12 );
		TS_ASSERT_DELTA( goldF2.x(), decompF2.x(), 1e-12 );
		TS_ASSERT_DELTA( goldF2.y(), decompF2.y(), 1e-12 );
		TS_ASSERT_DELTA( goldF2.z(), decompF2.z(), 1e-12 );
	}*/

	void test_start_func_matches_start_score_w_full_bbflex()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -25.94299600233435 );
	}
	void test_etab_numeric_deriv_check()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );

		adv.simple_deriv_check( true, 5e-3 );
	}

	/// @brief Note the derivative vectors computed below belong to the exact-derivative evaluation machinery
	/// and not the interpolated-derivative machinery.  Since we're continuing to use interpolated
	/// derivatives, these derivatives will all register as wrong.
	/*void dont_test_etable_derivs_w_full_bbflex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );

		//adv.compute_pose_atom_derivs(); // This line was used to generate the next set of lines:

		using namespace core;
		using namespace core::id;
		AtomDerivList adl;
		adl.add( AtomDeriv( AtomID( 1, 1), 0.08958673088646721,0.1821287577691115,-0.02467191655365314,-0.1267167031487318,0.05847556303004517,-0.02845557327366756));
		adl.add( AtomDeriv( AtomID( 2, 1), 0.01313409561604382,-0.1322652312771539,-0.351189772962506,-0.142248432900986,0.01142875745887183,-0.009624231762076366));
		adl.add( AtomDeriv( AtomID( 3, 1), -1.038253020435439,-1.029073919042521,2.323795191234375,0.2727793371592986,0.2135543166896935,0.2164464192970266));
		adl.add( AtomDeriv( AtomID( 4, 1), -0.3636883667504551,-0.5246389507286071,0.757449686465184,0.4617165298191914,-0.01105374770479745,0.2140362679278998));
		adl.add( AtomDeriv( AtomID( 5, 1), 0.2242177527630103,0.2972163599946891,-0.3939775038602061,-0.1329912946252242,0.006875939995982572,-0.07049988156333614));
		adl.add( AtomDeriv( AtomID( 6, 1), 0.03910495425028759,0.2146518889943101,0.100303935651133,-0.03267604302344655,0.01879322154501967,-0.02747853630590417));
		adl.add( AtomDeriv( AtomID( 7, 1), 0.01515065605532481,0.07148182000822254,0.01311383457570255,-0.01597204265231821,0.004942856882793981,-0.008490078218793412));
		adl.add( AtomDeriv( AtomID( 8, 1), 0.02643072988762418,-0.005423362043679335,-0.08997783287975986,-0.01899827313241162,-0.00209148614844731,-0.005454625024783314));
		adl.add( AtomDeriv( AtomID( 9, 1), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 10, 1), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 11, 1), 0.0884617028514188,0.1497562218571111,0.01698364456754008,0.01433747632246182,-0.006360521517302941,-0.01859376515734874));
		adl.add( AtomDeriv( AtomID( 12, 1), -0.1730499033115221,-1.005773689162697,-0.621886618238608,0.04730248334366361,-0.08121467605366189,0.1181853605208611));
		adl.add( AtomDeriv( AtomID( 13, 1), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 14, 1), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 15, 1), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 16, 1), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 2), -0.07364164620207672,0.5434463830781087,1.251728463819886,-0.04849740553418862,0.2256272229689221,-0.1008107833877967));
		adl.add( AtomDeriv( AtomID( 2, 2), -0.09329652394481736,0.2396968491813572,0.5777008452288551,-0.03800719642235138,0.1458384409903754,-0.06664860268770192));
		adl.add( AtomDeriv( AtomID( 3, 2), -0.3580926447053159,-2.336536678891831,-4.034561371479854,-0.7124526038851279,-0.4279958955096536,0.3111000254695575));
		adl.add( AtomDeriv( AtomID( 4, 2), 0.1187106427112947,-0.3885166006089103,-0.7102142389121755,-0.4504180551916486,0.08128820415114343,-0.1197543345473459));
		adl.add( AtomDeriv( AtomID( 5, 2), 0.3596652688014228,0.3929047070738018,0.2095125813531147,-0.03038462963918534,0.08477662477605546,-0.1068233649415618));
		adl.add( AtomDeriv( AtomID( 6, 2), -0.2328758051591025,-0.3002100934690699,-0.5443608067438241,-0.02514159551424499,-0.07944723628256704,0.05456993075291793));
		adl.add( AtomDeriv( AtomID( 7, 2), 1.297610392120572,0.5820229095396163,-1.45890363742245,-0.1106314100579755,-0.1949758465440468,-0.1761849585458727));
		adl.add( AtomDeriv( AtomID( 8, 2), -0.05976074192586889,0.2579177321651763,1.116372620579691,0.03093367267967449,0.1601170219792362,-0.03533623024507073));
		adl.add( AtomDeriv( AtomID( 9, 2), -0.06544218957276748,-0.1535033287041948,-0.1474223111839272,-0.01951010821788871,-0.007457896205870368,0.01642625240280308));
		adl.add( AtomDeriv( AtomID( 10, 2), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 11, 2), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 12, 2), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 13, 2), -0.3940882406701322,-0.491446404538916,-0.7875051236562226,0.08682411196222649,-0.2307006147206415,0.1005207759549681));
		adl.add( AtomDeriv( AtomID( 14, 2), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 15, 2), -0.05151338003006062,-0.02445991344277985,0.2031903691941558,0.1277159640931885,-0.1696798606352068,0.01195295967889673));
		adl.add( AtomDeriv( AtomID( 16, 2), 0.03694762070203393,0.06906384312019161,0.3302004468221989,0.04985845667475158,-0.0105636182178469,-0.003369429948796));
		adl.add( AtomDeriv( AtomID( 17, 2), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 18, 2), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 19, 2), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 3), 0.7519379953693255,1.024595359886998,1.137435672100456,-0.4291098518048824,0.5875900190516515,-0.2456209279606356));
		adl.add( AtomDeriv( AtomID( 2, 3), -0.2378890858842913,-0.4467457483551225,0.3732865046173544,-0.05383983282611252,0.1402309295731803,0.1335158446476045));
		adl.add( AtomDeriv( AtomID( 3, 3), 0.3942339791858468,1.758781495931073,-0.6209284140368777,-0.1832702892882811,-0.09359355479864458,-0.3814642435043402));
		adl.add( AtomDeriv( AtomID( 4, 3), -0.156740242792483,-0.1684204743874812,-0.3221645654111227,-0.0524563530085104,-0.09032509304114056,0.07274113618218339));
		adl.add( AtomDeriv( AtomID( 5, 3), -0.6552392170676312,-1.024986430760521,0.4167103036400238,0.01359483912047071,0.09331176225972224,0.2508965124566283));
		adl.add( AtomDeriv( AtomID( 6, 3), -0.172309583401373,0.01521716630406212,-0.1656733006596874,-0.1594008929393974,0.1300848479535027,0.1777342764578992));
		adl.add( AtomDeriv( AtomID( 7, 3), -0.09985248760372152,-0.07300096977112241,-0.02608838031842699,-0.1045757929129598,0.07719722053371708,0.1842460530023404));
		adl.add( AtomDeriv( AtomID( 8, 3), -0.9886381048589845,0.2638854119815348,-0.8947366142415994,-0.2174271974196162,0.02157374444327204,0.2466086726864921));
		adl.add( AtomDeriv( AtomID( 9, 3), -0.3415963929241937,-0.7004834990944748,0.1867397628298539,0.1549777447775477,-0.08745673909297993,-0.04456557025452165));
		adl.add( AtomDeriv( AtomID( 10, 3), -1.056972068530434,0.8169955857214577,-1.01726268845573,-0.2523460324586381,0.05404812062605101,0.3056042332125783));
		adl.add( AtomDeriv( AtomID( 11, 3), -0.9498302654413734,-0.2565189119939917,0.02758038003993056,0.04155825775409024,-0.1465178513781708,0.06847951942244258));
		adl.add( AtomDeriv( AtomID( 12, 3), -0.2366663506499247,0.305161342398713,-0.1416229370988149,-0.04797151217647502,-0.008578837929507976,0.06168007246651162));
		adl.add( AtomDeriv( AtomID( 13, 3), 0.6783145767522059,0.8144672572516209,1.040687209488278,-0.03784433154978174,0.2346375015857509,-0.1589663052634284));
		adl.add( AtomDeriv( AtomID( 14, 3), 0.006645983604119448,0.007370109243716431,-0.003876015435345033,-0.0008625398331536492,-0.0005648310491337137,-0.002552954778033968));
		adl.add( AtomDeriv( AtomID( 15, 3), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 16, 3), 0.0006862488284215601,0.01196709822340477,-0.01212025011066616,-0.0005557711177442126,-0.002151337525140892,-0.002155620924892293));
		adl.add( AtomDeriv( AtomID( 17, 3), -0.3200483705765229,-0.07614237831829591,-0.2467240062966292,0.0001479505454088344,-0.07037514276613813,0.02152680435698619));
		adl.add( AtomDeriv( AtomID( 18, 3), 0.0357052889256847,0.007689929330244131,0.0105290047386416,0.002023986638293168,6.726239583451689e-05,-0.006912739862810956));
		adl.add( AtomDeriv( AtomID( 19, 3), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 20, 3), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 21, 3), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 4), 0.6784503338621983,1.975818963869844,2.460325504994235,-0.02202501269603423,0.4617385542531395,-0.3647358501160426));
		adl.add( AtomDeriv( AtomID( 2, 4), -0.01516324492342256,0.1900898164690547,-0.02857292931096747,-0.2305371912395893,-0.01973098418047029,-0.008923385632629486));
		adl.add( AtomDeriv( AtomID( 3, 4), -0.2400933443008526,0.9165162234648185,0.3355012499883968,0.5160075654185959,0.1842711122764236,-0.1341201615580086));
		adl.add( AtomDeriv( AtomID( 4, 4), -0.08867545520749008,0.1109548041534378,0.702449372005472,-0.2335340655408258,0.0345857868642441,-0.03494372655814983));
		adl.add( AtomDeriv( AtomID( 5, 4), -0.0102389851154303,0.1475623770450024,-0.04091301652400484,-0.06452275744996168,-0.00937743556480946,-0.01767430492368309));
		adl.add( AtomDeriv( AtomID( 6, 4), -0.06365395979976039,0.2057744713636616,-0.5665091022093046,-0.1683400945418602,-0.05774163519125147,-0.002058644498295906));
		adl.add( AtomDeriv( AtomID( 7, 4), 0.09885430984905288,-0.5492698102071191,-0.4659952925563294,-0.1174577795783915,-0.08639868359475233,0.07692133670916332));
		adl.add( AtomDeriv( AtomID( 8, 4), -0.3835221067379756,-0.5602737043080308,-1.074926196966013,-0.2565800470808359,-0.1224868482585412,0.1553876730120113));
		adl.add( AtomDeriv( AtomID( 9, 4), 0.1445711324721675,0.3362774126049267,0.5852819667595099,0.05178822156061252,0.08220206585905163,-0.06002197547978379));
		adl.add( AtomDeriv( AtomID( 10, 4), 0.001549653786652389,-0.002251621990312821,0.004397695511994591,0.001897062810165875,0.001079878957981087,-0.0001155858465503708));
		adl.add( AtomDeriv( AtomID( 11, 4), -0.1686737255579723,1.002548786005925,1.455102981749616,0.0969040965354754,0.1957954414439529,-0.1236677468176789));
		adl.add( AtomDeriv( AtomID( 12, 4), 0.1650194580281748,0.6463460927361643,0.3304032870493718,-0.02114948058487586,0.0466278074790097,-0.0806518167273453));
		adl.add( AtomDeriv( AtomID( 13, 4), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 14, 4), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 15, 4), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 16, 4), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 17, 4), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 18, 4), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 19, 4), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 5), -0.2988708405541932,-0.7283233701863441,1.149243372249374,0.02384866769157635,0.2390753963115649,0.1577140874648646));
		adl.add( AtomDeriv( AtomID( 2, 5), -0.0623003241809924,0.7529145977124581,-0.2347140631372836,0.01081222970228794,-0.05365647857821006,-0.1749888815777048));
		adl.add( AtomDeriv( AtomID( 3, 5), 0.07666904411360723,-0.04807582189326871,-0.05897954830249846,-0.0174708591842168,-0.03120606397022877,0.002726082253134582));
		adl.add( AtomDeriv( AtomID( 4, 5), 2.571962483506588,-1.743405110082832,-0.4135859390250755,-0.5721503393751999,-0.8600143153835589,0.06722700595443273));
		adl.add( AtomDeriv( AtomID( 5, 5), 0.05864456835344859,1.176684327071584,-0.164056988213993,0.03738314443371783,-0.03571064218966413,-0.2427682906946846));
		adl.add( AtomDeriv( AtomID( 6, 5), -0.4420501733193337,0.4120186844782897,0.2622253940464486,-0.01381947337115303,0.05870414834392537,-0.1155346021431326));
		adl.add( AtomDeriv( AtomID( 7, 5), -0.345993736517278,0.4898358679179395,0.1809216605335777,-0.03863381381193592,0.01987733772882552,-0.1276999697375313));
		adl.add( AtomDeriv( AtomID( 8, 5), 2.884317601171079,2.89983345123384,-4.75135973877537,0.5000960272348105,-0.5693118208156183,-0.04387663742939671));
		adl.add( AtomDeriv( AtomID( 9, 5), -0.2374429504901756,-0.003520344342778162,0.2534700590166595,-0.02454285486559135,0.03502847496576013,-0.022504494618362));
		adl.add( AtomDeriv( AtomID( 10, 5), -0.0731819920318251,-0.1914149477803336,0.1756259448189903,0.006063337013865194,0.03328743929971754,0.03880657008482304));
		adl.add( AtomDeriv( AtomID( 11, 5), 1.259798887956005,-0.1118889928107435,-1.906164735084106,0.03270505994800633,-0.4028895828621262,0.04526403526692555));
		adl.add( AtomDeriv( AtomID( 12, 5), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 13, 5), 0.002248107352783552,-0.002973157295006521,-0.002201770328217693,0.001102220260134824,-0.0007634122850179242,0.002156289516054068));
		adl.add( AtomDeriv( AtomID( 14, 5), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 15, 5), 0.4036288110419861,-0.003602602777429414,-0.2884659379486441,0.001818296898500712,-0.08152649563277695,0.003562377597062628));
		adl.add( AtomDeriv( AtomID( 16, 5), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 17, 5), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 6), 0.176687080387397,-0.8602621133909742,0.04267314138858454,-0.444935314432684,-0.08684575640981101,0.09149098456721609));
		adl.add( AtomDeriv( AtomID( 2, 6), -0.06872949159302087,0.1828355302490002,-0.06362645494233111,0.1532262837096738,0.0450646163414779,-0.03601884712686175));
		adl.add( AtomDeriv( AtomID( 3, 6), -0.05992125308915881,0.0213353447164176,0.06138357123397745,-0.01366782327976952,0.05762806500153633,-0.03337224747423872));
		adl.add( AtomDeriv( AtomID( 4, 6), 0.08454076543947278,0.01187811857808684,-0.06693297964520997,-0.01540456228096821,0.313059452196372,0.03609945084611588));
		adl.add( AtomDeriv( AtomID( 5, 6), 0.1145823309588809,0.1630358010387771,0.196093527409052,0.325991498012873,-0.07366592824941548,-0.129237728692691));
		adl.add( AtomDeriv( AtomID( 6, 6), 0.1254389803863799,-0.1768862485975857,-0.1494682446454895,-0.0513790109575589,0.1101517588987655,-0.1734767289018781));
		adl.add( AtomDeriv( AtomID( 7, 6), 0.249998948768422,-0.5803659643613481,0.4631749217276971,0.06609113058204122,-0.1889403393621447,-0.2724181502712673));
		adl.add( AtomDeriv( AtomID( 8, 6), -0.3106678555072777,-0.3699213246961113,0.4068071149217193,0.1117943069337315,-0.2421862901192801,-0.1348525472119443));
		adl.add( AtomDeriv( AtomID( 9, 6), 0.1908622120124963,-1.619430040920232,0.4917457261387844,0.3016752645659259,-0.1026579311781775,-0.4551656152675946));
		adl.add( AtomDeriv( AtomID( 10, 6), -0.3369047443485477,0.6320501836430164,0.3918331863529977,-0.1954799952557722,-0.1716462611350711,0.1087986279108883));
		adl.add( AtomDeriv( AtomID( 11, 6), -0.3677771833775565,0.404712811926473,0.01398337496051153,-0.08631890789577107,-0.08669270645803109,0.2388210425196591));
		adl.add( AtomDeriv( AtomID( 12, 6), -2.584527726754407,0.4515798131083444,3.107362009062089,-0.05254522125641399,-0.8905309936104918,0.08571297378272261));
		adl.add( AtomDeriv( AtomID( 13, 6), 0.6402084585885791,1.923668787387808,-0.4472507739421661,-0.3710085165812317,0.1988702908774152,0.3242873779681085));
		adl.add( AtomDeriv( AtomID( 14, 6), -1.373671203329626,-0.2342268377733169,0.9606799215047135,0.1491790156641856,-0.3274753402622007,0.1334673512898704));
		adl.add( AtomDeriv( AtomID( 15, 6), 0.03152270172342855,-0.3082198182288761,-0.1415258239344526,-0.1741015966957925,-0.03535442190409301,0.03821762446094276));
		adl.add( AtomDeriv( AtomID( 16, 6), 0.008734201621918354,-0.01646617639733152,0.007368549773348499,-0.009201842506955094,-0.006941992296264667,-0.004605698949483233));
		adl.add( AtomDeriv( AtomID( 17, 6), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 18, 6), -0.03006165431926123,0.002692931368867108,-0.07262236707386779,-0.0622240388104351,-0.01382756418009669,0.02524457598396962));
		adl.add( AtomDeriv( AtomID( 19, 6), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 20, 6), 0.05759478591363058,-0.5600263551791395,0.1576034212863581,0.2095512011042681,-0.01374033253236479,-0.1254034002008357));
		adl.add( AtomDeriv( AtomID( 21, 6), -0.1103335290807858,0.2506752718779794,-0.06739899380323919,-0.1002079483743582,-0.02393636273449537,0.07501658503402729));
		adl.add( AtomDeriv( AtomID( 22, 6), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 23, 6), 0.2737695678819826,0.3366579190253379,-0.07579095166219332,-0.06817853131926735,0.05612219753981284,0.003018502496154915));
		adl.add( AtomDeriv( AtomID( 24, 6), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 7), -0.08852766696918825,0.2395888423574791,0.3632490494571925,-0.1499645599075818,0.04677655611580459,-0.06740046142590066));
		adl.add( AtomDeriv( AtomID( 2, 7), 0.3655205703575232,-0.09358365387747841,0.1633904104196398,-0.1058582749990939,-0.09025126440319364,0.1851230674256469));
		adl.add( AtomDeriv( AtomID( 3, 7), -10.66822183312729,2.066487571166231,-7.32526850128551,1.193739108852847,-1.228448398398061,-2.085063089298384));
		adl.add( AtomDeriv( AtomID( 4, 7), -1.173488426780546,0.04209220299027792,-0.5981449288191755,0.07238835322540878,-0.1966975001080868,-0.1558590925948189));
		adl.add( AtomDeriv( AtomID( 5, 7), 0.1441557468359077,-0.5700109060572955,-0.3737494769221433,-0.06978223524287153,-0.2064979590982761,0.2880180045033207));
		adl.add( AtomDeriv( AtomID( 6, 7), -0.1156083796460985,-0.1056182122929276,-0.1083163037439827,-0.04927400069972657,-0.0801920485595308,0.1307857423017759));
		adl.add( AtomDeriv( AtomID( 7, 7), -0.1138560063148491,-0.942957560598206,-0.3371846914914517,-0.01318973216115849,-0.1013428545935418,0.287865207493495));
		adl.add( AtomDeriv( AtomID( 8, 7), 0.2866540301163428,0.08743912718799489,0.1037050216584047,-0.0664879803652853,-0.03034567099901488,0.2093673591199095));
		adl.add( AtomDeriv( AtomID( 9, 7), -0.07244469257184404,0.1747603751906684,0.06436401726775573,-0.008115523337901482,0.02204249135810975,-0.06898389568636448));
		adl.add( AtomDeriv( AtomID( 10, 7), -0.01448279722967578,-0.001889061326849134,-0.006056889269676148,0.002285309520681788,0.0008543433120571461,-0.005730925523027964));
		adl.add( AtomDeriv( AtomID( 11, 7), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 12, 7), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 13, 7), 0.3002514037831872,-0.1914877682949534,0.1388525672050734,0.0784057691707406,0.09346413539155821,-0.04064889594826535));
		adl.add( AtomDeriv( AtomID( 14, 7), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 15, 7), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 16, 7), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 17, 7), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 18, 7), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 19, 7), 0.04156178356183329,-0.02541897679464268,-0.01086527884383126,0.006244124415389993,0.01066905455670129,-0.001074938147454922));
		adl.add( AtomDeriv( AtomID( 1, 8), -0.2113683888762971,0.06398159165225892,-0.7317499513931749,0.1343825136519483,-0.0500909836322884,-0.04319660862812487));
		adl.add( AtomDeriv( AtomID( 2, 8), 0.2597584698468707,-0.132809458124603,-0.01740536865076826,-0.08352337307614167,-0.1684977938459324,0.03919463717485544));
		adl.add( AtomDeriv( AtomID( 3, 8), 2.984294836341832,-0.2964388492995966,-3.713663071353501,0.5443563347707038,-0.1460149084291444,0.4490995166703807));
		adl.add( AtomDeriv( AtomID( 4, 8), -0.3412837228612389,0.02298581947175426,0.3533470691466471,-0.06594194686625189,-0.08704902694779394,-0.05802798916695967));
		adl.add( AtomDeriv( AtomID( 5, 8), 0.3347329088889932,-0.1946839926511182,-0.2320055811766531,-0.1624299034618587,-0.2793503271817621,6.208026428615888e-05));
		adl.add( AtomDeriv( AtomID( 6, 8), -0.747318772619459,0.2617053068115619,0.6127183789955598,0.03937655105701775,0.2142641835141192,-0.04349018895606925));
		adl.add( AtomDeriv( AtomID( 7, 8), 0.4234855930337771,0.3659221115748787,-1.009970077495869,0.03023510258263293,-0.2312858489067479,-0.07111940983419476));
		adl.add( AtomDeriv( AtomID( 8, 8), -1.050175743497831,1.73559480118695,-1.074057630572538,0.3027374038228277,0.1418253360870126,-0.06682710509376404));
		adl.add( AtomDeriv( AtomID( 9, 8), -0.9561454560523069,-1.071218715186675,1.733557057916432,-0.1385339534579621,0.2133152537013762,0.05540554979330006));
		adl.add( AtomDeriv( AtomID( 10, 8), 0.07665578466917866,-0.06593335957156471,-0.0004179083257564904,0.03698303767756603,0.04286557454305807,0.02079985608979888));
		adl.add( AtomDeriv( AtomID( 11, 8), 0.002256296697972348,-0.001371319172434929,-0.002612964376520623,0.0005640064119070987,0.0003169731870381028,0.0003206681284357159));
		adl.add( AtomDeriv( AtomID( 12, 8), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 13, 8), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 14, 8), -0.9147265091787038,0.1141672504626341,1.358986310576893,0.004198489586262316,0.3015101926525946,-0.02250362621059721));
		adl.add( AtomDeriv( AtomID( 15, 8), -2.38719600964327,1.570587920393812,-0.04297702675084447,0.5593810656782536,0.8476883151534953,-0.09268253574534571));
		adl.add( AtomDeriv( AtomID( 16, 8), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 17, 8), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 18, 8), -0.9462344042803468,-0.482039784343277,1.097053828301704,-0.1276596871314149,0.1666420276140676,-0.03688779885344796));
		adl.add( AtomDeriv( AtomID( 19, 8), 1.460977991871222,3.946294865473181,-4.64692836193422,0.9304672363923701,0.17174295756297,0.4383843148574998));
		adl.add( AtomDeriv( AtomID( 20, 8), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 21, 8), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 22, 8), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 9), 8.82921683561727,-2.265231816224718,4.751247106458741,-0.7852904609186561,0.5597687831392706,1.726179597045198));
		adl.add( AtomDeriv( AtomID( 2, 9), -0.1669532944435822,0.1138680953582734,-0.4151094224388727,0.0823964100060318,-0.1587292764779679,-0.07667990840574428));
		adl.add( AtomDeriv( AtomID( 3, 9), -1.625253609065312,-0.04315762786706984,-2.79044203116581,0.4810382921956427,-0.1821160303978925,-0.2773573204449247));
		adl.add( AtomDeriv( AtomID( 4, 9), -5.480016002876077,-0.302951949730942,-10.78458084477986,1.28167664115946,0.7533464802674221,-0.6724263458599243));
		adl.add( AtomDeriv( AtomID( 5, 9), -0.2346859123037976,0.2690794811221221,-0.571721900205883,0.167001339259625,-0.2262335232578935,-0.1750285596453961));
		adl.add( AtomDeriv( AtomID( 6, 9), -0.2592343069950416,-0.01734485967574118,-0.02130570530289078,0.02385891512770897,-0.1438235428985721,-0.173214174638975));
		adl.add( AtomDeriv( AtomID( 7, 9), -2.265351839665199,-7.015613190471083,8.224669389371947,-1.643706811561426,-0.2715035516623874,-0.6843239382991263));
		adl.add( AtomDeriv( AtomID( 8, 9), 0.2009324997203767,-0.9114410425950024,0.9326232682056854,-0.1639009282686617,-0.09541303139920829,-0.05793371389508979));
		adl.add( AtomDeriv( AtomID( 9, 9), 2.643524312263445,-0.8881193746912136,1.645882705774965,-0.3521378514009713,0.06006211194138365,0.5979935831959974));
		adl.add( AtomDeriv( AtomID( 10, 9), 0.1699559724789146,-0.06929371530988,0.214310343308506,-0.03087543510202386,-0.01149734872341349,0.02076787578473728));
		adl.add( AtomDeriv( AtomID( 11, 9), 0.007777247139242849,0.00103267408301424,0.001059054644578169,-0.0006532128444872727,0.001744540162448838,0.003095832996297256));
		adl.add( AtomDeriv( AtomID( 12, 9), 0.3394045726084037,0.6865825073707378,-1.05773495913434,0.23787105385564,-0.06811190337202463,0.03211587333619701));
		adl.add( AtomDeriv( AtomID( 1, 10), 0.5618156556977153,0.4021216025073708,3.661521792781784,-0.63430851506138,0.1779165678112358,0.07778742693938751));
		adl.add( AtomDeriv( AtomID( 2, 10), 0.213574151222596,0.2994320101287142,0.9646750845412697,-0.1149825200918477,-0.1041600336885479,0.05778748025092491));
		adl.add( AtomDeriv( AtomID( 3, 10), 0.5773323632800725,0.7140863521632415,0.9992628234940537,-0.1145556166011534,-0.1429113390058891,0.1683117771005821));
		adl.add( AtomDeriv( AtomID( 4, 10), 0.185463308052699,-0.3815953305391472,-0.8863537120160678,0.0796778537386148,0.1237335453679184,-0.03659805828168157));
		adl.add( AtomDeriv( AtomID( 5, 10), -0.9097893180839821,0.3883654962872808,3.617168480989811,-0.6641855498713666,0.2138507804943127,-0.1900163585437196));
		adl.add( AtomDeriv( AtomID( 6, 10), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 7, 10), -0.6356083876397565,-0.239265340364442,0.03002071846742068,-0.008901454362702553,0.01236097156043412,-0.08994744714102069));
		adl.add( AtomDeriv( AtomID( 1, 11), 0.03735601047464545,0.0661214321507146,0.09577941759440019,0.1521846535850826,-0.3356368092699237,0.1723520085278118));
		adl.add( AtomDeriv( AtomID( 2, 11), 1.556163999907484,0.1895839155051959,-1.299438777730494,0.1616696475007083,0.304624542632641,0.2380538461456995));
		adl.add( AtomDeriv( AtomID( 3, 11), 0.4662862429489061,0.6798142907790149,-0.002588221664967166,0.9730354705127721,-0.6644775604811665,0.7695292529638509));
		adl.add( AtomDeriv( AtomID( 4, 11), 0.4225221438709026,1.198943379052845,0.9420204235542478,0.02308903925713583,-0.2145621225604236,0.262724883329717));
		adl.add( AtomDeriv( AtomID( 5, 11), -1.009795855730008,-0.1304356349255242,0.7412928618483761,-0.05876346948363244,-0.3336800177164992,-0.1387614507297659));
		adl.add( AtomDeriv( AtomID( 6, 11), 0.9586347616270897,-0.2748532149454453,-1.139492328009496,0.05587633973843914,0.4706842811762688,-0.06652443760832308));
		adl.add( AtomDeriv( AtomID( 7, 11), 6.304459109559227,-1.68635133110381,-5.384968958152478,1.184122621462921,0.7036479669658233,1.165959360358113));
		adl.add( AtomDeriv( AtomID( 1, 12), 0.3004810571028704,0.3838634779063482,-0.109610667944847,0.1109000022472245,-0.04002033673618433,0.1638618265738445));
		adl.add( AtomDeriv( AtomID( 2, 12), 0.372625555808952,0.4488540767081955,-0.3245868222191327,0.09393745718956213,0.01353971409328065,0.1265635270831742));
		adl.add( AtomDeriv( AtomID( 3, 12), -22.13588131140835,-30.04831727356888,14.417083415325,-0.6609278677194499,-1.703842224183725,-4.565959055327549));
		adl.add( AtomDeriv( AtomID( 4, 12), 0.4845801061028043,1.133522339285006,0.03701628628880903,0.07367040287654293,-0.03679846204676496,0.1624330191277407));
		adl.add( AtomDeriv( AtomID( 5, 12), 0.5783840785352289,0.2659380609049521,-0.6177895286117125,0.09711522319236675,0.05721895846324837,0.1155516473369115));
		adl.add( AtomDeriv( AtomID( 6, 12), 0.3553243372273199,0.2005679805513944,-0.2583804215700096,0.09034802291386243,0.0140833104036895,0.1351786341882313));
		adl.add( AtomDeriv( AtomID( 7, 12), -1.050254700420957,1.609916550973376,1.740783455148325,-0.6146754466816351,-0.02082236264473519,-0.3515908363263134));
		adl.add( AtomDeriv( AtomID( 8, 12), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 9, 12), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 10, 12), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 11, 12), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 12, 12), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 13, 12), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 14, 12), -0.00285928672920752,-0.001507872247411754,0.001115047783243694,-0.0006375089652046594,0.0003036615141681561,-0.00122410732046718));
		adl.add( AtomDeriv( AtomID( 15, 12), -3.916444904004648,0.09607933158417897,2.626907976925555,-0.6011955700538771,-0.4025449536145331,-0.8815965754850035));
		adl.add( AtomDeriv( AtomID( 1, 13), 1.702471027354817,1.348418350601046,-1.391326221120487,-0.8846998271591618,0.8918300789943003,-0.2182204107861057));
		adl.add( AtomDeriv( AtomID( 2, 13), -0.1877329261203957,-0.1746541399893689,0.0966478658519852,0.1160746807940266,-0.1141949219636729,0.01910464934696973));
		adl.add( AtomDeriv( AtomID( 3, 13), -0.6918261011449622,-0.9232493815188453,-0.8813341775125112,0.301678330145334,-0.1470986842783492,-0.08271570026695398));
		adl.add( AtomDeriv( AtomID( 4, 13), -0.3306496499949613,-0.2984488384493626,0.1758478494993143,0.03666242689806703,-0.06359892150218459,-0.03900306797425013));
		adl.add( AtomDeriv( AtomID( 5, 13), -0.3484735155174486,-0.02900058005501442,0.6017873475963476,0.08488399179858255,-0.1886569249244194,0.04006176413528367));
		adl.add( AtomDeriv( AtomID( 6, 13), -0.215618113920491,-0.01831868577030926,0.2737309310173394,0.05539281453604421,-0.100847507449275,0.03688403191339629));
		adl.add( AtomDeriv( AtomID( 7, 13), 1.031684191821098,0.245428589538492,-1.278872003497666,-0.1974552788544389,0.4122449943433876,-0.08017595349598813));
		adl.add( AtomDeriv( AtomID( 8, 13), 0.008530205094986455,-0.006726853258468792,-0.03743588323244587,-1.231257045171725e-05,0.004624249673938021,-0.0008337369134448927));
		adl.add( AtomDeriv( AtomID( 9, 13), 0.2044862276199192,0.2238069211265432,0.1682237404046878,-0.07266679187912954,0.06762699824880285,-0.00164086304243195));
		adl.add( AtomDeriv( AtomID( 10, 13), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 11, 13), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 14), 7.724072417453533,9.122649101894833,-9.451531751392437,0.7273687409539027,0.8463132090413974,1.411291589190086));
		adl.add( AtomDeriv( AtomID( 2, 14), -0.2419412281651836,-0.3480223828482487,0.2800931755723721,0.1685360004267227,-0.1814376585457868,-0.07986113651827528));
		adl.add( AtomDeriv( AtomID( 3, 14), -1.461369171687849,0.6375779361412226,-7.680401398940963,0.2968157768014474,0.8770050450445707,0.01632748527102153));
		adl.add( AtomDeriv( AtomID( 4, 14), 0.1479981899955668,0.5630275924474173,-0.4076674397363539,0.07647658614010357,0.01507787527228874,0.04858778064876738));
		adl.add( AtomDeriv( AtomID( 5, 14), -0.5298567067400407,-1.266411127333827,2.936585900492932,-0.2340582601330186,-0.3577095436220801,-0.1964950813116119));
		adl.add( AtomDeriv( AtomID( 6, 14), 3.07079172600758,-0.03485926921433452,5.960906377802939,-0.7177230206124554,-0.4197911693099473,0.3672837922819491));
		adl.add( AtomDeriv( AtomID( 7, 14), 14.94663406082702,21.69867245709663,-7.469178137651745,0.4426798405969494,0.8293740941982518,3.29526086328623));
		adl.add( AtomDeriv( AtomID( 8, 14), 0.3456127915773505,0.3372095471710272,0.10388932275851,0.01228261523259937,-0.02769562461981064,0.04903487637605815));
		adl.add( AtomDeriv( AtomID( 9, 14), -0.07316685616636993,-0.09122094570604795,0.437720660745513,-0.04910854954714736,-0.04130052562916189,-0.01681572712324764));
		adl.add( AtomDeriv( AtomID( 10, 14), -1.280462449421133,-0.5893161311301182,-2.699730826973799,0.1227772933858897,0.457917398509037,-0.1581897792260356));
		adl.add( AtomDeriv( AtomID( 11, 14), 3.694247333610815,0.3793469278312429,5.678112256584719,-0.7547964679920298,-0.2480058436763308,0.5076484866671147));
		adl.add( AtomDeriv( AtomID( 1, 15), -0.3285682839427466,-0.5531025917441792,-0.5781785348402371,0.1405057739642983,0.004846436339559806,-0.08448310442574258));
		adl.add( AtomDeriv( AtomID( 2, 15), -0.165977178488858,-0.3186005369262037,-0.3286811969015422,0.1004983708253035,0.002160530277018161,-0.05284385691164753));
		adl.add( AtomDeriv( AtomID( 3, 15), 0.1837469117716549,-0.1889125949028607,0.8719725338940245,0.1314383425934477,-0.1232449353804314,-0.05439839931592199));
		adl.add( AtomDeriv( AtomID( 4, 15), -0.1678509834840849,0.06041668937731451,-0.5125210551525557,0.07757622611805659,0.04505835751977662,-0.02009472382692582));
		adl.add( AtomDeriv( AtomID( 5, 15), 0.32854245384353,0.6280399288822874,2.36020994835296,-0.2434312194427882,-0.1633389196373864,0.07734941282576484));
		adl.add( AtomDeriv( AtomID( 6, 15), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 7, 15), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 16), 1.091650517938819,-1.243512483269547,6.522487761611499,-0.0542582000351506,-0.8181076477092719,-0.1468912040054559));
		adl.add( AtomDeriv( AtomID( 2, 16), -0.1574114948246968,-0.2544404264783457,-0.5868163015103213,0.07837057144092992,0.08151051450254072,-0.05636516706872116));
		adl.add( AtomDeriv( AtomID( 3, 16), 0.01044517667414647,-0.6517764183401441,-0.6406749136051237,-0.06887033433630341,0.08599521661473689,-0.08860814725809218));
		adl.add( AtomDeriv( AtomID( 4, 16), 0.08917447359957791,-0.8955263712214487,-1.457696308859393,0.006723045358372876,0.2351244708535244,-0.1440359276907887));
		adl.add( AtomDeriv( AtomID( 5, 16), -0.207885018446635,-1.04264536210021,-0.3544797262823148,0.05340130497303724,0.05637079090035753,-0.1971229093759852));
		adl.add( AtomDeriv( AtomID( 6, 16), 0.0598895442479053,-0.8122533594819399,0.518919154875463,-0.07059980116205433,-0.08738048873272117,-0.1286267908664061));
		adl.add( AtomDeriv( AtomID( 7, 16), -0.8378980497349953,-0.9068508130120361,-0.5709637382719587,-0.0100865539064814,0.147430079984949,-0.219358210784196));
		adl.add( AtomDeriv( AtomID( 8, 16), -0.9365920117043919,-0.4660048563086669,-0.6541529360559812,0.001139772487312672,0.1536022020248832,-0.1110548770568554));
		adl.add( AtomDeriv( AtomID( 9, 16), -0.9164067048196118,-0.4156271627289856,-0.6144548709168114,0.05897040105972373,0.09130830895222476,-0.1497117016310243));
		adl.add( AtomDeriv( AtomID( 10, 16), -0.3389108771764102,-0.5702436377900009,-0.1220887492796988,-0.03333380945712048,0.03503748901647395,-0.07111805669452517));
		adl.add( AtomDeriv( AtomID( 11, 16), -1.119367626474311,-0.817508209333685,-0.3939233446358574,0.0474717463797891,0.03939541255586376,-0.2166523269219632));
		adl.add( AtomDeriv( AtomID( 12, 16), 0.05800492244556581,-0.3273999412570584,1.376185248622343,0.04235345930652788,-0.1910426425894021,-0.04723496284298469));
		adl.add( AtomDeriv( AtomID( 13, 16), -2.125642975007124,1.605867289242894,-6.308658684552529,0.002858474978252305,0.7956827096725591,0.2015776732281337));
		adl.add( AtomDeriv( AtomID( 14, 16), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 15, 16), 0.01079780261882736,0.0410512392470169,0.01148253942263518,0.005060434011321179,-0.002263655682321766,0.003334132107793554));
		adl.add( AtomDeriv( AtomID( 16, 16), 1.32812108612535,0.6146036908914012,2.949888421181757,0.008022747392112589,-0.541235459547202,0.1091533255247327));
		adl.add( AtomDeriv( AtomID( 17, 16), -0.239028658341125,1.582323003408824,-1.368914557038202,-0.4365993312362113,0.3696038072174517,0.503459003443445));
		adl.add( AtomDeriv( AtomID( 18, 16), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 19, 16), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 20, 16), 0.01617327185210707,0.01577340343416909,0.001581023041264311,0.0010776676710932,-0.001396511022349831,0.002908445623469984));
		adl.add( AtomDeriv( AtomID( 21, 16), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 22, 16), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 23, 16), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 24, 16), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 1, 17), 0.01709862927555459,-0.2174650862348108,-0.4460406952947044,0.1756378572166395,0.1167895081113127,-0.05020724360285407));
		adl.add( AtomDeriv( AtomID( 2, 17), -0.05036259379821408,0.1351583043801213,-0.2097079435549041,0.1116470391623495,0.08168839027897111,0.02583607346338844));
		adl.add( AtomDeriv( AtomID( 3, 17), -0.08148805877448446,0.03527833338573065,0.8704917618704033,0.1002631419451234,-0.1032140022066973,0.01356872896566318));
		adl.add( AtomDeriv( AtomID( 4, 17), 0.1484205374041597,-0.03129012065870886,0.8729691803802487,0.2057211303550763,-0.04958184684832191,-0.03675348845543634));
		adl.add( AtomDeriv( AtomID( 5, 17), -0.789107895006551,1.396733408729407,-0.1859945150132656,0.01353529305417538,0.03234259480484909,0.185452652064463));
		adl.add( AtomDeriv( AtomID( 6, 17), 0.2135994011911075,-0.3729500580990511,0.02012522494016553,0.1021753536541398,0.05492100976265485,-0.06667257529605421));
		adl.add( AtomDeriv( AtomID( 7, 17), 1.751858025271757,-1.707318410632396,4.746734641769034,0.07464076783143475,-0.5810452245319245,-0.236539289032784));
		adl.add( AtomDeriv( AtomID( 8, 17), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 9, 17), -0.369416604450204,1.878570541788247,-1.718658057399713,-0.2139744993114232,0.1696739764628876,0.2314539911840626));
		adl.add( AtomDeriv( AtomID( 10, 17), -0.1723027153926766,0.2340479442686705,0.004194780160536512,-0.02071358926237178,-0.01581916065398411,0.0318096182824763));
		adl.add( AtomDeriv( AtomID( 11, 17), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 12, 17), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 13, 17), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 14, 17), 0.2224614299797059,0.1901899531693124,0.7643851700706192,0.02109464908401795,-0.09619865247796699,0.01779635704906207));
		adl.add( AtomDeriv( AtomID( 15, 17), 0.2257716147056588,-0.3024790016031212,0.643341206979845,-0.04156429553573508,-0.08482755580680318,-0.02529686594780585));
		adl.add( AtomDeriv( AtomID( 1, 18), 0.17320165127825,0.1141400026536342,-0.832554033783761,0.1880424036256555,0.3280447183138392,0.08409337652111029));
		adl.add( AtomDeriv( AtomID( 2, 18), 0.02013713286105723,-0.1247586786469226,0.2664714292304207,0.06037209113783362,-0.006241414644017823,-0.007484447657212265));
		adl.add( AtomDeriv( AtomID( 3, 18), 0.5626677145660757,-0.2620531362719178,-0.243146289182517,0.0695079132064525,0.2121178821502959,-0.06776289970990443));
		adl.add( AtomDeriv( AtomID( 4, 18), -0.1743108114286235,0.09333387595518262,0.05494064654955707,0.003876850297193132,-0.008288419341230754,0.0263806000616075));
		adl.add( AtomDeriv( AtomID( 5, 18), -0.3706113524320216,-0.2963167168493726,0.7518046829436369,0.229412316594187,0.06703184535138951,0.1395115880972461));
		adl.add( AtomDeriv( AtomID( 6, 18), 0.08815957601846217,-0.08036964575382743,-0.005656760074487455,0.1914221052631812,0.2000409044439638,0.1411541242005381));
		adl.add( AtomDeriv( AtomID( 7, 18), 0.4706743568340843,-1.432094194670444,1.109074432796944,0.3457112582060988,0.04655527590341946,-0.08660003414890217));
		adl.add( AtomDeriv( AtomID( 8, 18), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 9, 18), 1.315560209908188,-0.05030362060269983,-2.715142272525945,-0.7263460385007268,0.1207141148796548,-0.3541708711336434));
		adl.add( AtomDeriv( AtomID( 10, 18), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 11, 18), -0.9065577152233103,0.01365789312930078,0.6873156917380642,0.005429115257047609,-0.240340334228602,0.01193679850655069));
		adl.add( AtomDeriv( AtomID( 12, 18), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 13, 18), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 14, 18), 0.8454310661127505,-1.53087811626669,0.6114099528492447,0.08764556156179537,-0.0334972487765234,-0.2050645153453606));
		adl.add( AtomDeriv( AtomID( 15, 18), 0.1622410848983626,-0.5225861717008106,0.1939736394383352,0.03983973897903773,-0.01314758474271754,-0.06874330187568836));
		adl.add( AtomDeriv( AtomID( 1, 19), 0.6570131896889005,-0.165992170232332,-0.6415114131763857,-0.02109791916531087,0.3092226542746896,-0.1016196271736994));
		adl.add( AtomDeriv( AtomID( 2, 19), 0.4205103598195471,-0.07478122377754427,0.05634627520905846,0.01454542838327634,0.03213190361108027,-0.06590746655504363));
		adl.add( AtomDeriv( AtomID( 3, 19), 0.6302588227532011,-0.08917681062756812,0.3216425182593681,0.05412214078924549,0.1203341180969858,-0.0726892203296621));
		adl.add( AtomDeriv( AtomID( 4, 19), 0.5596310815677805,0.009698841604767562,-0.01387556847691172,-0.00392195515285747,0.1656305769926307,-0.04240714705673399));
		adl.add( AtomDeriv( AtomID( 5, 19), 0.5233163275728089,0.1343556935371346,0.6994855748847056,0.1171871616702875,0.1017962747561063,-0.1072257482803585));
		adl.add( AtomDeriv( AtomID( 6, 19), 0.6436706845534792,0.07670783268966477,0.1722738241098364,-0.0002079078298411231,0.3734520717002289,-0.1655090377529036));
		adl.add( AtomDeriv( AtomID( 7, 19), -0.3367168316305637,-0.09969581751955761,2.863949776089468,0.6384417539392051,-0.2126796006654767,0.06765859496795104));
		adl.add( AtomDeriv( AtomID( 8, 19), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 9, 19), -0.07762275221781899,0.02627613306836571,0.03206035952427931,0.006323012073707895,0.006257884210846353,0.01018006678839385));
		adl.add( AtomDeriv( AtomID( 10, 19), -0.1245965489991851,-0.144625724761106,-0.6657352163432696,-0.09930579568909648,0.009972748683704843,0.01641920565999047));
		adl.add( AtomDeriv( AtomID( 11, 19), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 12, 19), 0.1105292852341068,0.0681345031287361,-0.6758559067274577,-0.1439790724347716,0.0558111559078983,-0.01791985313364115));
		adl.add( AtomDeriv( AtomID( 13, 19), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 14, 19), -0.06411249617706993,-0.1091380862699261,0.366814029924709,0.09396898232482366,-0.01521284012070141,0.01189781580452499));
		adl.add( AtomDeriv( AtomID( 15, 19), -0.318379527275172,0.1627265835171113,0.5017015570046214,0.08346695525239117,-0.1835833781911118,0.1125132359902238));
		adl.add( AtomDeriv( AtomID( 1, 20), 0.1786942482356135,0.001003592599970518,-0.4965641541405323,-0.03884994592812342,0.1160984019610949,-0.01374595070471386));
		adl.add( AtomDeriv( AtomID( 2, 20), -0.06342154628688898,0.02473704784709901,-0.1533998435922527,-0.008896473542104022,0.07463144409831264,0.01571311714310264));
		adl.add( AtomDeriv( AtomID( 3, 20), -0.303195721337871,-0.02286307191518862,0.1857421563079399,0.01458778178090551,0.07873486904536066,0.03350383196090913));
		adl.add( AtomDeriv( AtomID( 4, 20), -0.375307564203911,-0.09928961890496434,0.5155534902651573,0.0481849625654906,0.02119705915982887,0.03915952318975257));
		adl.add( AtomDeriv( AtomID( 5, 20), 0.07340700210240951,0.001989237441200961,-0.03116152168739558,-0.00363888982870768,0.03196980682571038,-0.00653127464770243));
		adl.add( AtomDeriv( AtomID( 6, 20), -0.01311420220394584,0.01917518754596183,-0.1418362866634252,-0.002101059609421389,0.06060454557258256,0.008387538027384987));
		adl.add( AtomDeriv( AtomID( 7, 20), -0.05192769155698423,0.01568744748132186,-0.00190804544986202,0.01407210676190731,0.04825698709498849,0.01378213038483194));
		adl.add( AtomDeriv( AtomID( 8, 20), 0.03473683150254361,-0.007941163958207436,0.01353367679280487,-0.0005288315151479559,-0.008700428057912475,-0.003747805955178999));
		adl.add( AtomDeriv( AtomID( 9, 20), 0.00589348375504029,0.002979401092312876,-0.01440876591044504,-0.001882850274506107,-0.00562030806940073,-0.001932275094211892));
		adl.add( AtomDeriv( AtomID( 10, 20), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 11, 20), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 12, 20), 0,0,0,0,0,0));

		adv.validate_atom_deriv_list( adl );
	}*/

	/// @brief Numeric deriv check
	void dont_test_atom_tree_minimize_with_etable_energy()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -25.94299600233435, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, true, false );
		min_options.deriv_check_to_stdout( false );

		minimizer.run( pose, movemap, sfxn, min_options );

		NumericalDerivCheckResultOP deriv_check_result = minimizer.deriv_check_result();
		for ( Size ii = 1, iiend = deriv_check_result->n_deriv_check_results(); ii <= iiend; ++ii ) {
			NumDerivCheckData const & iidata( deriv_check_result->deriv_check_result( ii ) );
			TS_ASSERT( iidata.nsteps() >= 1 );
			for ( Size jj = 1; jj <= iidata.nangles(); ++jj ) {
				/// 1e-1.  That's not very good.  That's the problem with interpolating derivatives instead of using the real
				/// derivatives.
				TS_ASSERT_DELTA( iidata.dof_step_data( jj, 1 ).num_deriv(), iidata.dof_step_data( jj, 1 ).ana_deriv(), 1e-1 );
			}
		}

		//Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		//TS_ASSERT_DELTA( -26.48575563731416, end_score, 1e-12 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_start_func_matches_start_score_w_partial_bbflex()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -25.94299600233435 );
	}


	void dont_test_atom_tree_minimize_with_etable_energy2()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -25.94299600233435, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -26.15005763826676, end_score, 1e-12 );
	}

	void dont_test_atom_tree_minimize_with_intrares_etable_energy()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -22.82696911663496, end_score, 1e-12 );
	}

	void dont_test_atom_tree_minimize_with_intrares_etable_energy2()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -22.49034811587382, end_score, 1e-12 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_start_func_matches_start_score_w_full_bbflex_and_intraresidue()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -22.28334391305056 );
	}

	/// @brief Make sure that the domain map logic inside the ScoreFunction
	/// operates correctly for the intraresidue portions of two-body energies.
	void test_start_func_matches_start_score_w_partial_bbflex_and_intraresidue()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -22.28334391305056 );
	}

	/*void dont_test_etable_derivatives_w_intraresidue_terms_and_full_bb_flex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );

		adv.add_res_for_deriv( 9 );
		adv.add_res_for_deriv( 10 );
		adv.add_res_for_deriv( 11 );

		//adv.compute_pose_atom_derivs();
		using namespace core;
		using namespace core::id;
		AtomDerivList adl;
		adl.add( AtomDeriv( AtomID( 1, 9), 8.870229654729471,-2.230899786804663,4.568061214041253,-0.7458987923648381,0.6037756366965,1.743246456980084));
		adl.add( AtomDeriv( AtomID( 2, 9), -0.9194829750509632,-0.01342076116459836,-0.1956472707984394,0.04251969824598709,-0.1209947218103617,-0.1915298752953238));
		adl.add( AtomDeriv( AtomID( 3, 9), -5.780766298022331,-0.7551964912427289,-4.974801883984178,0.717216633071126,0.2720736746243307,-0.8747143960212155));
		adl.add( AtomDeriv( AtomID( 4, 9), -7.038245358969004,-0.8584241219216738,-10.9061353063898,1.205515878812623,1.007445317287538,-0.8572726853677711));
		adl.add( AtomDeriv( AtomID( 5, 9), 0.4525475150942307,0.4327799782120202,-0.5970654298387187,0.1888594631229107,-0.3027068221211068,-0.07626897969322977));
		adl.add( AtomDeriv( AtomID( 6, 9), 0.9345010564613776,-0.2008774428784887,0.3896422196996073,-0.04311210834188221,-0.2587926038509379,-0.03002058071297762));
		adl.add( AtomDeriv( AtomID( 7, 9), -2.270610777728292,-6.899999576711538,8.093815861886982,-1.619062369783125,-0.2827566838621938,-0.6952569173340969));
		adl.add( AtomDeriv( AtomID( 8, 9), -0.005257280783539625,-0.8857584219791812,0.844962576430577,-0.1510101487017299,-0.07859238493299497,-0.08332649465748718));
		adl.add( AtomDeriv( AtomID( 9, 9), 2.497870843644781,-0.8194163669972006,1.466598503343748,-0.3137795281242261,0.05352129872726993,0.5643241561390084));
		adl.add( AtomDeriv( AtomID( 10, 9), -0.8850939766765746,-0.02143177917411314,-0.04341772278233186,0.001749226624900419,0.1103921361443422,-0.09015050958726513));
		adl.add( AtomDeriv( AtomID( 11, 9), 0.2580272027175695,-0.06666602938643937,0.1991272534364095,-0.04138573474506205,-0.04796844850768407,0.03756783283928272));
		adl.add( AtomDeriv( AtomID( 12, 9), 6.045596869061917,1.876719981417103,1.28375698969165,0.05566573369119911,-0.737902691073536,0.8165907938802381));
		adl.add( AtomDeriv( AtomID( 1, 10), 0.5066602428212209,0.3954511161799208,3.711577108686203,-0.642193853663245,0.1762894873724578,0.06888182888205918));
		adl.add( AtomDeriv( AtomID( 2, 10), 0.213574151222596,0.2994320101287142,0.9646750845412697,-0.1149825200918477,-0.1041600336885479,0.05778748025092491));
		adl.add( AtomDeriv( AtomID( 3, 10), 0.5436543430687971,0.7176719489125598,1.045702965011548,-0.1228751692583829,-0.1434294992976939,0.1623185105086581));
		adl.add( AtomDeriv( AtomID( 4, 10), 0.1975081726244682,-0.3869741033353062,-0.9094329280474898,0.08550297842543522,0.121017254358682,-0.03292492007638309));
		adl.add( AtomDeriv( AtomID( 5, 10), -0.8761112978727067,0.3847798995379625,3.570728339472317,-0.6558659972141372,0.2143689407861175,-0.1840230919517956));
		adl.add( AtomDeriv( AtomID( 6, 10), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 7, 10), -0.5924978393350311,-0.2272160812408331,0.00304461859442277,-0.006841240447657951,0.01670434300844861,-0.08471498728899086));
		adl.add( AtomDeriv( AtomID( 1, 11), 0.03960381876528875,0.03551869886455637,0.03053933495700529,0.1642270767585205,-0.3289754528862154,0.1696422226578495));
		adl.add( AtomDeriv( AtomID( 2, 11), 1.556163999907484,0.1895839155051959,-1.299438777730494,0.1616696475007083,0.304624542632641,0.2380538461456995));
		adl.add( AtomDeriv( AtomID( 3, 11), 0.4662862429489061,0.6798142907790149,-0.002588221664967166,0.9730354705127721,-0.6644775604811665,0.7695292529638509));
		adl.add( AtomDeriv( AtomID( 4, 11), 0.4202743355802593,1.229546112339003,1.007260506191643,0.01104661608369788,-0.2212234789441319,0.2654346691996792));
		adl.add( AtomDeriv( AtomID( 5, 11), -1.009795855730008,-0.1304356349255242,0.7412928618483761,-0.05876346948363244,-0.3336800177164992,-0.1387614507297659));
		adl.add( AtomDeriv( AtomID( 6, 11), 0.9586347616270897,-0.2748532149454453,-1.139492328009496,0.05587633973843914,0.4706842811762688,-0.06652443760832308));
		adl.add( AtomDeriv( AtomID( 7, 11), 6.304459109559227,-1.68635133110381,-5.384968958152478,1.184122621462921,0.7036479669658233,1.165959360358113));

		adv.validate_atom_deriv_list( adl );
	}*/

	/*void dont_test_etable_derivatives_w_intraresidue_terms_and_partial_bb_flex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );

		adv.add_res_for_deriv( 9 );
		adv.add_res_for_deriv( 10 );
		adv.add_res_for_deriv( 11 );

		//adv.compute_pose_atom_derivs();

		using namespace core;
		using namespace core::id;
		AtomDerivList adl;
		adl.add( AtomDeriv( AtomID( 1, 9), 1.189437192888129,-0.05915294988937453,-0.4856878332655125,0.08483931280373207,-0.02437570755682031,0.2107380956761624));
		adl.add( AtomDeriv( AtomID( 2, 9), -0.1143457467976166,-0.1307893868449699,0.3970805867101087,-0.06461247364021053,-0.0261948023788113,-0.02723417878994933));
		adl.add( AtomDeriv( AtomID( 3, 9), -1.220153558527872,-0.0754365423249008,-1.740839960141839,0.2721569522379748,0.005169885374109654,-0.1909786537764307));
		adl.add( AtomDeriv( AtomID( 4, 9), -5.416137163355387,-0.3245726381831474,-10.50104002784313,1.236263952555874,0.7652978423830882,-0.6612840117207842));
		adl.add( AtomDeriv( AtomID( 5, 9), 0.09809889676819149,-0.4790407769384676,0.8936189721289668,-0.1635990134816922,-0.09423234913715206,-0.03255554762805684));
		adl.add( AtomDeriv( AtomID( 6, 9), 0.194074702550198,-0.8984445775649909,1.098436461899415,-0.2104986779487496,-0.07594688003298976,-0.0249277907173227));
		adl.add( AtomDeriv( AtomID( 7, 9), -0.02557681891386881,0.001373208934434602,0.001252185887329203,-0.0002406511123768966,5.787597742414968e-05,-0.004978945853854591));
		adl.add( AtomDeriv( AtomID( 8, 9), 0.281644782975596,-1.065918870204982,1.104972755710165,-0.1992535214137278,-0.09723493945958744,-0.04301087223365593));
		adl.add( AtomDeriv( AtomID( 9, 9), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 10, 9), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 11, 9), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 12, 9), 0.3394045726084037,0.6865825073707378,-1.05773495913434,0.23787105385564,-0.06811190337202463,0.03211587333619701));
		adl.add( AtomDeriv( AtomID( 1, 10), 0.5066602428212209,0.3954511161799208,3.711577108686203,-0.642193853663245,0.1762894873724578,0.06888182888205918));
		adl.add( AtomDeriv( AtomID( 2, 10), 0.213574151222596,0.2994320101287142,0.9646750845412697,-0.1149825200918477,-0.1041600336885479,0.05778748025092491));
		adl.add( AtomDeriv( AtomID( 3, 10), 0.5436543430687971,0.7176719489125598,1.045702965011548,-0.1228751692583829,-0.1434294992976939,0.1623185105086581));
		adl.add( AtomDeriv( AtomID( 4, 10), 0.1975081726244682,-0.3869741033353062,-0.9094329280474898,0.08550297842543522,0.121017254358682,-0.03292492007638309));
		adl.add( AtomDeriv( AtomID( 5, 10), -0.8761112978727067,0.3847798995379625,3.570728339472317,-0.6558659972141372,0.2143689407861175,-0.1840230919517956));
		adl.add( AtomDeriv( AtomID( 6, 10), 0,0,0,0,0,0));
		adl.add( AtomDeriv( AtomID( 7, 10), -0.5924978393350311,-0.2272160812408331,0.00304461859442277,-0.006841240447657951,0.01670434300844861,-0.08471498728899086));
		adl.add( AtomDeriv( AtomID( 1, 11), 0.03960381876528875,0.03551869886455637,0.03053933495700529,0.1642270767585205,-0.3289754528862154,0.1696422226578495));
		adl.add( AtomDeriv( AtomID( 2, 11), 1.556163999907484,0.1895839155051959,-1.299438777730494,0.1616696475007083,0.304624542632641,0.2380538461456995));
		adl.add( AtomDeriv( AtomID( 3, 11), 0.4662862429489061,0.6798142907790149,-0.002588221664967166,0.9730354705127721,-0.6644775604811665,0.7695292529638509));
		adl.add( AtomDeriv( AtomID( 4, 11), 0.4202743355802593,1.229546112339003,1.007260506191643,0.01104661608369788,-0.2212234789441319,0.2654346691996792));
		adl.add( AtomDeriv( AtomID( 5, 11), -1.009795855730008,-0.1304356349255242,0.7412928618483761,-0.05876346948363244,-0.3336800177164992,-0.1387614507297659));
		adl.add( AtomDeriv( AtomID( 6, 11), 0.9586347616270897,-0.2748532149454453,-1.139492328009496,0.05587633973843914,0.4706842811762688,-0.06652443760832308));
		adl.add( AtomDeriv( AtomID( 7, 11), 6.304459109559227,-1.68635133110381,-5.384968958152478,1.184122621462921,0.7036479669658233,1.165959360358113));

		adv.validate_atom_deriv_list( adl );
	}*/

	void test_setup_for_minimizing_with_autoupdate_w_full_bb_flexibility()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );
		adv.validate_start_func_matches_start_score( -22.28334391305056 );
	}

	void test_setup_for_minimizing_with_autoupdate_w_partial_bb_flexibility()
	{

		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );
		adv.validate_start_func_matches_start_score( -22.28334391305056 );

		/*using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		/// BEGIN AtomTreeMinimizer setup block
		MinimizerMap min_map;
		min_map.setup( pose, movemap );

		pose.energies().set_use_nblist( pose, min_map.domain_map(), min_options.nblist_auto_update() );
		sfxn.setup_for_minimizing( pose, min_map );
		/// END AtomTreeMinimizer setup block

		Real start_func = sfxn(pose);
		//std::cout << "start_func: " << start_func << std::endl;
		TS_ASSERT_DELTA( -22.28334391305056, start_func, 1e-12 );*/
	}

	void test_etable_derivatives_w_autoupdate_intraresidue_terms_and_full_bb_flex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );

		adv.add_res_for_deriv( 9 );
		adv.add_res_for_deriv( 10 );
		adv.add_res_for_deriv( 11 );

		adv.simple_deriv_check( true, 5e-3 );

	}

	void test_etable_derivatives_w_autoupdate_intraresidue_terms_and_partial_bb_flex()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nblist_auto_update( true );

		adv.simple_deriv_check( true, 5e-3 );
	}

	void dont_test_atom_tree_minimize_with_autoupdate()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -22.82696911663496, end_score, 1e-12 );
	}

	void test_setup_for_minimizing_with_autoupdate2()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout.precision( 16 );
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		/// BEGIN AtomTreeMinimizer setup block
		MinimizerMap min_map;
		min_map.setup( pose, movemap );

		pose.energies().set_use_nblist( pose, min_map.domain_map(), min_options.nblist_auto_update() );
		sfxn.setup_for_minimizing( pose, min_map );
		/// END AtomTreeMinimizer setup block

		Real start_func = sfxn(pose);
		//std::cout << "start_func: " << start_func << std::endl;
		TS_ASSERT_DELTA( -22.28334391305056, start_func, 1e-12 );
	}

	void dont_test_atom_tree_minimize_with_autoupdate2()
	{
		using namespace core;
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.25 );
		sfxn.set_weight( fa_sol, 0.125 );
		sfxn.set_weight( fa_intra_atr, 0.5 );
		sfxn.set_weight( fa_intra_rep, 0.25 );
		sfxn.set_weight( fa_intra_sol, 0.125 );

		AtomTreeMinimizer minimizer;
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		Real start_score = sfxn(pose);
		TS_ASSERT_DELTA( -22.28334391305056, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );
		min_options.nblist_auto_update( true );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -22.49034811587382, end_score, 1e-12 );
	}

};
