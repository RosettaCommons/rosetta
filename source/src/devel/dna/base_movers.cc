// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


#include <devel/dna/base_movers.hh>
#include <devel/dna/relax_util.hh>

#include <devel/cartesian_frags/DNA_FragLib.hh>
#include <devel/cartesian_frags/dna_util.hh>

#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/conformation/Residue.hh>

#include <basic/prof.hh> // profiling
#include <basic/Tracer.hh>

#include <core/chemical/ResidueType.hh>

#include <numeric/random/random.hh>

//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>


namespace devel {
namespace dna {

using namespace core;
using namespace devel::cartesian_frags;
using utility::vector1;

static THREAD_LOCAL basic::Tracer tt( "devel.dna.base_movers", basic::t_trace );
static THREAD_LOCAL basic::Tracer td( "devel.dna.base_movers", basic::t_debug );
static THREAD_LOCAL basic::Tracer ti( "devel.dna.base_movers", basic::t_info );
static THREAD_LOCAL basic::Tracer tw( "devel.dna.base_movers", basic::t_warning );


/// @details  Setup a pointgraph for later use in dma calcs
conformation::PointGraphOP
setup_dme_point_graph( pose::Pose const & ref_pose, Real const threshold )
{
	conformation::PointGraphOP pg( new conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( ref_pose.conformation(), *pg );
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, threshold );
	return pg;
}

/// @details  Calculate the dme using a pointgraph
Real
point_graph_dme( conformation::PointGraph const & pg, pose::Pose const & pose )
{
	Size total(0);
	Real dme(0.0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & i_rsd( pose.residue(i) );
		for ( auto
				i_iter     = pg.get_vertex( i ).const_upper_edge_list_begin(),
				i_end_iter = pg.get_vertex( i ).const_upper_edge_list_end();
				i_iter != i_end_iter; ++i_iter ) {
			Size const j = i_iter->upper_vertex();
			Real const reference_distance( std::sqrt( i_iter->data().dsq() ) );
			Real const pose_distance( i_rsd.nbr_atom_xyz().distance( pose.residue(j).nbr_atom_xyz() ) );
			dme += ( reference_distance - pose_distance ) * ( reference_distance - pose_distance );
			++total;
		}
	}
	tt << "dme_nbrs: " << total << std::endl;
	return std::sqrt( dme / total );
}

///////////////////////////////////////////////////////////////////////////////
///
/// @details  Make a base-pair fragment insertion between residue seqpos and seqpos-partner
/// try to patch up the backbone after making the insertion

void
make_base_pair_move(
	Size const seqpos,
	DNA_FragLib const & lib,
	Real const frag_dev_threshold,
	Size const max_tries,
	Real const max_score_increase,
	scoring::ScoreFunction const & scorefxn,
	pose::Pose & pose_inout
)
{
	PROF_START( basic::MAKE_BASE_PAIR_MOVE );

	//// scorefxn is most likely a chainbreak scoring function.
	using namespace conformation;
	using namespace kinematics;
	using namespace scoring::dna;
	using namespace id;

	pose::Pose const start_pose( pose_inout );
	pose::Pose pose( pose_inout ); // our local pose

	// test the use of rsd-rsd distances to gauge the impact of a fragment insertion
	conformation::PointGraphOP pg( setup_dme_point_graph( pose, 10.0 ) );

	Real const start_score( scorefxn( pose ) );

	Size const seqpos_partner( retrieve_base_partner_from_pose( pose )[ seqpos ] );
	chemical::ResidueType const &  rsd_type( pose.residue_type( seqpos         ) );
	chemical::ResidueType const & prsd_type( pose.residue_type( seqpos_partner ) );
	std::string const bp( std::string() + rsd_type.name1() + prsd_type.name1() );
	vector1< CartesianFragment > const & bps( lib.base_pairs(bp) );
	vector1< std::pair< Real, Size > > frag_devs;

	Size top_nn( 0 );
	{ // compute deviations from current for all frags
		Stub const chi  ( torsion_stub( TorsionID( seqpos        , CHI, 1 ), Backward, pose.conformation() ) );
		Stub const chi_p( torsion_stub( TorsionID( seqpos_partner, CHI, 1 ), Backward, pose.conformation() ) );
		RT const current( RT( chi, chi_p ) );
		for ( Size i=1; i<= bps.size(); ++i ) {
			frag_devs.push_back( std::make_pair( current.distance_squared( bps[i].rt(1) ), i ) );
		}
		std::sort( frag_devs.begin(), frag_devs.end() );
		td << "bp_move: frag_devs: " << frag_devs[1].first << ' ' << frag_devs[2].first << std::endl;
		top_nn = std::min( Size(5), frag_devs.size() );
		while ( top_nn < frag_devs.size() && frag_devs[ top_nn ].first < frag_dev_threshold ) ++top_nn;
	}

	// need these for making the fragment insertion
	vector1< Size > offsets;
	offsets.push_back( seqpos );
	offsets.push_back( seqpos_partner );

	Size ntries( 0 );
	Real best_score( 999.9 ), dme(0.0), rmsd(0.0), frag_dev(0.0);
	while ( ntries < max_tries ) { // keep looping until we like the move
		++ntries;

		// make a random basepair fragment insertion
		Size const nn( static_cast< int >( top_nn * numeric::random::rg().uniform() ) + 1 );
		bps[ frag_devs[ nn ].second ].insert( pose.conformation(), offsets );


		dme = point_graph_dme( *pg, pose );
		rmsd = scoring::nbr_atom_rmsd( pose, start_pose );
		frag_dev = frag_devs[ nn ].first;

		tt << "bp_move: choose_frag: nn= " << nn << " top_nn= " <<top_nn << " max_nn= " << bps.size() <<
			" max_frag_dev= " << frag_dev_threshold << " frag_dev= " << frag_dev <<
			" dme= " << dme << " rmsd= " << rmsd << std::endl;

		// by inserting this fragment we have perturbed some of the backbone segments


		// did we perturb the seqpos-1 --> seqpos   connection?
		if ( !rsd_type.is_lower_terminus() && !pose.fold_tree().jump_exists( seqpos-1, seqpos ) ) {
			patch_up_backbone_link( seqpos-1, lib, scorefxn, pose );
		}
		// did we perturb the seqpos   --> seqpos+1 connection?
		if ( !rsd_type.is_upper_terminus() && !pose.fold_tree().jump_exists( seqpos, seqpos+1 ) ) {
			patch_up_backbone_link( seqpos  , lib, scorefxn, pose );
		}

		if ( !prsd_type.is_lower_terminus() && !pose.fold_tree().jump_exists( seqpos_partner-1, seqpos_partner ) ) {
			patch_up_backbone_link( seqpos_partner-1, lib, scorefxn, pose );
		}
		if ( !prsd_type.is_upper_terminus() && !pose.fold_tree().jump_exists( seqpos_partner, seqpos_partner+1 ) ) {
			patch_up_backbone_link( seqpos_partner  , lib, scorefxn, pose );
		}

		Real const final_score( scorefxn( pose ) );
		tt << "bp_move: ntries= " << ntries << " score_increase= " << final_score - start_score <<
			" max_score_increase= " << max_score_increase << " dme,rmsd " << dme << ' ' << rmsd << std::endl;
		scorefxn.show( tt, pose );

		if ( ntries == 1 || final_score < best_score ) {
			best_score = final_score;
			pose_inout = pose;
		}

		if ( final_score - start_score < max_score_increase ) break;

		pose = start_pose; // unmake the move
	}
	td << "make_base_pair_move: score_increase= " << best_score - start_score << " dme= " << dme << " rmsd= " <<
		rmsd << " frag_dev= " << frag_dev << " ntries= " << ntries << std::endl;
	PROF_STOP( basic::MAKE_BASE_PAIR_MOVE );

}


///////////////////////////////////////////////////////////////////////////////
///
/// @details  Make a base-step fragment insertion between residue seqpos and seqpos+1
/// try to patch up the backbone after making the insertion

void
make_base_step_move(
	Size const seqpos,
	DNA_FragLib const & lib,
	Real const frag_dev_threshold,
	Size const max_tries,
	Real const max_score_increase,
	scoring::ScoreFunction const & scorefxn,
	pose::Pose & pose_inout
)
{
	//// scorefxn is most likely a chainbreak scoring function.
	PROF_START( basic::MAKE_BASE_STEP_MOVE );

	using namespace conformation;
	using namespace kinematics;
	using namespace scoring::dna;
	using namespace id;

	pose::Pose const start_pose( pose_inout );
	pose::Pose pose( pose_inout ); // our local pose
	assert( pose.fold_tree().jump_exists( seqpos, seqpos+1 ) );

	// test the use of rsd-rsd distances to gauge the impact of a fragment insertion
	conformation::PointGraphOP pg( setup_dme_point_graph( pose, 10.0 ) );

	Real const start_score( scorefxn( pose ) );

	Size const seqpos_partner( retrieve_base_partner_from_pose( pose )[ seqpos ] );
	chemical::ResidueType const &      rsd_type( pose.residue_type( seqpos   ) );
	chemical::ResidueType const & next_rsd_type( pose.residue_type( seqpos+1 ) );
	//Conformation & conf( pose.conformation() );

	std::string const bs( std::string() + rsd_type.name1() + next_rsd_type.name1() );
	vector1< CartesianFragment > const & bss( lib.base_steps( bs ) );
	vector1< std::pair< Real, Size > > frag_devs;

	Size top_nn( 0 );
	{ // compute deviations from current for all frags
		Stub const chi  ( torsion_stub( TorsionID( seqpos  , CHI, 1 ), Backward, pose.conformation() ) );
		Stub const chi_p( torsion_stub( TorsionID( seqpos+1, CHI, 1 ), Backward, pose.conformation() ) );
		RT const current( RT( chi, chi_p ) );
		for ( Size i=1; i<= bss.size(); ++i ) {
			frag_devs.push_back( std::make_pair( current.distance_squared( bss[i].rt(1) ), i ) );
		}
		std::sort( frag_devs.begin(), frag_devs.end() );
		tt << "bs_move: frag_devs: " << frag_devs[1].first << ' ' << frag_devs[2].first << std::endl;
		top_nn = std::min( Size(5), frag_devs.size() );
		while ( top_nn < frag_devs.size() && frag_devs[ top_nn ].first < frag_dev_threshold ) ++top_nn;
	}

	// need these for making the fragment insertion
	vector1< Size > offsets;
	offsets.push_back( seqpos   );
	offsets.push_back( seqpos+1 );

	Size ntries( 0 );
	Real best_score( 999.9 ), dme(0.0), rmsd(0.0), frag_dev(0.0);

	while ( ntries < max_tries ) { // keep looping until we like the move
		++ntries;

		// make a random basepair fragment insertion
		Size const nn( static_cast< int >( top_nn * numeric::random::rg().uniform() ) + 1 );
		bss[ frag_devs[ nn ].second ].insert( pose.conformation(), offsets );


		dme = point_graph_dme( *pg, pose );
		rmsd = scoring::nbr_atom_rmsd( pose, start_pose );
		frag_dev = frag_devs[ nn ].first;

		tt << "bs_move: choose_frag: nn= " << nn << " top_nn= " <<top_nn << " max_nn= " << bss.size() <<
			" max_frag_dev= " << frag_dev_threshold << " frag_dev= " << frag_dev <<
			" dme= " << dme << " rmsd= " << rmsd << std::endl;

		// by inserting this fragment we have perturbed two of the backbone segments
		patch_up_backbone_link( seqpos, lib, scorefxn, pose );

		patch_up_backbone_link( seqpos_partner-1, lib, scorefxn, pose );

		Real const final_score( scorefxn( pose ) );
		tt << "bs_move: ntries= " << ntries << " score_increase= " << final_score - start_score <<
			" max_score_increase= " << max_score_increase << std::endl;

		scorefxn.show( tt, pose );

		if ( ntries == 1 || final_score < best_score ) {
			best_score = final_score;
			pose_inout = pose;
		}

		if ( final_score - start_score < max_score_increase ) break;

		pose = start_pose; // unmake the move
	}
	td << "make_base_step_move: score_increase= " << best_score - start_score << " dme= " << dme << " rmsd= " <<
		rmsd << " frag_dev= " << frag_dev << " ntries= " << ntries << std::endl;

	PROF_STOP ( basic::MAKE_BASE_STEP_MOVE );
}


/// @details  Make a basepair move using make_base_pair_move
void
BasePairMover::apply( pose::Pose & pose )
{
	make_base_pair_move(
		choose_random_base_pair( pose ), *lib_, frag_dev_threshold_, max_tries_,
		max_score_increase_, *scorefxn_, pose );
}

std::string
BasePairMover::get_name() const {
	return "BasePairMover";
}


/// @details  Make a basestep move using make_base_step_move
void
BaseStepMover::apply( core::pose::Pose & pose )
{
	make_base_step_move(
		choose_random_base_step_jump( pose ), *lib_, frag_dev_threshold_, max_tries_,
		max_score_increase_, *scorefxn_, pose );
}

std::string
BaseStepMover::get_name() const {
	return "BaseStepMover";
}

} // ns dna
} // ns devel
