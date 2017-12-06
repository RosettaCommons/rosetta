// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/RNA_Motif.cc
/// @brief
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/rna/RNA_Motif.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.rna.RNA_Motif" );

using namespace core::pose;

namespace core {
namespace scoring {
namespace rna {

RNA_Motifs
get_rna_motifs( pose::Pose const & pose,
	core::scoring::rna::RNA_LowResolutionPotential const & potential,
	core::pose::rna::RNA_FilteredBaseBaseInfo const & filtered_base_base_info )
{
	// TODO:
	// 1. tandem GA sheared and WC should probably be generalized to any purine
	// residues. The only reason they haven't been? We currently arbitrarily use
	// 'which residue is the G' to enforce uniqueness. (Otherwise we would find
	// G-A and A-G in both positions.)

	using namespace core::pose::rna;
	using namespace core::scoring::rna;
	using namespace core::chemical::rna;
	using namespace core::chemical;
	using namespace utility::tools;

	RNA_Motifs rna_motifs;

	// we should pre-cache this information.
	utility::vector1< core::pose::rna::BasePair > base_pairs;
	RNA_BasePairList const & base_pair_list  = filtered_base_base_info.base_pair_list();
	for ( auto const & base_pair : base_pair_list ) {
		base_pairs.push_back( base_pair );
		base_pairs.push_back( base_pair.flipped() );
	}
	utility::vector1< core::pose::rna::BaseStack> base_stacks;
	RNA_BaseStackList const & base_stack_list = filtered_base_base_info.base_stack_list();
	for ( auto const & base_stack : base_stack_list ) {
		base_stacks.push_back( base_stack );
		base_stacks.push_back( base_stack.flipped() );
	}

	/////////////////////////////////////////
	// Straight up W/C base pair step?
	/////////////////////////////////////////
	for ( auto const & base_pair1  : base_pairs ) {
		if ( base_pair1.edge1() != WATSON_CRICK )  continue;
		if ( base_pair1.edge2() != WATSON_CRICK )  continue;
		if ( base_pair1.orientation() != ANTIPARALLEL )  continue;
		if ( !check_watson_crick_sequence( pose, base_pair1.res1(), base_pair1.res2() ) ) continue;
		for ( auto const & base_pair2  : base_pairs ) {
			if ( base_pair2.res1() != base_pair1.res1() + 1 ) continue;
			if ( base_pair1.res2() != base_pair2.res2() + 1 ) continue;
			if ( pose.pdb_info() && pose.pdb_info()->number( base_pair2.res1() ) != pose.pdb_info()->number( base_pair1.res1() ) + 1 ) continue;
			if ( pose.pdb_info() && pose.pdb_info()->chain(  base_pair2.res1() ) != pose.pdb_info()->chain(  base_pair1.res1() ) ) continue;
			if ( pose.pdb_info() && pose.pdb_info()->number( base_pair1.res2() ) != pose.pdb_info()->number( base_pair2.res2() ) + 1 ) continue;
			if ( pose.pdb_info() && pose.pdb_info()->chain(  base_pair1.res2() ) != pose.pdb_info()->chain(  base_pair2.res2() ) ) continue;
			if ( base_pair2.edge1() != WATSON_CRICK )  continue;
			if ( base_pair2.edge2() != WATSON_CRICK )  continue;
			if ( base_pair2.orientation() != ANTIPARALLEL )  continue;
			if ( !check_watson_crick_sequence( pose, base_pair2.res1(), base_pair2.res2() ) ) continue;
			// Following ends up being too stringent:
			//  if ( !check_stack( base_pair2.res1(), base_pair1.res1(), base_stacks ) ) continue;
			// if ( !check_stack( base_pair2.res2(), base_pair1.res2(), base_stacks ) ) continue;
			rna_motifs.push_back( RNA_Motif( WC_STACKED_PAIR,
				{ base_pair1.res1(), base_pair2.res1(), base_pair2.res2(), base_pair1.res2() } ) );
		}
	}

	/////////////////////////
	//  U-turns
	//
	//  2
	//  N
	// /  R 3
	// 1 U .. p
	//
	// See, e.g.,
	//
	//   Moore, Annual Review Biochem., 1999.
	//
	/////////////////////////
	for ( Size n = 1; n <= pose.size(); n++ ) {
		utility::vector1< Size > u_turn;
		if ( n > pose.size() - 3 ) continue;
		if ( !check_rna_loop( pose, n, 4 ) ) continue;
		if ( potential.get_base_backbone( pose.residue( n ), pose.residue( n + 3 ), 1 /* O1P */) >= 0.0 &&
				potential.get_base_backbone( pose.residue( n ), pose.residue( n + 3 ), 2 /* O2P */) >= 0.0 ) continue;
		if ( potential.get_base_backbone( pose.residue( n+2 ), pose.residue( n ), 6 /* O2' */) >= 0.0 )  continue;
		// should we also enforce stacking n+1, n+2?
		for ( Size i = n; i <= n+2; i++ ) u_turn.push_back( i );
		rna_motifs.push_back( RNA_Motif( U_TURN, u_turn ) );
	}


	/////////////////////////////////////////
	//  UA-handles
	//
	// 2 U o-@ A 3
	//   |   N <-- [could be lots of N's]
	// 1 G -o- C 4
	//
	/////////////////////////////////////////
	for ( Size n = 2; n <= pose.size(); n++ ) {
		for ( auto const & base_pair : base_pairs ) {
			// Look for the U-A Watson-Crick/Hoogsteen pair
			if ( base_pair.res1() == n && base_pair.edge1() == WATSON_CRICK && base_pair.edge2() == HOOGSTEEN ) {
				utility::vector1< Size > ua_handle;
				// Look for W/C pair stacked immediately previous, connected to U
				if ( !check_stack( n-1, n, base_stacks ) ) continue;
				for ( auto const & base_pair_prev : base_pairs ) {
					if ( base_pair_prev.res1() == n-1 && base_pair_prev.edge1() == WATSON_CRICK && base_pair_prev.edge2() == WATSON_CRICK ) {
						if ( base_pair_prev.res2() == base_pair.res2()+1 ) continue;
						if ( !check_stack( base_pair.res2(), base_pair_prev.res2(), base_stacks ) ) continue;
						ua_handle.push_back( base_pair_prev.res1() );
						ua_handle.push_back( base_pair.res1() /* n */);
						ua_handle.push_back( base_pair.res2() );
						ua_handle.push_back( base_pair_prev.res2() );
						rna_motifs.push_back( RNA_Motif( UA_HANDLE, ua_handle ) );
						break;
					}
				}
			}
		}
	}

	/////////////////////////////////////////
	//  T-loop
	//
	//  4
	//  N
	// /  R 5
	// 3 U .. p
	//   |  \         //
	// 2 U o-@ A 6
	//   |   N <-- [could be lots of N's]
	// 1 G -o- C 7
	//
	/////////////////////////////////////////
	for ( auto const & u_turn : rna_motifs.get_motifs( U_TURN ) ) {
		for ( auto const & ua_handle : rna_motifs.get_motifs( UA_HANDLE ) ) {
			// right on top:
			if ( u_turn[1] != ua_handle[2]+1 ) continue;
			if ( u_turn[3] != ua_handle[3]-1 ) continue;
			// check for stacking
			if ( !check_stack( ua_handle[2], u_turn[1], base_stacks ) ) continue;
			// check for pocket (lack of stacking!) that would allow intercalation
			// (this should probably be a bonus feature, not obligate)
			if ( check_stack(  ua_handle[3], u_turn[3], base_stacks ) ) continue;
			utility::vector1< Size > t_loop = make_vector1( ua_handle[1],ua_handle[2],u_turn[1],u_turn[2],u_turn[3],
				ua_handle[3],ua_handle[4]);
			rna_motifs.push_back( RNA_Motif( T_LOOP, t_loop ) );
		}
	}

	/////////////////////////////////////////
	//  Intercalated T-loops
	//
	//  4
	//  N
	// /  R 5
	// 3 U X 8 <--- intercalated base from *outside* T-loop
	//   |  \         //
	// 2 U o-@ A 6
	//   |   N <-- [could be lots of N's]
	// 1 G -o- C 7
	//
	/////////////////////////////////////////
	for ( auto const & t_loop : rna_motifs.get_motifs( T_LOOP ) ) {
		for ( auto const & base_stack : base_stacks ) {
			if ( base_stack.res1() == t_loop[ 5 ] ) {
				Size const & intercalator = base_stack.res2();
				for ( auto const & base_stack2 : base_stacks ) {
					if ( base_stack2.res1() == intercalator &&
							base_stack2.res2() == t_loop[ 6 ] ) {
						utility::vector1< Size > intercalated_t_loop = t_loop.residues();
						intercalated_t_loop.push_back( intercalator );
						rna_motifs.push_back( RNA_Motif( INTERCALATED_T_LOOP, intercalated_t_loop ) );
						break;
					}
				}
			}
		}
	}


	/////////////////////////////////////////
	// GNRA tetraloop
	//
	//  2
	//  N
	// /  R 3
	// 1 G .. p A 4
	//   G -o- C
	//
	// Following generalizes to GNRA with bulges, e.g.,
	//  pentaloop in SARS virus, GAGUA. Bulge must be
	//  after U-turn (4) and before WC (6).
	/////////////////////////////////////////
	for ( auto const & u_turn : rna_motifs.get_motifs( U_TURN ) ) {
		for ( auto const & base_pair : base_pairs ) {
			if ( base_pair.edge1() != WATSON_CRICK ) continue;
			if ( base_pair.edge2() != WATSON_CRICK ) continue;
			if ( base_pair.res1() + 1 != u_turn[ 1 ] ) continue;
			if ( !check_rna_loop( pose, base_pair.res1(), base_pair.res2() - base_pair.res1() + 1 ) ) continue;
			if ( !check_stack( base_pair.res1(), u_turn[1], base_stacks ) ) continue;
			if ( u_turn[3] >= base_pair.res2() ) continue;
			// Now look for 'A' in GNRA:
			for ( Size i = u_turn[3]+1; i < base_pair.res2(); i++ ) {
				if ( !check_stack( u_turn[3],  i, base_stacks ) ) continue;
				if ( !check_stack( i, base_pair.res2(), base_stacks ) ) continue;
				rna_motifs.push_back( RNA_Motif( GNRA_TETRALOOP,
					{ u_turn[1], u_turn[2], u_turn[3], i } ) );
				break;
			}
		}
	}

	/////////////////////////////////////////
	// A-minor motif
	//
	//   1 A o-@ 2 C -o- G 3
	//  A    G -o- C
	//
	//  Note that 1 & 2 are point of WC-Sugar contact.
	// The other A may be 5' or 3' of the A
	// The stacked pair may be 5' or 3' of the other stacked pair.
	//
	// Could also set up Type 0, I, II, III ?
	//  [note this supercedes 'ribose zipper']
	//
	//  An alternative would be to use
	//   Grabow et al., 2013 WIRES nomenclature --
	//   A-planar & A-twisted.
	//  Note that they force stereotyped sugar-sugar base pair as anchor:
	//
	//   A > G -o- C
	//   A   X -o- X
	//   5'  5' 3'
	//   Not implemented yet.
	/////////////////////////////////////////
	for ( Size i = 1; i < pose.size(); i++ ) {
		if ( pose.residue( i ).aa() != na_rad &&  pose.residue( i ).na_analogue() != na_rad ) continue;
		if ( pose.residue( i+1 ).aa() != na_rad && pose.residue( i+1 ).na_analogue() != na_rad ) continue;
		if ( !check_rna_loop( pose, i, 2 ) ) continue;
		if ( !check_stack( i, i+1, base_stacks ) ) continue;
		bool found_a_minor( false );
		for ( Size a_offset = 0; a_offset <= 1; a_offset++ ) {
			for ( auto const & base_pair : base_pairs ) {
				if ( base_pair.res1() != i + a_offset ) continue;
				// if ( base_pair.edge1() != WATSON_CRICK ) continue; // should be SUGAR in Grabow.
				if ( base_pair.edge2() != SUGAR ) continue;
				for ( auto const & stacked_pair : rna_motifs.get_motifs( WC_STACKED_PAIR ) ) {
					for ( Size stacked_pair_offset = 0; stacked_pair_offset <= 1; stacked_pair_offset++ ) {
						if ( stacked_pair[ 1+stacked_pair_offset ] == base_pair.res2() ) {
							// utility::vector1< Size > a_minor_res = a_offset ? make_vector1(i+1,i) : make_vector1(i+1,i);
							// if ( stacked_pair_offset ) {
							//  a_minor_res.append( {stacked_pair[2],stacked_pair[1],stacked_pair[4],stacked_pair[3]} );
							// } else {
							//  a_minor_res.append( stacked_pair.residues() );
							// }
							utility::vector1< Size > a_minor_res = { i + a_offset };
							if ( stacked_pair_offset ) {
								a_minor_res.append( {stacked_pair[2],stacked_pair[3]} );
							} else {
								a_minor_res.append( {stacked_pair[1],stacked_pair[4]} );
							}
							rna_motifs.push_back( RNA_Motif( A_MINOR, a_minor_res ) );
							found_a_minor = true; break;
						}
					}
					if ( found_a_minor ) break;
				} // stacked_pair
				if ( found_a_minor ) break;
			} // base_pair
			if ( found_a_minor ) {
				if ( a_offset ) i += 1; // already covered i,i+1
				break;
			}
		} // a_offset
	}

	/////////////////////////////////////////
	// Platform
	//   1  2
	//   A <-@ A
	//   U -o- G 3
	//   4
	//
	//  1 -> 2 -> 3 should be sequence contiguous!
	//
	// See Cate et al., P4-P6 RNA
	//
	/////////////////////////////////////////
	// "platform"
	for ( auto const & base_pair : base_pairs ) {
		// sequence adjacent!
		if ( base_pair.res1()+1 != base_pair.res2() /* canonical */ &&
				base_pair.res1()+2 != base_pair.res2() /* C7.2 */ ) continue;
		if ( !check_rna_loop( pose, base_pair.res2(), base_pair.res2() - base_pair.res1() + 1 ) ) continue;
		if ( base_pair.edge1() != SUGAR ) continue;
		if ( base_pair.edge2() != HOOGSTEEN /* canonical */ &&  base_pair.edge2() != WATSON_CRICK /* C7.2 */ ) continue;
		for ( auto const & base_pair2 : base_pairs ) {
			if ( base_pair2.res1() != base_pair.res2() + 1 ) continue; // sequence adjacent
			if ( !check_rna_loop( pose, base_pair.res2(), 2 ) ) continue;
			if ( !check_stack( base_pair.res1(), base_pair2.res2(), base_stacks ) ) continue;
			if ( !check_stack( base_pair.res2(), base_pair2.res1(), base_stacks ) ) continue;
			// perhaps could relax this -- then we could include bulged-G motif
			// if ( base_pair2.edge1() != WATSON_CRICK ) continue;
			// if ( base_pair2.edge2() != WATSON_CRICK ) continue;
			rna_motifs.push_back( RNA_Motif( PLATFORM,
				{ base_pair.res1(), base_pair.res2(),
				base_pair2.res1(), base_pair2.res2() } ) );
		}
	}


	/////////////////////////////////////////
	// Tetraloop receptor
	//
	//  classic 11-nt:
	// 1 C -o- G 10 <-- stacked pair
	// 2 C -o- G 9 -[U] <-- bulge is optional.
	// 3 U o-@ A 8 /
	// 4 A  <-@  A 5   <-- PLATFORM! Note strand crossover
	// 7 U  -o-  G 6
	//
	//  R(1):
	// 1 C -o- G 10 <-- stacked pair
	// 2 C -o- G 9
	// 3 U -o- U 8  <-- does not have to be UA handle.
	// 4 G  <-@  U 5   <-- PLATFORM! Note strand crossover
	// 7 A  -o-  G 6
	//
	// Based on similarities across
	// models of C7.2, C7.10, R1 &
	// classic 11-nt tetraloop receptors
	//
	/////////////////////////////////////////
	for ( auto const & platform : rna_motifs.get_motifs( PLATFORM ) ) {
		for ( auto const & base_pair : base_pairs ) {
			if ( base_pair.res1() == platform[4] ) continue;
			if ( base_pair.res2() == platform[3] ) continue;

			// sequence adjacent on 5' strand
			if ( base_pair.res1() >= platform[1] ) continue;
			if ( platform[1] - base_pair.res1() > 3 ) continue;

			// sequence adjacent on 3' strand
			if ( base_pair.res2() <= platform[4] ) continue;
			if ( base_pair.res2() - platform[4] > 3 ) continue;

			if ( !check_stack( base_pair.res1(), platform[1], base_stacks ) ) continue;

			for ( auto const & stacked_pair : rna_motifs.get_motifs( WC_STACKED_PAIR ) ) {

				// sequence adjacent on 5' strand
				if ( stacked_pair[2] >= base_pair.res1() ) continue;
				if ( base_pair.res1() - stacked_pair[2] > 3 ) continue;

				// sequence adjacent on 3' strand
				if ( stacked_pair[3] <= base_pair.res2() ) continue;
				if ( stacked_pair[3] - base_pair.res2() > 3 ) continue;

				if ( !check_stack( stacked_pair[2], base_pair.res1(), base_stacks ) &&
						!check_stack( stacked_pair[3], base_pair.res2(), base_stacks ) ) continue;
				rna_motifs.push_back( RNA_Motif( TL_RECEPTOR,
					{ stacked_pair[1], stacked_pair[2], base_pair.res1(), platform[1], platform[2], platform[3], platform[4], base_pair.res2(), stacked_pair[3], stacked_pair[4] } ) );
			}
		}
	}



	///////////////////////////////////////////
	// Tetraloop/tetraloop-receptor (docked!)
	//
	// 1 C -o- G 10  _ 4A  G 1
	// 2 C -o- G 9 _/  3A /
	// 3 U o-@ A 8  2A
	// 4 A  <-@  A 5
	// 7 U  -o-  G 6
	//
	///////////////////////////////////////////
	for ( auto const & tetraloop : rna_motifs.get_motifs( GNRA_TETRALOOP ) ) {
		for ( auto const & receptor : rna_motifs.get_motifs( TL_RECEPTOR ) ) {
			if ( !check_stack( tetraloop[2], receptor[5], base_stacks ) ) continue;
			for ( auto const & a_minor : rna_motifs.get_motifs( A_MINOR ) ) {
				if ( a_minor[1] != tetraloop[ 4 ] ) continue;
				if ( a_minor[2] != receptor[ 9 ] ) continue;
				utility::vector1< Size > res( tetraloop.residues() );
				res.append( receptor.residues() );
				rna_motifs.push_back( RNA_Motif( TETRALOOP_TL_RECEPTOR, res ) );
			}
		}
	}


	/////////////////////////////////////////////////////////////////////////////////
	// loopE-submotif
	//
	// 3  A @-o U 4
	// 2  G <-@ A 5
	// 1  G -o- C 6
	//
	//  occurs in loopE, SRL, SRP -- note that
	//  sequence can shift in, e.g., SRP
	//
	// See, e.g., Leontis & Westhof, Comp Funct Genom 2002; 3: 518–524.
	//  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2448414/pdf/CFG-03-518.pdf
	//
	// Specific case of "cross-strand purine stack" from Moore, Ann. Rev. Biochem. 1999.
	// Also called GA_UA motif [see, e.g., Grabow et al., WIREs RNA 2013]
	//
	/////////////////////////////////////////////////////////////////////////////////
	for ( auto const & base_pair1 : base_pairs ) {
		if ( base_pair1.edge1() != WATSON_CRICK ) continue;
		if ( base_pair1.edge2() != WATSON_CRICK ) continue;
		if ( !check_rna_loop( pose, base_pair1.res1(),   3 ) ) continue;
		if ( !check_rna_loop( pose, base_pair1.res2()-2, 3 ) ) continue;
		for ( auto const & base_pair2 : base_pairs ) {
			// direct sequence adjacent
			if ( base_pair2.res1() != base_pair1.res1()+1 ) continue;
			if ( base_pair1.res2() != base_pair2.res2()+1 ) continue;
			// 'sheared' base pair on top of W/C pair.
			if ( base_pair2.edge1() != SUGAR &&
					potential.get_base_backbone( pose.residue( base_pair2.res2() ),
					pose.residue( base_pair2.res1() ), 6 /* O2' */) >= 0.0  /* for 1lnt edge case*/ )  continue;
			if ( base_pair2.edge2() != HOOGSTEEN ) continue;
			// directly stacking of 'sheared' base pair on W/C
			if ( !check_stack( base_pair1.res1(), base_pair2.res1(), base_stacks ) ) continue;
			if ( !check_stack( base_pair1.res2(), base_pair2.res2(), base_stacks ) ) continue;
			for ( auto const & base_pair3 : base_pairs ) {
				if ( base_pair3.edge1() != HOOGSTEEN ) continue;
				if ( base_pair3.edge2() != WATSON_CRICK ) continue;
				// direct sequence adjacent
				if ( base_pair3.res1() != base_pair2.res1()+1 ) continue;
				if ( base_pair2.res2() != base_pair3.res2()+1 ) continue;
				// actually, these two will not stack if in bulged-G
				// if ( !check_stack( base_pair2.res1(), base_pair3.res1(), base_stacks ) ) continue;
				// this is the "cross-strand" purine stack:
				if ( !check_stack( base_pair2.res2(), base_pair3.res1(), base_stacks ) ) continue;
				rna_motifs.push_back( RNA_Motif( LOOP_E_SUBMOTIF,
					{ base_pair1.res1(), base_pair2.res1(), base_pair3.res1(),
					base_pair3.res2(), base_pair2.res2(), base_pair1.res2() } ) );
			}
		}
	}

	/////////////////////////////////////////////////////
	// bulged-G motif
	//     5  4
	// 3  A @-o  U @-> G
	// 2  G <-@  A 6
	// 1  G -o-  C 7
	//
	//  bulged-G motif from, e.g., Sarcin-Ricin loop
	//
	// Note: should we consider the G4-U5 a 'platform'?
	// TODO: These motifs usually have an additional
	//   N and A 5' of G4 with a cool turn -- include here?
	//   Or define a more extended motif to include the S-turn?
	//
	/////////////////////////////////////////////////////
	for ( auto const & loop_e_submotif : rna_motifs.get_motifs( LOOP_E_SUBMOTIF ) ) {
		for ( auto const & base_pair : base_pairs ) {
			if ( base_pair.res2() == loop_e_submotif[ 4 ] ) {
				if ( base_pair.edge1() != SUGAR ) continue; /*bulged G*/
				if ( base_pair.edge2() != HOOGSTEEN ) continue;
				// "platform"
				if ( base_pair.res1()+1 == base_pair.res2() /* SRL*/ ||
						base_pair.res1()+2 == base_pair.res2() /* FMN aptamer, PDB 1FMN*/ ) {
					if ( !check_rna_loop( pose, base_pair.res1(), base_pair.res2() - base_pair.res1() ) ) continue;
					// check that bulged G makes contact with phosphate
					if ( potential.get_base_backbone( pose.residue(base_pair.res1()), pose.residue(loop_e_submotif[3]), 1 /* O1P */) >= 0.0 &&
							potential.get_base_backbone( pose.residue(base_pair.res1()), pose.residue(loop_e_submotif[3]), 2 /* O2P */) >= 0.0 ) continue;
					rna_motifs.push_back( RNA_Motif( BULGED_G,
						{ loop_e_submotif[1], loop_e_submotif[2], loop_e_submotif[3],
						base_pair.res1(),
						loop_e_submotif[4], loop_e_submotif[5], loop_e_submotif[6] } ) );
				}
			}
		}
	}


	/////////////////////////////////////////
	// Z-turn
	//  2
	//  N
	// /  C 3
	// 1 U <-o G 4
	//   G -o- C
	/////////////////////////////////////////

	// AMW: first implementing the classic Z-turn, UNCG.
	// Theoretically: could implement CCNG (note cWW CG bp; bulged U at 3' end)
	for ( auto const & base_pair : base_pairs ) {
		if ( base_pair.edge1() != WATSON_CRICK ) continue;
		if ( base_pair.edge2() != WATSON_CRICK ) continue;
		if ( base_pair.res1() != base_pair.res2() - 5 ) continue;
		if ( !check_rna_loop( pose, base_pair.res1(), base_pair.res2() - base_pair.res1() + 1 ) ) continue;
		if ( !check_stack( base_pair.res1(), base_pair.res1()+1, base_stacks ) ) continue;

		// Look for U-G pair
		for ( auto const & base_pair2 : base_pairs ) {
			if ( base_pair2.res1() != base_pair.res1() + 1 ) continue; // Should we accomodate bulges?
			if ( base_pair2.res2() != base_pair.res2() - 1 ) continue; // Should we accomodate bulges?
			if ( base_pair2.edge1() != SUGAR ) continue;
			if ( base_pair2.edge2() != WATSON_CRICK ) continue;

			// I'm struggling to define this orientation properly: 3U4M is, as desired, trans. But
			// the tetraloops in 1F7Y are not trans. Both 3U4M and 1F7Y are parallel.
			// Furthermore, in the CCNG variant this becomes a cWW bp, so maybe we shouldn't
			// be so anal about LW orientation?
			//if ( base_pair2.LW_orientation() != TRANS ) continue;

			// Stacked on prior bp
			if ( !check_stack( base_pair.res1(), base_pair2.res1(), base_stacks ) ) continue;
			if ( !check_stack( base_pair.res2(), base_pair2.res2(), base_stacks ) ) continue;

			// Must be U-G
			if ( pose.aa( base_pair2.res1() ) != na_ura && pose.residue_type( base_pair2.res1() ).na_analogue() != na_ura ) continue;
			if ( pose.aa( base_pair2.res2() ) != na_rgu && pose.residue_type( base_pair2.res2() ).na_analogue() != na_rgu ) continue;

			// Penultimate nt must be C
			if ( pose.aa( base_pair2.res2() - 1 ) != na_rcy && pose.residue_type( base_pair2.res2() - 1 ).na_analogue() != na_rcy ) continue;

			rna_motifs.push_back( RNA_Motif( Z_TURN, { base_pair2.res1(), base_pair2.res1() + 1, base_pair2.res2() - 1, base_pair2.res2() } ) );
			break;
		}
	}

	/////////////////////////////////////////
	// U-turn loop
	//  generalization of GNRA
	//  see Auffinger
	// AMW TODO: make GNRA loop come after this one,
	// so we can only look for GNRA within this class
	/////////////////////////////////////////

	/////////////////////////////////////////
	// G-quartet
	/////////////////////////////////////////

	/////////////////////////////////////////
	// G-quartet-stack
	/////////////////////////////////////////

	/////////////////////////////////////////
	// Tandem G/A (sheared)
	// a la Cruz/Westhof.
	// or cross-strand purine stack [Moore, 1999]
	// 1 G <-@ A 4
	// 2 A @-< G 3
	/////////////////////////////////////////

	// AMW: enforce uniqueness?
	for ( auto const & base_pair1 : base_pairs ) {
		if ( pose.aa( base_pair1.res1() ) != na_rgu && pose.residue_type( base_pair1.res1() ).na_analogue() != na_rgu ) continue;
		if ( pose.aa( base_pair1.res2() ) != na_rad && pose.residue_type( base_pair1.res2() ).na_analogue() != na_rad ) continue;

		if ( base_pair1.edge1() != SUGAR ) continue;
		if ( base_pair1.edge2() != HOOGSTEEN ) continue;
		if ( base_pair1.LW_orientation() != TRANS ) continue;

		for ( auto const & base_pair2 : base_pairs ) {

			if ( pose.aa( base_pair2.res1() ) != na_rgu && pose.residue_type( base_pair2.res1() ).na_analogue() != na_rgu ) continue;
			if ( pose.aa( base_pair2.res2() ) != na_rad && pose.residue_type( base_pair2.res2() ).na_analogue() != na_rad ) continue;

			if ( base_pair2.edge1() != SUGAR ) continue;
			if ( base_pair2.edge2() != HOOGSTEEN ) continue;
			if ( base_pair2.LW_orientation() != TRANS ) continue;

			// This base pair has to have its A adjacent to the prior G, and vice-versa
			if ( !( base_pair2.res2() == base_pair1.res1() + 1 && base_pair1.res1() == base_pair2.res2() - 1 ) ) continue;

			// Avoid duplication: bp1 r1 < bp2 r1
			if ( base_pair1.res1() > base_pair2.res1() ) continue;

			// AA stack; GG stack
			if ( !check_stack( base_pair2.res1(), base_pair1.res1(), base_stacks ) ) continue;
			if ( !check_stack( base_pair2.res2(), base_pair1.res2(), base_stacks ) ) continue;

			rna_motifs.push_back( RNA_Motif( TANDEM_GA_SHEARED, { base_pair1.res1(), base_pair2.res2(), base_pair2.res1(), base_pair1.res2() } ) );
			break;
		}
	}

	/////////////////////////////////////////
	// GA-minor
	// important for kink-turn -- tried to keep
	// this general but when necessary I ensured
	// this would capture the one in a KT
	//
	//    [s] N 5
	// 1 A -<- g 4
	// 2 G -o- C 3
	// here [s] represents base of A stacking on ribose of N
	// 23 is truly any WC pair.
	// Note: here we permit A to be any nt as long as its sugar
	// edge makes a sufficient interaction with 4's O2'
	/////////////////////////////////////////

	for ( auto const & base_pair1 : base_pairs ) {
		if ( base_pair1.edge1() != WATSON_CRICK ) continue;
		if ( base_pair1.edge2() != WATSON_CRICK ) continue;

		if ( base_pair1.res1() == pose.size() ) continue;
		if ( base_pair1.res2() == 1 ) continue;

		// Arbitrarily to disambiguate: bp1.res1 is next to the non-A, bp1.res2 is next to the A.

		// Two next residues to consider; they are not-quite base paired to each other.
		for ( Size ii = base_pair1.res1() + 1;
				ii <= ( base_pair1.res1() + 4 > pose.size() - 1 ? pose.size() - 1 : base_pair1.res1() + 4 );
				++ii ) {
			if ( !check_stack( ii, base_pair1.res1(), base_stacks ) ) continue;

			for ( Size jj = ( base_pair1.res2() - 4 < 1 ? 1 : base_pair1.res2() - 4 );
					jj <= base_pair1.res2() - 1;
					++jj ) {

				// The sugar edges of ii and jj are related iff there is a strong base-O2' interaction. Hooray!
				if ( potential.get_base_backbone( pose.residue( jj ), pose.residue( ii ), 6 /*O2'*/) > -1.0 ) continue;

				// SECOND: stack b/n jj and ii + 1 ribose.
				// Can't use base-backbone potential for this -- it won't pick up on anything!
				// Instead, measure base centroid to sugar distance.
				auto jj_centroid = get_rna_base_centroid( pose.residue( jj ), false );
				// AMW TODO: convert all these atom names to atom indices stemming from the RNA_info object
				// This will provide generality for some chemically modified residues. (though mods to these
				// sugar atoms are pretty rare)
				if ( jj_centroid.distance( pose.residue( ii + 1 ).xyz( "O4'" ) ) < 4.0
						|| jj_centroid.distance( pose.residue( ii + 1 ).xyz( "C4'" ) ) < 4.0
						|| jj_centroid.distance( pose.residue( ii + 1 ).xyz( "C3'" ) ) < 4.0
						|| jj_centroid.distance( pose.residue( ii + 1 ).xyz( "C2'" ) ) < 4.0
						|| jj_centroid.distance( pose.residue( ii + 1 ).xyz( "C1'" ) ) < 4.0 ) {
					rna_motifs.push_back( RNA_Motif( GA_MINOR, { jj, base_pair1.res2(), base_pair1.res1(), ii, ii + 1 } ) );
				}
			}
		}
	}


	/////////////////////////////////////////
	// Kink-turn
	// contains tandem G/A
	// Note: one generalization permits one of the GAs to be cWW
	/////////////////////////////////////////
	//for ( auto const & tandem_ga_sheared : rna_motifs.get_motifs( TANDEM_GA_SHEARED ) ) {
	// Look for nnn 5' of the g on the strand where g comes first.

	//}


	/////////////////////////////////////////
	// Tandem W/C
	//  A -o- A
	//  G -o- A
	// should this include 'normal' base pair
	//  steps?
	// may be needed to counter-balance
	//  sheared G/A
	// Happens in 1FMN (FMN aptamer)
	/////////////////////////////////////////

	for ( auto const & base_pair1 : base_pairs ) {
		if ( pose.aa( base_pair1.res1() ) != na_rgu && pose.residue_type( base_pair1.res1() ).na_analogue() != na_rgu ) continue;
		if ( pose.aa( base_pair1.res2() ) != na_rad && pose.residue_type( base_pair1.res2() ).na_analogue() != na_rad ) continue;

		if ( base_pair1.edge1() != WATSON_CRICK ) continue;
		if ( base_pair1.edge2() != WATSON_CRICK ) continue;
		if ( base_pair1.LW_orientation() != TRANS ) continue;

		for ( auto const & base_pair2 : base_pairs ) {

			if ( pose.aa( base_pair2.res1() ) != na_rgu && pose.residue_type( base_pair2.res1() ).na_analogue() != na_rgu ) continue;
			if ( pose.aa( base_pair2.res2() ) != na_rad && pose.residue_type( base_pair2.res2() ).na_analogue() != na_rad ) continue;

			if ( base_pair2.edge1() != WATSON_CRICK ) continue;
			if ( base_pair2.edge2() != WATSON_CRICK ) continue;
			if ( base_pair2.LW_orientation() != TRANS ) continue;

			// This base pair has to have its A adjacent to the prior G, and vice-versa
			if ( !( base_pair2.res2() == base_pair1.res1() + 1 && base_pair1.res1() == base_pair2.res2() - 1 ) ) continue;

			// Avoid duplication: bp1 r1 < bp2 r1
			if ( base_pair1.res1() > base_pair2.res1() ) continue;

			// AA stack; GG stack
			if ( !check_stack( base_pair2.res1(), base_pair1.res1(), base_stacks ) ) continue;
			if ( !check_stack( base_pair2.res2(), base_pair1.res2(), base_stacks ) ) continue;

			rna_motifs.push_back( RNA_Motif( TANDEM_GA_WATSON_CRICK, { base_pair1.res1(), base_pair2.res2(), base_pair2.res1(), base_pair1.res2() } ) );
			break;
		}
	}

	/////////////////////////////////////////
	// what is P4-P6 j5/5a? --
	//  W/C stacked pair with minor groove
	//  docked into non-W/C stacked-pair
	// May be describable as an 'off-label'
	// UA_h turn, since it has a UA_handle
	//  made with a CC instead of a UA.
	/////////////////////////////////////////

	/////////////////////////////////////////
	// G-ribo
	/////////////////////////////////////////

	/////////////////////////////////////////
	// Double T-loop
	//  a la T-box, RNase P
	//
	// Note: this is sometimes so separated from the
	// WC pair, and sometimes not stacked on them (see 3q1q)
	// so while maybe you COULD find one clean loop w/ WC and
	// all, just don't bother...
	//
	//  4
	//  N ----------- N 7
	// / /-- 12 A @-o N 8
	// 3 N  R 5  \    |
	//   |   \ 11 R   N 9
	// 2 N o-@ A 6 ---/ /
	// 1 N ----------- N
	//                 10
	//
	/////////////////////////////////////////

	for ( auto const & u_turn1 : rna_motifs.get_motifs( U_TURN ) ) {
		// First search for U-turns. Note that this RELIES on the definition of
		// a U turn continuing to permit starting on non-U.
		for ( auto const & base_pair1 : base_pairs ) {
			// convert to pdbinfo numbering? AMW TODO
			if ( base_pair1.res1() != u_turn1[1] - 1 ) continue;
			if ( base_pair1.res2() != u_turn1[3] + 1 ) continue;
			if ( base_pair1.edge1() != WATSON_CRICK ) continue;
			if ( base_pair1.edge2() != HOOGSTEEN ) continue;
			if ( base_pair1.LW_orientation() != TRANS ) continue;

			for ( auto const & u_turn2 : rna_motifs.get_motifs( U_TURN ) ) {
				for ( auto const & base_pair2 : base_pairs ) {
					if ( base_pair2.res1() != u_turn2[1] - 1 ) continue;
					if ( base_pair2.res2() != u_turn2[3] + 1 ) continue;
					if ( base_pair2.edge1() != WATSON_CRICK ) continue;
					if ( base_pair2.edge2() != HOOGSTEEN ) continue;
					if ( base_pair2.LW_orientation() != TRANS ) continue;

					// We have now identified bp1 and bp2, which are the UA-handle analogues.
					// Now we loop through BPs again, hoping to find:
					// u_turn1[2] to bp2.res1() - 1
					// bp1.res2() to u_turn2[1]
					// u_turn2[2] to bp1.res1() - 1
					// bp2.res2() to u_turn1[1]
					//
					// AMW TODO: figure out a way to loop only once, tracking for each
					// of these that are found.
					for ( auto const & base_pair3 : base_pairs ) {
						if ( base_pair3.res1() != u_turn1[2] ) continue;
						if ( base_pair3.res2() != base_pair2.res1() - 1 ) continue;
						if ( base_pair3.edge1() != HOOGSTEEN ) continue;
						if ( base_pair3.edge2() != WATSON_CRICK ) continue;
						if ( base_pair3.LW_orientation() != CIS ) continue;

						for ( auto const & base_pair4 : base_pairs ) {
							if ( base_pair4.res1() != base_pair1.res2() ) continue;
							if ( base_pair4.res2() != u_turn2[1] ) continue;
							if ( base_pair4.edge1() != WATSON_CRICK ) continue;
							if ( base_pair4.edge2() != SUGAR ) continue;
							if ( base_pair4.LW_orientation() != TRANS  ) continue;

							for ( auto const & base_pair5 : base_pairs ) {
								if ( base_pair5.res1() != u_turn2[2] ) continue;
								if ( base_pair5.res2() != base_pair1.res1() - 1 ) continue;
								if ( base_pair5.edge1() != HOOGSTEEN ) continue;
								if ( base_pair5.edge2() != WATSON_CRICK ) continue;
								if ( base_pair5.LW_orientation() != CIS  ) continue;

								// Sometimes this base pair won't form. This is despite two quite evident
								// H-bonds (as in, for example, 3q1q).
								//for ( auto const & base_pair6 : base_pairs ) {
								// if ( base_pair6.res1() != base_pair2.res2() ) continue;
								// if ( base_pair6.res2() != u_turn1[1] ) continue;
								// if ( base_pair6.edge1() != WATSON_CRICK ) continue;
								// if ( base_pair6.edge2() != SUGAR ) continue;
								// if ( base_pair6.LW_orientation() != TRANS  ) continue;

								rna_motifs.push_back( RNA_Motif( DOUBLE_T_LOOP, {
									base_pair1.res1() - 1, base_pair1.res1(), u_turn1[1], u_turn1[2], u_turn1[3], base_pair1.res2(),
									base_pair2.res1() - 1, base_pair2.res1(), u_turn2[1], u_turn2[2], u_turn2[3], base_pair2.res2(),
									} ) );
								//}
							}
						}
					}
				}
			}

		}
	}


	/////////////////////////////////////////
	// A-minor junctions
	//  a la Jaeger
	/////////////////////////////////////////

	/////////////////////////////////////////
	// T-loop/PK
	//  a la Jaeger
	/////////////////////////////////////////

	return rna_motifs;
}

/////////////////////////////////////////////////////////////////
/// @details used in RiboDraw (which takes output from rna_motif app)
void
output_rna_motifs(
	std::ostream & out,
	core::pose::Pose const & pose,
	RNA_Motifs const & rna_motifs,
	bool const output_WC_stacked_pair/* = false */ )
{
	using utility::vector1;
	for ( auto const & motif : rna_motifs ) {

		if ( !output_WC_stacked_pair ) {
			if ( rna_motif_bonus.count( motif.type() ) == 0 ) continue;
			if ( rna_motif_bonus.find( motif.type() )->second >= 0.0  ) continue;
		}

		vector1< int > res;
		vector1< char > chain;
		vector1< std::string > segid;
		for ( auto const & m : motif ) {
			res.push_back( pose.pdb_info()->number( m ) );
			chain.push_back( pose.pdb_info()->chain( m ) );
			segid.push_back( pose.pdb_info()->segmentID( m ) );
		}
		out << to_string( motif.type() ) << " " << make_tag_with_dashes( res, chain, segid ) << std::endl;
	}
}

// @brief super-simple helper function for PyMOL commands.
// @details for speed, could also define PyMOL boolean property and then color things at end based on property.
void
output_motifs_to_pymol( std::ostream & out,
	pose::Pose const & pose,
	RNA_Motifs const & rna_motifs ) {

	std::string tag( tag_from_pose( pose ) );
	tag = utility::replace_in( tag, ".pdb", "" );
	for ( auto const & motif : rna_motifs ) {
		for ( auto const & res : motif ) {
			out << "color " << motif_color.find( motif.type() )->second << ", " << tag << " and chain " << pose.pdb_info()->chain( res ) << " and resi " << pose.pdb_info()->number( res ) << std::endl;

		}
	}
}

} //rna
} //scoring
} //core
