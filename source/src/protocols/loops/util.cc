// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Oliver Lange
/// @author Mike Tyka

// Unit Headers
#include <protocols/loops/util.hh>

// Package Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Project Headers
#include <numeric/xyz.functions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/loops/loops_main.hh> //for getting ss from dssp

#include <core/fragment/SecondaryStructure.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>

#include <core/scoring/func/HarmonicFunc.hh>

#include <core/conformation/util.hh> //idealize
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <basic/Tracer.hh>

//numeric headers
#include <numeric/conversions.hh>

//// C++ headers
#include <list>

//Auto Headers
#include <utility/vector1.hh>

//Auto Headers
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS

static thread_local basic::Tracer TR( "protocols.loops.util" );

namespace protocols {
namespace loops {

/// TODO(someone) so bad!
using namespace core;
using namespace pose;
using namespace kinematics;

// Numeric constants
static const double EXT_PHI = -150;
static const double EXT_PSI = +150;
static const double EXT_OMG =  180;

void
fix_with_coord_cst( Loops const& rigid, core::pose::Pose& pose, bool bCstAllAtom, utility::vector1< core::Real > &weights ) {
	bool const bReadWeights = ( weights.size() >= pose.total_residue() );
	if ( !bReadWeights ) {
		weights.resize( pose.total_residue() );
	}

	for ( Loops::const_iterator it = rigid.begin(), eit = rigid.end();
			it!=eit; ++it ) {
		for ( Size pos = it->start(); pos <= it->stop(); ++pos ) {
			Size const seq_dist( std::min( (int) pos - it->start(), (int) it->stop() - pos ) + 1);
			Real coord_sdev;
			if ( bReadWeights ) {
				coord_sdev = weights[ pos ];
			} else {
				coord_sdev = ( 1.0/seq_dist ); //or something else?
				weights[ pos ] = coord_sdev;
			}
			core::scoring::func::FuncOP fx( new scoring::func::HarmonicFunc( 0.0, coord_sdev ) );
			conformation::Residue const & rsd( pose.residue( pos ) );
			if ( bCstAllAtom ) {
				for ( Size ii = 1; ii<= rsd.natoms(); ++ii ) {
					pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new scoring::constraints::CoordinateConstraint(
						id::AtomID( ii, pos),
						id::AtomID( 1, pos ) /*this is completely ignored! */,
						rsd.xyz( ii ),
						fx
						) ) ) );
				}
			} else {
				id::AtomID atomID( pose.residue_type(pos).atom_index("CA"), pos );
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new scoring::constraints::CoordinateConstraint(
					atomID,
					id::AtomID( 1, pos ) /*this is completely ignored! */,
					rsd.xyz( atomID.atomno() ),
					fx
					) ) ) );
			}
		}
	}
}

/// @brief get frags that are fully within the Loop --- shorten(=true/false) frags that are close to the end of loops.
void select_loop_frags(
	loops::Loops const& loops,
	core::fragment::FragSet& source,
	core::fragment::FragSet& loop_frags,
	Size shorten
) {
	using namespace core::fragment;
	//assuming backbone degrees of freedom are wanted by the loops.
	// Jumps are filtered out
	kinematics::MoveMap movemap;
	movemap.set_bb( false );
	loops.switch_movemap( movemap, id::BB, true );


	InsertMap insert_map;
	InsertSize insert_size;
	source.generate_insert_map( movemap, insert_map, insert_size);

	{ //debug
		Size const total_insert = insert_map.size();
		TR.Trace << "size of insertmap: " << total_insert << " -- ";
		for ( Size i = 1; i<=total_insert; i++ ) TR.Trace << " " << insert_map[ i ];
		TR.Trace << "insert_size: \nRESIDUES: ";
		for ( Size i = 1; i<=insert_map[ total_insert ]; i++ ) TR.Trace << " " << RJ(3,i);
		TR.Trace <<"\nINSERT:   ";
		for ( Size i = 1; i<=insert_map[ total_insert ]; i++ ) TR.Trace << " " << RJ(3,insert_size[ i ]);
		TR.Trace << std::endl;
	}

	Size const total_insert = insert_map.size();

	for ( Size i = 1; i<=total_insert; i++ ) {
		Size const pos ( insert_map[ i ] );
		Size const size ( insert_size[ pos ] );
		FrameList copy_frames;
		source.frames( pos, copy_frames );
		for ( FrameList::iterator it = copy_frames.begin(), eit = copy_frames.end();
				it!=eit; ++it ) {
			TR.Trace << "add frame at pos " << pos << " " << (*it)->length() << " insert_size " << size << std::endl;
			if ( (*it)->length() == size ) loop_frags.add( *it );
			else if ( shorten && size > shorten ) loop_frags.add( (*it)->generate_sub_frame( size ) );
		}
	}
} //select_loop_frags

void safe_set_extended_torsions_and_idealize_loops(const protocols::loops::Loops& loops,
	core::pose::Pose* pose) {
	if ( loops.empty() ) {
		return;
	}

	set_extended_torsions_and_idealize_loops(*pose, loops);
}

void set_extended_torsions_and_idealize_loops(core::pose::Pose& pose, loops::Loops loops) {
	using core::Size;
	using core::conformation::Conformation;
	using protocols::loops::Loop;
	using protocols::loops::Loops;

	// if no loops, we want to have extended structure everywhere.
	// it is a by-value parameter -- as intenden the change is kept local
	if ( loops.empty() ) {
		loops.add_loop(1, pose.total_residue(), 0);
	}

	TR.Debug << "extend structure for " << loops << std::endl;

	Conformation& conf = pose.conformation();
	for ( Loops::const_iterator i = loops.begin(); i != loops.end(); ++i ) {
		const Loop& loop = *i;

		for ( Size j = loop.start(); j <= loop.stop(); ++j ) {
			core::conformation::idealize_position(j, conf);
			pose.set_phi(j, EXT_PHI);
			pose.set_psi(j, EXT_PSI);
			pose.set_omega(j, EXT_OMG);
		}
	}
}

void addScoresForLoopParts(
	core::pose::Pose & pose,
	loops::Loops loops,
	const core::scoring::ScoreFunction & scorefxn,
	core::pose::Pose & native_pose,
	core::Size nloops
) {

	using namespace core;

	Size nres = pose.total_residue();
	utility::vector1< core::Size > all_loop_list;
	for ( Size i = 1; i < nres; i ++ ) {
		if ( loops.is_loop_residue(i) ) all_loop_list.push_back( i );
	}
	scorefxn(pose);
	setPoseExtraScore( pose, "ScoreCore",  scorefxn.get_sub_score_exclude_res( pose, all_loop_list ) );

	Real score =  scorefxn(pose);

	for ( Size l = 1; l <= nloops; l++ ) {
		if ( l > loops.size() ) {
			setPoseExtraScore( pose, "ScoreLoopI" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			setPoseExtraScore( pose, "ScoreLoopC" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			setPoseExtraScore( pose, "ScoreLoopL" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			continue;
		}
		utility::vector1< core::Size > loop_list;
		utility::vector1< core::Size > non_loop_list;
		for ( Size i = 1; i < nres; i ++ ) {
			if ( ( i < loops[l].start() ) || ( i > loops[l].stop() ) ) {
				loop_list.push_back( i );
			} else {
				non_loop_list.push_back( i );
			}
		}

		Real loopscore = scorefxn.get_sub_score_exclude_res( pose, loop_list );
		Real nonloopscore = scorefxn.get_sub_score_exclude_res( pose, non_loop_list );
		setPoseExtraScore( pose, "ScoreLoopI" + ObjexxFCL::right_string_of(l,3,'0'), loopscore );
		setPoseExtraScore( pose, "ScoreLoopC" + ObjexxFCL::right_string_of(l,3,'0'), score - loopscore - nonloopscore );
		setPoseExtraScore( pose, "ScoreLoopL" + ObjexxFCL::right_string_of(l,3,'0'), score - nonloopscore );
	}

	// Work out RMS values too

	core::pose::Pose native_pose_super = native_pose;
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, native_pose_super, core::id::BOGUS_ATOM_ID );
	for ( core::Size ir=1; ir <= native_pose.total_residue(); ++ir ) {
		runtime_assert( ir <=  pose.total_residue() );
		runtime_assert( ir <=  native_pose_super.total_residue() );
		if ( ( !loops.is_loop_residue( ir ) ) && pose.residue(ir).is_protein() ) {
			id::AtomID const id1( native_pose_super.residue(ir).atom_index("CA"), ir );
			id::AtomID const id2( pose.residue(ir).atom_index("CA"), ir );
			atom_map.set(id1, id2);
		}
	}
	core::scoring::superimpose_pose( native_pose_super, pose, atom_map );

	int corelength;
	setPoseExtraScore( pose, "corerms", native_loop_core_CA_rmsd( native_pose, pose, loops, corelength ) );

	for ( Size l = 1; l <= nloops; l++ ) {
		if ( l > loops.size() ) {
			setPoseExtraScore( pose, "RMSLoop" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			continue;
		}
		Loops temploops;
		temploops.add_loop( loops[l] );
		setPoseExtraScore( pose, "RMSLoop" + ObjexxFCL::right_string_of(l,3,'0'),
			loops::loop_rmsd( native_pose_super, pose, temploops, true ) );
	}
}

loops::Loops compute_ss_regions(
	core::Real max_loop_frac,
	core::Size min_length,
	core::fragment::SecondaryStructure const & ss
) {

	using core::Size;

	Size start( 0 );
	Size last( 0 );
	Size max_gap( 2 );
	loops::Loops ss_regions;
	for ( Size pos = 1; pos <= ss.total_residue(); ++pos ) {
		if ( ss.loop_fraction( pos ) <= max_loop_frac ) {
			if ( !start ) {
				start = pos;
				last = pos - 1;
			}
			if ( last + max_gap < pos ) {
				if ( last - start >= min_length ) {
					ss_regions.add_loop( start, last );
				}
				start=0;
			}
			last = pos;
		}
	}
	return ss_regions;
}

core::scoring::ScoreFunctionOP get_cen_scorefxn() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string const weights( option[ OptionKeys::loops::cen_weights ]() ),
		patch( option[ OptionKeys::loops::cen_patch ]() );
	return scoring::ScoreFunctionFactory::create_score_function(
		weights, patch
	);
}

core::scoring::ScoreFunctionOP get_fa_scorefxn() {
	return core::scoring::get_score_function();
}


void add_coordinate_constraints_to_pose( core::pose::Pose & pose, const core::pose::Pose &constraint_target_pose,  protocols::loops::Loops &exclude_regions ){
	using namespace conformation;
	using namespace core;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace id;
	using namespace scoring::constraints;

	core::Size nnonvrt_cst_target = constraint_target_pose.total_residue();
	core::Size nnonvrt_pose = pose.total_residue();

	while ( pose.residue( nnonvrt_pose ).aa() == core::chemical::aa_vrt ) { nnonvrt_pose--; }
	while ( constraint_target_pose.residue( nnonvrt_cst_target ).aa() == core::chemical::aa_vrt ) { nnonvrt_cst_target--; }

	protocols::loops::Loops coordconstraint_segments;
	coordconstraint_segments = exclude_regions.invert( nnonvrt_cst_target );

	//TR.Info << coordconstraint_segments << std::endl;

	if ( nnonvrt_pose != nnonvrt_cst_target ) {
		TR.Error << "ERROR coord constraint pose length mismatch with input pose: " << nnonvrt_cst_target << " vs. " << nnonvrt_pose << std::endl;
		utility_exit();
	}

	if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		pose.append_residue_by_jump
			( *ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ),
			pose.total_residue()/2 );
	}


	Size nres = pose.total_residue();
	Real const coord_sdev( 0.5 );
	core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, coord_sdev ) );

	for ( Size i = 1; i<= (Size)nres; ++i ) {
		if ( i==(Size)pose.fold_tree().root() ) continue;
		if ( coordconstraint_segments.is_loop_residue( i ) ) {
			Residue const & nat_i_rsd( pose.residue(i) );
			for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint(
					AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
					fx ) ) ) );
			}
		}
	}

}

LoopsOP
loops_from_string( std::string const & loop_str, core::pose::Pose const & pose ){
	utility::vector1< std::string > const loops_vec( utility::string_split( loop_str, ',' ) );

	// each loop should have the format loop_start:loop_end:cut
	// if cut is not set then it's taken to be 0. Residue numbering can follow the
	// pdb numbering
	LoopsOP loops_from_tag( new Loops() );
	BOOST_FOREACH ( std::string const residue_pair, loops_vec ) {
		utility::vector1< std::string > const residues( utility::string_split( residue_pair, ':' ) );
		if ( residues.size() != 2 && residues.size() != 3 ) {
			utility_exit_with_message(
				"To specify a loops string it must have the format "
				"\"loop_start:loop_end:cut[,loop_start:loop_end:cut...]\". "
				"If the cut is not set, then it is taken to be zero. "
				"The residue numbering can use the pdb numbering.");
		}
		core::Size const loop_start( core::pose::parse_resnum( residues[ 1 ], pose ) );
		core::Size const loop_stop( core::pose::parse_resnum( residues[ 2 ], pose ) );
		core::Size loop_cut( 0 );
		if ( residues.size() == 3 ) {
			loop_cut = core::pose::parse_resnum( residues[ 3 ], pose );
		}
		runtime_assert( loop_start <= loop_stop );
		runtime_assert( loop_start >= 1 );
		runtime_assert( loop_stop <= pose.total_residue() );
		loops_from_tag->add_loop( loop_start, loop_stop, loop_cut );
	}
	return( loops_from_tag );
}

void define_scorable_core_from_secondary_structure(
	core::fragment::SecondaryStructure const& ss_def,
	protocols::loops::Loops& scored_core )
{
	using namespace core;
	using namespace basic::options;
	// Size const max_loop_size( 3 );
	// Size const max_short_helix( 5 );
	Size const max_loop_size( option[ OptionKeys::evaluation::score_sscore_maxloop ]() );
	Size const max_short_helix( option[ OptionKeys::evaluation::score_sscore_short_helix ]() );

	//find residues that are part of a short helix -- less than or equal to 5 residues
	utility::vector1< bool > short_helix( ss_def.total_residue(), false );

	//selection of loop definitions...
	//these loops define regions that are scored. Add all loops that are 4 residues or longer.
	//subsequently add also helices that have fewer than 6 residues if they terminated a long loop (>=4)
	loops::Loops unscored_loops;

	for ( Size pos=1; pos <= ss_def.total_residue(); pos++ ) {

		//detect loops
		if ( ss_def.loop_fraction( pos ) > 0.1 ) {
			//go to end of loop
			Size lpos = 1;
			for ( ; ( lpos+pos <= ss_def.total_residue() ) && ( ss_def.loop_fraction( pos+lpos ) > 0.1); ++lpos ) {}
			if ( lpos > max_loop_size ) { //this loop has 4 or more residues
				unscored_loops.add_loop( pos, pos+lpos-1 );
			}
			pos+=lpos-1;
		} // have found a loop

		// look for short helices and store in short_helix
		if ( ss_def.helix_fraction( pos ) > 0.1 ) {
			Size hpos = 1;
			for ( ; ( hpos+pos <= ss_def.total_residue() ) && ( ss_def.helix_fraction( pos+hpos ) > 0.1); ++hpos ) {}
			if ( hpos <= max_short_helix   ) { //this helix has 5 or fewer residues
				for ( Size ipos = 0; ipos < hpos; ++ipos ) {
					short_helix[ pos+ipos] = true;
				}
			}
		}

		//finished parsing secondary structure definition
	}

	//elongate loops if they are terminated by a short helix
	loops::Loops removed_short_helices( unscored_loops );
	for ( loops::Loops::const_iterator it=unscored_loops.begin(); it != unscored_loops.end(); ++it ) {
		Size npos( it->stop() + 1 );
		while ( short_helix[ npos ] ) {
			removed_short_helices.add_loop( npos-1, npos );
			npos++;
		}
	}

	scored_core = removed_short_helices.invert( ss_def.total_residue() );
}

/// @brief Extract secondary structure chunks from the secondary structure
protocols::loops::Loops extract_secondary_structure_chunks(core::pose::Pose const & pose,
	char const extracted_ss_type) {
	using namespace core;
	using protocols::loops::Loops;
	Loops secondary_structure_chunks;

	bool ss_chunk_started = false;
	Size chunk_start_seqpos(0);
	Size chunk_end_seqpos(0);

	for ( core::Size ires = 1; ires <= pose.total_residue(); ++ires ) {
		if ( !ss_chunk_started ) {
			if ( pose.secstruct(ires) == extracted_ss_type ) {
				ss_chunk_started = true;
				chunk_start_seqpos = ires;
			}
		} else {
			if ( pose.secstruct(ires) != extracted_ss_type ) {
				ss_chunk_started = false;
				chunk_end_seqpos = ires - 1;
				secondary_structure_chunks.add_loop( chunk_start_seqpos, chunk_end_seqpos);
			} else if ( ! pose.residue_type(ires).is_protein() ) {
				ss_chunk_started = false;
				chunk_end_seqpos = ires - 1;
				secondary_structure_chunks.add_loop( chunk_start_seqpos, chunk_end_seqpos);
			}
		}
	}

	// if the input sequence ends with the last ss chunk
	if ( ss_chunk_started ) {
		chunk_end_seqpos = pose.total_residue();
		secondary_structure_chunks.add_loop( chunk_start_seqpos, chunk_end_seqpos);
	}
	return secondary_structure_chunks;
}

protocols::loops::Loops split_by_resSeq(core::pose::Pose const & pose) {
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	Loops chunks;

	//fpd don't include ligand
	core::Size nres = pose.total_residue();
	while ( nres>1 && !pose.residue(nres).is_protein() ) nres--;

	Loop new_loop;
	new_loop.set_start(1);
	new_loop.set_stop(nres);
	chunks.add_loop(new_loop);

	chunks = split_by_resSeq(pose, chunks);
	return chunks;
}

protocols::loops::Loops find_non_protein_chunks(core::pose::Pose const & pose) {
	Loops chunks;
	Loop new_loop;
	bool chunk_started = false;

	for ( core::Size ires = 1; ires < pose.total_residue(); ++ires ) {
		if ( pose.residue_type(ires).is_protein() ) continue;
		if ( !chunk_started ) {
			new_loop.set_start(ires);
			chunk_started = true;
		}

		if ( !pose.residue_type(ires).is_polymer() || pose.residue_type(ires).is_upper_terminus() ) {
			chunk_started = false;
			new_loop.set_stop(ires);
			chunks.add_loop(new_loop);
		}
	}
	return chunks;
}

protocols::loops::Loops split_by_resSeq(core::pose::Pose const & pose,
	protocols::loops::Loops const & input_chunks) {
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	Loops chunks;

	Loops::LoopList::const_iterator eit, it;
	for (  it = input_chunks.begin(), eit = input_chunks.end();
			it != eit; ++it ) {

		Loop new_loop(*it);

		for ( core::Size ires = it->start(); ires < it->stop(); ++ires ) {
			if ( pose.pdb_info()->number(ires+1) - pose.pdb_info()->number(ires) != 1 ||
					pose.pdb_info()->chain(ires+1) != pose.pdb_info()->chain(ires) ) {
				new_loop.set_stop(ires);
				chunks.add_loop(new_loop);

				new_loop.set_start(ires+1);
				new_loop.set_stop(it->stop());
			}
		}
		chunks.add_loop( new_loop );
	}
	return chunks;

}

protocols::loops::Loops split_by_ca_ca_dist(core::pose::Pose const & pose,
	protocols::loops::Loops const & input_chunks,
	core::Real const CA_CA_distance_cutoff) {
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	Loops continuous_chunks;

	Loops::LoopList::const_iterator eit, it;
	for (  it = input_chunks.begin(), eit = input_chunks.end();
			it != eit; ++it ) {

		Loop new_loop(*it);

		for ( core::Size ires = it->start(); ires < it->stop(); ++ires ) {
			if ( ! pose.residue(ires).is_protein() ) continue;

			if ( pose.residue(ires+1).is_protein() ) {
				if ( pose.residue(ires).xyz("CA").distance(pose.residue(ires+1).xyz("CA")) > CA_CA_distance_cutoff ) {
					new_loop.set_stop(ires);
					continuous_chunks.add_loop(new_loop);

					new_loop.set_start(ires+1);
					new_loop.set_stop(it->stop());
				}
			} else {
				new_loop.set_stop(ires);
				continuous_chunks.add_loop(new_loop);

				new_loop.set_start(ires+1);
				new_loop.set_stop(it->stop());
			}
		}
		continuous_chunks.add_loop( new_loop );
	}
	return continuous_chunks;
}

// gap_size- if two chunks are seperated by a gap of this size (or less), consider them one big chunk
protocols::loops::Loops remove_small_gaps(protocols::loops::Loops const & input_chunks, core::Size gap_size) {
	using protocols::loops::Loops;
	using protocols::loops::Loop;
	Loops secondary_structure_chunks;

	// join ss_chunks seperated by a small gap
	core::Size i_chunk = 1;
	while ( i_chunk <= input_chunks.num_loop() ) {
		Loop new_loop(input_chunks[i_chunk]);
		while ( i_chunk < input_chunks.num_loop() ) {
			const core::Size gap_length = input_chunks[i_chunk+1].start() - input_chunks[i_chunk].stop() - 1;
			if ( gap_length <= gap_size ) {
				new_loop.set_stop(input_chunks[i_chunk+1].stop());
				++i_chunk;
			} else {
				break;
			}
		}
		secondary_structure_chunks.add_loop(new_loop);
		++i_chunk;
	}
	return secondary_structure_chunks;
}

protocols::loops::Loops remove_short_chunks(protocols::loops::Loops const & input_chunks, core::Size minimum_length_of_chunk) {
	using protocols::loops::Loops;
	Loops secondary_structure_chunks;

	//remove short ss_chunks
	Loops::LoopList::const_iterator eit, it;
	for (  it = input_chunks.begin(), eit = input_chunks.end();
			it != eit; ++it ) {
		if ( it->size() >= minimum_length_of_chunk ) {
			secondary_structure_chunks.add_loop(*it);
		}
	}
	return secondary_structure_chunks;
}

// gap_size- if two chunks are seperated by a gap of this size (or less), consider them one big chunk
protocols::loops::Loops extract_secondary_structure_chunks(core::pose::Pose const & pose,
	std::string extracted_ss_types,
	core::Size gap_size,
	core::Size minimum_length_of_chunk_helix,
	core::Size minimum_length_of_chunk_strand,
	core::Real CA_CA_distance_cutoff) {
	using protocols::loops::Loops;
	Loops secondary_structure_chunks;

	for ( core::Size i_ss = 0; i_ss < extracted_ss_types.size(); ++i_ss ) {
		char ss = extracted_ss_types[i_ss];
		Loops secondary_structure_chunks_this_ss;

		// this order might be the best to deal with chain breaks in the middle of secondary structure chunk
		secondary_structure_chunks_this_ss = extract_secondary_structure_chunks(pose, ss);
		secondary_structure_chunks_this_ss = remove_small_gaps(secondary_structure_chunks_this_ss, gap_size);
		secondary_structure_chunks_this_ss = split_by_resSeq(pose, secondary_structure_chunks_this_ss);
		secondary_structure_chunks_this_ss = split_by_ca_ca_dist(pose, secondary_structure_chunks_this_ss, CA_CA_distance_cutoff);
		if ( ss == 'H' ) {
			secondary_structure_chunks_this_ss = remove_short_chunks(secondary_structure_chunks_this_ss, minimum_length_of_chunk_helix);
		}
		if ( ss == 'E' ) {
			secondary_structure_chunks_this_ss = remove_short_chunks(secondary_structure_chunks_this_ss, minimum_length_of_chunk_strand);
		}

		Loops::LoopList::const_iterator eit, it;
		for (  it = secondary_structure_chunks_this_ss.begin(), eit = secondary_structure_chunks_this_ss.end();
				it != eit; ++it ) {
			secondary_structure_chunks.add_loop(*it);
		}
	}

	// insert non protein chunks
	/*
	Loops non_prot_chunks = find_non_protein_chunks(pose);
	Loops::LoopList::const_iterator eit, it;
	for (  it = non_prot_chunks.begin(), eit = non_prot_chunks.end();
	it != eit; ++it ) {
	secondary_structure_chunks.add_loop(*it);
	}
	*/
	return secondary_structure_chunks;
}

protocols::loops::Loops extract_continuous_chunks(core::pose::Pose const & pose,
	core::Size const minimum_size,
	core::Real const CA_CA_distance_cutoff) {

	using protocols::loops::Loops;
	Loops chunks;
	chunks = split_by_resSeq(pose);
	chunks = split_by_ca_ca_dist(pose, chunks, CA_CA_distance_cutoff);
	chunks = remove_short_chunks(chunks, minimum_size);
	return chunks;
}

std::pair<bool, core::Size>
has_severe_pep_bond_geom_issues(
	const core::pose::Pose & pose,
	const Loop & loop,
	bool check_bonds, /*true*/
	bool check_angles, /*true*/
	core::Real max_c_n_dis, /*2.0*/
	core::Real allowed_ca_c_n_deviation, /* 25.0 */
	core::Real allowed_c_n_ca_deviation /* 25.0*/ ) {

	//JAB - these values are based on the CDL.  No peptide bond without severe chainbreaks or missing residues should have values
	// out of this range.  This can also indicate verypoorly solved crystal structures.
	// Berkholz DS, Shapovalov MV, Dunbrack RL Jr, Karplus PA (2009)
	// Conformation dependence of backbone geometry in proteins. Structure 17: 1316–1325.
	//

	TR << "checking peptide bond geometry: " << std::endl;
	for ( core::Size i =  loop.start(); i <= loop.stop() - 1; ++i ) {
		std::pair<bool, core::Size> checked = has_severe_pep_bond_geom_issues(
			pose, i, check_bonds, check_angles, max_c_n_dis, allowed_ca_c_n_deviation, allowed_c_n_ca_deviation);

		if ( checked.first ) return checked;
		else {
			continue;
		}

	}

	return std::make_pair(false, 0);

}

std::pair<bool, core::Size>
has_severe_pep_bond_geom_issues(
	core::pose::Pose const & pose,
	core::Size resnum,
	bool check_bonds,
	bool check_angles,
	core::Real max_c_n_dis,
	core::Real allowed_ca_c_n_deviation,
	core::Real allowed_c_n_ca_deviation)
{

	//Skip measuring virtual residue geometry or non-protein residues.
	if ( pose.residue(resnum).is_virtual_residue() || pose.residue(resnum+1).is_virtual_residue() ) {
		TR << "Cannot measure geometry of virtual residue" << std::endl;
		return std::make_pair(false, 0);
	}
	if ( ! pose.residue(resnum).is_protein() || ! pose.residue(resnum+1).is_protein() ) {
		TR << "Cannot measure geometry of non-protein residue" << std::endl;
		return std::make_pair(false, 0);
	}

	//25 degree dif default from  min/max in  CDL. 1.0 A crystal structures. This should cover NMR and MD structures as well
	core::Real min_ca_c_n_ang = 114.5 - allowed_ca_c_n_deviation;
	core::Real max_ca_c_n_ang = 119.5 + allowed_ca_c_n_deviation;

	core::Real min_c_n_ca_ang =  120.0 - allowed_c_n_ca_deviation;
	core::Real max_c_n_ca_ang = 126.0 + allowed_c_n_ca_deviation;

	//core::id::AtomID n_0 = core::id::AtomID(pose.residue(i).atom_index("N"), i);
	numeric::xyzVector<Real> const ca_0 = pose.residue(resnum).atom("CA").xyz();
	numeric::xyzVector<Real> const c_0 = pose.residue(resnum).atom("C").xyz();


	numeric::xyzVector<Real> const n_1 = pose.residue(resnum +1).atom("N").xyz();
	numeric::xyzVector<Real> const ca_1 = pose.residue(resnum +1).atom("CA").xyz();
	//core::id::AtomID c_1 = core::id::AtomID(pose.residue(i +1).atom_index("C"), i +1);


	if ( check_bonds ) {
		core::Real c_n_dis =c_0.distance(n_1);
		if ( c_n_dis > max_c_n_dis ) {
			TR<< "PepBondGeom: "<< pose.pdb_info()->pose2pdb(resnum)<<" C-N --  " << c_n_dis << " max: " << max_c_n_dis << std::endl;
			return std::make_pair(true, resnum);
		}
	}

	if ( check_angles ) {

		core::Real ca_c_n_ang = numeric::angle_degrees(ca_0, c_0, n_1);
		core::Real c_n_ca_ang = numeric::angle_degrees(c_0, n_1, ca_1);
		//core::Real ca_c_n_ang = numeric::conversions::degrees(pose.conformation().bond_angle(ca_0, c_0, n_1));
		//core::Real c_n_ca_ang = numeric::conversions::degrees(pose.conformation().bond_angle(c_0, n_1, ca_1));

		if ( ca_c_n_ang < min_ca_c_n_ang || ca_c_n_ang > max_ca_c_n_ang ) {
			TR<< "PepBondGeom: " << resnum <<" - " << pose.pdb_info()->pose2pdb(resnum)
				<< " Ca-C-N --  " << min_ca_c_n_ang <<" (min): "<< ca_c_n_ang <<" : "<< max_ca_c_n_ang << " (max)"<< std::endl;
			return std::make_pair(true, resnum);

		} else if ( c_n_ca_ang < min_c_n_ca_ang || c_n_ca_ang > max_c_n_ca_ang ) {
			TR << "PepBondGeom: " << resnum << " - " << pose.pdb_info()->pose2pdb(resnum)
				<< " C-N-Ca --  "<< min_c_n_ca_ang <<" (min) : " << c_n_ca_ang << " : "<< max_c_n_ca_ang << " (max)" << std::endl;
			return std::make_pair(true, resnum);

		}


	}

	return std::make_pair(false, 0);
}


} // loops
} // protocols


