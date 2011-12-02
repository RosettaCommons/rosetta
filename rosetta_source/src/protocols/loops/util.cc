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
#include <protocols/loops/Loops.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/loops/loops_main.hh> //for getting ss from dssp

#include <core/fragment/SecondaryStructure.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>

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

//// C++ headers
#include <list>

//Auto Headers
#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <protocols/rosetta_scripts/util.hh>
#include <utility/string_util.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS

static basic::Tracer TR("protocols.loops.util");

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
      conformation::Residue const & rsd( pose.residue( pos ) );
			if ( bCstAllAtom ) {
				for ( Size ii = 1; ii<= rsd.natoms(); ++ii ) {
					pose.add_constraint( new scoring::constraints::CoordinateConstraint(
	      id::AtomID( ii, pos),
	      id::AtomID( 1, pos ) /*this is completely ignored! */,
	      rsd.xyz( ii ),
	      new scoring::constraints::HarmonicFunc( 0.0, coord_sdev )
						) );
				}
      } else {
				id::AtomID atomID( pose.residue_type(pos).atom_index("CA"), pos );
				pose.add_constraint( new scoring::constraints::CoordinateConstraint(
						atomID,
						id::AtomID( 1, pos ) /*this is completely ignored! */,
						rsd.xyz( atomID.atomno() ),
						new scoring::constraints::HarmonicFunc( 0.0, coord_sdev )
					) );
      }
    }
  }
}

///@brief get frags that are fully within the Loop --- shorten(=true/false) frags that are close to the end of loops.
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
  if (loops.empty())
    return;

  set_extended_torsions_and_idealize_loops(*pose, loops);
}

void set_extended_torsions_and_idealize_loops(core::pose::Pose& pose, loops::Loops loops) {
  using core::Size;
  using core::conformation::Conformation;
  using protocols::loops::Loop;
  using protocols::loops::Loops;

	// if no loops, we want to have extended structure everywhere.
	// it is a by-value parameter -- as intenden the change is kept local
	if (loops.empty())
		loops.add_loop(1, pose.total_residue(), 0);

	TR.Debug << "extend structure for " << loops << std::endl;

  Conformation& conf = pose.conformation();
  for (Loops::const_iterator i = loops.begin(); i != loops.end(); ++i) {
    const Loop& loop = *i;

    for (Size j = loop.start(); j <= loop.stop(); ++j) {
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
	for( Size i = 1; i < nres; i ++ ){
		if( loops.is_loop_residue(i) ) all_loop_list.push_back( i );
	}
	scorefxn(pose);
	setPoseExtraScores( pose, "ScoreCore", 	scorefxn.get_sub_score_exclude_res( pose, all_loop_list ) );

	Real score =  scorefxn(pose);

	for( Size l = 1; l <= nloops; l++ ){
		if( l > loops.size() ){
			setPoseExtraScores( pose, "ScoreLoopI" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			setPoseExtraScores( pose, "ScoreLoopC" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			setPoseExtraScores( pose, "ScoreLoopL" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			continue;
		}
		utility::vector1< core::Size > loop_list;
		utility::vector1< core::Size > non_loop_list;
		for( Size i = 1; i < nres; i ++ ){
			if( ( i < loops[l].start() ) || ( i > loops[l].stop() ) ){
				loop_list.push_back( i );
			}else{
				non_loop_list.push_back( i );
			}
		}

		Real loopscore = scorefxn.get_sub_score_exclude_res( pose, loop_list );
		Real nonloopscore = scorefxn.get_sub_score_exclude_res( pose, non_loop_list );
		setPoseExtraScores( pose, "ScoreLoopI" + ObjexxFCL::right_string_of(l,3,'0'), loopscore );
		setPoseExtraScores( pose, "ScoreLoopC" + ObjexxFCL::right_string_of(l,3,'0'), score - loopscore - nonloopscore );
		setPoseExtraScores( pose, "ScoreLoopL" + ObjexxFCL::right_string_of(l,3,'0'), score - nonloopscore );
	}

	// Work out RMS values too

	core::pose::Pose native_pose_super = native_pose;
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, native_pose_super, core::id::BOGUS_ATOM_ID );
	for ( core::Size ir=1; ir <= native_pose.total_residue(); ++ir ) {
		runtime_assert( ir <=  pose.total_residue() );
		runtime_assert( ir <=  native_pose_super.total_residue() );
		if( ( !loops.is_loop_residue( ir ) ) && pose.residue(ir).is_protein() ) {
			id::AtomID const id1( native_pose_super.residue(ir).atom_index("CA"), ir );
			id::AtomID const id2( pose.residue(ir).atom_index("CA"), ir );
			atom_map.set(id1, id2);
		}
	}
	core::scoring::superimpose_pose( native_pose_super, pose, atom_map );

	int corelength;
	setPoseExtraScores(	pose, "corerms", native_loop_core_CA_rmsd( native_pose, pose, loops, corelength )	);

	for( Size l = 1; l <= nloops; l++ ){
		if( l > loops.size() ){
			setPoseExtraScores( pose, "RMSLoop" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			continue;
		}
		Loops temploops;
		temploops.add_loop( loops[l] );
		setPoseExtraScores( pose, "RMSLoop" + ObjexxFCL::right_string_of(l,3,'0'),
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
	return core::scoring::getScoreFunction();
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
  for ( Size i = 1; i<= (Size)nres; ++i ) {
    if ( i==(Size)pose.fold_tree().root() ) continue;
    if( coordconstraint_segments.is_loop_residue( i ) ) {
      Residue const & nat_i_rsd( pose.residue(i) );
      for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
        pose.add_constraint( new CoordinateConstraint(
          AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
          new HarmonicFunc( 0.0, coord_sdev ) ) );
      }
    }
  }

}

Loops
loops_from_string( std::string const loop_str, core::pose::Pose const & pose ){
  utility::vector1< std::string > const loops_vec( utility::string_split( loop_str, ',' ) );

// each loop should have the format loop_start:loop_end:cut
// if cut is not set then it's taken to be 0. Residue numbering can follow the
// pdb numbering
  Loops loops_from_tag;
 	foreach( std::string const residue_pair, loops_vec ){
    utility::vector1< std::string > const residues( utility::string_split( residue_pair, ':' ) );
    runtime_assert( residues.size() == 2 || residues.size() == 3 );
    core::Size const loop_start( protocols::rosetta_scripts::parse_resnum( residues[ 1 ], pose ) );
    core::Size const loop_stop( protocols::rosetta_scripts::parse_resnum( residues[ 2 ], pose ) );
    core::Size loop_cut( 0 );
	  if( residues.size() == 3 )
      loop_cut = protocols::rosetta_scripts::parse_resnum( residues[ 3 ], pose );
    runtime_assert( loop_start <= loop_stop );
	  runtime_assert( loop_start >= 1 );
    runtime_assert( loop_stop <= pose.total_residue() );
    loops_from_tag.add_loop( loop_start, loop_stop, loop_cut );
	}
	return( loops_from_tag );
}

void define_scorable_core_from_secondary_structure(
   core::fragment::SecondaryStructure const& ss_def,
	 protocols::loops::Loops& scored_core )
{
	using namespace core;
	using namespace basic::options;
	//	Size const max_loop_size( 3 );
	//	Size const max_short_helix( 5 );
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

	scored_core =	removed_short_helices.invert( ss_def.total_residue() );
}



} // loops
} // protocols
