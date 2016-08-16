// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilot/kblacklock/loodo.cc
///
/// @brief Pilot Application for LooDo Algorithm.
/// @details This algorithm, named LooDo, utilizes large libraries of peptide fragments to perform loop-directed domain placement searches of an insert domain with respect to a parent domain.
///
/// @author Kristin Blacklock (kristin.blacklock@rutgers.edu)
/// @note Last Modified: 9/23/2015

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loodo.OptionKeys.gen.hh>
#include <basic/prof.hh>

#include <core/conformation/Residue.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/MatcherTask.hh>
#include <protocols/match/MatchSet.hh>
#include <protocols/match/OccupiedSpaceHash.hh>
#include <protocols/match/downstream/ActiveSiteGrid.hh>
#include <protocols/sic_dock/util.hh>

// Numeric headers
#include <numeric/geometry/BoundingBox.fwd.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/LexicographicalIterator.hh>
#include <utility/OrderedTuple.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <map>
#include <iostream>
#include <string>
#include <fstream>
#include <iostream>

#include <ObjexxFCL/string.functions.hh>

static THREAD_LOCAL basic::Tracer TR( "app.loodo" );

//protocols::forge::components::VarLengthBuild::FrameList
//pick_fragments(
// std::string const & complete_ss,
// std::string const & complete_aa,
// utility::vector1< std::string > const & complete_abego,
// protocols::forge::build::Interval const & interval,
// core::Size const frag_length,
// core::Size const n_frags
//)
//{
//   using core::fragment::Frame;
//  using core::fragment::FrameOP;
//  using core::fragment::IndependentBBTorsionSRFD;
//  using core::fragment::picking_old::vall::pick_fragments;
//
//  FrameList frames;
//
//  for ( Size j = 0, je = interval.length(); j < je; ++j ) {
//    TR << "picking " << n_frags << " " << frag_length << "-mers for position " << ( interval.left + j ) << std::endl;
//    String ss_sub = complete_ss.substr( interval.left + j - 1, frag_length );
//    if ( ss_sub.length() < frag_length ) {
//      ss_sub.append( frag_length - ss_sub.length(), 'D' );
//    }
//    String aa_sub;
//    if ( !complete_aa.empty() ) {
//      aa_sub = complete_aa.substr( interval.left + j - 1, frag_length );
//      if ( aa_sub.length() < frag_length ) {
//        aa_sub.append( frag_length - aa_sub.length(), '.' );
//      }
//    } else {
//      aa_sub = "";
//    }
//    utility::vector1< String > abego_sub;
//    if ( complete_abego.size() > 0 ) {
//      runtime_assert( complete_ss.length() == complete_abego.size() );
//      Size pos( 1 );
//      abego_sub.resize( frag_length );
//      for( Size ii = interval.left + j; ii <= interval.left + j + frag_length - 1; ++ii, ++pos ) {
//        if ( ii > complete_abego.size() ) {
//          abego_sub[ pos ] = "X";
//        } else {
//          abego_sub[ pos ] = complete_abego[ ii ];
//        }
//      }
//    } else {
//      abego_sub.clear();
//    }
//    FrameOP frame( new Frame( interval.left + j, frag_length ) );
//    frame->add_fragment( pick_fragments( ss_sub, aa_sub, abego_sub, n_frags, true, IndependentBBTorsionSRFD() ) );
//    frames.push_back( frame );
//  }
//  return frames;
//}

class VarLengthBuildKristin : public protocols::forge::components::VarLengthBuild {
public:
	VarLengthBuild::FrameList
	public_pick_fragments(
		String const & complete_ss,
		std::string const & complete_aa,
		utility::vector1< String > const & complete_abego,
		Interval const & interval,
		core::Size const frag_length,
		core::Size const n_frags
	)
	{
		protocols::forge::components::VarLengthBuild::FrameList frames = pick_fragments( complete_ss, complete_aa, complete_abego, interval, frag_length, n_frags );
		return frames;
	}
};


///// generate fragment library.
core::fragment::ConstantLengthFragSetOP
get_vallfrags( int round,
	core::Size fraglength,
	std::string cap_pose_ss,
	std::string cap_pose_aa,
	std::string bot_pose_ss,
	std::string bot_pose_aa,
	int ins_begin,
	int ins_end )
{

	utility::vector1< std::string > abego;
	std::string ss;
	std::string aa;

	// round 1 = linkerA
	if ( round == 1 ) {
		//TR << "Round " << round << std::endl;

		// Generate SS sequence of linker, where two flanking residues at either end
		// have native SS identities.
		ss += bot_pose_ss.at(ins_begin-2);
		TR << "SS (+bot @ ins_begin-2)= " << ss << std::endl;
		ss += bot_pose_ss.at(ins_begin-1);
		TR << "SS (+bot @ ins_begin-1)= " << ss << std::endl;
		for ( core::Size j=1; j<=fraglength; j++ ) {
			ss += "L";
		}
		ss += cap_pose_ss.at(0);
		ss += cap_pose_ss.at(1);


		// Generate AA sequence of linker, where two flanking residues at either end
		// have native AA identity.
		aa += bot_pose_aa.at(ins_begin-2);
		TR << "AA (+bot @ ins_begin-2)= " << aa << std::endl;
		aa += bot_pose_aa.at(ins_begin-1);
		TR << "AA (+bot @ ins_begin-1)= " << aa << std::endl;
		for ( core::Size k = 1; k <= fraglength; k++ ) {
			aa += "B";
		}
		TR << "AA (+B x fraglength)= " << aa << std::endl;
		aa += cap_pose_aa.at(0);
		TR << "AA (+cap @ 0)= " << aa << std::endl;
		aa += cap_pose_aa.at(1);
		TR << "AA (+cap @ 1)= " << aa << std::endl;

		// Generate ABEGO sequence based on SS sequence.
		for ( unsigned s = 0; s <= ss.size()-1; s++ ) {
			if ( ss.at(s) == 'H' ) {
				abego.push_back("A");
			} else if ( ss.at(s) == 'E' ) {
				abego.push_back("B");
			} else {
				abego.push_back("X");
			}
		}
	}

	// round 2 = linker b
	if ( round == 2 ) {

		// Generate SS sequence of linker, where two flanking residues at either end
		// have native SS identities.
		ss += cap_pose_ss.at(cap_pose_ss.size()-2);
		ss += cap_pose_ss.at(cap_pose_ss.size()-1);
		for ( core::Size j = 1; j <= fraglength; j++ ) {
			ss += "L";
		}
		ss += bot_pose_ss.at(ins_end-1);
		ss += bot_pose_ss.at(ins_end);


		// Generate AA sequence of linker, where two flanking residues at either end
		// have native AA identity.
		aa += cap_pose_aa.at(cap_pose_aa.size()-2);
		aa += cap_pose_aa.at(cap_pose_aa.size()-1);
		for ( core::Size k = 1; k <= fraglength; k++ ) {
			aa += "B";
		}
		aa += bot_pose_aa.at(ins_end-1);
		aa += bot_pose_aa.at(ins_end);


		// Generate ABEGO sequence based on SS sequence.
		for ( unsigned t = 0; t <= ss.size()-1; t++ ) {
			if ( ss.at(t) == 'H' ) {
				abego.push_back("A");
			} else {
				abego.push_back("B");
			}
		}
	}

	protocols::forge::build::Interval interval(1,fraglength+4);

	core::Size num_frags = basic::options::option[ basic::options::OptionKeys::loodo::num_frags ]();

	core::fragment::ConstantLengthFragSetOP fraglib;
	fraglib = core::fragment::ConstantLengthFragSetOP ( new core::fragment::ConstantLengthFragSet( fraglength+4 ) );

	// If the user has specified a fragment library, use this instead of the vall.
	if ( basic::options::option[basic::options::OptionKeys::loodo::use_fraglib ].user() ) {
		TR << "Using user-defined fragment library file: " << basic::options::option[basic::options::OptionKeys::loodo::use_fraglib ]() << std::endl;
		fraglib->read_fragment_file( basic::options::option[basic::options::OptionKeys::loodo::use_fraglib ](), num_frags, 1, false );
	} else {
		//else if ( basic::options::Option[basic::options::OptionKeys::loodo::use_fraglibsc ].user() ){
		//  TR << "using user-defined fragment library file with sidechain torsions: " << basic::options::Option[basic::options::OptionKeys::loodo::use_fraglibsc ]() << std::endl;
		//  //fraglib->read_fragment_file( basic::options::Option[basic::options::OptionKeys::loodo::use_fraglibsc ](), num_frags, 1, false, true );
		//}
		//fraglib->add( VarLengthBuildKristin::public_pick_fragments( ss, aa, abego, interval, fraglength+4, num_frags) );
		VarLengthBuildKristin vlb;
		fraglib->add( vlb.public_pick_fragments( ss, aa, abego, interval, fraglength+4, num_frags) );

	}
	//TR << "FragLib has " << fraglib->size() << " fragments." << std::endl;
	//if (debug) { TR << "Returning fraglib." << std::endl; }

	return fraglib;
}

core::fragment::ConstantLengthFragSetOP
molten_get_vallfrags( core::Size fraglength,
	std::string frag_ss,
	std::string frag_aa,
	int known)
{
	TR << "Native fragment secondary structure: " << frag_ss << std::endl;
	TR << "Native fragment AA sequence: " << frag_aa << std::endl;

	utility::vector1< std::string > abego;
	std::string ss;
	std::string aa;

	core::Size nonmolten = known/2;

	//SS SEQUENCE //
	ss += frag_ss.at(0);
	ss += frag_ss.at(1);
	for ( core::Size i = 1; i <= nonmolten; i++ ) {
		ss += frag_ss.at(1+i);
	}
	for ( core::Size j=1; j <= fraglength-known; j++ ) {
		ss += "L";
	}
	for ( core::Size k=nonmolten; k >= 1; k-- ) {
		ss += frag_ss.at(6-k);
	}
	ss+= frag_ss.at(6);
	ss += frag_ss.at(7);

	TR << "Secondary Structure: " << ss << std::endl;

	// AA SEQUENCE //
	aa += frag_aa.at(0);
	aa += frag_aa.at(1);
	for ( core::Size i = 1; i <= nonmolten; i++ ) {
		aa += frag_aa.at(1+i);
	}
	for ( core::Size j=1; j <= fraglength-known; j++ ) {
		aa += "B";
	}
	for ( core::Size k=nonmolten; k >= 1; k-- ) {
		aa += frag_aa.at(6-k);
	}
	aa += frag_aa.at(6);
	aa += frag_aa.at(7);

	TR << "Amino Acid Seq: " << aa << std::endl;

	for ( unsigned s = 0; s <= ss.size()-1; s++ ) {
		if ( ss.at(s) == 'H' ) {
			abego.push_back("A");
		} else if ( ss.at(s) == 'E' ) {
			abego.push_back("B");
		} else {
			abego.push_back("X");
		}
	}

	protocols::forge::build::Interval interval(1,fraglength+4);

	core::Size num_frags = basic::options::option[ basic::options::OptionKeys::loodo::num_frags ]();

	core::fragment::ConstantLengthFragSetOP fraglib;
	fraglib = core::fragment::ConstantLengthFragSetOP ( new core::fragment::ConstantLengthFragSet( fraglength+4 ) );

	if ( basic::options::option[basic::options::OptionKeys::loodo::use_fraglib ].user() ) {
		fraglib->read_fragment_file( basic::options::option[basic::options::OptionKeys::loodo::use_fraglib ](), num_frags, 1, false );
	} else {
		//else if ( basic::options::option[basic::options::OptionKeys::loodo::use_fraglibsc ].user() ) {
		//  fraglib->read_fragment_file( basic::options::option[basic::options::OptionKeys::loodo::use_fraglibsc ](), num_frags, 1, false, true );
		//}
		//protocols::forge::components::VarLengthBuild vlb;
		//fraglib->add(vlb.pick_fragments( ss, aa, abego, interval, fraglength+4, num_frags));
		VarLengthBuildKristin vlb;
		fraglib->add( vlb.public_pick_fragments( ss, aa, abego, interval, fraglength+4, num_frags) );
	}

	return fraglib;
}

///// Terminal-residue iterative superimposition.
core::pose::Pose
perturb( std::string round,
	std::string pair,
	core::pose::Pose & stationary_pose,
	core::pose::Pose & moving_pose,
	int ins_begin,
	int ins_end )
{

	if ( round == "A" ) {
		if ( pair == "ftb" ) {
			std::map< core::id::AtomID, core::id::AtomID > res_map_ftb;

			//res_map_ftb[1] = ins_begin-1; //[frag] = bot
			core::id::AtomID f1_bb1(  1 , 1 );
			core::id::AtomID b1_bb1(  1 , ins_begin-1 );
			core::id::AtomID f1_bb2(  2 , 1 );
			core::id::AtomID b1_bb2(  2 , ins_begin-1 );
			core::id::AtomID f1_bb3(  3 , 1 );
			core::id::AtomID b1_bb3(  3 , ins_begin-1 );
			res_map_ftb[ f1_bb1 ] = b1_bb1;
			res_map_ftb[ f1_bb2 ] = b1_bb2;
			res_map_ftb[ f1_bb3 ] = b1_bb3;

			//res_map_ftb[2] = ins_begin;
			core::id::AtomID f2_bb1( 1 , 2 );
			core::id::AtomID b2_bb1( 1 , ins_begin );
			core::id::AtomID f2_bb2( 2 , 2 );
			core::id::AtomID b2_bb2( 2 , ins_begin );
			core::id::AtomID f2_bb3( 3 , 2 );
			core::id::AtomID b2_bb3( 3 , ins_begin );
			res_map_ftb[ f2_bb1 ] = b2_bb1;
			res_map_ftb[ f2_bb2 ] = b2_bb2;
			res_map_ftb[ f2_bb3 ] = b2_bb3;

			core::scoring::superimpose_pose( moving_pose, stationary_pose, res_map_ftb);
			return moving_pose;
		}
		if ( pair == "ctf" ) {
			std::map< core::id::AtomID, core::id::AtomID > res_map_ctf;

			//res_map_ctf[2] = stationary_pose.total_residue(); //[cap] = frag
			core::id::AtomID c1_bb1( 1 , 2 );
			core::id::AtomID f1_bb1( 1 , stationary_pose.total_residue() );
			core::id::AtomID c1_bb2( 2 , 2 );
			core::id::AtomID f1_bb2( 2 , stationary_pose.total_residue() );
			core::id::AtomID c1_bb3( 3 , 2 );
			core::id::AtomID f1_bb3( 3 , stationary_pose.total_residue() );
			res_map_ctf[ c1_bb1 ] = f1_bb1;
			res_map_ctf[ c1_bb2 ] = f1_bb2;
			res_map_ctf[ c1_bb3 ] = f1_bb3;

			//res_map_ctf[1] = stationary_pose.total_residue()-1;
			core::id::AtomID c2_bb1( 1 , 1 );
			core::id::AtomID f2_bb1( 1 , stationary_pose.total_residue()-1 );
			core::id::AtomID c2_bb2( 2 , 1 );
			core::id::AtomID f2_bb2( 2 , stationary_pose.total_residue()-1 );
			core::id::AtomID c2_bb3( 3 , 1 );
			core::id::AtomID f2_bb3( 3 , stationary_pose.total_residue()-1 );
			res_map_ctf[ c2_bb1 ] = f2_bb1;
			res_map_ctf[ c2_bb2 ] = f2_bb2;
			res_map_ctf[ c2_bb3 ] = f2_bb3;

			core::scoring::superimpose_pose( moving_pose, stationary_pose, res_map_ctf);
			return moving_pose;
		}

	} else {
		if ( pair == "ftb" ) {
			std::map< core::id::AtomID, core::id::AtomID > res_map_ftb;

			//res_map_ftb[moving_pose.total_residue()] = ins_end+1; //[frag] = bot
			core::id::AtomID f1_bb1( 1 , moving_pose.total_residue() );
			core::id::AtomID b1_bb1( 1 , ins_end+1 );
			core::id::AtomID f1_bb2( 2 , moving_pose.total_residue() );
			core::id::AtomID b1_bb2( 2 , ins_end+1 );
			core::id::AtomID f1_bb3( 3 , moving_pose.total_residue() );
			core::id::AtomID b1_bb3( 3 , ins_end+1 );
			res_map_ftb[ f1_bb1 ] = b1_bb1;
			res_map_ftb[ f1_bb2 ] = b1_bb2;
			res_map_ftb[ f1_bb3 ] = b1_bb3;

			//res_map_ftb[moving_pose.total_residue()-1] = ins_end;
			core::id::AtomID f2_bb1( 1 , moving_pose.total_residue()-1 );
			core::id::AtomID b2_bb1( 1 , ins_end );
			core::id::AtomID f2_bb2( 2 , moving_pose.total_residue()-1 );
			core::id::AtomID b2_bb2( 3 , ins_end );
			core::id::AtomID f2_bb3( 3 , moving_pose.total_residue()-1 );
			core::id::AtomID b2_bb3( 3 , ins_end );
			res_map_ftb[ f2_bb1 ] = b2_bb1;
			res_map_ftb[ f2_bb2 ] = b2_bb2;
			res_map_ftb[ f2_bb3 ] = b2_bb3;

			core::scoring::superimpose_pose( moving_pose, stationary_pose, res_map_ftb);
			return moving_pose;
		}
		if ( pair == "ctf" ) {
			std::map< core::id::AtomID, core::id::AtomID > res_map_ctf;

			//res_map_ctf[moving_pose.total_residue()] = 2; //[cap] = frag
			core::id::AtomID c1_bb1( 1 , moving_pose.total_residue() );
			core::id::AtomID f1_bb1( 1 , 2 );
			core::id::AtomID c1_bb2( 2 , moving_pose.total_residue() );
			core::id::AtomID f1_bb2( 2 , 2 );
			core::id::AtomID c1_bb3( 3 , moving_pose.total_residue() );
			core::id::AtomID f1_bb3( 3 , 2 );
			res_map_ctf[ c1_bb1 ] = f1_bb1;
			res_map_ctf[ c1_bb2 ] = f1_bb2;
			res_map_ctf[ c1_bb3 ] = f1_bb3;

			//res_map_ctf[moving_pose.total_residue()-1] = 1;
			core::id::AtomID c2_bb1( 1 , moving_pose.total_residue()-1 );
			core::id::AtomID f2_bb1( 1 , 1 );
			core::id::AtomID c2_bb2( 2 , moving_pose.total_residue()-1 );
			core::id::AtomID f2_bb2( 2 , 1 );
			core::id::AtomID c2_bb3( 3 , moving_pose.total_residue()-1 );
			core::id::AtomID f2_bb3( 3 , 1 );
			res_map_ctf[ c2_bb1 ] = f2_bb1;
			res_map_ctf[ c2_bb2 ] = f2_bb2;
			res_map_ctf[ c2_bb3 ] = f2_bb3;

			core::scoring::superimpose_pose( moving_pose, stationary_pose, res_map_ctf);
			return moving_pose;
		}
	}
	return moving_pose; //Won't get here.
}

///// Make BumpGrid around parent domain.
protocols::match::BumpGrid
make_bb_grid( core::pose::Pose & bot_pose)
{
	protocols::match::BumpGrid bb_grid = *protocols::match::bump_grid_to_enclose_pose( bot_pose );
	utility::vector1< protocols::match::BumpGrid > original_scaffold_residue_bump_grids;
	original_scaffold_residue_bump_grids.resize( bot_pose.total_residue() );

	for ( core::Size ii = 1; ii <= bot_pose.total_residue(); ii++ ) {
		protocols::match::BumpGrid resbg = *bump_grid_to_enclose_residue_backbone( bot_pose.residue( ii ), bb_grid );
		fill_grid_with_backbone_heavyatom_spheres( bot_pose.residue( ii ), resbg ); //Changed All Radii to H_ARO (1.0)
		//fill_grid_with_residue_spheres( bot_pose.residue( ii ), resbg); //Takes Liggy Into Account
		bb_grid.or_with( resbg );
		original_scaffold_residue_bump_grids[ ii ] = resbg;
	}

	return bb_grid;
}

///// Make Placement Grid from gridlig file.
protocols::match::downstream::ActiveSiteGrid
make_active_site_grid()
{
	protocols::match::downstream::ActiveSiteGrid active_site_grid;
	std::string gridligpath = basic::options::option[ basic::options::OptionKeys::loodo::gridligpath ]();
	TR << "Making Active Site Grid from " << gridligpath << std::endl;
	active_site_grid.initialize_from_gridlig_file( gridligpath );

	return active_site_grid;
}

///// Find the pose-number of the residue in the pose that is closest to the center of mass.
core::Size
find_centerest_residue( core::pose::Pose & cap_pose )
{
	// find pose1 COM
	core::Real cap_pose_atom_count = 0; // for division later
	numeric::xyzVector<core::Real> cap_pose_COM(0,0,0);
	for ( core::Size ii = 1; ii <= cap_pose.total_residue(); ii++ ) {
		for ( core::Size jj = 1; jj <= cap_pose.residue_type(ii).natoms(); jj++ ) {
			cap_pose_COM += cap_pose.residue(ii).atom(jj).xyz();
			cap_pose_atom_count++;
		}
	}

	if ( cap_pose_atom_count < 1 ) {
		std::cerr << "\n" << std::endl;
		std::cerr << "WARNING!" << std::endl;
		std::cerr << "There are no CA atoms in the insert domain pose." << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	} else {
		cap_pose_COM /= cap_pose_atom_count;
	}

	core::Real min_dist = 10.0;
	core::Real new_dist = 100.0;
	core::Size res_i = 1;
	for ( core::Size ii = 1; ii <= cap_pose.total_residue(); ii++ ) {
		if ( cap_pose.residue(ii).is_protein() ) {
			new_dist = cap_pose_COM.distance( cap_pose.residue(ii).atom("CA").xyz() );
			if ( new_dist < min_dist ) {
				min_dist = new_dist;
				res_i = ii;
			}
		}
	}

	return res_i;
}

///// Check pose CA's for containment within the placement grid and occupancy within the bumpgrid.
bool
containment_and_clash_checks( protocols::match::BumpGrid &  bb_grid,
	protocols::match::downstream::ActiveSiteGrid & active_site_grid,
	core::pose::Pose &  cap_pose,
	core::Size & COM,
	core::Size & clash_fail,
	core::Size & cont_fail
	/*core::pose::Pose frag_pose*/ )
{
	// Bump check
	protocols::match::ProbeRadius radius_type = protocols::match::probe_radius_for_atom_type(1);
	for ( core::Size ii = 1; ii <= cap_pose.total_residue(); ++ii ) {
		if ( bb_grid.occupied( radius_type, cap_pose.residue(ii).atom("CA").xyz() ) ) {
			++clash_fail;
			return false;
		}
	}


	///// Containment Check (at least x% of Ca atoms lying in grid) /////
	if ( basic::options::option[ basic::options::OptionKeys::loodo::com_in_grid ]() ) {
		if ( !active_site_grid.occupied( cap_pose.residue( COM ).atom("CA").xyz() ) ) {
			++cont_fail;
			return false;
		}
	}

	///// Quantify ratio of CA atoms lying in placement grid.
	core::Size count = 0;
	for ( core::Size cur_res = 1; cur_res <= cap_pose.total_residue(); ++cur_res ) {
		if ( active_site_grid.occupied( cap_pose.residue( cur_res ).atom("CA").xyz() ) ) {
			count++;
		}
	}

	///// If ratio is less than threshold, this linker does not pass.
	core::Size ratio = count/cap_pose.total_residue();
	core::Real user_ratio = basic::options::option[ basic::options::OptionKeys::loodo::ca_ratio ]();
	if ( ratio >= user_ratio ) {
		return true;
	} else {
		++cont_fail;
		return false;
	}
}

numeric::geometry::hashing::Real6
real6_from_rt( core::kinematics::RT  & rt )
{
	numeric::geometry::hashing::Real6 rt_6;

	numeric::HomogeneousTransform<core::Real>  ht( rt.get_rotation() , rt.get_translation() );
	numeric::xyzVector<core::Real> euler_angles =  ht.euler_angles_deg();

	rt_6[1] = rt.get_translation().x();
	rt_6[2] = rt.get_translation().y();
	rt_6[3] = rt.get_translation().z();

	/// Euler angles lie in (0.0,360.0]
	rt_6[4] = fmod(euler_angles.x(),360.0);
	rt_6[5] = fmod(euler_angles.y(),360.0);
	rt_6[6] = fmod(euler_angles.z(),360.0);
	rt_6[4] = rt_6[4]<0.0 ? rt_6[4]+360.0 : rt_6[4];
	rt_6[5] = rt_6[5]<0.0 ? rt_6[5]+360.0 : rt_6[5];
	rt_6[6] = rt_6[6]<0.0 ? rt_6[6]+360.0 : rt_6[6];

	return rt_6;
}

numeric::geometry::hashing::Real6
six_coords_to_real6(numeric::xyzVector<core::Real> bot_ca_1,
	numeric::xyzVector<core::Real> bot_ca_2,
	numeric::xyzVector<core::Real> bot_ca_3,
	numeric::xyzVector<core::Real> cap_ca_1,
	numeric::xyzVector<core::Real> cap_ca_2,
	numeric::xyzVector<core::Real> cap_ca_3 )
{

	core::kinematics::Stub  bot_stub( bot_ca_1, bot_ca_2, bot_ca_3 );
	core::kinematics::Stub  cap_stub( cap_ca_1, cap_ca_2, cap_ca_3 );

	//TR << "Bot stub: " << bot_stub << std::endl;
	//TR << "Cap stub: " << cap_stub << std::endl;

	core::kinematics::RT rt;
	rt.from_stubs( cap_stub, bot_stub );

	return real6_from_rt( rt );
}


bool
real6_similarity( numeric::geometry::hashing::Real6 & first_real6,
	numeric::geometry::hashing::Real6 & second_real6,
	core::Real dist_tolerance,
	core::Real euler_tolerance )
{

	if ( std::abs( first_real6[1] - second_real6[1] ) > dist_tolerance ||
			std::abs( first_real6[2] - second_real6[2] ) > dist_tolerance ||
			std::abs( first_real6[3] - second_real6[3] ) > dist_tolerance ||
			( std::abs( first_real6[4] - second_real6[4] ) > euler_tolerance && std::abs( first_real6[4] - second_real6[4] ) < (360.0 - euler_tolerance) ) ||
			( std::abs( first_real6[5] - second_real6[5] ) > euler_tolerance && std::abs( first_real6[5] - second_real6[5] ) < (360.0 - euler_tolerance) ) ||
			( std::abs( first_real6[6] - second_real6[6] ) > euler_tolerance && std::abs( first_real6[6] - second_real6[6] ) < (360.0 - euler_tolerance) )  ) {
		//std::abs( first_real6[4] - second_real6[4] ) > euler_tolerance ||
		//std::abs( first_real6[5] - second_real6[5] ) > euler_tolerance ||
		//std::abs( first_real6[6] - second_real6[6] ) > euler_tolerance ) {

		return false;
	} else {
		return true;
	}
}

/// @brief Main method
int main(int argc, char *argv[])
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		typedef numeric::geometry::hashing::Real6 Real6;
		typedef std::pair< std::pair< core::Size, core::Size >, Real6  > Key;
		typedef core::fragment::Frame Frame;
		typedef utility::vector1< std::pair<Key, Frame> >::iterator it_type;

		using namespace ObjexxFCL;

		devel::init(argc, argv);

		/*  Make cap_pose and bot_pose */
		TR << "Initializing LooDo." << std::endl;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////// Gathering User Input //////////

		// The ins_begin input variable denote the beginning of the insertion site (pose numbering) in the parent domain).
		int ins_begin = option[ OptionKeys::loodo::ins_begin ]();
		int ins_end = ins_begin+1; //Defaults to next residue for now.

		// Input cap pdb
		std::string cap = option[ OptionKeys::loodo::cap ]();
		core::pose::Pose cap_pose;
		core::import_pose::pose_from_file( cap_pose, cap , core::import_pose::PDB_file);

		// Input bottom pdb
		std::string bot = option[ OptionKeys::loodo::bot ]();
		core::pose::Pose bot_pose;
		core::import_pose::pose_from_file( bot_pose, bot , core::import_pose::PDB_file);

		// The fragA and fragB vectors contain lengths of loop fragments to insert.
		utility::vector1< int > fragAvector;
		if ( option[ OptionKeys::loodo::fragAlength ].user() ) {
			fragAvector = option[ OptionKeys::loodo::fragAlength ]();
		}
		utility::vector1< int > fragBvector;
		if ( option[ OptionKeys::loodo::fragBlength ].user() ) {
			fragBvector = option[ OptionKeys::loodo::fragBlength ]();
		}

		// Initializing Pose Data
		core::scoring::dssp::Dssp DSSP_cap(cap_pose);
		DSSP_cap.insert_ss_into_pose(cap_pose);
		std::string capss = cap_pose.secstruct();
		std::string capaa = cap_pose.sequence();

		core::scoring::dssp::Dssp DSSP_bot(bot_pose);
		DSSP_bot.insert_ss_into_pose(bot_pose);
		std::string botss = bot_pose.secstruct();
		std::string botaa = bot_pose.sequence();

		// Finding Center of Mass Residues
		core::Size botCOM = find_centerest_residue(bot_pose);
		core::Size capCOM = find_centerest_residue(cap_pose);

		// User Sanity Check
		TR << "Cap PDB: " << cap << std::endl;
		TR << "Number of residues: " << cap_pose.total_residue() << std::endl;
		TR << "Amino Acid sequence: " << cap_pose.sequence() << std::endl;
		TR << "Secondary structure sequence: " << cap_pose.secstruct() << std::endl;
		TR << "Center of Mass residue: " << capCOM << std::endl;
		if ( basic::options::option[ basic::options::OptionKeys::loodo::com_in_grid ]() ) {
			TR << "Requiring cap COM to lie in active site grid." << std::endl;
		}
		TR << "Ratio of CA atoms required to lie in active site grid: " << basic::options::option[ basic::options::OptionKeys::loodo::ca_ratio ]() << std::endl;

		TR << "Bottom PDB: " << bot << std::endl;
		TR << "Number of residues: " << bot_pose.total_residue() << std::endl;
		TR << "Amino Acid sequence: " << bot_pose.sequence() << std::endl;
		TR << "Secondary structure sequence: " << bot_pose.secstruct() << std::endl;
		TR << "Center of Mass residue: " << botCOM << std::endl;

		TR << "Insertion point begins at residue " << ins_begin << " [" << bot_pose.residue(ins_begin).name3() << "] and ends at residue " << ins_end << " [" << bot_pose.residue(ins_end).name3() << "]." << std::endl;

		bool loud = basic::options::option[ basic::options::OptionKeys::loodo::loud ]();
		core::Size clash_fail = 0;
		core::Size cont_fail = 0;
		if ( loud ) {
			TR << "Loud output." << std::endl;
		}

		core::Real dist_tolerance = basic::options::option[ basic::options::OptionKeys::loodo::distance_tolerance ]();
		core::Real euler_tolerance = basic::options::option[ basic::options::OptionKeys::loodo::euler_tolerance ]();
		TR << "Distance/Euler Angle tolerance for Real6 similarity: " << dist_tolerance << "/" << euler_tolerance << std::endl;

		// This block of code will warn you if you have specified a different length of fragments to be used compared to what is in the library you've specified.
		if ( basic::options::option[basic::options::OptionKeys::loodo::use_fraglib ].user() /* || basic::options::option[basic::options::OptionKeys::loodo::use_fraglibsc ].user() */  ) {

			std::string filename;
			//if ( basic::options::option[basic::options::OptionKeys::loodo::use_fraglib ].user() ){
			filename = basic::options::option[basic::options::OptionKeys::loodo::use_fraglib ]();
			//}
			//else { filename = basic::options::option[basic::options::OptionKeys::loodo::use_fraglibsc ](); }
			utility::io::izstream data( filename );

			if ( !data.good() ) {
				std::cerr << "Open failed for file: " << data.filename() << std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}

			std::string line;
			core::Size blank_counter = 0;
			core::Size frag_length = 0;
			while ( getline( data, line ) ) {
				if ( line == "" || line == " " ) {
					++blank_counter;
					if ( blank_counter == 1 ) {
						frag_length = -1;
					}
					if ( blank_counter == 2 ) {
						break;
					}
				}
				++frag_length;
			}
			for ( utility::vector1< int >::iterator fragAlength_iterator = fragAvector.begin(); fragAlength_iterator != fragAvector.end(); ++fragAlength_iterator ) {
				core::Size A_length = *fragAlength_iterator;
				if ( A_length+4 != frag_length ) {
					std::cerr << "\n" << std::endl;
					std::cerr << "WARNING!" << std::endl;
					std::cerr << "The first fragment in the library has " << frag_length << " residues, but -loodo:fragAlength specifies fragments of length " << A_length+4 << "." << std::endl;
					utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
				}
			}
			for ( utility::vector1< int >::iterator fragBlength_iterator = fragBvector.begin(); fragBlength_iterator != fragBvector.end(); ++fragBlength_iterator ) {
				core::Size B_length = *fragBlength_iterator;
				if ( B_length+4 != frag_length ) {
					std::cerr << "\n" << std::endl;
					std::cerr << "WARNING!" << std::endl;
					std::cerr << "The first fragment in the library has " << frag_length << " residues, but -loodo:fragBlength specifies fragments of length " << B_length+4 << "." << std::endl;
					utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
				}
			}
		}
		// End warning

		// Center the Native Cap for STUB calculations.
		core::pose::Pose centered_cap;
		core::import_pose::pose_from_file( centered_cap, cap , core::import_pose::PDB_file);
		centered_cap.center();
		//centered_cap.dump_pdb("Centered_Native.pdb");

		// Open Output Files
		std::ofstream caphitRTs;
		caphitRTs.open("CapHit_RT.txt");
		std::ofstream DebugRTs;
		DebugRTs.open("Debug_RT.txt");

		// Initialize Bump and AS Grids.
		TR << "Initializing BumpGrid." << std::endl;
		protocols::match::BumpGrid bb_grid = make_bb_grid(bot_pose);
		TR << "Initializing Placement Grid." << std::endl;
		protocols::match::downstream::ActiveSiteGrid active_site_grid = make_active_site_grid();
		TR << "Done." << std::endl;

		// Required by fragment_as_pose method:
		core::chemical::ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard");

		// Store the native cap CA coordinates for residues involved in STUB calculation
		numeric::xyzVector<core::Real> cap1_ca = cap_pose.residue(capCOM).atom("CA").xyz();
		numeric::xyzVector<core::Real> cap2_ca = cap_pose.residue(1).atom("CA").xyz();
		numeric::xyzVector<core::Real> cap3_ca = cap_pose.residue(cap_pose.n_residue()).atom("CA").xyz();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// end setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		// For every vector length of FragA:
		for ( utility::vector1< int >::iterator fragAlength_iterator = fragAvector.begin(); fragAlength_iterator != fragAvector.end(); ++fragAlength_iterator ) {

			std::stringstream VV;
			VV << *fragAlength_iterator;
			std::string fragAlength = VV.str();
			TR << "Fragment A Length: " << fragAlength << std::endl;

			// For every vector length of FragB:
			for ( utility::vector1< int >::iterator fragBlength_iterator = fragBvector.begin(); fragBlength_iterator != fragBvector.end(); ++fragBlength_iterator ) {

				std::stringstream VW;
				VW << *fragBlength_iterator;
				std::string fragBlength = VW.str();
				TR << "Fragment B Length: " << fragBlength << std::endl;

				utility::vector1< std::pair<Key, Frame> > library_A;

				core::fragment::ConstantLengthFragSetOP fraglibA;
				core::fragment::ConstantLengthFragSetOP fraglibB;

				// Round 1 = Linker A
				int round = 1;

				TR << "Making FragA Library." << std::endl;

				// If the user has specified a PDB corresponding to the native linkerA, then bias the fragment picker
				// to pick fragments closer to the native linker ss and aa.
				if ( option[ OptionKeys::loodo::known ].user() && option[ OptionKeys::loodo::fragAnative ].user() ) {
					int known = option[ OptionKeys::loodo::known ]();
					std::string fragA_nat = option[ OptionKeys::loodo::fragAnative ]();
					core::pose::Pose fragA_native;
					core::import_pose::pose_from_file( fragA_native, fragA_nat, core::import_pose::PDB_file);
					core::scoring::dssp::Dssp DSSP_fragA(fragA_native);
					DSSP_fragA.insert_ss_into_pose(fragA_native);
					std::string fragAnat_ss = fragA_native.secstruct();
					std::string fragAnat_aa = fragA_native.sequence();

					fraglibA = molten_get_vallfrags(*fragAlength_iterator, fragAnat_ss, fragAnat_aa, known);

				} else {
					fraglibA = get_vallfrags(round, *fragAlength_iterator, capss, capaa, botss, botaa, ins_begin, ins_end);
				}
				TR << "Made FragA library containing " << fraglibA->size() << " fragments." << std::endl;

				// For each fragment in library A, superimpose the fragment to the parent domain, then
				// superimpose the insert domain to the fragment. Check the insert domain for clashes and
				// containment within the placement grid.
				TR << "Beginning containment and clash check for each loop in FragA library." << std::endl;
				for ( core::fragment::ConstFrameIterator i = fraglibA->begin(); i != fraglibA->end(); ++i ) {
					//TR << "Getting position... " << std::endl;
					core::Size positionA = (*i)->start();
					//TR << "Getting frame..." << std::endl;
					core::fragment::Frame frameA = **i;
					//TR << "Frame has " << frameA.nr_frags() << " frags." << std::endl;
					for ( core::Size h = 1; h<=frameA.nr_frags(); h++ ) {
						//TR << "For each frag in this frame..." << std::endl;

						core::pose::Pose frag_pose;
						frameA.fragment_as_pose (h, frag_pose, rsd_set);

						// For visualizing translated caps and frags:
						std::stringstream HH;
						HH << h;
						std::string H = HH.str();

						std::stringstream II;
						II << positionA;
						std::string I = II.str();

						// Superimpose fragA to parent domain
						core::pose::Pose trans_frag_pose = perturb("A", "ftb", bot_pose, frag_pose, ins_begin, ins_end);
						// Superimpose insert domain to fragA
						core::pose::Pose trans_cap_pose = perturb("A", "ctf", trans_frag_pose, cap_pose, ins_begin, ins_end);

						if ( containment_and_clash_checks( bb_grid, active_site_grid, trans_cap_pose, capCOM, clash_fail, cont_fail /*trans_frag_pose*/) ) {

							// Get real6
							Real6 real6_A = six_coords_to_real6( bot_pose.residue(botCOM).atom("CA").xyz(),
								bot_pose.residue(ins_begin).atom("CA").xyz(),
								bot_pose.residue(ins_end).atom("CA").xyz(),
								trans_cap_pose.residue(capCOM).atom("CA").xyz(),
								trans_cap_pose.residue(1).atom("CA").xyz(),
								trans_cap_pose.residue(trans_cap_pose.n_residue()).atom("CA").xyz() );

							if ( option[ OptionKeys::loodo::debug ]() ) {
								DebugRTs << "CapHitA_fl"+H+"_"+I+" || "
									<< real6_A[4] << " " << real6_A[5] << " " << real6_A[6] << " "
									<< real6_A[1] << " " << real6_A[2] << " " << real6_A[3]
									<< std::endl;
							}
							if ( option[ OptionKeys::loodo::dump_all_As ]() ) {
								trans_cap_pose.dump_pdb("CapHitA_fl"+H+"_"+I+".pdb");
							}

							std::pair<core::Size, core::Size> PH;
							PH = std::make_pair(positionA, h);

							Key mykey;
							mykey = std::make_pair(PH, real6_A);

							Frame frame = **i;

							std::pair<Key,Frame> mypair = std::make_pair(mykey, frame);

							// Keep information for generating fragment in library_A.
							library_A.push_back(mypair);

						}//end: if containment_and_clash_checks
					}//end: for each fragment in frame
				}//end: for each frame

				TR << "FragA: Number of fragments passing containment/clash checks: " << library_A.size() << std::endl;
				if ( loud ) {
					TR << "Failed Clash Check: " << clash_fail << std::endl;
					TR << "Failed Containment Check: " << cont_fail << std::endl;
				}

				// Repeat above for FragB.
				// Round 2 = LinkerB
				round++; //round = 2

				TR << "Making FragB Library." << std::endl;

				if ( option[ OptionKeys::loodo::known ].user() && option[ OptionKeys::loodo::fragBnative ].user() ) {
					int known = option[ OptionKeys::loodo::known ]();
					std::string fragB_nat = option[ OptionKeys::loodo::fragBnative ]();
					core::pose::Pose fragB_native;
					core::import_pose::pose_from_file( fragB_native, fragB_nat, core::import_pose::PDB_file);
					core::scoring::dssp::Dssp DSSP_fragB( fragB_native );
					DSSP_fragB.insert_ss_into_pose( fragB_native );
					std::string fragBnat_ss = fragB_native.secstruct();
					std::string fragBnat_aa = fragB_native.sequence();

					fraglibB = molten_get_vallfrags(*fragBlength_iterator, fragBnat_ss, fragBnat_aa, known);
				} else {
					fraglibB = get_vallfrags(round, *fragBlength_iterator, capss, capaa, botss, botaa, ins_begin, ins_end);
				}

				TR << "Made FragB library containing " << fraglibB->size() << " fragments." << std::endl;

				TR << "Beginning containment and clash check for each loop in FragB library." << std::endl;
				core::Size b_count(0);
				cont_fail = 0;
				clash_fail = 0;
				for ( core::fragment::ConstFrameIterator j = fraglibB->begin(); j != fraglibB->end(); ++j ) {
					core::Size positionB = (*j)->start();
					core::fragment::Frame frameB= **j;

					for ( core::Size num_frag_B = 1; num_frag_B <= frameB.nr_frags(); num_frag_B++ ) {

						core::pose::Pose frag_pose;
						frameB.fragment_as_pose (num_frag_B, frag_pose, rsd_set);

						// For assigning unique names to linkerBs and CapHitBs
						std::stringstream KK;
						KK << num_frag_B;
						std::string linkerB_fragnum = KK.str();

						std::stringstream JJ;
						JJ << positionB;
						std::string linkerB_position = JJ.str();

						// Superimpose fragB to bot (ftb)
						core::pose::Pose trans_frag_pose = perturb("B", "ftb", bot_pose, frag_pose, ins_begin, ins_end);
						// Superimpose cap to fragB (ctf)
						core::pose::Pose trans_cap_poseB = perturb("B", "ctf", trans_frag_pose, cap_pose, ins_begin, ins_end);

						if ( containment_and_clash_checks( bb_grid, active_site_grid, trans_cap_poseB, capCOM, clash_fail, cont_fail /*trans_frag_pose*/) ) {
							++b_count;
							int LAnum = 0;

							// Calculate Real6 for CapHitB
							Real6 real6_B = six_coords_to_real6( bot_pose.residue(botCOM).atom("CA").xyz(),
								bot_pose.residue(ins_begin).atom("CA").xyz(),
								bot_pose.residue(ins_end).atom("CA").xyz(),
								trans_cap_poseB.residue(capCOM).atom("CA").xyz(),
								trans_cap_poseB.residue(1).atom("CA").xyz(),
								trans_cap_poseB.residue(trans_cap_poseB.n_residue()).atom("CA").xyz() );

							if ( option[ OptionKeys::loodo::debug ]() ) {
								DebugRTs << "CapHitB_fl" + linkerB_position + "_" + linkerB_fragnum + " || "
									<< real6_B[4] << " " << real6_B[5] << " " << real6_B[6] << " "
									<< real6_B[1] << " " << real6_B[2] << " " << real6_B[3]
									<< std::endl;
							}
							if ( option[ OptionKeys::loodo::dump_all_Bs ]() ) {
								trans_cap_poseB.dump_pdb("CapHitB_fl" + linkerB_position + "_" + linkerB_fragnum + ".pdb");
							}


							for ( it_type it = library_A.begin(); it != library_A.end(); it++ ) {

								// Get FragA information again.
								core::pose::Pose libpose;
								Key fraginfo = it->first;
								Real6 real6_A = fraginfo.second;

								core::fragment::Frame fragframe = it->second;
								fragframe.fragment_as_pose(fraginfo.first.second, libpose, rsd_set);

								std::stringstream HHH;
								HHH << fraginfo.first.second;
								std::string linkerA_fragnum = HHH.str();

								std::stringstream III;
								III << fraginfo.first.first;
								std::string linkerA_position = III.str();

								if ( real6_similarity( real6_A, real6_B, dist_tolerance, euler_tolerance) ) {
									LAnum++;
									std::stringstream L;
									L << LAnum;
									std::string linkerAnumber = L.str();


									//Get CapHitB Stub
									centered_cap.center();
									core::kinematics::Stub stub_b = protocols::sic_dock::getxform( trans_cap_poseB.residue( capCOM ), centered_cap.residue( capCOM) );

									//Print to Output File
									caphitRTs
										<< "CapHitB_fl" + fragAlength + fragBlength + "_" + linkerB_position + "_" + linkerB_fragnum + "_" + linkerAnumber + "_" + linkerA_position + "_" + linkerA_fragnum + ".pdb"
										<< " || "
										<< " REAL6 "
										<< real6_B[4] << " " << real6_B[5] << " " << real6_B[6] << " "
										<< real6_B[1] << " " << real6_B[2] << " " << real6_B[3]
										<< " || "
										<< stub_b
										<< std::endl;

									if ( option[ OptionKeys::loodo::debug ]() ) {
										trans_cap_poseB.dump_pdb( "CapHitB_fl" + fragAlength + fragBlength + "_" + linkerB_position + "_" + linkerB_fragnum + "_" + linkerAnumber + "_" + linkerA_position + "_" + linkerA_fragnum + ".pdb");
									}

									// Superimpose FragA to Bottom
									core::pose::Pose LinkerA = perturb("A", "ftb", bot_pose, libpose, ins_begin, ins_end);

									// Superimpose Cap to FragA

									core::pose::Pose trans_cap_poseA = perturb("A", "ctf", LinkerA, cap_pose, ins_begin, ins_end);
									if ( option[ OptionKeys::loodo::debug ]() ) {
										trans_cap_poseA.dump_pdb( "CapHitA_fl" + fragAlength + fragBlength + "_" + linkerB_position + "_" + linkerB_fragnum + "_" + linkerAnumber + "_" + linkerA_position + "_" + linkerA_fragnum + ".pdb");
									}

									// Get CapHitA Stub
									centered_cap.center();
									core::kinematics::Stub stub_a = protocols::sic_dock::getxform( trans_cap_poseA.residue( capCOM ), centered_cap.residue( capCOM ));

									//Print to Output File
									caphitRTs
										<< "CapHitA_fl" + fragAlength + fragBlength + "_" + linkerB_position + "_" + linkerB_fragnum + "_" + linkerAnumber + "_" + linkerA_position + "_" + linkerA_fragnum + ".pdb"
										<< " || "
										<< " REAL6 "
										<< real6_A[4] << " " << real6_A[5] << " " << real6_A[6] << " "
										<< real6_A[1] << " " << real6_A[2] << " " << real6_A[3]
										<< " || "
										<< stub_a
										<< std::endl;

									// Output the Linker poses.
									trans_frag_pose.dump_pdb( "LinkerB_fl" + fragAlength + fragBlength + "_" + linkerB_position + "_" + linkerB_fragnum + ".pdb" );
									LinkerA.dump_pdb( "LinkerA_fl" + fragAlength + fragBlength + "_" + linkerA_position + "_" + linkerA_fragnum + ".pdb" );

								} // end: for all similar & passing A frags
							} // end: for all passing A frags
						} // end: B Containment/clash check
					} // end: for all frags in frame
				} // end: for all frames in fragment library B

				TR << "FragB: Number of fragments passing containment/clash checks: " << b_count << std::endl;
				if ( loud ) {
					TR << "Failed Clash Check: " << clash_fail << std::endl;
					TR << "Failed Containment Check: " << cont_fail << std::endl;
				}
				TR << "FINISHED: FragA = "+ fragAlength+", FragB = "+ fragBlength << std::endl;

			} // end: for every fragment B length
		} // end: for every fragment A length

		// Close Output Files
		caphitRTs.close();
		DebugRTs.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
