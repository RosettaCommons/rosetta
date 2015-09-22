// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file r_trjconv.cc
/// @brief tool to handle and convert pdb/silent rosetta structure output
/// @author Oliver Lange
// libRosetta headers


//#include <protocols/jobdist/JobDistributors.hh>
//#include <protocols/jobdist/Jobs.hh>

#include <core/types.hh>
#include <devel/init.hh>

#include <core/conformation/Conformation.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/kinematics/FoldTree.hh>

#include <protocols/abinitio/ClassicAbinitio.hh>

#include <protocols/evaluation/PoseEvaluator.hh>


#include <core/pose/Pose.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>


#include <basic/options/option.hh>
#include <utility/excn/Exceptions.hh>
//#include <basic/options/OptionKeys.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
// C++ headers
#include <iostream>
#include <string>

#include <basic/Tracer.hh>

#include <core/kinematics/Stub.hh>
#include <core/id/types.hh>


static THREAD_LOCAL basic::Tracer tracer( "r_trjconv" );

using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;


#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>


class ThisApplication {
public:
	static void register_options();
	// void run();
	void pick_frags();
	// void pick_chunks();
};

using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( Boolean, check_simple_foldtree )
OPT_KEY( File, pick_frags )
OPT_KEY( File, pick_chunks )
OPT_KEY( Integer, size )
OPT_KEY( FileVector, rmsd_target )
OPT_KEY( StringVector, rmsd_column_name )
OPT_KEY( Boolean, dump_struct )
//OPT_1GRP_KEY( File, evaluate, cs_file )
void ThisApplication::register_options() {
	OPT( in::file::silent ); // input silent file
	OPT( out::file::silent );
	OPT( in::file::psipred_ss2 );
	NEW_OPT( check_simple_foldtree, "read structures change foldtree to simple and compute RMSD between complex/simple foldtree",false );
	NEW_OPT( pick_frags, "pick fragments from each decoy", "" );
	NEW_OPT( pick_chunks, "pick non-local fragments from each decoy", "" );
	NEW_OPT( size, "length of fragment to pick", 9 );
	NEW_OPT( rmsd_target,"[vector] determine rmsd against this structure ","" );
	NEW_OPT( rmsd_column_name, "[vector] use xxx as column name: rms_xxx ","");
	NEW_OPT( dump_struct,"write also structures to silent:out",false);
	// NEW_OPT( evaluate::cs_file, "compute a SPARTA chemical shift score for each model", false);
	//NEW_OPT( tag_selector, "a list of tags to be extracted","");
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// why was this function duplicated here? another copy lives in core/chemical/util.cc/hh, which seems like a better
/// place for it.
// void
// protocols::toolbox::switch_to_residue_type_set(
//  pose::Pose & pose,
//  std::string const & type_set_name
// )
// {
//  protocols::toolbox::switch_to_residue_type_set( pose, type_set_name );
// }

// void compute_chi( pose::Pose &pose ) {
//  //this little method demonstrates access to CHI-angles, atom-coordinates, Atom-Tree, Stubs
//  // what we are going to do is to compute the position of SG on CYS outside of the pose
//  // we use the same mechanism as deep inside the atomtree, i.e., we get the stub on CB and transform this stub to get SG.
//  // we copy the internal coordinates d, theta and the chi-angle from the pose and then compare with the pose-computed position of SG.
//  using namespace core;
//  using namespace id;
//  using namespace numeric;
//  using namespace kinematics;
//  // using namespace conformation;
//  // ATTENTION no error checking in this code. If you have an invalid AtomID you are screwed.

//  core::Size const seqpos( 286 );//a CYS
//  id::NamedAtomID sulfur_id ("SG", seqpos );
//  core::kinematics::tree::Atom const& sulfur( pose.conformation().atom_tree().atom( id::AtomID( sulfur_id, pose ) ) );

//  id::DOF_Type dof_id_theta( THETA ); //the angle SG-CB-CA
//  id::DOF_Type dof_id_d( D ); //the distance SG-CB

//  core::Real const theta( //getting the angle of SG
//           sulfur.dof( dof_id_theta )
//  );

//  //three different ways of getting the S - CB distance.
//  core::Real const S_CB_dist( sulfur.distance ( *sulfur.parent() )); //I know CB is parent in Atom-Tree ( BAD BEHAVIOUR )
//  core::Real const S_CB_dist2( sulfur.dof( dof_id_d ) ); // Thats the better way of doing it. BUT STILL BAD !
//  core::Real const S_CB_dist3( distance( pose.xyz( sulfur_id ), pose.xyz( NamedAtomID( "CB", seqpos ) )) ); //GOOD BEHAVIOUR

//  tracer.Info << "theta " << theta << " dist1 " << S_CB_dist << " dist2 " << S_CB_dist2 << " " << S_CB_dist3 << std::endl;

//  //the STUB of CB is where we start
//  id::NamedStubID CB_stub_id( "CB", "CA", "N", seqpos ); //sequence center-atom, parent, grand-parent
//  kinematics::Stub manual_stub( pose.stub_from_id( CB_stub_id ) ); //GOOD WAY OF GETTING IT
//  //just to demonstrate how the atom-tree works the BAD WAY OF GETTING THE STUB
//  kinematics::Stub atree_stub( //get the same stub directly from the AtmoTree --- asuming N2C folding direction
//                pose.conformation().atom_tree().atom(
//                                   id::AtomID( id::NamedAtomID( "CB", seqpos ), pose )
//                ).get_stub() );

//  tracer.Info << "\natree Stub " << atree_stub << "\n" << "manual     " << manual_stub << std::endl;
//  //get the current chi-angle that determines the position of SG ( chi1 ).
//  TorsionID chi1_id( seqpos, CHI, 1);
//  core::Real const chi( pose.torsion( chi1_id ) );
//  tracer.Info << "chi1 " << chi << std::endl;

//  //now compute the new STUB for SG. the center of the new stub is the position of SG.
//  manual_stub.M *= x_rotation_matrix_degrees( chi );
//  Stub new_stub( manual_stub.M * z_rotation_matrix_radians( theta ), manual_stub.v );
//  new_stub.v += S_CB_dist * new_stub.M.col_x(); //new center of stub -- should coincide with sulfur

//  Vector sulfur_xyz( pose.xyz( sulfur_id ) );
//  tracer.Info << " compare to pose.xyz " << distance(sulfur_xyz,new_stub.v) << std::endl;

// }

void ThisApplication::pick_frags()
{
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace io::silent;
	using namespace pose;
	using namespace fragment;
	SilentFileData sfd;
	sfd.read_file( *(option[ in::file::silent ]().begin()) );
	ConstantLengthFragSet frags( option[ OptionKeys::size ] ); //steal good old 9mers
	for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
		Pose pose;
		std::string tag = it->decoy_tag();
		it->fill_pose( pose );
		steal_constant_length_frag_set_from_pose( pose, frags );
	}
	FragmentIO().write_data( option[ OptionKeys::pick_frags ](), frags );
}


// void generate_frames_based_on_ss(
//  FrameList& frames,
//  core::conformation::SecondaryStructure const& ss
// ) {
//  using namespace loops;
//  //Loops ss_region = ss.compute_ss_regions( 0.3, 2 );
//  Loops ss_region = compute_ss_regions( 0.3, 2, ss );

//  for ( Loops::const_iterator it1 = ss_region.begin(); it1 != ss_region.end(); ++it1 ) {
//   for ( Loops::const_iterator it2 = it1 + 1 ; it2 != ss_region.end(); ++it2 ) {
//    JumpingFrameOP frame = new JumpingFrame( it1->start(), it2->stop(), it1->size() + it2->size() + 2 /*for jumps*/);
//    FragDataOP frag_data = new FragData;
//    Size pos = 1;
//    //set torsion residues
//    for ( Size i = it1->start(); i<= it1->stop(); i++ ) {
//     frame->set_pos( pos++, i );
//     frag_data->add_residue( new BBTorsionSRFD( 3, 'E', 'X' ) );
//    }
//    //now set two jump residues
//    Size up( it1->start() + it1->size()/2 );
//    Size down( it2->start() + it2->size()/2 );
//    tracer.Trace << "setting up FRAME for " << *it1  << " to " << *it2 << " jumping via " << up << " " << down << std::endl;
//    frame->set_pos( pos++, up);
//    frame->set_pos( pos++, down);
//    frag_data->add_residue( new UpJumpSRFD );
//    frag_data->add_residue( new DownJumpSRFD );
//    for ( Size i = it2->start(); i<= it2->stop(); i++ ) {
//     frame->set_pos( pos++, i );
//     frag_data->add_residue( new BBTorsionSRFD( 3, 'E', 'X' ) );
//    }
//    frame->add_fragment( frag_data );
//    frames.push_back( frame );
//    tracer.Trace << *frame << std::endl;
//   } //it2
//  } //it1
// }

//make sure that no cut-points within frame region
void checked_steal_fragment( Frame& frame, pose::Pose const& pose ) {
	for ( Size i = 1; i<=frame.length(); ++i ) {
		Size pos( frame.seqpos( i ) );
		//don't want to use fragments that go across a cutpoint... certainly non-ideal geometry there
		if ( i < frame.length() && ( pos + 1 == frame.seqpos( i+1 ) ) //consecutive
				&& pose.fold_tree().is_cutpoint( pos )  ) return;
	}
	//okay all consecutive parts in fragment are also consecutive in pose
	frame.steal( pose );
}

// void ThisApplication::pick_chunks()
// {
//  using namespace core;
//  using namespace basic::options;
//  using namespace basic::options::OptionKeys;
//  using namespace io::silent;
//  using namespace pose;
//  using namespace fragment;

//  core::conformation::SecondaryStructure ss_def;
//  ss_def.read_psipred_ss2( option[ in::file::psipred_ss2 ] );

//  FrameList frames;
//  generate_frames_based_on_ss( frames, ss_def );

//  SilentFileData sfd;
//  sfd.read_file( *(option[ in::file::silent ]().begin()) );

//  for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
//   Pose pose;
//   std::string tag = it->decoy_tag();
//   it->fill_pose( pose );
//   for ( FrameList::iterator frame_it = frames.begin(); frame_it != frames.end(); ++frame_it ) {
//    checked_steal_fragment( **frame_it, pose );
//   }
//  }
//  OrderedFragSet frags;
//  frags.add( frames );
//  FragmentIO().write_data( option[ OptionKeys::pick_chunks ](), frags );
// }


// void ThisApplication::run()
// {
//  using namespace core;
//  using namespace basic::options;

//  using namespace io::silent;
//  using namespace pose;

//  // a bunch of PoseEvaluators for process_decoy() --- if available
//  protocols::evaluation::MetaPoseEvaluator evaluator;


//  SilentFileData sfd;
//  sfd.read_file( *(option[ in::file::silent ]().begin()) );
//  SilentFileData sfd_out;

//  if ( option[ rmsd_target ].user() ) {
//   utility::vector1< std::string > const& rmsd_target( option[ OptionKeys::rmsd_target ]() );
//   utility::vector1< std::string > const& rmsd_col_name( option[ OptionKeys::rmsd_column_name ]() );
//   if ( rmsd_target.size() > 1 && rmsd_col_name.size() != rmsd_target.size() ) {
//    utility_exit_with_message("if you specify multiple rmsd_targets you have to specify the same number of column names");
//   }
//   for ( Size ct = 1; ct <= rmsd_target.size(); ct ++ ) {
//    pose::PoseOP rmsd_pose = new pose::Pose;
//    core::import_pose::pose_from_pdb( *rmsd_pose, rmsd_target[ ct ] );
//    std::string tag("");
//    if ( rmsd_col_name.size() >= ct ) tag = rmsd_col_name[ ct ];
//    evaluator.add_evaluation( new protocols::simple_filters::SelectRmsdEvaluator( rmsd_pose, tag ) );
//   }
//  }
//  if ( option[ OptionKeys::evaluation::chemical_shifts ].user() ) {
//   evaluator.add_evaluation( new protocols::simple_filters::ChemicalShiftEvaluator( "cs_score", option[  OptionKeys::evaluation::chemical_shifts ]() ) );
//  }


//  using namespace basic::options::OptionKeys;
//  Size ct( 0 );
//  for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
//   Pose pose;
//   std::string tag = it->decoy_tag();
//   if ( 1 ) { //  if ( tag == opt.tag_selected || opt.tag_selected.size()==0 ) {
//    tracer.Info << ++ct << " " <<  tag << std::endl;
//    it->fill_pose( pose );

//    //  it->second->print_conformation( std::cout );
//    if ( option [ out::file::silent ].user() ) {
//     // run PoseEvaluators
//     SilentStructOP ss = *it;
//     ProteinSilentStructOP pss= dynamic_cast< ProteinSilentStruct* > ( ss.get() );
//     evaluator.apply( pose, tag, *pss );
//     sfd_out.add_structure ( *it );

//    } else if ( option[ check_simple_foldtree ] ) {
//     it->fill_pose( pose );
//     pose.dump_pdb( tag + ".pdb");
//     pose::Pose pose_sf;
//     core::pose::make_pose_from_sequence(
//        pose_sf,
//        pose.sequence(),
//        *( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ))
//     );
//     // steal torsions
//     // determine length of segment to copy from native
//     Size seg_len = pose.total_residue();
//     fragment::Frame long_frame(1, seg_len);
//     fragment::FragData frag( new fragment::BBTorsionSRFD, seg_len);

//     // get torsion angles from complex fold-tree pose
//     frag.steal( pose, long_frame);
//     // apply torsions to extended structue
//     tracer.Debug << "stolen torsions: " << frag << std::endl;
//     frag.apply( pose_sf, long_frame);
//     pose_sf.dump_pdb( tag + ".sf.pdb");

//     core::Real rmsd = core::scoring::CA_rmsd( pose_sf, pose );
//     tracer.Info << tag << " rmsd for simple fold-tree " << rmsd << " chainbreak:" << it->get_energy("chainbreak") << std::endl;

//    } else {
//     //    pose.dump_pdb( tag + ".1.pdb");
//     it->fill_pose( pose );
//     pose.dump_pdb( tag + ".pdb");
//    }

//   }
//  }
//  if ( option[ out::file::silent ].user() ) {
//   sfd_out.write_all( option[ out::file::silent], !option[ dump_struct ] );
//  }
// }

int
main( int argc, char * argv [] )
{
	try{
		ThisApplication::register_options();
		devel::init( argc, argv );
		try {
			//   if ( option[ OptionKeys::in::file::s ].user() ) {
			//    core::pose::Pose pose;
			//    core::import_pose::pose_from_pdb( pose, *core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
			//     option[ OptionKeys::in::file::s ]()[ 1 ] );
			//    compute_chi( pose );
			//   }

			ThisApplication app;
			if ( option[ pick_frags ].user() ) {
				app.pick_frags();
			}
			//  else if ( option[ pick_chunks ].user() ) {
			//   app.pick_chunks();
			// } else app.run();
		} catch ( utility::excn::EXCN_Base& excn ) {
			std::cerr << "Exception : " << std::endl;
			excn.show( std::cerr );
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

