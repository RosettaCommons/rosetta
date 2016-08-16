// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  read and write XEASY format peak-lists
/// @author Oliver Lange
#include <protocols/noesy_assign/CrossPeakList.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/NoesyModule.hh>
#include <protocols/noesy_assign/PeakFileFormat.hh>
#include <utility/excn/Exceptions.hh>
//#include <devel/NoesyAssign/NoeNetwork.hh>
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>


#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

// for switching residue type set to centroid

#include <core/chemical/ChemicalManager.fwd.hh> //required for type-set-name FA_STANDARD


#include <basic/prof.hh>

#include <basic/Tracer.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>


static THREAD_LOCAL basic::Tracer tr( "main" );

using namespace core;
//using namespace protocols;

// OPT_KEY( FileVector, i )
// OPT_KEY( File, o )
// OPT_KEY( File, r )
// OPT_KEY( File, ro )
// OPT_KEY( File, simple_cst )
// //OPT_KEY( File, s )
OPT_2GRP_KEY( File, noesy, out, cst )
// OPT_KEY( String, assign )
// OPT_KEY( Boolean, ignore_assignments )
// OPT_KEY( Boolean, force_symmetry )
// OPT_KEY( Boolean, sequential )

// OPT_KEY( File, frags_for_ss )
// OPT_KEY( Real, hop_penalty )
// OPT_KEY( Integer, max_hops )
// OPT_KEY( Integer, prune_cutoff )
// OPT_KEY( Real, start_weight )
// OPT_KEY( Integer, edge_redundancy )
// OPT_KEY( Boolean, no_decoys )
// OPT_KEY( Boolean, no_network )
// OPT_KEY( Boolean, no_symm )
// OPT_KEY( Boolean, no_cs )
// OPT_KEY( Boolean, no_upper )
// OPT_KEY( Boolean, no_remove_diagonal )
// OPT_KEY( Boolean, no_calibrate )


// Had these "iterative" options in here to enable testing of NOE in script with same flags as archive run...

//OPT_1GRP_KEY( Boolean, iterative, assign_noes )
OPT_2GRP_KEY( Integer, noesy, out, min_seq_sep )
OPT_2GRP_KEY( Boolean, noesy, out, split )
OPT_2GRP_KEY( Integer, noesy, out, worst_prob_class )

//OPT_1GRP_KEY( Real, iterative, centroid_before_quickrelax_weight )
//OPT_1GRP_KEY( Real, iterative, fullatom_after_quickrelax_weight )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	protocols::noesy_assign::NoesyModule::register_options();
	//  Templates::register_options();
	// NEW_OPT( i, "read NOEs from this file", "NOE_in.dat" );
	//   NEW_OPT( o, "write NOEs to this file", "NOE_out.dat" );
	//   NEW_OPT( r, "read resonance assignments from file", "cs_in.dat" );
	//   NEW_OPT( ro, "write resonance assignments to file", "cs_out.dat" );
	//   NEW_OPT( simple_cst, "make simple constraints from assigned peaks", "NOE_simple.cst" );
	// OPT( in::file::s );
	//  OPT( in::file::silent );
	OPT( in::file::fasta );
	OPT( in::file::native );
	//  //NEW_OPT( s, "structure", "native.pdb" );
	NEW_OPT( noesy::out::cst, "write constraints", "assigned.cst" );
	NEW_OPT( noesy::out::min_seq_sep, "do not write constraints from peaks that have any residue pair less than sep. assigned", 2 );
	//  NEW_OPT( assign, "use an assignment strategy [simple, cyana]", "cyana" );
	//  NEW_OPT( ignore_assignments, "don't read existing assignments from peak-file", true );
	//  NEW_OPT( force_symmetry, "assign only symmetric peaks (makes sense, really)", false );
	//  NEW_OPT( sequential, "always assign cross-peak as sequential NOE if possible", false );
	//  NEW_OPT( frags_for_ss, "fragments that define SS structure for NOE Network analysis", "frags3.dat" );
	//  NEW_OPT( hop_penalty, "extra resistance for each hop", 2.0 );
	//  NEW_OPT( max_hops, "how many hops do we evaluate at most", 5);
	//  NEW_OPT( prune_cutoff, "if expected R_inv is at least X times larger, then total R, don't evaluated path", 10 );
	//  NEW_OPT( start_weight, "start with this weight", 0.1 );
	//  NEW_OPT( edge_redundancy, "need for direct NOEs to keep weight ", 4 );
	//  NEW_OPT( no_decoys, "check comp. with decoys", false );
	//  NEW_OPT( no_network, "check comp. with decoys", false );
	//  NEW_OPT( no_symm, "check comp. with decoys", false );
	//  NEW_OPT( no_cs, "check comp. with decoys", false );
	//  NEW_OPT( no_upper, "check upper", false );
	//  NEW_OPT( no_remove_diagonal, "", false );
	//  NEW_OPT( no_calibrate, "don't calibrate the NOE distance bound", false );
	// NEW_OPT( iterative::fullatom_after_quickrelax_weight, "[IGNORED] just to make cyana_test happy", 1.0);
	// NEW_OPT( iterative::centroid_before_quickrelax_weight, "[IGNORED] just to make cyana_test happy", 1.0);
	NEW_OPT( noesy::out::worst_prob_class, "write only restraints with class X and better [default: write all]", 5 );
	NEW_OPT( noesy::out::split, "write a separate file with the good (classes 0, 1, 2) restraints", false  );
}

void run_old() {
	using namespace protocols::noesy_assign;


	//  std::string fasta_sequence;
	//  pose::PoseOP native_pose = NULL;
	//  if ( basic::options::option[ basic::options::OptionKeys::in::file::fasta ].user() ) {
	//   fasta_sequence = core::sequence::read_fasta_file( basic::options::option[ basic::options::OptionKeys::in::file::fasta ]()[1] )[1]->sequence();
	//   tr.Info << "read fasta sequence: " << fasta_sequence.size() << " residues\n"  << fasta_sequence << std::endl;
	//  } else if ( native_pose ) {
	//   fasta_sequence = native_pose->sequence();
	//   tr.Info << "take sequence from native : " << fasta_sequence << std::endl;
	//  }
	//  ResonanceListOP resonances( new ResonanceList( fasta_sequence ) );
	//   { // read Resonances
	//     utility::io::izstream input_file(basic::options::option[ basic::options::OptionKeys::r ]() );
	//     utility::io::ozstream output_file( basic::options::option[ basic::options::OptionKeys::ro ]() );
	//     if ( input_file.good() ) {
	//       resonances->read_from_stream( input_file );
	//       resonances->write_to_stream( output_file );
	//     } else {
	//       tr.Error << "cannot read " << basic::options::option[ basic::options::OptionKeys::r ]() << std::endl;
	//     }
	//   }

	//  CrossPeakList cpl( resonances );
	//   { // read CrossPeaks
	//   Size nfiles( basic::options::option[ basic::options::OptionKeys::i ]().size() );
	//   for ( core::Size ifile = 1; ifile <= nfiles; ++ifile ) {
	//    std::string file( basic::options::option[ basic::options::OptionKeys::i ]()[ ifile ] );
	//    utility::io::izstream input_file( file );
	//    if ( input_file.good() ) {
	//     PeakFileFormat_xeasy format;
	//     format.set_ignore_assignments( basic::options::option[ basic::options::OptionKeys::ignore_assignments ]() );
	//     cpl.read_from_stream( input_file, format );
	//    } else {
	//     tr.Error << "cannot read " << file << std::endl;
	//    }
	//   }
	//   }


	//  if ( basic::options::option[ basic::options::OptionKeys::assign ].user() ) {
	//   std::string assign_strategy( basic::options::option[ basic::options::OptionKeys::assign ]() );
	//   if ( assign_strategy == "simple" ) {
	//    cpl.find_assignments();
	//   } else if ( assign_strategy == "cyana" ) {
	//    using namespace basic::options;
	//    using namespace basic::options::OptionKeys;
	//    cpl.find_assignments();
	//    if ( !option[ no_remove_diagonal ]() ) cpl.delete_diagonal_peaks();
	//    if ( !option[ no_cs ]() ) cpl.update_chemshiftscore();
	//    if ( !option[ no_symm ]() ) cpl.update_symmetry_score();
	//    if ( !option[ no_upper ]() ) cpl.update_upperdistance_score();
	//    core::io::silent::SilentFileData sfd;
	//    if ( !option[ no_decoys ]() && option[ in::file::silent ].user() ) {
	//     sfd.read_file( basic::options::option[ basic::options::OptionKeys::in::file::silent ]()[ 1 ] );
	//     cpl.update_decoy_compatibility_score( sfd.begin(), sfd.end() );
	//    } else {
	//     cpl.set_trivial_decoy_compatibility_score();
	//    }
	//    if ( !option[ no_network ]() ) cpl.network_analysis();
	//    cpl.update_peak_volumina();
	//    if ( !option[  no_calibrate ] () ) cpl.calibrate( sfd.begin(), sfd.end() );
	//    cpl.eliminate_spurious_peaks();
	//   } else {
	//    tr.Error << "assing strategy " << assign_strategy << " not recognized " << std::endl;
	//   }
	//  }

	//  if ( basic::options::option[ basic::options::OptionKeys::force_symmetry ]() ) {
	//   PeakAssignmentList assignments( resonances );
	//   assignments.add( cpl );
	//   assignments.check_for_symmetric_peaks( cpl );
	//  }

	//  if ( basic::options::option[ basic::options::OptionKeys::sequential ]() ) {
	//   PeakAssignmentList assignments( resonances );
	//   assignments.add( cpl );
	//   assignments.invalidate_competitors_to_sequential_NOE( cpl );
	//  }

	//  core::pose::Pose pose;
	//  core::import_pose::pose_from_file( pose, basic::core::options::option[ options::OptionKeys::in::file::s ]()[ 1 ] , core::import_pose::PDB_file);

	//  core::scoring::constraints::ConstraintSetOP cstset = cpl.generate_constraints( pose );
	//  core::scoring::constraints::ConstraintIO::write_constraints(  basic::options::option[ basic::options::OptionKeys::cst_out ](), *cstset, pose );

	//  core::scoring::constraints::ConstraintSetOP centroid_cstset = cpl.generate_constraints( pose, true );
	//  core::pose::Pose centroid_pose = pose;
	//  core::util::switch_to_residue_type_set( centroid_pose, core::chemical::CENTROID );
	//  core::scoring::constraints::ConstraintIO::write_constraints(  basic::options::option[ basic::options::OptionKeys::cst_out ]().name()
	//   + ".centroid", *centroid_cstset, centroid_pose );


	//  utility::io::ozstream output_file( basic::options::option[ basic::options::OptionKeys::o ]() );
	//  PeakFileFormat_xeasy format;
	//  cpl.write_to_stream( output_file, format );

	//  fragment::FragSetOP ss_frags = fragment::FragmentIO().read_data( basic::options::option[ basic::options::OptionKeys::frags_for_ss ]() );
	//  core::conformation::SecondaryStructure ss_def( *ss_frags, false /*no JustUseCentralResidue */ );

	//  utility::vector1< bool > beta_ss( ss_def.total_residue(), false );
	//  utility::vector1< bool > helix_ss( ss_def.total_residue(), false );
	//  utility::vector1< bool > any_ss( ss_def.total_residue(), false );
	//  for ( Size i=1; i<=ss_def.total_residue(); ++i ) {
	//   if ( ss_def.strand_fraction( i )>0.7 ) beta_ss[ i ]=true;
	//   if ( ss_def.helix_fraction( i )>0.7 ) helix_ss[ i ] = true;
	//   if ( helix_ss[ i ] ) tr.Debug << "helix " << i << std::endl;
	//   if ( beta_ss[ i ] ) tr.Debug << "beta: " << i << std::endl;
	//   if ( helix_ss[ i ] || beta_ss[ i ] ) any_ss[ i ]=true;
	//  }

	//  NoeNetwork network(
	//                  resonances,
	//                  basic::options::option[ basic::options::OptionKeys::max_hops ](),
	//                  basic::options::option[ basic::options::OptionKeys::hop_penalty ](),
	//                  basic::options::option[ basic::options::OptionKeys::prune_cutoff ]()
	//  );
	//  network.set_starting_weight( basic::options::option[ basic::options::OptionKeys::start_weight ]() );
	//  network.set_edge_redundancy( basic::options::option[ basic::options::OptionKeys::edge_redundancy ]() );
	//  network.add( cpl );
	//  network.add( ss_def );
	//  Real delta( 1.0 );
	//  while( delta > 10 ) {
	//   tr.Debug << network << std::endl;
	//   delta = network.evaluate_edges();
	//  }
	//  tr.Debug << network << std::endl;
	//pose.dump_pdb("the_reference_pose.pdb" );
	/*
	utility::io::ozstream cst_output_file( basic::options::option[ basic::options::OptionKeys::cst_out ]() );
	//  core::scoring::constraints::ConstraintSetOP cstset( new core::scoring::constraints::ConstraintSet );
	{//make simple constraint list from assigned peaks
	using namespace core::scoring::constraints;
	Size ct( 1 );
	for ( CrossPeakList::CrossPeaks::const_iterator it = cpl.peaks().begin(); it != cpl.peaks().end(); ++it, ++ct ) {
	tr.Debug << " find assignments for peak " << ct << std::endl;
	//   (*it)->find_assignments( resonances );
	if ( (*it)->assigned() && !(*it)->ambiguous() ) {
	PeakAssignment const& assignment( **((*it)->assignments().begin()));
	id::NamedAtomID const& atom1( assignment.atom( *resonances, 1 ) ); //[ (*it)->proton( 1 ).assignment( ind_assigned ) ].atom() );
	id::NamedAtomID const& atom2( assignment.atom( *resonances, 2 ) ); //[ (*it)->proton( 2 ).assignment( ind_assigned ) ].atom() );
	if ( atom1.rsd() != atom2.rsd() ) {
	//   if ( beta_ss[ atom1.rsd() ] && beta_ss[ atom2.rsd() ]
	//       && atom1.rsd() > 1 && atom1.rsd() < ss_def.total_residue() && ( beta_ss[ atom1.rsd()-1 ] || beta_ss[ atom1.rsd()+1] )
	//       && atom2.rsd() > 1 && atom2.rsd() < ss_def.total_residue() && ( beta_ss[ atom2.rsd()-1 ] || beta_ss[ atom2.rsd()+1] )
	//      ) {
	cst_output_file << "AmbiguousNMRDistance " << atom1 << " " << atom2 << " BOUNDED 1.5 6 0.5 NOE ; unambiguous NOE " << ct << std::endl;     //}
	}
	tr.Debug << " write assignments as constraints for peak " << ct << std::endl;
	// cstset->add_constraint( new AmbiguousNMRDistanceConstraint( atom1, atom2, pose, new BoundFunc( 1.5, 5.5, 0.5, "NOE" ) ) );
	} else if ( (*it)->assigned() ) {
	cst_output_file << "AmbiguousConstraint" << std::endl;
	//    Size ind_assigned( 0 );
	for ( CrossPeak::PeakAssignments::const_iterator ait=(*it)->assignments().begin(); ait!=(*it)->assignments().end(); ++ait ) {
	PeakAssignment const& assignment( **ait );
	id::NamedAtomID const& atom1( assignment.atom( *resonances, 1 ) ); //[ (*it)->proton( 1 ).assignment( ind_assigned ) ].atom() );
	id::NamedAtomID const& atom2( assignment.atom( *resonances, 2 ) ); //[ (*it)->proton( 2 ).assignment( ind_assigned ) ].atom() );
	if ( atom1.rsd() != atom2.rsd() ) {
	//     if ( beta_ss[ atom1.rsd() ] && beta_ss[ atom2.rsd() ]
	//       && atom1.rsd() > 1 && atom1.rsd() < ss_def.total_residue() && ( beta_ss[ atom1.rsd()-1 ] || beta_ss[ atom1.rsd()+1] )
	//       && atom2.rsd() > 1 && atom2.rsd() < ss_def.total_residue() && ( beta_ss[ atom2.rsd()-1 ] || beta_ss[ atom2.rsd()+1] )
	// ) {
	cst_output_file << "AmbiguousNMRDistance " << atom1 << " " << atom2 << " BOUNDED 1.5 6 0.5 NOE ; ambiguous NOE " << ct << std::endl;
	}
	}
	cst_output_file << "END_Ambiguous" << std::endl;
	}
	}
	}
	*/
}

void run() {
	using namespace protocols::noesy_assign;
	using namespace basic::options;

	std::string fasta_sequence;
	if ( option[ OptionKeys::in::file::fasta ].user() ) {
		fasta_sequence = core::sequence::read_fasta_file( option[ OptionKeys::in::file::fasta ]()[1] )[1]->sequence();
		tr.Info << "read fasta sequence: " << fasta_sequence.size() << " residues\n"  << fasta_sequence << std::endl;
	}

	PROF_START( basic::NOESY_ASSIGN_TOTAL );

	NoesyModule nm( fasta_sequence );

	nm.assign();

	core::pose::Pose pose;
	//  core::import_pose::pose_from_file( pose, option[ options::OptionKeys::in::file::s ]()[ 1 ] , core::import_pose::PDB_file);

	core::pose::make_pose_from_sequence(
		pose,
		fasta_sequence,
		chemical::FA_STANDARD
	);

	std::string cst_file( option[ OptionKeys::noesy::out::cst ]() );
	std::string cst_centroid_file( cst_file + ".centroid");
	std::string best_cst_file( cst_file + ".good");
	std::string best_cst_centroid_file( best_cst_file + ".centroid");

	if ( option[ OptionKeys::noesy::out::split ] ) {
		nm.generate_constraint_files( pose, best_cst_file, best_cst_centroid_file,
			option[ OptionKeys::noesy::out::min_seq_sep ](), 0, 2 );
		nm.generate_constraint_files( pose, cst_file, cst_centroid_file,
			option[ OptionKeys::noesy::out::min_seq_sep ](), 3, option[ OptionKeys::noesy::out::worst_prob_class ]() );
	} else {
		nm.generate_constraint_files( pose, cst_file, cst_centroid_file,
			option[ OptionKeys::noesy::out::min_seq_sep ](), 0, option[ OptionKeys::noesy::out::worst_prob_class ]() );
	}


	if ( option[ OptionKeys::in::file::native ].user() ) {
		core::pose::Pose native_pose;
		core::import_pose::pose_from_file( native_pose, option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		nm.add_dist_viol_to_assignments(native_pose);
	}

	nm.write_assignments();

	PROF_STOP( basic::NOESY_ASSIGN_TOTAL );
	basic::prof_show();
}

int
main( int argc, char * argv [] )
{
	try{
		register_options();
		protocols::noesy_assign::PeakAssignmentParameters::register_options();
		protocols::noesy_assign::PeakFileFormat_xeasy::register_options();
		devel::init( argc, argv );
		protocols::noesy_assign::PeakAssignmentParameters::get_instance();
		try{
			run();
		} catch ( utility::excn::EXCN_Base& anExcn ) {
			anExcn.show( std::cerr );
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
