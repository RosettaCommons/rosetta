// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file NoesyModule.cc
/// @brief main hook-up for the automatic NOESY assignment module
/// @details
///   handling of input-output options
///   class NoesyModule:
///       read input files
///       write assignments, constraints
///       run assignment
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/NoesyModule.hh>

// Package Headers
#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/PeakFileFormat.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/PeakAssignment.hh>
#include <protocols/noesy_assign/PeakAssignmentResidueMap.hh>
#include <protocols/noesy_assign/PeakFileFormat_Sparky.hh>
#include <protocols/noesy_assign/CovalentCompliance.hh>

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

// for switching residue type set to centroid
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

// Utility headers
#include <utility>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/prof.hh>

#include <ObjexxFCL/string.functions.hh>
//Auto Headers
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <protocols/noesy_assign/CrossPeakList.impl.hh>
#include <protocols/noesy_assign/NoesyModule.impl.hh>
#include <utility/vector1.hh>

OPT_2GRP_KEY( File, noesy, in, resonances )
OPT_2GRP_KEY( FileVector, noesy, in, peaks )
OPT_2GRP_KEY( FileVector, noesy, in, peak_resonance_pairs )
OPT_2GRP_KEY( Boolean, noesy, in, use_assignments )
OPT_2GRP_KEY( File, noesy, in, decoys )
OPT_2GRP_KEY( File, noesy, in, local_dist_table ) /// WARNING WARNING WARNING THIS CODE IS THREAD-UNSAFE!

OPT_2GRP_KEY( File, noesy, out, resonances )
OPT_2GRP_KEY( File, noesy, out, peaks )
OPT_2GRP_KEY( File, noesy, out, talos )
OPT_2GRP_KEY( Boolean, noesy, out, names )
OPT_2GRP_KEY( Boolean, noesy, out, separate_peak_files )
OPT_2GRP_KEY( Boolean, noesy, out, unambiguous )
OPT_2GRP_KEY( String, noesy, out, format )
OPT_2GRP_KEY( Real, noesy, out, minVC )
OPT_2GRP_KEY( Real, noesy, out, padding )

OPT_1GRP_KEY( Boolean, noesy, no_decoys )

// OPT_1GRP_KEY( Boolean, noesy, no_symm )
// OPT_1GRP_KEY( Boolean, noesy, no_cs )
// OPT_1GRP_KEY( Boolean, noesy, no_upper )
// OPT_1GRP_KEY( Boolean, noesy, no_remove_diagonal )
// OPT_1GRP_KEY( Boolean, noesy, no_calibrate )

bool protocols::noesy_assign::NoesyModule::options_registered_( false );


bool protocols::noesy_assign::NoesyModule::cmdline_options_activated() {
	using namespace basic::options;
	using namespace OptionKeys;
	return ( option[ noesy::in::resonances ].user() && option[ noesy::in::peaks ].user() )
		|| option[ noesy::in::peak_resonance_pairs ].user();
}

void protocols::noesy_assign::NoesyModule::register_options() {
	using namespace basic::options;
	using namespace OptionKeys;

	PeakAssignmentParameters::register_options();
	PeakFileFormat_xeasy::register_options();
	if ( options_registered_ ) return;
	options_registered_ = true;

	//....
	NEW_OPT( noesy::in::resonances, "file with assigned chemical shifts", "" );

	NEW_OPT( noesy::in::peaks, "file with noesy peaks", "" );
	NEW_OPT( noesy::in::peak_resonance_pairs, "pairs of files that belong together: cc.peaks cc.prot ilv.peaks ilv.prot", "" );
	NEW_OPT( noesy::in::use_assignments, "when reading peaks the already existing assignments are not ignored", false );
	NEW_OPT( noesy::in::decoys, "silent file with decoys used for 3D structural compatibility test", "" );
	NEW_OPT( noesy::in::local_dist_table, "file with distances determined for instances with FragToAtomDist ", "" );

	NEW_OPT( noesy::out::resonances, "the parsed resonances file with translated atom names etc.", "cs_out.dat" );
	NEW_OPT( noesy::out::peaks, "the parsed peaks file with assignments", "NOE_out.dat" );
	NEW_OPT( noesy::out::talos, "write the resonances also as talos file", "cs_prot.tab" );
	NEW_OPT( noesy::out::names, "write atom-names rather then resonance ID for assignments", true );
	NEW_OPT( noesy::out::separate_peak_files, "write peaks to independent files with out::peaks as prefix", false );
	NEW_OPT( noesy::out::unambiguous, "write only the assignment with hightest VC", false );
	NEW_OPT( noesy::out::minVC, "write only assignments that contribute more than X to the peak-volume", 0.0 );
	NEW_OPT( noesy::out::format, "write as xeasy or sparky file", "xeasy" );
	NEW_OPT( noesy::out::padding, "add padding of X Angstrom to final constraints", 0.0 );
	NEW_OPT( noesy::no_decoys, "check comp. with decoys", false );

}

static basic::Tracer tr( "protocols.noesy_assign.NoesyModule" );

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
//using namespace OptionKeys;
///// templates

namespace protocols {
namespace noesy_assign {

/// @details Auto-generated virtual destructor
NoesyModule::~NoesyModule() = default;


///Constructor   - read input files / requires options to be initialized
NoesyModule::NoesyModule( std::string  fasta_sequence ) :
	crosspeaks_( /* NULL */ ),
	sequence_(std::move( fasta_sequence ))
	//  main_resonances_( new ResonanceList( fasta_sequence ) )
{
	read_input_files();
	runtime_assert( options_registered_ );
	// moved this option into PeakAssignmentParameters.cc
	// skip_network_analysis_ = basic::options::option[  options::OptionKeys::noesy::no_network ]();
}


///delete all data and read input files again...  aka fresh_instance()
void NoesyModule::reset() {
	crosspeaks_ = nullptr;
	// main_resonances_ = new ResonanceList( main_resonances_->sequence() );
	read_input_files();
}

///read peak- and resonance files
void NoesyModule::read_input_files() {
	ResonanceListOP main_resonances( new ResonanceList( sequence() ) );
	basic::ProfileThis doit( basic::NOESY_ASSIGN_READ_INPUT );
	using namespace basic::options;
	//read resonances
	if ( option[ OptionKeys::noesy::in::resonances ].user() ) {
		utility::io::izstream input_file( option[ OptionKeys::noesy::in::resonances ]() );
		utility::io::ozstream output_file( option[ OptionKeys::noesy::out::resonances ]() );
		if ( input_file.good() ) {
			main_resonances->read_from_stream( input_file );
			main_resonances->write_to_stream( output_file );
			if ( option[ OptionKeys::noesy::out::talos ].user() ) {
				utility::io::ozstream talos_file( option[ OptionKeys::noesy::out::talos ]() );
				main_resonances->write_talos_format( talos_file, true /*backbone only*/ );
			}
		} else {
			tr.Error << "cannot read " << input_file << std::endl;
			throw utility::excn::EXCN_FileNotFound( option[ OptionKeys::noesy::in::resonances ]() );
		}
	}

	{ //scope
		//read peak lists
		crosspeaks_ = CrossPeakListOP( new CrossPeakList() );
		Size nfiles( option[ OptionKeys::noesy::in::peaks ]().size() );

		//loop-all files: read one file at a time
		for ( core::Size ifile = 1; ifile <= nfiles; ++ifile ) {  //as many as there are files
			std::string file( option[ OptionKeys::noesy::in::peaks ]()[ ifile ] );
			utility::io::izstream input_file( file );
			if ( input_file.good() ) {
				if ( main_resonances->size() < 1 ) {
					throw utility::excn::EXCN_BadInput( "attempt to read Peak-Files with option -noesy:in:peaks without a global resonance file: -noesy:in:resonances" );
				}
				PeakFileFormat_xeasy format;
				format.set_filename( option[ OptionKeys::noesy::in::peaks ]()[ ifile ].base() );
				format.set_ignore_assignments( !option[ OptionKeys::noesy::in::use_assignments ]() );
				tr.Info << "reading " << file << "... " << std::endl;
				crosspeaks_->read_from_stream( input_file, format, main_resonances );
			} else {
				tr.Error << "cannot read " << file << std::endl;
			}
		}
	} // scope end

	//specific resonances for a given peak-file
	if ( option[ OptionKeys::noesy::in::peak_resonance_pairs ].user() ) { //scope
		Size n_pair_files(  option[ OptionKeys::noesy::in::peak_resonance_pairs ]().size() );
		if ( n_pair_files % 2 != 0 ) {
			throw utility::excn::EXCN_BadInput( "odd number of entries in option -noesy:in:peak_resonance_pairs, always provide pairs of files  <*.peaks> <*.prot>" );
		}
		for ( core::Size ifile = 1; ifile <= n_pair_files; ifile += 2 ) {
			ResonanceListOP resonances( new ResonanceList( sequence() ) );
			utility::io::izstream res_input_file( option[ OptionKeys::noesy::in::peak_resonance_pairs ]()[ ifile+1 ] );
			if ( res_input_file.good() ) {
				resonances->read_from_stream( res_input_file );
				if ( tr.Debug.visible() ) {
					utility::io::ozstream output_file( "cs_"+ObjexxFCL::string_of( ceil(ifile/2. ))+".dat" );
					resonances->write_to_stream( output_file );
				}
			} else {
				tr.Error << "cannot read " << res_input_file << std::endl;
				throw utility::excn::EXCN_FileNotFound( option[ OptionKeys::noesy::in::resonances ]() );
			}
			std::string file( option[ OptionKeys::noesy::in::peak_resonance_pairs ]()[ ifile ] );
			utility::io::izstream input_file( file );
			if ( input_file.good() ) {
				PeakFileFormat_xeasy format;
				format.set_filename( option[ OptionKeys::noesy::in::peak_resonance_pairs ]()[ ifile ].base() );
				format.set_ignore_assignments( ! option[ OptionKeys::noesy::in::use_assignments ]() );
				tr.Info << "reading " << file << "... " << std::endl;
				crosspeaks_->read_from_stream( input_file, format, resonances );
			} else {
				tr.Error << "cannot read " << file << std::endl;
			}
		}
	} // scope end

	if ( option[ OptionKeys::noesy::in::local_dist_table ].user() ) {
		// WARNING WARNING WARNING! THREAD UNSAFE!
		CovalentCompliance::get_instance()->load_dist_table( option[ OptionKeys::noesy::in::local_dist_table ]() );
	}
}

void NoesyModule::add_dist_viol_to_assignments( core::pose::Pose native_pose) {
#ifndef WIN32
	using namespace core::scoring::constraints;
	for ( auto it = crosspeaks_->begin(); it != crosspeaks_->end(); ++it ) {
		if ( (*it)->eliminated() ) continue;
		if ( (*it)->min_seq_separation_residue_assignment( 0.1 ) < 1 ) continue;
		for ( auto ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			core::scoring::func::FuncOP func( new BoundFunc(1.5, (*it)->distance_bound(), 1, "NOE Peak " ) );
			ConstraintOP new_cst( (*ait)->create_constraint( native_pose, func ) );
			(*ait)->set_native_distance_viol( new_cst->score( native_pose ) );
		}
	}
#endif
}

/// @brief write peak assignments into peak-file (sparky, cyana)
void NoesyModule::write_assignments( std::string file_name ) {
	using namespace basic::options;

	ProfileThis doit( NOESY_ASSIGN_WRITE_ASSIGNMENTS );
	if ( file_name == "use_cmd_line" ) {
		file_name = option[ OptionKeys::noesy::out::peaks ]();
	}
	PeakFileFormatOP format;
	std::string const format_str( option[ OptionKeys::noesy::out::format ]() );
	if ( format_str == "xeasy" ) {
		format = PeakFileFormatOP( new PeakFileFormat_xeasy() );
	} else if ( format_str == "sparky" ) {
		format = PeakFileFormatOP( new PeakFileFormat_Sparky() );
	} else utility_exit_with_message( "NOE_data output format "+format_str+" is not known! ");
	format->set_write_atom_names( option[ OptionKeys::noesy::out::names ]() );
	format->set_write_only_highest_VC( option[ OptionKeys::noesy::out::unambiguous ]() );
	format->set_min_VC_to_write( option[ OptionKeys::noesy::out::minVC ]() );
	if ( option[ OptionKeys::noesy::out::separate_peak_files ]() ) {
		crosspeaks_->write_peak_files( file_name, *format );
	} else {
		utility::io::ozstream output_file( file_name );
		crosspeaks_->write_to_stream( output_file, *format );
	}
}

/// @brief assign peaks ( no explicit decoys - wrapper )
void NoesyModule::assign() {
	using namespace basic::options;
	using namespace OptionKeys;

	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData sfd( opts );
	std::string file_name("NO_FILE");
	if ( !option[ noesy::no_decoys ]() ) {
		if ( option[ noesy::in::decoys ].user() ) {
			file_name = option[ noesy::in::decoys ]();
			sfd.read_file( file_name );
		}
		if ( sfd.size() == 0 && option[ OptionKeys::in::file::silent ].user() ) {
			file_name = option[ OptionKeys::in::file::silent ]()[ 1 ];
			sfd.read_file( file_name );
		}
	}
	if ( sfd.begin() != sfd.end() ) {
		if ( sfd.begin()->sequence().one_letter_sequence() != sequence() ) {
			std::ostringstream str;
			str << "Sequence of input structures in file " << file_name << " does not match the sequence in resonance file and fasta-file\n";
			str << "Resonances: " << sequence() << "\n";
			str << "Ensemble: " << sfd.begin()->sequence().one_letter_sequence();
			utility_exit_with_message( str.str() );
		}
	}
	assign( sfd.begin(), sfd.end() );
}

/// @brief generate constraint files from assignments
void NoesyModule::generate_constraint_files(
	core::pose::Pose const& pose,
	std::string const& cst_fa_file,
	std::string const& cst_centroid_file,
	core::Size min_seq_separation,
	core::Size min_quali,
	core::Size max_quali,
	bool ignore_elimination_candidates, /*default = true*/
	bool elimination_candidates /*default = false */
) const {

	PROF_START( NOESY_ASSIGN_GEN_CST );

	using namespace core::scoring::constraints;
	using namespace basic::options;
	using namespace OptionKeys;

	core::pose::Pose centroid_pose = pose;
	core::util::switch_to_residue_type_set( centroid_pose, core::chemical::CENTROID_t );

	ConstraintSetOP cstset( new ConstraintSet );
	ConstraintSetOP centroid_cstset( new ConstraintSet );
	tr.Info << "generate constraints..." << std::endl;
	crosspeaks_->generate_fa_and_cen_constraints(
		cstset,
		centroid_cstset,
		pose,
		centroid_pose,
		min_seq_separation,
		min_quali,
		max_quali,
		option[ noesy::out::padding ](),
		ignore_elimination_candidates,
		elimination_candidates
	);

	PROF_STOP( NOESY_ASSIGN_GEN_CST );

	tr.Info << "write constraints..." << std::endl;
	PROF_START( NOESY_ASSIGN_WRITE_CST );

	ConstraintIO::write_constraints( cst_fa_file, *cstset, pose );
	ConstraintIO::write_constraints( cst_centroid_file, *centroid_cstset, centroid_pose );

	PROF_STOP( NOESY_ASSIGN_WRITE_CST );
}

}
}
