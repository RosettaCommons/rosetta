// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/devel/denovo_design/filters/PsiPredInterface.cc
/// @brief  Interface to PsiPred for secondary structure prediction
/// @author Tom Linsky

// unit headers
#include <devel/denovo_design/filters/PsiPredInterface.hh>

// package headers

// project headers
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

// boost
#include <boost/lexical_cast.hpp>

// C++ headers

static basic::Tracer TR("devel.denovo_design.filters.psipredinterface");

namespace devel {
namespace denovo_design {
namespace filters {

/// @brief default constructor
PsiPredInterface::PsiPredInterface( std::string const & cmd )
	: ReferenceCount(),
		cmd_( cmd )
{}

/// @brief copy constructor
PsiPredInterface::PsiPredInterface( PsiPredInterface const & rval )
	: ReferenceCount(),
		psipred_output_( rval.psipred_output_ ),
		cmd_( rval.cmd_ )
{}

/// @brief given a pose, generates a fasta string
std::string
PsiPredInterface::convert_to_fasta( core::pose::Pose const & pose ) const{
	std::string fasta_str( ">" + pose.pdb_info()->name() + "\n" );
	for ( core::Size i=1; i<=pose.total_residue(); i++ ){
		if ( pose.residue( i ).is_protein() ){
			fasta_str += pose.residue( i ).name1();
			if ( i % 80 == 0 ){
				fasta_str += "\n";
			} //if 80 chars
		} //if is protein
	} // for each residue
	return fasta_str;
}

/// @brief creates a fasta file
/// returns filename if successful
/// exits rosetta if not successful
std::string
PsiPredInterface::create_fasta_file( core::pose::Pose const & pose ) const{
	std::string fasta_filename( pose.pdb_info()->name() + ".fasta" );
	// get rid of the path and just use the current directory
	if ( fasta_filename.find('/') != std::string::npos ) {
		fasta_filename = fasta_filename.substr( fasta_filename.find_last_of( '/' )+1, std::string::npos );
	}
	std::string fasta( convert_to_fasta( pose ) );
	TR.Debug << "Fasta: " << fasta << std::endl;
	TR << "Fasta filename: " << fasta_filename << std::endl;

	// Write out a fasta file
	std::ofstream fastafile( fasta_filename.c_str() );
	if ( !fastafile.is_open() ) {
		utility_exit_with_message( "Error: The psipred filter could not open " + fasta_filename + " for writing." );
	}
	fastafile << fasta << std::endl;
	if ( fastafile.bad() ) {
		utility_exit_with_message( "The psipred filter encountered an error writing to " + fasta_filename );
	}
	fastafile.close();
	return fasta_filename;
}

/// @brief Parses the psipred pass2 horiz file and returns predicted secondary structure and confidence for each position.
std::pair< std::string, utility::vector1< core::Size > >
PsiPredInterface::parse_psipred_horiz_output( std::string const & psipred_horiz_filename ) const
{
	// open and read psipred output into a stringstream buffer
	utility::io::izstream t( psipred_horiz_filename );
	std::stringstream buffer;
	buffer << t.rdbuf();
	t.close();
	std::string line;
	std::string ss;
	utility::vector1< core::Size > conf;
	while ( getline( buffer, line ) ) {
		// skip comments and blank lines
		if ( line[0] == '#' || line == "" )
			continue;
		// we will basically read a block at a time in this while loop
		// first 6 characters are header, strip them
		std::string const data( line.substr( 6, std::string::npos ) );
		for ( core::Size i=0; i<data.size(); ++i ) {
			core::Size const num( data[i] - '0' );
			if ( num > 9 ) {
				utility_exit_with_message( std::string( "Invalid psipred integer conversion!, char=" ) + data[i] + ", line=" + line );
			}
			conf.push_back( num );
		}
		if ( ! getline( buffer, line ) ) {
			TR.Error << "Error parsing psipred output -- last line=" + line << std::endl;
		}
		std::string data2( line.substr( 6, std::string::npos ) );
		// convert C to L
		for ( core::Size i=0; i<data2.size(); ++i ) {
			if ( data2[i] == 'C' ) {
				data2[i] = 'L';
			}
		}
		ss += data2;
		//ignore the next two lines -- they are for Amino acid sequence, and position number, respectively
		for ( core::Size i=1; i<=2; ++i ) {
			if ( ! getline( buffer, line ) ) {
				TR.Error << "Error parsing psipred output -- last line=" + line << std::endl;
			}
		}
	}
	TR.Debug << "ss=" << ss << std::endl;
	return std::pair< std::string, utility::vector1< core::Size > >( ss, conf );
}

/// @brief Parses the psipred output and returns the predicted secondary structure and likelihoods of the blueprint secondary structure being present on the pose at each position.
PsiPredResult
PsiPredInterface::parse_psipred_output( core::pose::Pose const & /*pose*/,
		std::string const & psipred_str,
		std::string const & blueprint_ss,
		std::string const & psipred_horiz_filename ) const {
	std::istringstream psipred( psipred_str );
	utility::vector1< core::Real > probabilities;
	std::string pred_ss( "" );
	std::string line;
	core::Size count( 0 );
	while ( getline( psipred, line ) ){
		if( line[0] == '#' || line == "" )
			continue;
		std::istringstream line_stream( line );
		core::Size resi( 0 );
		char resn( 'X' ), ss( 'Y' );
		core::Real loop( 0.0 ), helix( 0.0 ), sheet( 0.0 );
		line_stream >> resi >> resn >> ss >> loop >> helix >> sheet;
		//TR << "Read line, resi=" << resi << " resn=" << resn << " ss=" << ss << " loop=" << loop << " helix=" << helix << " casheet=" << sheet << std::endl;

		if ( blueprint_ss[count] == 'H' ) { probabilities.push_back( helix ); }
		else if ( blueprint_ss[count] == 'E' ) { probabilities.push_back( sheet ); }
		else if ( blueprint_ss[count] == 'L' ) { probabilities.push_back( loop ); }
		else {
			utility_exit_with_message( "Unknown blueprint type " + std::string(1, blueprint_ss[count]) + " in blueprint string " + blueprint_ss  );
		}
		if ( ss == 'C' ) {
			pred_ss += 'L';
		} else {
			pred_ss += ss;
		}
		++count;
	} // while file is good
	TR << "Predicted SS = " << pred_ss << std::endl << "Wanted SS    = " << blueprint_ss << std::endl;
	runtime_assert( blueprint_ss.size() == probabilities.size() );

	PsiPredResult result;
	result.psipred_prob = probabilities;
	std::pair< std::string, utility::vector1< core::Size > > const & horiz_output( parse_psipred_horiz_output( psipred_horiz_filename ) );
	result.pred_ss = horiz_output.first;
	result.psipred2_confidence = horiz_output.second;
	for ( core::Size i=0; i<blueprint_ss.size(); ++i ) {
		if ( horiz_output.first[i] == blueprint_ss[i] ) {
			++(result.nres);
		}
	}

	return result;
}

/// @brief create fasta file, run psipred, returns psipred output if successful
PsiPredResult
PsiPredInterface::run_psipred( core::pose::Pose const & pose, std::string const & blueprint_ss ){
	/*// first, try to lookup data from the sequence of pose in the cache
	PsiPredResultMap::const_iterator it( psipred_output_.find( pose.sequence() ) );
	// if the lookup succeeded and we found something, return it
	if ( it != psipred_output_.end() ) {
		TR << "found output for " << pose.sequence() << " in the cache; returning it; " << std::endl;
		return it;
		}*/

	utility::vector1< core::pose::PoseOP > poses( pose.split_by_chain() );
	PsiPredResult result;
	core::Size start_residue = 1;

	for (core::Size i=1; i<=poses.size(); ++i) {
		std::string const fasta_filename( create_fasta_file( *(poses[i]) ) );
		core::Size replace_pos( fasta_filename.find( ".fasta" ) );
		runtime_assert( replace_pos != std::string::npos ); //sanity check
		std::string psipred_filename( fasta_filename );
		std::string psipred_horiz_filename( fasta_filename );
		psipred_filename.replace( replace_pos, 6, ".ss2" );
		psipred_horiz_filename.replace( replace_pos, 6, ".horiz" );

		//TR << "Psipred filename: " << psipred_filename << std::endl;
		// ensure we can get a shell
		runtime_assert( cmd_ != "" );
		// Call psipred on the fasta file
		std::string command = cmd_ + " " + fasta_filename;
		#ifndef WIN32
			command += " > /dev/null";
		#endif
			TR.Debug << "Trying to run " << command << std::endl;

		#ifdef __native_client__
	  		core::Size retval = 1;
		#else
	  		core::Size retval = system( command.c_str() );
		#endif
	  	if ( retval != 0 ){
			utility_exit_with_message( "Failed to run the psipred command, which was \"" + command + "\". Something went wrong. Make sure you specified the full path to the psipred command in your XML file. Return code=" + boost::lexical_cast<std::string>( retval ) );
		}

		// open and read psipred output into a stringstream buffer
		utility::io::izstream t( psipred_filename );
		std::stringstream buffer;
		buffer << t.rdbuf();
		t.close();

		std::string chain_blueprint(blueprint_ss.substr(start_residue-1, poses[i]->total_residue()));

		PsiPredResult chain_result( parse_psipred_output( *(poses[i]), buffer.str(), chain_blueprint, psipred_horiz_filename ) );

		result.psipred_prob.insert(result.psipred_prob.end(), chain_result.psipred_prob.begin(), chain_result.psipred_prob.end());
		result.psipred2_confidence.insert(result.psipred2_confidence.end(), chain_result.psipred2_confidence.begin(), chain_result.psipred2_confidence.end());
		result.pred_ss += chain_result.pred_ss;
		result.nres += chain_result.nres;

		// clean up psipred files
		cleanup_after_psipred( psipred_filename );

		start_residue += poses[i]->total_residue();

	// parse and save result in the cache
	}

	return result;
}

void
PsiPredInterface::cleanup_after_psipred( std::string const & psipred_filename ) const {
	// clean up the files created by psipred
	// psipred also creates .ss and .horiz files
	std::string ss_filename( psipred_filename.substr( 0, psipred_filename.find( ".ss2" ) ) + ".ss" );
	std::string horiz_filename( psipred_filename.substr( 0, psipred_filename.find( ".ss2" ) ) + ".horiz" );
	std::string fasta_filename( psipred_filename.substr( 0, psipred_filename.find( ".ss2" ) ) + ".fasta" );
	utility::vector1< std::string > files_to_remove;
	files_to_remove.push_back( fasta_filename );
	files_to_remove.push_back( psipred_filename );
	files_to_remove.push_back( ss_filename );
	files_to_remove.push_back( horiz_filename );
	//TR << "Before removing files! len=" << files_to_remove.size() << std::endl;
	for ( core::Size i = 1; i <= files_to_remove.size(); ++i ) {
		if ( remove( files_to_remove[i].c_str() ) ) {
			TR.Debug << "Removing " << files_to_remove[i] << std::endl;
			TR.Warning << "failed to remove " << files_to_remove[i] << std::endl;
		}
	}
}

/// @brief generates scores for each residue based on psipred confidence and the desired secondary structure
utility::vector1< core::Real >
generate_prob( PsiPredResult const & psipred_result, std::string desired_ss )
{
	utility::vector1< core::Real > retval;
	for ( core::Size i=1; i<=psipred_result.psipred2_confidence.size(); ++i ) {
		// these weightings are copied from the NNMAKE source code
		core::Real pct_h( 0.3 - 0.0333*psipred_result.psipred2_confidence[i] );
		core::Real pct_e( 0.3 - 0.0333*psipred_result.psipred2_confidence[i] );
		core::Real pct_l( 0.3 - 0.0333*psipred_result.psipred2_confidence[i] );
		if ( psipred_result.pred_ss[i-1] == 'H' ) {
			pct_h = 0.3 + (0.0667*psipred_result.psipred2_confidence[i]);
		} else if ( psipred_result.pred_ss[i-1] == 'E' ) {
			pct_e = 0.3 + (0.0667*psipred_result.psipred2_confidence[i]);
		} else if ( psipred_result.pred_ss[i-1] == 'L' ) {
			pct_l = 0.3 + (0.0667*psipred_result.psipred2_confidence[i]);
		} else {
			utility_exit_with_message( std::string( "Invalid secondary structure character in psipred result: " ) + psipred_result.pred_ss[i-1] );
		}
		core::Real resi_value;
		if ( desired_ss[i-1] == 'H' ) {
			resi_value = pct_h;
		} else if ( desired_ss[i-1] == 'E' ) {
			resi_value = pct_e;
		} else if ( desired_ss[i-1] == 'L' ) {
			resi_value = pct_l;
		} else {
			utility_exit_with_message( std::string( "Invalid secondary structure character in blueprint file: " ) + desired_ss[i-1] );
		}
		retval.push_back( resi_value );
	}
	return retval;
}

// returns the residue numbers of the two SS strings that don't match up
utility::vector1< core::Size >
nonmatching_residues( std::string const & blueprint_ss, std::string const & pred_ss ) {
	utility::vector1< core::Size > residues;
	runtime_assert( blueprint_ss.size() == pred_ss.size() );
	for ( core::Size i=1; i<=blueprint_ss.size(); i++ ) {
			if ( blueprint_ss[i-1] != pred_ss[i-1] ) {
				residues.push_back( i );
			}
	}
	return residues;
}

//namespaces
}
}
}
