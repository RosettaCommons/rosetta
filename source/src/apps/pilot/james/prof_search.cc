// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file prof_align.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

#include <basic/options/option.hh>

// AUTO-REMOVED #include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
// AUTO-REMOVED #include <core/sequence/L1ScoringScheme.hh>
// AUTO-REMOVED #include <core/sequence/MatrixScoringScheme.hh>
// AUTO-REMOVED #include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

// AUTO-REMOVED #include <core/scoring/rms_util.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructFactory.hh>

// AUTO-REMOVED #include <protocols/comparative_modeling/ThreadingMover.hh>
// AUTO-REMOVED #include <protocols/jobdist/standard_mains.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

// C++ headers
#include <map>
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/sequence/ScoringScheme.hh>
#include <protocols/jobdist/Jobs.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end



///////////////////////////////////////////////////////////////////////////////

using std::map;
using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using namespace basic;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;

class SequenceProfileDB : public utility::pointer::ReferenceCount {
public:
	SequenceProfileDB( utility::file::FileName const & fn ) {
		init_private_data_();
		open_file( fn );
	}

	~SequenceProfileDB() {
		input_.close();
	}


	SequenceProfileOP get_next_profile() {
		if ( !input_ ) {
			std::string const msg(
				"ERROR: Trying to create a SequenceProfile, but I don't have an open istream!"
			);
			utility_exit_with_message(msg);
		}

		bool finished( false );

		utility::vector1< utility::vector1< core::Real > > prof_rows;
		std::string sequence, id;
		while ( !finished ) {
			std::string line, tag;
			getline( input_, line );
			if ( input_.eof() ) {
				break; // end of file!
			}

			std::istringstream line_stream( line );
			line_stream >> tag;
			if ( tag == "ENTRY:" ) {
				//std::cout << "read id from line " << line << std::endl;
				line_stream >> id;
			} else if ( tag == "SEQUENCE:" ) {
				//std::cout << "read sequence from line " << line << std::endl;
				line_stream >> sequence;
			} else if ( tag == "--" ) {
				//std::cout << "finished with line " << line << std::endl;
				finished = true;
			} else { // we should be on a profile line
				//std::cout << "making profile from line " << line << std::endl;
				//std::cout << "read tag " << tag << std::endl;
				utility::vector1< core::Real > profile_row( 20, 0.0 );
				//std::cout << "read tag " << tag << std::endl;
				profile_row[ order_[1] ] = float_of(tag);

				for ( Size idx = 2; idx <= 20; ++idx ) {
					line_stream >> tag;
					//std::cout << "read tag " << tag << std::endl;
					profile_row[ order_[idx] ] = float_of(tag);
				}

				prof_rows.push_back( profile_row );
			}
		} // while (!finished)


		SequenceProfileOP next_prof = new SequenceProfile( prof_rows, sequence, id );
		next_prof->id( id );
		next_prof->sequence( sequence );
		return next_prof;
	} // get_next_profile

	bool has_another_profile() {
		input_.peek();
		if ( !input_.eof() ) {
			return true;
		} else {
			return false;
		}
	}

private:
	void open_file( FileName const & fn ) {
		input_.open( fn );
		if ( !input_ ) {
			utility_exit_with_message(
				"ERROR: Unable to open file (" + fn.name() + ")"
			);
		}
	}


	void init_private_data_() {
		order_.resize( 20 );

		utility::vector1< char > aa_names;
		aa_names.resize( 20 );
		aa_names[ 1] = 'A';
		aa_names[ 2] = 'R';
		aa_names[ 3] = 'N';
		aa_names[ 4] = 'D';
		aa_names[ 5] = 'C';
		aa_names[ 6] = 'Q';
		aa_names[ 7] = 'E';
		aa_names[ 8] = 'G';
		aa_names[ 9] = 'H';
		aa_names[10] = 'I';
		aa_names[11] = 'L';
		aa_names[12] = 'K';
		aa_names[13] = 'M';
		aa_names[14] = 'F';
		aa_names[15] = 'P';
		aa_names[16] = 'S';
		aa_names[17] = 'T';
		aa_names[18] = 'W';
		aa_names[19] = 'Y';
		aa_names[20] = 'V';

		for ( Size i = 1; i <= aa_names.size(); ++i ) {
			order_[i] = core::chemical::aa_from_oneletter_code( aa_names[i] );
		}

		has_another_prof_ = true;
	} // init_private_data_

	bool has_another_prof_;
	utility::vector1< core::chemical::AA > order_;
	utility::io::izstream input_;
}; // SequenceProfileDB

int
main( int argc, char* argv [] ) {
	try {

	using core::Real;
	using core::Size;

	// options, random initialization
	devel::init( argc, argv );

	// query and template profiles
	FileName pssm_file( option[ in::file::pssm ]()[1] );
	SequenceProfileOP query_prof( new SequenceProfile );
	query_prof->read_from_file( pssm_file );
	query_prof->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()

	FileName pssm_db  ( option[ in::file::pssm ]()[2] );
	SequenceProfileDB prof_db( pssm_db );

	// scoring scheme for aligning profiles
	std::string const scoring_scheme_type( option[ frags::scoring::profile_score ]() );
	ScoringSchemeFactory ssf;
	ScoringSchemeOP ss( ssf.get_scoring_scheme( scoring_scheme_type ) );
	//ScoringSchemeOP ss( new MatrixScoringScheme( -4, -1, FileName("BLOSUM62") );

	// for the ProfSim scoring scheme, the optimal opening and extension
	// penalties were 2 and 0.2, with a scoring shift of -0.45 applied to
	// all ungapped and aligned pairs.
	Real const gap_open    (  -2   );
	Real const gap_extend  (  -0.5 );
	ss->gap_open  ( gap_open   );
	ss->gap_extend( gap_extend );

	// construct alignments
	SWAligner sw_aligner;
	NWAligner nw_aligner;

	std::string output_fn = option[ out::file::alignment ]();
	utility::io::ozstream output( output_fn );

	query_prof->convert_profile_to_probs();
	while ( prof_db.has_another_profile() ) {
		SequenceProfileOP db_prof = prof_db.get_next_profile();
		db_prof->convert_profile_to_probs();
		std::cout << "generating alignment of " << query_prof->id() << " with "
			<< db_prof->id() << std::endl;
		SequenceAlignment local_align = sw_aligner.align( query_prof, db_prof, ss );
		output << local_align << "--" << std::endl;
		//SequenceAlignment global_align = nw_aligner.align( query_prof, db_prof, ss );
		//output << global_align << "--" << std::endl;
		//output << *db_prof << std::endl;
	}

	output.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
