// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    mp_mutate_relax.cc
/// @brief   Mutate a residue, then do range relax for a membrane protein
/// @author  JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_MPMutateRelaxMover_cc
#define INCLUDED_protocols_membrane_MPMutateRelaxMover_cc

// Unit Headers
#include <protocols/relax/membrane/MPMutateRelaxMover.hh>
#include <protocols/relax/membrane/MPMutateRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/relax/membrane/MPRangeRelaxMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/membrane/util.hh>
#include <utility/io/util.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <core/pose/util.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.MPMutateRelaxMover" );

namespace protocols {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Mutations will be read in by commandline
MPMutateRelaxMover::MPMutateRelaxMover() : protocols::moves::Mover()
{
	register_options();
	init_from_cmd();
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
MPMutateRelaxMover::MPMutateRelaxMover( MPMutateRelaxMover const & ) = default;

/// @brief Assignment Operator
MPMutateRelaxMover & MPMutateRelaxMover::operator = ( MPMutateRelaxMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new MPMutateRelaxMover( *this ) );
}

/// @brief Destructor
MPMutateRelaxMover::~MPMutateRelaxMover() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPMutateRelaxMover::clone() const {
	return ( protocols::moves::MoverOP( new MPMutateRelaxMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPMutateRelaxMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MPMutateRelaxMover() );
}

/// @brief Parse Rosetta Scripts Options for this Mover
void
MPMutateRelaxMover::parse_my_tag(
	utility::tag::TagCOP /*tag*/,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// TODO: implement this

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
MPMutateRelaxMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MPMutateRelaxMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MPMutateRelaxMoverCreator::keyname() const {
	return MPMutateRelaxMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MPMutateRelaxMoverCreator::mover_name() {
	return "MPMutateRelaxMover";
}

/// @brief Get the name of this Mover (MPMutateRelaxMover)
std::string
MPMutateRelaxMover::get_name() const {
	return "MPMutateRelaxMover";
}

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Mutate residue and then range relax the membrane protein
void MPMutateRelaxMover::apply( core::pose::Pose & pose ) {

	using namespace utility;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::membrane;
	using namespace protocols::relax::membrane;
	using namespace protocols::simple_moves;
	using namespace core::scoring;
	using namespace core::pack::task;

	// finalize setup
	finalize_setup( pose );

	TR << "Running MPMutateRelax protocol with SCOREFUNCTION " << utility::CSI_Green() << sfxn_->get_name() << utility::CSI_Reset() << std::endl;

	// final foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// WHAT IS HAPPENING HERE:
	// I am dumping the poses outside of JD2
	// Why would I do that?
	// because otherwise I would need to put the read_mutant_file function into
	//   the pilot app (maybe people want to use it?) and I want the mutant to
	//   be part of the output filename

	// make a deepcopy of the pose to start from
	Pose original_pose = Pose( pose );
	Pose working_pose;

	// construct counter
	Size counter(0);

	// go through each construct (i.e. outer vector)
	for ( core::Size c = 1; c <= wt_res_.size(); ++c ) {

		TR << "going through construct " << c << std::endl;

		// counter
		counter = 1;
		std::string mutations;

		// iterate over nstruct
		while ( counter <= nstruct_ ) {

			TR << "working on nstruct " << counter << std::endl;

			// get original starting pose
			working_pose = Pose( original_pose );

			// make mutations
			std::string mutations = make_mutations( working_pose, c );

			// if repacking
			if ( repack_mutation_only_ == true || repack_radius_ > 0 ) {

				// create task factory and repack
				TR << "Repacking only..." << std::endl;
				PackerTaskOP repack = TaskFactory::create_packer_task( working_pose );
				TR.Debug << "repacking vector size: " << repack_residues_[c].size() << std::endl;
				TR.Debug << "pose length: " << working_pose.size() << std::endl;
				repack->restrict_to_residues( repack_residues_[c] );
				repack->restrict_to_repacking();
				core::pack::pack_rotamers( working_pose, *sfxn_, repack );

				// if relaxing
			} else if ( relax_ == true ) {

				// do range relax
				TR << "Running MPRangeRelax..." << std::endl;
				MPRangeRelaxMoverOP relax( new MPRangeRelaxMover() );
				//    relax->optimize_membrane( false );
				relax->apply( working_pose );
			}

			// create output filename
			std::string output;

			// if model exists, increment counter
			// this means that the app should be able to start from an existing file number
			Size a = counter;
			while ( a <= nstruct_ ) {
				output = output_filename( mutations, a );
				if ( utility::file::file_exists( output ) ) {
					continue;
				} else {
					break;
				}
				++a;
			}

			// dump pose
			working_pose.dump_scored_pdb( output, *sfxn_ );

			// increment counter
			++counter;

		}// nstruct
	}// iterate over construct

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void MPMutateRelaxMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::in::file::s );
	option.add_relevant( OptionKeys::mp::mutate_relax::mutation );
	option.add_relevant( OptionKeys::mp::mutate_relax::mutant_file );
	option.add_relevant( OptionKeys::mp::mutate_relax::nmodels );
	option.add_relevant( OptionKeys::mp::mutate_relax::repack_mutation_only );
	option.add_relevant( OptionKeys::mp::mutate_relax::repack_radius );
	option.add_relevant( OptionKeys::mp::mutate_relax::relax );

}

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize from commandline
void MPMutateRelaxMover::init_from_cmd() {

	using namespace basic::options;

	// making sure models won't be overwritten
	if ( option[ OptionKeys::out::nstruct ]() > 1 ) {
		utility_exit_with_message( "Set '-out:nstruct 1'! Otherwise models will be overwritten!" );
	}

	// input checking
	if ( ! option[ OptionKeys::mp::mutate_relax::mutant_file ].user() &&
			! option[ OptionKeys::mp::mutate_relax::mutation ].user() ) {

		utility_exit_with_message( "Too few inputs: You must EITHER specify the mutant file with -mp:mutate_relax:mutant_file OR specify a single mutation with -mp:mutate_relax:mutation! Quitting..." );
	}

	// more input checking
	if ( option[ OptionKeys::mp::mutate_relax::mutant_file ].user() &&
			option[ OptionKeys::mp::mutate_relax::mutation ].user() ) {

		utility_exit_with_message( "Too many inputs: You must EITHER specify the mutant file with -mp:mutate_relax:mutant_file OR specify a single mutation with -mp:mutate_relax:mutation! Quitting..." );
	}

	// get mutant file
	if ( option[ OptionKeys::mp::mutate_relax::mutant_file ].user() ) {
		mutant_file_ = option[ OptionKeys::mp::mutate_relax::mutant_file ]();
	}

	// Number of iterations to run
	nstruct_ = 1;
	if ( option[ OptionKeys::mp::mutate_relax::nmodels ].user() ) {
		nstruct_ = option[ OptionKeys::mp::mutate_relax::nmodels ]();
	}

	// get protein name for dumping PDBs
	if ( option[ OptionKeys::in::file::s ].user() ) {
		protein_ = option[ OptionKeys::in::file::s ](1);
	} else {
		utility_exit_with_message("No PDB given, please use -in:file:s to provide PDB.");
	}

	// repack options, overwrite iterations to 1
	repack_mutation_only_ = false;
	if ( option[ OptionKeys::mp::mutate_relax::repack_mutation_only ].user() ) {
		repack_mutation_only_ = option[ OptionKeys::mp::mutate_relax::repack_mutation_only ]();
		nstruct_ = 1;
		TR << "Repacking only: setting number of iterations to 1." << std::endl;
	}

	// repacking within radius, overwrite iterations to 1
	repack_radius_ = 0;
	if ( option[ OptionKeys::mp::mutate_relax::repack_radius ].user() ) {
		repack_radius_ = option[ OptionKeys::mp::mutate_relax::repack_radius ]();
		nstruct_ = 1;
		TR << "Repacking only: setting number of iterations to 1." << std::endl;
	}

	// only running relax
	relax_ = false;
	if ( option[ OptionKeys::mp::mutate_relax::relax ].user() ) {
		relax_ = option[ OptionKeys::mp::mutate_relax::relax ]();
		TR << "Setting relax option to " << relax_ << std::endl;
	}

	// check number of iterations
	if ( relax_ == true && nstruct_ == 0 ) {
		nstruct_ = 100;
		TR << "Relaxing structures: setting number of iterations to 100." << std::endl;
	} else if ( relax_ == true ) {
		TR << "Relaxing structures: number of iterations is " << nstruct_ << " as per user-input." << std::endl;
		TR << "I assume you know what you are doing, a good number of iterations is 100." << std::endl;
	} else if ( relax_ == false && repack_mutation_only_ == false && repack_radius_ == 0 ) {
		nstruct_ = 1;
		TR << "Neither repacking nor relaxing structures: setting number of iterations to 1." << std::endl;
	}

	// checking inputs
	if ( repack_mutation_only_ == true && repack_radius_ > 0 ) {
		utility_exit_with_message("Set EITHER repack_mutation_only OR repack_radius, not both!");
	}
	if ( ( repack_mutation_only_ == true || repack_radius_ > 0 ) && relax_ == true ) {
		utility_exit_with_message("Set EITHER repacking option OR relax option, not both!");
	}

}// init from cmdline

////////////////////////////////////////////////////////////////////////////////

/// @brief Get repack residues
void MPMutateRelaxMover::get_repack_residues( Pose & pose ) {

	using namespace core::pose;

	// go through number of constructs
	for ( Size i = 1; i <= resn_.size(); ++i ) {

		// initialize boolean vector with false
		utility::vector1< bool > repack_res( pose.size(), false );

		// if repacking
		if ( repack_mutation_only_ == true ) {

			// go through number of mutations within construct
			for ( Size j = 1; j <= resn_[i].size(); ++j ) {

				// set residue number of the mutation to true
				core::Size resn = resn_[i][j];
				repack_res[ resn ] = true;
			}
		} else if ( repack_radius_ > 0 ) {

			// go through number of mutations within construct
			for ( Size j = 1; j <= resn_[i].size(); ++j ) {

				core::Size resn = resn_[i][j];
				core::Vector xyz_mut = pose.residue( resn ).xyz( "CA" );

				// go through residues in pose
				for ( Size k = 1; k <= nres_protein( pose ); ++k ) {

					// check whether the residue is within repack radius
					core::Vector xyz_k = pose.residue( k ).xyz( "CA" );
					core::Real dist = ( xyz_k - xyz_mut ).length();

					// if yes, set repack flag of this residue to true
					if ( dist <= repack_radius_ ) {
						repack_res[ k ] = true;
					}
				}
			}
		}

		// add vector describing construct to total constructs
		repack_residues_.push_back( repack_res );
	}

} // get repack residues

////////////////////////////////////////////////////////////////////////////////

/// @brief Finalize setup
void MPMutateRelaxMover::finalize_setup( Pose & pose ){

	using namespace basic::options;

	// get mutation
	if ( option[ OptionKeys::mp::mutate_relax::mutation ].user() ) {

		// input format A163F
		std::string mutation = option[ OptionKeys::mp::mutate_relax::mutation ]();

		TR << "Looking at mutation " << mutation << std::endl;

		// add string to private data
		add_mutant_to_vectors( mutation, pose );
	}

	// read mutant file
	if ( mutant_file_.size() > 0 ) {
		read_mutant_file( pose );
	}

	// error checking
	check_mutant_file( pose );

	// call AddMembraneMover
	if ( ! pose.conformation().is_membrane() ) {
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );
	}

	// get repack residues
	get_repack_residues( pose );

	// create scorefunction
	sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );

}// finalize setup

////////////////////////////////////////////////////////////////////////////////

/// @brief Make mutations, returns string of output file
std::string MPMutateRelaxMover::make_mutations( Pose & pose, core::Size num_construct ) {

	using namespace protocols::simple_moves;
	using namespace utility;

	// mutations
	std::string mutations = "";
	core::Size c = num_construct;

	// go through each mutation for this construct
	for ( Size m = 1; m <= wt_res_[ c ].size(); ++m ) {

		// mutate residue
		TR << "Mutating residue " << wt_res_[c][m] << resn_[c][m] << " to " << new_res_[c][m] << std::endl;
		MutateResidueOP mutate( new MutateResidue( resn_[c][m], one2three( new_res_[c][m] ) ) );
		mutate->apply( pose );

		// get mutation as output tag
		mutations += "_" + wt_res_[c][m] + to_string( resn_[c][m] ) + new_res_[c][m];

	}// iterate over mutations in construct

	return mutations;
} // make mutations

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize from commandline
void MPMutateRelaxMover::read_mutant_file( core::pose::Pose & pose ) {

	using namespace utility;
	using namespace utility::io;
	TR << "reading mutant file" << std::endl;

	// get file content
	utility::vector1< std::string > lines = get_lines_from_file_data( mutant_file_ );

	// go through lines
	for ( core::Size i = 1; i <= lines.size(); ++i ) {

		// add mutants in line to private data vectors
		add_mutant_to_vectors( lines[ i ], pose );
	}

} // read mutant file

////////////////////////////////////////////////////////////////////////////////

/// @brief Add mutants to private data: list of A163F into vectors
void MPMutateRelaxMover::add_mutant_to_vectors( std::string mutations, core::pose::Pose & pose ) {

	using namespace utility;
	TR << "adding mutants to vectors" << std::endl;

	// initialize line vectors
	utility::vector1< std::string > wildtypes;
	utility::vector1< core::Size > seqids;
	utility::vector1< std::string > mutants;

	// initialize single point variables
	std::string wt, mut;
	core::Size seqid;

	// PDB numbering?
	bool pdb_numbering( false );

	// split string by whitespace, write output into vector
	utility::vector1< std::string > all_mutants = split_whitespace( mutations );

	// iterate over mutants
	for ( core::Size col = 1; col <= all_mutants.size(); ++col ) {

		// get single mutant
		std::string wt_id_mut = all_mutants[ col ];

		// PDB or pose numbering?
		// PDB numbering: A_E15L, chain, underscore, wt, resn, mutant
		// pose numbering: E15L, wt, resn, mutant
		// if string contains underscore, then PDB numbering
		if ( col == 1 ) {
			if ( wt_id_mut.find( "_" ) != std::string::npos ) {
				TR << "Mutations are in PDB numbering scheme." << std::endl;
				pdb_numbering = true;
			} else {
				TR << "Mutations are in pose numbering scheme." << std::endl;
			}
		}

		// check input file format
		if ( ! pdb_numbering && wt_id_mut.find( "_" ) != std::string::npos ) {
			utility_exit_with_message( "Mutations should either be in PDB numbering (A_E15L) where A is the chain, or in pose numbering (E15L) where the residues are renumbered without gaps, starting from 1. Quitting." );
		} else if ( pdb_numbering && wt_id_mut.find( "_" ) == std::string::npos ) {
			utility_exit_with_message( "Mutations should either be in PDB numbering (A_E15L) where A is the chain, or in pose numbering (E15L) where the residues are renumbered without gaps, starting from 1. Quitting." );
		}

		// pose numbering
		if ( pdb_numbering == false ) {

			// get wt and mutant from vector: split string by character
			wt = wt_id_mut[ 0 ];
			mut = wt_id_mut[ wt_id_mut.size()-1 ];

			// get residue number: get each digit
			utility::vector1< std::string > tmp;
			for ( core::Size i = 1; i <= wt_id_mut.size()-2; ++i ) {
				tmp.push_back( to_string( wt_id_mut[ i ] ) );
			}

			// join digits together to the residue number
			seqid = string2Size( join( tmp, "" ) );
		} else {
			// PDB numbering

			// get chain
			char chain = wt_id_mut[ 0 ];

			// remove the chain and the underscore from the string
			wt_id_mut.erase( 0, 2 );

			// get wt and mutant from vector: split string by character
			wt = wt_id_mut[ 0 ];
			mut = wt_id_mut[ wt_id_mut.size()-1 ];

			// get residue number: get each digit
			utility::vector1< std::string > tmp;
			for ( core::Size i = 1; i <= wt_id_mut.size()-2; ++i ) {
				tmp.push_back( to_string( wt_id_mut[ i ] ) );
			}

			// join digits together to the residue number
			int seqid_pdb = string2int( join( tmp, "" ) );

			// get pose resnumber from pdb numbering
			seqid = pose.pdb_info()->pdb2pose( chain, seqid_pdb );
		}

		TR << "wt " << wt << ", seqid " << seqid << ", mut " << mut << std::endl;

		// put all of the wt / seqid / mut info into the line vectors
		wildtypes.push_back( wt );
		seqids.push_back( seqid );
		mutants.push_back( mut );

	} // iterate over columns in line

	// add the line to the file (or outer vector in this case)
	wt_res_.push_back( wildtypes );
	resn_.push_back( seqids );
	new_res_.push_back( mutants );

} // add mutants to private data

////////////////////////////////////////////////////////////////////////////////

/// @brief Check mutant file for errors
/// @details If Rosetta doesn't start crying, you're good to go
void MPMutateRelaxMover::check_mutant_file( core::pose::Pose & pose ) {

	using namespace utility;
	TR << "checking mutant file" << std::endl;

	// go through outer loop
	for ( core::Size i = 1; i <= wt_res_.size(); ++i ) {

		// go through inner loop
		for ( core::Size j = 1; j <= wt_res_[ i ].size(); ++j ) {

			// check whether wt residue has the same identity in the pose
			if ( wt_res_[ i ][ j ] != to_string( pose.residue_type( resn_[ i ][ j ] ).name1() ) ) {
				TR << "i " << i << ", j " << j << std::endl;
				TR << "residue identity in PDB " << to_string( pose.residue_type( j ).name1() ) << ", " << wt_res_[ i ][ j ] << std::endl;
				utility_exit_with_message( "Residue identity in input file doesn't match the pose!" );
			}

		} // inner loop
	} // outer loop
} // check mutant file

////////////////////////////////////////////////////////////////////////////////

/// @brief Convert one to three letter code
std::string MPMutateRelaxMover::one2three( std::string one ) {

	// error checking
	if ( one == "B" || one == "J" || one == "O" || one == "U" || one == "X" || one == "Z" ) {
		utility_exit_with_message( "One letter code doesn't belong to any of the 20 natural amino acids!" );
	}

	// create the vectors
	utility::vector1< std::string > olc;
	olc.push_back( "A" );
	olc.push_back( "C" );
	olc.push_back( "D" );
	olc.push_back( "E" );
	olc.push_back( "F" );
	olc.push_back( "G" );
	olc.push_back( "H" );
	olc.push_back( "I" );
	olc.push_back( "K" );
	olc.push_back( "L" );
	olc.push_back( "M" );
	olc.push_back( "N" );
	olc.push_back( "P" );
	olc.push_back( "Q" );
	olc.push_back( "R" );
	olc.push_back( "S" );
	olc.push_back( "T" );
	olc.push_back( "V" );
	olc.push_back( "W" );
	olc.push_back( "Y" );

	utility::vector1< std::string > tlc;
	tlc.push_back( "ALA" );
	tlc.push_back( "CYS" );
	tlc.push_back( "ASP" );
	tlc.push_back( "GLU" );
	tlc.push_back( "PHE" );
	tlc.push_back( "GLY" );
	tlc.push_back( "HIS" );
	tlc.push_back( "ILE" );
	tlc.push_back( "LYS" );
	tlc.push_back( "LEU" );
	tlc.push_back( "MET" );
	tlc.push_back( "ASN" );
	tlc.push_back( "PRO" );
	tlc.push_back( "GLN" );
	tlc.push_back( "ARG" );
	tlc.push_back( "SER" );
	tlc.push_back( "THR" );
	tlc.push_back( "VAL" );
	tlc.push_back( "TRP" );
	tlc.push_back( "TYR" );

	// do the conversion
	for ( core::Size i = 1; i <= 20; ++i ) {
		if ( olc[ i ] == one ) {
			return tlc[ i ];
		}
	}

	return "WARNING: Your mutant seems to be non-canonical. Please use a different application!";

} // one to three letter code

////////////////////////////////////////////////////////////////////////////////

/// @brief Create output filename
std::string MPMutateRelaxMover::output_filename( std::string mutation_tag, core::Size counter ) {

	using namespace basic::options;
	using namespace utility;

	// create output filename to /my/path/myprotein.pdb to myprotein
	const std::string tmp( file_basename( protein_ ) );
	std::string output;

	// get pathname
	if ( option[ OptionKeys::out::path::pdb ].user() ) {
		output = option[ OptionKeys::out::path::pdb ]().path() + trim( tmp, ".pdb");
	} else {
		output = trim( tmp, ".pdb");
	}

	// add mutation to output filename
	output += mutation_tag;

	// add counter to output filename
	std::string cnt;
	if ( counter < 10 ) {
		cnt = "000" + to_string( counter );
	} else if ( counter > 9 && counter < 100 ) {
		cnt = "00" + to_string( counter );
	} else if ( counter > 99 && counter < 1000 ) {
		cnt = "0" + to_string( counter );
	} else if ( counter > 999 && counter < 10000 ) {
		cnt = to_string( counter );
	} else {
		utility_exit_with_message( "Please choose an -nstruct < 9999." );
	}
	output += "_" + cnt + ".pdb";

	return output;

}
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MPMutateRelaxMover_cc
