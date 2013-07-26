// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/BinaryRNASilentStruct.cc
///
/// @brief
/// @author Rhiju Das

// C++ Headers
#include <cmath>
#include <cstdlib>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <utility>
#include <vector>
// AUTO-REMOVED #include <list>
#include <string>
#include <map>
#include <sstream>

// mini headers
#include <ObjexxFCL/string.functions.hh>

#include <utility/Binary_Util.hh>

#include <basic/Tracer.hh>

// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>

// AUTO-REMOVED #include <basic/options/option.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/NamedStubID.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

#include <numeric/model_quality/rms.hh>


// option key includes
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <basic/options/keys/OptionKeys.hh>



namespace core {
namespace io {
namespace silent {

static basic::Tracer tr("core.io.silent");


/// @brief Constructors.
BinaryRNASilentStruct::BinaryRNASilentStruct( Size const nres_in )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	fullatom_    = true; //option[ in::file::fullatom ]();
	bJumps_use_IntraResStub_ = false;
	nres  ( nres_in );
	resize( nres_in );
}

BinaryRNASilentStruct::BinaryRNASilentStruct()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	fullatom_    = true; //ption[ in::file::fullatom ]();
	bJumps_use_IntraResStub_ = false;
	nres( 0 );
	decoy_tag( "empty" );
}


BinaryRNASilentStruct::BinaryRNASilentStruct(
	core::pose::Pose const & pose,
	std::string tag
) {
	fullatom_ = true;
	bJumps_use_IntraResStub_ = false;
	fill_struct( pose, tag );
} // BinaryRNASilentStruct

void
BinaryRNASilentStruct::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {
	tr.Trace << "binary:fill_struct... " << std::endl;
	decoy_tag( tag );

	if ( tag == "empty_tag" ) set_tag_from_pose( pose );

	fullatom_ = !pose.residue(1).is_coarse();
	//	using namespace core::chemical;
	//	if ( pose.residue(1).residue_type_set().name() == ChemicalManager::get_instance()->residue_type_set(FA_STANDARD)->name() ) {
	//		fullatom_ = true;
	//	} else {
	//		fullatom_ = false;
	//	}
	tr.Trace << "get energies from pose..." << std::endl;
	energies_from_pose( pose );

	// conformation information
	//sequence_ = pose.annotated_sequence();
	sequence( pose.annotated_sequence() );
	resize( pose.total_residue() );

	tr.Trace << "read coords..." << std::endl;
	for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
		core::conformation::Residue const& resi = pose.residue(i);

		atm_coords_[i].resize( resi.natoms() );
		for (unsigned int j = 1; j <= resi.natoms(); ++j) {
			atm_coords_[i][j] = resi.atom(j).xyz();
		}
		secstruct_[i] = pose.secstruct(i);
	} // for ( unsigned int i = 1; i <= pose.total_residue(); ++i )

	fold_tree_ = pose.fold_tree();
	jumps_.clear();
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
		add_jump( pose.jump(nr) );
	}

	fill_struct_with_residue_numbers( pose ); // grabs residue numbers from pose PDBInfo object.

} // BinaryRNASilentStruct

bool BinaryRNASilentStruct::init_from_lines(
	utility::vector1< std::string > const & lines,
	SilentFileData & container
) {
	utility::vector1< std::string > energy_names_;
	utility::vector1< std::string >::const_iterator iter = lines.begin();
	if ( iter->substr(0,9) != "SEQUENCE:" ) {
		// get sequence and scorename data from the silent-file data object, because I don't have it!
		EnergyNamesOP enames = EnergyNamesOP(
			static_cast< EnergyNames * > ( container.get_shared_silent_data( energynames )() )
		);

		SimpleSequenceDataOP seqdata = SimpleSequenceDataOP(
			static_cast< SimpleSequenceData * > ( container.get_shared_silent_data( simplesequencedata )() )
		);

		sequence      ( seqdata->sequence()   );
		energy_names_ = enames ->energy_names();
	} else {
		// get sequence and scorename data from the first two lines provided,
		// put them into container for further use by other SilentStruct
		// objects.

		// first line is SEQUENCE:
		std::istringstream line_stream( *iter );
		std::string tag;
		tr.Debug << "reading sequence from " << *iter << std::endl;
		++iter;

		std::string temp_seq;
		line_stream >> tag >> temp_seq;
		if ( line_stream.fail() || tag != "SEQUENCE:" ) {
			tr.Error << "bad format in sequence line of silent file" << std::endl;
			tr.Error << "line = " << *iter << std::endl;
			tr.Error << "tag = " << tag << std::endl;
			return false;
		}
		sequence( temp_seq );

		// second line is a list of score names
		std::istringstream score_line_stream( *iter );
		tr.Debug << "reading score names from " << *iter << std::endl;
		++iter;

		score_line_stream >> tag; // SCORE:
		if ( score_line_stream.fail() || tag != "SCORE:" ) {
			tr.Error << "bad format in second line of silent file" << std::endl;
			tr.Error << "tag = "  << tag << std::endl;
			tr.Error << "line = " << *iter << std::endl;
		}

		score_line_stream >> tag; // first score name
		while ( ! score_line_stream.fail() ) {
			energy_names_.push_back( tag );
			score_line_stream >> tag; // try to get next score name
		}

		EnergyNamesOP enames( new EnergyNames() );
		SimpleSequenceDataOP seqdata( new SimpleSequenceData() );

		enames ->energy_names( energy_names_ );
		seqdata->set_sequence( sequence()    );

		container.set_shared_silent_data( energynames       , enames  );
		container.set_shared_silent_data( simplesequencedata, seqdata );
	} // get header information

	int currpos = 1;
	bool bitflip = false;
	for ( utility::vector1< std::string >::const_iterator end = lines.end();
				iter != end; ++iter
	) {
		std::string tag;
		std::istringstream line_stream( *iter );

		//		std::cout << (*iter) << std::endl;

		if ( iter->substr(0,6) == "REMARK" ){
			get_parent_remark_from_line( line_stream.str() );
			continue;  // skip comments if record_old_remarks==false
		}

		if ( iter->substr(0,7) == "SCORE: " ) {
			// SCORE: line with values from this structure.
			Size nres = one_letter_sequence().length();
			resize( nres );

			std::string tag;
			line_stream >> tag;
			if ( line_stream.fail() || tag != "SCORE:" ) {
				tr.Error << "bad format in first score line of silent file" << std::endl;
				tr.Error << "line = " << *iter << std::endl;
				tr.Error << "tag = " << tag << std::endl;
			}

			parse_energies( line_stream, energy_names_ );

		} else { // conformation lines

			if ( Size(currpos) > nres() ) continue;

			// parse fold_tree and jump lines
			if ( iter->substr(0,10) == "FOLD_TREE " ) {
				kinematics::FoldTree f;
				line_stream >> f;
				//NOTE!!!!!!!!!
				// In BinaryProteinSilentStruct, used the function fold_tree( f ),
				// which does something weird (e.g., doesn't use a const fold_tree as input).
				set_fold_tree( f ); // add fold-tree to this SilentStruct
				tr.Debug << "read fold-tree " << f; //"\n" is in fold-tree output
				tr.Debug << "reading " << f.num_jump() << " jumps " << std::endl;
				continue;
			} else if ( iter->substr(0,2) == "RT" ) {
				kinematics::Jump jump;
				line_stream >> jump;
				tr.Debug << "read jump " << jump << std::endl;
				add_jump( jump );
				// modern style jumps, defined completely with the FoldTree
				bJumps_use_IntraResStub_ = false;
				continue;
			} else if ( iter->substr(0,9) == "SEQUENCE:" ) {
				tr.Debug << "Skipping duplicate sequence declaration " << std::endl;
				continue;
			} else if ( iter->substr(0,19) == "ANNOTATED_SEQUENCE:" ) {
				std::string annotated_seq;
				line_stream >> tag; //ANNOTATED_SEQUENCE
				line_stream >> annotated_seq;
				sequence( annotated_seq );
				tr.Debug << "read annotated sequence as: " << sequence() << std::endl;
				// resize pose according to number of resiudes in annotated sequence
				resize( one_letter_sequence().length() );
				continue;
			} else if ( iter->substr(0,4) == "JUMP" ) {
				// support for rosetta++ silent files
				std::string tag;
				Size nr;
				line_stream >> tag; //JUMP
				line_stream >> nr;
				if ( nr != fold_tree().num_jump() ) {
					tr.Warning
						<< "WARNING: corrupted silent file read line JUMP X -- X should match number of jumps in FOLD_TREE " << std::endl;
				}
				for ( Size i = 1; i<= nr; i++ ) {
					kinematics::Jump jump;
					line_stream >> jump;
					add_jump( jump );
				}
				bJumps_use_IntraResStub_ = true;// jump is defined via N-C-CA rosetta++ style
				continue;
			} else if ( iter->substr(0,7) == "RES_NUM" ) {
				figure_out_residue_numbers_from_line( line_stream );
				continue;
			}

			// parse coords
			line_stream >> tag;

			if (tag.length() < 1) {
				tr.Warning << "WARNING:  read blank line in decoy tag " << decoy_tag() << std::endl;
				continue;
			}
			secstruct_[currpos] = tag[0];  // first char is sec struct

			int natoms = (tag.length()-1) / 16;
			utility::vector1< numeric::xyzVector <float> > atm_buff( natoms+1 );
			utility::decode6bit( (unsigned char*)&(atm_buff[1]) , tag.substr(1) );

			// endianness check ...
			//   check the dist between atoms 1 and 2 is unreasonable .. and flipping fixes then turn bitflip
			//   on

			if (currpos == 1) {
				core::Real len_check12 = (atm_buff[1]-atm_buff[2]).length();
				if ( len_check12 < 0.5 || len_check12 > 2.0 ) {
					utility::swap4_aligned ( (void*) &(atm_buff[1][0]) , 3*natoms );
						// recheck; if not better flip back
					len_check12 = (atm_buff[1]-atm_buff[2]).length();
					if ( len_check12 < 0.5 || len_check12 > 2.0 ) {
						utility::swap4_aligned ( (void*) &(atm_buff[1][0]) , 3*natoms );
					} else {
						tr.Warning << "reading big-endian binary silent file! " << decoy_tag() << std::endl;
						bitflip = true;
					}
				}
			} else {
				if (bitflip ) {
					utility::swap4_aligned ( (void*) &(atm_buff[1][0]) , 3*natoms );
				}
			}

			atm_coords_[currpos].resize( natoms ); // allocate space for coords
			for (int j=1; j<=natoms; ++j) {
				atm_coords_[currpos][j] = atm_buff[j];
			}
			currpos++;
			//tr.Debug << "processing line " << *iter << std::endl;
		} // conformation lines
	} // for ( iter ... )

	if ( fold_tree().num_jump() != jumps_.size() ) {
		tr.Warning	<< "parse error:  found " << jumps_.size()
								<< " RT lines for a fold-tree with " << fold_tree().num_jump()
								<< " for decoy tag " << decoy_tag() << std::endl;
		return false;
	}

	//if ( (unsigned int) currpos != total_residue + 1 ) {
	if ( atm_coords_.size() != nres() ) {
		tr.Error << "ERROR: didn't find coordinates for all sequence positions of "
							<< decoy_tag() << std::endl;
		tr.Error << "       expected " << nres()
							<< ", found " << currpos-1 << std::endl;
		return false; //no success
	}

	if ( fold_tree().size() < 1 ) {
		fold_tree_.simple_tree( nres() );
		tr.Debug << " generating simple fold-tree " << fold_tree();
	}

	if ( bJumps_use_IntraResStub_ ) { //for rosetta++ file-format
		//prepares of setting RT via N, CA, C
		fold_tree_.put_jump_stubs_intra_residue();
		//on could also think of making this a temporary change after read is
		//finished return to a standard fold_tree...
	}

	tr.Debug << "(TEX) FOLD TREE: " << fold_tree();

	return true;
} // init_from_lines

/// @brief Resize this silent-struct to the appropriate number of residues.
void
BinaryRNASilentStruct::resize(
	Size const nres_in
) {
	//nres_ = nres_in;
	nres( nres_in );
	secstruct_.resize( nres() );
	atm_coords_.resize( nres() );

	// make a new FoldTree if we're just replacing a simple FoldTree, otherwise
	// trust the user to have provided us a reasonable tree in a FOLD_TREE line.
	if ( fold_tree_.is_simple_tree() ) {
		fold_tree_.simple_tree( nres() );
	}
}

void BinaryRNASilentStruct::fill_pose (
	core::pose::Pose & pose
) const {
	using namespace core::chemical;
	ResidueTypeSetCAP residue_set;
	//	std::cout << "RESIDUE TYPE SET RNA " << std::endl;
	if ( one_letter_sequence()[0] != 'Z' /* Mg(2+) */ && atm_coords_[1].size() < 8 ) { //hmm, may be dangerous.
		residue_set = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );
	} else {
		residue_set = ChemicalManager::get_instance()->residue_type_set( RNA );
	}
	fill_pose( pose, *residue_set );
} // fill_pose

void BinaryRNASilentStruct::fill_pose (
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set
) const {
	core::pose::make_pose_from_sequence( pose, sequence(), residue_set );
	tr.Debug << "SEQUENCE: " << sequence();
	tr.Debug << "FOLD TREE: " << fold_tree();

	// set fold_tree
	pose.fold_tree( fold_tree() );

	// WE DON'T NEED THIS, SINCE XYZ OF ALL ATOMS IS DEFINED... RIGHT?
	// set jumps
	//	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
	//		if ( !bJumps_use_IntraResStub_ ) { //default modern file-format
	//			pose.set_jump( nr, jump( nr ) );
	//		}
	//	}

	tr.Debug << "nres = " << nres() << std::endl;
	tr.Debug << "one_letter_sequence() = " << one_letter_sequence().length() << std::endl;

	if( nres() != one_letter_sequence().length() ){
		utility_exit_with_message( "RuntimeAssert failed: nres() == one_letter_sequence().length()" );
	}

	// coords
	for ( Size seqpos = 1; seqpos <= nres(); ++seqpos ) {
		int natoms_pose = pose.residue_type(seqpos).natoms() ,
		    natoms_struct = atm_coords_[seqpos].size();
		int natoms_total = std::min( natoms_pose , natoms_struct );

		if ( natoms_pose != natoms_struct) {
			tr.Warning 	<< "[ WARNING ] Number of atoms in pose and silent file disagree! ";
			tr.Warning 	<< "Attempting to continue ..." << std::endl;
			tr.Warning 	<< "[ WARNING ]    (in residue "
									<< seqpos << "  natoms_pose=" << natoms_pose
									<< "  natoms_struct=" << natoms_struct << ")" << std::endl;
		}

		natoms_pose  = pose.residue_type(seqpos).natoms();
		natoms_total = std::min( natoms_pose, natoms_struct );

		for ( int j = 1; j <= natoms_total; ++j ){
			id::AtomID id( j, seqpos );
			numeric::xyzVector< core::Real> atom_i(atm_coords_[seqpos][j][0], atm_coords_[seqpos][j][1], atm_coords_[seqpos][j][2]);
			pose.set_xyz( id, atom_i );
		}
		pose.set_secstruct( seqpos, secstruct_[seqpos] );

	} // for ( seqpos )

	tr.Debug << "Hallelujah! " << pose.total_residue() << std::endl;

	finish_pose( pose );

} // fill_pose

void
BinaryRNASilentStruct::print_header( std::ostream& out ) const
{
	SilentStruct::print_header( out );
	if ( fullatom_ ) {
		out << "REMARK BINARY_SILENTFILE RNA \n";
	} else {
		out << "REMARK BINARY_SILENTFILE RNA COARSE\n";
	}
}


void BinaryRNASilentStruct::print_conformation(
	std::ostream & output
) const {
	if ( fold_tree().size() > 1 || fold_tree().num_jump() > 0 ) { //assume non-trivial fold_tree only if more than one edge, i.e., EDGE 1 <nres> -1
		output << "FOLD_TREE ";
		for ( kinematics::FoldTree::const_iterator
					it = fold_tree().begin(), it_end = fold_tree().end();
					it != it_end; ++it
		) {
			output << *it;
		}
		//		output << fold_tree(); this produces a new-line --- wrong behaviour
		//		of fold_tree but I don't want to fix 1000 u-tracer unit-tests!
		output << ' ' << decoy_tag() << "\n";
	}
	for ( Size i = 1; i <= fold_tree().num_jump(); i++ ) {
		output << jump( i ) << ' ' << decoy_tag() << "\n";
	}
	output << "ANNOTATED_SEQUENCE: " << sequence() << " " << decoy_tag() << "\n"; //chu print annotated_sequence per decoy
	//tr.Debug << "FOLD_TREE Size: " << fold_tree().size() << " " << fold_tree() << std::endl;

	// fullatom flag
	//int fullatom_flag = (fullatom_? 1 : 2);
	std::string resline;
	//encode6bit(  (unsigned char*)&fullatom_flag, 4, resline );
	//output << resline << "\n";

	for ( Size i = 1; i <= nres(); ++i ) {
		// make sure secstruct is valid
		char this_secstr = secstruct_[i];
		if (this_secstr < 'A' || this_secstr > 'Z') this_secstr = 'L';
		utility::encode6bit(  (unsigned char*)&atm_coords_[i][1][0], atm_coords_[i].size()*12, resline );  // ASSUMES FLOAT == 4 BYTES!!! (eep!)
		output << this_secstr << resline << ' ' << decoy_tag() << "\n";
	} // for ( Size i = 1; i <= nres; ++i )
} // print_conformation


Real BinaryRNASilentStruct::get_debug_rmsd() {
	pose::Pose temp_pose;
	ObjexxFCL::FArray2D< Real > rebuilt_coords ( 3, atm_coords_.size() ),
		original_coords( 3, atm_coords_.size() );

	// build temp_pose from coordinates
	fill_pose( temp_pose );

	Size const c4star_index = temp_pose.residue(1).atom_index( " C4'" );

	for ( Size i = 1; i <= temp_pose.total_residue(); ++i ) {
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			rebuilt_coords (k,i) = temp_pose.residue(i).xyz( " C4'" )[k-1];
			original_coords(k,i) = atm_coords_[i][c4star_index][k-1];
		}
	}

	Real rmsd = numeric::model_quality::rms_wrapper( temp_pose.total_residue(), rebuilt_coords, original_coords );
	return rmsd;
}

ObjexxFCL::FArray2D< Real >
BinaryRNASilentStruct::get_CA_xyz() const {
	core::Size n_residues = nres();
	ObjexxFCL::FArray2D< Real > my_coords( 3, n_residues );
	for ( Size i = 1; i <= n_residues; ++i ) { // i = n_residues
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			my_coords(k,i) = atm_coords_[i][2][k-1];
		} // k
	} // i

	return my_coords;
} // get_CA_positions

Real BinaryRNASilentStruct::CA_rmsd( RNA_SilentStruct other_pss ) {
	ObjexxFCL::FArray2D< Real > my_coords    = get_CA_xyz();
	ObjexxFCL::FArray2D< Real > other_coords = other_pss.get_CA_xyz();
	Real rmsd = numeric::model_quality::rms_wrapper( nres(), my_coords, other_coords );

	return rmsd;
} // RNA_SilentStruct::CA_rmsd


} // namespace silent
} // namespace io
} // namespace core
