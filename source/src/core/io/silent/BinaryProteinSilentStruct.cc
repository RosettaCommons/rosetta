// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/BinaryProteinSilentStruct.cc
///
/// @brief
/// @author Frank DiMaio
/// @author Mike Tyka

// C++ Headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <sstream>

// mini headers
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>
// AUTO-REMOVED #include <utility/string_util.hh>

#include <core/scoring/Energies.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/io/raw_data/DisulfideFile.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Conformation.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedStubID.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

#include <numeric/model_quality/rms.hh>

#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <boost/lexical_cast.hpp>

#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <utility/Binary_Util.hh>


static basic::Tracer tr("core.io.silent");

namespace core {
namespace io {
namespace silent {



/// @brief Constructors.
BinaryProteinSilentStruct::BinaryProteinSilentStruct( Size const nres_in )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	fullatom_    = false;
	bJumps_use_IntraResStub_ = false;
	nres  ( nres_in );
	resize( nres_in );
	symminfo_ = new core::conformation::symmetry::SymmetryInfo();
	symminfo_->set_use_symmetry(false);
}

BinaryProteinSilentStruct::BinaryProteinSilentStruct()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	fullatom_    = false;
	bJumps_use_IntraResStub_ = false;
	nres( 0 );
	decoy_tag( "empty" );
	symminfo_ = new core::conformation::symmetry::SymmetryInfo();
	symminfo_->set_use_symmetry(false);
}


BinaryProteinSilentStruct::BinaryProteinSilentStruct(
	core::pose::Pose const & pose,
	std::string tag
) {
	bJumps_use_IntraResStub_ = false;
	symminfo_ = new core::conformation::symmetry::SymmetryInfo();
	symminfo_->set_use_symmetry(false);
	fill_struct( pose, tag );
} // BinaryProteinSilentStruct

void
BinaryProteinSilentStruct::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {
	tr.Trace << "binary:fill_struct... " << std::endl;
	decoy_tag( tag );

	if ( tag == "empty_tag" ) set_tag_from_pose( pose );

	fullatom( pose.is_fullatom() );
	tr.Trace << "get energies from pose..." << std::endl;
	energies_from_pose( pose );

	// conformation information
	sequence( pose.annotated_sequence( true ) );

	if( !core::pose::symmetry::is_symmetric(pose) ) {
		resize( pose.total_residue() );
	} else {    // core::pose::symmetry::is_symmetric(pose)
		//fpd previous implementation stored all atom coords and only wrote asymm unit coords
		//fpd however, if this struct is doubling as temporary storage for poses (as in batch relax)
		//fpd    we only want to store asymm unit coords
		//fpd in other words, the code should work if we call fill_struct()->fill_pose() without an intervening print_conformation()
		symmetry_info( *core::pose::symmetry::symmetry_info(pose)->clone() );
		core::Size nres_effective = symmetry_info()->num_virtuals() + symmetry_info()->num_independent_residues();
		resize( nres_effective );
	}

	tr.Trace << "read coords..." << std::endl;
	for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
		core::conformation::Residue const& resi = pose.residue(i);
		//int natoms = pose.residue(i).natoms();
		if( is_symmetric() && !symmetry_info()->bb_is_independent( i ) ) continue;
		int i_asymm =  symmetry_info()->get_asymmetric_seqpos( i ); // remaps virtual ids

		atm_coords_[i_asymm].resize( resi.natoms() );
		for (unsigned int j = 1; j <= resi.natoms(); ++j) {
			atm_coords_[i_asymm][j] = resi.atom(j).xyz();
		}
		secstruct_[i_asymm] = pose.secstruct(i);
	} // for ( unsigned int i = 1; i <= pose.total_residue(); ++i )

	fold_tree_ = pose.fold_tree();
	jumps_.clear();
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
		add_jump( pose.jump(nr) );
	}

	chain_endings( pose.conformation().chain_endings() );

	fill_struct_with_residue_numbers( pose ); // grabs residue numbers from pose PDBInfo object.

} // BinaryProteinSilentStruct


void BinaryProteinSilentStruct::add_chain_ending( Size const seqpos ) {
	core::Size nres_pose = nres();
	if ( is_symmetric() )
		 nres_pose = symmetry_info()->num_total_residues_with_pseudo() ;
	if ( seqpos < 1 || seqpos >= nres_pose ) {
		tr.Fatal << "ERROR: add_chain_ending() invalid chain ending " << seqpos << std::endl;
		utility_exit();
	}

	chain_endings_.push_back( seqpos );
	std::sort( chain_endings_.begin(), chain_endings_.end() ); // keep the list sorted
}

void BinaryProteinSilentStruct::parse_chain_endings( std::istream & stream ) {
	std::string s;
	stream >> s; // first column is "CHAIN_ENDINGS" tag, skip it

	utility::vector1< std::string > v;
	while ( stream.good() ) {
		stream >> s;
		v.push_back( s );
	}
	// remember to skip the last entry in the vector, which is the structure's nametag
	for ( Size i = 1, ie = v.size(); i < ie; ++i ) {
		add_chain_ending( boost::lexical_cast< Size >( v[ i ] ) );
	}
}

std::string BinaryProteinSilentStruct::chain_endings_str() const {
	std::ostringstream ss;
	ss << "CHAIN_ENDINGS ";

	for ( utility::vector1< Size >::const_iterator i = chain_endings().begin(), ie = chain_endings().end(); i != ie; ++i ) {
		ss << ' ' << (*i);
	}

	return ss.str();
}

void BinaryProteinSilentStruct::chain_endings( utility::vector1< Size > const & endings ) {
	core::Size nres_pose = nres();
	if ( is_symmetric() )
		 nres_pose = symmetry_info()->num_total_residues_with_pseudo() ;
	for ( utility::vector1< Size >::const_iterator i = endings.begin(), ie = endings.end(); i != ie; ++i ) {
		if ( (*i) < 1 || (*i) > nres_pose ) {  //fpd if symmetric, chainendings may be > nres (in asu)
			tr.Fatal << "ERROR: chain_endings() invalid chain ending " << (*i) << std::endl;
			utility_exit();
		}
	}

	chain_endings_ = endings;
	std::sort( chain_endings_.begin(), chain_endings_.end() ); // keep the list sorted
}

bool BinaryProteinSilentStruct::init_from_lines(
	utility::vector1< std::string > const & lines,
	SilentFileData & container
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< std::string > energy_names_;
	utility::vector1< std::string >::const_iterator iter = lines.begin();

	if ( iter->substr(0,9) != "SEQUENCE:" ) {
		// get sequence and scorename data from the silent-file data object,
		// because I don't have it!
		EnergyNamesOP enames = EnergyNamesOP(
			static_cast< EnergyNames * >
			( container.get_shared_silent_data( energynames )() )
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
		++iter;
		if ( iter == lines.end() ) {
			utility_exit_with_message( "While reading binary silent structure, encountered end of structure too early after reading the 'SEQUENCE' line" );
		}
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


	core::Size currpos = 1;
	bool bitflip = false;
	fullatom_ = false; //start with fullatom_ false and update to true as soon as a residue with too many atoms is read...
	bool fullatom_well_defined = false;
	for ( utility::vector1< std::string >::const_iterator end = lines.end();
				iter != end; ++iter ) {
		std::string tag;
		std::istringstream line_stream( *iter );

		if ( iter->substr(0,6) == "REMARK" ) {
			std::string tag;
			std::string comment;
			std::string value;
			line_stream >> tag >> comment >> value;
			add_comment( comment, value );
			continue;  // skip comments
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
			// parse fold_tree and jump lines
			if ( iter->substr(0,10) == "FOLD_TREE " ) {
				kinematics::FoldTree f;
				line_stream >> f;
				fold_tree( f ); // add fold-tree to this SilentStruct
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
				tr.Warning << "Skipping duplicate sequence declaration " << std::endl;
				continue;
			} else if ( iter->substr(0,19) == "ANNOTATED_SEQUENCE:" ) {
				std::string annotated_seq;
				line_stream >> tag; //ANNOTATED_SEQUENCE
				line_stream >> annotated_seq;
				sequence( annotated_seq );
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
			} else if ( iter->substr(0,13) == "SYMMETRY_INFO" ) { // symmetry information
				core::conformation::symmetry::SymmetryInfo s;
				line_stream >> s;
				symmetry_info( s );
				core::Size nres_effective = symmetry_info()->num_virtuals() + symmetry_info()->num_independent_residues();
				resize( nres_effective );   //fpd resize containers to ASU
				continue;
			} else if ( iter->substr( 0, 13 ) == "CHAIN_ENDINGS" ) {
				chain_endings_.clear();
				parse_chain_endings( line_stream );
				continue;
			} else if ( iter->substr(0,7) == "RES_NUM" ) {
				figure_out_residue_numbers_from_line( line_stream );
				continue;
			}

			// parse coords
			line_stream >> tag;

			if ( tag.length() < 1 ) {
				tr.Warning << "WARNING:  read blank line in decoy tag " << decoy_tag()
					<< std::endl;
				continue;
			}
			if ( static_cast< core::Size > (currpos) > secstruct_.size() ) {
				tr.Error << "ERROR: trying to index off the end of the secstruct array"
					<< ", idx is " << currpos << " secstruct_.size() is " << secstruct_.size()
					<< std::endl;
				return false;
			}
			secstruct_[currpos] = tag[0];  // first char is sec struct

			int natoms = (tag.length()-1) / 16;
			utility::vector1< numeric::xyzVector <float> > atm_buff( natoms+1 );
			utility::decode6bit( (unsigned char*)&(atm_buff[1]) , tag.substr(1) );

			// endianness check ...
			//   check the dist between atoms 1 and 2 as well as atoms 2 and 3 if
			//   EITHER is unreasonable .. and flipping fixes BOTH then turn bitflip on
			if (currpos == 1) {
				core::Real len_check12 = (atm_buff[1]-atm_buff[2]).length();
				core::Real len_check23 = (atm_buff[3]-atm_buff[2]).length();
				if ( len_check12 < 0.5 || len_check12 > 2.0 || len_check23 < 0.5 ||
							len_check23 > 2.0
				) {
					utility::swap4_aligned ( (void*) &(atm_buff[1][0]) , 3*natoms );
						// recheck; if not better flip back
					len_check12 = (atm_buff[1]-atm_buff[2]).length();
					len_check23 = (atm_buff[3]-atm_buff[2]).length();
					if ( len_check12 < 0.5 || len_check12 > 2.0 || len_check23 < 0.5 || len_check23 > 2.0 ) {
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
			if ( !symmetry_info()->get_use_symmetry() || currpos <=symmetry_info()->num_independent_residues()  ) {
				// always run this if we're not dealing with a symmetric pose, or, if we're dealing with a symmetric pose
				// and we're still reading in data for the assymetric unit.  But if we're reading in a symmetric pose
				// and currpos > the number of residues in the asymmetric unit (i.e. we're reading in a virtual residue),
				// then DON'T use a count of the number of atoms in the given residue as an indication of whether we're dealing
				// with a fullatom structure.
				detect_fullatom( currpos, natoms, fullatom_, fullatom_well_defined );
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

	//	if ( (unsigned int) currpos != nres() + 1 ) 		return false;
	//	if ( nres() == 0 ) 		return false;


	if ( atm_coords_.size() != nres() ) {
		tr << "ERROR: didn't find coordinates for all sequence positions of "
							<< decoy_tag() << std::endl;
		tr << "       expected " << nres()
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

	//tr.Debug << "(TEX) FOLD TREE: " << fold_tree();
	if ( !fullatom_well_defined ) fullatom_ = option[ in::file::fullatom ]();
	return true;
} // init_from_lines

/// @brief Resize this silent-struct to the appropriate number of residues.
void
BinaryProteinSilentStruct::resize(
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

void BinaryProteinSilentStruct::fill_pose (
	core::pose::Pose & pose
) const {
	using namespace core::chemical;
	ResidueTypeSetCAP residue_set;
	tr.Debug << "fill_pose: SilentStruct is " << ( fullatom() ? "fullatom" : "centroid" ) << std::endl;
	if ( fullatom() ) {
		residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	} else {
		residue_set = ChemicalManager::get_instance()->residue_type_set( CENTROID );
	}
	fill_pose( pose, *residue_set );
} // fill_pose

void BinaryProteinSilentStruct::fill_pose (
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set
) const {

	runtime_assert( nres() != 0 );
	runtime_assert( sequence() != "" );
	if ( pose.annotated_sequence() != sequence()
			|| fullatom() != pose.is_fullatom() ) {  //fpd
		//fpd  this function assumes an asymmetric conformation to start (which will later be symmetrized)
		if (core::pose::symmetry::is_symmetric(pose))
			core::pose::symmetry::make_asymmetric_pose( pose );
		core::pose::make_pose_from_sequence( pose, sequence(), residue_set, false /*auto_termini*/ );
	}

	//fpd ???
	pose.energies().clear();

	// coords
	utility::vector1< id::AtomID > atm_ids;
	utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;
	core::Size nres_pose = nres();
	if ( is_symmetric() )
		 nres_pose = symmetry_info()->num_total_residues_with_pseudo() ;
	for ( Size seqpos = 1; seqpos <= nres_pose; ++seqpos ) {
		int natoms_pose = pose.residue(seqpos).natoms() ;
		int atm_seqpos =  symmetry_info()->get_asymmetric_seqpos( seqpos ) ;
		int natoms_struct = atm_coords_[atm_seqpos].size();
		int natoms_total = std::min( natoms_pose , natoms_struct );

		if ( natoms_pose != natoms_struct) {
			tr.Warning 	<< "[ WARNING ] "
				<< "Number of atoms in pose and silent file disagree! "
				<< "Attempting to continue ..." << std::endl
				<< "[ WARNING ]    (in residue " << seqpos
				<< "  natoms_pose=" << natoms_pose
				<< "atm_seqpos " << atm_seqpos << "  natoms_struct="
				<< natoms_struct << ")" << std::endl;
			tr.flush();
		}

		//natoms_pose  = pose.residue(seqpos).natoms();
		//natoms_total = std::min( natoms_pose, natoms_struct );

		for ( int j = 1; j <= natoms_total; ++j ){
			id::AtomID id( j, seqpos );
			//numeric::xyzVector< core::Real> atom_i(atm_coords_[atm_seqpos][j][0], atm_coords_[atm_seqpos][j][1], atm_coords_[atm_seqpos][j][2]);
			//pose.set_xyz( id, atom_i );
			atm_ids.push_back( id );
			atm_xyzs.push_back(
			     numeric::xyzVector< core::Real>(atm_coords_[atm_seqpos][j][0], atm_coords_[atm_seqpos][j][1], atm_coords_[atm_seqpos][j][2]) );
		}
		pose.set_secstruct( seqpos, secstruct_[atm_seqpos] );
	} // for ( seqpos )
	pose.batch_set_xyz( atm_ids, atm_xyzs );

	tr.Debug << "FOLD TREE: " << fold_tree();
	// set fold_tree
	pose.fold_tree( fold_tree() );

	// SYMMETRY - setup up Symmetry stuff here.
	//fpd  As a note, this doesn't symmetrize the conformation like calling the constructor with symmdata would
	//fpd  Instead, this just adds symminfo to an already-symmetrized pose
	//fpd  When jumps are symmetrized the pose will also be symmetric
	if( is_symmetric() ) {
		core::pose::symmetry::make_symmetric_pose( pose, *symmetry_info() );
	}

	// set jumps
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
		if( is_symmetric() && !symmetry_info()->jump_is_independent(nr) ) continue;
		if ( !bJumps_use_IntraResStub_ ) { //default modern file-format
			pose.set_jump( nr, jump( nr ) );
		} else { //support for rosetta++ format
			Size start = fold_tree().jump_edge( nr ).start();
			Size stop = fold_tree().jump_edge( nr ).stop();
			id::StubID up_stub  ( core::pose::named_stub_id_to_stub_id(id::NamedStubID( "CA", "N", "CA", "C", start < stop ? start : stop ), pose ) );
			id::StubID down_stub( core::pose::named_stub_id_to_stub_id(id::NamedStubID( "CA", "N", "CA", "C", start < stop ? stop : start ), pose ) );
			pose.conformation().set_stub_transform( up_stub, down_stub, jump( nr ) );
		}
	}


	// fpd if symmetric 'nres' covers the ASU while 'sequence' covers the symmetric pose
	tr.Debug << "nres = " << nres_pose << std::endl;
	tr.Debug << "one_letter_sequence() = " << one_letter_sequence().length() << std::endl;

	if( nres_pose != one_letter_sequence().length() ){
		utility_exit_with_message( "RuntimeAssert failed: nres() == one_letter_sequence().length()" );
	}


  if ( !chain_endings().empty() ) {
    pose.conformation().chain_endings( chain_endings() );
  }


	core::pose::initialize_disulfide_bonds(pose);

	finish_pose( pose );
} // fill_pose

void
BinaryProteinSilentStruct::print_header( std::ostream & out ) const
{
	SilentStruct::print_header( out );
}


void BinaryProteinSilentStruct::print_conformation(
	std::ostream & output
) const {
	output << "REMARK BINARY SILENTFILE\n";

	// fold tree
	// assume non-trivial fold_tree only if more than one edge, i.e., EDGE 1 <nres> -1?
	// no -- can have a fold tree with a single jump, actually.
	if ( fold_tree().size() > 1 || fold_tree().num_jump() > 1 ) {
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

	// sequence
	output << "ANNOTATED_SEQUENCE: " << sequence() << " " << decoy_tag() << "\n"; //chu print annotated_sequence per decoy

	//lin print out the SYMMETRY_INFO line here
	if( is_symmetric() ) {
		output << *symmetry_info() << ' ' << decoy_tag() << "\n";
	}

	//tr.Debug << "FOLD_TREE Size: " << fold_tree().size() << " " << fold_tree() << std::endl;

	// chain endings
	if ( !chain_endings().empty() ) {
		output << chain_endings_str() << ' ' << decoy_tag() << '\n';
	}

	// fullatom flag
	//int fullatom_flag = (fullatom_? 1 : 2);
	std::string resline;
	//encode6bit(  (unsigned char*)&fullatom_flag, 4, resline );
	//output << resline << "\n";

	using namespace basic::options;
	if( option[ OptionKeys::run::write_failures ].user() && option[ OptionKeys::run::write_failures ] == false  && decoy_tag().substr( 0, 8 ) == "FAILURE_" ) return;
	// the encoded coordinates
	for ( Size i = 1; i <= nres(); ++i ) {
		//fpd  symmetric residues are now trimmed in fill_struct()
		//if( !symmetry_info()->is_asymmetric_seqpos(i) ) continue; //Symmetry only output asymmetric unit

		char this_secstr = secstruct_[i];
		if (this_secstr < 'A' || this_secstr > 'Z') this_secstr = 'L';

		utility::encode6bit(  (unsigned char*)&atm_coords_[i][1][0], atm_coords_[i].size()*12, resline );  // ASSUMES FLOAT == 4 BYTES!!! (eep!)
		output << this_secstr << resline << ' ' << decoy_tag() << "\n";
	} // for ( Size i = 1; i <= nres; ++i )
} // print_conformation


Real BinaryProteinSilentStruct::get_debug_rmsd() {
	pose::Pose temp_pose;
	ObjexxFCL::FArray2D< Real > rebuilt_coords ( 3, atm_coords_.size() ),
		original_coords( 3, atm_coords_.size() );

	// build temp_pose from coordinates
	fill_pose( temp_pose );

	for ( Size i = 1; i <= temp_pose.total_residue(); ++i ) {
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			rebuilt_coords (k,i) = temp_pose.residue(i).xyz( "CA" )[k-1];
			original_coords(k,i) = atm_coords_[i][2][k-1];
		}
	}

	Real rmsd = numeric::model_quality::rms_wrapper( temp_pose.total_residue(), rebuilt_coords, original_coords );
	return rmsd;
}

ObjexxFCL::FArray2D< Real >
BinaryProteinSilentStruct::get_CA_xyz() const {
	core::Size n_residues = nres();
	ObjexxFCL::FArray2D< Real > my_coords( 3, n_residues );
	for ( Size i = 1; i <= n_residues; ++i ) { // i = n_residues
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			my_coords(k,i) = atm_coords_[i][2][k-1];
		} // k
	} // i

	return my_coords;
} // get_CA_positions

Real BinaryProteinSilentStruct::CA_rmsd( BinaryProteinSilentStruct other_pss ) {
	ObjexxFCL::FArray2D< Real > my_coords    = get_CA_xyz();
	ObjexxFCL::FArray2D< Real > other_coords = other_pss.get_CA_xyz();
	Real rmsd = numeric::model_quality::rms_wrapper( nres(), my_coords, other_coords );

	return rmsd;
}

BinaryProteinSilentStruct & BinaryProteinSilentStruct::operator= (
	BinaryProteinSilentStruct const &
)
{
	utility_exit_with_message( "called BinaryProteinSilentStruct::operator=)" );
	exit(1);
}

} // namespace silent
} // namespace io
} // namespace core
