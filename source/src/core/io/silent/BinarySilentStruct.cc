// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/BinarySilentStruct.cc
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

#include <core/scoring/Energies.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/raw_data/DisulfideFile.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedStubID.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <core/conformation/Residue.hh>

#include <numeric/model_quality/rms.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <boost/lexical_cast.hpp>

#include <utility/vector1.hh>
#include <utility/Binary_Util.hh>
#include <utility/pointer/owning_ptr.hh>


static THREAD_LOCAL basic::Tracer tr( "core.io.silent" );

namespace core {
namespace io {
namespace silent {

/// @brief Constructors.
BinarySilentStruct::BinarySilentStruct( Size const nres_in )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	fullatom_    = false;
	bJumps_use_IntraResStub_ = false;
	nres  ( nres_in );
	resize( nres_in );
	symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo() );
	symminfo_->set_use_symmetry(false);
}

BinarySilentStruct::BinarySilentStruct()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	fullatom_    = false;
	bJumps_use_IntraResStub_ = false;
	nres( 0 );
	decoy_tag( "empty" );
	symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo() );
	symminfo_->set_use_symmetry(false);
}


BinarySilentStruct::BinarySilentStruct(
	core::pose::Pose const & pose,
	std::string tag
) {
	bJumps_use_IntraResStub_ = false;
	symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo() );
	symminfo_->set_use_symmetry(false);
	fill_struct( pose, tag );
} // BinarySilentStruct

BinarySilentStruct::~BinarySilentStruct()
{
}

void
BinarySilentStruct::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {

	Parent::fill_struct( pose, tag );

	fullatom( pose.is_fullatom() );

	if ( !core::pose::symmetry::is_symmetric(pose) ) {
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
		if ( is_symmetric() && !symmetry_info()->bb_is_independent( i ) ) continue;
		int i_asymm =  symmetry_info()->get_asymmetric_seqpos( i ); // remaps virtual ids

		atm_coords_[i_asymm].resize( resi.natoms() );
		for ( unsigned int j = 1; j <= resi.natoms(); ++j ) {
			atm_coords_[i_asymm][j] = resi.atom(j).xyz();
		}
		secstruct_[i_asymm] = pose.secstruct(i);
	} // for ( unsigned int i = 1; i <= pose.total_residue(); ++i )

	noncanonical_residue_connections_.clear();
	for ( Size ires=1; ires<=pose.total_residue(); ++ires ) {
		for ( int icon=1; icon<=(int)pose.residue_type(ires).n_possible_residue_connections(); ++icon ) {
			Size atom_num = pose.residue(ires).residue_connection( icon ).atomno();
			core::id::NamedAtomID atom1(pose.residue_type(ires).atom_name(atom_num), ires);
			if ( pose.residue(ires).connected_residue_at_resconn(icon) != 0 ) {
				Size anchor_rsd = pose.residue(ires).connected_residue_at_resconn(icon);
				if ( anchor_rsd > ires ) continue;
				int anchor_conid = (int) pose.residue(ires).connect_map(icon).connid();
				Size anchor_atom_num = pose.residue(anchor_rsd).residue_connection( anchor_conid ).atomno();
				core::id::NamedAtomID atom2(pose.residue_type(anchor_rsd).atom_name(anchor_atom_num), anchor_rsd);

				if ( ires == anchor_rsd+1 ) {
					if ( icon == (int) pose.residue_type(ires).lower_connect_id() ) {
						if ( anchor_conid == (int) pose.residue_type(anchor_rsd).upper_connect_id() ) {
							// canonical polymer connection, skip
							continue;
						}
					}
				}
				noncanonical_residue_connections_.push_back(std::pair<core::id::NamedAtomID, core::id::NamedAtomID> (atom1, atom2) );
			}
		}
	}

	fold_tree_ = pose.fold_tree();
	jumps_.clear();
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++ )  {
		add_jump( pose.jump(nr) );
	}

	chain_endings( pose.conformation().chain_endings() );

	fill_struct_with_residue_numbers( pose ); // grabs residue numbers from pose PDBInfo object.
	fill_struct_with_submotif_info_list( pose );
	fill_other_struct_list( pose );

} // BinarySilentStruct

void BinarySilentStruct::add_chain_ending( Size const seqpos ) {
	core::Size nres_pose = nres();
	if ( is_symmetric() ) {
		nres_pose = symmetry_info()->num_total_residues_with_pseudo() ;
	}
	if ( seqpos < 1 || seqpos >= nres_pose ) {
		tr.Fatal << "ERROR: add_chain_ending() invalid chain ending " << seqpos << std::endl;
		utility_exit();
	}

	chain_endings_.push_back( seqpos );
	std::sort( chain_endings_.begin(), chain_endings_.end() ); // keep the list sorted
}

void BinarySilentStruct::parse_chain_endings( std::istream & stream ) {
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

std::string BinarySilentStruct::chain_endings_str() const {
	std::ostringstream ss;
	ss << "CHAIN_ENDINGS ";

	for ( utility::vector1< Size >::const_iterator i = chain_endings().begin(), ie = chain_endings().end(); i != ie; ++i ) {
		ss << ' ' << (*i);
	}

	return ss.str();
}

void BinarySilentStruct::chain_endings( utility::vector1< Size > const & endings ) {
	core::Size nres_pose = nres();
	if ( is_symmetric() ) {
		nres_pose = symmetry_info()->num_total_residues_with_pseudo() ;
	}
	for ( utility::vector1< Size >::const_iterator i = endings.begin(), ie = endings.end(); i != ie; ++i ) {
		if ( (*i) < 1 || (*i) > nres_pose ) {  //fpd if symmetric, chainendings may be > nres (in asu)
			tr.Fatal << "ERROR: chain_endings() invalid chain ending " << (*i) << std::endl;
			utility_exit();
		}
	}

	chain_endings_ = endings;
	std::sort( chain_endings_.begin(), chain_endings_.end() ); // keep the list sorted
}

bool BinarySilentStruct::init_from_lines(
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
			utility::pointer::static_pointer_cast< core::io::silent::EnergyNames > ( container.get_shared_silent_data( energynames ) )
		);

		SimpleSequenceDataOP seqdata = SimpleSequenceDataOP(
			utility::pointer::static_pointer_cast< core::io::silent::SimpleSequenceData > ( container.get_shared_silent_data( simplesequencedata ) )
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
		resize( temp_seq.size() );

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
		if ( iter->substr(0,7) == "SCORE: " || iter->substr(0,7) == "OTHER: " ) {
			// SCORE: line with values from this structure.
			Size nres = one_letter_sequence().length();
			resize( nres );

			std::string tag;
			line_stream >> tag;
			if ( line_stream.fail() || ( tag != "SCORE:" && tag != "OTHER:" ) ) {
				tr.Error << "bad format in first score line of silent file" << std::endl;
				tr.Error << "line = " << *iter << std::endl;
				tr.Error << "tag = " << tag << std::endl;
			}
			scoreline_prefix( tag );

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
			} else if ( iter->substr(0,24) == "NONCANONICAL_CONNECTION:" ) {
				line_stream >> tag;
				int res1, res2;
				std::string atom1, atom2;
				line_stream >> res1 >> atom1 >> res2 >> atom2 ;
				core::id::NamedAtomID atom_id1(atom1, res1);
				core::id::NamedAtomID atom_id2(atom2, res2);
				noncanonical_residue_connections_.push_back( std::pair< core::id::NamedAtomID, core::id::NamedAtomID > (atom_id1, atom_id2) ) ;
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
			} else if ( iter->substr(0,11) == "SEGMENT_IDS" ) {
				figure_out_segment_ids_from_line( line_stream );
				continue;
			} else if ( iter->substr(0,13) == "SUBMOTIF_INFO" ) {
				add_submotif_info_from_line( line_stream );
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

			// option to force bit flip:
			if ( force_bitflip() ) {
				tr.Warning << "Forcing binary silent file big-endian/little-endian flipping!  Tag: " << decoy_tag() << std::endl;
				bitflip=true;
			}

			// endianness check ...
			//   check the dist between atoms 1 and 2 as well as atoms 2 and 3 if
			//   EITHER is unreasonable .. and flipping fixes BOTH then turn bitflip on
			// [ hey wait, len23 assumed protein -- only look at atoms 1 and 2 -- rhiju ]
			if ( currpos == 1 && !bitflip ) {
				core::Real len_check12 = (atm_buff[1]-atm_buff[2]).length();
				if ( len_check12 < 0.5 || len_check12 > 2.0 || !numeric::is_a_finitenumber( len_check12, 1.0, 0.0 ) ) {
					utility::swap4_aligned ( (void*) &(atm_buff[1][0]) , 3*natoms );

					// recheck; if not better flip back
					len_check12 = (atm_buff[1]-atm_buff[2]).length();

					if ( len_check12 < 0.5 || len_check12 > 2.0 || !numeric::is_a_finitenumber( len_check12, 1.0, 0.0 ) ) { //|| len_check23 < 0.5 || len_check23 > 2.0 ) {
						utility::swap4_aligned ( (void*) &(atm_buff[1][0]) , 3*natoms );
					} else {
						tr.Warning << "Reading binary silent file with inverted endian-ness!  Will attempt to flip automatically. Tag: " << decoy_tag() << std::endl;
						bitflip = true;
					}
				}
			} else {
				if ( bitflip ) {
					utility::swap4_aligned ( (void*) &(atm_buff[1][0]) , 3*natoms );
				}
			}
			if ( !symmetry_info()->get_use_symmetry() || currpos <=symmetry_info()->num_independent_residues()  ) {
				// always run this if we're not dealing with a symmetric pose,
				// or, if we're dealing with a symmetric pose and we're still
				// reading in data for the assymetric unit.  But if we're
				// reading in a symmetric pose and currpos > the number of
				// residues in the asymmetric unit (i.e. we're reading in a
				// virtual residue), then DON'T use a count of the number of
				// atoms in the given residue as an indication of whether we're
				// dealing with a fullatom structure.
				detect_fullatom( currpos, natoms, fullatom_, fullatom_well_defined );
			}

			atm_coords_[currpos].resize( natoms ); // allocate space for coords
			for ( int j=1; j<=natoms; ++j ) {
				atm_coords_[currpos][j] = atm_buff[j];
			}
			currpos++;
			//tr.Debug << "processing line " << *iter << std::endl;
		} // conformation lines
	} // for ( iter ... )

	if ( fold_tree().num_jump() != jumps_.size() ) {
		tr.Warning << "parse error:  found " << jumps_.size()
			<< " RT lines for a fold-tree with " << fold_tree().num_jump()
			<< " for decoy tag " << decoy_tag() << std::endl;
		return false;
	}

	// if ( (unsigned int) currpos != nres() + 1 )   return false;
	// if ( nres() == 0 )   return false;


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
BinarySilentStruct::resize(
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

void BinarySilentStruct::fill_pose (
	core::pose::Pose & pose
) const {
	using namespace core::chemical;
	ResidueTypeSetCOP residue_set;
	tr.Debug << "fill_pose: SilentStruct is " << ( fullatom() ? "fullatom" : "centroid" ) << std::endl;
	if ( fullatom() ) {
		residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	} else {
		residue_set = ChemicalManager::get_instance()->residue_type_set( CENTROID );
	}
	fill_pose( pose, *residue_set );
} // fill_pose

void BinarySilentStruct::fill_pose (
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set
) const {

	runtime_assert( nres() != 0 );
	runtime_assert( sequence() != "" );
	if ( pose.annotated_sequence() != sequence()
			|| fullatom() != pose.is_fullatom() ) {  //fpd
		//fpd  this function assumes an asymmetric conformation to start (which will later be symmetrized)
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			core::pose::symmetry::make_asymmetric_pose( pose );
		}
		core::pose::make_pose_from_sequence( pose, sequence(), residue_set, false /*auto_termini*/ );
	}

	for ( Size i_conn=1; i_conn <= noncanonical_residue_connections_.size(); ++i_conn ) {
		Size res1 = noncanonical_residue_connections_[i_conn].first.rsd();
		std::string atom1 = noncanonical_residue_connections_[i_conn].first.atom();
		Size res2 = noncanonical_residue_connections_[i_conn].second.rsd();
		std::string atom2 = noncanonical_residue_connections_[i_conn].second.atom();

		pose.conformation().declare_chemical_bond(res1, atom1, res2, atom2);
	}

	//fpd ???
	pose.energies().clear();

	// coords
	utility::vector1< id::AtomID > atm_ids;
	utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;
	core::Size nres_pose = nres();
	if ( is_symmetric() ) {
		nres_pose = symmetry_info()->num_total_residues_with_pseudo() ;
	}
	for ( Size seqpos = 1; seqpos <= nres_pose; ++seqpos ) {
		core::conformation::Residue const & rsd = pose.residue(seqpos);
		Size natoms_pose = rsd.natoms() ;
		Size atm_seqpos =  symmetry_info()->get_asymmetric_seqpos( seqpos ) ;
		Size natoms_struct = atm_coords_[atm_seqpos].size();
		bool fixup( false );
		if ( natoms_pose != natoms_struct ) {

			// for backwards compatibility -- sorry for hack. If this
			// hack needs to be expanded, would be best to move into separate function.
			// N_ACETYLATION & C_METHYLAMIDATION lost 3 hydrogen atoms.
			if ( (natoms_pose + 3) == natoms_struct && rsd.has_variant_type( chemical::C_METHYLAMIDATION ) ) fixup = true;
			if ( (natoms_pose + 3) == natoms_struct && rsd.has_variant_type( chemical::N_ACETYLATION ) ) fixup = true;

			if ( !fixup )  {
				tr.Warning  << tr.Red << "[ WARNING ] "
					<< "Number of atoms in pose and silent file disagree! "
					<< "Attempting to continue ..." << std::endl
					<< "[ WARNING ]    (in residue " << pose.residue(seqpos).name() << " at " << seqpos
					<< "  natoms_pose=" << natoms_pose
					<< "atm_seqpos " << atm_seqpos << "  natoms_struct="
					<< natoms_struct << ")" << tr.Reset << std::endl;
				tr.flush();
			}
		}

		Size count( 0 );

		// fpd
		Size natoms_total = natoms_struct;
		if ( !fixup ) natoms_total = std::min( natoms_struct, natoms_pose );

		for ( Size j = 1; j <= natoms_total; ++j ) {

			// hydrogens that were present in old rosetta, now removed.
			if ( fixup && rsd.has_variant_type( chemical::C_METHYLAMIDATION ) && (j > rsd.nheavyatoms()+3) && (j < rsd.nheavyatoms()+7 ) ) continue;
			if ( fixup && rsd.has_variant_type( chemical::N_ACETYLATION ) && (j > rsd.nheavyatoms()+2) && (j < rsd.nheavyatoms()+6 ) ) continue;

			id::AtomID id( ++count, seqpos );
			atm_ids.push_back( id );
			atm_xyzs.push_back(
				numeric::xyzVector< core::Real>(atm_coords_[atm_seqpos][j][0],
				atm_coords_[atm_seqpos][j][1],
				atm_coords_[atm_seqpos][j][2]) );
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
	if ( is_symmetric() ) {
		core::pose::symmetry::make_symmetric_pose( pose, *symmetry_info() );
	}

	// set jumps
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++ )  {
		if ( is_symmetric() && !symmetry_info()->jump_is_independent(nr) ) continue;
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

	// additional setup for mirrored poses
	if ( core::pose::symmetry::is_mirror_symmetric( pose ) ) {
		conformation::symmetry::MirrorSymmetricConformation & mirror_conf(
			dynamic_cast< conformation::symmetry::MirrorSymmetricConformation & >( pose.conformation() ) );
		mirror_conf.update_residue_identities();
	}

	// fpd if symmetric 'nres' covers the ASU while 'sequence' covers the symmetric pose
	tr.Debug << "nres = " << nres_pose << std::endl;
	tr.Debug << "one_letter_sequence() = " << one_letter_sequence().length() << std::endl;

	if ( nres_pose != one_letter_sequence().length() ) {
		utility_exit_with_message( "RuntimeAssert failed: nres() == one_letter_sequence().length()" );
	}

	if ( !chain_endings().empty() ) {
		pose.conformation().chain_endings( chain_endings() );
	}

	core::pose::initialize_disulfide_bonds(pose);

	finish_pose( pose );

	setup_other_poses( pose, residue_set );

} // fill_pose

/// @details  for stepwise modeling, setup other_poses inside full_model_info. this could go into SilentStruct parent class, if needed.
void
BinarySilentStruct::setup_other_poses( pose::Pose & pose, core::chemical::ResidueTypeSet const & residue_set ) const {
	for ( Size n = 1; n <= other_struct_list().size(); n++ ) {
		core::pose::PoseOP other_pose( new core::pose::Pose );
		other_struct_list()[n]->fill_pose( *other_pose, residue_set );
		core::pose::full_model_info::nonconst_full_model_info( pose ).add_other_pose( other_pose );
	}
}

void
BinarySilentStruct::print_header( std::ostream & out ) const
{
	SilentStruct::print_header( out );
	out << "REMARK BINARY SILENTFILE";
	if ( full_model_parameters() != 0 ) out << "    " << *full_model_parameters();
	out << std::endl;
}


void BinarySilentStruct::print_conformation(
	std::ostream & output
) const {

	// fold tree
	// assume non-trivial fold_tree only if more than one edge, i.e., EDGE 1 <nres> -1?
	// no -- can have a fold tree with a single jump, actually.
	if ( fold_tree().size() > 1 || fold_tree().num_jump() > 0 ) {
		output << "FOLD_TREE ";
		for ( kinematics::FoldTree::const_iterator
				it = fold_tree().begin(), it_end = fold_tree().end();
				it != it_end; ++it
				) {
			output << *it;
		}
		//  output << fold_tree(); this produces a new-line --- wrong behaviour
		//  of fold_tree but I don't want to fix 1000 u-tracer unit-tests!
		output << ' ' << decoy_tag() << "\n";
	}
	for ( Size i = 1; i <= fold_tree().num_jump(); i++ ) {
		output << jump( i ) << ' ' << decoy_tag() << "\n";
	}

	// sequence
	output << "ANNOTATED_SEQUENCE: " << sequence() << " " << decoy_tag() << "\n"; //chu print annotated_sequence per decoy

	//lin print out the SYMMETRY_INFO line here
	if ( is_symmetric() ) {
		output << *symmetry_info() << ' ' << decoy_tag() << "\n";
	}

	//tr.Debug << "FOLD_TREE Size: " << fold_tree().size() << " " << fold_tree() << std::endl;

	for ( Size i_conn=1; i_conn<=noncanonical_residue_connections_.size(); ++i_conn ) {
		output << "NONCANONICAL_CONNECTION: " << noncanonical_residue_connections_[i_conn].first.rsd() << " " << noncanonical_residue_connections_[i_conn].first.atom() << " " << noncanonical_residue_connections_[i_conn].second.rsd() << " " << noncanonical_residue_connections_[i_conn].second.atom() << std::endl;
	}

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
	if ( option[ OptionKeys::run::write_failures ].user() && option[ OptionKeys::run::write_failures ] == false  && decoy_tag().substr( 0, 8 ) == "FAILURE_" ) return;
	// the encoded coordinates
	for ( Size i = 1; i <= nres(); ++i ) {
		//fpd  symmetric residues are now trimmed in fill_struct()
		//if( !symmetry_info()->is_asymmetric_seqpos(i) ) continue; //Symmetry only output asymmetric unit

		char this_secstr = secstruct_[i];
		if ( this_secstr < 'A' || this_secstr > 'Z' ) this_secstr = 'L';

		utility::encode6bit(  (unsigned char*)&atm_coords_[i][1][0], atm_coords_[i].size()*12, resline );  // ASSUMES FLOAT == 4 BYTES!!! (eep!)
		output << this_secstr << resline << ' ' << decoy_tag() << "\n";
	} // for ( Size i = 1; i <= nres; ++i )
} // print_conformation


Real BinarySilentStruct::get_debug_rmsd() {
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
BinarySilentStruct::get_CA_xyz() const {
	core::Size n_residues = nres();
	ObjexxFCL::FArray2D< Real > my_coords( 3, n_residues );
	for ( Size i = 1; i <= n_residues; ++i ) { // i = n_residues
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			my_coords(k,i) = atm_coords_[i][2][k-1];
		} // k
	} // i

	return my_coords;
} // get_CA_positions

Real BinarySilentStruct::CA_rmsd( BinarySilentStruct other_pss ) {
	ObjexxFCL::FArray2D< Real > my_coords    = get_CA_xyz();
	ObjexxFCL::FArray2D< Real > other_coords = other_pss.get_CA_xyz();
	Real rmsd = numeric::model_quality::rms_wrapper( nres(), my_coords, other_coords );

	return rmsd;
}

BinarySilentStruct & BinarySilentStruct::operator= (
	BinarySilentStruct const &
)
{
	utility_exit_with_message( "called BinarySilentStruct::operator=)" );
	exit(1);
}

} // namespace silent
} // namespace io
} // namespace core
