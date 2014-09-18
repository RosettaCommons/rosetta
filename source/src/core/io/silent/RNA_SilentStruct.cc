// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/RNA_SilentStruct.cc
///
/// @brief Representation of rosetta++ protein silent-file structures.
/// @author James Thompson

// C++ Headers
#include <cmath>
#include <cstdlib>
// AUTO-REMOVED #include <fstream>
#include <iostream>
// AUTO-REMOVED #include <utility>
#include <vector>
// AUTO-REMOVED #include <list>
#include <string>
#include <map>
#include <sstream>

// mini headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/id/TorsionID.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/rna/util.hh>

#include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <numeric/model_quality/rms.hh>

#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>


namespace core {
namespace io {
namespace silent {

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

static thread_local basic::Tracer tr( "core.io.silent" );

/////////////////////////////////////////////////////////////////////////
// Following should be easy to generalize for protein, RNA, DNA.
// This may eventually be critical as we start to look at mixed systems.
// For now, just for safety and (perhaps) to avoid confusion,
// we'll go ahead and make this a separate class.
//   -- Rhiju, April 2008
/////////////////////////////////////////////////////////////////////////

RNA_SilentStruct::RNA_SilentStruct(
	core::pose::Pose const & pose,
	std::string tag,
	bool fa
) :
	fullatom_( fa )
{
	fill_struct( pose, tag );
} // RNA_SilentStruct

// RNA_SilentStruct::RNA_SilentStruct( core::io::silent::RNA_SilentStruct const & src )
// {
// 	nres( src.nres_ );
// 	fullatom_    = src.fullatom_;
// 	non_main_chain_sugar_coords_defined_ = src.non_main_chain_sugar_coords_defined_;
// 	resize( nres );
// 	for (Size n = 1; n <= nres; n++ ) {
// 		secstruct_[ n ] = src.secstruct_[ n ];
// 		for (Size k = 1; k <= src.mainchain_torsions_[n].size();k++ ) {
// 			mainchain_torsions_[n][k] = src.mainchain_torsions_[n][k];
// 		}
// 	}
// }

void
RNA_SilentStruct::fill_struct( core::pose::Pose const & pose, std::string tag ) {
	decoy_tag( tag );

	energies_from_pose( pose );

	// conformation information
	sequence( pose.sequence() );
	resize( pose.total_residue() );
	static const std::string important_atom = "C4'";
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		core::conformation::Residue resi = pose.residue(i);

		secstruct_[i] = pose.secstruct(i);
		// I wonder if we can just grab these torsions...
		mainchain_torsions_ [i] = resi.mainchain_torsions();
		coords_[i]    = resi.xyz( important_atom );
		if ( fullatom() ) {
			chi_torsions_[i] = resi.chi();
		} // if ( fullatom )

		//New (Feb. 2009)...
		// x-y-z of coordinates of C2', C1', and O4', in a local coordinate system defined
		// by C3', C4', and C5' (as "stub" atoms).
		{
			non_main_chain_sugar_coords_defined_ = true;
			kinematics::Stub const input_stub( resi.xyz( " C3'" ), resi.xyz( " C3'" ), resi.xyz( " C4'" ), resi.xyz( " C5'" ) );
			utility::vector1< Vector > vecs;
			for (Size n = 1; n <= chemical::rna::non_main_chain_sugar_atoms.size(); n++  ) {
				Vector v = input_stub.global2local( resi.xyz( chemical::rna::non_main_chain_sugar_atoms[ n ] ) );
				vecs.push_back( v );
			}
			non_main_chain_sugar_coords_[i] = vecs;
		}

	} // for ( Size i = 1; i <= pose.total_residue(); ++i )

	fold_tree_ = pose.fold_tree();
	jumps_.clear();
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
		add_jump( pose.jump(nr) );
	}

	fill_struct_with_residue_numbers( pose ); // grabs residue numbers from pose PDBInfo object.

} // RNA_SilentStruct

	//Following should be easy to generalize for protein vs. RNA.

bool RNA_SilentStruct::init_from_lines(
	utility::vector1< std::string > const & lines,
	SilentFileData & container
) {
	bool success( false );

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

		sequence( seqdata->sequence() );
		energy_names_ = enames ->energy_names();
	} else {
		// get sequence and scorename data from the first two lines provided, put them into container for further use
		// by other RNA_SilentStruct objects.

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
			return success;
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

	// resize myself appropriately, according to length of sequence
	resize( sequence().length() );

	for ( utility::vector1< std::string >::const_iterator end = lines.end(); iter != end;	++iter ) {
		std::string tag;
		std::istringstream line_stream( *iter );

		if ( iter->substr(0,7) == "SCORE: " ) { // SCORE: line with values from this structure.
			resize( sequence().length() ); // sequence_ should be defined by now.

			std::string tag;
			line_stream >> tag;
			if ( line_stream.fail() || tag != "SCORE:" ) {
				tr.Error << "bad format in first score line of silent file" << std::endl;
				tr.Error << "line = " << *iter << std::endl;
				tr.Error << "tag = " << tag << std::endl;
			}

			utility::vector1< std::string >::const_iterator energy_iter;
			for ( energy_iter = energy_names_.begin(); energy_iter != energy_names_.end(); ++energy_iter ) {
				line_stream >> tag;
				if ( *energy_iter != "description" ) { // currently the only text-based field, might change in future.
					Real score_val = (Real) float_of( tag );
					add_energy( *energy_iter, score_val );
				} else {
					line_stream >> tag;
				}
			} // for ( energy_iter ... )
			decoy_tag( tag ); // decoy_tag should be last column of this line.
		} else { // conformation lines
			// parse fold_tree and jump lines
			if ( iter->substr(0,10) == "FOLD_TREE " ) {
				kinematics::FoldTree f;
				line_stream >> f;
				set_fold_tree( f ); // add fold-tree to this SilentStruct
				tr.Debug << "read fold-tree " << f; //"\n" is in fold-tree output
				tr.Debug << "reading " << f.num_jump() << " jumps " << std::endl;
				continue;
			} else if ( iter->substr(0,2) == "RT" ) {
				kinematics::Jump jump;
				line_stream >> jump;
				tr.Debug << "read jump " << jump << std::endl;
				add_jump( jump );
				continue;
			} else if ( iter->substr(0,9) == "SEQUENCE:" ) {
				//tr.Warning << "skipping duplicate sequence declaration " << std::endl;
				continue;
			} else if ( iter->substr(0,6) == "REMARK" ) {
				//tr.Warning << "skipping duplicate sequence declaration " << std::endl;
				continue;
			} else if ( iter->substr(0,7) == "RES_NUM" ) {
				figure_out_residue_numbers_from_line( line_stream );
				continue;
			}


			// parse ss,torsions, and c-alpha coords
			int seqpos;
			Real x, y, z, torsion_value;
			char ss;
			utility::vector1< Real > temp_mainchain_torsions, temp_chi_torsions;

			line_stream >> tag;
			if ( !is_int( tag ) ) {
				tr.Error << "ERROR:  !is_int( " << tag << " ) from line (" << *iter << ")\n";
			}
			runtime_assert( is_int( tag ) ); // this tag should represent the sequence position within the silent-file
			seqpos = int_of( tag );
			line_stream >> ss;

			// It would be nice to not hard-wire these values --
			//  since we'll eventually need to use RNA (or protein) .params files
			//  to create the pose, couldn't we look up the residue
			//  and figure out how many torsions are required?
			for ( Size n = 1; n <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; n++ ){
				line_stream	>> torsion_value;
				temp_mainchain_torsions.push_back( torsion_value  );
			}
			if (fullatom_) {
				for ( Size n = 1; n <= chemical::rna::NUM_RNA_CHI_TORSIONS; n++ ){
					line_stream	>> torsion_value;
					temp_chi_torsions.push_back( torsion_value  );
				}
			}

			//Added Feb. 2009... information on C2', C1', and O4' -- in general they have
			// varying bond lengths and bond angles to keep sugar ring closed.
			line_stream >> x >> y >> z;
			Vector temp_vec( x, y, z );

			line_stream >> tag;

			if ( is_float( tag )  /* New silent format with extra info on sugar atoms*/){
				utility::vector1< Vector> vecs ;
				vecs.push_back( temp_vec );

				x = float_of( tag );
				line_stream >> y >> z;
				temp_vec = Vector( x, y, z );
				vecs.push_back( temp_vec );

				line_stream >> x >> y >> z;
				temp_vec = Vector( x, y, z );
				vecs.push_back( temp_vec );

				set_non_main_chain_sugar_coords( seqpos, vecs );

				line_stream >> x >> y >> z;
				temp_vec = Vector( x, y, z );

				line_stream >> tag;
			}

			set_secstruct( seqpos, ss       );
			set_coords   ( seqpos, temp_vec );
			set_mainchain_torsions( seqpos, temp_mainchain_torsions );
			set_chi_torsions     ( seqpos, temp_chi_torsions );


			if ( tag != decoy_tag() ) { // decoy_tag should be last tag.
				tr.Warning 	<< "parse error(" << *iter << ") " << tag << " != " << decoy_tag() << std::endl;
				success = false;
				break;
			}
		} // conformation lines
	} // for ( iter ... )
	// if no fold-tree available generate a standard tree
	if ( fold_tree().size() < 1 ) {
		fold_tree_.simple_tree( nres() );
		tr.Debug << " generating simple fold-tree " << fold_tree();
	}

	success = true;
	return success;
} // init_from_lines

/// @brief Resize this silent-struct to the appropriate number of residues.
void
RNA_SilentStruct::resize(
	Size const nres_in
) {
	nres( nres_in );
	secstruct_.resize( nres() );
	coords_   .resize( nres() );
	mainchain_torsions_.resize( nres() );
	chi_torsions_     .resize( nres() );
	non_main_chain_sugar_coords_.resize( nres() );
	fold_tree_.simple_tree( nres() );
}

// @brief Fill a Pose with the data in this RNA_SilentStruct.
void RNA_SilentStruct::fill_pose(
	core::pose::Pose & pose
) const {
	using namespace core::chemical;
	ResidueTypeSetCAP residue_set;
	if ( fullatom() ) {
		residue_set = ChemicalManager::get_instance()->residue_type_set( FA_RNA );
	}
	fill_pose( pose, *residue_set );
} // fill_pose

void RNA_SilentStruct::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & /*residue_set*/
) const {
	using namespace core::chemical;

	bool const use_input_pose( false ); // tex hack for refactoring!
	if (use_input_pose)	{
		tr.Info << "Using bond lengths and angles from an input pose." << std::endl;
	}	else {
		tr.Info << "Using ideal geometry from params files..." << std::endl;
		//RHIJU HACK!
		//tr.Info << "USING RNA PARAMS FILES " << std::endl;
		static const ResidueTypeSetCAP rna_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_RNA );
		core::pose::make_pose_from_sequence( pose, sequence(), *rna_residue_set );
	}
	tr.Debug << "FOLD TREE: " << fold_tree();


	// set fold_tree
	pose.fold_tree( fold_tree() );

	// set jumps
	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
		pose.set_jump( nr, jump( nr ) );
	}

	assert( nres() == sequence().length() );

	for ( Size seqpos = 1; seqpos <= nres(); ++seqpos ) {

		// It would be nice to not hard-wire these values --
		//  since we'll eventually need to use RNA (or protein) .params files
		//  to create the pose, couldn't we look up the residue
		//  and figure out how many torsions are required?
		for ( Size n = 1; n <= core::chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; n++ ){
			id::TorsionID rna_torsion_id( seqpos, id::BB, n );
			//			std::cout << rna_torsion_id << " " << mainchain_torsions_[ seqpos ][n ] << std::endl;
			pose.set_torsion( rna_torsion_id,
												mainchain_torsions_[seqpos][n] );
		}


		if (fullatom_) {
			for ( Size n = 1; n <= core::chemical::rna::NUM_RNA_CHI_TORSIONS; n++ ){
				id::TorsionID rna_torsion_id( seqpos, id::CHI, n );
				pose.set_torsion( rna_torsion_id,
													chi_torsions_[seqpos][n] );
			}
		}

		pose.set_secstruct( seqpos, secstruct_[seqpos] );
	}


	if ( non_main_chain_sugar_coords_defined_ ) {
		//Force one refold.
		pose.residue(1).xyz( 1 );

		pose::Pose const & reference_pose = pose; /*try to avoid refolds*/
		for ( Size seqpos = 1; seqpos <= nres(); ++seqpos ) {
			pose::rna::apply_non_main_chain_sugar_coords( non_main_chain_sugar_coords_[ seqpos ], pose, reference_pose, seqpos );
		}
	}

	finish_pose( pose );
} // fill_pose


void
RNA_SilentStruct::print_header( std::ostream& out ) const
{
	SilentStruct::print_header( out );
	out << "REMARK RNA \n";
}

void RNA_SilentStruct::print_conformation( std::ostream & output ) const {

	if ( fold_tree().size() > 1 ) { //assume non-trivial fold_tree only if more than one edge, i.e., EDGE 1 <nres> -1
		output << fold_tree();
	}
	for ( Size i = 1; i <= fold_tree().num_jump(); i++ ) {
		output << jump( i ) << "\n";
	}

	tr.Debug << "FOLD_TREE Size: " << fold_tree().size() << " " << fold_tree()
		<< std::endl;
	for ( Size i = 1; i <= nres(); ++i ) {
		// make sure secstruct is valid
		char this_secstr = secstruct_[i];
		if (this_secstr < 'A' || this_secstr > 'Z')
			this_secstr = 'L';

		output << I( 4, i ) << ' '
					 << this_secstr << ' ';

		for ( Size n = 1; n <= core::chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; n++ ){
			output << F( 9, 3, mainchain_torsions_[i][n] );
		}
		if ( fullatom_ ) {
			for ( Size n = 1; n <= core::chemical::rna::NUM_RNA_CHI_TORSIONS; n++ ){
				output << F( 9, 3, chi_torsions_[i][n] );
			}
		}

		//New, Feb. 2009
		if ( non_main_chain_sugar_coords_defined_ ) {
			for (Size n = 1; n <= non_main_chain_sugar_coords_[i].size(); n++ ) {
				output << F( 12, 6, non_main_chain_sugar_coords_[i][n].x() )
							 << F( 12, 6, non_main_chain_sugar_coords_[i][n].y() )
							 << F( 12, 6, non_main_chain_sugar_coords_[i][n].z() );
			}
		}

		output << F( 9, 3, coords_[i].x() )
					 << F( 9, 3, coords_[i].y() )
					 << F( 9, 3, coords_[i].z() );

		output << ' ' << decoy_tag();
		output << "\n";
	} // for ( Size i = 1; i <= nres; ++i )
} // print_conformation

Real RNA_SilentStruct::get_debug_rmsd() {
	pose::Pose temp_pose;
	FArray2D< Real > rebuilt_coords (3, coords_.size() ), original_coords( 3, coords_.size() );
	static std::string atom_name = "C4'";

	// build temp_pose from coordinates
	fill_pose( temp_pose );

	for ( Size i = 1; i <= temp_pose.total_residue(); ++i ) {
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			rebuilt_coords (k,i) = temp_pose.residue(i).xyz( atom_name )[k-1];
			original_coords(k,i) = coords_[i][k-1];
		}
	}

	Real rmsd = numeric::model_quality::rms_wrapper( temp_pose.total_residue(), rebuilt_coords, original_coords );
	return rmsd;
}

ObjexxFCL::FArray2D< Real >
RNA_SilentStruct::get_CA_xyz() const {
	core::Size n_residues = nres();
	FArray2D< Real > my_coords( 3, n_residues );
	for ( Size i = 1; i <= n_residues; ++i ) { // i = n_residues
		for ( Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
			my_coords(k,i) = coords_[i][k-1];
		} // k
	} // i

	return my_coords;
} // get_CA_positions

Real RNA_SilentStruct::CA_rmsd( RNA_SilentStruct other_pss ) {
	FArray2D< Real > my_coords    = get_CA_xyz();
	FArray2D< Real > other_coords = other_pss.get_CA_xyz();
	Real rmsd = numeric::model_quality::rms_wrapper( nres(), my_coords, other_coords );

	return rmsd;
} // RNA_SilentStruct::CA_rmsd

RNA_SilentStruct & RNA_SilentStruct::operator= (
	RNA_SilentStruct const &
)
{
	utility_exit_with_message( "called ProteinSilentStruct::operator=)" );
	exit(0); //  just to keep the compiler happy
}

} // namespace silent
} // namespace io
} // namespace core
