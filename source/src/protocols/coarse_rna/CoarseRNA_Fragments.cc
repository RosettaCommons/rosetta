// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_Fragments
/// @details
/// @author Rhiju Das


#include <protocols/coarse_rna/CoarseRNA_Fragments.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/toolbox/AtomID_Mapper.hh>
#include <protocols/farna/util.hh>
#include <protocols/farna/secstruct/RNA_SecStructInfo.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/util.hh>
#include <core/pose/copydofs/util.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// External library headers


//C++ headers
#include <vector>
#include <string>
#include <sstream>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>


using namespace core;
using basic::T;

static THREAD_LOCAL basic::Tracer TR( "protocols.coarse_rna.coarse_rna_fragment_mover" );


namespace protocols {
namespace coarse_rna {

SourcePositions::~SourcePositions() {}

CoarseRNA_Fragments::~CoarseRNA_Fragments() {}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This sort of repeats a lot of stuff in protocols/farna/fragments/RNA_Fragments
//
//  Not quite sure whether we should unify, or make subclasses of a Fragments class...
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
CoarseRNA_Fragments::CoarseRNA_Fragments( std::string const & frag_source_file ):
	RNA_Fragments(),
	frag_source_file_( frag_source_file )
{
	initialize_frag_source_pose();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_Fragments::initialize_frag_source_pose(){

	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::kinematics;

	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( "coarse_rna" ) );

	if ( frag_source_file_.substr( frag_source_file_.size()-4, frag_source_file_.size() ) == ".pdb" ) {
		Pose pose;
		import_pose::pose_from_pdb( pose, *rsd_set, frag_source_file_ );
		protocols::farna::figure_out_secstruct( pose );
		frag_source_secstruct_ = protocols::farna::secstruct::get_rna_secstruct( pose );
		frag_source_pose_ = core::pose::MiniPoseOP( new MiniPose( pose ) );
	} else {

		// Figure out correspondence: P,S,CEN,X,Y --> 1,2,3,4,5. Hopefully!
		ResidueTypeCOP const & rsd_type = rsd_set->get_representative_type_aa( na_rad );
		for ( Size i = 1; i <= rsd_type->natoms(); i++ ) {
			coarse_rna_name_to_num_[ rsd_type->atom_name( i ) ] = i;
		}

		std::string line;
		char dummy_char;
		Real x,y,z;
		std::string sequence = "";
		frag_source_secstruct_ = "";
		utility::vector1< utility::vector1< PointPosition > > all_res_coords;
		utility::io::izstream coords_in( frag_source_file_.c_str() );
		while (  getline( coords_in, line) ) {

			std::istringstream line_stream( line );
			line_stream >> dummy_char;
			sequence += dummy_char;

			line_stream >> dummy_char;
			frag_source_secstruct_ += dummy_char;

			utility::vector1< PointPosition > res_coords;
			while ( line_stream.good() ) {
				line_stream >> x >> y >> z;
				res_coords.push_back( Vector( x,y,z) );
			}
			if ( res_coords.size() != coarse_rna_name_to_num_.size() ) utility_exit_with_message( "Should only have 5? atoms in mini ppose file.");
			all_res_coords.push_back( res_coords );
		}

		//  std::cout << "SIZE     ! " << all_res_coords.size() << std::endl;
		//  std::cout << "SEQUENCE ! " << sequence << std::endl;
		//  std::cout << "SECSTRUCT! " << frag_source_secstruct_ << std::endl;
		frag_source_pose_ = core::pose::MiniPoseOP( new MiniPose( all_res_coords, FoldTree( all_res_coords.size() ), sequence ) );

	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_Fragments::insert_fragment(
	core::pose::Pose & pose,
	Size const & insert_res,
	Size const & source_res,
	Size const & frag_size,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) const
{

	using namespace core::id;

	std::map< Size, Size > res_map;

	for ( Size offset = 0; offset < frag_size; offset ++ ) {

		Size const insert_offset = insert_res + offset;
		Size const source_offset = source_res + offset;

		if ( (insert_offset) > pose.total_residue() ) continue;
		if ( (source_offset) > frag_source_pose_->total_residue() ) continue;

		res_map[ insert_offset ] = source_offset;

	}
	//  std::cout << std::endl;

	std::map< AtomID, AtomID > atom_id_map = atom_level_domain_map->atom_id_mapper()->calculate_atom_id_map( pose, res_map, frag_source_pose_->fold_tree() );
	core::pose::copydofs::copy_dofs(  pose, *frag_source_pose_, atom_id_map, atom_level_domain_map->calculate_atom_id_domain_map( pose ) );

}

///////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_Fragments::find_source_positions( SequenceSecStructPair const & key ) const
{

	using namespace protocols::farna;

	SourcePositionsOP source_positions( new SourcePositions );

	std::string const RNA_string = key.first;
	std::string const RNA_secstruct_string = key.second;

	Size const size = RNA_string.length();
	static Distance const DIST_CUTOFF( 10.0 ); /*distance of S to next P*/

	runtime_assert( RNA_string.length() == RNA_secstruct_string.length() );

	// dummy initialization.
	std::string vall_current_sequence ( RNA_string );
	std::string vall_current_secstruct( RNA_secstruct_string );

	std::string const & source_secstruct( frag_source_secstruct_ );
	std::string const & source_sequence( frag_source_pose_->sequence() );

	for ( Size i = 1; i <= source_sequence.size() - size + 1; i++ ) {

		bool match( true );

		for ( Size offset = 0; offset < size; offset++ ) {
			vall_current_sequence [offset] = source_sequence[ i - 1 + offset ];
			vall_current_secstruct[offset] = source_secstruct[ i - 1 + offset ];

			if (
					!compare_RNA_char( vall_current_sequence[offset], RNA_string[ offset ] ) ||
					!compare_RNA_secstruct( vall_current_secstruct[offset], RNA_secstruct_string[ offset ] ) ) {
				match = false;
				break;
			}
			// check for chainbreak
			if ( offset > 1 &&
					( ( frag_source_pose_->xyz( id::AtomID( 2, i+offset-1 ) ) -
					frag_source_pose_->xyz( id::AtomID( 1, i+offset   ) ) ).length() > DIST_CUTOFF ) ) continue;

		}

		if ( match ) {
			source_positions->push_back( i );
		}

	}

	std::cout << "Picked Fragment Library for sequence " << RNA_string << " " <<
		" and sec. struct " << RNA_secstruct_string << " ... found " <<
		source_positions->size() << " potential fragments" << std::endl;

	source_positions_map_[ key ] = source_positions;

}

///////////////////////////////////////////////////////////////////////////////////////
Size
CoarseRNA_Fragments::pick_random_fragment(
	std::string const RNA_string,
	std::string const RNA_secstruct_string,
	Size const type ) const
{

	std::string const RNA_string_local = protocols::farna::convert_based_on_match_type( RNA_string, type );

	SequenceSecStructPair const key( std::make_pair( RNA_string_local, RNA_secstruct_string ) );

	if ( ! source_positions_map_.count( key ) ) {
		find_source_positions( key );
	}

	SourcePositionsOP source_positions = source_positions_map_[ key ];

	Size const num_frags = source_positions->size();

	if ( num_frags == 0 ) { //trouble.
		std::cout << "Fragment Library: zero fragments found for " << RNA_string_local << std::endl;
		std::cerr << "Fragment Library: zero fragments found for " << RNA_string_local << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	Size const which_frag = static_cast <Size> ( numeric::random::rg().uniform() * num_frags) + 1;

	return (*source_positions)[ which_frag ];

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
Size
CoarseRNA_Fragments::pick_random_fragment(
	core::pose::Pose & pose,
	Size const position,
	Size const size,
	Size const type ) const
{

	std::string const & RNA_sequence( pose.sequence() );
	std::string const & RNA_string = RNA_sequence.substr( position - 1, size );

	std::string const & RNA_secstruct( protocols::farna::secstruct::get_rna_secstruct( pose ) );
	std::string const & RNA_secstruct_string = RNA_secstruct.substr( position - 1, size );

	return pick_random_fragment( RNA_string, RNA_secstruct_string, type );

}

////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_Fragments::apply_random_fragment(
	core::pose::Pose & pose,
	core::Size const position,
	core::Size const size,
	core::Size const type,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) const
{
	Size const source_res = pick_random_fragment( pose, position, size, type );
	//  std::cout << "applying to fragment position " << position << " from source position " << source_res << std::endl;
	insert_fragment( pose, position, source_res, size, atom_level_domain_map );
}

////////////////////////////////////////////////////////////////////////////////////////
bool
CoarseRNA_Fragments::is_fullatom(){ return false; }


} // namespace rna
} // namespace protocols
