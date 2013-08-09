// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA de novo fragment assembly Structure Parameters
/// @brief User input parameters for RNA structure inference
/// @detailed
/// @author Rhiju Das


// Unit Headers
#include <protocols/rna/RNA_StructureParameters.hh>
#include <protocols/rna/RNA_JumpLibrary.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/rna/RNA_SecStructInfo.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh> // for  Correctly_position_cutpoint_phosphate_torsions
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_DataInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/ObjexxFCL.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <numeric/random/random.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <list>

#include <core/chemical/VariantType.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/tree/Atom.hh>
//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


namespace protocols{
namespace rna{

static numeric::random::RandomGenerator RG(144620);  // <- Magic number, do not change it!

static basic::Tracer tr( "protocols.rna.rna_structure_parameters" ) ;

using namespace core;

RNA_StructureParameters::RNA_StructureParameters():
	secstruct_defined_( false ),
	assume_non_stem_is_loop( false ),
	root_at_first_rigid_body_( false )
	{
	}

RNA_StructureParameters::~RNA_StructureParameters() {}

//////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::initialize(
	core::pose::Pose & pose,
	std::string const rna_params_file,
	std::string const jump_library_file,
	bool const ignore_secstruct )
{


	if ( rna_params_file.length() > 0 ) read_parameters_from_file( rna_params_file );

	if ( rna_pairing_list_.size() > 0 || chain_connections_.size() > 0 ) {
		rna_jump_library_ = RNA_JumpLibraryOP( new RNA_JumpLibrary( jump_library_file) );
	}

	initialize_allow_insert( pose );

	initialize_secstruct( pose );

	if  ( ignore_secstruct ) override_secstruct( pose );

	if ( virtual_anchor_attachment_points_.size() > 0 ) append_virtual_anchor( pose );

	setup_virtual_phosphate_variants( pose );

}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::override_secstruct( core::pose::Pose & pose ){
	rna_secstruct_ = std::string( pose.total_residue(), 'X' );
	std::cout << "OVER-RIDING SECONDARY STRUCTURE WITH:   " << rna_secstruct_ << std::endl;
	set_rna_secstruct( pose, rna_secstruct_ );
}
/////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::append_virtual_anchor( pose::Pose & pose )
{

	if ( virtual_anchor_attachment_points_.size() == 0 ) return;

	std::cout << "Current last residue is type: " << pose.residue( pose.total_residue() ).name3()  << std::endl;
	std::cout << pose.annotated_sequence() << std::endl;
	if ( pose.residue( pose.total_residue() ).name3() == "XXX" ) return; //already did virtual residue attachment.

	// Fix up the pose.
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type(1).residue_type_set();

	//	std::cout << " CHECK XXX " << residue_set.name3_map("XXX").size() << std::endl;
	//	std::cout << " CHECK YYY " << residue_set.name3_map("YYY").size() << std::endl;
	//	std::cout << " CHECK VRT " << residue_set.name3_map("VRT").size() << std::endl;

	core::chemical::ResidueTypeCOPs const & rsd_type_list( residue_set.name3_map("XXX") );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type_list[1] ) );
	pose.append_residue_by_jump( *new_res, virtual_anchor_attachment_points_[1] );

	Size const virt_res = pose.total_residue();

	// allow insert.
	allow_insert_->append_residue( pose, virt_res, false );

	// Info on pairings and cutpoints.
	cutpoints_open_.push_back( (virt_res - 1) );

	for ( Size n = 1; n <= virtual_anchor_attachment_points_.size(); n++ ){

		RNA_Pairing p;
		p.pos1 = virtual_anchor_attachment_points_[ n ];
		p.pos2 = virt_res;
		p.edge1 = 'X';
		p.edge2 = 'X';
		p.orientation = 'X';
		rna_pairing_list_.push_back( p );

		utility::vector1< Size > obligate_pair;
		obligate_pair.push_back( rna_pairing_list_.size() );
		obligate_pairing_sets_.push_back( obligate_pair );
	}


}


/////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::initialize_secstruct( core::pose::Pose & pose  )
{

	if ( !secstruct_defined_ ) {

		rna_secstruct_ = std::string( pose.total_residue(), 'X' );

		if ( rna_pairing_list_.size() > 0 && assume_non_stem_is_loop ){
			rna_secstruct_ = std::string( pose.total_residue(), 'L' );
		}

		for ( Size n = 1; n <= rna_pairing_list_.size(); n++ ){
			RNA_Pairing const & rna_pairing( rna_pairing_list_[ n ] );
			if (  rna_pairing.edge1 == 'W' &&
						rna_pairing.edge2 == 'W' &&
						rna_pairing.orientation == 'A' &&
						core::chemical::rna::possibly_canonical( pose.residue( rna_pairing.pos1 ).aa(),
																										pose.residue( rna_pairing.pos2 ).aa() ) )	 {
				rna_secstruct_[ rna_pairing.pos1 - 1 ] = 'H';
				rna_secstruct_[ rna_pairing.pos2 - 1 ] = 'H';
			}
		}
	}

	std::cout << "Setting desired secondary structure to: " << rna_secstruct_ << std::endl;

	set_rna_secstruct( pose, rna_secstruct_ );
}

//////////////////////////////////////////////////////////////////////////////
std::list< Size >
RNA_StructureParameters::get_stem_residues( core::pose::Pose const & pose ) const
{
	std::list< Size > in_stem;

	for ( Size n = 1; n <= rna_pairing_list_.size(); n++ ){
		RNA_Pairing const & rna_pairing( rna_pairing_list_[ n ] );
		if (  rna_pairing.edge1 == 'W' &&
					rna_pairing.edge2 == 'W' &&
					rna_pairing.orientation == 'A' &&
					core::chemical::rna::possibly_canonical( pose.residue( rna_pairing.pos1 ).aa(),
																									pose.residue( rna_pairing.pos2 ).aa() ) )	 {
			in_stem.push_back( rna_pairing.pos1 );
			in_stem.push_back( rna_pairing.pos2 );
		}
	}

	in_stem.sort();
	in_stem.unique();

	return in_stem;
}

 /////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::initialize_allow_insert( core::pose::Pose & pose  )
{

 	allow_insert_ = new toolbox::AllowInsert( pose );

 	if (allow_insert_res_.size() > 0 ) {
 		allow_insert_->set( false );
 		for (Size n = 1; n <= allow_insert_res_.size(); n++ ) {

			Size const i = allow_insert_res_[ n ];
			allow_insert_->set( i, true );
			// new -- make sure loops are moveable (& closeable!) at the 3'-endpoint
			if ( i < pose.total_residue() ) allow_insert_->set_phosphate( i+1, pose, true );

 		}
 	} else {
 		allow_insert_->set( true );
 	}

	//std::cout << "ALLOW_INSERT after INITIALIZE_ALLOW_INSERT " << std::endl;
	//allow_insert_->show();

	// We don't trust phosphates at the beginning of chains!
	// Wait, this it total nonsense... captured by virtual phosphate stuff later on!
	// 	allow_insert_->set_phosphate( 1, pose, false );
	//	for ( Size i = 1; i <= cutpoints_open_.size(); i++ ) {
	//		allow_insert_->set_phosphate( cutpoints_open_[i]+1, pose, false );
	//	}

	// std::cout << "ALLOW_INSERT! ALLOW_INSERT! ALLOW_INSERT!" << std::endl;
	// for (Size i = 1; i <= pose.total_residue(); i++ ){
	// 	std::cout << allow_insert_->get( i );
	// }
	// std::cout << std::endl;
	// std::cout << "ALLOW_INSERT! ALLOW_INSERT! ALLOW_INSERT!"<< std::endl;


}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::get_pairings_from_line(
	std::istringstream & line_stream,
	bool const obligate_pairing_set )
{

	Size a,b;
	char e1,e2,o;
	std::string tag;

	utility::vector1< Size > line_pairings;

	while ( !line_stream.fail() ) {

		line_stream >> a >> b ;

		if ( line_stream.fail() ) {
			std::cout << "Parse error!!!" << a << ' ' << b << std::endl;
		}

		line_stream >> tag;
		if ( line_stream.fail() || tag == "PAIR" ) {
			e1 = 'W';
			e2 = 'W';
			o = 'A';
		} else {
			e1 = tag[0];
			line_stream >> e2 >> o;
			if ( line_stream.fail() )  utility_exit_with_message(  "Problem with PAIR readin: "+tag );
		}

		RNA_Pairing p;

		if ( a < b ) {
			p.pos1 = a;
			p.pos2 = b;
		} else {
			p.pos1 = b;
			p.pos2 = a;
		}

		if (a == b) {
			utility_exit_with_message(   "Can't base pair a residue with itself: "+I(3,a)+" "+I(3,b)  );
		}

		p.edge1 = e1;
		p.edge2 = e2;
		p.orientation = o;

		rna_pairing_list_.push_back( p );

		line_pairings.push_back( rna_pairing_list_.size() );

		line_stream >> tag;
		if ( !line_stream.fail() && tag != "PAIR" )  utility_exit_with_message(  "Problem with PAIR readin: " + tag );
	}

	if ( obligate_pairing_set ) {
		obligate_pairing_sets_.push_back( line_pairings );
	} else {
		possible_pairing_sets_.push_back( line_pairings );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::read_chain_connection( std::istringstream & line_stream ) {

	std::string tag;
	Size pos1( 0 ), pos2( 0 ), which_segment( 1 );

	utility::vector1< Size > res_list1, res_list2;

	while ( !line_stream.fail() ) {

		//Find the next non-whitespace character
		char checkchar( ' ' );
		while ( checkchar == ' ' && !line_stream.fail()) checkchar = line_stream.get();

		if (line_stream.fail()  ) break;
		line_stream.putback( checkchar );

		if ( checkchar == 'S' ) {
			line_stream >> tag;

			if ( tag == "SEGMENT1" ) {
				which_segment = 1;
			} else if ( tag == "SEGMENT2" ) {
				which_segment = 2;
			} else {
				utility_exit_with_message(  "Looking for SEGMENT1 or SEGMENT2 in CHAIN_CONNECTION line, but got " + tag );
			}

		} else { //Better be two numbers.

			if (!line_stream.fail()  ){
				line_stream >> pos1 >> pos2;
				runtime_assert( pos2 > pos1 );
				for (Size i = pos1; i <= pos2; i++ ){
					if (which_segment == 1 ) res_list1.push_back( i );
					if (which_segment == 2 ) res_list2.push_back( i );
				}
			}
		}
	}

	if ( res_list1.size() == 0 || res_list2.size() == 0 ){
		utility_exit_with_message(  "Did not specify SEGMENT1 or SEGMENT2 in CHAIN_CONNECTION line?" );
	}

	chain_connections_.push_back( std::make_pair(res_list1, res_list2) );

}


/////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::read_parameters_from_file( std::string const & filename ) {

	tr << "Reading RNA parameters file: " << filename << std::endl;

	Size a;

	utility::io::izstream pairing_stream( filename );
	if ( !pairing_stream ) {
		pairing_stream.close();
		utility_exit_with_message(  "Pairing file? " + filename );
	}

	std::string line, tag;
	//	Size num_stems( 0 );

	while ( getline( pairing_stream, line ) ) {

		std::istringstream line_stream( line );
		line_stream >> tag;

		if (line_stream.fail() ) continue; //Probably a blank line.

		if (tag == "OBLIGATE" ) {

			line_stream >> tag;
			if ( line_stream.fail() || tag != "PAIR" )  utility_exit_with_message(  "Problem with STEM readin: " + tag );
			get_pairings_from_line( line_stream, true /*obligate jump*/ );

		} else if (tag == "STEM"  || tag == "POSSIBLE" ) {

			line_stream >> tag;
			if ( line_stream.fail() || tag != "PAIR" )  utility_exit_with_message(  "Problem with STEM readin: " + tag );
			get_pairings_from_line( line_stream, false /*obligate jump*/ );

		} else if (tag == "ALLOW_INSERT" ) {
			// deprecated!!! switch to ALLOW_INSERT_RES!!
			Size pos1, pos2;
			while ( !line_stream.fail() ) {
				line_stream >> pos1 >> pos2;
				runtime_assert( pos2 >= pos1 );
				for (Size i = pos1; i <= pos2; i++ ) allow_insert_res_.push_back( i );
			}
			//utility_exit_with_message( "No longer reading in ALLOW_INSERT from command line. Try using -s <pdb> instead." );

		} else if (tag == "ALLOW_INSERT_RES" ) {

			Size pos;
			while ( !line_stream.fail() ) {
				line_stream >> pos;
				allow_insert_res_.push_back( pos );
			}
			//utility_exit_with_message( "No longer reading in ALLOW_INSERT from command line. Try using -s <pdb> instead." );

		} else if (tag == "VIRTUAL_ANCHOR" ) {

			while ( !line_stream.fail() ) {
				line_stream >> a;
				if (!line_stream.fail()  ){
					virtual_anchor_attachment_points_.push_back( a );
				}
			}


		} else if (tag == "CUTPOINT_CLOSED" ) {
			while ( !line_stream.fail() ) {
				line_stream >> a;
				if (!line_stream.fail()  ){
					cutpoints_closed_.push_back( a );
				}
			}

		} else if (tag == "CUTPOINT_OPEN" ) {
			while ( !line_stream.fail() ) {
				line_stream >> a;
				if (!line_stream.fail()  ){
					cutpoints_open_.push_back( a );
				}
			}

		} else if (tag == "CHAIN_CONNECTION" ) {

			read_chain_connection( line_stream );

		} else if (tag == "SECSTRUCT" ) {
			line_stream >> rna_secstruct_;
			secstruct_defined_ = true;

		} else if  ( tag[0] == '#' ) {
			continue;

		} else {
			utility_exit_with_message(   "Unrecognized tag in pairing file: " + tag );
		}
	}
	pairing_stream.close();

}


/////////////////////////////////////////////////////////////////////////////////////////////////
bool
in_list( Size const & pos, utility::vector1 < Size > const & res_list )
{
	for (Size n = 1; n <= res_list.size(); n++ ) {
		if ( res_list[n] == pos ) return true;
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_StructureParameters::check_in_chain_connections( Size const & pos1, Size const & pos2 ) const
{
	Size const num_chain_connections( chain_connections_.size() );

	for (Size n = 1; n <= num_chain_connections ; n++ ){
		utility::vector1 < Size > const & res_list1( chain_connections_[n].first );
		utility::vector1 < Size > const & res_list2( chain_connections_[n].second);

		if ( in_list( pos1, res_list1 ) && in_list( pos2, res_list2) ) return n;
		if ( in_list( pos2, res_list1 ) && in_list( pos1, res_list2) ) return n;

	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_StructureParameters::check_forward_backward(
		 pose::Pose & pose,
		 Size const jump_pos ) const
{
	using namespace core::kinematics::tree;

	if ( pose.residue( jump_pos ).is_coarse() || !pose.residue( jump_pos ).is_RNA() )  return true;

	Size atomno = pose.residue( jump_pos ).atom_index( " O5'" );
	AtomCOP current_atom ( & pose.atom_tree().atom( id::AtomID( atomno,jump_pos) ) );
	id::AtomID const parent_id( current_atom->parent()->id() );
	std::string const & parent_name( pose.residue(parent_id.rsd()).atom_name(parent_id.atomno()) );
	//	std::cout << "HELLO ==> " <<  parent_name << std::endl;
	if ( parent_name == " P  " ) {
		return true;
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::add_new_RNA_jump(
  pose::Pose & pose,
	Size const & which_jump,
	bool & success ) const
{
	kinematics::FoldTree fold_tree( pose.fold_tree() ); //Make a copy.

	Size const jump_pos1( fold_tree.upstream_jump_residue( which_jump ) );
	Size const jump_pos2( fold_tree.downstream_jump_residue( which_jump ) );

	// Later can be smart about virtual residue. (rigid body jumps)
	if ( !pose.residue( jump_pos1 ).is_RNA() || !pose.residue( jump_pos2 ).is_RNA() ) return;

	char e1('W') ,e2('W'), o('A');
	bool found_pairing( false );
	//bool flip( false ); // wait, why do we need flip?  // We don't. ~Labonte
	for ( Size n = 1;  n <= rna_pairing_list_.size(); n++ ){
		RNA_Pairing pairing = rna_pairing_list_[n];
		if (pairing.pos1 == jump_pos1 && pairing.pos2 == jump_pos2 ){
			e1 = pairing.edge1;
			e2 = pairing.edge2;
			o  = pairing.orientation;
			found_pairing = true;
			//flip = false;  // set but never used ~Labonte
			break;
		}
		if (pairing.pos1 == jump_pos2 && pairing.pos2 == jump_pos1 ){
			e1 = pairing.edge2;
			e2 = pairing.edge1;
			o  = pairing.orientation;
			found_pairing = true;
			//flip = true;  // set but never used ~Labonte
			break;
		}
	}

	//The jump may have come from a "chain_connection", which
	// doesn't know parallel vs. antiparallel..
	if (!found_pairing) {
		e1 = 'X';
		e2 = 'X';
		o  = 'X';

		// To save time following could be as runtime_assert statement, only
		// active in debug builds.
		if (check_in_chain_connections( jump_pos1, jump_pos2 ) > 0) found_pairing = true;

	}

	if (!found_pairing){
		utility_exit_with_message(  "Trouble finding foldtree pairing in input pairing list? " + I(3,jump_pos1)+' '+I(3,jump_pos2) );
	}

	bool const forward1 = check_forward_backward( pose, jump_pos1 );
	bool const forward2 = check_forward_backward( pose, jump_pos2 );
	//	std::cout << "FORWARD? " << forward1 << " " << forward2 << std::endl;

	std::string atom_name1, atom_name2;

	kinematics::Jump const new_jump = rna_jump_library_->get_random_base_pair_jump(
																	     pose.residue(jump_pos1).name1(),
																			 pose.residue(jump_pos2).name1(),
																			 e1, e2, o,
																			 atom_name1, atom_name2,
																			 success,
																			 forward1, forward2 );

	//	std::cout << "ATOM NAME? " << atom_name1 << " " << atom_name2 << " /  " << jump_pos1 << " " << jump_pos2 << " / " << flip << "  . " << success << std::endl;

	if (!success) return; //shh don't do anything.

	fold_tree.set_jump_atoms( which_jump, atom_name1, atom_name2 );

	pose.fold_tree( fold_tree );

	//	std::cout << " PUTTING IN NEW JUMP: " << which_jump << "    " << atom_name1 << " " << atom_name2 << " " << e1 << " " << e2 << " " << o << "-->" << new_jump << std::endl;
	pose.set_jump( which_jump, new_jump );


}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::sample_alternative_chain_connection( pose::Pose & pose, Size const & which_jump ) const
{
	kinematics::FoldTree fold_tree( pose.fold_tree() ); //Make a copy.

	Size const jump_pos1( fold_tree.upstream_jump_residue( which_jump ) );
	Size const jump_pos2( fold_tree.downstream_jump_residue( which_jump ) );

	Size const chain_connection_number = check_in_chain_connections( jump_pos1, jump_pos2 );

	// Is this jump even part of a chain connection?
	if ( chain_connection_number < 1 ) return;

	utility::vector1 < Size > const & res_list1( chain_connections_[chain_connection_number].first );
	utility::vector1 < Size > const & res_list2( chain_connections_[chain_connection_number].second);

	//Now shift around the jump positions. Hmm, this is potentially dangerous.
	Size ntries( 0 );
	Size const MAX_TRIES( 1000 );
	bool success( false );

	// Get ready for a new fold-tree
	Size const num_jumps( fold_tree.num_jump() );
	ObjexxFCL::FArray2D <int> jump_points( 2, num_jumps );
	ObjexxFCL::FArray1D <int> cuts( num_jumps );
	for (Size n = 1; n <= num_jumps; n++ ){
		jump_points( 1, n ) = fold_tree.jump_point( 1, n );
		jump_points( 2, n ) = fold_tree.jump_point( 2, n );
		cuts( n ) = fold_tree.cutpoint( n );
	}

	// IS this necessary? Or does which_jump == which_jump_in_list?
	Size which_jump_in_list( 0 );
	for ( which_jump_in_list = 1; which_jump_in_list <= num_jumps; which_jump_in_list++ ){
		if ( Size( jump_points( 1, which_jump_in_list ) ) == jump_pos1 &&
				 Size( jump_points( 2, which_jump_in_list ) ) == jump_pos2 ) break;
		if ( Size( jump_points( 1, which_jump_in_list ) ) == jump_pos2 &&
				 Size( jump_points( 2, which_jump_in_list ) ) == jump_pos1 ) break;
	}
	if (which_jump_in_list > num_jumps ) {
		utility_exit_with_message( "Problem with fold tree change --> Jump " + I(3, jump_pos1) + " " + I(3, jump_pos2) );
	}

	while( !success && ntries++ < MAX_TRIES ){
			Size const new1( static_cast<Size>( RG.uniform() * res_list1.size() )  + 1 );
			Size const new2( static_cast<Size>( RG.uniform() * res_list2.size() )  + 1 );
			jump_points( 1, which_jump_in_list ) = res_list1[ new1 ];
			jump_points( 2, which_jump_in_list ) = res_list2[ new2 ];
			success = fold_tree.tree_from_jumps_and_cuts( pose.total_residue(), num_jumps,
																										jump_points, cuts, 1, false /*verbose*/ );
	}

	fill_in_default_jump_atoms( fold_tree, pose );

	//	tr << "Changing fold_tree ==> " << fold_tree << std::endl;

	if (success) pose.fold_tree( fold_tree );


}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::insert_base_pair_jumps( pose::Pose & pose, bool & success ) const
{

	Size const num_jump( pose.num_jump() );

	for (Size i = 1; i <= num_jump; i++ ){

		// Check that we can actually insert here. At least one of the jump partners
		// should allow moves. (I guess the other one can stay fixed).
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( i ) );
		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( i ) );

		if (allow_insert_->get( jump_pos1 ) || allow_insert_->get( jump_pos2 ) ) {
			add_new_RNA_jump( pose, i, success );
		}

		//std::cout << "INSERT RANDOM JUMP: " << i << " " << jump_pos1 << " " << jump_pos2 << std::endl;

		if (!success) return;

	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::setup_fold_tree_and_jumps_and_variants( pose::Pose & pose )
{
	//	std::cout << "ABOUT TO SET UP JUMPS " << std::endl;
	setup_jumps( pose );
	//	std::cout << "ABOUT TO SET UP CHAINBREAKS " << std::endl;
	setup_chainbreak_variants( pose );

	tr << pose.annotated_sequence() << std::endl;
	tr << pose.fold_tree() << std::endl;


}

///////////////////////////////////////////////////////////////
void
RNA_StructureParameters::setup_jumps( pose::Pose & pose )
{

	///////////////////////////////////////////////////////////
	// Basic setup ==> How many jumps? cuts?
	///////////////////////////////////////////////////////////
	Size const nres = pose.total_residue();
	kinematics::FoldTree f( nres );

	Size const num_cuts_closed( cutpoints_closed_.size() );
	Size const num_cuts_open  ( cutpoints_open_.size() );
	Size const num_cuts_total ( num_cuts_closed + num_cuts_open );

	Size const num_pairings( rna_pairing_list_.size() );
	Size const num_obligate_pairing_sets( obligate_pairing_sets_.size() );
	Size const num_possible_pairing_sets( possible_pairing_sets_.size() );
	Size const num_chain_connections( chain_connections_.size() );
	runtime_assert( num_possible_pairing_sets + num_obligate_pairing_sets <= num_pairings );

	////////////////////////////////////////////////////////////////////////
	// Two possibilities for desired fold tree topology:
	//   Jumps dominate, or cuts dominate.
	Size const num_pairings_to_force = std::max( num_obligate_pairing_sets + num_chain_connections,
																							 num_cuts_total );

	////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray2D <int> jump_points( 2, num_pairings_to_force );
	ObjexxFCL::FArray1D <int> cuts( num_pairings_to_force );

	//////////////////////////////////////////////////////////////////////
	// Cuts.
	//////////////////////////////////////////////////////////////////////
	std::vector< int > obligate_cut_points; //switch this to utility::vector1?
	for (Size n = 1; n<= num_cuts_closed; n++ ) 	  obligate_cut_points.push_back( cutpoints_closed_[ n ] );
	for (Size n = 1; n<= num_cuts_open  ; n++ ) 		obligate_cut_points.push_back( cutpoints_open_[n] );
	// for (Size n = 1; n <= num_cuts_total; n++ ) std::cout << "CUT " << obligate_cut_points[ n-1 ]  << std::endl;

	//////////////////////////////////////////////////////////////////////
	// If a cut needs to be randomly chosen, will generally try to
	// place it in a loopy region.
	FArray1D_float cut_bias( nres, 1.0 );
	std::string const & rna_secstruct( get_rna_secstruct( pose ) );
	for ( Size i = 1; i < nres; i++ ) {
		if ( rna_secstruct[i] == 'H' && rna_secstruct[i+1] == 'H' ) {
			cut_bias( i )   = 0.1;
		}
		if ( !pose.residue(i+1).is_RNA() ||
				 !allow_insert_->get( core::id::AtomID( named_atom_id_to_atom_id( core::id::NamedAtomID( " P  ", i+1 ), pose )  ) ) ){
			cut_bias( i ) = 0.0;
		}
		//std::cout << "CUT_BIAS " << i << " " << cut_bias( i ) << std::endl;
	}

	//std::cout << "ALLOW INSERT IN CHECKING CUT_BIAS" << std::endl;
	//	allow_insert_->show();

	//////////////////////////////////////////////////////////////////////
	// Jump residues.
	//////////////////////////////////////////////////////////////////////
	Size ntries( 0 );
	Size const MAX_TRIES( 1000 );
	bool success( false );

	while (!success && ++ntries < MAX_TRIES ){

		// First, put in a jump from each obligate pairing set.
		Size count( 0 );

		for (Size n = 1; n <= num_obligate_pairing_sets ; n++ ){
			Size const pairing_index_in_stem( static_cast<Size>( RG.uniform() * obligate_pairing_sets_[n].size() )  + 1 );
			Size const which_pairing = obligate_pairing_sets_[n][pairing_index_in_stem];
			count++;
			jump_points(1, count) = rna_pairing_list_[which_pairing].pos1;
			jump_points(2, count) = rna_pairing_list_[which_pairing].pos2;
			//			std::cout << "JUMPS1 " <<  jump_points(1,count) << ' ' << jump_points(2,count ) << std::endl;
		}


		// "Chain connections" provide less information about specific residues to pair --
		//  but they're assumed to be obligatory.
		for (Size n = 1; n <= num_chain_connections ; n++ ){
			utility::vector1 < Size > const & res_list1( chain_connections_[n].first );
			utility::vector1 < Size > const & res_list2( chain_connections_[n].second);
			Size const pairing_index_in_list1( static_cast<Size>( RG.uniform() * res_list1.size() )  + 1 );
			Size const pairing_index_in_list2( static_cast<Size>( RG.uniform() * res_list2.size() )  + 1 );
			count++;
			jump_points(1, count) = res_list1[ pairing_index_in_list1 ];
			jump_points(2, count) = res_list2[ pairing_index_in_list2 ];
			//			std::cout << "JUMPS2 " <<  jump_points(1,count) << ' ' << jump_points(2,count ) << std::endl;
		}
		//		std::cout << std::endl;

		// Then, to fill out pairings, look at remaining possible pairing sets (these
		// should typically be Watson-Crick stems, but this setup is general )
		// Note that there might be three stems defined, but we only want two --
		//  following picks a random set of two.
		FArray1D < bool > used_set( num_possible_pairing_sets, false );
		Size num_sets_left( num_possible_pairing_sets );

		while ( count < num_pairings_to_force ) {

			// Find a random pairing among what's remaining.
			Size const set_index( static_cast<Size>( RG.uniform() * num_sets_left )  + 1 );
			Size m( 0 ), set_count( 0 );
			for (m = 1; m <= num_possible_pairing_sets; m++ ){
				if ( !used_set( m ) ) {
					set_count++;
					if (set_count == set_index ) break;
				}
			}

			if (m > num_possible_pairing_sets ) {
				utility_exit_with_message( "Problem with pairing search "+I(3,num_possible_pairing_sets)+" "+I(3,m) );
			}

			Size const pairing_index_in_set( static_cast<Size>( RG.uniform() * possible_pairing_sets_[m].size() )  + 1 );
			Size const which_pairing = possible_pairing_sets_[m][pairing_index_in_set];

			//			std::cout << "USING SET: " << m  << " ==> " << pairing_index_in_set << std::endl;

			count++;
			jump_points(1, count) = rna_pairing_list_[which_pairing].pos1;
			jump_points(2, count) = rna_pairing_list_[which_pairing].pos2;

			used_set( m ) = true;
			num_sets_left--;
		}

		////////////////////////////////////////////////////////////////////////////////
		// Do it, get the fold tree. and set up the pose.
		////////////////////////////////////////////////////////////////////////////////
		//f.tree_from_jumps_and_cuts( nres, num_pairings_to_force, jump_points, cuts, true /* verbose */);

		//		std::cout << "Making attempt " << ntries << std::endl;
		//		for (Size n = 1; n <= num_pairings_to_force; n++ ){
		//			std::cout << "JUMPS " << jump_points(1, n) <<
		//				" " << 	jump_points(2, n)  <<  std::endl;
		//		}

		success = f.random_tree_from_jump_points( nres, num_pairings_to_force, jump_points, obligate_cut_points, cut_bias, 1, true /*enable 1 or NRES jumps*/ );
	}

	if (!success)  utility_exit_with_message( "Couldn't find a freaking tree!" );

	// Hold on to torsion angles in case we need to set up chainbreak residues...
	pose::Pose pose_copy = pose;

	fill_in_default_jump_atoms( f, pose );

	if ( virtual_anchor_attachment_points_.size() > 0 ){

		f.reorder( pose.total_residue() ); //reroot so that virtual residue is fixed.

		if ( root_at_first_rigid_body_ ) {
			utility::vector1< Size > rigid_body_jumps = get_rigid_body_jumps( pose );
			runtime_assert( rigid_body_jumps.size() > 0 );
			Size const anchor_rsd = pose.fold_tree().upstream_jump_residue( rigid_body_jumps[1] );
			std::cout << "ROOTING AT RSD" << anchor_rsd << std::endl;
			f.reorder( anchor_rsd ); //reroot so that partner of virtual residue is fixed
		}

	} else {

		// also useful -- if user has an input pdb, put the root in there, if possible.
		for (Size n = pose.total_residue(); n >= 1; n-- ){ // not sure why I did this backwards...
			if ( pose.residue(n).is_RNA() &&
					 allow_insert_->get_domain( named_atom_id_to_atom_id( id::NamedAtomID( " C1'", n ), pose ) ) == 1 /*1 means the first inputted pose*/ &&
					 possible_root(f,n) ) {
				f.reorder( n );
				break;
			}
		}

	}

	pose.fold_tree( f );

	bool const random_jumps( true ); // For now this is true... perhaps should also have a more deterministic procedure.
	if (random_jumps) insert_base_pair_jumps( pose, success );

	//	pose.dump_pdb( "insert_jumps.pdb" );

	if (!success) {
		utility_exit_with_message( "Trouble inserting base pair jumps into pose -- check residues and edges." );
	}
}



///////////////////////////////////////////////////////////////
void
RNA_StructureParameters::fill_in_default_jump_atoms( kinematics::FoldTree & f, pose::Pose const & pose ) const
{

	// default atoms for jump connections.
	for (Size i = 1; i <= f.num_jump(); i++ ){
		Size const jump_pos1( f.upstream_jump_residue( i ) );
		Size const jump_pos2( f.downstream_jump_residue( i ) );
		f.set_jump_atoms( i,
											chemical::rna::default_jump_atom( pose.residue( jump_pos1 ) ),
											chemical::rna::default_jump_atom( pose.residue( jump_pos2) ) );
	}

}

///////////////////////////////////////////////////////////////
void
RNA_StructureParameters::setup_chainbreak_variants( pose::Pose & pose )
{

	pose::Pose pose_copy = pose;

	// Create cutpoint variants to force chainbreak score computation.
	for ( Size cutpos = 1; cutpos < pose.total_residue(); cutpos++ ){

		if ( ! pose.fold_tree().is_cutpoint( cutpos ) ) continue;

		// Don't assign a chainbreak penalty if user said this was an "open" cutpoint.
		if ( std::find( cutpoints_open_.begin(), cutpoints_open_.end(), cutpos) != cutpoints_open_.end() ) continue;

		//		if ( !pose.residue( cutpos ).is_coarse() ) {//its regular RNA

		std::cout << "Adding chainbreak variants to " << cutpos << std::endl;

		// important! Taken from SWA code.
		protocols::swa::rna::Correctly_position_cutpoint_phosphate_torsions( pose, cutpos, false /*verbose*/ );

		pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpos   );
		pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpos+1 );

		for (Size i = cutpos; i <= cutpos + 1; i++ ){
			for (Size j = 1; j <= pose.residue( i ).mainchain_torsions().size(); j++ ) {
				id::TorsionID torsion_id( i, id::BB, j );
				pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
			} // j
		} // i
		//	} else {
		//		setup_coarse_chainbreak_constraints( pose, cutpos );
		//	}
	}

	allow_insert_->renumber_after_variant_changes( pose );

}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::setup_virtual_phosphate_variants( pose::Pose & pose )
{
	using namespace id;
	using namespace chemical;

	if ( pose.residue( 1 ).is_RNA() ) {
		if ( ! pose.residue_type( 1 ).has_variant_type( VIRTUAL_PHOSPHATE) ) pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, 1  );
		allow_insert_->set_phosphate( 1, pose, false );
	}

	for ( Size i = 1; i <= cutpoints_open_.size(); i++ ){

		Size n = cutpoints_open_[ i ];

		if ( n == pose.total_residue() ){
			utility_exit_with_message( "Do not specify cutpoint_open at last residue of model" );
		}

		if ( pose.residue_type( n   ).has_variant_type( CUTPOINT_LOWER ) ||
				 pose.residue_type( n+1 ).has_variant_type( CUTPOINT_UPPER ) ){
			utility_exit_with_message( "conflicting cutpoint_open & cutpoint_closed" );
		}

		if ( pose.residue_type( n+1 ).is_RNA() ){
					if ( ! pose.residue_type( n+1 ).has_variant_type( VIRTUAL_PHOSPHATE) ) pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, n+1  );
			allow_insert_->set_phosphate( n+1, pose, false );
		}

	}


}

/////////////////////////////////////////////////////////
bool
RNA_StructureParameters::random_jump_change( pose::Pose & pose ) const
{

	using namespace core::conformation;
	using namespace core::id;

	// Any jumps in here to tweak?
	Size const num_jump = pose.num_jump();
	if (num_jump == 0) return false;

	Size ntries( 0 );
	Size const MAX_TRIES( 1000 );

	Size which_jump( 0 );

	while ( ntries++ < MAX_TRIES ){
		// Randomly pick one.
		which_jump = static_cast<Size>( RG.uniform() * num_jump ) + 1 ;

		// Check that we can actually insert here. At least one of the jump partners
		// should allow moves. (I guess the other one can stay fixed).
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( which_jump ) );
		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( which_jump ) ); // Unused variable causes warning.

		Residue const & rsd1 = pose.residue( jump_pos1 );
		AtomID jump_atom_id1( rsd1.atom_index( default_jump_atom( rsd1 ) ), jump_pos1 );

		Residue const & rsd2 = pose.residue( jump_pos2 ); // Unused variable causes warning.
		AtomID jump_atom_id2( rsd2.atom_index( default_jump_atom( rsd2 ) ), jump_pos2 ); // Unused variable causes warning.

		if ( allow_insert_->get( jump_atom_id1 ) ||
				 allow_insert_->get( jump_atom_id2 ) ) break; // found an OK jump.

	}
	if (ntries >= MAX_TRIES ) return false;

	//	std::cout << "will try to change jump " << which_jump << std::endl;

	// Tweak it -- for connections between chains for which the residue is unknown
	//  consider alternative connection -- if there aren't any jumps of these types,
	//  following does nothing.
	sample_alternative_chain_connection( pose, which_jump );

	// Shift jump.
	bool success( false );
	add_new_RNA_jump( pose, which_jump, success );

	return success;

}

/////////////////////////////////////////////////////////
bool
RNA_StructureParameters::check_base_pairs( pose::Pose & pose ) const
{
	using namespace core::pose::rna;
	static scoring::rna::RNA_LowResolutionPotential const rna_low_resolution_potential;

	for (Size n = 1; n <= rna_pairing_list_.size(); n++ ){

		RNA_Pairing const & rna_pairing( rna_pairing_list_[ n ] );
		Size i( rna_pairing.pos1 );
		Size j( rna_pairing.pos2 );
		if (i > j ) {
			i = rna_pairing.pos2;
			j = rna_pairing.pos1;
		}

		// are these part of the pose that is actually being moved?
		if  ( !allow_insert_->get( named_atom_id_to_atom_id( id::NamedAtomID( "C1'", i ), pose ) ) ) continue;
		if  ( !allow_insert_->get( named_atom_id_to_atom_id( id::NamedAtomID( "C1'", j ), pose ) ) ) continue;

		if ( !rna_low_resolution_potential.check_forming_base_pair(pose,i,j) ) {
			tr << "MISSING BASE PAIR " << i << " " << j << std::endl;
			return false;
		}

		if (rna_pairing.edge1 == 'W' && rna_pairing.edge2 == 'W' && rna_pairing.orientation=='A' ) {

			if ( is_cutpoint_open(pose, i) && is_cutpoint_open( pose, j-1) ) {

				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, i, +1 /* sign */)) return false;
				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, j, -1 /* sign*/ )) return false;

			} else if ( is_cutpoint_open( pose, i-1 ) && is_cutpoint_open( pose, j) ){

				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, i, -1 /* sign */)) return false;
				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, j, +1 /* sign*/ )) return false;

			}
		}


	}

	return true;

}

///////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::setup_base_pair_constraints( core::pose::Pose & pose )
{

	using namespace core::id;

	utility::vector1< std::pair< Size, Size > > pairings;

	// Go through all pairings and define atom pair constraints that will bring together
	//  appropriate atoms on bases.
	for (Size n = 1; n <= rna_pairing_list_.size(); n++ ){

		RNA_Pairing const & rna_pairing( rna_pairing_list_[ n ] );
		Size const & i( rna_pairing.pos1 );
		Size const & j( rna_pairing.pos2 );

		//Basic check that its canonical...
		if ( !( rna_pairing.edge1 == 'W' && rna_pairing.edge2 == 'W' && rna_pairing.orientation == 'A' ) ) {
			tr <<  "skipping constraints for non-canonical base pair: " << I(3,i) << " " << I(3,j) << " " << rna_pairing.edge1 << " " << rna_pairing.edge2 << " " << rna_pairing.orientation << std::endl;
			continue;
		}

		if ( !pose.residue(i).is_RNA() ) continue;
		if ( !pose.residue(j).is_RNA() ) continue;

		if ( !allow_insert_->get( named_atom_id_to_atom_id( NamedAtomID( "C1'", i ), pose ) ) &&
				 !allow_insert_->get( named_atom_id_to_atom_id( NamedAtomID( "C1'", j ), pose ) ) ) continue; //assumed to be frozen, so no need to set up constraints?

		pairings.push_back( std::make_pair( i, j ) ) ;
	}

	// In RNA_ProtocolUtil.cc :
	protocols::rna::setup_base_pair_constraints( pose, pairings, suppress_bp_constraint_ );
}

/////////////////////////////////////////////////////////////////////
std::map< Size, Size >
RNA_StructureParameters::connections() const
{

	std::map< Size, Size > connections_local;
	for ( Size n = 1; n <= rna_pairing_list_.size(); n++ ) {
		RNA_Pairing pairing = rna_pairing_list_[n];
		connections_local[ pairing.pos1 ] = pairing.pos2;
		connections_local[ pairing.pos2 ] = pairing.pos1;
	}
	return connections_local;
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::set_allow_insert(toolbox::AllowInsertOP allow_insert_in )
{
	allow_insert_ = allow_insert_in;
}

////////////////////////////////////////////////////////////////////////////////////////
toolbox::AllowInsertOP
RNA_StructureParameters::allow_insert(){ return allow_insert_; }

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::set_jump_library( RNA_JumpLibraryOP rna_jump_library )
{
	rna_jump_library_ = rna_jump_library;
}

}
}



