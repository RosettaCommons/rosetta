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
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/rna/RNA_SecStructInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_DataInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>

#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/id/AtomID.hh>

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

//Auto Headers
#include <core/chemical/VariantType.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


namespace protocols{
namespace rna{

static numeric::random::RandomGenerator RG(144620);  // <- Magic number, do not change it!

static basic::Tracer tr( "protocols.rna.rna_structure_parameters" ) ;

using namespace core;

RNA_StructureParameters::RNA_StructureParameters(){
	secstruct_defined_ = false;
	assume_non_stem_is_loop = false;
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

	initialize_secstruct( pose );
	if  ( ignore_secstruct ) override_secstruct( pose );

	initialize_allow_insert( pose );

	if ( rna_pairing_list_.size() > 0 || chain_connections_.size() > 0 ) {
		rna_jump_library_ = RNA_JumpLibraryOP( new RNA_JumpLibrary( jump_library_file) );
	}

	//setup_base_pair_constraints( pose );

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
						core::scoring::rna::possibly_canonical( pose.residue( rna_pairing.pos1 ).aa(),
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
					core::scoring::rna::possibly_canonical( pose.residue( rna_pairing.pos1 ).aa(),
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

	allow_insert_.dimension( pose.total_residue() );

	if (allow_insert_segments_.size() > 0 ) {
		allow_insert_ = false;
		for (Size n = 1; n <= allow_insert_segments_.size(); n++ ) {
			for (Size i = allow_insert_segments_[n].first;
					 i <= allow_insert_segments_[n].second;
					 i++ ){
				allow_insert_( i ) = true;
			}
		}
	} else {
		allow_insert_ = true;
	}

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
			Size pos1, pos2;
			while ( !line_stream.fail() ) {
				line_stream >> pos1 >> pos2;
				allow_insert_segments_.push_back( std::make_pair( pos1, pos2 ) );
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

	Size atomno = pose.residue( jump_pos ).atom_index( " O5*" );
	Atom const * current_atom ( & pose.atom_tree().atom( id::AtomID( atomno,jump_pos) ) );
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


	char e1('W') ,e2('W'), o('A');
	bool found_pairing( false );
	bool flip( false ); // wait, why do we need flip?
	for ( Size n = 1;  n <= rna_pairing_list_.size(); n++ ){
		RNA_Pairing pairing = rna_pairing_list_[n];
		if (pairing.pos1 == jump_pos1 && pairing.pos2 == jump_pos2 ){
			e1 = pairing.edge1;
			e2 = pairing.edge2;
			o  = pairing.orientation;
			found_pairing = true;
			flip = false;
			break;
		}
		if (pairing.pos1 == jump_pos2 && pairing.pos2 == jump_pos1 ){
			e1 = pairing.edge2;
			e2 = pairing.edge1;
			o  = pairing.orientation;
			found_pairing = true;
			flip = true;
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

	if (!success) return; //sshh... don't do anything.

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
		if (allow_insert_( jump_pos1 ) || allow_insert_( jump_pos2 ) ) {
			add_new_RNA_jump( pose, i, success );
		}

		//std::cout << "INSERT RANDOM JUMP: " << i << " " << jump_pos1 << " " << jump_pos2 << std::endl;

		if (!success) return;

	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StructureParameters::setup_jumps( pose::Pose & pose, bool const random_jumps /* = true */ )
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


	//////////////////////////////////////////////////////////////////////
	// If a cut needs to be randomly chosen, will generally try to
	// place it in a loopy region.
	FArray1D_float cut_bias( nres, 0.1 );
	std::string const & rna_secstruct( get_rna_secstruct( pose ) );
	for ( Size i = 2; i <= nres; i++ ) {

		if ( rna_secstruct[i-1] != 'H' ) {
			cut_bias( i )   = 1.0;
			cut_bias( i-1 ) = 1.0;
		}

		if ( !allow_insert_(i-1) && !allow_insert_(i) ){
			cut_bias( i - 1) = 0.0;
		}

	}


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
		}


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

	pose.fold_tree( f );

	if (random_jumps) insert_base_pair_jumps( pose, success );

	pose.dump_pdb( "insert_jumps.pdb" );

	if (!success) {
		utility_exit_with_message( "Trouble inserting base pair jumps into pose -- check residues and edges." );
	}

	// Create cutpoint variants to force chainbreak score computation.
	for ( Size cutpos = 1; cutpos < nres; cutpos++ ){

		if ( ! f.is_cutpoint( cutpos ) ) continue;

		// Don't assign a chainbreak penalty if user said this was an "open" cutpoint.
		if ( std::find( cutpoints_open_.begin(), cutpoints_open_.end(), cutpos) != cutpoints_open_.end() ) continue;

		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpos   );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpos+1 );

		for (Size i = cutpos; i <= cutpos + 1; i++ ){
			for (Size j = 1; j <= scoring::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
				id::TorsionID torsion_id( i, id::BB, j );
				pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
			} // j
		} // i

	} // n

}


/////////////////////////////////////////////////////////
bool
RNA_StructureParameters::random_jump_change( pose::Pose & pose ) const
{
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
		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( which_jump ) );
		if (allow_insert_( jump_pos1 ) || allow_insert_( jump_pos2 ) ) break;
	}
	if (ntries >= MAX_TRIES ) return false;

	//std::cout << "will try to change jump " << which_jump << std::endl;

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
	using namespace scoring::rna;
	static RNA_LowResolutionPotential const rna_low_resolution_potential;

	for (Size n = 1; n <= rna_pairing_list_.size(); n++ ){

		RNA_Pairing const & rna_pairing( rna_pairing_list_[ n ] );
		Size i( rna_pairing.pos1 );
		Size j( rna_pairing.pos2 );
		if (i > j ) {
			i = rna_pairing.pos2;
			j = rna_pairing.pos1;
		}

		if ( !rna_low_resolution_potential.check_forming_base_pair(pose,i,j) ) {
			std::cout << "MISSING BASE PAIR " << i << " " << j << std::endl;
			return false;
		}

		if (rna_pairing.edge1 == 'W' && rna_pairing.edge2 == 'W' && rna_pairing.orientation=='A' ) {
			if ( is_rna_chainbreak(pose, i) || is_rna_chainbreak( pose, j-1) ) {
				//std::cout << "CHECK1: " << i << " " << j << std::endl;
				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, i, +1 /* sign */)) return false;
				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, j, -1 /* sign*/ )) return false;
			} else if ( is_rna_chainbreak( pose, i-1 ) || is_rna_chainbreak( pose, j) ){
				//std::cout << "CHECK2: " << i << " " << j << std::endl;
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
		pairings.push_back( std::make_pair( i, j ) ) ;
	}

	// In RNA_ProtocolUtil.cc :
	protocols::rna::setup_base_pair_constraints( pose, pairings );
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
RNA_StructureParameters::and_allow_insert( FArray1D< bool > const & allow_insert_in )
{
	assert( allow_insert_in.size() == allow_insert_.size() );
	for (Size n = 1; n <= allow_insert_.size(); n++ ) {
		allow_insert_( n ) = allow_insert_( n ) && allow_insert_in( n );
		//		std::cout << "ALLOW_INSERT " << n << " " << allow_insert_( n ) << std::endl;
	}
}

}
}



