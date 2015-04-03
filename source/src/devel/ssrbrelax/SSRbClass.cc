// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Vatsan Raman

// Unit Headers
#include <devel/ssrbrelax/SSRbClass.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <iostream>
#include <sstream>
#include <map>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/xyzVector.hh>

// option key includes

#include <basic/options/keys/SSrbrelax.OptionKeys.gen.hh>


using namespace core;
using basic::T;
using basic::Error;
using basic::Warning;

namespace devel{
namespace ssrbrelax {

	std::ostream & operator<<( std::ostream & os, const Rb & rb ) {
		os << rb.seg_begin_ << " " << rb.seg_end_ << " " << rb.anchor_pos_ << " "
			 << std::endl;
		return os;
	}

	std::ostream & operator<<( std::ostream & /*os*/, const RbSegments & /*rbsegments*/ ) {
	  std::cerr << "operator<< for RbSegments not implemented in mini " << std::endl;
    exit(-1);
	}


	/////////////////////////////////////////////////////////

	void
	RbSegments::one_segment_fold_tree(
						kinematics::FoldTree & f,
						int const total_residue
						)
	{
		using namespace kinematics;
		f.clear();

		if ( num_rb() > 1 ) {
			std::cout << "This is a single segment fold tree"
		<< "Multiple segments are present in the RbSegments object"
		<< std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}

		Rb one_segment_rb( rb_list.at( 0 ) ); //segment currently passed

		int seg_begin( one_segment_rb.seg_begin() );
		int seg_end( one_segment_rb.seg_end() );
		int anchor_pos( one_segment_rb.anchor_pos() );

		assert( seg_begin < seg_end );
		assert( anchor_pos - 1 != total_residue - anchor_pos );//unlikely
		assert( anchor_pos < seg_begin || anchor_pos > seg_end );

		//Is the anchor position close to n or c terminus.That will determine how the fold tree is set up
		//    bool const close_to_n ( ( anchor_pos - 1 < total_residue - anchor_pos ) ? true : false );
		//now a different definition of close_to_n
		bool const close_to_n( ( anchor_pos < seg_begin ) ? true : false );
		int const flex_jump ( seg_begin + ( seg_end - seg_begin )/2 );
		int const fixed_jump( close_to_n ? seg_end + ( total_residue - seg_end )/2 : ( seg_begin - 1 )/2 );

		std::cout << seg_begin << " " << seg_end << " " << anchor_pos << " " <<  flex_jump << " " << fixed_jump << " " << close_to_n << std::endl;

		FlexJumpRes[ std::make_pair( seg_begin, seg_end )] = flex_jump;


		if ( close_to_n ) {

			f.add_edge( 1, anchor_pos, Edge::PEPTIDE );//EDGES
			f.add_edge( anchor_pos, seg_begin - 1, Edge::PEPTIDE );
			f.add_edge( seg_end + 1, fixed_jump, Edge::PEPTIDE );
			f.add_edge( fixed_jump, total_residue, Edge::PEPTIDE );
			f.add_edge( seg_begin, flex_jump, Edge::PEPTIDE );
			f.add_edge( flex_jump, seg_end, Edge::PEPTIDE );

			//set jump 1 as flexible
			f.add_edge( anchor_pos, flex_jump, 1 ); //JUMPS
			f.add_edge( anchor_pos, fixed_jump, 2 );

		} else { ///THERE IS SOME PROBLEM HERE : CHECK THIS AGAIN ******************


			f.add_edge( 1, fixed_jump, Edge::PEPTIDE ); //EDGES
			f.add_edge( fixed_jump, seg_begin - 1, Edge::PEPTIDE );
			f.add_edge( seg_begin, flex_jump, Edge::PEPTIDE );
			f.add_edge( flex_jump, seg_end, Edge::PEPTIDE );
			f.add_edge( seg_end + 1, anchor_pos, Edge::PEPTIDE );
			f.add_edge( anchor_pos, total_residue, Edge::PEPTIDE );

			//set jump 1 as flexible
			f.add_edge( anchor_pos, flex_jump, 1 ); //JUMPS
			f.add_edge( anchor_pos, fixed_jump, 2 );
		}
		std::cout << "Fold tree " << f << std::endl;
		f.reorder( 1 );
		std::cout << "Fold tree reordereed " << f << std::endl;

	}

	///////////////////////////////////////////////////////////
	void
	RbSegments::read_segments_from_file()
	{

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::string filename( option[ OptionKeys::SSrbrelax::rb_file ]().name() );

		utility::io::izstream data( filename.c_str() );
		if( !data ) {
			std::cout << "Couldn't open rb file " << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		std::string line;
		while( getline( data, line ) ) {
			std::istringstream line_stream( line );
			int index, seg_begin, seg_end, anchor_pos,loop_stub_n, loop_stub_c; //the loop stub is the start position for loop modeling
			line_stream >> index >> seg_begin >> seg_end >> anchor_pos >> loop_stub_n >> loop_stub_c;
//There are two loop stubs, one for loop n-term of segment and other c-term of segment
// The loops are stored as a loop object
// Note that the cut point is always the seg_begin or seg_end positions
			if ( !line_stream.fail() ) {
	add_segment( index, seg_begin, seg_end, anchor_pos, loop_stub_n, loop_stub_c );
			}
		}
		data.close();
		data.clear();

	}

	///////////////////////////////////////////////////////////

	void
	RbSegments::read_param_file()
	{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::string param_file( option[ OptionKeys::SSrbrelax::rb_param_file ]().name() );

		utility::io::izstream data( param_file.c_str() );
		if( !data ) {
			std::cout << "Couldn't open rb param file " << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		std::string line;
		while( getline( data, line ) ) {
			std::istringstream line_stream( line );
			int index, dof;
			float stddev;
			line_stream >> index >> dof >> stddev;
			if( !line_stream.fail() ) {
				std::cout << "index dof " << index << " " << dof << " " << stddev << std::endl;
				rbfunc[ std::make_pair( index, dof ) ]  = stddev;
			}
		}

	}


	///////////////////////////////////////////////////////////////////

	void
	RbSegments::add_segment(
				int const index,
				int const seg_begin,
				int const seg_end,
				int const anchor_pos,
				int const loop_stub_n,
				int const loop_stub_c
				)
	{

		bool segment_ok = verify_segment( seg_begin, seg_end, anchor_pos, loop_stub_n, loop_stub_c );

		if ( !segment_ok ) {
			std::cout << "Error adding segment"
		<< std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		rb_list.push_back( Rb( index, seg_begin, seg_end, anchor_pos, loop_stub_n, loop_stub_c ) );
	}

	//////////////////////////////////////////////////////////////////
	void
	RbSegments::add_segment(
				const RbSegments::iterator & it
				)
	{
		add_segment( it->index(), it->seg_begin(), it->seg_end(), it->anchor_pos(), it->loop_stub_n(), it->loop_stub_c() );
	}


	///////////////////////////////////////////////////////////////////////

	void
	RbSegments::delete_segment(
					 int const seg_begin,
					 int const seg_end
					 )
	{

		assert( seg_begin < seg_end );
		bool segment_deleted( false );

		for ( iterator it=rb_list.begin(), it_end=rb_list.end(); it!=it_end; ++it ){
			if ( seg_begin == it->seg_begin() && seg_end == it->seg_end() ) {
	rb_list.erase( it );
	segment_deleted = true;
			}
		}
		if ( !segment_deleted ) {
			std::cout << "Segment not present in rb list "
		<< seg_begin << " " << seg_end
		<< std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}

	////////////////////////////////////////////////////////////////////
	void
	RbSegments::delete_segment(
					 const RbSegments::iterator & it
					 )
	{
		delete_segment( it->seg_begin(), it->seg_end() );
	}


	/////////////////////////////////////////////////////////////////////

	RbSegments::iterator
	RbSegments::one_random_segment()

	{
		int size = num_rb();
		assert( size > 0 );
		int index = 0;
		int const end = int ( numeric::random::uniform()*size );
		iterator it = rb_list.begin();
		while( index != end ) {
			++index;
			++it;
		}
		return it;

	}
	///////////////////////////////////////////////////////////////////////


	bool
	RbSegments::verify_segment(
					 int const seg_begin,
					 int const seg_end,
					 int const anchor_pos,
					 int const loop_stub_n,
					 int const loop_stub_c
					 )
	{

		//need to put some check for jump anchor positions as well.
		if ( seg_begin < seg_end ) {
			for ( const_iterator it = rb_list.begin(), it_end = rb_list.end(); it != it_end; ++it ) {
	if ( seg_end < it->seg_begin() || seg_begin > it->seg_end() ) {
		//do nothing, just makes sure that there is no overlap between multiple segments
	} else {
		std::cout << "Rb segment error\n "
				<< "Overlapping regions\n"
				<< "Existing segment definition " << it->seg_begin() << " " << it->seg_end()
				<< "New segment added " << seg_begin << " " << seg_end
				<< std::endl;
		return false;

	}
			}
			if( seg_begin < anchor_pos && anchor_pos < seg_end ) {
	std::cout << "Anchor position cannot lie between seg_begin and seg_end "
			<< seg_begin << " " << seg_end << " " << anchor_pos
			<< std::endl;
	return false;
			} else if ( ( seg_begin < loop_stub_n && loop_stub_n < seg_end ) || ( seg_begin < loop_stub_c  && loop_stub_c < seg_end ) ) {
	std::cout << "Loop stub cannot lie between seg_begin and seg_end"
			<< seg_begin << " " << seg_end << " " << loop_stub_n << " " << loop_stub_c
			<< std::endl;
	return false;
			} else if ( ( anchor_pos >= loop_stub_n && anchor_pos <= seg_begin ) && ( anchor_pos >= seg_end ) && ( anchor_pos <= loop_stub_c ) ) {
	std::cout << "Anchor position has to be in a fixed regions of the protein "
			<< "Here the anchor position is the loop_region between the loop_stub and the segment tip "
			<< seg_begin << " " << seg_end << " " << anchor_pos << " " << loop_stub_n << " " << loop_stub_c
			<< std::endl;
	return false;
			}

		} else {
			std::cout << "Rb segment definition error\n"
		<< "seg_begin, seg_end, anchor_pos " << seg_begin << " " << seg_end << " " << anchor_pos
		<< std::endl;
			return false;
		}
		return true;
	}
	/////////////////////////////////////////////////////////////////////////////////

	float RbSegments::get_gaussian_parameters(
						 int const & index,
						 int const & dof
						 )

	{
		std::cout << "RbSegments::get_gaussian_parameters " << index << " " << dof << std::endl;
		for ( RbSegments::RbFunc_const_it it=rbfunc.begin(), it_end = rbfunc.end();  it !=it_end; ++it ) {
			std::cout << "get_gaussian index " << it->first.first << " " << it->first.second << std::endl;
			if ( it->first.first == index && it->first.second == dof ) {
				std::cout << "great ! stddev " << it->second << std::endl;
				return it->second;
			}
		}
		return 0.0;
	}

	/////////////////////////////////////////////////////////////////////////////////

	bool RbSegments::distribution_exists(
							 int const & index,
							 int const & dof
							 )
	{
		for ( RbSegments::RbFunc_const_it it=rbfunc.begin(), it_end = rbfunc.end(); it != it_end; ++it ) {
			if ( it->first.first == index && it->first.second == dof ) {
	return true;
			}
		}
		return false;

	}

	/////////////////////////////////////////////////////////////////////////////////
	int RbSegments::get_flex_jump_pos(
																		int const & seg_begin,
																		int const & seg_end
																		)
	{
		for ( std::map< std::pair< int, int >, int >::const_iterator it = FlexJumpRes.begin(), it_end = FlexJumpRes.end(); it != it_end; ++it ) {
			if ( it->first.first == seg_begin && it->first.second == seg_end ) {
				return it->second;
		}
		}
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////
	Vector const RbSegments::X_axis(
																	core::pose::Pose & pose
																	)
	{

		using namespace conformation;

		assert( this->num_rb() == 1 );
		//The X-axis has to be chosen in a secondary structure dependent way
		//If the anchor_pos_res is in a strand, the X-axis is the vector
		// passing through i-1 and i+1 C-alphas
		//If the fixed_jum_pos_res is in a helix, the X-axis is the vector
		//passing through i and i+3 or i-3
		int const anchor_pos( this->anchor_res() );

		if ( pose.secstruct( anchor_pos ) != 'H' && pose.secstruct( anchor_pos ) != 'E' ) {
			Error() <<"The anchor pos has to be in a strand or helix " << anchor_pos << " SS type " << pose.secstruct( anchor_pos ) << "\n";
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}

		if ( pose.secstruct( anchor_pos ) == 'E' ) {
			Residue const & upstream_rsd( pose.residue( anchor_pos - 1 ) );
			Residue const & downstream_rsd( pose.residue( anchor_pos + 1 ) );
			if ( pose.secstruct( anchor_pos - 1 ) != 'E' && pose.secstruct( anchor_pos + 1 ) != 'E' ) {
				Error() << "The anchor_pos is in a strand; the positions up and downstream must be E also " << anchor_pos - 1 << " " << anchor_pos + 1 << "\n";
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}

			return ( upstream_rsd.xyz("CA") - downstream_rsd.xyz("CA") ).normalized();
		}

		//******* X-axis for HELIX undefined yet !!!! ***********
		//the following is just to suppress compiler warning
			Residue const & upstream_rsd( pose.residue( anchor_pos - 1 ) );
			Residue const & downstream_rsd( pose.residue( anchor_pos + 1 ) );
			return ( upstream_rsd.xyz("CA") - downstream_rsd.xyz("CA") ).normalized();
	}


	/////////////////////////////////////////////////////////////////////////////////
	Vector const RbSegments::Z_axis(
																	core::pose::Pose & pose
																	)

	{
		using namespace conformation;
		assert( this->num_rb() == 1 );

		//The Z-axis is chosen as a vector connecting CAs anchor_pos and flex_jump_pos
		int const anchor_pos( this->anchor_res() );
		int const flex_jump_pos( get_flex_jump_pos( this->seg_start(), this->seg_stop() ) );
		assert( flex_jump_pos != 0 );

		Residue const & rsd_anchor( pose.residue( anchor_pos ) );
		Residue const & rsd_flex_jump( pose.residue( flex_jump_pos ) );

		return ( rsd_anchor.xyz("CA") - rsd_flex_jump.xyz("CA") ).normalized();

	}

	/////////////////////////////////////////////////////////////////////////////////
	Vector const RbSegments::Y_axis(
																	core::pose::Pose & pose
																	)
	{
		assert( this->num_rb() == 1 );

		//The Y-axis is orthogonal X and Z axes
		numeric::xyzVector < Real > Z_axis;
		Z_axis = cross( this->X_axis( pose ), this->Z_axis( pose ) );
		return Z_axis;
	}


	/////////////////////////////////////////////////////////////////////////////////
Vector const RbSegments::alt_X_axis(
	core::pose::Pose & pose
)

{
	using namespace conformation;

	assert( this->num_rb() == 1 );
	//The X-axis is assigned as the helix axis vector. This is approximated as the vector
	//connecting i - 3 and i + 3 residues. Where i is the flex_jump_pos
	//The Z-axis is set as the line connecting the anchor_pos and flex_jump_pos
	//The Y-axis is calculated as X cross Z.

	int const flex_jump_pos( get_flex_jump_pos( this->seg_start(), this->seg_stop() ) );

	if ( pose.secstruct( flex_jump_pos ) != 'H' && pose.secstruct( flex_jump_pos ) != 'E' ) {
		Error() <<"The flex_jump_pos has to be in a strand or helix " << flex_jump_pos << " SS type " << pose.secstruct( flex_jump_pos ) << "\n";
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	//Now define X-axis for helix
	if( pose.secstruct( flex_jump_pos ) == 'H' ) {
		if( pose.secstruct( flex_jump_pos - 3 ) == 'H' && pose.secstruct( flex_jump_pos + 3 ) == 'H' ) {
			Residue const & upstream_rsd( pose.residue( flex_jump_pos - 3 ) );
			Residue const & downstream_rsd( pose.residue( flex_jump_pos + 3 ) );
			return ( upstream_rsd.xyz("CA") - downstream_rsd.xyz("CA") ).normalized();
		} else if ( pose.secstruct( flex_jump_pos - 3 ) == 'H' && pose.secstruct( flex_jump_pos + 3 ) != 'H' ) {
			Warning() << "The flex_jump_pos may be close to the edge of the helix"
								<< "This might not truly represent the axis of the helix " << "\n";
			Residue const & upstream_rsd( pose.residue( flex_jump_pos - 3) );
			Residue const & downstream_rsd( pose.residue( flex_jump_pos ));
			return ( upstream_rsd.xyz("CA") - downstream_rsd.xyz("CA") ).normalized();
		} else if ( pose.secstruct( flex_jump_pos - 3 ) != 'H' && pose.secstruct( flex_jump_pos + 3 ) == 'H' ) {
			Warning() << "The flex_jump_pos may be close to the edge of the helix"
								<< "This might not truly represent the axis of the helix " << "\n";
			Residue const & upstream_rsd( pose.residue( flex_jump_pos ) );
			Residue const & downstream_rsd( pose.residue( flex_jump_pos + 3 ));
			return ( upstream_rsd.xyz("CA") - downstream_rsd.xyz("CA") ).normalized();
		} else {
			Error() << "This is a single turn helix. Why do you want to rigid body perturb this"
							<< "You are better off modeling this as a loop" << "\n";
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}
	assert( pose.secstruct( flex_jump_pos ) == 'E' );
	//Now define X-axis for strand
	if( pose.secstruct( flex_jump_pos - 1 ) != 'E' && pose.secstruct( flex_jump_pos + 1 ) != 'E' ) {
		Error() << "This is a very short strand. You are better off modeling this as loop " << "\n";
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}
	Residue const & upstream_rsd( pose.residue( flex_jump_pos - 1 ) );
	Residue const & downstream_rsd( pose.residue( flex_jump_pos + 1 ) );
	return ( upstream_rsd.xyz("CA") - downstream_rsd.xyz("CA") ).normalized();
}

	/////////////////////////////////////////////////////////////////////////////////

Vector const RbSegments::alt_Y_axis(
	core::pose::Pose & pose
)
{
	assert( this->num_rb() == 1 );
	//The Y-axis is orthogonal to X and Z axes
	numeric::xyzVector < Real > Y_axis( cross( this->alt_X_axis( pose ),
																			this->alt_Z_axis( pose ) ) );
	return Y_axis;
}

	/////////////////////////////////////////////////////////////////////////////////

Vector const RbSegments::alt_Z_axis(
	core::pose::Pose & pose
)
{
	using namespace conformation;
	assert( this->num_rb() );

	//The Z-axis is chosen as the vector connecting CAs of anchor_pos and flex_jump_pos
	int const anchor_pos( this-> anchor_res() );
	int const flex_jump_pos( get_flex_jump_pos( this->seg_start(), this->seg_stop() ) );

	Residue const & rsd_anchor( pose.residue( anchor_pos ) );
	Residue const & rsd_flex_jump( pose.residue( flex_jump_pos ) );

	return ( rsd_anchor.xyz("CA") - rsd_flex_jump.xyz("CA") ).normalized();
}

	/////////////////////////////////////////////////////////////////////////////////
}//SSrbrelax
}//devel
