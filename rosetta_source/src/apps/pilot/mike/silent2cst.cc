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


// libRosetta headers

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>

using namespace core;
using namespace fragment;
using namespace protocols;
using namespace utility;
using namespace core;

namespace protocols {
namespace moves {

class SimpleCstMover : public moves::Mover
{
public:
	/// @brief
	/// 	empty constructor fills values with the values
	///		read in from the commandline
	SimpleCstMover( Size length) :
		Mover()
	{
		Mover::type( "SimpleCstMover" );
		dist_matrix = new Real[length * length * 100];
		for( Size i=0; i < length*length; i ++ ) dist_matrix[i] = 0;
	}


	~SimpleCstMover(){
		delete [] dist_matrix;
		delete [] dist2_matrix;

	}

	virtual void apply( pose::Pose & pose );
	virtual void test_move( pose::Pose & pose ){apply(pose);}

	void write( const std::string &filename );


private:  //  Data

	inline void set_dist(Size i, Size j, Real value ){    dist_matrix[j  *length + i ] = value; }
	inline void add_dist(Size i, Size j, Real value ){    dist_matrix[j  *length + i ] += value; }
	Real        get_dist(Size i, Size j ) const { return  dist_matrix[j  *length + i ]; }
	inline void add_dist2(Size i, Size j, Real value ){   dist2_matrix[j *length + i ] += value; }
	Real        get_dist2(Size i, Size j ) const { return dist2_matrix[j *length + i ]; }

	Size length;

	Real   *dist_matrix;
	Real  *dist2_matrix;

};

typedef utility::pointer::owning_ptr< SimpleCstMover > SimpleCstMoverOP;
typedef utility::pointer::owning_ptr< SimpleCstMover const > SimpleCstMoverCOP;




void SimpleCstMover::apply( pose::Pose & pose )
{
	std::cout << "Stealing .. " << std::endl;
	static int count=0;

	for( Size i = 1;  i < pose.total_residue(); i ++ ){
		for( Size j = i+3;  j < pose.total_residue(); j ++ ){
			Real dist;
			std::cout << "C " << count << " I " << i << " J " << j << " D " << dist << std::endl;
 			dist = numeric::distance( pose.xyz(  core::id::AtomID( 3, i ) ) , pose.xyz(  core::id::AtomID( 3, j  ) ) );
			add_dist(i,j,dist);
			add_dist2(i,j,dist*dist );
			std::cout << "C " << count << " I " << i << " J " << j << " D " << dist << std::endl;
		}
	}
 count++;



}

void SimpleCstMover::write( const std::string &filename )
{

	// first calculate the mean.





}

} // moves
} // protocols



///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	devel::init(argc, argv);
	using namespace protocols::moves;
	SimpleCstMoverOP capture3mers =  new SimpleCstMover( 300 );
	SequenceMoverOP seqmov = new SequenceMover;
	seqmov->add_mover( capture3mers );
	MoverOP mover = seqmov;
 	protocols::jd2::JobDistributor::get_instance()->go( *seqmov );



	return 0;
}

