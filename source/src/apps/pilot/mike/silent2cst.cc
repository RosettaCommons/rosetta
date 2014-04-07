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
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyzVector.hh>

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
		Mover(),
		length_( length )
	{
		Mover::type( "SimpleCstMover" );
		distmin_matrix = new Real[length * length ];
		distmax_matrix = new Real[length * length ];
		for( Size i=0; i < length*length; i ++ ) distmin_matrix[i] = 1000000.0;
		for( Size i=0; i < length*length; i ++ ) distmax_matrix[i] = 0;
	}


	~SimpleCstMover(){
		delete [] distmin_matrix;
		delete [] distmax_matrix;
	}

	virtual void apply( pose::Pose & pose );
	virtual void test_move( pose::Pose & pose ){apply(pose);}
 	virtual std::string get_name() const { return "MyMover"; }

	void write( const std::string &filename );


private:  //  Data

	inline void set_distmin(Size i, Size j, Real value ){    distmin_matrix[j  *length + i ] = value; }
	inline void add_distmin(Size i, Size j, Real value ){    distmin_matrix[j  *length + i ] += value; }
	Real        get_distmin(Size i, Size j ) const { return  distmin_matrix[j  *length + i ]; }

	inline void set_distmax(Size i, Size j, Real value ){    distmax_matrix[j  *length + i ] = value; }
	inline void add_distmax(Size i, Size j, Real value ){   distmax_matrix[j *length + i ] += value; }
	Real        get_distmax(Size i, Size j ) const { return distmax_matrix[j *length + i ]; }

	Size length_;

	Real  *distmin_matrix;
	Real  *distmax_matrix;

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
 			dist = pose.xyz(  core::id::AtomID( 3, i ) ).distance( pose.xyz(  core::id::AtomID( 3, j  ) ) );
			if ( dist < get_distmin( i, j ) ) set_distmin( i, j, dist );
			if ( dist > get_distmax( i, j ) ) set_distmax( i, j, dist );
			std::cout << "C " << count << " I " << i << " J " << j << " D " << distmin << std::endl;
		}
	}
 count++;

}

void SimpleCstMover::write( const std::string &filename )
{
  Real dist_limit = 12;
	for( Size i = 1;  i < pose.total_residue(); i ++ ){
		for( Size j = i+3;  j < pose.total_residue(); j ++ ){

			if( get_distmax( i, j ) < 12 ){
				if( get_distmax( i, j ) > get_distmin( i, j ) ) { std::cerr "ASSERTION"; }


			}

			if ( dist < get_distmin( i, j ) ) set_distmin( i, j, dist );
			if ( dist > get_distmax( i, j ) ) set_distmax( i, j, dist );
			std::cout << "C " << count << " I " << i << " J " << j << " D " << distmin << std::endl;
		}
	}


}

} // moves
} // protocols



///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
    try {
    	devel::init(argc, argv);
    	using namespace protocols::moves;
    	SimpleCstMoverOP capture3mers =  new SimpleCstMover( 300 );
    	SequenceMoverOP seqmov = new SequenceMover;
    	seqmov->add_mover( capture3mers );
    	MoverOP mover = seqmov;
    	using namespace protocols::jd2;
    	try{
    		JobDistributor::get_instance()->go( mover );
    	} catch ( utility::excn::EXCN_Base& excn ) {
    		std::cerr << "Exception: " << std::endl;
    		excn.show( std::cerr );
    		std::cout << "Exception: " << std::endl;
    		excn.show( std::cout ); //so its also seen in a >LOG file
    	}
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}

