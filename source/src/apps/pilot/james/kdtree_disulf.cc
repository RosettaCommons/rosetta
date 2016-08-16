// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file kdtree_disulf.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/RT.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/robert.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

#include <basic/database/open.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <ObjexxFCL/string.functions.hh>

#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>

#include <numeric/random/random.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzMatrix.fwd.hh>

// C++ headers
#include <iostream>
#include <string>

//Auto Headers
#include <protocols/jumping/StrandPairing.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>

#include <numeric/kdtree/util.hh>
#include <numeric/kdtree/calc_distance.hh>
#include <numeric/kdtree/nearest_neighbors.hh>
#include <numeric/kdtree/KDTree.hh>
#include <numeric/kdtree/KDPoint.hh>
#include <numeric/kdtree/KDPointList.hh>
#include <numeric/kdtree/WrappedPrimitive.hh>

#include <utility/excn/Exceptions.hh>


///////////////////////////////////////////////////////////////////////////////

using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using ObjexxFCL::FArray2A_float;
using ObjexxFCL::FArray2D_float;
using namespace core;
using namespace pose;
using namespace conformation;

typedef numeric::kdtree::WrappedPrimitive< core::kinematics::RT > WrappedRT;
typedef utility::pointer::owning_ptr< WrappedRT > WrappedRTOP;

int
main( int argc, char* argv [] ) {
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// options, random initialization
	devel::init( argc, argv );

	utility::io::izstream disulfide_database;
	basic::database::open(disulfide_database, "sampling/disulfide_jump_database_wip.dat");

	using std::string;
	using utility::vector1;
	using core::kinematics::RT;
	using core::kinematics::RTOP;
	using utility::pointer::ReferenceCountOP;

	string line;
	vector1< ReferenceCountOP > rts;
	getline(disulfide_database, line);
	while ( !disulfide_database.eof() ) {
		std::istringstream line_stream(line);
		string junk;
		//ReferenceCountOP rt_op( new RT );

		RT rt;
		line_stream >> junk >> junk >> junk >> junk >> junk
			>> junk >> junk >> junk >> junk >> rt;

		//ReferenceCountOP rt_op( rt );
		rts.push_back( new WrappedRT( rt ) );

		getline(disulfide_database, line);
	}

	// setup kdtree. db_points is the vector in k-space that will be searched,
	// each point is an (x,y,z) translation associated with an RT.
	using namespace numeric::kdtree;
	vector1< vector1< Real > > db_points;

	typedef vector1< ReferenceCountOP >::const_iterator iter;
	for ( iter it = rts.begin(), end = rts.end(); it != end; ++it ) {
		WrappedRTOP wrt = dynamic_cast< WrappedRT * > ( (*it)() );
		RT rt = wrt->val();
		numeric::xyzVector< Real > translation = rt.get_translation();
		vector1< Real > pt( 3, 0.0 );
		pt[1] = translation.x();
		pt[2] = translation.y();
		pt[3] = translation.z();
		db_points.push_back( pt );
	}
	vector1< Real > search_pt = db_points.front();
	search_pt[1] =  0;
	search_pt[2] = -2.5;
	search_pt[3] = -2.3;

	KDTree tree( db_points, rts );

	Size const max_nn( 2000 );
	Real const max_dist( option[ james::real ]() );
	Real const max_dist_sq( max_dist * max_dist );

	// naive search for nn
	utility::io::ozstream out_naive( "naive.txt" );
	clock_t start_naive = clock();
	vector1< ReferenceCountOP > found_pts;
	typedef vector1< vector1< Real > >::const_iterator pt_iter;
	pt_iter pt_it, pt_end;
	for ( iter it = rts.begin(), end = rts.end(); it != end; ++it ) {
		WrappedRTOP wrt = dynamic_cast< WrappedRT * > ( (*it)() );
		RT rt = wrt->val();
		numeric::xyzVector< Real > translation = rt.get_translation();
		vector1< Real > pt( 3, 0.0 );
		pt[1] = translation.x();
		pt[2] = translation.y();
		pt[3] = translation.z();

		Real const dist_sq( sq_vec_distance( search_pt, pt ) );
		if ( dist_sq < max_dist_sq ) {
			found_pts.push_back( wrt );
			out_naive << "found " << rt << " at distance " << dist_sq << std::endl;
		} else {
			out_naive << "rejected " << rt << " at distance " << dist_sq << std::endl;
		}
	}
	clock_t finish_naive = clock();

	for ( iter it = found_pts.begin(), end = found_pts.end(); it != end; ++it ) {
		WrappedRTOP wrt = dynamic_cast< WrappedRT * > ( (*it)() );
		RT rt = wrt->val();
		//std::cout << "found " << rt << std::endl;
		//out_naive << "found " << rt << std::endl;
	}
	out_naive.close();

	clock_t start_kdtree  = clock();
	KDPointList neighbors = nearest_neighbors(
		tree, search_pt, max_nn, max_dist
	);
	clock_t finish_kdtree = clock();

	utility::io::ozstream out_kdtree( "kdtree.txt" );
	for ( vector1< KDPointOP >::const_iterator
			it = neighbors.begin(), end = neighbors.end(); it != end; ++it
	) {
		WrappedRTOP wrt = dynamic_cast< WrappedRT * > ( (*it)->data()() );
		RT rt = wrt->val();
		out_kdtree << "found " << rt << " at distance " << (*it)->distance()
			<< std::endl;
	}
	out_kdtree.close();

	double kdtree_time = ( double(finish_kdtree) - double(start_kdtree) ) / CLOCKS_PER_SEC;
	double naive_time  = ( double(finish_naive)  - double(start_naive) )  / CLOCKS_PER_SEC;

	//std::cout << "kdtree: " << kdtree_time << std::endl;
	//std::cout << "naive: "  << naive_time << std::endl;
	std::cout << kdtree_time << ' ' << naive_time << ' ' << max_dist << std::endl;
	//std::cout << kdtree_time << ' ' << max_dist << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
