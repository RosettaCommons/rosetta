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
/// @author

#include <core/conformation/Conformation.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FragData.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>


#include <core/fragment/util.hh>
#include <basic/Tracer.hh>

#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <core/fragment/SecondaryStructure.hh>

#include <core/scoring/rms_util.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>


#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>

#include <basic/options/option.hh>

#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option_macros.hh>


#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


#include <protocols/toolbox/Cluster.hh>

// option key includes

#include <basic/options/keys/cluster.OptionKeys.gen.hh>
//#include <basic/options/keys/OptionKeys.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <utility/excn/Exceptions.hh>


static THREAD_LOCAL basic::Tracer tr( "main" );


class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}

OPT_KEY( File, f )
OPT_KEY( Boolean, torsion )
OPT_KEY( Integer, chop )
OPT_KEY( File, jumpss )
OPT_KEY( File, fold_tree )
OPT_KEY( Integer, ref_size )
OPT_KEY( File, exclude )
OPT_KEY( File, filter )
OPT_KEY( File, write )
OPT_KEY( Boolean, intrinsic )
OPT_KEY( File, write_big_cluster)
OPT_KEY( Integer, min_chop_in_quality_check )
OPT_KEY( File, fill_frags)
OPT_KEY( Integer, cluster_size)
OPT_1GRP_KEY( File, cluster, out )
OPT_1GRP_KEY( File, cluster, in )
OPT_1GRP_KEY( Real, cluster, acc_size )
OPT_1GRP_KEY( Real, cluster, acc_sum )
OPT_1GRP_KEY( Real, cluster, nmax )
OPT_1GRP_KEY( IntegerVector, cluster, range )
OPT_1GRP_KEY( File, out, qual )
OPT_KEY( File, ss_content )

using namespace core;
using namespace fragment;
using namespace pose;
using namespace kinematics;


using namespace protocols::jumping;
using namespace protocols::toolbox;


using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace ObjexxFCL::format;

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::native );
	NEW_OPT( f, "fragment file", "" );
	NEW_OPT( out::qual, "sequence_quality", "frag_sequence_quality.dat" );
}


int main( int argc, char** argv ) {
	try {

	ThisApplication::register_options();
	devel::init( argc, argv );
	Pose native;
	std::string const native_pdb ( option[ in::file::native ]() );
	core::import_pose::pose_from_file( native, native_pdb , core::import_pose::PDB_file);
	core::util::switch_to_residue_type_set( native, chemical::CENTROID );

	Pose test_pose(native);
    utility::io::ozstream out( option[ out::qual ] );
    FragSetOP orig_frags;
    orig_frags = FragmentIO().read_data( option[ f ]() );
    Size max_pdb( 0 );
    typedef boost::tuple<Real, std::string, Size, std::string, std::string, char, Size> FragRMSDAndSequence;
    for ( FrameIterator frame = orig_frags->begin(), eframe=orig_frags->end(); frame != eframe; ++frame ) {
        assert( frame->length() == 9);
        out << "\nstart: " << LJ(5,frame->start())  << "end: " << LJ(5,frame->end())<< std::endl;
				out << RJ(16,"rmsd") << RJ(5,"sequence") <<RJ(15,"secstruct") <<RJ(15,"pdb") <<RJ(8,"chain")<< RJ(6,"start")<< std::endl;
        std::set<FragRMSDAndSequence> data;
        for ( Size i=1; i<=frame->nr_frags(); i++ ) {
            frame->apply( i, test_pose );
            data.insert(boost::make_tuple(scoring::CA_rmsd( native, test_pose, frame->start(), frame->end()),frame->fragment_ptr(i)->sequence(),i,frame->fragment_ptr(i)->secstruct(), frame->fragment_ptr(i)->pdbid(),frame->fragment_ptr(i)->chain(), frame->fragment_ptr(i)->pdbpos()  ));
        }
				for(std::set<FragRMSDAndSequence>::iterator fd = data.begin(); fd!=data.end(); fd++){
						out << RJ(16,boost::get<0>(*fd)) << RJ(5,boost::get<2>(*fd)) <<RJ(15,boost::get<1>(*fd)) <<RJ(15,boost::get<3>(*fd)) <<RJ(8,boost::get<4>(*fd))<< RJ(6,boost::get<5>(*fd))<< std::endl;
				}
    }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}
