// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/types.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobInputter.hh>

// AUTO-REMOVED #include <protocols/jd2/util.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>
// AUTO-REMOVED #include <protocols/toolbox/DecoySetEvaluation.impl.hh>
#include <protocols/toolbox/InteratomicVarianceMatrix.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/toolbox/Cluster.hh>
#include <protocols/toolbox/Cluster.impl.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
// AUTO-REMOVED #include <core/scoring/constraints/BoundConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/NamedAtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/LocalCoordinateConstraint.hh>


// AUTO-REMOVED #include <core/id/NamedStubID.hh>
#include <core/id/AtomID.hh>

#include <core/chemical/ChemicalManager.hh>


// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/excn/Exceptions.hh>

// option key includes
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

// ObjexxFCL includes
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <ostream>
// AUTO-REMOVED #include <iterator>
#include <algorithm>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <utility/excn/EXCN_Base.hh>


static thread_local basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace toolbox;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

OPT_KEY( Real, wRMSD )
OPT_KEY( Real, tolerance )
OPT_KEY( Boolean, dump_fit )
OPT_1GRP_KEY( File, rmsf, out )
OPT_1GRP_KEY( File, rigid, out )
OPT_1GRP_KEY( Real, rigid, cutoff )
OPT_1GRP_KEY( Integer, rigid, min_gap )
OPT_1GRP_KEY( File, ivm, out )
OPT_1GRP_KEY( Integer, ivm, grid )
OPT_1GRP_KEY( RealVector, ivm, bounds )
OPT_1GRP_KEY( File, dist_cst, out )
OPT_1GRP_KEY( File, coord_cst, out )
OPT_2GRP_KEY( File, out, file, cluster )
OPT_KEY( Boolean, clustering )
OPT_1GRP_KEY( File, rigid, in )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
  OPT( in::file::s );
	OPT( in::file::silent );
	OPT( in::file::native );
	OPT( cluster::limit_cluster_size );
	NEW_OPT( wRMSD, "compute wRMSD with this sigma", 1 );
	NEW_OPT( tolerance, "stop wRMSD iteration if <tolearance change in wSum", 0.00001);
	NEW_OPT( dump_fit, "outptu pdbs of fitted structures ", false );
	NEW_OPT( ivm::out, "compute ivm", "ivm.dat" );
	NEW_OPT( ivm::grid, "number of grid points", 100 );
	utility::vector1< core::Real > bounds;
	bounds.push_back( 0.01 );
	bounds.push_back( 2 );
	NEW_OPT( rigid::in, "residues that are considered for wRMSD or clustering", "rigid.loop");
	NEW_OPT( ivm::bounds, "lower and upper bound for epsilon", bounds );
	NEW_OPT( rmsf::out, "write rmsf into this file", "rmsf.dat" );
	NEW_OPT( rigid::out, "write a RIGID definition", "rigid.loops" );
	NEW_OPT( rigid::cutoff, "residues with less then loop_cutoff will go into RIGID definition", 1.0 );
	NEW_OPT( rigid::min_gap, "don't have less then min_gap residues between rigid regions", 5 );
	NEW_OPT( dist_cst::out, "write constraint set with distances derived from decoy set", "dist.cst" );
	NEW_OPT( coord_cst::out, "write constraint set with bounds derived from superimposed structures","coord.cst" );
	NEW_OPT( out::file::cluster, "write clustered structures to silent file with this name", "cluster.out" );
	NEW_OPT( clustering, "compute RMSD distance matrix and cluster", false );
	DecoySetEvaluation::register_options();
}

// Forward
class RmsfMover;

// Types
typedef  utility::pointer::shared_ptr< RmsfMover >  RmsfMoverOP;
typedef  utility::pointer::shared_ptr< RmsfMover const >  RmsfMoverCOP;

class RmsfMover : public moves::Mover {
public:
	virtual void apply( core::pose::Pose& );
	std::string get_name() const { return "RmsfMover"; }

	DecoySetEvaluation eval_;
};

void RmsfMover::apply( core::pose::Pose &pose ) {
	if ( eval_.n_decoys_max() < eval_.n_decoys() + 1 ) {
		eval_.reserve( eval_.n_decoys() + 100 );
	}
	eval_.push_back( pose );
}

class FitMover : public moves::Mover {
public:
	virtual void apply( core::pose::Pose& );
	std::string get_name() const { return "FitMover"; }
	FitMover() : first( true ), iref( 1 ) {};
	ObjexxFCL::FArray1D_double weights_;
	core::pose::Pose ref_pose_;
	bool first;
	Size iref;
};

typedef utility::pointer::shared_ptr< FitMover > FitMoverOP;

void FitMover::apply( core::pose::Pose &pose ) {
	if ( iref && first ) {
		--iref;
		return;
	}
	if ( first ) {
		ref_pose_ = pose;
		first = false;
		return;
	}
	runtime_assert( !first );
	CA_superimpose( weights_, ref_pose_, pose );
	//pose.dump_pdb( "FIT_"+get_current_tag()+".pdb" );
}

void run() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
		basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/olange/minirosetta_database");
	}

	bool store_energies ( option[ clustering ] );

	RmsfMoverOP rmsf_tool( new RmsfMover );
	protocols::jd2::SilentFileJobInputterOP sfd_inputter ( utility::pointer::dynamic_pointer_cast< protocols::jd2::SilentFileJobInputter > ( protocols::jd2::JobDistributor::get_instance()->job_inputter() ) );
	if ( sfd_inputter ) { //option[ in::file::silent ].user() ) {
		io::silent::SilentFileData const& sfd( sfd_inputter->silent_file_data() );
		//		sfd.read_file( option[ in::file::silent ]
		rmsf_tool->eval_.push_back_CA_xyz_from_silent_file( sfd, store_energies );
	} else {
		protocols::jd2::JobDistributor::get_instance()->go( rmsf_tool, jd2::JobOutputterOP( new jd2::NoOutputJobOutputter ) );
	}

	FArray1D_double weights( rmsf_tool->eval_.n_atoms(), 1.0 );
	FArray1D_double input_weights( rmsf_tool->eval_.n_atoms(), 1.0 );

	if ( option[ rigid::in ].user() ) {
			std::ifstream is( option[ rigid::in ]().name().c_str() );

			if (!is.good()) {
				utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + option[ rigid::in ]().name() + "'" );
			}

			loops::PoseNumberedLoopFileReader reader;
			reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
			loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file(is, option[ rigid::in ](), false );
			loops::Loops rigid = loops::Loops( loops );

			for ( Size i=1;i<=rmsf_tool->eval_.n_atoms(); ++i ) {
				if (rigid.is_loop_residue( i ) ) weights( i )=1.0;
				else weights( i )=0.0;
			}
	}
	input_weights=weights;
	rmsf_tool->eval_.set_weights( weights );

	utility::vector1< Real > result;
	Size icenter = 1;
	if ( option[ rmsf::out ].user() || option[ coord_cst::out ].user() || option[ rigid::out ].user() || option[ dump_fit ] )	{//output rmsf file
		rmsf_tool->eval_.superimpose();
		if ( option[ wRMSD ].user() ) {
			icenter = rmsf_tool->eval_.wRMSD( option[ wRMSD ], option[ tolerance ](), weights );
		}
		rmsf_tool->eval_.rmsf( result );
		if ( option[ rmsf::out ].user() ) {
			utility::io::ozstream os_rmsf( option[ rmsf::out ]()); //output rmsf file
			Size ct = 1;
			for ( utility::vector1<Real>::const_iterator it = result.begin(); it != result.end(); ++it, ++ct ) {
				os_rmsf << ct << " " << *it << " " << weights( ct ) << std::endl;
			}
		}
	}

	if ( option[ rigid::out ].user() ) {//create RIGID output
		//		utility::ozstream os_core( option[ out::rigid ]() );
		Real const cutoff ( option[ rigid::cutoff ] );
		tr.Info << "make rigid with cutoff " << cutoff << " and write to file... " << option[ rigid::out ]() << std::endl;
		loops::Loops rigid;
		for ( Size i=1; i<=result.size(); ++i ) {
			if ( result[ i ]<cutoff && input_weights( i ) > 0 ) rigid.add_loop( loops::Loop(  i, i ), option[ rigid::min_gap ]() );
		}
		tr.Info << rigid << std::endl;
		rigid.write_loops_to_file( option[ rigid::out ](), "RIGID" );
		std::string fname =  option[ rigid::out ]();
		loops::Loops loops = rigid.invert( result.size() );
		loops.write_loops_to_file( fname + ".loopfile" , "LOOP" );


	}


	if ( option[ OptionKeys::clustering ]() ) {

		protocols::jd2::SilentFileJobInputterOP sfd_inputter ( utility::pointer::dynamic_pointer_cast< protocols::jd2::SilentFileJobInputter > ( protocols::jd2::JobDistributor::get_instance()->job_inputter() ) );
		if ( sfd_inputter ) { //option[ in::file::silent ].user() ) {
			io::silent::SilentFileData const& sfd( sfd_inputter->silent_file_data() );
			core::io::silent::SilentFileData kept_decoys;
			cluster_silent_structs( rmsf_tool->eval_, sfd.begin(), sfd.end(), kept_decoys, ClusterOptions( true ) );
			kept_decoys.write_all( option[ out::file::cluster ]() );
		}
	}

	if ( option[ ivm::out ].user() ) {
		InteratomicVarianceMatrix ivm;
		ivm.init( rmsf_tool->eval_.n_atoms(), rmsf_tool->eval_.n_decoys(), rmsf_tool->eval_.coords() );
		ivm.optimize_kurtosis( option[ ivm::grid ], option[ ivm::bounds ]()[ 1 ], option[ ivm::bounds ]()[ 2 ]);
	}

	if ( option[ dist_cst::out ].user() ) {
		scoring::constraints::ConstraintSet cst_set;
		rmsf_tool->eval_.create_dist_constraints(	cst_set );
		scoring::constraints::ConstraintIO::write_constraints( option[ dist_cst::out ], cst_set, rmsf_tool->eval_.ref_pose() );
	}

	if ( option[ coord_cst::out ].user() ) {
		scoring::constraints::ConstraintSet cst_set;
		rmsf_tool->eval_.create_xyz_constraints_median(
					cst_set,
					rmsf_tool->eval_.ref_pose(),
					1
		);
		scoring::constraints::ConstraintIO::write_constraints( option[ coord_cst::out ], cst_set, rmsf_tool->eval_.ref_pose() );
	}

	if ( option[ dump_fit ]() ) {
		FitMoverOP fit_tool( new FitMover );
		fit_tool->weights_ = weights;
		if ( icenter > 1 ) {
			fit_tool->iref = icenter - 1; //ignores the first iref structures before it takes the pose as reference pose.
			protocols::jd2::JobDistributor::get_instance()->restart();
			protocols::jd2::JobDistributor::get_instance()->go( fit_tool, jd2::JobOutputterOP( new jd2::NoOutputJobOutputter ) );
		}
		protocols::jd2::JobDistributor::get_instance()->restart();
		protocols::jd2::JobDistributor::get_instance()->go( fit_tool );
		if ( option[ in::file::native ].user() ) {
				pose::Pose native_pose;
				core::import_pose::pose_from_pdb( native_pose,
					*core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ), option[ in::file::native ]() );
				fit_tool->apply( native_pose );
				native_pose.dump_pdb( "fit_native.pdb");
		}
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	register_options();
	devel::init( argc, argv );

	try{
		run();
	} catch ( utility::excn::EXCN_Base& excn ) {
		excn.show( std::cerr );
	}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


