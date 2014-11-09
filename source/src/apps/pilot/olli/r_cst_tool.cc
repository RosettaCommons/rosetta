// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

#include <protocols/abinitio/Templates.hh>
// AUTO-REMOVED #include <protocols/abinitio/TemplateJumpSetup.hh>
#include <protocols/abinitio/PairingStatistics.hh>
#include <protocols/abinitio/StrandConstraints.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/constraints_additional/AdditionalConstraintCreators.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>

#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/util.hh>
#include <protocols/loops/Loops.hh>
#include <core/chemical/ResidueType.hh>
// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <utility/excn/EXCN_Base.hh>


#include <numeric/random/random.hh>



static thread_local basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
using namespace abinitio;

using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace scoring;
using namespace scoring::constraints;


OPT_KEY( Integer, level ) //viol_level already in AbrelaxApplication definier
OPT_1GRP_KEY( File, in, top )
OPT_KEY( Real, threshold )
OPT_KEY( Boolean, no_out )
OPT_KEY( Boolean, clean )
OPT_1GRP_KEY( File, out, cst )
OPT_KEY( Real, skip )
OPT_KEY( Boolean, show_atom_numbers )
OPT_KEY( Boolean, skip_redundant )
OPT_KEY( Integer, skip_redundant_width )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
  Templates::register_options();
  OPT( in::file::s );
  OPT( constraints::cst_file );
  OPT( out::prefix );
	NEW_OPT( threshold, "what makes a violation", 1 );
  NEW_OPT( level, "how much detail for violation output", 71 );
	NEW_OPT( in::top, "read topology from this file for checking", "");
	NEW_OPT( no_out, "don't write individual .viol files", false);
  //NEW_OPT( viol_type, "work only on these types of constraints", "");
	NEW_OPT( out::cst, "write the resulting constraint set to this file", "scanned.cst");
	NEW_OPT( clean, "write only PASSED constraints to file ( out::cst )", false );
	NEW_OPT( skip, "write only X-th part of constraints (default 100percent: 1.0)", 1.0 );
	NEW_OPT( show_atom_numbers, "dump information on all residues ", false );
	NEW_OPT( skip_redundant, "skip redundant constraints ( util.cc ) ", false );
	NEW_OPT( skip_redundant_width, "skip residue pairs in X neighbourhood ", 1);
	//	NEW_OPT( combine, "combine constraints randomly with OR ", 1 );
	OPT( constraints::combine );
}

// Forward
class ConstraintToolMover;

// Types
typedef  utility::pointer::shared_ptr< ConstraintToolMover >  ConstraintToolMoverOP;
typedef  utility::pointer::shared_ptr< ConstraintToolMover const >  ConstraintToolMoverCOP;

class ConstraintToolMover : public moves::Mover {
public:
	virtual void apply( core::pose::Pose& );
	ConstraintSet const& cstset() const {
		return *cstset_;
	}
	void show_cstset( std::ostream& os ) const {

		ConstraintSetOP output_cst = cstset_;
		if ( option[ skip ] < 1.0 ) {
			ConstraintCOPs all_constraints = cstset_->get_all_constraints();
			ConstraintSetOP filtered_cstset( new ConstraintSet );
			for ( ConstraintCOPs::const_iterator it = all_constraints.begin(); it != all_constraints.end(); ++it ) {
				Real r = numeric::random::rg().uniform();
				if ( r < option[ skip ] ) {
					filtered_cstset->add_constraint( *it );
				}
			}
			output_cst=filtered_cstset;
		}

		output_cst->show_definition( os, reference_pose_ );
	}
	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return protocols::moves::MoverOP( new ConstraintToolMover() );
	}

	std::string get_name() const { return "ConstraintToolMover"; }

private:
	ConstraintSetOP cstset_;
	pose::Pose reference_pose_; //this pose has been used to read constraint set
};

void ConstraintToolMover::apply( core::pose::Pose &pose ) {
	if ( option[ OptionKeys::show_atom_numbers ]() ) {
		tr.Info << pose.annotated_sequence() << std::endl;
		for ( Size i =1; i<=pose.total_residue(); ++i ) {
			tr.Info << " residue " << i << std::endl;
			pose.residue_type( i ).show_all_atom_names( tr.Info );
		}
	}
	if ( !cstset_ && option[ OptionKeys::constraints::cst_file ].user() ) {
		// reads and sets constraints
		tr.Info << "read constraints... : " << std::endl;
		cstset_ = ConstraintIO::get_instance()->read_constraints( option[ OptionKeys::constraints::cst_file ]()[1], ConstraintSetOP( new ConstraintSet ), pose );
		ConstraintCOPs added_constraints = cstset_->get_all_constraints();
		utility::vector1< bool > exclude_res;
		//if ( option[ OptionKeys::constraints::combine_exclude_region ].user() ) {
		//	std::string const file( option[ OptionKeys::constraints::combine_exclude_region ]() );
		//	loops::Loops rigid_core;
		//	rigid_core.read_loop_file( file, false /*no strict looprlx checking*/, "RIGID" );
		//	exclude_res.resize( pose.total_residue(), false );
		//	rigid_core.transfer_to_residue_vector( exclude_res, true );
		//}
		if ( option[ OptionKeys::skip_redundant ]() ) skip_redundant_constraints( added_constraints, pose.total_residue(), option[ OptionKeys::skip_redundant_width ]() );

		core::kinematics::ShortestPathInFoldTree sp( pose.fold_tree() );
		combine_constraints( added_constraints, option[ OptionKeys::constraints::combine ](), exclude_res, sp ); // if combine_ratio_ > 1 this will randomly combine constraints into multipletts w
		cstset_ = ConstraintSetOP( new ConstraintSet );
		cstset_->add_constraints( added_constraints );
		reference_pose_ = pose;
	}
	if ( !cstset_ && option[ OptionKeys::in::top ] ) {
		ConstraintCOPs my_strand_cst;
		utility::io::izstream is( option[ OptionKeys::in::top ] );
		PairingStatisticsOP ps( new PairingStatistics );
		is >> *ps;
		tr.Info << *ps << std::endl;
		my_strand_cst = StrandConstraints( *ps ).build_constraints( pose );
		cstset_ = ConstraintSetOP( new ConstraintSet );
		cstset_->add_constraints( my_strand_cst );
		reference_pose_ = pose;
	}
	cstset_->show_definition( tr.Debug, pose );
	pose.constraint_set( cstset_ );
	scoring::ScoreFunction score;
	score.set_weight( scoring::atom_pair_constraint, 1.0 );
	score.set_weight( scoring::coordinate_constraint, 1.0 );

	std::string fname( jd2::current_output_name() );
	tr.Info << fname << " " << score( pose ) << std::endl;
	if ( !option[ no_out ]() ) {
		std::string out_name( option[ out::prefix ]() + fname + ".viol" );
		utility::io::ozstream vfile( out_name );
		cstset_->show_violations( vfile, pose, option[ level ](), option[ threshold ] );
	} else {
		cstset_->show_violations( tr.Debug, pose, option[ level ](), option[ threshold ] );
	}

	if ( option[ clean ]() ) {  //cleaning
		ConstraintCOPs all_constraints = cstset_->get_all_constraints();
		ConstraintSetOP filtered_cstset( new ConstraintSet );
		for ( ConstraintCOPs::const_iterator it = all_constraints.begin(); it != all_constraints.end(); ++it ) {
			if ( (*it)->show_violations( tr.Debug, pose, option[ level ](), option[ threshold ] ) == 0 ) {
				filtered_cstset->add_constraint( *it );
			}
		}
		cstset_=filtered_cstset;
	}
}



void run() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
		basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/olange/minirosetta_database");
	}

	if ( option[ OptionKeys::constraints::no_linearize_bounded ] ) {
		tr.Info << "use fully harmonic potential for BOUNDED " << std::endl;
		ConstraintIO::get_func_factory().add_type("BOUNDED", scoring::func::FuncOP( new scoring::constraints::BoundFunc(0,0,0,1000,"dummy") ) );
	}
	if ( option[ OptionKeys::constraints::named ] ) {
		tr.Info << "use named constraints in AtomPairConstraint to avoid problems with cutpoint-variants " << std::endl;
		//ConstraintIO::get_cst_factory().add_type( new scoring::constraints::NamedAtomPairConstraint( id::NamedAtomID(), id::NamedAtomID(), NULL) );
		/// WARNING WARNING WARNING. THREAD UNSAFE. DO NOT USE SINGLETONS THIS WAY.
		core::scoring::constraints::ConstraintFactory::get_instance()->replace_creator( ConstraintCreatorCOP( ConstraintCreatorOP( new constraints_additional::NamedAtomPairConstraintCreator ) ) );
	}

	ConstraintToolMoverOP cst_tool( new ConstraintToolMover );
	protocols::jd2::JobDistributor::get_instance()->go( cst_tool, jd2::JobOutputterOP( new jd2::NoOutputJobOutputter ) );


	if ( option[ out::cst ].user() ) {
		utility::io::ozstream vfile( option[ out::cst ]() );
		cst_tool->show_cstset( vfile );
	} else {
		cst_tool->show_cstset( std::cout );
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


