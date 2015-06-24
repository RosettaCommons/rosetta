// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilot/doug/PeptoidRotamerRecoverer.cc
/// @brief Simply prints backbone and side chain dihedral angles of a peptoid
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// protocols header
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/CyclizationMover.hh>

// core headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/conformation/Residue.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>

#include <core/graph/Graph.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResiduePeptoidLibrary.hh>
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.tmpl.hh>

// devel headers
#include <devel/init.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/cyclization.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/database/open.hh>

// numeric headers
#include <numeric/angle.functions.hh>

// c++
#include <sstream>
#include <iomanip>
#include <map>

// tracer
static thread_local basic::Tracer TR( "PeptoidRotamerRecoverer" );

// local options
basic::options::BooleanOptionKey const cyclic( "cyclic" );
basic::options::IntegerVectorOptionKey const chains_to_compare( "chains_to_compare" );

// a few utility functions
core::Real
get_symm_corrected_angle( core::Size chi_num, core::conformation::Residue res )
{
	using namespace core;
	using namespace chemical;

	std::string res_type_name3( res.type().name3() );
	Real temp_chi( numeric::nonnegative_principal_angle_degrees( res.chi( chi_num ) ) );

	// if else chain of all symm side chains that we can model in the peptoid databank
	if ( res_type_name3 == "601" && chi_num == 2 && temp_chi >= 135 && temp_chi <= 315 ) {
		return temp_chi - 180;
	}	else if ( res_type_name3 == "101" && chi_num == 2 && temp_chi >= 135 && temp_chi <= 315 ) {
		return temp_chi - 180;
	}	else if ( res_type_name3 == "401" && chi_num == 3 && temp_chi >= 135 && temp_chi <= 315 ) {
		return temp_chi - 180;
	} else {
		return res.chi( chi_num );
	}
}

core::Real
angle_diff( core::Real a1, core::Real a2 )
{
	using namespace core;

	Real pad1( numeric::principal_angle_degrees( a1 ) );
	Real pad2( numeric::principal_angle_degrees( a2 ) );

	Real t1( fabs( pad1 - pad2 ) );
	Real t2( 360.00 - t1 );

	return( t1 <= t2 ? t1 : t2 );
}

core::Real
calc_dist( core::conformation::Residue res1, core::conformation::Residue res2 )
{
	using namespace core;
	using namespace conformation;

	Size nchi( res1.type().nchi() );
	Real sd( 0 );

	for( Size i( 1 ); i <= nchi; ++i ) {
		sd += pow( angle_diff( get_symm_corrected_angle( i, res1 ), get_symm_corrected_angle( i, res2 ) ), 2 );
	}

	return sqrt( sd/nchi );
}


// super simple class to grab and print stuff
class PeptoidRotamerRecoverer : public protocols::moves::Mover {
public:
// ctor
PeptoidRotamerRecoverer( bool cyclic );

//dtor
virtual ~PeptoidRotamerRecoverer(){}

// mover interface
virtual void apply( core::pose::Pose & pose );
virtual std::string get_name() const { return "PeptoidRotamerRecoverer"; }
virtual protocols::moves::MoverOP clone() const { return new PeptoidRotamerRecoverer( *this ); }
virtual protocols::moves::MoverOP fresh_instance() const { return clone(); }

private:
bool cyclic_;

};

PeptoidRotamerRecoverer::PeptoidRotamerRecoverer( bool cyclic ) :
  cyclic_( cyclic )
{}

void
PeptoidRotamerRecoverer::apply( core::pose::Pose & pose )
{
	using namespace core;
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
  using namespace basic::options;
  using namespace basic::options::OptionKeys::packing;

	// first copy the original pose
	Pose orig_pose( pose );

	// setup task factory
	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory );
	tf->push_back( new core::pack::task::operation::RestrictResidueToRepacking() );
	tf->push_back( new core::pack::task::operation::InitializeFromCommandline );
	if ( option[ OptionKeys::packing::resfile ].user() ) {
		tf->push_back( new core::pack::task::operation::ReadResfile );
	}


	// setup packer task
	pack::task::PackerTaskOP pt( tf->create_task_and_apply_taskoperations( pose ) );
	pt->set_bump_check( false ); // no bump check, use all possible rotamers

	// setup score function
	core::scoring::ScoreFunctionOP scrfxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::MM_STD_WTS ) );

	// setup pack rotamers mover
	protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover( scrfxn, pt ) );

	// repack the pose
	pack_mover->apply( pose );

	// compare orig and packed poses
	Size x1_total( 0 );
	Size x1_correct( 0 );
	Size x1x2_total( 0 );
	Size x1x2_correct( 0 );
	Real const max_diff( 30 );

	for( Size i( 1 ); i <= pose.total_residue(); ++i ) { // total residue might not be compatible with symetric poses

		Size chain_id( orig_pose.residue( i ).chain() );
		bool chain_match( false );
		if ( option[chains_to_compare].user() ) {
			for ( Size i(1); i <= option[chains_to_compare].value().size(); ++i ) {
				if ( chain_id == core::Size( option[chains_to_compare].value()[ i ] ) ) chain_match = true;
			}
		} else {
			chain_match = true;
		}

		if ( chain_match ) {

			conformation::Residue res_orig( orig_pose.residue( i ) ) ;
			conformation::Residue res_pack( pose.residue( i ) ) ;

			if( res_orig.type().nchi() >= 1 ) {
				x1_total++;
				Real x1_diff( angle_diff( get_symm_corrected_angle( 1, res_orig ), get_symm_corrected_angle( 1, res_pack ) ) );
				if( x1_diff <= max_diff ) x1_correct++;
			}

			if( res_orig.type().nchi() >= 2 ) {
				x1x2_total++;
				Real x1_diff( angle_diff( get_symm_corrected_angle( 1, res_orig ), get_symm_corrected_angle( 1, res_pack ) ) );
				Real x2_diff( angle_diff( get_symm_corrected_angle( 2, res_orig ), get_symm_corrected_angle( 2, res_pack ) ) );
				if( ( x1_diff <= max_diff ) && (x2_diff <= max_diff ) ) x1x2_correct++;
			}
		}
	}

	// report results
	TR << std::fixed << std::setprecision(3) << "{ 'pdb_name': \"" << pose.pdb_info()->name() << "\", 'x1 correct/total': \"" << x1_correct << "/" <<  x1_total << "\", 'x1x2 correct/total': \"" << x1x2_correct << "/" <<  x1x2_total << std::endl;
}

// typedefs
typedef utility::pointer::owning_ptr< PeptoidRotamerRecoverer > PeptoidRotamerRecovererOP;

int
main( int argc, char * argv [] )
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys::cyclization;
  using namespace protocols::simple_moves;
  using namespace protocols::moves;

  // add local options
 	option.add( cyclic, "cyclic" ).def("False");
	option.add( chains_to_compare, "chains_to_compare");

  // init
  devel::init( argc, argv );

	// setup sequence mover
	SequenceMoverOP sm( new SequenceMover() );

  // setup the cyclization mover(s) ( just add patches and constraints don't minimize )
	if ( option[chains_to_cyclize].user() && option[cyclic].value() == true ) {
		core::Size num_cyclic_chains( option[chains_to_cyclize].value().size() );
		for ( core::Size i(1); i <= num_cyclic_chains; ++i ) {
			sm->add_mover( new CyclizationMover( option[chains_to_cyclize].value()[i], true, false, 0 ) );
		}
	}

	// setup peptoid dihedral grabber mover
	PeptoidRotamerRecovererOP prr( new PeptoidRotamerRecoverer( option[cyclic].value() ) );
	sm->add_mover( prr );

  // go go go
  protocols::jd2::JobDistributor::get_instance()->go( sm );

  TR << "\n+-----------------------------------------------------------------+\n"
     <<   "|                              DONE                               |\n"
     <<   "+-----------------------------------------------------------------+" << std::endl;

  return 0;
}
