// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilot/doug/PeptoidDihedralGrabber.cc
/// @brief Simply prints backbone and side chain dihedral angles of a peptoid
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// protocols header
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
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
#include <core/pack/dunbrack/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResiduePeptoidLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResiduePeptoidLibrary.tmpl.hh>

// devel headers
#include <devel/init.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/cyclization.OptionKeys.gen.hh>
#include <basic/database/open.hh>

// numeric headers
#include <numeric/angle.functions.hh>

// c++
#include <sstream>
#include <iomanip>
#include <map>

// tracer
static basic::Tracer TR("PeptoidDihedralGrabber");
static basic::Tracer TR_EX("ExperimentalRotamers");
static basic::Tracer TR_IR("InterpolatedRotamers");

// local options
basic::options::BooleanOptionKey const cyclic( "cyclic" );
basic::options::StringOptionKey const tlc( "tlc" );

// a few utility functions
core::Real
get_symm_corrected_angle( core::Size chi_num, std::string tlc, core::Real chi )
{
	using namespace core;
	using namespace chemical;

	std::string res_type_name3( tlc );
	Real temp_chi( numeric::nonnegative_principal_angle_degrees( chi ) );

	//need to add 001

	// if else chain of all symm side chains that we can model in the peptoid databank
	if ( res_type_name3 == "PHE" && chi_num == 2 && temp_chi >= 135 && temp_chi <= 315 ) {
		return temp_chi - 180;
	} else {
		return chi;
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
		sd += pow( angle_diff( get_symm_corrected_angle( i, res1.type().name3(), res1.chi( i ) ), get_symm_corrected_angle( i, res2.type().name3(), res2.chi( i ) ) ), 2 );
	}

	return sqrt( sd/nchi );
}


// super simple class to grab and print stuff
class PeptoidDihedralGrabber : public protocols::moves::Mover {
public:
// ctor
PeptoidDihedralGrabber( bool cyclic );

//dtor
virtual ~PeptoidDihedralGrabber(){}

// mover interface
virtual void apply( core::pose::Pose & pose );
virtual std::string get_name() const { return "PeptoidDihedralGrabber"; }
virtual protocols::moves::MoverOP clone() const { return new PeptoidDihedralGrabber( *this ); }
virtual protocols::moves::MoverOP fresh_instance() const { return clone(); }

private:
bool cyclic_;

};

PeptoidDihedralGrabber::PeptoidDihedralGrabber( bool cyclic ) :
  cyclic_( cyclic )
{}

/*
rot_data_61_cis = [
    { 'aa':  "601", 'omg': -0.206, 'phi': -72.082, 'psi': 162.446, 'x1': 66.409, 'x2': -158.748 },
    { 'aa':  "601", 'omg': -12.364, 'phi': 79.333, 'psi': 176.360, 'x1': 37.311, 'x2': 83.005 }]
 */
void
PeptoidDihedralGrabber::apply( core::pose::Pose & pose )
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace pose;
	using namespace conformation;
	using namespace chemical;

	// get a rotlib
	core::pack::dunbrack::RotamerLibrary const & rl( core::pack::dunbrack::RotamerLibrary::get_instance() );

	// setup patcker task
	core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory;
	tf->push_back( new core::pack::task::operation::RestrictResidueToRepacking() );
	pack::task::PackerTaskOP pt( tf->create_task_and_apply_taskoperations( pose ) );
	pt->set_bump_check( false );

	//TR_IR << *pt << std::endl;

	// setup score function
	core::scoring::ScoreFunctionOP scrfxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::MM_STD_WTS ) );

	// buried?
	bool buried( true );


	// setup neighbor graph
	(*scrfxn)( pose ); // if you don'r score the pose before getting the png you get an assertion failure
	graph::GraphOP packer_neighbor_graph( new graph::Graph( pose.energies().energy_graph() ) );

	for( Size i(1); i <= pose.total_residue(); ++i ) {

		if ( pose.residue( i ).type().name3() == option[ tlc ].value() ) {

			// print out experimental rotamer data
			//TR << "resnum: " << i << " " << pose.residue( i ).type().name()  << " " << pose.residue( i ).has_variant_type( chemical::ACETYLATED_NTERMINUS ) << std::endl;
			// print name
			TR_EX << std::fixed << std::setprecision(3)
			<<   "{ 'pdb_name': \"" << pose.pdb_info()->name()
			<< "\", 'aa': \"" << pose.residue( i ).type().name3()
			<< "\", 'chain': " << pose.residue( i ).chain()
			<<   ", 'num_chi': " << pose.residue( i ).type().nchi() -  pose.residue( i ).type().n_proton_chi()
			<<   ", 'res': " << std::setw( 3 ) << i << ", " ;

			// print preceding omg, phi psi
			Real phi, psi;

			phi = numeric::principal_angle_degrees( pose.phi( i ) );
			psi = numeric::principal_angle_degrees( pose.psi( i ) );

			TR_EX << "'phi': " << std::setw( 9 ) << numeric::principal_angle_degrees( phi )
			<< ", 'psi': " << std::setw( 9 ) << numeric::principal_angle_degrees( psi ) << ", ";

			// print sidechain info
			for( Size j(1); j <= pose.residue( i ).type().nchi(); ++j ) {
				std::stringstream chi_string;
				chi_string << "'x" << j << "': ";
				if( j == pose.residue( i ).type().nchi() ) { // if it is the last don't print the ", "
					TR_EX << chi_string.str() << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( j, pose.residue( i ).type().name3(), pose.residue( i ).chi( j ) ) );
				} else {
					TR_EX << chi_string.str() << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( j, pose.residue( i ).type().name3(), pose.residue( i ).chi( j ) ) ) << ", ";
				}
			}
			TR_EX << " }," << std::endl;


			// print out interpolated rotamers if rotamer lib exists

			// need this for each res
			chemical::ResidueTypeCOP concrete_residue( pose.residue( i ).type() );
			conformation::Residue existing_residue( pose.residue( i ) );
			core::pack::dunbrack::RotamerVector rotamers;
			utility::vector1< utility::vector1< Real > > extra_chi_steps( concrete_residue->nchi() );

			std::cout << "DEBUG BEFORE YUCK" << std::endl;
			// get srdl (yuck)
			core::pack::dunbrack::SingleResidueRotamerLibraryCAP      blah1      = core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( *concrete_residue );
			const core::pack::dunbrack::SingleResidueRotamerLibrary*  blah2      = blah1.get();
			const core::pack::dunbrack::SingleResidueDunbrackLibrary* peptide_rl = static_cast<const core::pack::dunbrack::SingleResidueDunbrackLibrary*>(blah2);
			std::cout << "AFTER BEFORE YUCK" << std::endl;

			// fill rotamer vector
			peptide_rl->fill_rotamer_vector( pose, *scrfxn, *pt, packer_neighbor_graph, concrete_residue, existing_residue, extra_chi_steps, buried, rotamers );
			std::cout << "DEBUG AFTER FRV" << std::endl;

			Real phi_rl( pose.phi(i) );
			Real psi_rl( pose.psi(i) );
			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > drsd( peptide_rl->get_all_rotamer_samples( phi_rl, psi_rl ) );
			std::cout << "DEBUG AFTER RSD" << std::endl;



			for ( Size j(1); j <= drsd.size(); ++j ) {

				core::conformation::Residue resj( *rotamers[j] );
				core::conformation::Residue resi( pose.residue(i) );
				std::stringstream stupid;
				stupid << "\"" << pose.pdb_info()->name() << "\"";
				TR_IR << std::fixed << std::setprecision(3)
				<< "{ 'pdb_name': "	<< std::setw(21) << stupid.str()
				<<  ", 'aa': \"" << pose.residue( i ).type().name3()
				<< "\", 'chain': " << pose.residue( i ).chain()
				<<   ", 'num_chi': " << pose.residue( i ).type().nchi() -  pose.residue( i ).type().n_proton_chi()
				<< ", 'res': " << std::setw( 3 ) << i	<< ", 'rot': " << std::setw(3) << j << ", ";

				TR_IR << "'phi': " << std::setw( 9 ) << numeric::principal_angle_degrees( phi ) << ", "	<< "'psi': " << std::setw( 9 ) << numeric::principal_angle_degrees( psi ) << ", "
				<< "'x1': "  << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( 1, pose.residue( i ).type().name3(), drsd[j].chi_mean()[1] ) ) << ", "
				<< "'x2': "  << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( 2, pose.residue( i ).type().name3(), drsd[j].chi_mean()[2] ) ) << ", "
				<< "'x3': "  << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( 3, pose.residue( i ).type().name3(), drsd[j].chi_mean()[3] ) ) << ", "
				<< "'x4': "  << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( 4, pose.residue( i ).type().name3(), drsd[j].chi_mean()[4] ) ) << ", "
				<< "'sd1': " << std::setw( 9 ) <<  drsd[j].chi_sd()[1] << ", "
				<< "'sd2': " << std::setw( 9 ) <<  drsd[j].chi_sd()[2] << ", "
				<< "'sd3': " << std::setw( 9 ) <<  drsd[j].chi_sd()[3] << ", "
				<< "'sd4': " << std::setw( 9 ) <<  drsd[j].chi_sd()[4] << ", "
				<< "'prob': " << std::setw( 9 ) << drsd[j].probability() << ", "
				<< "'rms_dist': " << std::setw( 9 ) << calc_dist( resj, resi )
				<< " }," << std::endl;

			}
		}


	}

}

// typedefs
typedef utility::pointer::owning_ptr< PeptoidDihedralGrabber > PeptoidDihedralGrabberOP;

int
main( int argc, char * argv [] )
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys::cyclization;
  using namespace protocols::simple_moves;
  using namespace protocols::moves;

  // add local options
 	option.add( cyclic, "cyclic" ).def("False");
 	option.add( tlc, "tlc" ).def("PHE");

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
	PeptoidDihedralGrabberOP pdg( new PeptoidDihedralGrabber( option[cyclic].value() ) );
	sm->add_mover( pdg );

  // go go go
	protocols::jd2::JobDistributor::get_instance()->go( sm );

  TR << "\n+-----------------------------------------------------------------+\n"
     <<   "|                              DONE                               |\n"
     <<   "+-----------------------------------------------------------------+" << std::endl;

  return 0;
}
