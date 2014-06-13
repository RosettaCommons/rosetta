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
basic::options::BooleanOptionKey const qms( "qms" );
basic::options::BooleanOptionKey const kmc( "kmc" );

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
	if ( res_type_name3 == "601" && chi_num == 2 && temp_chi >= 135 && temp_chi <= 315 ) {
		return temp_chi - 180;
	} else if ( res_type_name3 == "602" && chi_num == 2 && temp_chi >= 135 && temp_chi <= 315 ) {
		return temp_chi - 180;
	}	else if ( res_type_name3 == "101" && chi_num == 2 && temp_chi >= 135 && temp_chi <= 315 ) {
		return temp_chi - 180;
	}	else if ( res_type_name3 == "401" && chi_num == 3 && temp_chi >= 135 && temp_chi <= 315 ) {
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


// some ugly functions to get rotamer sample data from rotlibs
std::map< std::string, core::pack::dunbrack::SingleResiduePeptoidLibraryCOP > peptoid_rotlibs_;
core::pack::dunbrack::SingleResiduePeptoidLibraryCAP
get_rsrpl(core::chemical::ResidueType const & rsd_type )
{
	using namespace core;
	using namespace pack;
	using namespace dunbrack;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

		// get some info about amino acid type
	std::string aa_name3( rsd_type.name3() );
	Size n_rotlib_chi( rsd_type.nchi() - rsd_type.n_proton_chi() );

	if ( peptoid_rotlibs_.find( aa_name3 ) == peptoid_rotlibs_.end() ) {

		// create izstream from path
		std::string dir_name;
		if ( option[ qms ].value() == true ) {
			dir_name = basic::database::full_name( "/rotamer/peptoid_rotlibs_qms/" );
		} else {
			dir_name = basic::database::full_name( "/rotamer/peptoid_rotlibs_kmc/" );
		}
		std::string file_name = rsd_type.get_peptoid_rotlib_path();
		utility::io::izstream rotlib_in( dir_name + file_name );
		std::cout << "Reading in rot lib " << dir_name + file_name << "...";

		// get an instance of RotamericSingleResiduePeptoidLibrary, but need a RotamerLibrary to do it
		// this means that when ever you read in the PEPTOID libraries you will also read in the Dunbrack libraries
		// this may need to be a pointer to the full type and not just a SRRLOP

		// this comes almost directally from RotmerLibrary.cc::create_rotameric_dunlib()
		SingleResiduePeptoidLibraryOP peptoid_rotlib;

		switch ( n_rotlib_chi ) {
		case 1: {
			RotamericSingleResiduePeptoidLibrary< ONE > * r1 =
				new RotamericSingleResiduePeptoidLibrary< ONE >();
			r1->set_n_chi_bins( rsd_type.get_peptoid_rotlib_n_bin_per_rot() );
			r1->read_from_file( rotlib_in );
			peptoid_rotlib = r1;
			break;
		}
		case 2: {
			RotamericSingleResiduePeptoidLibrary< TWO > * r2 =
				new RotamericSingleResiduePeptoidLibrary< TWO >();
			r2->set_n_chi_bins( rsd_type.get_peptoid_rotlib_n_bin_per_rot() );
			r2->read_from_file( rotlib_in );
			peptoid_rotlib = r2;
			break;
		}
		case 3: {
			RotamericSingleResiduePeptoidLibrary< THREE > * r3 =
				new RotamericSingleResiduePeptoidLibrary< THREE >();
			r3->set_n_chi_bins( rsd_type.get_peptoid_rotlib_n_bin_per_rot() );
			r3->read_from_file( rotlib_in );
			peptoid_rotlib = r3;
			break;
		}
		case 4: {
			RotamericSingleResiduePeptoidLibrary< FOUR > * r4 =
				new RotamericSingleResiduePeptoidLibrary< FOUR >();
			r4->set_n_chi_bins( rsd_type.get_peptoid_rotlib_n_bin_per_rot() );
			r4->read_from_file( rotlib_in );
			peptoid_rotlib = r4;
			break;
		}
		default:
			utility_exit_with_message( "ERROR: too many chi angles desired for peptoid library: " + n_rotlib_chi );
			break;
		}

		// add new rotamer library to map
		peptoid_rotlibs_[ aa_name3 ] = peptoid_rotlib;
		std::cout << "done!" << std::endl;
	}
	return ( peptoid_rotlibs_.find( aa_name3 )->second)();
}

/// @details Returns the preceeding omg.
core::Real
get_omg_from_rsd(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose
)
{
	using namespace core;
	assert( rsd.is_peptoid() );

	if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && pose.residue( pose.conformation().chain_end( rsd.chain() ) ).has_variant_type( chemical::CTERM_CONNECT ) ) {
		assert( pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_protein() || pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_peptoid() );
		return pose.residue( pose.conformation().chain_end( rsd.chain() ) ).mainchain_torsion( 3 );

	} else if ( rsd.is_lower_terminus() ) {
		return 180.0;

	}	else {
		assert( pose.residue( rsd.seqpos() - 1 ).is_protein() || pose.residue( rsd.seqpos() - 1 ).is_peptoid() );
		return pose.residue( rsd.seqpos() - 1 ).mainchain_torsion( 3 );
	}

}

/// @details Handle lower-term residues by returning a "neutral" phi value
core::Real
get_phi_from_rsd(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose
)
{
	using namespace core;
	assert( rsd.is_peptoid() );

	if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && pose.residue( pose.conformation().chain_end( rsd.chain() ) ).has_variant_type( chemical::CTERM_CONNECT ) ) {
		assert( pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_protein() || pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_peptoid() );
		return rsd.mainchain_torsion( 1 );

	} else if ( rsd.has_variant_type( chemical::ACETYLATED_NTERMINUS ) ) {
		return rsd.mainchain_torsion( 1 );

	} else if ( rsd.is_lower_terminus() ) {
		return -90.0;

	}	else {
		return rsd.mainchain_torsion( 1 );
	}

}

core::Real
get_psi_from_rsd(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose
)
{
	using namespace core;
	assert( rsd.is_peptoid() );

	if ( rsd.has_variant_type( chemical::CTERM_CONNECT ) && pose.residue( pose.conformation().chain_begin( rsd.chain() ) ).has_variant_type( chemical::NTERM_CONNECT ) ) {
		assert( pose.residue( pose.conformation().chain_begin( rsd.chain() ) ).is_protein() || pose.residue( pose.conformation().chain_begin( rsd.chain() ) ).is_peptoid() );
		return rsd.mainchain_torsion( 2 );

	} else if ( rsd.has_variant_type( chemical::METHYLATED_CTERMINUS ) ) {
		return rsd.mainchain_torsion( 2 );

	} else if ( rsd.is_upper_terminus() ) {
		return 180.0;

	}	else {
		return rsd.mainchain_torsion( 2 );
	}

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
		Real omg, phi, psi;

		// first residue
		if ( i == 1 ) {
			if ( cyclic_ ) {
			  omg = numeric::principal_angle_degrees( pose.omega( pose.total_residue() ) );
			} else {
				omg = numeric::principal_angle_degrees( 0.0 );
			}
		} else {
		  omg = numeric::principal_angle_degrees( pose.omega( i - 1 ) );
		}

		phi = numeric::principal_angle_degrees( pose.phi( i ) );
		psi = numeric::principal_angle_degrees( pose.psi( i ) );

		TR_EX << "'omg': " << std::setw( 9 ) << numeric::principal_angle_degrees( omg )
		<< ", 'phi': " << std::setw( 9 ) << numeric::principal_angle_degrees( phi )
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

		if( !concrete_residue->get_peptoid_rotlib_path().empty() ) {

			// get srrl
			core::pack::dunbrack::SingleResiduePeptoidLibraryCAP peptoid_rl;
			peptoid_rl = get_rsrpl( pose.residue( i ).type() );
			//peptoid_rl =  rl.get_peptoid_rotamer_library( pose.residue( i ).type() );

			// fill rotamer vector
			peptoid_rl->fill_rotamer_vector( pose, *scrfxn, *pt, packer_neighbor_graph, concrete_residue, existing_residue, extra_chi_steps, buried, rotamers );

			Real omg_rl( get_omg_from_rsd( pose.residue( i ), pose ) );
			Real phi_rl( get_phi_from_rsd( pose.residue( i ), pose ) );
			Real psi_rl( get_psi_from_rsd( pose.residue( i ), pose ) );
			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > drsd( peptoid_rl->get_all_rotamer_samples( omg_rl, phi_rl, psi_rl ) );

			for ( Size j(1); j <= drsd.size(); ++j ) {
				std::stringstream stupid;
				stupid << "\"" << pose.pdb_info()->name() << "\"";
				TR_IR << std::fixed << std::setprecision(3)
				<< "{ 'pdb_name': "	<< std::setw(21) << stupid.str()
				<<  ", 'aa': \"" << pose.residue( i ).type().name3()
				<< "\", 'chain': " << pose.residue( i ).chain()
				<<   ", 'num_chi': " << pose.residue( i ).type().nchi() -  pose.residue( i ).type().n_proton_chi()
				<< ", 'res': " << std::setw( 3 ) << i	<< ", 'rot': " << std::setw(3) << j << ", ";

				TR_IR << "'omg': " << std::setw( 9 ) << numeric::principal_angle_degrees( omg ) << ", "	<< "'phi': " << std::setw( 9 ) << numeric::principal_angle_degrees( phi ) << ", "	<< "'psi': " << std::setw( 9 ) << numeric::principal_angle_degrees( psi ) << ", "
				<< "'x1': "  << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( 1, pose.residue( i ).type().name3(), drsd[j].chi_mean()[1] ) ) << ", "
				<< "'x2': "  << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( 2, pose.residue( i ).type().name3(), drsd[j].chi_mean()[2] ) ) << ", "
				<< "'x3': "  << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( 3, pose.residue( i ).type().name3(), drsd[j].chi_mean()[3] ) ) << ", "
				<< "'x4': "  << std::setw( 9 ) << numeric::principal_angle_degrees( get_symm_corrected_angle( 4, pose.residue( i ).type().name3(), drsd[j].chi_mean()[4] ) ) << ", "
				<< "'sd1': " << std::setw( 9 ) <<  drsd[j].chi_sd()[1] << ", "
				<< "'sd2': " << std::setw( 9 ) <<  drsd[j].chi_sd()[2] << ", "
				<< "'sd3': " << std::setw( 9 ) <<  drsd[j].chi_sd()[3] << ", "
				<< "'sd4': " << std::setw( 9 ) <<  drsd[j].chi_sd()[4] << ", "
				<< "'prob': " << std::setw( 9 ) << drsd[j].probability() << ", "
				<< "'rms_dist': " << std::setw( 9 ) << calc_dist( *rotamers[j],	pose.residue( i ) )
				<< " }," << std::endl;

			}
		}
		/*
		// now print out the RotamerVector
		for ( Size j(1); j <= rotamers.size(); ++j ) {
			TR_IR << std::fixed << std::setprecision(3)	<< "{ 'pdb_name': \""	<< pose.pdb_info()->name()<< "\", 'aa': \"" << pose.residue( i ).type().name3() << "\", 'res': " << std::setw( 3 ) << i	<< ", rot: " << j << " ";

			for( Size k(1); k <= rotamers[j]->type().nchi(); ++k ) {
				std::stringstream chi_string;
				chi_string << "'x" << k << "': ";
				if( k == pose.residue( i ).type().nchi() ) { // if it is the last don't print the ", "
					TR_IR << chi_string.str() << std::setw( 9 ) << numeric::principal_angle_degrees( rotamers[ j ]->chi( k ) );
				} else {
					TR_IR << chi_string.str() << std::setw( 9 ) << numeric::principal_angle_degrees( rotamers[ j ]->chi( k ) ) << ", ";
				}
			}

			// distance from native rotamer
			TR_IR << " rms_dist: " << calc_dist( *rotamers[j],	pose.residue( i ) ) << std::endl;
		}
		*/

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
 	option.add( qms, "qms" ).def("False");
 	option.add( kmc, "kmc" ).def("True");

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


	// derp derp
	//for ( core::Real i( -180.00 ); i <= 360.00; i += 5.00 ) {
	//	core::Real nnpa( numeric::nonnegative_principal_angle_degrees( i ) );
	//	core::Real symc( nnpa );
	//	if ( nnpa >= 135 && nnpa <= 315 ) { symc = symc - 180.00; }
	//	core::Real sypa( numeric::principal_angle_degrees( symc ) );
	//	std::cout << "DEBUG: " << i << "\t" << nnpa << "\t" << symc << "\t" << sypa << std::endl;
	//}

  // go go go
	protocols::jd2::JobDistributor::get_instance()->go( sm );

  TR << "\n+-----------------------------------------------------------------+\n"
     <<   "|                              DONE                               |\n"
     <<   "+-----------------------------------------------------------------+" << std::endl;

  return 0;
}
