// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    interface_statistics.cc
/// @brief   Computing interface statistics for membrane protein interfaces
/// @details last Modified: 5/14/15
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/scoring/Interface.fwd.hh>
#include <protocols/scoring/Interface.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/PiPiCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// Package Headers
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <core/pose/util.hh>
#include <utility/string_util.hh>


static basic::Tracer TR( "apps.pilot.jkleman.interface_statistics" );

using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;
using namespace protocols::moves;
using namespace protocols::scoring;
using namespace protocols::toolbox::pose_metric_calculators;

class MPInterfaceStatistics : public Mover {
	
public:
	
	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Construct a Default Membrane Position Mover
	MPInterfaceStatistics( Size jump );
	
	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	MPInterfaceStatistics( MPInterfaceStatistics const & src );
	
	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	MPInterfaceStatistics &
	operator=( MPInterfaceStatistics const & src );
	
	/// @brief Destructor
	~MPInterfaceStatistics();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	virtual std::string get_name() const;

	/// @brief Get Membrane protein interface statistics
	virtual void apply( Pose & pose );

	/// @brief Register Options with JD2
	void register_options();
	
	/// @brief Initialize Mover options from the comandline
	void init_from_cmd();
	
	/// @brief Get size of the interface
	Real get_size( Pose & pose );
	
	/// @brief Get charge in the interface
	SSize get_charge( Pose & pose );
	
	/// @brief Get number of charged residues in the interface
	Size get_number_charges( Pose & pose );
	
	/// @brief Get hydrophobicity of the interface (UHS, Koehler & Meiler, Proteins, 2008)
	Real get_hydrophobicity( Pose & pose );
	
	/// @brief Get number of contacts
	Size get_number_of_contacts();
	
	/// @brief Get surface complementarity
	Real get_surface_complementarity( Pose & pose );

	/// @brief Get size of residues in the interface
	Real get_residue_size( Pose & pose );
	
	/// @brief Get number of Hbonds across interface
	Size get_number_hbonds( Pose & pose );
	
	/// @brief Get number of pi-stacking across interface
	Size get_number_pi_stacking( Pose & pose );
	
	/// @brief Compute a ruggedness factor of the interface???
	Real get_ruggedness_factor( Pose & pose );
	
	/// @brief Compute an average deviation from a plane along the interface
	Real get_plane_deviation( Pose & pose );
	
	/// @brief Get contact preferences for AAs: Yan & Jernigan, 2008, ProteinJ
	
	/// @brief Get in/out orientation
	Real get_inout_orientation( Pose & pose );

	/*
	size
	charge
	hydrophobicity
	size of residues
	surface complementarity
	number of contacts
	Hbonds across interface
	pi-stacking across interface?
	ruggedness factor
	deviation from a plane
	contact preferences of AAs: Yan, 2008, ProteinJ
	NEED TO SPLIT INTERFACE INTO TM AND SOL
	*/

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

//	/// @brief Create a Clone of this mover
//	virtual protocols::moves::MoverOP clone() const;
//
//	/// @brief Create a Fresh Instance of this Mover
//	virtual protocols::moves::MoverOP fresh_instance() const;
//
//	/// @brief Pase Rosetta Scripts Options for this Mover
//	void parse_my_tag(
//					  utility::tag::TagCOP tag,
//					  basic::datacache::DataMap &,
//					  protocols::filters::Filters_map const &,
//					  protocols::moves::Movers_map const &,
//					  core::pose::Pose const &
//					  );
	
private:
	
	/// @brief jump
	Size jump_;
	
	/// @brief Interface
	Interface interface_;
	
	/// @brief Partners;
	std::string partners_;
	
};

//////////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
MPInterfaceStatistics::MPInterfaceStatistics( Size jump ) :
	Mover(),
	jump_( jump ),
	interface_(),
	partners_()
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
MPInterfaceStatistics::MPInterfaceStatistics( MPInterfaceStatistics const & src ) :
	Mover( src ),
	jump_( src.jump_ ),
	interface_( src.interface_ ),
	partners_( src.partners_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
MPInterfaceStatistics &
MPInterfaceStatistics::operator=( MPInterfaceStatistics const & src )
{
	// Abort self-assignment.
	if (this == &src) {
	return *this;
	}

	// Otherwise, create a new object
	return *( new MPInterfaceStatistics( *this ) );
}

/// @brief Destructor
MPInterfaceStatistics::~MPInterfaceStatistics() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
MPInterfaceStatistics::get_name() const {
	return "MPInterfaceStatistics";
}

/// @brief Get Membrane protein interface statistics
void MPInterfaceStatistics::apply( Pose & pose ) {
	
	using namespace protocols::membrane;
	
	TR << "Calling MPInterfaceStatistics" << std::endl;

	register_options();
	init_from_cmd();

	// add membrane
	AddMembraneMoverOP addmem( new AddMembraneMover() );
	addmem->apply( pose );

	// Check the pose is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}

	// create scorefunction and score the pose
	core::scoring::ScoreFunctionOP highres_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
	( *highres_sfxn )( pose );

	// calculate interface
	interface_ = Interface( jump_ );
	interface_.calculate( pose );

	// get size of the interface
	Real size = get_size( pose );
	TR << "sasa " << size << std::endl;
	
	// Get charge in the interface
	SSize charge = get_charge( pose );
	TR << "charge " << charge << std::endl;
	
	// Get number of charged residues in the interface
	Size ncharges = get_number_charges( pose );
	TR << "ncharges " << ncharges << std::endl;
	
	// Get hydrophobicity of the interface
	Real hydro = get_hydrophobicity( pose );
	TR << "hydro " << hydro << std::endl;
	
	// Get number of contacts
	Size ncontacts = get_number_of_contacts();
	TR << "ncontacts " << ncontacts << std::endl;
	
	// Get surface complementarity
	Real shape = get_surface_complementarity( pose );
	TR << "shape_compl " << shape << std::endl;
	
	// Get size of residues in the interface
	Real res_size = get_residue_size( pose );
	TR << "avg_res_size " << res_size << std::endl;
	
	// Get number of Hbonds across interface
	Size nhbonds = get_number_hbonds( pose );
	TR << "nhbonds " << nhbonds << std::endl;
	
	// Get number of pi-stacking across interface
//	Size npi = get_number_pi_stacking( pose );
//	TR << "pi stacking " << npi << std::endl;
	
//	// Compute a ruggedness factor of the interface???
//	Real rugged = get_ruggedness_factor( pose );
//	TR << "rugged " << rugged << std::endl;
//	
//	// Compute an average deviation from a plane along the interface
//	Real plane = get_plane_deviation( pose );
//	TR << "plane dev " << plane << std::endl;
	
	// Get contact preferences for AAs: Yan & Jernigan, 2008, ProteinJ
	
	// Get in/out orientation
	Real inout = get_inout_orientation( pose );
	TR << "inout " << inout << std::endl;
	
} // apply

/// @brief Register Options with JD2
void MPInterfaceStatistics::register_options() {
	
	using namespace basic::options;
	option.add_relevant( OptionKeys::docking::partners );

} // register options

/// @brief Initialize Mover options from the comandline
void MPInterfaceStatistics::init_from_cmd() {
	
	using namespace basic::options;
	
	// docking partners
	if( option[ OptionKeys::docking::partners ].user() ) {
		partners_ = option[ OptionKeys::docking::partners ]();
	}
} // init from commandline

/////////////////////////////////
///// Rosetta Scripts Methods ///
/////////////////////////////////
//
///// @brief Create a Clone of this mover
//protocols::moves::MoverOP
//MPInterfaceStatistics::clone() const {
//	return ( protocols::moves::MoverOP( new MPInterfaceStatistics( *this ) ) );
//}
//
///// @brief Create a Fresh Instance of this Mover
//protocols::moves::MoverOP
//MPInterfaceStatistics::fresh_instance() const {
//	return protocols::moves::MoverOP( new MPInterfaceStatistics() );
//}
//
///// @brief Pase Rosetta Scripts Options for this Mover
//void
//MPInterfaceStatistics::parse_my_tag(
//									   utility::tag::TagCOP tag,
//									   basic::datacache::DataMap &,
//									   protocols::filters::Filters_map const &,
//									   protocols::moves::Movers_map const &,
//									   core::pose::Pose const &
//									   ) {
//		// TODO
//
//}
//
///// @brief Create a new copy of this mover
//protocols::moves::MoverOP
//MPInterfaceStatisticsCreator::create_mover() const {
//	return protocols::moves::MoverOP( new MPInterfaceStatistics );
//}
//
///// @brief Return the Name of this mover (as seen by Rscripts)
//std::string
//MPInterfaceStatisticsCreator::keyname() const {
//	return MPInterfaceStatisticsCreator::mover_name();
//}
//
///// @brief Mover name for Rosetta Scripts
//std::string
//MPInterfaceStatisticsCreator::mover_name() {
//	return "MPInterfaceStatistics";
//}

//////////////////////////////////////////////////////////////////////

typedef utility::pointer::shared_ptr< MPInterfaceStatistics > MPInterfaceStatisticsOP;

//////////////////////////////////////////////////////////////////////

/// @brief Get size of the interface
Real MPInterfaceStatistics::get_size( Pose & pose ) {

	using namespace numeric;
	using namespace core::pose;
	using namespace core::scoring::sasa;
	using namespace protocols::scoring;
	
	// get per-residue SASA
	SasaCalc calc = SasaCalc();
	Real total_sasa( calc.calculate( pose ) );
	
	// partition pose by jump
	Pose partner1, partner2;
	partition_pose_by_jump( pose, jump_, partner1, partner2 );

	// calculate SASAs of individual partners
	Real p1_sasa( calc.calculate( partner1 ) );
	Real p2_sasa( calc.calculate( partner2 ) );
	
	// calculate interface
	Real intf_sasa = - ( total_sasa - p1_sasa - p2_sasa) / 2;

	return intf_sasa;

}// get_size

	//////////////////////////////////////////////////////////////////////

/// @brief Get charge in the interface
SSize MPInterfaceStatistics::get_charge( Pose & pose ) {

	SSize charge( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
	
		std::string res = pose.residue( i ).name3();
		
		if ( interface_.is_interface( i ) && ( res == "ARG" || res == "LYS" ) &&
			 ( pose.residue( i ).atom( "CA" ).xyz().z() >= -15.0 &&
			   pose.residue( i ).atom( "CA" ).xyz().z() <= 15.0 ) ) {
			++charge;
		}
		if ( interface_.is_interface( i ) && ( res == "ASP" || res == "GLU" ) &&
			 ( pose.residue( i ).atom( "CA" ).xyz().z() >= -15.0 &&
			   pose.residue( i ).atom( "CA" ).xyz().z() <= 15.0 ) ) {
			--charge;
		}
	}
	return charge;

} // get charge

	//////////////////////////////////////////////////////////////////////

/// @brief Get number of charged residues in the interface
Size MPInterfaceStatistics::get_number_charges( Pose & pose ) {

	Size charge( 0 );
	
	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		
		std::string res = pose.residue( i ).name3();
		
		if ( interface_.is_interface( i ) &&
			( res == "ARG" || res == "LYS" || res == "ASP" || res == "GLU" ) &&
			( pose.residue( i ).atom( "CA" ).xyz().z() >= -15.0 &&
			 pose.residue( i ).atom( "CA" ).xyz().z() <= 15.0 ) ) {
			++charge;
			TR << "charge: " << res << i << std::endl;
		}
	}
	return charge;

} // get number of charges

	//////////////////////////////////////////////////////////////////////

/// @brief Get hydrophobicity of the interface (UHS, Koehler & Meiler, Proteins, 2008)
Real MPInterfaceStatistics::get_hydrophobicity( Pose & pose ) {

	// hydrophobicity scale
	std::map< char, Real > hydro;
	
	// polar
	hydro[ 'C' ] = 0.01;
	hydro[ 'N' ] = 0.50;
	hydro[ 'Q' ] = 0.46;
	hydro[ 'S' ] = 0.06;
	hydro[ 'T' ] = -0.01;
	
	// charged
	hydro[ 'D' ] = 0.73;
	hydro[ 'E' ] = 0.70;
	hydro[ 'K' ] = 0.90;
	hydro[ 'R' ] = 0.55;
	
	// hydrophobic
	hydro[ 'A' ] = -0.16;
	hydro[ 'G' ] = -0.20;
	hydro[ 'I' ] = -0.39;
	hydro[ 'L' ] = -0.30;
	hydro[ 'M' ] = -0.20;
	hydro[ 'P' ] = 0.50;
	hydro[ 'V' ] = -0.25;
	
	// aromatic
	hydro[ 'F' ] = -0.46;
	hydro[ 'H' ] = 0.38;
	hydro[ 'W' ] = -0.03;
	hydro[ 'Y' ] = -0.12;

	Real total_hydro( 0 );
	
	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		
		if ( interface_.is_interface( i ) ) {
			total_hydro += hydro[ pose.residue( i ).name1() ];
		}
	}
	return total_hydro;

} // get hydrophobicity

	//////////////////////////////////////////////////////////////////////

/// @brief Get number of contacts
Size MPInterfaceStatistics::get_number_of_contacts() {

	return interface_.interface_nres();
	
} // get number of contacts in interface

	//////////////////////////////////////////////////////////////////////

/// @brief Get surface complementarity
Real MPInterfaceStatistics::get_surface_complementarity( Pose & pose ) {

	using namespace core::scoring::sc;
	ShapeComplementarityCalculator shape;
	shape.Init();
	
	return shape.CalcSc( pose, jump_, 1 );

} // get surface complementarity of interface

	//////////////////////////////////////////////////////////////////////

/// @brief Get avg residue in the interface
Real MPInterfaceStatistics::get_residue_size( Pose & pose ) {

	Real avg_size( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		
		if ( interface_.is_interface( i ) ) {
			avg_size += static_cast< Real >( pose.residue( i ).natoms() );
		}
	}
	avg_size /= interface_.interface_nres();
	
	return avg_size;
	
} // get residue size

	//////////////////////////////////////////////////////////////////////

/// @brief Get number of Hbonds across interface
Size MPInterfaceStatistics::get_number_hbonds( Pose & pose ) {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;
	
	// get the partners
	utility::vector1< std::string > partners( utility::string_split( partners_, '_' ) );
	utility::vector1< std::string > partner1( utility::split( partners[1] ) );
	utility::vector1< std::string > partner2( utility::split( partners[2] ) );

	// get the Hbonds
	hbonds::HBondDatabaseCOP hb_database(hbonds::HBondDatabase::get_database());
	hbonds::HBondSet hb_set;
	pose.update_residue_neighbors();
	hbonds::fill_hbond_set( pose, false, hb_set, false );
//	hb_set.show(pose, true, std::cout);
	
	Size nhbonds( 0 );
	
	// go through Hbonds
	for ( Size i = 1; i <= hb_set.nhbonds(); ++i ) {

		HBond hbond = hb_set.hbond( i );
		Size don = hbond.don_res();
		Size acc = hbond.acc_res();

		// check whether donor/acc is in partner1
		// doesn't check whether the residues are defined as interface (they should anyway)
		
		bool don_in_1( false );
		bool acc_in_1( false );
		
		// donor in partner1?
		for ( Size j = 1; j <= partner1.size(); ++j ) {
			don_in_1 = res_in_chain( pose, don, partner1[ j ] );
		}

		// acceptor in partner1?
		for ( Size j = 1; j <= partner1.size(); ++j ) {
			acc_in_1 = res_in_chain( pose, acc, partner1[ j ] );
		}
		
//		TR << "hbond: don " << don << ", acc " << acc;
//		TR << ", don: " << don_in_1 << ", acc: " << acc_in_1 << std::endl;

		// if donor in 1 and acceptor not, then we have Hbond across interface
		if ( don_in_1 == true && acc_in_1 == false ) {
			++nhbonds;
		}
		// similarly, if donor not in 1 but acc is, we also have Hbond across interface
		else if ( don_in_1 == false && acc_in_1 == true ) {
			++nhbonds;
		}
	} // go through Hbonds
	

	return nhbonds;
} // number hbonds across interface

	//////////////////////////////////////////////////////////////////////

/// @brief Get number of pi-stacking across interface
Size MPInterfaceStatistics::get_number_pi_stacking( Pose & pose ) {
	TR << pose.total_residue() << std::endl;

	core::pose::metrics::PoseMetricCalculatorOP pi_pi_calculator ( new protocols::toolbox::pose_metric_calculators::PiPiCalculator() );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "pi_pi", pi_pi_calculator );

	TR << "Pi-stacking " << pose.print_metric( "pi_pi", "pi_pi" ) << std::endl;

	return 0;
}

	//////////////////////////////////////////////////////////////////////

/// @brief Compute a ruggedness factor of the interface???
Real MPInterfaceStatistics::get_ruggedness_factor( Pose & pose ) {
	TR << pose.total_residue() << std::endl;
	return 0.0;
}

	//////////////////////////////////////////////////////////////////////

/// @brief Compute an average deviation from a plane along the interface
Real MPInterfaceStatistics::get_plane_deviation( Pose & pose ) {
	TR << pose.total_residue() << std::endl;
	return 0.0;
}
	//////////////////////////////////////////////////////////////////////


/// @brief Get contact preferences for AAs: Yan & Jernigan, 2008, ProteinJ

	//////////////////////////////////////////////////////////////////////

/// @brief Get in/out orientation
Real MPInterfaceStatistics::get_inout_orientation( Pose & pose ) {
	
	if ( pose.residue( nres_protein( pose ) ).atom("CA").xyz().z() > 0.0 ){
		return 0.0;
	}
	return 1.0;

}// get in/out orientation


////////////////////////////// MAIN ////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::jd2;
		
		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);
		
		MPInterfaceStatisticsOP mpis( new MPInterfaceStatistics( 1 ) );
		JobDistributor::get_instance()->go(mpis);
	}
	catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
