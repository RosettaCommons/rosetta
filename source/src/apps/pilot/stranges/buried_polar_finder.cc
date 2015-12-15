// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Ben Stranges (stranges@unc.edu)
/// @brief  runs calculators for some random optimization purposes

// Unit headers
#include <devel/init.hh>


//project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/moves/Mover.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
// //vishal's stuff
// #include <core/scoring/cont_sasa.hh>

// Andrew's stuff
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

// Utility Headers
#include <basic/MetricValue.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/exit.hh>


// jd2 headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <basic/Tracer.hh>

// option key includes
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.stranges_test.sasa_calcs_test" );
static THREAD_LOCAL basic::Tracer TR2( "CalcOutput" );
static THREAD_LOCAL basic::Tracer TRdebug( "sasa_calc_debug" );


using namespace core;
using namespace utility;
using namespace protocols;
using namespace protocols::moves;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::jd2;
using namespace core::pose::metrics;
using namespace protocols::toolbox::pose_metric_calculators;


basic::options::BooleanOptionKey const print_only_unsat( "print_only_unsat" );

OPT_1GRP_KEY( Boolean, bpf, variable_sasa_radii )


// mover deffinition
class CalcsTestMover : public Mover {
	public:

	CalcsTestMover();

	virtual void apply( core::pose::Pose& pose );

	virtual std::string get_name() const{
		return "CalcsTestMover";
}

virtual void register_calculators( );

virtual void calc_stuff(
	pose::Pose & pose,
	core::id::AtomID_Map< core::Real > & atom_sasa,
	core::id::AtomID_Map< core::Size > & atom_hbonds
);

virtual void pretty_print(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Size resnum,
	Size atom_number,
	Real sasa,
	Size numHbonds,
	bool unsat
) const;

virtual core::Size satisfaction_cutoff( std::string atom_type );

virtual MoverOP clone() const {
	return MoverOP( new CalcsTestMover( *this ) );
}

virtual MoverOP fresh_instance() const {
	return clone();
}

private:
core::scoring::ScoreFunctionOP scorefxn_;

utility::file::FileName posename_;

//names of calculators;
std::string Sasa_, NumberHBonds_, BuriedUnsatisfiedPolars_;
bool calcs_ready_;
//probe radius
core::Real probe_radius_;
//numbers of hbonds sat/unsat
Real n_burried_unsat_, n_exposed_sat_;
//cutoff to be considered burried...
core::Real burial_cutoff_;
//print only burried unsat atoms
bool print_unsat_only_;
//number of atoms of different types in a protein
Real n_heavy_atoms_, n_burried_N_, n_burried_O, n_burried_polars_,n_exposed_polars_, n_polars_;
};

CalcsTestMover::CalcsTestMover() {
	// variable definitions
	scorefxn_ = core::scoring::get_score_function();
	calcs_ready_=false;

	burial_cutoff_ = basic::options::option[basic::options::OptionKeys::pose_metrics::atomic_burial_cutoff];
	print_unsat_only_ = option[ print_only_unsat ];
	probe_radius_ = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius];
}

//APPLY
void CalcsTestMover::apply( core::pose::Pose & pose ){
	using namespace std;
	posename_ = pose.pdb_info()->name();

	if ( !calcs_ready_ ) {
		register_calculators();
		calcs_ready_ = true;
	}

	//score for safety
	(*scorefxn_)(pose);

	//metric values for things
	basic::MetricValue< core::id::AtomID_Map< core::Real > > mv_atom_sasa;
	basic::MetricValue< core::id::AtomID_Map< core::Size > > mv_atom_hbonds;
	basic::MetricValue< core::id::AtomID_Map< bool > > mv_unsat_map;

	//run calculators
	pose.metric(Sasa_, "atom_sasa", mv_atom_sasa);
	pose.metric(NumberHBonds_, "atom_Hbonds", mv_atom_hbonds);
	pose.metric(BuriedUnsatisfiedPolars_, "atom_bur_unsat", mv_unsat_map);

	//set values
	core::id::AtomID_Map< core::Real > atom_sasa ( mv_atom_sasa.value() );
	core::id::AtomID_Map< core::Size > atom_hbonds ( mv_atom_hbonds.value() );
	core::id::AtomID_Map< bool > unsat_map (mv_unsat_map.value() ) ;

	//continuous_sasa calculation
	//total_sasa_ = core::scoring::calc_per_atom_cont_sasa( this_pose, atom_sasa_, residue_sasa_, probe_radius_ );


	//figure out unsat stuff calculated values
	calc_stuff( pose, atom_sasa, atom_hbonds);

	//figure out some percentages
	//core::Real burried_fraction_unsat ( n_burried_unsat_ / n_burried_polars_);  // unused ~Labonte
	//core::Real total_fraction_unsat ( n_burried_unsat_ / n_polars_ );  // unused ~Labonte
	//core::Real exposed_fraction_sat ( n_exposed_sat_ / n_exposed_polars_);  // unused ~Labonte
	//core::Real total_exposed_sat ( n_exposed_sat_ / n_polars_);  // unused ~Labonte
	//  //debugging
	//  cout << "Total burried:  "<<  n_burried_unsat_ << endl;
	//  cout << "Burried Polars:  " << n_burried_polars_ <<  " fraction unsat:  "<< burried_fraction_unsat << endl;
	//  cout << "Total Polars:  " << n_polars_ <<" n_polars_ << total faction unsat: " << total_fraction_unsat << endl;

	//  //print totals
	//  TR << "INPUT:  " << setw(25)  << posename_.base()
	//    << "   total_burried_unsat: "<< setprecision(0) << fixed << n_burried_unsat_
	//    << "   faction_burried_polars_unsat: "<< setprecision(3) << fixed << burried_fraction_unsat
	//    << "   faction_total_polars_unsat: "<< setprecision(3) << fixed << total_fraction_unsat
	//    << endl;

	//  TR << "INPUT:  " << setw(25)  << posename_.base()
	//    << "   total_exposed_sat: "<< setprecision(0) << fixed << n_exposed_sat_
	//    << "   faction_exposed_polars_sat: "<< setprecision(3) << fixed << exposed_fraction_sat
	//    << "   faction_total_exposed_sat: "<< setprecision(3) << fixed << total_exposed_sat
	//    << endl;

	//print totals
	TR << "INPUT:  " << setw(25)  << posename_.base()
		<< "   burried_unsat: "<< setprecision(0) << fixed << n_burried_unsat_
		<< "   burried_polar: "<< setprecision(0) << fixed << n_burried_polars_
		<< "   exposed_sat: "<< setprecision(0) << fixed << n_exposed_sat_
		<< "   exposed_polar: " << setprecision(0) << fixed << n_exposed_polars_
		<< "   total_polar: " << setprecision(0) << fixed << n_polars_
		<< "   probe_radius:  "<< setprecision(2) << fixed << probe_radius_
		<< endl;

	set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY); //hacky way to prevent output

	return;
}//end apply


//Prints out info about all don/accp atoms
//Borrowed heavily from Florian's BuriedUnsatPolarsCalculator
void CalcsTestMover::calc_stuff(
	pose::Pose & pose,
	core::id::AtomID_Map< core::Real > & atom_sasa,
	core::id::AtomID_Map< core::Size > & atom_hbonds
){
	using namespace protocols::toolbox::pose_metric_calculators;

	//setup initial values
	n_burried_unsat_ = 0;
	n_exposed_sat_ = 0;

	n_heavy_atoms_ = 0;
	n_burried_N_ = 0;
	n_burried_O = 0;
	n_burried_polars_ = 0;
	n_polars_ = 0;
	n_exposed_polars_ = 0;

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue rsd = pose.residue( ii );
		//for all atoms in res
		for ( Size at = 1; at <= rsd.nheavyatoms(); ++at ) {
			core::id::AtomID atid( at, ii );
			bool this_atom_bur_unsat(false);
			bool this_atom_exposed_sat(false);
			++n_heavy_atoms_;

			//is it a don/accp?
			if ( rsd.atom_type( at ).is_acceptor() || rsd.atom_type( at ).is_donor() ) {
				++n_polars_;
				//how many hbonds to atom
				Size nhbonds ( atom_hbonds[ atid ] );
				Size bonded_heavyatoms = rsd.n_bonded_neighbor_all_res( at ) - rsd.type().number_bonded_hydrogens( at );
				Size num_bonds (nhbonds + bonded_heavyatoms);
				//find sasa by adding in H's
				Real cursasa =  atom_sasa[ atid ];
				TRdebug << "Residue: "<< ii << "  Atom: "<< rsd.type().atom_name( at ) << " SASA: " << atom_sasa[ atid ]
					<< "\n";
				for ( Size hcount = rsd.type().attached_H_begin( at ); hcount<= rsd.type().attached_H_end( at ); hcount++ ) {
					cursasa = cursasa + atom_sasa[ core::id::AtomID ( hcount, ii ) ];
					TRdebug << "Hydrogen: " << rsd.type().atom_name( hcount ) << "  SASA: "
						<<  atom_sasa[ core::id::AtomID ( hcount, ii ) ] << "\n";
				}
				TRdebug<< std::endl;
				//how many does this require?
				Size satisfaction_cut = CalcsTestMover::satisfaction_cutoff( rsd.type().atom_type( at ).name() );
				//What is this atom's situation?
				if ( cursasa < burial_cutoff_ ) {
					++n_burried_polars_;
					if ( num_bonds < satisfaction_cut ) {
						this_atom_bur_unsat = true;
					}
				}
				if ( cursasa > burial_cutoff_ ) {
					++n_exposed_polars_;
					if ( num_bonds >= satisfaction_cut )  {
						this_atom_exposed_sat = true;
					}
				}
				//now print info about all polar atoms
				if ( this_atom_bur_unsat ) {
					CalcsTestMover::pretty_print( pose, rsd, ii, at, cursasa, nhbonds, this_atom_bur_unsat );
				} else if ( !print_unsat_only_ ) {
					CalcsTestMover::pretty_print( pose, rsd, ii, at, cursasa, nhbonds, this_atom_bur_unsat );
				}

			}//if is don/accp
			if ( this_atom_bur_unsat ) {
				++n_burried_unsat_;
			}
			if ( this_atom_exposed_sat ) {
				++n_exposed_sat_;
			}
		}//for atoms in a residue
	}//for all residues
	return;
}//end calc_stuff

//print out a prettier list of the outputs
void CalcsTestMover::pretty_print(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Size resnum,
	Size atom_number,
	Real sasa,
	Size numHbonds,
	bool unsat
) const
{
	using namespace std;
	utility::file::FileName filename( pose.pdb_info()->name() );

	TR2 << "STRCTURE:  " << setw(25) << filename.base() << "   "
		<< "CHAIN:  " << setw(3) << pose.pdb_info()->chain(resnum) << "   "
		<< "RESIDUE:  " << setw(4) << pose.pdb_info()->number(resnum) << "   "
		<< "ATOM:  " << setw(6)  << rsd.type().atom_name( atom_number )  << "   "
		<< "SASA:  " << setw(8)  << setprecision(3) << fixed<< sasa << "   "
		<< "HBONDS:  " << setw(3) << numHbonds << "   "
		<< "BUR_UNSAT:  " << unsat
		<< endl;

	return;
}//end pretty print

//register the calculators used here...
void CalcsTestMover::register_calculators(){
	using namespace core::pose::metrics;
	using namespace protocols::toolbox::pose_metric_calculators;

	//find name of job
	// protocols::jd2::JobOP const job_me( JobDistributor::get_instance()->current_job() );
	std::string name( posename_.base() );


	//now check on calculators
	Sasa_ = "Sasa_" + name;
	if ( CalculatorFactory::Instance().check_calculator_exists( Sasa_ ) ) {
		Warning() << "In InterfaceAnalyzerMover, calculator " << Sasa_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		if ( basic::options::option[ basic::options::OptionKeys::bpf::variable_sasa_radii ] ) {
			using namespace devel::vardist_solaccess;
			CalculatorFactory::Instance().register_calculator( Sasa_, VarSolDistSasaCalculatorOP( new VarSolDistSasaCalculator ) );
		} else {
			CalculatorFactory::Instance().register_calculator( Sasa_, PoseMetricCalculatorOP( new simple_calculators::SasaCalculatorLegacy ) );
		}

	}

	NumberHBonds_ = "NumberHBonds_" + name;
	if ( CalculatorFactory::Instance().check_calculator_exists( NumberHBonds_ ) ) {
		Warning() << "In InterfaceAnalyzerMover, calculator " << NumberHBonds_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator( NumberHBonds_, PoseMetricCalculatorOP( new NumberHBondsCalculator ));
	}

	BuriedUnsatisfiedPolars_ = "BuriedUnsatisfiedPolars_" + name;
	if ( CalculatorFactory::Instance().check_calculator_exists( BuriedUnsatisfiedPolars_ ) ) {
		Warning() << "In InterfaceAnalyzerMover, calculator " << BuriedUnsatisfiedPolars_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		CalculatorFactory::Instance().register_calculator(  BuriedUnsatisfiedPolars_, PoseMetricCalculatorOP( new BuriedUnsatisfiedPolarsCalculator(Sasa_, NumberHBonds_) ));
	}

	return;
}


core::Size
CalcsTestMover::satisfaction_cutoff( std::string atom_type )
{

	//according to jk, buried hydroxyls are often seen making only one hydrogen bond. also, ether oxygens often are bad h-bond acceptors
	if ( atom_type == "OH" ) return 2;

	//backbone oxygens also only have one h-bbond in most secondary structure elements
	else if ( atom_type == "OCbb" ) return 2;

	else if ( atom_type ==  "S" ) return 2;

	//everything else we expect to have 3 bonded/h-bonded neighbours to count as satisfied
	else return 3;

}


////////////////////////////////////////////////
// Main
////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


		option.add( print_only_unsat, "Print only unsatisfied atoms" ).def(false);
		NEW_OPT( bpf::variable_sasa_radii, "Use the ", false );

		// init
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new CalcsTestMover ) );

		std::cout << "Done! -------------------------------"<< std::endl;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

} //end main
