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
/// @details last Modified: 8/28/15
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <core/conformation/membrane/MembraneInfo.hh>
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
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/intf.OptionKeys.gen.hh>

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

// External
#include <ObjexxFCL/format.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa/util.hh>
#include <protocols/docking/util.hh>
#include <protocols/membrane/util.hh>


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.jkleman.interface_statistics" );

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
	MPInterfaceStatistics();

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

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	// /// @brief Create a Clone of this mover
	// virtual protocols::moves::MoverOP clone() const;
	//
	// /// @brief Create a Fresh Instance of this Mover
	// virtual protocols::moves::MoverOP fresh_instance() const;
	//
	// /// @brief Pase Rosetta Scripts Options for this Mover
	// void parse_my_tag(
	//       utility::tag::TagCOP tag,
	//       basic::datacache::DataMap &,
	//       protocols::filters::Filters_map const &,
	//       protocols::moves::Movers_map const &,
	//       core::pose::Pose const &
	//       );

private: // methods

	/// @brief Register Options with JD2
	void register_options();

	/// @brief Initialize Mover options from the comandline
	void init_from_cmd();

	/// @brief Get the chains from the commandline
	utility::vector1< bool > get_chains_from_cmd( Pose & pose );

	/// @brief calculate interfaces
	void calculate_interfaces( Pose & pose );

	/// @brief fill vectors with data
	void fill_vectors_with_data( Pose & pose );

	/// @brief Add output to scorefile
	void add_to_scorefile( Pose & pose );

	/// @brief Print amino acid statistics
	/// @details Doesn't print to scorefile, too many columns; prints to out
	void print_aa_statistics( Pose & pose );

	/// @brief Get size of the interface
	Real get_size( Pose & pose );

	/// @brief Get charge for all states
	utility::vector1< SSize > get_charge( Pose & pose );

	/// @brief Get avg charge for all states
	utility::vector1< Real > get_avg_charge( Pose & pose );

	/// @brief Get number of charged residues for all states
	utility::vector1< Size > get_number_charges( Pose & pose );

	/// @brief Get avg number of charged residues for all states
	utility::vector1< Real > get_avg_number_charges( Pose & pose );

	/// @brief Get hydrophobicity for all states (UHS, Koehler & Meiler, Proteins, 2008)
	utility::vector1< Real > get_hydrophobicity( Pose & pose );

	/// @brief Get avg hydrophobicity for all states (UHS, Koehler & Meiler, Proteins, 2008)
	utility::vector1< Real > get_avg_hydrophobicity( Pose & pose );

	/// @brief Get number of residues in each state, this is what we are normalizing with
	utility::vector1< Size > get_number_of_residues( Pose & pose );

	/// @brief Get surface complementarity
	Real get_surface_complementarity( Pose & pose );

	/// @brief Get avg residue size (in number of atoms) for all states
	utility::vector1< Real > get_residue_size( Pose & pose );

	/// @brief Get number of Hbonds across interface
	Size get_number_hbonds( Pose & pose );

	/// @brief Get number of pi-stacking across interface
	Size get_number_pi_stacking( Pose & pose );

	/// @brief Get in/out orientation
	Real get_inout_orientation( Pose & pose );

	/// @brief Get AA statistics
	/// @details These are counts, outer vector is over 6 bins, inner vector
	///   over 20 AAs
	utility::vector1< utility::vector1< Size > > get_aa_statistics( Pose & pose );

private: // data

	/// @brief Membrane protein?
	bool membrane_;

	/// @brief scorefunction
	core::scoring::ScoreFunctionOP sfxn_;

	/// @brief Partners;
	std::string partner_;

	/// @brief Partners split up
	utility::vector1< std::string > partners_;

	/// @brief jumps
	utility::vector1< core::Size > jumps_;

	/// @brief Interfaces
	utility::vector1< Interface > interfaces_;

	/// @brief Chains to consider for computing statistics
	std::string chains_;

	/// @brief Unique chains to minimize bias from multiple counts
	utility::vector1< bool > uniq_;

	/// @brief Are residues in the membrane?
	utility::vector1< bool > in_mem_;

	/// @brief Are the residues in the interface?
	utility::vector1< bool > intf_;

	/// @brief Are the residues exposed?
	utility::vector1< bool > exposed_;
};

//////////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
MPInterfaceStatistics::MPInterfaceStatistics() :
	Mover(),
	membrane_(),
	sfxn_(),
	partner_(),
	partners_(),
	jumps_(),
	interfaces_(),
	chains_(),
	uniq_(),
	in_mem_(),
	intf_(),
	exposed_()
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
MPInterfaceStatistics::MPInterfaceStatistics( MPInterfaceStatistics const & src ) :
	Mover( src ),
	membrane_( src.membrane_ ),
	sfxn_( src.sfxn_ ),
	partner_( src.partner_ ),
	partners_( src.partners_ ),
	jumps_( src.jumps_ ),
	interfaces_( src.interfaces_ ),
	chains_( src.chains_ ),
	uniq_( src.uniq_ ),
	in_mem_( src.in_mem_ ),
	intf_( src.intf_ ),
	exposed_( src.exposed_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
MPInterfaceStatistics &
MPInterfaceStatistics::operator=( MPInterfaceStatistics const & src )
{
	// Abort self-assignment.
	if ( this == &src ) {
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
	using namespace protocols::docking;

	TR << "Calling MPInterfaceStatistics" << std::endl;

	register_options();
	init_from_cmd();

	// if membrane protein
	if ( membrane_ == true ) {

		// add membrane
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );

		// create scorefunction and score the pose
		// this isn't really the main objective of the application, but if we are
		// writing a scorefile anyway, we might as well...
		sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );

		// get foldtree from partners and movable jumps
		// see ASCII art in protocols/membrane/util for how the foldtree is set up
		jumps_ = create_membrane_multi_partner_foldtree_anchor_tmcom( pose, partner_ );
		TR << "membrane interface jumps_ from foldtree: " << std::endl;
		for ( core::Size i = 1; i <= jumps_.size(); ++i ) {
			TR << jumps_[ i ] << " ";
		}
		TR << std::endl;

	} else {
		// for soluble protein
		// create scorefunction and score the pose
		// this isn't really the main objective of the application, but if we are
		// writing a scorefile anyway, we might as well...
		sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "talaris2014.wts" );

		// set up soluble foldtree and get interface jumps
		utility::vector1< int > jumps;
		setup_foldtree( pose, partner_, jumps );

		// cast ints in vector to core::Sizes
		for ( core::Size i = 1; i <= jumps.size(); ++i ) {
			core::Size thisjump = static_cast< core::Size >( jumps[ i ] );
			jumps_.push_back( thisjump );
		}

		TR << "soluble protein interface jumps_ from foldtree: " << std::endl;
		for ( core::Size i = 1; i <= jumps_.size(); ++i ) {
			TR << jumps_[ i ] << " ";
		}
		TR << std::endl;

	}

	( *sfxn_ )( pose );

	// get unique chains, either from commandline or unique chains from pose
	if ( chains_.size() != 0 ) {
		uniq_ = get_chains_from_cmd( pose );
	} else {
		uniq_ = compute_unique_chains( pose );
	}

	// check uniq chains for correctness: output to std::out
	utility::vector1< core::Size > chains = get_chains( pose );
	utility::vector1< core::Size > chain_last_res = chain_end_res( pose );
	for ( core::Size i = 1; i < chain_last_res.size(); ++i ) {
		if ( uniq_[ chain_last_res[ i ] ] == true ) {
			TR << "unique: " << get_chain_from_chain_id( chains[ i ], pose) << std::endl;
		} else {
			TR << "NOT unique: " << get_chain_from_chain_id( chains[ i ], pose) << std::endl;
		}
	}

	// calculate interface
	calculate_interfaces( pose );

	// fill boolean vectors with data
	fill_vectors_with_data( pose );

	// add measures to scorefile
	add_to_scorefile( pose );

	// print AA statistics
	print_aa_statistics( pose );

	// other ideas:
	// Compute a ruggedness factor of the interface???
	// Compute an average deviation from a plane along the interface
	// Get contact preferences for AAs: Yan & Jernigan, 2008, ProteinJ

} // apply

/// @brief Register Options with JD2
void MPInterfaceStatistics::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::mp::setup::spanfiles );
	option.add_relevant( OptionKeys::docking::partners );
	option.add_relevant( OptionKeys::intf::chains );

} // register options

/// @brief Initialize Mover options from the comandline
void MPInterfaceStatistics::init_from_cmd() {

	using namespace basic::options;

	membrane_ = false;

	// docking partners
	// requires all chains in the PDB file and is only for setting up the foldtree!!!
	if ( option[ OptionKeys::docking::partners ].user() ) {
		partner_ = option[ OptionKeys::docking::partners ]();
	} else {
		utility_exit_with_message( "Please provide -docking:partners flag!" );
	}

	// membrane protein?
	if ( option[ OptionKeys::mp::setup::spanfiles ].user() ) {
		membrane_ = true;
	}

	// which chains to consider for getting statistics
	if ( option[ OptionKeys::intf::chains ].user() ) {
		chains_ = option[ OptionKeys::intf::chains ]();
	}

} // init from commandline

/// @brief Get the chains from the commandline
utility::vector1< bool > MPInterfaceStatistics::get_chains_from_cmd( Pose & pose ) {

	// set residue vector to false, i.e. not unique
	utility::vector1< bool > tmp ( pose.total_residue(), false );

	// go through pose residues
	for ( core::Size r = 1; r <= pose.total_residue(); ++r ) {

		// get chain
		core::Size res_chain = pose.residue( r ).chain();

		// go through chains in string, 0-indexed
		for ( core::Size c = 0; c < chains_.size(); ++c ) {

			// get chain
			core::Size chain = get_chain_id_from_chain( chains_[ c ], pose );

			// if residue chain is the same as in user-provided string, set to true
			if ( res_chain == chain ) {
				tmp[ r ] = true;
			}
		} // chains
	} // residues

	return tmp;
} // get chains from commandline

/////////////////////////////////
///// Rosetta Scripts Methods ///
/////////////////////////////////
//
///// @brief Create a Clone of this mover
//protocols::moves::MoverOP
//MPInterfaceStatistics::clone() const {
// return ( protocols::moves::MoverOP( new MPInterfaceStatistics( *this ) ) );
//}
//
///// @brief Create a Fresh Instance of this Mover
//protocols::moves::MoverOP
//MPInterfaceStatistics::fresh_instance() const {
// return protocols::moves::MoverOP( new MPInterfaceStatistics() );
//}
//
///// @brief Pase Rosetta Scripts Options for this Mover
//void
//MPInterfaceStatistics::parse_my_tag(
//            utility::tag::TagCOP tag,
//            basic::datacache::DataMap &,
//            protocols::filters::Filters_map const &,
//            protocols::moves::Movers_map const &,
//            core::pose::Pose const &
//            ) {
//  // TODO
//
//}
//
///// @brief Create a new copy of this mover
//protocols::moves::MoverOP
//MPInterfaceStatisticsCreator::create_mover() const {
// return protocols::moves::MoverOP( new MPInterfaceStatistics );
//}
//
///// @brief Return the Name of this mover (as seen by Rscripts)
//std::string
//MPInterfaceStatisticsCreator::keyname() const {
// return MPInterfaceStatisticsCreator::mover_name();
//}
//
///// @brief Mover name for Rosetta Scripts
//std::string
//MPInterfaceStatisticsCreator::mover_name() {
// return "MPInterfaceStatistics";
//}

//////////////////////////////////////////////////////////////////////

typedef utility::pointer::shared_ptr< MPInterfaceStatistics > MPInterfaceStatisticsOP;

//////////////////////////////////////////////////////////////////////

/// @brief Get size of the interface in square Angstrom
/// @details Supports multiple partners; the interface is sum of all interfaces,
///   independent of uniq partners
Real MPInterfaceStatistics::get_size( Pose & pose ) {

	using namespace numeric;
	using namespace core::pose;
	using namespace core::scoring::sasa;
	using namespace protocols::scoring;

	// split pose by chains
	utility::vector1< PoseOP > pose_chains( pose.split_by_chain() );
	utility::vector1< core::Size > chains( get_chains( pose ) );
	utility::vector1< Pose > pose_partners;

	// go through partners and append the chains into subposes to fit the partners
	for ( core::Size p = 1; p <= partners_.size(); ++p ) {

		Pose subpose;

		// go through chains in partner, string indexes from 0
		for ( core::Size c = 0; c < partners_[ p ].size(); ++c ) {

			// get chainID from chain
			core::Size chainid( get_chain_id_from_chain( partners_[ p ][ c ], pose ) );

			// go through pose chains to find pose chain in vector
			for ( core::Size i = 1; i <= pose_chains.size(); ++i ) {

				// if the chainids match, append chain to subpose
				if ( chainid == chains[ i ] ) {
					append_pose_to_pose( subpose, *pose_chains[ i ], true );
				}

			} // chains
		} // chains in partner

		// add subpose to vector
		pose_partners.push_back( subpose );

	} // partners

	// get total SASA
	SasaCalc calc = SasaCalc();
	core::Real total_sasa( calc.calculate( pose ) );

	// start with total SASA for interface
	core::Real intf_sasa = total_sasa;

	// go through partners and subtract from SASA
	for ( core::Size p = 1; p <= pose_partners.size(); ++p ) {

		// compute SASA for that partner
		core::Real p_sasa( calc.calculate( pose_partners[ p ] ) );

		// subtract from total SASA
		intf_sasa -= p_sasa;
	}

	// compute: intf_sasa = - ( total_sasa - sum( p_sasa ) ) / 2;
	intf_sasa /= -2;

	return intf_sasa;

}// get_size

//////////////////////////////////////////////////////////////////////

/// @brief Get charge of the interface
utility::vector1< SSize > MPInterfaceStatistics::get_charge( Pose & pose ) {

	SSize tm_int( 0 );
	SSize tm_nonint( 0 );
	SSize tm_core( 0 );
	SSize sol_int( 0 );
	SSize sol_nonint( 0 );
	SSize sol_core( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		std::string res = pose.residue( i ).name3();

		// is residue positive?
		bool pos( false );
		if ( res == "ARG" || res == "LYS" ) {
			pos = true;
		}

		// is residue negative?
		bool neg( false );
		if ( res == "ASP" || res == "GLU" ) {
			neg = true;
		}

		// interface, membrane, positive / negative
		if ( uniq_[i] && intf_[i] && in_mem_[i] && pos ) {
			++tm_int;
		} else if ( uniq_[i] && intf_[i] && in_mem_[i] && neg ) {
			--tm_int;
		} else if ( uniq_[i] && intf_[i] && !in_mem_[i] && pos ) {
			// interface, sol, positive / negative
			++sol_int;
		} else if ( uniq_[i] && intf_[i] && !in_mem_[i] && neg ) {
			--sol_int;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && in_mem_[i] && pos ) {
			// not interface, membrane, positive / negative
			++tm_nonint;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && in_mem_[i] && neg ) {
			--tm_nonint;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && !in_mem_[i] && pos ) {
			// not interface, sol, positive / negative
			++sol_nonint;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && !in_mem_[i] && neg ) {
			--sol_nonint;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && in_mem_[i] && pos ) {
			// if core, membrane, positive / negative
			++tm_core;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && in_mem_[i] && neg ) {
			--tm_core;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && !in_mem_[i] && pos ) {
			// if core, solution, positive / negative
			++sol_core;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && !in_mem_[i] && neg ) {
			--sol_core;
		}

	}

	// write counts into vector
	utility::vector1< SSize > charge;
	charge.push_back( tm_int );
	charge.push_back( tm_nonint );
	charge.push_back( tm_core );
	charge.push_back( sol_int );
	charge.push_back( sol_nonint );
	charge.push_back( sol_core );

	return charge;

} // get charge

//////////////////////////////////////////////////////////////////////

/// @brief Get average charge of the interface, meaning averaged over the
///   number of residues
utility::vector1< Real > MPInterfaceStatistics::get_avg_charge( Pose & pose ) {

	utility::vector1< SSize > charge( get_charge( pose ) );
	utility::vector1< Size > nres( get_number_of_residues( pose ) );
	utility::vector1< Real > avg_charge;

	// get average
	for ( Size i = 1; i <= nres.size(); ++i ) {
		avg_charge.push_back( static_cast< Real >( charge[ i ] ) / static_cast< Real >( nres[ i ] ) );
	}

	return avg_charge;

} // get avg charge

//////////////////////////////////////////////////////////////////////

/// @brief Get number of charged residues in the interface
utility::vector1< Size > MPInterfaceStatistics::get_number_charges( Pose & pose ) {

	Size tm_int( 0 );
	Size tm_nonint( 0 );
	Size tm_core( 0 );
	Size sol_int( 0 );
	Size sol_nonint( 0 );
	Size sol_core( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		std::string res = pose.residue( i ).name3();

		// is residue charged?
		bool charged( false );
		if ( res == "ARG" || res == "LYS" || res == "ASP" || res == "GLU" ) {
			charged = true;
		}

		// interface, membrane, charged
		if ( uniq_[i] && intf_[i] && in_mem_[i] && charged ) {
			++tm_int;
		} else if ( uniq_[i] && intf_[i] && !in_mem_[i] && charged ) {
			// interface, sol, charged
			++sol_int;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && in_mem_[i] && charged ) {
			// not interface, membrane, charged
			++tm_nonint;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && !in_mem_[i] && charged ) {
			// not interface, sol, charged
			++sol_nonint;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && in_mem_[i] && charged ) {
			// core, membrane
			++tm_core;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && !in_mem_[i] && charged ) {
			// core, membrane
			++sol_core;
		}

	}

	// write counts into vector
	utility::vector1< Size > ncharges;
	ncharges.push_back( tm_int );
	ncharges.push_back( tm_nonint );
	ncharges.push_back( tm_core );
	ncharges.push_back( sol_int );
	ncharges.push_back( sol_nonint );
	ncharges.push_back( sol_core );

	return ncharges;

} // get number of charges

//////////////////////////////////////////////////////////////////////

/// @brief Get average number of charges, meaning averaged over the
///   number of residues
utility::vector1< Real > MPInterfaceStatistics::get_avg_number_charges( Pose & pose ) {

	utility::vector1< Size > ncharges( get_number_charges( pose ) );
	utility::vector1< Size > nres( get_number_of_residues( pose ) );
	utility::vector1< Real > avg_ncharges;

	// get average
	for ( Size i = 1; i <= nres.size(); ++i ) {
		avg_ncharges.push_back( static_cast< Real >( ncharges[ i ] ) / static_cast< Real >( nres[ i ] ) );
	}

	return avg_ncharges;

} // get avg charge

//////////////////////////////////////////////////////////////////////

/// @brief Get hydrophobicity of the interface (UHS, Koehler & Meiler, Proteins, 2008)
utility::vector1< Real > MPInterfaceStatistics::get_hydrophobicity( Pose & pose ) {

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

	Real tm_int( 0 );
	Real tm_nonint( 0 );
	Real tm_core( 0 );
	Real sol_int( 0 );
	Real sol_nonint( 0 );
	Real sol_core( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		// interface, mem
		if ( uniq_[i] && intf_[i] && in_mem_[i] ) {
			tm_int += hydro[ pose.residue( i ).name1() ];
		} else if ( uniq_[i] && intf_[i] && !in_mem_[i] ) {
			// interface, sol
			sol_int += hydro[ pose.residue( i ).name1() ];
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && in_mem_[i] ) {
			// not interface, mem
			tm_nonint += hydro[ pose.residue( i ).name1() ];
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && !in_mem_[i] ) {
			// not interface, sol
			sol_nonint += hydro[ pose.residue( i ).name1() ];
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && in_mem_[i] ) {
			// core, mem
			tm_core += hydro[ pose.residue( i ).name1() ];
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && !in_mem_[i] ) {
			// core, sol
			sol_core += hydro[ pose.residue( i ).name1() ];
		}
	}

	// write counts into vector
	utility::vector1< Real > total_hydro;
	total_hydro.push_back( tm_int );
	total_hydro.push_back( tm_nonint );
	total_hydro.push_back( tm_core );
	total_hydro.push_back( sol_int );
	total_hydro.push_back( sol_nonint );
	total_hydro.push_back( sol_core );

	return total_hydro;

} // get hydrophobicity

//////////////////////////////////////////////////////////////////////

/// @brief Get avg hydrophobicity for all states (UHS, Koehler & Meiler, Proteins, 2008)
utility::vector1< Real > MPInterfaceStatistics::get_avg_hydrophobicity( Pose & pose ) {

	utility::vector1< Real > hydro( get_hydrophobicity( pose ) );
	utility::vector1< Size > nres( get_number_of_residues( pose ) );
	utility::vector1< Real > avg_hydro;

	// get average
	for ( Size i = 1; i <= nres.size(); ++i ) {
		avg_hydro.push_back( hydro[ i ] / nres[ i ] );
	}

	return avg_hydro;

} // average hydrophobicity

//////////////////////////////////////////////////////////////////////


/// @brief Get number of residues in each bin, this is what we are normalizing by
utility::vector1< Size > MPInterfaceStatistics::get_number_of_residues( Pose & pose ) {

	Size tm_int( 0 );
	Size tm_nonint( 0 );
	Size tm_core( 0 );
	Size sol_int( 0 );
	Size sol_nonint( 0 );
	Size sol_core( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		// interface, mem
		if ( uniq_[i] && intf_[i] && in_mem_[i] ) {
			++tm_int;
		} else if ( uniq_[i] && intf_[i] && !in_mem_[i] ) {
			// interface, sol
			++sol_int;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && in_mem_[i] ) {
			// non interface, mem
			++tm_nonint;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && !in_mem_[i] ) {
			// not interface, sol
			++sol_nonint;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && in_mem_[i] ) {
			// core, mem
			++tm_core;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && !in_mem_[i] ) {
			// core, sol
			++sol_core;
		}
	}

	// write counts into vector
	utility::vector1< Size > nresidues;
	nresidues.push_back( tm_int );
	nresidues.push_back( tm_nonint );
	nresidues.push_back( tm_core );
	nresidues.push_back( sol_int );
	nresidues.push_back( sol_nonint );
	nresidues.push_back( sol_core );

	return nresidues;

} // get number of contacts in interface

//////////////////////////////////////////////////////////////////////

/// @brief Get surface complementarity
Real MPInterfaceStatistics::get_surface_complementarity( Pose & pose ) {

	using namespace core::scoring::sc;
	ShapeComplementarityCalculator shape;
	shape.Init();

	// takes jumpID and whether a quick computation or not
	return shape.CalcSc( pose, jumps_[1], 1 );

} // get surface complementarity of interface

//////////////////////////////////////////////////////////////////////

/// @brief Get avg residue in the interface
utility::vector1 < Real > MPInterfaceStatistics::get_residue_size( Pose & pose ) {

	Real tm_int( 0 );
	Real tm_nonint( 0 );
	Real tm_core( 0 );
	Real sol_int( 0 );
	Real sol_nonint( 0 );
	Real sol_core( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		// interface, mem
		if ( uniq_[i] && intf_[i] && in_mem_[i] ) {
			tm_int += static_cast< Real >( pose.residue( i ).natoms() );
		} else if ( uniq_[i] && intf_[i] && !in_mem_[i] ) {
			// interface, sol
			sol_int += static_cast< Real >( pose.residue( i ).natoms() );
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && in_mem_[i] ) {
			// not interface, mem
			tm_nonint += static_cast< Real >( pose.residue( i ).natoms() );
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && !in_mem_[i] ) {
			// not interface, sol
			sol_nonint += static_cast< Real >( pose.residue( i ).natoms() );
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && in_mem_[i] ) {
			// core, mem
			tm_core += static_cast< Real >( pose.residue( i ).natoms() );
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && !in_mem_[i] ) {
			// core, sol
			sol_core += static_cast< Real >( pose.residue( i ).natoms() );
		}
	}

	utility::vector1< Size > nres( get_number_of_residues( pose ) );

	TR << "number of residues:" << std::endl;
	TR << "tm int " << nres[1] << std::endl;
	TR << "tm nonint " << nres[2] << std::endl;
	TR << "tm core " << nres[3] << std::endl;
	TR << "sol int " << nres[4] << std::endl;
	TR << "sol nonint " << nres[5] << std::endl;
	TR << "sol core " << nres[6] << std::endl;

	TR << "residue size:" << std::endl;
	TR << "tm int " << tm_int << std::endl;
	TR << "tm nonint " << tm_nonint << std::endl;
	TR << "tm core " << tm_core << std::endl;
	TR << "sol int " << sol_int << std::endl;
	TR << "sol nonint " << sol_nonint << std::endl;
	TR << "sol core " << sol_core << std::endl;

	// get average residue size
	tm_int /= static_cast< Real > ( nres[1] );
	tm_nonint /= static_cast< Real > ( nres[2] );
	tm_core /= static_cast< Real > ( nres[3] );
	sol_int /= static_cast< Real > ( nres[4] );
	sol_nonint /= static_cast< Real > ( nres[5] );
	sol_core /= static_cast< Real > ( nres[6] );

	// write counts into vector
	utility::vector1< Real > avg_size;
	avg_size.push_back( tm_int );
	avg_size.push_back( tm_nonint );
	avg_size.push_back( tm_core );
	avg_size.push_back( sol_int );
	avg_size.push_back( sol_nonint );
	avg_size.push_back( sol_core );

	return avg_size;

} // get residue size

//////////////////////////////////////////////////////////////////////

/// @brief Get number of Hbonds across interface
Size MPInterfaceStatistics::get_number_hbonds( Pose & pose ) {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	// get the partners
	utility::vector1< std::string > partners( utility::string_split( partner_, '_' ) );
	utility::vector1< std::string > partner1( utility::split( partners[1] ) );
	utility::vector1< std::string > partner2( utility::split( partners[2] ) );

	// get the Hbonds
	hbonds::HBondDatabaseCOP hb_database(hbonds::HBondDatabase::get_database());
	hbonds::HBondSet hb_set;
	pose.update_residue_neighbors();
	hbonds::fill_hbond_set( pose, false, hb_set, false );
	// hb_set.show(pose, true, std::cout);

	Size nhbonds( 0 );

	// go through Hbonds
	for ( Size i = 1; i <= hb_set.nhbonds(); ++i ) {

		HBond hbond = hb_set.hbond( i );
		Size don = hbond.don_res();
		Size acc = hbond.acc_res();

		// check whether donor/acc is in partner1
		// doesn't check whether the residues are defined as interface (they should anyway because they are between the two partners)

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

		//  TR << "hbond: don " << don << ", acc " << acc;
		//  TR << ", don: " << don_in_1 << ", acc: " << acc_in_1 << std::endl;

		// if donor in 1 and acceptor not, then we have Hbond across interface
		if ( don_in_1 == true && acc_in_1 == false ) {
			++nhbonds;
		} else if ( don_in_1 == false && acc_in_1 == true ) {
			// similarly, if donor not in 1 but acc is, we also have Hbond across interface
			++nhbonds;
		}
	} // go through Hbonds


	return nhbonds;
} // number hbonds across interface

//////////////////////////////////////////////////////////////////////

/// @brief Get number of pi-stacking across interface
Size MPInterfaceStatistics::get_number_pi_stacking( Pose & pose ) {

	core::pose::metrics::PoseMetricCalculatorOP pi_pi_calculator ( new protocols::toolbox::pose_metric_calculators::PiPiCalculator() );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "pi_pi", pi_pi_calculator );

	TR << "WARNING: PLEASE CHECK THE NUMBER OF YOUR Pi-stacking INTERACTIONS VISUALLY TO CONFIRM!" << std::endl;

	return utility::string2Size( pose.print_metric( "pi_pi", "pi_pi" ) );
}

//////////////////////////////////////////////////////////////////////

/// @brief Get in/out orientation
Real MPInterfaceStatistics::get_inout_orientation( Pose & pose ) {

	if ( pose.residue( nres_protein( pose ) ).atom("CA").xyz().z() > 0.0 ) {
		return 0.0;
	}
	return 1.0;

}// get in/out orientation

//////////////////////////////////////////////////////////////////////

/// @brief Get AA statistics
/// @details These are counts, outer vector is over 6 bins, inner vector
///   over 20 AAs
utility::vector1< utility::vector1< Size > >
MPInterfaceStatistics::get_aa_statistics( Pose & pose ) {

	// initialize vectors with 20 zeros
	utility::vector1< Size > tm_int( 20, 0 );
	utility::vector1< Size > tm_nonint( 20, 0 );
	utility::vector1< Size > tm_core( 20, 0 );
	utility::vector1< Size > sol_int( 20, 0 );
	utility::vector1< Size > sol_nonint( 20, 0 );
	utility::vector1< Size > sol_core( 20, 0 );

	// create a max of the amino acids and it's index in the vector
	std::map< char, Size > aa_index;

	// polar
	aa_index[ 'C' ] = 1;
	aa_index[ 'N' ] = 2;
	aa_index[ 'Q' ] = 3;
	aa_index[ 'S' ] = 4;
	aa_index[ 'T' ] = 5;

	// charged
	aa_index[ 'D' ] = 6;
	aa_index[ 'E' ] = 7;
	aa_index[ 'K' ] = 8;
	aa_index[ 'R' ] = 9;

	// hydrophobic
	aa_index[ 'A' ] = 10;
	aa_index[ 'G' ] = 11;
	aa_index[ 'I' ] = 12;
	aa_index[ 'L' ] = 13;
	aa_index[ 'M' ] = 14;
	aa_index[ 'P' ] = 15;
	aa_index[ 'V' ] = 16;

	// aromatic
	aa_index[ 'F' ] = 17;
	aa_index[ 'H' ] = 18;
	aa_index[ 'W' ] = 19;
	aa_index[ 'Y' ] = 20;

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		// increase counter at the correct place
		if ( uniq_[i] && intf_[i] && in_mem_[i] ) {
			tm_int[ aa_index[ pose.residue( i ).name1() ] ] += 1;
		} else if ( uniq_[i] && intf_[i] && !in_mem_[i] ) {
			// interface, sol
			sol_int[ aa_index[ pose.residue( i ).name1() ] ] += 1;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && in_mem_[i] ) {
			// not interface, mem
			tm_nonint[ aa_index[ pose.residue( i ).name1() ] ] += 1;
		} else if ( uniq_[i] && !intf_[i] && exposed_[i] && !in_mem_[i] ) {
			// not interface, sol
			sol_nonint[ aa_index[ pose.residue( i ).name1() ] ] += 1;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && in_mem_[i] ) {
			// core, mem
			tm_core[ aa_index[ pose.residue( i ).name1() ] ] += 1;
		} else if ( uniq_[i] && !intf_[i] && !exposed_[i] && !in_mem_[i] ) {
			// core, sol
			sol_core[ aa_index[ pose.residue( i ).name1() ] ] += 1;
		}
	}

	// define overall vector
	utility::vector1< utility::vector1< Size > > aa_stats;

	// write counts into vector
	aa_stats.push_back( tm_int );
	aa_stats.push_back( tm_nonint );
	aa_stats.push_back( tm_core );
	aa_stats.push_back( sol_int );
	aa_stats.push_back( sol_nonint );
	aa_stats.push_back( sol_core );

	return aa_stats;

} // get AA statistics

//////////////////////////////////////////////////////////////////////

/// @brief Add output to scorefile
void MPInterfaceStatistics::add_to_scorefile( Pose & pose ) {

	// get job
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	// get size of the interface
	Real size = get_size( pose );
	job->add_string_real_pair( "Isize", size );
	TR << "intf size " << size << std::endl;

	// Get surface complementarity
	Real shape = get_surface_complementarity( pose );
	job->add_string_real_pair( "shape", shape );
	TR << "shape_compl " << shape << std::endl;

	// Get number of Hbonds across interface
	Size nhbonds = get_number_hbonds( pose );
	job->add_string_real_pair( "nhbonds", nhbonds );
	TR << "nhbonds " << nhbonds << std::endl;

	// Get number of pi-stacking across interface
	Size npi = get_number_pi_stacking( pose );
	job->add_string_real_pair( "npi", npi );
	TR << "pi stacking " << npi << std::endl;

	// Get in/out orientation
	Real inout = get_inout_orientation( pose );
	job->add_string_real_pair( "inout", inout );
	TR << "inout " << inout << std::endl;

	// Get charge in the interface
	utility::vector1< SSize > charge = get_charge( pose );
	job->add_string_real_pair( "chrg_tm_int", charge[1] );
	job->add_string_real_pair( "chrg_tm_nint", charge[2] );
	job->add_string_real_pair( "chrg_tm_core", charge[3] );
	job->add_string_real_pair( "chrg_so_int", charge[4] );
	job->add_string_real_pair( "chrg_so_nint", charge[5] );
	job->add_string_real_pair( "chrg_so_core", charge[6] );
	TR << "charge tm int " << charge[1] << std::endl;
	TR << "charge tm nonint " << charge[2] << std::endl;
	TR << "charge tm core " << charge[3] << std::endl;
	TR << "charge sol int " << charge[4] << std::endl;
	TR << "charge sol nonint " << charge[5] << std::endl;
	TR << "charge sol core " << charge[6] << std::endl;

	// Get avg charge in the interface
	utility::vector1< Real > avg_charge = get_avg_charge( pose );
	job->add_string_real_pair( "avg_chrg_tm_int", avg_charge[1] );
	job->add_string_real_pair( "avg_chrg_tm_nint", avg_charge[2] );
	job->add_string_real_pair( "avg_chrg_tm_core", avg_charge[3] );
	job->add_string_real_pair( "avg_chrg_so_int", avg_charge[4] );
	job->add_string_real_pair( "avg_chrg_so_nint", avg_charge[5] );
	job->add_string_real_pair( "avg_chrg_so_core", avg_charge[6] );
	TR << "avg charge tm int " << avg_charge[1] << std::endl;
	TR << "avg charge tm nonint " << avg_charge[2] << std::endl;
	TR << "avg charge tm core " << avg_charge[3] << std::endl;
	TR << "avg charge sol int " << avg_charge[4] << std::endl;
	TR << "avg charge sol nonint " << avg_charge[5] << std::endl;
	TR << "avg charge sol core " << avg_charge[6] << std::endl;

	// Get number of charged residues in the interface
	utility::vector1< Size > ncharges = get_number_charges( pose );
	job->add_string_real_pair( "nchrg_tm_int", ncharges[1] );
	job->add_string_real_pair( "nchrg_tm_nint", ncharges[2] );
	job->add_string_real_pair( "nchrg_tm_core", ncharges[3] );
	job->add_string_real_pair( "nchrg_so_int", ncharges[4] );
	job->add_string_real_pair( "nchrg_so_nint", ncharges[5] );
	job->add_string_real_pair( "nchrg_so_core", ncharges[6] );
	TR << "ncharges tm int " << ncharges[1] << std::endl;
	TR << "ncharges tm nonint " << ncharges[2] << std::endl;
	TR << "ncharges tm core " << ncharges[3] << std::endl;
	TR << "ncharges sol int " << ncharges[4] << std::endl;
	TR << "ncharges sol nonint " << ncharges[5] << std::endl;
	TR << "ncharges sol core " << ncharges[6] << std::endl;

	// Get avg number of charged residues in the interface
	utility::vector1< Real > avg_ncharges = get_avg_number_charges( pose );
	job->add_string_real_pair( "avg_nchrg_tm_int", avg_ncharges[1] );
	job->add_string_real_pair( "avg_nchrg_tm_nint", avg_ncharges[2] );
	job->add_string_real_pair( "avg_nchrg_tm_core", avg_ncharges[3] );
	job->add_string_real_pair( "avg_nchrg_so_int", avg_ncharges[4] );
	job->add_string_real_pair( "avg_nchrg_so_nint", avg_ncharges[5] );
	job->add_string_real_pair( "avg_nchrg_so_core", avg_ncharges[6] );
	TR << "avg ncharges tm int " << avg_ncharges[1] << std::endl;
	TR << "avg ncharges tm nonint " << avg_ncharges[2] << std::endl;
	TR << "avg ncharges tm core " << avg_ncharges[3] << std::endl;
	TR << "avg ncharges sol int " << avg_ncharges[4] << std::endl;
	TR << "avg ncharges sol nonint " << avg_ncharges[5] << std::endl;
	TR << "avg ncharges sol core " << avg_ncharges[6] << std::endl;

	// Get hydrophobicity of the interface
	utility::vector1< Real > hydro = get_hydrophobicity( pose );
	job->add_string_real_pair( "hydro_tm_int", hydro[1] );
	job->add_string_real_pair( "hydro_tm_nint", hydro[2] );
	job->add_string_real_pair( "hydro_tm_core", hydro[3] );
	job->add_string_real_pair( "hydro_so_int", hydro[4] );
	job->add_string_real_pair( "hydro_so_nint", hydro[5] );
	job->add_string_real_pair( "hydro_so_core", hydro[6] );
	TR << "hydro tm int " << hydro[1] << std::endl;
	TR << "hydro tm nonint " << hydro[2] << std::endl;
	TR << "hydro tm core " << hydro[3] << std::endl;
	TR << "hydro sol int " << hydro[4] << std::endl;
	TR << "hydro sol nonint " << hydro[5] << std::endl;
	TR << "hydro sol core " << hydro[6] << std::endl;

	// Get average hydrophobicity of the interface
	utility::vector1< Real > avg_hydro = get_avg_hydrophobicity( pose );
	job->add_string_real_pair( "avg_hydro_tm_int", avg_hydro[1] );
	job->add_string_real_pair( "avg_hydro_tm_nint", avg_hydro[2] );
	job->add_string_real_pair( "avg_hydro_tm_core", avg_hydro[3] );
	job->add_string_real_pair( "avg_hydro_so_int", avg_hydro[4] );
	job->add_string_real_pair( "avg_hydro_so_nint", avg_hydro[5] );
	job->add_string_real_pair( "avg_hydro_so_core", avg_hydro[6] );
	TR << "avg hydro tm int " << avg_hydro[1] << std::endl;
	TR << "avg hydro tm nonint " << avg_hydro[2] << std::endl;
	TR << "avg hydro tm core " << avg_hydro[3] << std::endl;
	TR << "avg hydro sol int " << avg_hydro[4] << std::endl;
	TR << "avg hydro sol nonint " << hydro[5] << std::endl;
	TR << "avg hydro sol core " << avg_hydro[6] << std::endl;

	// Get number of contacts
	utility::vector1< Size > nresidues = get_number_of_residues( pose );
	job->add_string_real_pair( "nres_tm_int", nresidues[1] );
	job->add_string_real_pair( "nres_tm_nint", nresidues[2] );
	job->add_string_real_pair( "nres_tm_core", nresidues[3] );
	job->add_string_real_pair( "nres_so_int", nresidues[4] );
	job->add_string_real_pair( "nres_so_nint", nresidues[5] );
	job->add_string_real_pair( "nres_so_core", nresidues[6] );
	TR << "nresidues tm int " << nresidues[1] << std::endl;
	TR << "nresidues tm nonint " << nresidues[2] << std::endl;
	TR << "nresidues tm core " << nresidues[3] << std::endl;
	TR << "nresidues sol int " << nresidues[4] << std::endl;
	TR << "nresidues sol nonint " << nresidues[5] << std::endl;
	TR << "nresidues sol core " << nresidues[6] << std::endl;

	// Get size of residues in the interface
	utility::vector1< Real > res_size = get_residue_size( pose );
	job->add_string_real_pair( "ressize_tm_int", res_size[1] );
	job->add_string_real_pair( "ressize_tm_nint", res_size[2] );
	job->add_string_real_pair( "ressize_tm_core", res_size[3] );
	job->add_string_real_pair( "ressize_so_int", res_size[4] );
	job->add_string_real_pair( "ressize_so_nint", res_size[5] );
	job->add_string_real_pair( "ressize_so_core", res_size[6] );
	TR << "avg res size tm int " << res_size[1] << std::endl;
	TR << "avg res size tm nonint " << res_size[2] << std::endl;
	TR << "avg res size tm core " << res_size[3] << std::endl;
	TR << "avg res size sol int " << res_size[4] << std::endl;
	TR << "avg res size sol nonint " << res_size[5] << std::endl;
	TR << "avg res size sol core " << res_size[6] << std::endl;

} // add_to_scorefile

/// @brief Print amino acid statistics
/// @details Doesn't print to scorefile, too many columns; prints to out
void MPInterfaceStatistics::print_aa_statistics( Pose & pose ) {

	using namespace ObjexxFCL::format;

	// get statistics
	utility::vector1< utility::vector1< Size > > aa_cnts = get_aa_statistics( pose );

	TR << ">region      C    N    Q    S    T    D    E    K    R    A    G    I    L    M    P    V    F    H    W    Y" << std::endl;

	// iterate over regions
	for ( Size i = 1; i <= aa_cnts.size(); ++i ) {

		// print regions
		i == 1 ? TR << ">tm-int  " :
			i == 2 ? TR << ">tm-nint " :
			i == 3 ? TR << ">tm-core " :
			i == 4 ? TR << ">so-int  " :
			i == 5 ? TR << ">so-nint " : TR << ">so-core ";

		// iterate over amino acids
		for ( Size aa = 1; aa <= aa_cnts[i].size(); ++aa ) {
			TR << I( 5, aa_cnts[ i ][ aa ] );
		}
		TR << std::endl;
	}

} // print AA statistics

/// @brief Fill vectors with data
void MPInterfaceStatistics::fill_vectors_with_data( Pose & pose ) {

	using namespace core::scoring::sasa;

	// get relative per residue SASA
	utility::vector1< Real > sasa( rel_per_res_sc_sasa( pose ) );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		// is residue in the membrane?
		// default for soluble proteins, special for membrane protein
		bool in_mem( false );
		if ( membrane_ == true ) {
			in_mem = pose.conformation().membrane_info()->in_membrane( i );
		}
		in_mem_.push_back( in_mem );

		// is the residue exposed?
		bool exposed( false );
		if ( sasa[i] >= 0.5 ) {
			exposed = true;
		}
		exposed_.push_back( exposed );

		TR << "res " << i << " " << pose.residue( i ).name3() << " uniq " << uniq_[ i ] << " mem " << in_mem << " exposed " << exposed << " sasa " << sasa[ i ] << " intf " << intf_[ i ] << std::endl;
	}

} // fill vectors with data

/// @brief Calculate interfaces
/// @detail This is for 2-body docking
void MPInterfaceStatistics::calculate_interfaces( Pose & pose ) {

	// fill partners_ vector
	partners_ = utility::string_split( partner_, '_' );

	// initialize intf_ vector with false for each residue
	utility::vector1< bool > tmp( pose.total_residue(), false );
	intf_ = tmp;

	// go through jumps and compute interfaces
	for ( core::Size i = 1; i <= jumps_.size(); ++i ) {

		// compute interface
		Interface intf = Interface( jumps_[ i ] );
		intf.calculate( pose );
		interfaces_.push_back( intf );

		// go through residues and set interface residue to true
		for ( core::Size r = 1; r <= nres_protein( pose ); ++r ) {
			bool is_intf = intf.is_interface( r );
			if ( is_intf == true ) {
				intf_[ r ] = true;
			}
		}
	}

} // calculate interfaces

////////////////////////////// MAIN ////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::jd2;

		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		MPInterfaceStatisticsOP mpis( new MPInterfaceStatistics() );
		JobDistributor::get_instance()->go(mpis);
	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;
}
