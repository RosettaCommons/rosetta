// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Package Headers
//#include <devel/noesy_assign/Exceptions.hh>

// Project Headers

// Utility headers
#include <basic/Tracer.hh>

// Utility headers
#include <basic/options/option_macros.hh>


//// C++ headers
#include <cstdlib>
#include <string>





OPT_1GRP_KEY( Real, noesy_weights, chemshift )
//OPT_1GRP_KEY( Real, noesy_weights, network_high )

OPT_1GRP_KEY( Real, noesy_weights, Vmin )
OPT_1GRP_KEY( Real, noesy_weights, symmetry )
OPT_1GRP_KEY( Real, noesy_weights, covalent )
OPT_1GRP_KEY( Real, noesy_weights, decoys )
OPT_1GRP_KEY( Real, noesy_weights, Smax )
OPT_1GRP_KEY( Real, noesy_weights, dcut )
OPT_1GRP_KEY( Real, noesy_weights, network_min )
OPT_1GRP_KEY( Real, noesy_weights, network_atom_min )
OPT_1GRP_KEY( Real, noesy_weights, dcalibrate )
OPT_1GRP_KEY( Real, noesy_weights, calibration_target )

OPT_1GRP_KEY( Boolean, noesy, atom_dependent_calibration )
OPT_1GRP_KEY( Boolean, noesy, ignore_resonancefile_tolerances )

OPT_2GRP_KEY( RealVector, noesy_weights, defaults, Vmin )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, symmetry )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, covalent )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, decoys )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, Smax )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, dcut )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, network_min )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, network_atom_min )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, dcalibrate )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, calibration_target )


OPT_1GRP_KEY( Real, noesy_weights, elim_dist_viol )
OPT_1GRP_KEY( Real, noesy_weights, centroid_padding )

OPT_1GRP_KEY( Real, noesy_weights, cst_strength )
OPT_1GRP_KEY( Integer, noesy_weights, cycle )

OPT_1GRP_KEY( Boolean, noesy, no_network )
OPT_2GRP_KEY( Boolean, noesy, network, include_reverse_dir )
OPT_2GRP_KEY( Boolean, noesy, network, allow_same_residue_connect )
OPT_2GRP_KEY( Boolean, noesy, network, use_all_covalent_atoms )

OPT_1GRP_KEY( Boolean, noesy, use_local_distviol )
OPT_2GRP_KEY( Real, noesy, local_distviol, range )
OPT_2GRP_KEY( Real, noesy, local_distviol, global_buffer )

bool protocols::noesy_assign::PeakAssignmentParameters::options_registered_( false );
core::Size protocols::noesy_assign::PeakAssignmentParameters::cycle_selector_( 0 );

protocols::noesy_assign::PeakAssignmentParameters* protocols::noesy_assign::PeakAssignmentParameters::instance_( NULL );

protocols::noesy_assign::PeakAssignmentParameters const*
protocols::noesy_assign::PeakAssignmentParameters::get_instance() {
  if ( instance_ ) return instance_;
  instance_ = new PeakAssignmentParameters;
  instance_->set_options_from_cmdline();
  return instance_;
}

protocols::noesy_assign::PeakAssignmentParameters*
protocols::noesy_assign::PeakAssignmentParameters::get_nonconst_instance() {
  if ( instance_ ) return instance_;
  instance_ = new PeakAssignmentParameters;
  instance_->set_options_from_cmdline();
  return instance_;
}

void protocols::noesy_assign::PeakAssignmentParameters::set_cycle( core::Size cycle ) {
	cycle_selector_ = cycle;
	if ( instance_ ) delete instance_;
	instance_ = NULL;
}

void protocols::noesy_assign::PeakAssignmentParameters::register_options() {
  using namespace basic::options;
  using namespace OptionKeys;
  if ( options_registered_ ) return;

	//cycle independent
  NEW_OPT( noesy_weights::chemshift, "contribution of chem. shift overlap to peak volume", 0.5 );
	//	NEW_OPT( noesy_weights::network_high, "contribution of network anchoring to peak volume", 10.0 );
	NEW_OPT( noesy::ignore_resonancefile_tolerances, "ignore the tolerances in the resonance file", false );

//cycle dependent -- cycle defaults
	NEW_OPT7( noesy_weights::defaults::Vmin, "acceptable minimial volume contribution per assignment", 0.01, 0.001, 0.005, 0.01, 0.025, 0.05, 0.501 );
	NEW_OPT7( noesy_weights::defaults::symmetry, "contribution of symmetry compliance to peak volume", 10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0 );
  NEW_OPT7( noesy_weights::defaults::covalent, "contribution of local covalent compliance to peak volume", 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 1.0 );
  NEW_OPT7( noesy_weights::defaults::decoys, "exponent controlling contribution of decoy compliance to peak volume", 3, 3, 6, 6, 6, 6, 6 );
	NEW_OPT7( noesy_weights::defaults::Smax,"maximum cumulative contribution of symmetry, covalent and network to peak volume", 20, 20, 20, 20, 10, 10, 10 );
  NEW_OPT7( noesy_weights::defaults::dcut, "upper limit on acceptable distance violation for elimination of spurious NOESY cross peaks (A)", -1, 1.5, 0.9, 0.6, 0.3, 0.1, 0.1 );
	NEW_OPT7( noesy_weights::defaults::network_min, "Threshold for acceptable lower limit of network-anchoring per residue", 1.0, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5 );
  NEW_OPT7( noesy_weights::defaults::network_atom_min, "contribution of network anchoring to peak volume", 0.25, 0.25, 0.4, 0.4, 0.4, 0.4, 0.4 );
  NEW_OPT7( noesy_weights::defaults::dcalibrate, "upper limit on acceptable distance violation for structure dependent calibration (A)", -1, 1.5, 0.9, 0.6, 0.3, 0.1, 0.1 );
  NEW_OPT7( noesy_weights::defaults::calibration_target, "target for NOE calibration > 1 -- (A) (structure independent ) <1 (%) of models violated", 5, 0.15, 0.15, 0.1, 0.1, 0.1, 0.1 );


	//cycle dependent
	NEW_OPT( noesy_weights::Vmin, "acceptable minimial volume contribution per assignment", 0.01 );
	NEW_OPT( noesy_weights::symmetry, "contribution of symmetry compliance to peak volume", 10.0 );
  NEW_OPT( noesy_weights::covalent, "contribution of local covalent compliance to peak volume", 10.0 );
  NEW_OPT( noesy_weights::decoys, "exponent controlling contribution of decoy compliance to peak volume", 3 );
	NEW_OPT( noesy_weights::Smax,"maximum cumulative contribution of symmetry, covalent and network to peak volume", 20 );
  NEW_OPT( noesy_weights::dcut, "upper limit on acceptable distance violation for elimination of spurious NOESY cross peaks (A)", -1 );

	NEW_OPT( noesy_weights::network_min, "Threshold for acceptable lower limit of network-anchoring per residue", 1.0 );
  NEW_OPT( noesy_weights::network_atom_min, "contribution of network anchoring to peak volume", 0.25 );
  NEW_OPT( noesy_weights::dcalibrate, "upper limit on acceptable distance violation for structure dependent calibration (A)", -1 );
  NEW_OPT( noesy_weights::calibration_target, "target for NOE calibration > 1 -- (A) (structure independent ) <1 (%) of models violated", 5 );
	NEW_OPT( noesy::atom_dependent_calibration, "individual calibration constants per atom-group: backbone, side-chain, methyl", false );



	NEW_OPT( noesy_weights::elim_dist_viol, "percentage of decoys that can be violated by distance constraints before it is eliminated ->set to 1 to switch this feature off", 0.5 );


	NEW_OPT( noesy_weights::cst_strength, "curvature of the constraint potential after upper distance is violated", 4 );
	NEW_OPT( noesy_weights::cycle, "set to cycle 1-7 to select a set of cycle-dependent options", 1 );
	NEW_OPT( noesy_weights::centroid_padding, "if sidechain NOESY are mapped to CB add some padding", 0.5 );

  NEW_OPT( noesy::no_network, "skip network analysis", false );
	NEW_OPT( noesy::network::include_reverse_dir, "scan a->y->b and also b->y->a noesy connections to validate a->b", true );
	NEW_OPT( noesy::network::allow_same_residue_connect, "scan also a->y->b', a'->y->b and a'->y-b' with a' and b' in same residues as a and b, respectively", true);
	NEW_OPT( noesy::network::use_all_covalent_atoms, "resonance list contains atoms (e.g., CA, CB) that are never part of an assignment (e.g., H, QD1, etc.) should CA, CB be included in network analysis? ", false );

	NEW_OPT( noesy::use_local_distviol, "don't use global distance-violation parameter dcut, but make-up local one using the distribution of distances", false );
	NEW_OPT( noesy::local_distviol::range, "get distance-difference between lowest X% and hightest X% of decoys", 0.00 );
	NEW_OPT( noesy::local_distviol::global_buffer, "add an extra grace-buffer (A) to the local range before you evaluate dist-cutoff", 0.5 );



  options_registered_ = true;
}


static basic::Tracer tr("protocols.noesy_assign.parameters");


void protocols::noesy_assign::PeakAssignmentParameters::set_options_from_cmdline() {
  using namespace basic::options::OptionKeys;
  using namespace basic::options;

	runtime_assert( options_registered_ );

	//set cycle from cmdline if not specified
	if ( cycle_selector_ == 0 ) cycle_selector_ = option[ noesy_weights::cycle ]();

	//cycle independent
  chemshift_overlap_weight_ = option[ noesy_weights::chemshift ](); //Gamma, eq. (4)
  ignore_resonancefile_tolerances_ = option[ noesy::ignore_resonancefile_tolerances ]();
	//  dmax_ = 5.5;
  vmin_ = 0.1; //minimum peak-volume contribution to network anchoring
  vmax_ = 1.0; //maximum peak-volume contribution to network anchoring
  nmax_ = 20; //maximum number of assignments
  nr_conformers_violatable_ = option[ noesy_weights::elim_dist_viol ](); // M/2
  network_reswise_high_ = 4.0; //used in network-based elimination
	centroid_mapping_distance_padding_ = option[ noesy_weights::centroid_padding ]();

	//cycle dependent
#define opt_read_out_macro( VAR, OPT )															\
{ VAR = option[ noesy_weights::defaults::OPT ]()[ cycle_selector_ ]; \
	if ( option [ noesy_weights::OPT ].user() ) {  \
		VAR = option[ noesy_weights::OPT ](); \
	}\
}\


	//read cycle-dependent default value -- overwrite if option is selected
	opt_read_out_macro( min_volume_, Vmin );
	opt_read_out_macro( symmetry_compliance_weight_, symmetry );
	opt_read_out_macro( covalent_compliance_weight_, covalent );
	opt_read_out_macro( decoy_compatibility_exponent_, decoys );
	opt_read_out_macro( smax_, Smax );
	opt_read_out_macro( dcut_, dcut );
	opt_read_out_macro( dcalibrate_, dcalibrate );
 	opt_read_out_macro( network_reswise_min_, network_min );
	opt_read_out_macro( network_atom_min_, network_atom_min );
	opt_read_out_macro( calibration_target_, calibration_target );

	atom_dependent_calibration_ = option[ noesy::atom_dependent_calibration ];

	no_network_ = option[ noesy::no_network ]();
	network_include_reverse_dir_= option[ noesy::network::include_reverse_dir ];
  network_allow_same_residue_connect_ = option[ noesy::network::allow_same_residue_connect ];
	network_use_all_covalent_atoms_ = option[ noesy::network::use_all_covalent_atoms ];
	cst_strength_ = option[ noesy_weights::cst_strength ]();

	use_local_distviol_ = option[ noesy::use_local_distviol ]();
	local_distviol_range_ = option[ noesy::local_distviol::range ]();
	local_distviol_global_buffer_ = option[ noesy::local_distviol::global_buffer ]();
	if ( local_distviol_range_ >= 0.5 ) {
		utility_exit_with_message( "local_distviol::range must be in the range 0.0...0.5" );
	}
	show_on_tracer();
}

void protocols::noesy_assign::PeakAssignmentParameters::show( std::ostream& os ) const {
	os << " ============== cycle-dependent parameter selection for Noesy peak assignment =============== " << std::endl;
	os << "   cycle: " << cycle_selector_ << std::endl;
	os << "----------------------------------------------------------------------------------------------" << std::endl;
	os << "        Vmin: " << min_volume_ << std::endl;
	os << "        T:    " << symmetry_compliance_weight_ << std::endl;
	os << "        O:    " << covalent_compliance_weight_ << std::endl;
	os << "        eta:  " << decoy_compatibility_exponent_ << std::endl;
	os << "        Smax: " << smax_ << std::endl;
	os << "        dcut: " << dcut_ << std::endl;
	os << "   Nmin(res): " << network_reswise_min_ << std::endl;
	os << "  Nmin(atom): " << network_atom_min_ << std::endl;
	os << " Calibration Target: " << calibration_target_ << ( calibration_target_ >= 1 ? " (A) average distance " : " fract. of constraints violate input structures " ) << std::endl;
	os << " Calibration distance: " << dcalibrate_ << " A above upper bound is a violation" << std::endl;
	os << "==============================================================================================" << std::endl;

	os << " ============== cycle-independent parameter for Noesy peak assignment =============== " << std::endl;
	os << "      chemical shift agreement-weight: " << chemshift_overlap_weight_ << std::endl;
	os << "      elimination of peak if distance-bound violates more than: " << nr_conformers_violatable_ << std::endl;
	os << "      elimination MAXASSIGN if assignments count is above: " << nmax_ << std::endl;
	//	os << "      contribution of network anchoring to peak volume XXX" << std::endl;
	os << "      strength (curvature) of the constraint potential after upper distance is violated " << cst_strength_ << std::endl;
	os << "      padding for CB-mapped protons in centroid mode: " << centroid_mapping_distance_padding_ << std::endl;
	os << "      atom_dependent_calibration: " << (atom_dependent_calibration_ ? "yes" : "no" ) << std::endl;
	os << "----------------------------------------------------------------------------------------------" << std::endl;
	if ( !no_network_ ) {
		os << "      network-anchoring options: " << std::endl;
		os << "                minimum contribution to peak-volume " << vmin_ << std::endl;
		os << "                maximum contribution to peak-volume " << vmax_ << std::endl;
		os << "                do not eliminate if network score per residue is above " << network_reswise_high_ << std::endl;
		os << "                network_include_reverse_dir:           " << (network_include_reverse_dir_ ? "yes ": "no " ) << std::endl;
		os << "                network_allow_same_residue_connection: " << (network_allow_same_residue_connect_ ? "yes ": "no " ) << std::endl;
		os << "                network_use_all_covalent_atoms:        " << (network_use_all_covalent_atoms_ ? "yes" : "no" ) << std::endl;
	} else {
		os << "      no network-anchoring applied! " << std::endl;
	}
	if ( use_local_distviol_ ) {
		os << "     use local distance violation " << std::endl;
		os << "                local distance violation range " << local_distviol_range_*100 << std::endl;
		os << "                add the following constant distance " << local_distviol_global_buffer_ << std::endl;
	}
	os << "==============================================================================================" << std::endl;
}

void protocols::noesy_assign::PeakAssignmentParameters::show_on_tracer() const {
	show( tr.Info );
}


namespace protocols {
namespace noesy_assign {
using namespace core;

}
}
