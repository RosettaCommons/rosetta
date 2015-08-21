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
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#define DEFINE_OPTIONS_NOW
#include <protocols/noesy_assign/PeakAssignmentOptionKeys.hh>
#undef DEFINE_OPTIONS_NOW

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

#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace noesy_assign {

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex PeakAssignmentParameters::singleton_mutex_;

std::mutex & PeakAssignmentParameters::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
PeakAssignmentParameters const * PeakAssignmentParameters::get_instance()
{
	boost::function< PeakAssignmentParameters * () > creator = boost::bind( &PeakAssignmentParameters::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

PeakAssignmentParameters *
PeakAssignmentParameters::create_singleton_instance()
{
	PeakAssignmentParameters * inst = new PeakAssignmentParameters;
	inst->set_options_from_cmdline();
	return inst;
}

bool PeakAssignmentParameters::options_registered_( false );

#if defined MULTI_THREADED && defined CXX11
std::atomic< PeakAssignmentParameters * > PeakAssignmentParameters::instance_( 0 );
#else
PeakAssignmentParameters * PeakAssignmentParameters::instance_( 0 );
#endif

/// @details DANGER DANGER DANGER! NOT IN ANY WAY THREADSAFE!!!
void
PeakAssignmentParameters::reset() {
	if ( instance_ ) delete instance_;
	instance_ = NULL;
}

/// @details DANGER DANGER DANGER! NOT IN ANY WAY THREADSAFE!!!
PeakAssignmentParameters*
PeakAssignmentParameters::get_nonconst_instance() {
	if ( instance_ ) return instance_;
	instance_ = new PeakAssignmentParameters;
#if defined CXX11 && defined MULTI_THREADED
	instance_.load()->set_options_from_cmdline();
#else
  instance_->set_options_from_cmdline();
#endif
	return instance_;
}

/// @details DANGER DANGER DNAGER: Completely thread unsafe.
void PeakAssignmentParameters::set_cycle( core::Size cycle ) {
	if ( instance_ ) delete instance_;
	instance_ = new PeakAssignmentParameters;
#if defined CXX11 && defined MULTI_THREADED
	instance_.load()->set_options_from_cmdline( cycle );
#else
  instance_->set_options_from_cmdline( cycle );
#endif
}

void PeakAssignmentParameters::register_options() {
	using namespace basic::options;
	using namespace OptionKeys;
	if ( options_registered_ ) return;

	//cycle independent
	NEW_OPT( noesy_weights::chemshift, "contribution of chem. shift overlap to peak volume", 0.5 );
	// NEW_OPT( noesy_weights::network_high, "contribution of network anchoring to peak volume", 10.0 );
	NEW_OPT( noesy::ignore_resonancefile_tolerances, "ignore the tolerances in the resonance file", false );
	NEW_OPT( noesy::ignore_resonancefile_intensities, "ignore the tolerances in the resonance file", false );

	//cycle dependent -- cycle defaults
	NEW_OPT7( noesy_weights::defaults::Vmin, "acceptable minimial volume contribution per assignment", 0.01, 0.001, 0.005, 0.01, 0.025, 0.05, 0.501 );
	NEW_OPT7( noesy_weights::defaults::symmetry, "contribution of symmetry compliance to peak volume", 10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0 );
	NEW_OPT7( noesy_weights::defaults::covalent, "contribution of local covalent compliance to peak volume", 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 1.0 );
	//  NEW_OPT7( noesy_weights::defaults::decoys, "exponent controlling contribution of decoy compliance to peak volume", 3, 3, 6, 6, 6, 6, 6 );
	NEW_OPT7( noesy_weights::defaults::Smax,"maximum cumulative contribution of symmetry, covalent and network to peak volume", 20, 20, 20, 20, 10, 10, 10 );
	NEW_OPT7( noesy_weights::defaults::dcut, "upper limit on acceptable distance violation for elimination of spurious NOESY cross peaks (A)", -1, 1.5, 0.9, 0.6, 0.3, 0.1, 0.1 );
	// NEW_OPT7( noesy_weights::defaults::reswise_min, "Threshold for acceptable lower limit of network-anchoring per residue", 1.0, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5 );
	//  NEW_OPT7( noesy_weights::defaults::atomwise_min, "contribution of network anchoring to peak volume", 0.25, 0.25, 0.4, 0.4, 0.4, 0.4, 0.4 );
	NEW_OPT7( noesy_weights::defaults::dcalibrate, "upper limit on acceptable distance violation for structure dependent calibration (A)", -1, 1.5, 0.9, 0.6, 0.3, 0.1, 0.1 );
	NEW_OPT7( noesy_weights::defaults::calibration_target, "target for NOE calibration > 1 -- (A) (structure independent ) <1 (%) of models violated", 3.8, 0.15, 0.15, 0.1, 0.1, 0.1, 0.1 );

	//cycle dependent

	NEW_OPT( noesy_weights::symmetry, "contribution of symmetry compliance to peak volume", 10.0 );
	NEW_OPT( noesy_weights::covalent, "contribution of local covalent compliance to peak volume", 10.0 );
	//  NEW_OPT( noesy_weights::decoys, "exponent controlling contribution of decoy compliance to peak volume", 3 );
	NEW_OPT( noesy_weights::Smax,"maximum cumulative contribution of symmetry, covalent and network to peak volume", 20 );


	// NEW_OPT( noesy_weights::network_min, "Threshold for acceptable lower limit of network-anchoring per residue", 1.0 );
	//NEW_OPT( noesy_weights::network_atom_min, "contribution of network anchoring to peak volume", 0.25 );
	NEW_OPT( noesy_weights::dcalibrate, "upper limit on acceptable distance violation for structure dependent calibration (A)", -1 );
	NEW_OPT( noesy_weights::calibration_target, "target for NOE calibration > 1 -- (A) (structure independent ) <1 (%) of models violated", 5 );
	NEW_OPT( noesy::atom_dependent_calibration, "individual calibration constants per atom-group: backbone, side-chain, methyl", false );

	NEW_OPT( noesy::elim::max_assign, "maximum number of assignments to a peak before it is eliminated", 20 );
	NEW_OPT( noesy::elim::vmin, "acceptable minimial volume contribution per assignment", 0.01 );
	NEW_OPT( noesy::elim::dist_viol, "percentage of decoys that can be violated by distance constraints before it is eliminated ->set to 1 to switch this feature off", 0.5 );
	NEW_OPT( noesy::elim::dcut, "upper limit on acceptable distance violation for elimination of spurious NOESY cross peaks (A)", -1 );

	NEW_OPT( noesy::map_to_cen_atom,"map the centroid restraints to CEN atom", false );

	NEW_OPT( noesy_weights::cst_strength, "curvature of the constraint potential after upper distance is violated", 4 );
	NEW_OPT( noesy_weights::cycle, "set to cycle 1-7 to select a set of cycle-dependent options", 1 );
	NEW_OPT( noesy_weights::centroid_padding, "if sidechain NOESY are mapped to CB add some padding", 0.5 );
	NEW_OPT( noesy_weights::min_symmetry_reinforcement, "minimum contribution to symmetry score from a symmetric peak", 0.0 );
	NEW_OPT( noesy::no_network, "skip network analysis", false );
	NEW_OPT( noesy::network::reswise_min, "Threshold for acceptable lower limit of network-anchoring per residue",  1.0 );
	NEW_OPT( noesy::network::reswise_high, "Threshold for always accepted (disregard atomwise_min) network anchoring per residue", 4.0 );
	NEW_OPT( noesy::network::atomwise_min, "Threshold for acceptable network-anchoring per assignment", 0.25 );
	NEW_OPT( noesy::network::include_reverse_dir, "scan a->y->b and also b->y->a noesy connections to validate a->b", true );
	NEW_OPT( noesy::network::allow_same_residue_connect, "scan also a->y->b', a'->y->b and a'->y-b' with a' and b' in same residues as a and b, respectively", true);
	NEW_OPT( noesy::network::use_all_covalent_atoms, "resonance list contains atoms (e.g., CA, CB) that are never part of an assignment (e.g., H, QD1, etc.) should CA, CB be included in network analysis? ", false );
	NEW_OPT( noesy::network::mode, "choose implementation of network anchoring to be used [orig,clean]", "clean" );
	NEW_OPT( noesy::network::vmin, "minimum contribution to peak-volume to be considered for network anchoring", 0.1 );
	NEW_OPT( noesy::network::vmax, "minimum contribution to peak-volume to be considered for network anchoring", 1.0 );
	NEW_OPT( noesy::use_local_distviol, "don't use global distance-violation parameter dcut, but make-up local one using the distribution of distances", false );
	NEW_OPT( noesy::local_distviol::range, "get distance-difference between lowest X% and hightest X% of decoys", 0.99 );
	NEW_OPT( noesy::local_distviol::global_buffer, "add an extra grace-buffer (A) to the local range before you evaluate dist-cutoff", 0.5 );
	NEW_OPT( noesy::local_distviol::global_factor, "allow violation of X*interval beofre eliminated", 1.0 );
	NEW_OPT( noesy::local_distviol::cutoff, "peaks with a low-quartil distance of more than X are always removed", 8.0 );
	NEW_OPT( noesy::local_distviol::cutoff_buffer, "peaks with a low-quartil distance of more than X+distance_bound are always removed", 2.0 );

	NEW_OPT( noesy::calibration::convergence, "use only distance with stddev (over ensemble) of less than X Angstrom for calibration (0=all)", 0 );
	NEW_OPT( noesy::calibration::max_noe_dist, "use max_noe_dist to cap upper distance bounds, this value is overwritten by content in file", 0 );
	NEW_OPT( noesy::calibration::max_nudging, "individual restraint can be extended by maximum of X% of original calibrated distance", 1.0 );
	NEW_OPT( noesy::calibration::start_nudging,"start nudging if X% or more conformers violate original calibrated distance", 0.1 );
	NEW_OPT( noesy::calibration::stop_nudging,"stop nudging if only X% of conformers violated original calibrated distance", 0.1 );
	NEW_OPT( noesy::calibration::eliminate, "enable elimination by distance violation and nudging even when dcalibrate -1", false );
	NEW_OPT( noesy::calibration::use_median, "use median and not average for calibration", false );
	NEW_OPT( noesy::calibration::cycles, "how many cycles of calibration and elimination by distance-violation", 1 );
	NEW_OPT( noesy::calibration::ignore_eliminated_peaks, "ignore peaks that are already eliminated", false );

	NEW_OPT5( noesy::prob::sigmoid::tau, "tau and m parameter for sigmoids in crosspeak-probability score", 0.3, 1, .2,  0.2, 2);
	NEW_OPT5( noesy::prob::sigmoid::m, "tau and m parameter for sigmoids in crosspeak-probability score", 0.5, 2, .4, 0.5, 4);
	NEW_OPT2( noesy::prob::level, "selection levels for HI, MED, LOW probability cross-peaks", 0.7, 0.45 );
	NEW_OPT5( noesy::prob::sigmoid::w, "sigmoid contribution", 0.2000, 0.6000, 0.2000, 0.4000, 0.2000 );
	options_registered_ = true;
}


static thread_local basic::Tracer tr( "protocols.noesy_assign.parameters" );


void PeakAssignmentParameters::set_options_from_cmdline( core::Size cycle_selector ) {
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

	runtime_assert( options_registered_ );

	//set cycle from cmdline if not specified
	//if ( cycle_selector_ == 0 )
	if ( !cycle_selector ) {
		cycle_selector_ = option[ noesy_weights::cycle ]();
	} else {
		cycle_selector_ = cycle_selector;
	}

	//cycle independent
	chemshift_overlap_weight_ = option[ noesy_weights::chemshift ](); //Gamma, eq. (4)
	ignore_resonancefile_tolerances_ = option[ noesy::ignore_resonancefile_tolerances ]();
	ignore_resonancefile_intensities_ = option[ noesy::ignore_resonancefile_intensities ]();
	//  dmax_ = 5.5;
	vmin_ = option[ noesy::network::vmin ](); //previously 0.1 minimum peak-volume contribution to network anchoring
	vmax_ = option[ noesy::network::vmax ](); //previously 1.0; //maximum peak-volume contribution to network anchoring
	nmax_ = option[ noesy::elim::max_assign ](); //20 maximum number of assignments
	nr_conformers_violatable_ = option[ noesy::elim::dist_viol ](); // M/2
	network_reswise_min_ = option[ noesy::network::reswise_min ]();
	network_atom_min_ = option[ noesy::network::atomwise_min ]();
	// opt_read_out_macro( network_reswise_min_, network_min );
	// opt_read_out_macro( network_atom_min_, network_atom_min );
	network_reswise_high_ = option[ noesy::network::reswise_high ]();// 4.0; //used in network-based elimination
	centroid_mapping_distance_padding_ = option[ noesy_weights::centroid_padding ]();

	//cycle dependent
	#define opt_read_out_macro( VAR, OPT )               \
{ VAR = option[ noesy_weights::defaults::OPT ]()[ cycle_selector_ ]; \
	if ( option [ noesy_weights::OPT ].user() ) {  \
		VAR = option[ noesy_weights::OPT ](); \
	}\
}\

	#define opt_read_out_macro2( VAR, OPT, DIROPT )               \
{ VAR = option[ noesy_weights::defaults::OPT ]()[ cycle_selector_ ]; \
	if ( option [ noesy::DIROPT ].user() ) {  \
		VAR = option[ noesy::DIROPT ](); \
	}\
}\


	//read cycle-dependent default value -- overwrite if option is selected
	opt_read_out_macro2( min_volume_, Vmin, elim::vmin );
	opt_read_out_macro( symmetry_compliance_weight_, symmetry );
	opt_read_out_macro( covalent_compliance_weight_, covalent );
	// opt_read_out_macro( decoy_compatibility_exponent_, decoys );
	opt_read_out_macro( smax_, Smax );
	opt_read_out_macro2( dcut_, dcut, elim::dcut );
	opt_read_out_macro( dcalibrate_, dcalibrate );

	opt_read_out_macro( calibration_target_, calibration_target );

	atom_dependent_calibration_ = option[ noesy::atom_dependent_calibration ];
	map_to_cen_atom_ = option[ noesy::map_to_cen_atom ];

	no_network_ = option[ noesy::no_network ]();
	network_include_reverse_dir_= option[ noesy::network::include_reverse_dir ];
	network_allow_same_residue_connect_ = option[ noesy::network::allow_same_residue_connect ];
	network_use_all_covalent_atoms_ = option[ noesy::network::use_all_covalent_atoms ];
	network_mode_ = option[ noesy::network::mode ]();
	cst_strength_ = option[ noesy_weights::cst_strength ]();
	min_contribution_symmetric_peaks_ = option[ noesy_weights::min_symmetry_reinforcement ];
	use_local_distviol_ = option[ noesy::use_local_distviol ]();
	local_distviol_range_ = option[ noesy::local_distviol::range ]();
	local_distviol_global_buffer_ = option[ noesy::local_distviol::global_buffer ]();
	local_distviol_global_factor_ = option[ noesy::local_distviol::global_factor ]();
	local_distviol_cutoff_ = option[ noesy::local_distviol::cutoff ]();
	local_distviol_cutoff_buffer_ = option[ noesy::local_distviol::cutoff_buffer ]();

	calibration_convergence_ = option[ noesy::calibration::convergence ]();
	calibration_max_noe_dist_ = option[ noesy::calibration::max_noe_dist ]();
	calibration_stop_nudging_ = option[ noesy::calibration::stop_nudging ]();
	calibration_start_nudging_ = option[ noesy::calibration::start_nudging ]();
	calibration_max_nudging_ = option[ noesy::calibration::max_nudging ]();
	calibration_eliminate_ = option[ noesy::calibration::eliminate ]();
	calibration_use_median_ = option[ noesy::calibration::use_median ]();
	calibration_ignore_eliminated_peaks_ = option[ noesy::calibration::ignore_eliminated_peaks ]();
	calibration_cycles_ = option[ noesy::calibration::cycles ]();

	if ( local_distviol_range_ <= 0.5 ) {
		utility_exit_with_message( "local_distviol::range must be in the range 0.5...1.0" );
	}

	prob_level_ = option[ noesy::prob::level ]();
	prob_sigmoid_tau_ = option[ noesy::prob::sigmoid::tau ]();
	prob_sigmoid_m_ = option[ noesy::prob::sigmoid::m ]();
	prob_sigmoid_w_ = option[ noesy::prob::sigmoid::w ]();
	show_on_tracer();
}

void PeakAssignmentParameters::show( std::ostream& os ) const {
	os << " ============== cycle-dependent parameter selection for Noesy peak assignment =============== " << std::endl;
	os << "   cycle: " << cycle_selector_ << std::endl;
	os << "----------------------------------------------------------------------------------------------" << std::endl;
	os << "        Vmin: " << min_volume_ << std::endl;
	os << "        T:    " << symmetry_compliance_weight_ << std::endl;
	os << "        V:    " << covalent_compliance_weight_ << std::endl;
	// os << "        eta:  " << decoy_compatibility_exponent_ << std::endl;
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
	// os << "      contribution of network anchoring to peak volume XXX" << std::endl;
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
		os << "                local distance violation range " << local_distviol_range_*100 << "%" << std::endl;
		os << "                add the following constant distance " << local_distviol_global_buffer_ << std::endl;
		os << "                multiply max-extension by " << local_distviol_global_factor_ << std::endl;
	}
	os << "==============================================================================================" << std::endl;
}

void PeakAssignmentParameters::show_on_tracer() const {
	show( tr.Info );
}


}
}

namespace protocols {
namespace noesy_assign {

using namespace core; //???

}
}
