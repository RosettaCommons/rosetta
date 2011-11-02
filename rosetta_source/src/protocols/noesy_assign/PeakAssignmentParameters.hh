// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PeakAssignmentParametersList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_PeakAssignmentParameters_HH
#define INCLUDED_protocols_noesy_assign_PeakAssignmentParameters_HH


// Unit Header

// Package Headers

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.fwd.hh>

// Utility headers
// AUTO-REMOVED #include <utility/pointer/ReferenceCount.hh>

#include <iosfwd>


//// C++ headers

namespace protocols {
namespace noesy_assign {

class PeakAssignmentParameters { //: public utility::pointer::ReferenceCount {

private:
  static bool options_registered_;
  PeakAssignmentParameters() {}; //private constructor

public:
  static void register_options();
  void set_options_from_cmdline();
  static void set_cycle( core::Size );
  void show( std::ostream& ) const;
  void show_on_tracer() const;
  static PeakAssignmentParameters const* get_instance();
  static PeakAssignmentParameters * get_nonconst_instance();
private:
  static PeakAssignmentParameters* instance_;
  static core::Size cycle_selector_;

public:
  /* maybe make all options const
     make cmd-line options vectors for 7 cycles and have static function that
     switches cycles by throwing away one instance and creating a new one with the given
     cycle number and setting all const-params in the constructor from the cmd-line.

     all good, unless I want to change the things from somewhere else than cmdline...
  */

  //cycle independent
  core::Real chemshift_overlap_weight_; //Gamma, eq. (4)
  bool ignore_resonancefile_tolerances_;
  //  core::Real dmax_; //unused
  core::Real vmin_;
  core::Real vmax_;
  core::Real nmax_;
  core::Real nr_conformers_violatable_; //Mvio in fraction of nr_conformers
  core::Real network_reswise_high_;
  core::Real centroid_mapping_distance_padding_;
  //cycle dependent
  core::Real symmetry_compliance_weight_; //T, eq. (5)
  core::Real covalent_compliance_weight_; //O
  core::Real decoy_compatibility_exponent_; //eta, eq. (6)

  core::Real smax_;
  core::Real dcut_;
  core::Real dcalibrate_;

  bool use_local_distviol_;
  core::Real local_distviol_range_; //how many (in percent) decoys at both ends of the range are ignored to calculate max_extension
  core::Real local_distviol_global_buffer_;
  core::Real network_reswise_min_; //N_bar_min
  core::Real network_atom_min_; //N_min per atom
  core::Real calibration_target_;
  bool atom_dependent_calibration_;
  core::Real min_volume_; //minimum volume contribution

  core::Real cst_strength_; //for  1/cst_strength ->sigma for BoundFunc

  bool no_network_;
  bool network_use_all_covalent_atoms_;
  bool network_include_reverse_dir_;
  bool network_allow_same_residue_connect_;
};

}
}

#endif
