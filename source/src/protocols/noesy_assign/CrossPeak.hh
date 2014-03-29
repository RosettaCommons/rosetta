// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_CrossPeak_HH
#define INCLUDED_protocols_noesy_assign_CrossPeak_HH


// Unit Header
#include <protocols/noesy_assign/CrossPeak.fwd.hh>

// Package Headers
#include <protocols/noesy_assign/CrossPeakInfo.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakAssignment.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>
//#include <list>
//#include <map>
#include <set>

#ifndef WIN32
#include <protocols/noesy_assign/PeakAssignment.fwd.hh>
#endif

namespace protocols {
namespace noesy_assign {

class CrossPeak : public utility::pointer::ReferenceCount {
public:
  ///@brief a single dimension of a multidimensional peak, stores putative assignments as indices into Resonance table
  class Spin {
  public:
    Spin( core::Real freq );
    Spin();
    ~Spin();

    ///@brief resonance frequency of peak
    core::Real freq() const { return freq_; }
    void set_freq( core::Real freq ) { freq_ = freq; }

    ///@brief add index-pointer into Resonance table to assign peak to putative spin
    void add_assignment( core::Size assignment ) {
      assignments_.push_back( assignment );
    }

    ///@brief return index-pointer into Resonance table
    core::Size assignment( core::Size nr ) const {
      return assignments_[ nr ];
    }

    core::Size n_assigned() const {
      return assignments_.size();
    }

    core::Size assignment_index( core::Size assignment ) const; //return 0 if not found

  private:
    core::Real freq_;
    typedef utility::vector1< core::Size > SpinAssignments;
    SpinAssignments assignments_; //index into resonance table...
  };

public:
#ifndef WIN32
  typedef utility::vector1< PeakAssignmentOP > PeakAssignments;
  typedef PeakAssignments::const_iterator const_iterator;
  typedef PeakAssignments::iterator iterator;
#endif
  enum EliminationReason {
    NOT_ELIMINATED = 0,
    EL_DISTVIOL,
    EL_NETWORK,
    EL_MINPEAKVOL,
    EL_MAXASSIGN
  };

  enum QualityClass {
    HI_UNAMBIG = 0,
    HI_NEAR_UNAMBIG,
    HI_AMBIG,
    MED_AMBIG,
    UNAMBIG_MED_PROB,
    BAD_LOW_PROB,
    MAX_CLASS
  };

  CrossPeak( Spin const&, Spin const&, core::Real strength );
  CrossPeak();
  virtual ~CrossPeak();
  virtual CrossPeakOP empty_clone() {
    return new CrossPeak;
  }

  virtual bool has_label( core::Size ) const { return false; }
  virtual Spin& label( core::Size ) { runtime_assert( false ); return proton1_; }
  virtual Spin const& label( core::Size ) const { runtime_assert( false ); return proton2_; }

  Spin const& proton( core::Size i ) const { return i>=2 ? proton2_ : proton1_ ; }
  Spin& proton( core::Size i ) { return i>=2 ? proton2_ : proton1_ ; }
  bool has_proton( core::Size ) const;

  ///@brief flat access to spins in [ a, b, h(a), h(b) ] order
  Spin const& spin( core::Size i ) const { return i>2 ? label( i-2 ) : proton( i ); }
  Spin& spin( core::Size i )  { return i>2 ? label( i-2 ) : proton( i ); }

  core::Size dimension() const { return has_label( 1 ) ? ( has_label( 2 ) ? 4 : 3 ) : 2; }
  ///@brief expect res_ids in order: spin1, spin2, label1, label2
  virtual void add_full_assignment( core::Size res_ids[] );

  FoldResonance const& folder( core::Size i ) {
    if ( i==1 ) {
      return info1_->fold_proton_resonance();
    } else if ( i==2 ) {
      return info2_->fold_proton_resonance();
    } else if ( i==3 ) {
      return info1_->fold_label_resonance();
    } else if ( i==4 ) {
      return info2_->fold_label_resonance();
    };
    runtime_assert( false );
    return UNFOLDED_; //to make compiler happy
  }

  core::Real fold_resonance( core::Real freq, core::Size i /*dim*/ ) {
    // dim = 1, 2 are the protons
    // dim = 3, 4 are the labels

    /// for now only labels can be folded
    if ( i==1 ) {
      return info1_->fold_proton_resonance()( freq );
    } else if ( i==2 ) {
      return info2_->fold_proton_resonance()( freq );
    } else if ( i==3 ) {
      return info1_->fold_label_resonance()( freq );
    } else if ( i==4 ) {
      return info2_->fold_label_resonance()( freq );
    };
    runtime_assert( false ); //should never get here...
    return 0.0; //to make compiler happy
  }
//  virtual void read_from_stream( std::istream& );
//   virtual void add_assignment_from_stream( std::istream& );
//   virtual void write_to_stream( std::ostream& ) const;

  virtual void find_assignments();
  virtual void assign_spin( core::Size spin_id );

  //provide array with spin 1 spin 2 label 1 label 2 -- returns index in assignment list
  virtual core::Size assign_spin( core::Size spin_id, core::Size res_ids[] );

  ///@brief the measured integral of the peak
  core::Real volume() const { return volume_; }
  void set_volume( core::Real val ) { volume_ = val; }

  void set_resonances( ResonanceListOP res_in ) {
    resonances_=res_in;
  }

  ResonanceList const& resonances() const { return *resonances_; }

  ///@brief the cumulative cyana-type weights for all assignments Vk
  core::Real cumulative_peak_volume() const {
    return cumulative_peak_volume_;
  }

  void set_cumulative_peak_volume( core::Real setting ) {
    cumulative_peak_volume_ = setting;
  }

  ///@brief the largest volume contribution (normalized) any of the assignments has
  core::Real max_volume_contribution() const;

  core::Real probability() const;

  void set_eliminated_due_to_dist_violations( bool setting ) {
    if ( eliminated_due_to_dist_violations_ && !setting ) eliminated_ = NOT_ELIMINATED; //remove elimination if this is not longer violated
    eliminated_due_to_dist_violations_ = setting;
    eliminated_ = ( setting ) ? EL_DISTVIOL : eliminated_;
  }

  void set_elimination_candidate( bool setting = true ) {
    elimination_candidate_ = setting;
  }

  void set_elimination_comment( std::string const& str ) {
    elimination_comment_ = str;
  }

  bool is_elimination_candidate() {
    return elimination_candidate_;
  }

  ///@brief returns true if this peak is to be ignored due to points (i)-(iv) on p215 of JMB 2002, 319,209-227
  /// do_not_compute ... for outputter that does not want to change state...
  bool eliminated( bool recompute = false, bool do_not_compute = false ) const;
  std::string elimination_reason() const;

  core::Size min_seq_separation_residue_assignment( core::Real volume_threshold ) const;

  void
  create_fa_and_cen_constraint(
       core::scoring::constraints::ConstraintOP& fa_cst,
       core::scoring::constraints::ConstraintOP& cen_cst,
       core::pose::Pose const& pose,
       core::pose::Pose const& centroid_pose,
       core::Size normalization,
       core::Real padding = 0.0,
       bool fa_only = false
  ) const;

#if 0
  core::scoring::constraints::ConstraintOP
  create_constraint(
    core::pose::Pose const& pose,
    core::Size normalization = 1
  ) const;

  core::scoring::constraints::ConstraintOP
  create_centroid_constraint(
    core::pose::Pose const& pose,
    core::pose::Pose const& centroid_pose,
    core::Size normalization = 1
  ) const;
#endif
//   ///@brief number of assigned protons
//   core::Size n_assigned() const { runtime_assert( proton1_.n_assigned() == proton2_.n_assigned() ); return proton1_.n_assigned(); }

//   ///@brief number of valid assignments
//   core::Size n_valid_assigned() const { return proton1_.n_valid_assigned(); }

//   ///@brief return the index for the next valid assigned spin.  -- return 0 for last
//   core::Size next_assigned( core::Size last = 0 ) const {
//     return proton1_.next_assigned( last );
//   }

#ifndef WIN32
  bool assigned() const { return assignments_.size(); }
  bool ambiguous() const { return assignments_.size() > 1; }
  core::Size n_assigned() const { return assignments_.size(); }
#endif
  ///@brief number of assignments with peak volume higher than params.min_volume_
  core::Size n_Vmin_assignments();

  //void invalidate_assignment( Size assignment );

#ifndef WIN32
  PeakAssignments const& assignments() const { return assignments_; }
  const_iterator begin() const { return assignments_.begin(); }
  iterator begin() { return assignments_.begin(); }
  const_iterator end() const { return assignments_.end(); }
#endif

  core::Real distance_bound() const { return distance_bound_; }
  void nudge_distance_bound( core::Real offset );
  core::Size peak_id() const { return peak_id_; };

  std::string const& filename() const { return info1_->filename(); }
  std::size_t exp_hash() const { return info1_->exp_hash(); }

  bool same_peak( CrossPeak const & other ) const {
    return peak_id() == other.peak_id() && exp_hash() == other.exp_hash();
  }

  void set_peak_id( core::Size val ) {
    peak_id_ = val;
  }

  //  CrossPeakInfo const& info_struct( core::Size i ) const { return i>=2 ? *info2_ : *info1_; }
  CrossPeakInfo const& info( core::Size i ) const { return i>=2 ? *info2_ : *info1_; }

  void set_info( core::Size i, CrossPeakInfoCOP info  ) { if ( i<2 ) info1_=info; else info2_=info; }

  void calibrate( PeakCalibrator const&, CalibrationTypeCumulator& calibration_types );

  core::Real tolerance( core::Size d ) {
    runtime_assert( info1_ && info2_ );
    if ( d == 1 ) {
      return info1_->proton_tolerance();
    } else if ( d == 2 ) {
      return info2_->proton_tolerance();
    } else if ( d == 3 ) {
      return info1_->label_tolerance();
    } else if ( d == 4 ) {
      return info1_->label_tolerance();
    }
    runtime_assert( false );
    return 0.0;
  }

  bool is4D() const {
    return has_label( 1 ) && has_label( 2 ) && (info1_->label_tolerance() < 10) && (info1_->proton_tolerance() <10)
      && (info2_->label_tolerance() <10) && (info2_->proton_tolerance() <10);
  }

  void print_peak_info( std::ostream& ) const;
  //  core::Size assignment( core::Size nr ) const;
  QualityClass quality_class() const;
  std::string quality_class_str() const;
  core::Real smallest_native_violation() const;
private:

  FoldResonance UNFOLDED_; //dummy folder
#ifndef WIN32
  PeakAssignments assignments_;
#endif
  ResonanceListOP resonances_;
  CrossPeakInfoCOP info1_, info2_;
  Spin proton1_;
  Spin proton2_;
  core::Real volume_; //exp intensity
  core::Size peak_id_;
  core::Real cumulative_peak_volume_; //sum Vk
  core::Real distance_bound_; //b -- computed after calibration...
  mutable EliminationReason eliminated_;
  bool eliminated_due_to_dist_violations_;//Mviol, dcut
  bool elimination_candidate_;
  std::string elimination_comment_;

protected:
  static std::set< core::id::NamedAtomID > unknown_resonances_;
};

inline std::ostream& operator<< ( std::ostream& os, CrossPeak const& peak ) {
  peak.print_peak_info( os );
  return os;
}

class CrossPeak3D : public CrossPeak {
public:
  CrossPeak3D( Spin const& sp1, Spin const& sp2, Spin const& label1, core::Real strength );
  CrossPeak3D();
  ~CrossPeak3D();
  virtual CrossPeakOP empty_clone() {
    return new CrossPeak3D;
  }
  virtual void assign_spin( Size spin_id );
  virtual core::Size assign_spin( core::Size spin_id, core::Size res_ids[] ); //provide array with spin 1 spin 2 label 1 label 2
  virtual void assign_labelled_spin( Size proton );
  virtual bool has_label( core::Size i ) const { return i==1; }
  virtual Spin& label( core::Size ) { return label1_; }
  virtual Spin const& label( core::Size ) const { return label1_; }
private:
  Spin label1_;
};

class CrossPeak4D : public CrossPeak3D {
public:
  CrossPeak4D( Spin const& sp1,  Spin const& sp2, Spin const& label1, Spin const& label2, core::Real strength );
  CrossPeak4D();
  ~CrossPeak4D();
  virtual CrossPeakOP empty_clone() {
    return new CrossPeak4D;
  }
  virtual void assign_spin( Size proton );
  virtual core::Size assign_spin( core::Size spin_id, core::Size res_ids[] ); //provide array with spin 1 spin 2 label 1 label 2
  virtual bool has_label( core::Size ) const { return true; }
  Spin& label( core::Size i ) { return i==1 ? CrossPeak3D::label( 1 ) : label2_; }
  Spin const& label( core::Size i ) const { return i==1 ? CrossPeak3D::label( 1 ) : label2_; }
private:
  Spin label2_;
};

}
}

#endif
