// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_PeakFileFormat_hh
#define INCLUDED_protocols_noesy_assign_PeakFileFormat_hh


// Unit Headers
#include <protocols/noesy_assign/CrossPeakInfo.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
//#include <core/id/NamedAtomID.fwd.hh>
//#include <core/chemical/AA.hh>

// Utility headers
// #include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
//#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers

#include <protocols/noesy_assign/CrossPeak.fwd.hh>
#include <protocols/noesy_assign/PeakAssignment.fwd.hh>
#include <protocols/noesy_assign/ResonanceList.fwd.hh>


namespace protocols {
namespace noesy_assign {

class PeakFileFormat_Base : public utility::pointer::ReferenceCount {
public:
  PeakFileFormat_Base()
    : ignore_assignments_( false ),
      min_VC_to_write_( 0.0 ),
      write_only_highest_VC_( false ),
      ignore_negative_intensity_( false ),
      minimum_peak_intensity_( 0.0 )
  {};

  virtual ~PeakFileFormat_Base() {
    if (info1_ && info2_ ) { //should this really be in the destructor ?
      info1_->set_filename( filename() );
      info2_->set_filename( filename() );
    }
  };

  virtual void write_peak( std::ostream&, core::Size ct, CrossPeak const& ) const = 0;
  virtual void write_resonances( std::ostream&, CrossPeak const& ) const = 0;
  virtual void write_strength( std::ostream&, CrossPeak const& ) const = 0;
  virtual void write_assignments( std::ostream&, CrossPeak const&, std::string const& first_line_end ) const = 0;

  virtual void read_resonances( std::istream&, CrossPeak& ) const = 0;
  virtual void read_assignments( std::istream& is,  std::istream& rest_line, CrossPeak&, std::string& new_peak_line ) const = 0;
  virtual void read_strength( std::istream&, CrossPeak& ) const = 0;

  void set_filename( std::string str ) {
    filename_ = str;
  }
  std::string const& filename() const {
    return filename_;
  }

  virtual void set_format_from_peak( CrossPeak const& ) = 0;
  virtual void write_header( std::ostream& ) = 0;
  virtual bool compatible_with_current_format( CrossPeak const& ) const = 0;
  virtual CrossPeakOP read_peak( std::istream&, std::string &next_line ) const = 0;
  virtual void output_diagnosis( std::ostream& ) const {};
  virtual void read_header( std::istream&, std::string& next_line ) = 0;

  bool ignore_assignments() const {
    return ignore_assignments_;
  }
  void set_ignore_assignments( bool setting = true ) {
    ignore_assignments_ = setting;
  }

  bool write_only_highest_VC() const {
    return write_only_highest_VC_;
  }

  void set_write_only_highest_VC( bool setting = true ) {
    write_only_highest_VC_ = setting;
  }

  core::Real min_VC_to_write() const {
    return min_VC_to_write_;
  }

  void set_min_VC_to_write( core::Real setting ) {
    min_VC_to_write_=setting;
  }

  bool ignore_negative_intensity() const {
    return ignore_negative_intensity_;
  }

  void set_ignore_negative_intensity( bool setting = true ) {
     ignore_negative_intensity_ = setting;
  }

  void set_minimum_peak_intensity( core::Real setting ) {
    minimum_peak_intensity_ = setting;
  }

  core::Real minimum_peak_intensity() const {
    return minimum_peak_intensity_;
  }

protected:
  CrossPeakInfoOP info1_;
  CrossPeakInfoOP info2_;
  std::string filename_;

private:
  bool ignore_assignments_;
  core::Real min_VC_to_write_;
  bool write_only_highest_VC_;
  bool ignore_negative_intensity_;
  core::Real minimum_peak_intensity_;
  //  virtual void write_header( std::ostream& );
};

class PeakFileFormat : public PeakFileFormat_Base {
public:
  PeakFileFormat();
  //  PeakFileFormat( ResonanceListOP const& );
  virtual ~PeakFileFormat();
  virtual void write_peak( std::ostream&, core::Size ct, CrossPeak const& ) const;
  virtual void write_resonances( std::ostream&, CrossPeak const& ) const;
  virtual void write_strength( std::ostream&, CrossPeak const& ) const;
  virtual void write_assignments( std::ostream&, CrossPeak const&, std::string const& first_line_end ) const;
  virtual void write_assignment( std::ostream&, PeakAssignment const& ) const;
  virtual void write_assignment_indent( std::ostream&, CrossPeak const& ) const;
  virtual void write_assignment_stats( std::ostream& os, PeakAssignment& pa ) const;
  virtual void write_nil_assignment( std::ostream&) const {};
  virtual void read_resonances( std::istream&, CrossPeak& ) const;
  virtual void read_assignments( std::istream& is, std::istream& rest_line, CrossPeak&, std::string& next_line ) const;
  virtual void read_strength( std::istream&, CrossPeak& ) const;

  virtual CrossPeakOP read_peak( std::istream&, std::string& next_line ) const;
  virtual void read_header( std::istream&, std::string& next_line );
  //  virtual void write_header( std::ostream& );
  virtual void output_diagnosis( std::ostream& ) const;

  virtual void set_format_from_peak( CrossPeak const& );
  virtual void write_header( std::ostream& );
  virtual bool compatible_with_current_format( CrossPeak const& ) const;

  static void register_options();

  void set_write_atom_names( bool setting = true ) {
    write_atom_names_ = setting;
  }
  bool write_atom_names() const {
    return write_atom_names_;
  }

  core::Size ncol() const { return col2proton_.size(); }
  //  ResonanceList const& resonances() const { return *resonances_; }
protected:
  utility::vector1< std::string > column_labels_; //eg. NhH
  utility::vector1< core::Size > col2proton_; //up to 4 columns may read 1 1 2 2 or 1 2 1 (NhH) or 1 2 2 for hHN
  utility::vector1< bool > col2islabel_;  // woulde be 0 0 1 for hHN and 1 0 0 for NhH

private:
  static bool options_registered_;
  //  ResonanceListOP resonances_;
  bool write_atom_names_;
};

class PeakFileFormat_xeasy : public PeakFileFormat {
public:
  PeakFileFormat_xeasy() : PeakFileFormat() {};
  //  PeakFileFormat_xeasy( ResonanceListOP const& rop ) : PeakFileFormat( rop ) {};
};

}
}

#endif
