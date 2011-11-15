// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_CrossPeakInfo_hh
#define INCLUDED_protocols_noesy_assign_CrossPeakInfo_hh


// Unit Headers
#include <protocols/noesy_assign/CrossPeakInfo.fwd.hh>

// Package Headers
// AUTO-REMOVED #include <protocols/noesy_assign/ResonanceList.fwd.hh>

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.fwd.hh>
#include <core/chemical/AA.hh>

// Utility headers
// AUTO-REMOVED #include <utility/exit.hh>
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
// AUTO-REMOVED #include <cstdlib>
#include <string>
// AUTO-REMOVED #include <list>
// AUTO-REMOVED #include <map>

namespace protocols {
namespace noesy_assign {


///@brief shared information about CrossPeak (one for each proton dimension)
class CrossPeakInfo : public utility::pointer::ReferenceCount {
public:
  CrossPeakInfo( std::string const& proton, std::string const& label, core::Real proton_tolerance, core::Real label_tolerance = 0.0 ) :
    proton_atom_name_( proton ),
    label_atom_type_( label ),
    proton_tolerance_( proton_tolerance ),
    label_tolerance_( label_tolerance )
  {}

  std::string const& main_atom() const { return proton_atom_name_; }
  std::string const& label_atom_type() const { return label_atom_type_; }

  core::Real const& proton_tolerance() const { return proton_tolerance_; }
  core::Real const& label_tolerance() const { return label_tolerance_; }
  bool has_label() const { return label_atom_type_ != "" && label_tolerance_ < 200; } //tolerance is sometimes set to 999

  //returns the corresponding label-atom depending on type: (HN) H --> N, (CH) QD1-->CD1 etc.
  std::string label_atom_name( std::string const& proton_name, core::chemical::AA aa ) const;

  void set_label( std::string label, core::Real tolerance ) {
    label_atom_type_ = label;
    label_tolerance_ = tolerance;
  }
  void set_proton( std::string name, core::Real tolerance ) {
    proton_atom_name_ = name;
    proton_tolerance_ = tolerance;
  }

  bool operator ==( CrossPeakInfo const& cpi ) const {
    return cpi.proton_atom_name_ == proton_atom_name_
      && cpi.label_atom_type_ == label_atom_type_
      && cpi.proton_tolerance_ == proton_tolerance_
      && cpi.label_tolerance_ == label_tolerance_;
  }

  void set_filename( std::string filename ) {
    filename_ = filename;
  }

  std::string const& filename() const {
    return filename_;
  }


private:
  std::string proton_atom_name_;
  std::string label_atom_type_;
  core::Real proton_tolerance_;
  core::Real label_tolerance_;
  std::string filename_;
};

}
}

#endif
