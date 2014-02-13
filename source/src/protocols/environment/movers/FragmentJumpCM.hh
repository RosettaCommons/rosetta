// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/FragmentJumpCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_FragmentJumpCM_hh
#define INCLUDED_protocols_environment_FragmentJumpCM_hh

// Unit Headers
#include <protocols/environment/movers/FragmentJumpCM.fwd.hh>

// Package headers
#include <protocols/environment/movers/FragmentCM.hh>

// Project headers
#include <protocols/jumping/JumpSample.hh>
#include <protocols/abinitio/TemplateJumpSetup.hh>


// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class FragmentJumpCM : public FragmentCM {
  typedef FragmentCM Parent;
public:
  FragmentJumpCM();

  FragmentJumpCM( std::string const& topol_filename,
                  std::string const& label,
                  std::string const& moverkey );

  virtual
  ~FragmentJumpCM() {};

  virtual
  claims::EnvClaims yield_claims( core::pose::Pose& );

  virtual std::string get_name() const;

  /// @brief Set the topology using a topology file.
  void set_topology( std::string const& topology_file );

  /// @brief Set the topology using a psipred_ss2 file, a pairing
  ///        list file, and a number of random sheets to pick.
  void set_topology( std::string const& ss_info_file,
                     std::string const& pairing_file,
                     core::Size const& n_sheets,
                     bool bRandomSheets );

  void set_moverkey( std::string const& key ) { moverkey_ = key; }

  std::string const& moverkey() { return moverkey_; }

  virtual void
  parse_my_tag( utility::tag::TagCOP const tag,
                basic::datacache::DataMap & data,
                protocols::filters::Filters_map const & filters,
                protocols::moves::Movers_map const & movers,
                core::pose::Pose const & pose );

  virtual
  moves::MoverOP fresh_instance() const;

  virtual
  moves::MoverOP clone() const;

protected:
  jumping::JumpSample calculate_jump_sample() const;

  jumping::JumpSample setup_fragments();

  claims::EnvClaims build_claims( jumping::JumpSample const& );

private:
  jumping::BaseJumpSetupOP jump_def_;

  // Hash is used as a quick and dirty check in the JumpSampleData
  // to make sure that we're not loading in total garbage
  std::string moverkey_;

}; // end FragmentJumpCM base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_FragmentJumpCM_hh
