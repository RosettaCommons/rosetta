// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragData.hh
/// @brief  A fragment as list of SingleResidue Data
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_SecstructSRFD_HH
#define INCLUDED_core_fragment_SecstructSRFD_HH

// Unit Headers
#include <core/fragment/SingleResidueFragData.hh>

// Package Headers
#include <core/fragment/Frame.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>

// AUTO-REMOVED #include <core/kinematics/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

//
class SecstructSRFD : public SingleResidueFragData {
  typedef SingleResidueFragData Parent;
public:
  SecstructSRFD( char secstruct = 'X', char sequence = 'X')
    : SingleResidueFragData( sequence ), secstruct_ ( secstruct )
  {};

  //
  SingleResidueFragDataOP clone() const {
    return SingleResidueFragDataOP( new SecstructSRFD( *this ) );
  };

	/// @brief create a new instance of this object
	virtual
	SingleResidueFragDataOP create() const {
		return SingleResidueFragDataOP( new SecstructSRFD() );
	}

  /// get secstruct for this position
  virtual char
  secstruct() const
  {
    return secstruct_;
  }

  void set_secstruct( char const ss ) {
    secstruct_ = ss;
  }

  virtual bool apply( pose::Pose&, Size seq_pos ) const;

  /// @brief apply secondary structure fragment data to the pose, movemap has no effect
  /// @remarks In this version of apply(), by convention MoveMap has no effect
  ///  because a setting for sec.struct currently does not exist within the map.
  /// @return always true
  virtual bool apply( kinematics::MoveMap const &, pose::Pose & pose, Size const seqpos ) const;

  virtual bool apply_ss( std::string&, Size seq_pos) const;
  virtual bool steal( pose::Pose const&, Size seq_pos );
  virtual bool is_compatible( SingleResidueFragData const& ) const;
  virtual bool is_applicable( kinematics::MoveMap const&, Size seq_pos ) const;

  virtual
  void show( std::ostream &out ) const;

	virtual
	void read_data( std::istream &in );

	virtual
	std::string type() const {
		return "Secstruct";
	}

	static std::string _static_type_name() {
		return "Secstruct";
	}

private:
  char secstruct_;
};

} //fragment
} //core

#endif
