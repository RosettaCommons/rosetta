// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_SecStructInfo.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_protocols_rna_RNA_SecStructInfo_hh
#define INCLUDED_protocols_rna_RNA_SecStructInfo_hh


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// Utility headers

// Numceric Headers

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace farna {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
//// Rhiju move this to its own namespace!
	class RNA_SecStructInfo: public basic::datacache::CacheableData  {

public:

RNA_SecStructInfo(): initialized_(false) {};

  RNA_SecStructInfo( std::string const rna_secstruct_string ) { rna_secstruct_  = rna_secstruct_string; }

  RNA_SecStructInfo( RNA_SecStructInfo const & src );

 basic::datacache::CacheableDataOP
  clone() const
  {
    return basic::datacache::CacheableDataOP( new RNA_SecStructInfo( *this ) );
  }

 	// Undefinded, comented out to make python bindings complile
	//void
	//update( core::pose::Pose const & pose );

  Size
  size() const {
    return rna_secstruct_.size();
  }

  bool
  initialized() const
  {
    return initialized_;
  }

  bool &
  initialized()
  {
    return initialized_;
  }

  void
  set_initialized( bool const & setting)
  {
    initialized_ = setting;
  }

	void
	set_secstruct( std::string const secstruct ){ rna_secstruct_ = secstruct; initialized_ = true; }

	std::string const &
	get_secstruct() const { return rna_secstruct_; }

private:

	std::string rna_secstruct_;

  bool initialized_;

};

std::string const &
get_rna_secstruct( core::pose::Pose & pose );

void
set_rna_secstruct( core::pose::Pose & pose, std::string const & rna_secstruct_string ); //By default, unknown, actually.


} //farna
} //protocols

#endif
