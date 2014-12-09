// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragData.cc
/// @brief  a collection classes of the FragData and SingleResidueFragData class hirachy
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///

// Unit Headers
#include <core/fragment/SecstructSRFD.hh>

// Package Headers
// AUTO-REMOVED #include <core/fragment/FragData.hh>
// AUTO-REMOVED #include <core/fragment/Frame.hh>

// Project Headers
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace core {
namespace fragment {

static thread_local basic::Tracer tr( "core.fragment" );


bool SecstructSRFD::apply( pose::Pose& pose, Size seqpos ) const {
  //  Parent::apply( pose, seqpos );
  pose.set_secstruct( seqpos, secstruct() );
  return true; //can something go wrong ? check has_torsion() ?
}


/// @brief apply secondary structure fragment data to the pose, movemap has no effect
/// @remarks In this version of apply(), by convention MoveMap has no effect
///  because a setting for sec.struct currently does not exist within the map.
/// @return always true
bool SecstructSRFD::apply( kinematics::MoveMap const &, pose::Pose & pose, Size const seqpos ) const {
	// see is_applicable() for additional comments
	pose.set_secstruct( seqpos, secstruct() );
	return true;
}


bool SecstructSRFD::apply_ss( std::string& ss, Size seqpos ) const {
  // Parent::apply_ss( ss, seqpos );
 ss[ seqpos - 1 ] = secstruct();
 return true; //can something go wrong ? check has_torsion() ?
}


bool SecstructSRFD::steal( pose::Pose const& pose, Size seqpos ) {
  Parent::steal( pose, seqpos );
  secstruct_ = pose.secstruct( seqpos );
	tr.Trace << "steal secstructur " << secstruct_ << " at position " << seqpos  << std::endl;
  return true;
}

bool SecstructSRFD::is_compatible( SingleResidueFragData const& aSRFD) const {
  //SecstructSRFD const* ptr = dynamic_cast< SecstructSRFD const* > ( & aSRFD );
  return dynamic_cast< SecstructSRFD const* > ( & aSRFD ); //cast succesful same type
}

bool SecstructSRFD::is_applicable( kinematics::MoveMap const&, Size) const {
  //movemap always allows changes of the secstruct id ( it doesn't move anything but changes the energy -- strange dof )
  return true;
}

void SecstructSRFD::show( std::ostream &out ) const {
  using namespace ObjexxFCL::format;
	Parent::show( out );
  out << sequence() << ' ' << secstruct() << ' ';
}

void SecstructSRFD::read_data( std::istream &in ) {
	Parent::read_data( in );
	char c;
	in >> c >> secstruct_;
	set_sequence( c );
}

} //fragment
} //core
