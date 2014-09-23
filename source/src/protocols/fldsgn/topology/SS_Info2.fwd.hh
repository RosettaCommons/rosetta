// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/SS_Info2.fwd.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )


#ifndef INCLUDED_protocols_fldsgn_topology_SS_Info2_fwd_hh
#define INCLUDED_protocols_fldsgn_topology_SS_Info2_fwd_hh

/// Utility headers
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace fldsgn {
namespace topology {

	class SS_Base;
	class Strand;
	class Helix;
	class Loop;
	class SS_Info2;

	typedef utility::pointer::shared_ptr< SS_Base > SS_BaseOP;
	typedef utility::pointer::shared_ptr< Helix > HelixOP;
	typedef utility::pointer::shared_ptr< Strand > StrandOP;
	typedef utility::pointer::shared_ptr< Loop > LoopOP;

	typedef utility::pointer::shared_ptr< SS_Base const > SS_BaseCOP;
	typedef utility::pointer::shared_ptr< Helix const > HelixCOP;
	typedef utility::pointer::shared_ptr< Strand const > StrandCOP;
	typedef utility::pointer::shared_ptr< Loop const > LoopCOP;


	typedef utility::vector1< HelixOP > Helices;
	typedef utility::vector1< StrandOP > Strands;
	typedef utility::vector1< LoopOP > Loops;

	typedef utility::pointer::shared_ptr< SS_Info2 > SS_Info2_OP;
	typedef utility::pointer::shared_ptr< SS_Info2 const > SS_Info2_COP;

} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
