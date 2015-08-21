// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/WorkUnit.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_MoverList_hh
#define INCLUDED_protocols_wum_MoverList_hh


#include <core/types.hh>
#include <protocols/moves/Mover.hh>

#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace wum {

class MoverList;
//typedef  utility::pointer::access_ptr< MoverList const >  MoverListCAP;
typedef  const MoverList* MoverListCAP;

class MoverList{
public:
	MoverList(){}

	void register_mover( const std::string &name, moves::MoverCOP the_mover );

	moves::MoverCOP get_mover( const std::string &name ) const;
protected:
	std::map< std::string, moves::MoverCOP > mover_list_;
};


}
}

#endif

