// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///


#ifndef INCLUDED_protocols_jd2_MpiFileBuffer_fwd_hh
#define INCLUDED_protocols_jd2_MpiFileBuffer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class MpiFileBuffer;
typedef utility::pointer::shared_ptr< MpiFileBuffer > MpiFileBufferOP;
typedef utility::pointer::shared_ptr< MpiFileBuffer const > MpiFileBufferCOP;


}
}

#endif
