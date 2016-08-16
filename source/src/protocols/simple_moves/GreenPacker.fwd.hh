// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/GreenPacker.fwd.hh
/// @brief  packing mover that makes extensive reuse of rotamer pair energies class forward declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_GreenPacker_fwd_hh
#define INCLUDED_protocols_simple_moves_GreenPacker_fwd_hh


/// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

class MinimalRotamer;

class GroupDiscriminator;
class ChainGroupDiscriminator;
class UserDefinedGroupDiscriminator;

class GreenPacker;

typedef utility::pointer::shared_ptr< MinimalRotamer > MinimalRotamerOP;
typedef utility::pointer::shared_ptr< MinimalRotamer const > MinimalRotamerCOP;
typedef utility::pointer::shared_ptr< GroupDiscriminator > GroupDiscriminatorOP;
typedef utility::pointer::shared_ptr< GroupDiscriminator const > GroupDiscriminatorCOP;
typedef utility::pointer::shared_ptr< UserDefinedGroupDiscriminator > UserDefinedGroupDiscriminatorOP;
typedef utility::pointer::shared_ptr< GreenPacker > GreenPackerOP;
typedef utility::pointer::shared_ptr< GreenPacker const > GreenPackerCOP;

} // namespace moves
} // namespace protocols

#endif // INCLUDED_protocols_simple_moves_GreenPacker_FWD_HH
