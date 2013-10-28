// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RemodelRotamerLinks.hh
/// @brief When operated on a PackerTask, adds RotamerLinks to it that will
//trigger the use of a FixbbLinkingRotamerAnnealer by pack_rotamers, which will
//multiple sidechain substitutions in the sim annealing process

#ifndef INCLUDED_protocols_forge_remodel_RemodelRotamerLinks_hh
#define INCLUDED_protocols_forge_remodel_RemodelRotamerLinks_hh

#include <protocols/forge/remodel/RemodelRotamerLinks.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

namespace protocols {
namespace forge {
namespace remodel {

class RemodelRotamerLinks : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pose::Pose Pose;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pack::task::PackerTask PackerTask;

public:
	RemodelRotamerLinks();
	virtual ~RemodelRotamerLinks();
	virtual TaskOperationOP clone() const;
	virtual	void apply( Pose const & pose, PackerTask & ptask ) const;
	virtual void parse_tag( TagCOP, DataMap & );

};

} // namespace remodel
} // namespace forge
} // namespace protocols

#endif
