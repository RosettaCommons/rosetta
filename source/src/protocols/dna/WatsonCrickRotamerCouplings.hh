// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file WatsonCrickRotamerCouplings.hh
/// @brief When operated on a PackerTask, adds RotamerCouplings to it that will trigger the use of a FixbbCoupledRotamerAnnealer by pack_rotamers, which will make base-paired DNA rotamer substitutions

#ifndef INCLUDED_protocols_dna_WatsonCrickRotamerCouplings_hh
#define INCLUDED_protocols_dna_WatsonCrickRotamerCouplings_hh

#include <protocols/dna/WatsonCrickRotamerCouplings.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

#include <utility/vector1.hh>

#include <utility/tag/Tag.fwd.hh>


namespace protocols {
namespace dna {

class WatsonCrickRotamerCouplings : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pose::Pose Pose;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pack::task::PackerTask PackerTask;

public:
	virtual ~WatsonCrickRotamerCouplings();
	virtual TaskOperationOP clone() const;
	virtual	void apply( Pose const & pose, PackerTask & ptask ) const;
	virtual void parse_tag( TagCOP, DataMap & );
};

} // namespace dna
} // namespace protocols

#endif
