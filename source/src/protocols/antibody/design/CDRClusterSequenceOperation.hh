// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file /CDRClusterSequenceOperation.hh
/// @brief
/// @author jadolfbr ()


#ifndef INCLUDED_protocols_antibody_design_CDRClusterSequenceOperation_hh
#define INCLUDED_protocols_antibody_design_CDRClusterSequenceOperation_hh

#include <protocols/antibody/design/CDRClusterSequenceOperation.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
namespace design {
			

	
	
class CDRClusterSequenceOperation : public utility::pointer::ReferenceCount {
public:

	CDRClusterSequenceOperation();
	CDRClusterSequenceOperation(CDRClusterSequenceOperation const & src);

	virtual ~CDRClusterSequenceOperation();

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;
	

};

} //design
} //antibody
} //protocols

#endif	//INCLUDED_protocols_antibody_design_CDRClusterSequenceOperation_hh



