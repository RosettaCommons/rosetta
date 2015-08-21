// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RestrictDesignToProteinDNAInterface.hh
/// @brief When operated on a PackerTask, pack/design is limited to the protein-DNA interface
/// @author ashworth

#ifndef INCLUDED_protocols_dna_RestrictDesignToProteinDNAInterface_hh
#define INCLUDED_protocols_dna_RestrictDesignToProteinDNAInterface_hh

#include <protocols/dna/RestrictDesignToProteinDNAInterface.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

#include <protocols/dna/DnaChains.fwd.hh>
#include <protocols/dna/DnaInterfaceFinder.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <protocols/dna/DnaDesignDef.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace dna {

class RestrictDesignToProteinDNAInterface : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseCOP PoseCOP;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef TaskOperation parent;

public:
	RestrictDesignToProteinDNAInterface();
	virtual ~RestrictDesignToProteinDNAInterface();

	virtual TaskOperationOP clone() const;

	virtual void apply( Pose const & pose, PackerTask & ptask ) const;

	void copy_dna_chains( DnaChainsCOP chains );
	DnaChainsCOP dna_chains() const;

	void copy_targeted_dna( DnaDesignDefOPs const & targeted_dna );
	DnaDesignDefOPs const & targeted_dna() const;

	void copy_interface( DnaInterfaceFinderCOP interface );
	DnaInterfaceFinderCOP interface() const;

	void set_reference_pose( PoseCOP pose );
	PoseCOP reference_pose() const;

	void set_base_only( bool value ) { base_only_ = value; }
	bool base_only() const { return base_only_; }

	void set_forget_chains_and_interface( bool value ) { forget_chains_and_interface_ = value; }
	bool forget_chains_and_interface() const { return forget_chains_and_interface_; }

	virtual void parse_tag( TagCOP, DataMap & );

private:
	// reference pose: sometimes the input pose is not actually 'native,' and it is useful to make the Task aware of an alternative 'native' pose (such as for reverting mutations that did not occur in a local context)
	mutable DnaChainsOP dna_chains_; // mutable to allow on-the-fly generation of default, and nulling to prevent accumulation of state
	DnaDesignDefOPs targeted_dna_;
	mutable DnaInterfaceFinderOP interface_; // mutable to allow on-the-fly generation of default, and nulling to prevent accumulation of state
	PoseCOP reference_pose_;
	bool base_only_;
	bool forget_chains_and_interface_; // WARNING: default true, false at own risk
	mutable core::Real z_cutoff_; // mutable to allow lazy default initialization
	core::Real close_threshold_, contact_threshold_;
};

} // namespace dna
} // namespace protocols

#endif
