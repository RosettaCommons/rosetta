// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_assembly/RedesignDomainAssemblyMover.hh
/// @brief Declaration of a class derived from devel::domain_assembly::DomainAssemblyMover which redesigns the domain interface
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

#ifndef INCLUDED_devel_domain_assembly_RedesignDomainAssemblyMover_hh
#define INCLUDED_devel_domain_assembly_RedesignDomainAssemblyMover_hh

// Unit Headers
#include <devel/domain_assembly/RedesignDomainAssemblyMover.fwd.hh>
#include <devel/domain_assembly/DomainAssemblyMover.hh>

namespace devel {
namespace domain_assembly {

class RedesignDomainAssemblyMover : public DomainAssemblyMover
{
public:
	RedesignDomainAssemblyMover();
	~RedesignDomainAssemblyMover();
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
private: // protocol stages
	virtual void run_fullatom_stage( core::pose::Pose & pose );
	virtual void run_fullatom_relax( core::pose::Pose & pose );
private: // initializers
	void initialize();
	void initialize_repack_only();
private: // data
	std::string residues_to_repack_only_;
};

}
}


#endif
