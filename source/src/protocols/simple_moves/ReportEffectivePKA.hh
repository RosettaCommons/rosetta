// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
///@author Yifan Song

#ifndef INCLUDED_protocols_simple_moves_ReportEffectivePKA_hh
#define INCLUDED_protocols_simple_moves_ReportEffectivePKA_hh

#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/optimization/Multifunc.hh>
#include <core/optimization/Minimizer.hh>

namespace protocols {
namespace simple_moves {

class IonizableResidue {
public:
    IonizableResidue(
                     std::string resname,
                     core::Real ref_pKa,
                     core::Real acid_base
                     ) {
        resname_ = resname;
        ref_pKa_ = ref_pKa;
        acid_base_ = acid_base;
    }
    void add_neutral_restype(std::string restype) {
        neutral_restype_.push_back(restype);
    }
    void add_ionized_restype(std::string restype) {
        ionized_restype_.push_back(restype);
    }
    core::Real ref_pKa() {
        return ref_pKa_;
    }
    core::Real acid_base_coefficient() {
        return acid_base_;
    }
    utility::vector1 < std::string > & neutral_restypes() {
        return neutral_restype_;
    }
    utility::vector1 < std::string > & ionized_restypes() {
        return ionized_restype_;
    }
    std::string name3() {
        return resname_;
    }
    
private:
    std::string resname_;
    utility::vector1<std::string> neutral_restype_;
    utility::vector1<std::string> ionized_restype_;
    core::Real ref_pKa_;
    core::Real acid_base_;
};

/// scale density map intensities to match a pose's
class ReportEffectivePKA : public moves::Mover {
public:
	ReportEffectivePKA();
	virtual ~ReportEffectivePKA() {}

	void init();

	virtual void apply( core::pose::Pose & );

	std::string get_name() const { return "ReportEffectivePKA"; }

	moves::MoverOP clone() const { return new ReportEffectivePKA( *this ); }
	moves::MoverOP fresh_instance() const { return new ReportEffectivePKA; }

    void task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
    core::pack::task::TaskFactoryOP task_factory() const { return task_factory_; }

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & );

private:
    core::scoring::ScoreFunctionOP scorefxn_; //dflt NULL
	core::pack::task::TaskFactoryOP task_factory_;
	core::pack::task::PackerTaskOP task_;
    utility::vector1<IonizableResidue> ionizables_;
};

} // moves
} // protocols

#endif
