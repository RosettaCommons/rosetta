// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/FabConstraint.hh
/// @brief This class is specific to antibodies and penalizes presence of residues flanking
/// @brief antibody cdr residues at Antigen-Antibody interfaces (ported from Fab constraint
/// @brief in rosetta++ which uses a constant constraint score of 0.5/flanking residue)
/// @author Krishna Kilambi (kkpraneeth@jhu.edu, April 2012)

#ifndef INCLUDED_core_scoring_constraints_FabConstraint_hh
#define INCLUDED_core_scoring_constraints_FabConstraint_hh

// Unit header

#include <core/scoring/constraints/FabConstraint.fwd.hh>
#include <core/scoring/constraints/MultiConstraint.hh>


#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/SequenceMapping.fwd.hh>


#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>


//Utility Headers
#include <numeric/xyzVector.fwd.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace constraints {


class FabConstraint : public MultiConstraint {
public:

    /// @brief Constructor
    FabConstraint();

    /// @brief Constructor
    FabConstraint(ConstraintCOPs & cst_in) ;

    ///
    virtual
    ConstraintOP clone() const {
        return ConstraintOP( new FabConstraint(*this) );
    }

    std::string type() const {
        return "FabConstraint";
    }

    void
    show(std::ostream& out) const;

    /// @brief read in constraint definition
    void
    read_def(std::istream& data, pose::Pose const& pose,func::FuncFactory const& func_factory);

    Size
    pose_res_no(core::pose::Pose const & pose, std::string tempres);

    utility::vector1<Real>
    calc_penalty_vec(Size start_res, Size stop_res, utility::vector1<Size> res1, utility::vector1<Size> res2);

    void
    setup_csts(core::pose::Pose const & pose, utility::vector1<Size> res1, utility::vector1<Size> res2, std::string antchains);

private:



}; //FabConstraint

} //constraints
} //scoring
} //core

#endif
