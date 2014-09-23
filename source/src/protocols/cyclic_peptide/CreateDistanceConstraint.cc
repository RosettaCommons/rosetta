// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
///
/// @file   protocols/cyclic_peptide/CreateDistanceConstraint.cc
/// @brief  Add distance constraints to the current pose conformation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author Yifan Song

#include <protocols/cyclic_peptide/CreateDistanceConstraint.hh>
#include <protocols/cyclic_peptide/CreateDistanceConstraintCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/id/AtomID.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <istream>

static thread_local basic::Tracer TR( "protocols.cyclic_peptide.CreateDistanceConstraint" );

namespace protocols {
namespace cyclic_peptide {

CreateDistanceConstraint::CreateDistanceConstraint() //:
{}
CreateDistanceConstraint::~CreateDistanceConstraint(){}

void CreateDistanceConstraint::apply( core::pose::Pose & pose )
{
    for (Size i_cst=1; i_cst<=cst_func_.size(); ++i_cst) {
        if (cst_func_[i_cst] == "") {
            ///TODO: use ideal bond length as default parameter
            
            //Real length, deviation;
            //new core::scoring::func::HarmonicFunc(length, deviation);
            /*
        pose.add_constraint(
        new core::scoring::constraints::AtomPairConstraint( core::id::AtomID(iatom,ires), core::id::AtomID(jatom,jres),
        new core::scoring::func::ScalarWeightedFunc( cst_weight_, new core::scoring::func::Ha( dist, coord_dev_ ) ) ) );
             */
            
        }
        else {
            std::istringstream data(cst_func_[i_cst]);
            std::string func_type;
            data >> func_type;
            core::scoring::func::FuncFactory func_factory;
            core::scoring::func::FuncOP func = func_factory.new_func( func_type );
            func->read_data( data );
            Size atomno1 = pose.residue_type(res1_[i_cst]).atom_index(atom1_[i_cst]);
            Size atomno2 = pose.residue_type(res2_[i_cst]).atom_index(atom2_[i_cst]);
            pose.add_constraint(
                core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AtomPairConstraint( core::id::AtomID(atomno1,res1_[i_cst]), core::id::AtomID(atomno2,res2_[i_cst]), func ) )
                );
        }
    }
	return;
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
CreateDistanceConstraint::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
    utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for (tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it) {
        if ( (*tag_it)->getName() == "Add" ) {
            res1_.push_back( (*tag_it)->getOption< Size >( "res1" ) );
            atom1_.push_back( (*tag_it)->getOption< std::string >( "atom1" ) );
            res2_.push_back( (*tag_it)->getOption< Size >( "res2" ) );
            atom2_.push_back( (*tag_it)->getOption< std::string >( "atom2" ) );
            cst_func_.push_back( (*tag_it)->getOption< std::string >( "cst_func", "" ) );
        }
    }
	return;
}
	
moves::MoverOP CreateDistanceConstraint::clone() const { return moves::MoverOP( new CreateDistanceConstraint( *this ) ); }
moves::MoverOP CreateDistanceConstraint::fresh_instance() const { return moves::MoverOP( new CreateDistanceConstraint ); }

protocols::moves::MoverOP
CreateDistanceConstraintCreator::create_mover() const {
	return protocols::moves::MoverOP( new CreateDistanceConstraint );
}

std::string
CreateDistanceConstraintCreator::keyname() const
{
	return CreateDistanceConstraintCreator::mover_name();
}

std::string
CreateDistanceConstraintCreator::mover_name()
{
	return "CreateDistanceConstraint";
}

std::string
CreateDistanceConstraint::get_name() const {
	return "CreateDistanceConstraint";
}
	
} // moves
} // protocols
