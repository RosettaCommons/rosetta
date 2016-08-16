// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/PyEnergy.hh
/// @brief  Various Energy classes for subclassing in PyRosetta.
/// @author Sergey Lyskov


#ifndef INCLUDED_core_scoring_methods_PyEnergy_hh
#define INCLUDED_core_scoring_methods_PyEnergy_hh

#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

/*
template< class T > utility::pointer::access_ptr<T>         _AP( T & o)        { return utility::pointer::access_ptr<T>        ( & o ); }
template< class T > utility::pointer::access_ptr<T const > _CAP( T const & o ) { return utility::pointer::access_ptr<T const > ( & o ); }


class PyContextIndependentOneBodyEnergy : public ContextIndependentOneBodyEnergy
{
public:
    PyContextIndependentOneBodyEnergy( EnergyMethodCreatorOP op) : ContextIndependentOneBodyEnergy(op) {}


    virtual void residue_energy(conformation::Residue const & rsd, pose::Pose const & pose, EnergyMap & emap) const {
        Py_residue_energy( _AP(rsd), _AP(pose), _AP(emap) ) ;
    }

    virtual void Py_residue_energy( utility::pointer::access_ptr< conformation::Residue const > const & ,
        pose::PoseCAP const &,
        utility::pointer::access_ptr< EnergyMap > const &
    ) const = 0;
};


class PyContextIndependentTwoBodyEnergy  : public ContextIndependentTwoBodyEnergy
{
public:
    PyContextIndependentTwoBodyEnergy( EnergyMethodCreatorOP op) : ContextIndependentTwoBodyEnergy(op) {}


    virtual void residue_pair_energy(conformation::Residue const & rsd1, conformation::Residue const & rsd2, pose::Pose const & pose, ScoreFunction const &sfn, EnergyMap & emap) const
    { Py_residue_pair_energy(_CAP(rsd1), _CAP(rsd2), _CAP(pose), _CAP(sfn), _AP(emap) ); }

    virtual void Py_residue_pair_energy(utility::pointer::access_ptr< conformation::Residue const > const & rsd1,
                                       utility::pointer::access_ptr< conformation::Residue const > const & rsd2,
                                       utility::pointer::access_ptr< pose::Pose const> const & pose,
                                       utility::pointer::access_ptr< ScoreFunction const > const & sfn,
                                       utility::pointer::access_ptr< EnergyMap > const & emap) const = 0;


    virtual bool defines_intrares_energy( EnergyMap const & emap ) const
    { return Py_defines_intrares_energy( _CAP(emap) ); }

    virtual bool Py_defines_intrares_energy( utility::pointer::access_ptr< EnergyMap const > const & emap ) const = 0;


    virtual void eval_intrares_energy(conformation::Residue const & rsd, pose::Pose const & pose, ScoreFunction const &sfn, EnergyMap & emap) const
    { Py_eval_intrares_energy( _CAP(rsd), _CAP(pose), _CAP(sfn), _AP(emap)); }

    virtual void Py_eval_intrares_energy(utility::pointer::access_ptr< conformation::Residue const> const & rsd,
                                         utility::pointer::access_ptr< pose::Pose const > const & pose,
                                         utility::pointer::access_ptr< ScoreFunction const > const & sfn,
                                         utility::pointer::access_ptr< EnergyMap > const & emap) const = 0;


    virtual Real eval_dof_derivative(id::DOF_ID const & dof_id, id::TorsionID const & tor_id, pose::Pose const & pose, ScoreFunction const & sfxn, EnergyMap const & weights) const
    { return Py_eval_dof_derivative(_CAP(dof_id), _CAP(tor_id), _CAP(pose), _CAP(sfxn), _CAP(weights)); }

    virtual Real Py_eval_dof_derivative(utility::pointer::access_ptr< id::DOF_ID const > const & dof_id,
                                        utility::pointer::access_ptr< id::TorsionID const > const & tor_id,
                                        utility::pointer::access_ptr< pose::Pose const > const & pose,
                                        utility::pointer::access_ptr< ScoreFunction const > const & sfxn,
                                        utility::pointer::access_ptr< EnergyMap const > const & weights) const = 0;


    virtual void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
    { Py_indicate_required_context_graphs( _AP(context_graphs_required) ); }

    virtual void Py_indicate_required_context_graphs( utility::pointer::access_ptr<  utility::vector1< bool > > const & context_graphs_required ) const = 0;

};
*/

}
}
}

#endif // INCLUDED_core_scoring_methods_PyEnergy_hh
