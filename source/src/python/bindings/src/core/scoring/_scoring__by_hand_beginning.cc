// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#include "boost/python.hpp"

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/SecondaryStructureWeights.hh>
#include <core/scoring/Energies.hh>


#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunctionInfo.hh>


namespace bp = boost::python;

// We no exactly need this function anymore - we can just call ScoringFunctionFactory from Python, so this is only for backward compatibility...
::core::scoring::ScoreFunctionOP create_scoring_function(::std::string const & s)
{
	return ::core::scoring::ScoreFunctionFactory::create_score_function(s);
}

::core::scoring::ScoreFunctionOP create_scoring_function_ws_patch(::std::string const & s, ::std::string const & p)
{
	return ::core::scoring::ScoreFunctionFactory::create_score_function(s, p);
}


//BOOST_PYTHON_FUNCTION_OVERLOADS(create_scoring_function_overloads, create_scoring_function, 1, 2)


void __scoring_by_hand_beginning__()
{
    // bp::def("create_score_function", create_scoring_function);
    // bp::def("create_score_function_ws_patch", create_scoring_function_ws_patch);


    // bp::implicitly_convertible< core::scoring::ScoreFunction const &, core::scoring::ScoreFunctionInfo >();

    // bp::implicitly_convertible< utility::pointer::owning_ptr< ::core::scoring::ScoreFunction >
    //                           , utility::pointer::owning_ptr< ::core::scoring::ScoreFunction const > >();
    // bp::implicitly_convertible< utility::pointer::owning_ptr< ::core::scoring::NeighborList >
    //                           , utility::pointer::owning_ptr< ::core::scoring::NeighborList const > >();
    // bp::implicitly_convertible< utility::pointer::owning_ptr< ::core::scoring::ScoreFunctionInfo >
    //                           , utility::pointer::owning_ptr< ::core::scoring::ScoreFunctionInfo const > >();
    // bp::implicitly_convertible< utility::pointer::owning_ptr< ::core::scoring::SecondaryStructureWeights >
    //                           , utility::pointer::owning_ptr< ::core::scoring::SecondaryStructureWeights const > >();

}
