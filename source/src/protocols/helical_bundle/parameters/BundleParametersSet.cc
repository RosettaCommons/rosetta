// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/helical_bundle/parameters/BundleParametersSet.cc
/// @brief  A class for holding sets of parameters for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>

// Package headers
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>

// Project headers

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/assert.hh>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>


namespace protocols {
	namespace helical_bundle {
		namespace parameters {

			static thread_local basic::Tracer TR( "protocols.helical_bundle.parameters.BundleParametersSet" );

			/// @brief Constructor.
			///
			BundleParametersSet::BundleParametersSet() :
				bundle_symmetry_(0),
				bundle_symmetry_copies_(0),
				n_helices_(0)
			{
			}

			BundleParametersSet::BundleParametersSet( BundleParametersSet const & src ) :
				core::conformation::parametric::ParametersSet(src),
				bundle_symmetry_( src.bundle_symmetry_ ),
				bundle_symmetry_copies_( src.bundle_symmetry_copies_ ),
				n_helices_(src.n_helices_)
			{
			}

			BundleParametersSet::~BundleParametersSet() {}


			/// @brief make a copy of this residue( allocate actual memory for it )
			///
			core::conformation::parametric::ParametersSetOP
			BundleParametersSet::clone() const
			{
				return ParametersSetOP( new BundleParametersSet( *this ) );
			}

			/// @brief Get a summary of this ParametersSet object, for output to remark lines of a PDB file.
			/// @details Default function can be overridden by derived classes.  This version actually outputs
			/// Crick parameter information.
			void BundleParametersSet::get_pdb_remark(std::stringstream &remark) const {
				//                                                                    X Anything from the X onward gets truncated
				remark << "----------------------------------------"                  << std::endl;
				remark << "CRICK PARAMETER INFORMATION"                               << std::endl;
				remark << "Bundle symmetry: " << bundle_symmetry() << (bundle_symmetry()==0 ? " (not symmetric)" : "") << std::endl;
				remark << "Bundle symmetry copies: " << bundle_symmetry_copies()      << std::endl;
				remark << "Number of helices per copy: " << n_helices()               << std::endl;
				if(n_helices() > 0) {
					for(core::Size i=1, imax=n_helices(); i<=imax; ++i) {
						BundleParametersCOP cur_helix( utility::pointer::dynamic_pointer_cast<BundleParameters const>( parameters(i) ) );
						remark << "---HELIX " << i << " PARAMETERS:---" << std::endl;
						cur_helix->get_pdb_remark(remark);
					}
				}
				remark << "----------------------------------------"                  << std::endl;
				return;
			}

		} // namespace parameters
	} // namespace helical_bundle
} // namespace protocols

