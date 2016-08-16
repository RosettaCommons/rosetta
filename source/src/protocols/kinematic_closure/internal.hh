// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @internal
/// @brief  Header file for internal helpers.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_kinematic_closure_internal_HH
#define INCLUDED_protocols_kinematic_closure_internal_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>

// Utility headers
#include <boost/iterator/iterator_facade.hpp>

namespace protocols {
namespace kinematic_closure {

/// @cond DOCUMENT_KIC_INTERNALS

/// @brief Catalog ideal geometries for a number of different lengths, angles,
/// and torsions.
/// @details This class should definitely not exist.  It is basically just a
/// bunch of magic numbers that must exist somewhere in the database already.
/// Unfortunately, I can't find them and so we're stuck with this nonsense.

struct IdealParameters {

	// Rosetta-derived parameters.
	static const Real c_n_ca_angle;
	static const Real n_ca_c_angle;
	static const Real ca_c_n_angle;
	static const Real c_n_length;
	static const Real n_ca_length;
	static const Real ca_c_length;
	static const Real omega_dihedral;

	// PDB-derived parameters.
	static const Real mean_n_ca_c_angle;
	static const Real std_n_ca_c_angle;
	static const Real min_n_ca_c_angle;
	static const Real max_n_ca_c_angle;

};

/// @brief Two global parameters that are sometimes useful for debugging.  They
/// may be removed at any time.

extern Size num_rama_filter_fails;
extern Size num_bump_filter_fails;

/// @brief Iterate sequentially through two solution lists.
/// @details This functionality is very important to the balanced solution
/// picker.  While it could be done without this class (basically by doubling
/// the number of for-loops) I think that this class meaningfully improves the
/// readability of the code.

class ChainedSolutionList {

public:

	ChainedSolutionList(SolutionList const &a, SolutionList const &b):
		first_(a), second_(b) {}

	class Iterator : public boost::iterator_facade // {{{1
		<Iterator, ClosureSolutionCOP,
		boost::forward_traversal_tag, ClosureSolutionCOP> {

	public:
		//Iterator();

		Iterator(ChainedSolutionList const *parent, SolutionList::const_iterator bookmark);

	private:
		friend class boost::iterator_core_access;

		void increment() {

			++bookmark_;

			if ( state_ == FIRST && bookmark_ == parent_->first_.end() ) {
				bookmark_ = parent_->second_.begin();
				state_ = SECOND;
			}
		}

		bool equal (Iterator const &other) const {
			return bookmark_ == other.bookmark_;
		}

		ClosureSolutionCOP const dereference() const {
			return *bookmark_;
		}

		ChainedSolutionList const *parent_;
		SolutionList::const_iterator bookmark_;
		enum {FIRST, SECOND} state_;

	}; // }}}1

	// These typedefs allow chained solution lists within boost foreach-loops.
	typedef Iterator iterator;
	typedef Iterator const_iterator;

	Iterator begin() { return Iterator(this, first_.begin()); }
	Iterator end() { return Iterator(this, second_.end()); }

	ClosureSolutionCOP const front() { return first_.front(); }
	ClosureSolutionCOP const back() { return second_.back(); }

	bool size() { return first_.size() + second_.size(); }
	bool empty() { return first_.empty() && second_.empty(); }
	bool not_empty() { return !empty(); }

private:
	SolutionList const &first_;
	SolutionList const &second_;

};

} // end namespace kinematic_closure
} // end namespace protocols

/// @endcond DOCUMENT_KIC_INTERNALS

#endif
