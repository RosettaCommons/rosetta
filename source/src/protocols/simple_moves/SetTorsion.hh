// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#ifndef INCLUDED_protocols_simple_moves_SetTorsion_hh
#define INCLUDED_protocols_simple_moves_SetTorsion_hh

#include <protocols/simple_moves/SetTorsion.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

#include <utility/vector1.hh>


// Utility headers

// C++ headers

// Unit headers

namespace protocols {
namespace simple_moves {

/// @brief A mover to change one torsion angle
class SetTorsion : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	///@brief default ctor
	SetTorsion();
	virtual ~SetTorsion();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return (protocols::moves::MoverOP( new protocols::simple_moves::SetTorsion( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new SetTorsion );
	}
    
    core::Size n_torsion_sets() const {
        return residues_.size();
    }
    
    utility::vector1<core::Size> residue_list(core::Size iset, core::pose::Pose const & pose);
    core::Real angle(core::Size iset) const;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
    
	//core::Real angle(core::Size const iset) const;
    //std::list <core::Size> residue(core::Size const iset) const;
	std::string torsion_name(core::Size const iset) {
        return torsion_name_[iset];
    }

	/// @brief Sets a residue index that will serve as the root of the FoldTree for the SetTorsion operation.
	/// @details The FoldTree is reset afterwards (i.e. the mover does not permanently change the FoldTree).
	void set_fold_tree_root(core::Size const root) { fold_tree_root_ = root; return; }

	///
	/// @brief Returns the residue index that will serve as the root of the FoldTree for the SetTorsion operation.
	core::Size get_fold_tree_root() const { return fold_tree_root_; }

private:
    bool random_set_;
	utility::vector1<std::string> angle_;
	utility::vector1<std::string> residues_;
	utility::vector1<std::string> torsion_name_; // phi/psi etc.
    utility::vector1<Size> extending_;
    utility::vector1< utility::vector1< core::id::NamedAtomID > > torsion_atoms_;

		///
		/// @brief The root for the FoldTree during the operation.  If 0, the default FoldTree is used.  The FoldTree is reset to the input FoldTree after the operation.
		core::Size fold_tree_root_;

};

} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_SetTorsion_HH_
