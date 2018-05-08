// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/data_storage/HashedSmartAssembly.hh
/// @brief an Assembly that makes use of the Hasher
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_data_storage_HashedSmartAssembly_hh
#define INCLUDED_protocols_sewing_data_storage_HashedSmartAssembly_hh

#include <protocols/sewing/data_storage/HashedSmartAssembly.fwd.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sewing {
namespace data_storage {

///@brief an Assembly that makes use of the Hasher
class HashedSmartAssembly : public SmartAssembly {

public:

	HashedSmartAssembly();
	HashedSmartAssembly(hashing::SegmentVectorCOP);
	HashedSmartAssembly(HashedSmartAssembly const & src);

	virtual ~HashedSmartAssembly();

	HashedSmartAssemblyOP
	clone() const;

	void
	set_basis_map_generator(hashing::BasisMapGeneratorOP new_bmg);

	bool
	add_segment(bool n_terminus, core::Size seg_ID_to_add, core::Size res_ID_1, core::Size res_ID_2 ) override;

	core::Size
	get_window_width() const;

private:
	hashing::BasisMapGeneratorOP existing_alignments_;
	//core::Size segID_1_;
	//core::Size segID_2_;
	SmartSegmentOP first_segment_;
};


} //protocols
} //sewing
} //data_storage



#endif //INCLUDED_protocols_sewing_data_storage_HashedSmartAssembly_hh





