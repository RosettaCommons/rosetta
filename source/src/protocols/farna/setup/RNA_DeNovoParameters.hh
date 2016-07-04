// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/setup/RNA_DeNovoParameters.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_DeNovoParameters_HH
#define INCLUDED_protocols_farna_RNA_DeNovoParameters_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/farna/setup/RNA_DeNovoPoseInitializer.fwd.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.fwd.hh>
#include <core/pose/rna/BasePair.hh>

namespace protocols {
namespace farna {
namespace setup {

class RNA_DeNovoParameters: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_DeNovoParameters( std::string const & filename );

	RNA_DeNovoParameters();

	//destructor
	~RNA_DeNovoParameters();

public:

	friend class RNA_DeNovoPoseInitializer;


	void set_rna_pairing_list( core::pose::rna::RNA_BasePairList const & setting ) { rna_pairing_list_ = setting; }
	core::pose::rna::RNA_BasePairList const & rna_pairing_list() const {
		return rna_pairing_list_;
	}

	void set_stem_pairing_sets( utility::vector1 < utility::vector1 <core::Size > > const & setting ) { stem_pairing_sets_ = setting; }
	utility::vector1 < utility::vector1 <core::Size > > const & stem_pairing_sets() const {
		return stem_pairing_sets_;
	}

	void set_obligate_pairing_sets( utility::vector1 < utility::vector1 <core::Size > > const & setting ) { obligate_pairing_sets_ = setting; }
	utility::vector1 < utility::vector1 <core::Size > > const & obligate_pairing_sets() const {
		return obligate_pairing_sets_;
	}

	void set_cutpoints_open( utility::vector1 <core::Size > const & setting ){ cutpoints_open_ = setting; }
	utility::vector1 <core::Size > cutpoints_open() const { return cutpoints_open_; }

	void set_cutpoints_closed( utility::vector1 <core::Size > const & setting ){ cutpoints_closed_ = setting; }
	utility::vector1 <core::Size > cutpoints_closed() const { return cutpoints_closed_; }

	void set_virtual_anchor_attachment_points( utility::vector1 <core::Size > const & setting ){ virtual_anchor_attachment_points_ = setting; }
	utility::vector1 <core::Size > virtual_anchor_attachment_points() const { return virtual_anchor_attachment_points_; }

	void set_allow_insert_res( utility::vector1 <core::Size > const & setting ){ allow_insert_res_ = setting; }
	utility::vector1 <core::Size > allow_insert_res() const { return allow_insert_res_; }

	void set_chain_connections( 	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > >  const & setting ) { chain_connections_ = setting; }
	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > >  const &
	chain_connections() const { return chain_connections_; }

private:

	void
	read_parameters_from_file( std::string const & filename );

	void
	save_res_lists_to_chain_connections_and_clear( utility::vector1< Size > & res_list1,
		utility::vector1< Size > & res_list2 );


	void
	read_chain_connection( std::istringstream & line_stream );

	void
	get_pairings_from_line(
		std::istringstream & line_stream,
		bool const in_stem );

	Size
	check_in_pairing_sets( utility::vector1 < utility::vector1 <core::Size > > pairing_sets,
		core::pose::rna::BasePair const & rna_pairing_check ) const;

private:

	std::string const filename_;

	core::pose::rna::RNA_BasePairList rna_pairing_list_;

	utility::vector1 < utility::vector1 <core::Size > > obligate_pairing_sets_;
	utility::vector1 < utility::vector1 <core::Size > > stem_pairing_sets_;

	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > > chain_connections_;

	utility::vector1 <core::Size > cutpoints_open_;
	utility::vector1 <core::Size > cutpoints_closed_;
	utility::vector1 <core::Size > virtual_anchor_attachment_points_;

	utility::vector1 < core::Size > allow_insert_res_;
	std::string rna_secstruct_legacy_;
	bool secstruct_defined_;

};

} //setup
} //farna
} //protocols

#endif
