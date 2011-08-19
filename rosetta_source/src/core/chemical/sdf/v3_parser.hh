// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /grid/src/core/chemical/sdf/v3_parser.hh
/// @author Sam DeLuca

#ifndef INCLUDED_core_chemical_sdf_v3_parser_hh
#define INCLUDED_core_chemical_sdf_v3_parser_hh

#include <core/chemical/sdf/v3_parser.fwd.hh>
#include <core/chemical/sdf/ctab_base.hh>

namespace core {
namespace chemical {
namespace sdf{

class V3Parser :public CtabBase
{
public:
	V3Parser(utility::vector1<std::string> const & connection_table_lines, core::chemical::ResidueTypeOP molecule_container, MolData const & mol_data);
	virtual void ParseTable();


private:
	virtual void ParseAtom(std::string const atom_line, core::Size const );
	virtual void ParseBond(std::string const bond_line);
	core::Real FindExtraParameter( std::vector<std::string> const extra_parameters,  std::string const query);

};

}
}
}

#endif /* V3_PARSER_HH_ */
