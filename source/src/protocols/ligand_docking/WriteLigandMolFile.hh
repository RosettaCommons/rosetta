// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/WriteLigandMolFile.hh
/// @brief  header class for Ligand MolFile writer
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_ligand_docking_WriteLigandMolFile_hh
#define INCLUDED_protocols_ligand_docking_WriteLigandMolFile_hh

#include <protocols/moves/Mover.hh>
#include <protocols/ligand_docking/WriteLigandMolFile.fwd.hh>

namespace protocols {
namespace ligand_docking {

class WriteLigandMolFile: public protocols::moves::Mover
{
public:
	WriteLigandMolFile();
	virtual ~WriteLigandMolFile();
	WriteLigandMolFile(WriteLigandMolFile const & that);

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	virtual void apply(core::pose::Pose & pose);

private:
	std::string chain_;
	std::string directory_;
    std::string prefix_;
	bool hash_file_names_;
};

}
}


#endif /* WRITELIGANDMOLFILE_HH_ */
