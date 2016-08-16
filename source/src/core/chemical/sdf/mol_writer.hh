// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/sdf/mol_writer.hh
/// @author Sam DeLuca

#ifndef INCLUDED_core_chemical_sdf_mol_writer_hh
#define INCLUDED_core_chemical_sdf_mol_writer_hh

#include <core/chemical/sdf/mol_writer.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <list>

//#include <utility/io/ozstream.fwd.hh>
#include <string>
#include <map>
namespace core {
namespace chemical {
namespace sdf {


class MolWriter
{
public:
	MolWriter();
	MolWriter(std::string const & ctab_mode);
	void output_residue(std::ostream & output_stream, core::conformation::Residue const & residue);
	void output_residue(std::ostream & output_stream, core::chemical::ResidueType const & residue_type);

	void output_residue(std::string const & file_name,core::conformation::Residue const & residue);
	void output_residue(std::string const & file_name, core::chemical::ResidueType const & residue_type);

	void output_residue(std::ostream & output_stream, core::conformation::ResidueCOP residue);
	void output_residue(std::ostream & output_stream, core::chemical::ResidueTypeCOP residue_type);

	void output_residue(std::string const & file_name,core::conformation::ResidueCOP residue);
	void output_residue(std::string const & file_name, core::chemical::ResidueTypeCOP residue_type);

	inline void set_job_data(std::map<std::string,std::string> const & job_data)
	{
		job_data_ = job_data;
	}

private:

	enum CtabMode {V2000,V3000};

	std::list<std::string> compose_metadata(core::conformation::Residue const & residue);
	std::list<std::string> compose_ctab(core::conformation::Residue const & residue);
	std::list<std::string> compose_atoms(core::conformation::Residue const & residue);
	std::list<std::string> compose_bonds(core::conformation::Residue const & residue);
	std::list<std::string> compose_properties(core::conformation::Residue const & residue);
	std::list<std::string> compose_typeinfo(core::conformation::Residue const & residue);
	std::list<std::string> compose_nbr_atom(core::conformation::Residue const & residue);
	std::list<std::string> compose_naming(core::conformation::Residue const & residue);
	std::list<std::string> compose_rosetta_properties(core::conformation::Residue const &  residue);
	std::list<std::string> compose_job_info();

	//utility::io::ozstream output_stream_;
	std::string const line_header_;
	std::map<std::string,std::string> job_data_;
	CtabMode ctab_mode_;
};

template< class Input_t >
void
output_residue(std::ostream & output_stream, Input_t residue ) {
	MolWriter writer;
	writer.output_residue( output_stream, residue );
}

}
}
}

#endif /* MOL_WRITER_HH_ */
