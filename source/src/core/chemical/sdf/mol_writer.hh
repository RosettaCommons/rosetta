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
#include <core/chemical/MutableResidueType.fwd.hh>
#include <utility/excn/Exceptions.hh>
#include <list>

#include <utility/io/ozstream.hh>
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

	void output_residue(std::ostream & output_stream, core::conformation::ResidueCOP residue);
	void output_residue(std::ostream & output_stream, core::chemical::ResidueTypeCOP residue_type);
	void output_residue(std::ostream & output_stream, core::chemical::MutableResidueTypeCOP residue_type);

	void output_residue(std::ostream & output_stream, core::conformation::Residue const & residue);
	void output_residue(std::ostream & output_stream, core::chemical::ResidueType const & residue_type);
	void output_residue(std::ostream & output_stream, core::chemical::MutableResidueType const & residue_type);

	template< class Input_t >
	void output_residue(std::string const & file_name, Input_t & residue ) {
		utility::io::ozstream outfile;
		outfile.open(file_name.c_str(), std::ios::out | std::ios::binary);
		if ( !outfile ) {
			throw CREATE_EXCEPTION(utility::excn::FileNotFound, "Cannot open file"+file_name);
		}
		output_residue(outfile,residue);
		outfile.close();
	}

	inline void set_job_data(std::map<std::string,std::string> const & job_data)
	{
		job_data_ = job_data;
	}

private:

	enum CtabMode {V2000,V3000};

	std::ostream & open_file( std::string const & file_name ) const;

	void output_residue_impl(std::ostream & output_stream, core::chemical::MutableResidueType const & residue_type, std::map< std::string, core::Vector > const & coords = {} );

	std::list<std::string> compose_metadata(core::chemical::MutableResidueType const & residue);
	std::list<std::string> compose_ctab(core::chemical::MutableResidueType const & residue, std::map< std::string, core::Vector > const & coords);
	std::list<std::string> compose_atoms(core::chemical::MutableResidueType const & residue, std::map< std::string, core::Vector > const & coords);
	std::list<std::string> compose_bonds(core::chemical::MutableResidueType const & residue);
	std::list<std::string> compose_properties(core::chemical::MutableResidueType const & residue);
	std::list<std::string> compose_atomnames(core::chemical::MutableResidueType const & residue);
	std::list<std::string> compose_typeinfo(core::chemical::MutableResidueType const & residue);
	std::list<std::string> compose_nbr_atom(core::chemical::MutableResidueType const & residue);
	std::list<std::string> compose_naming(core::chemical::MutableResidueType const & residue);
	std::list<std::string> compose_rosetta_properties(core::chemical::MutableResidueType const &  residue);
	std::list<std::string> compose_job_info();

private:

	//utility::io::ozstream output_stream_;
	std::string const line_header_;
	std::map<std::string,std::string> job_data_;
	CtabMode ctab_mode_;
};

template< class Input_t >
void
output_residue(std::string const & outfile, Input_t & residue ) {
	MolWriter writer;
	writer.output_residue( outfile, residue );
}

template< class Input_t >
void
output_residue(std::ostream & output_stream, Input_t & residue ) {
	MolWriter writer;
	writer.output_residue( output_stream, residue );
}

}
}
}

#endif /* MOL_WRITER_HH_ */
