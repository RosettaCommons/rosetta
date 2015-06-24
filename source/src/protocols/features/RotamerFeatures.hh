// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/RotamerFeatures.hh
/// @brief  report idealized torsional DOFs Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_RotamerFeatures_hh
#define INCLUDED_protocols_features_RotamerFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/RotamerFeatures.fwd.hh>
#include <basic/database/schema_generator/Schema.hh>


#include <protocols/features/RotamerFeatures.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

//External

// Project Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols{
namespace features{

//@brief Extract from the dunbrack Energy term the model for the
// rotamer conformation.
template < core::Size T, core::Size N >
class RotamerInitializer {

public:
	static
	bool
	initialize_rotamer(
		core::conformation::Residue const & residue,
		core::pack::dunbrack::RotamerLibraryScratchSpace & scratch,
		core::Size & rotamer_bin
	) {
		using namespace core::pack::dunbrack;
		using namespace core::pack::rotamers;

		SingleResidueRotamerLibraryCOP generic_rotlib =
			core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( residue.type() );

		if(!generic_rotlib){
			return false;
		}

		// This will throw a std::bad_cast if the residue has a non Dunbrack rotamer library.
		RotamericSingleResidueDunbrackLibrary< T, N > const & rotlib(
			dynamic_cast< RotamericSingleResidueDunbrackLibrary< T, N > const & >(
				* generic_rotlib));

		RotVector rotamer_vector;
		Size4 rotamer_fixed_vector;
		core::Size packed_rotno;

		// can't use get_rotamer_from_chi_static because it's private. Perhaps it should be made public?
		rotlib.get_rotamer_from_chi(residue.chi(), rotamer_vector);

		if(rotamer_vector.size() > 4){
			//eg LYS_p:dimethylated, perhaps there is a more direct way to
			//detect cases like this?
			return false;
		}

		copy(rotamer_vector.begin(), rotamer_vector.end(), rotamer_fixed_vector.begin());

		packed_rotno = rotlib.rotwell_2_packed_rotno(rotamer_vector);
		if(packed_rotno == 0){
			packed_rotno = rotlib.find_another_representative_for_unlikely_rotamer(
				residue, rotamer_fixed_vector);
			rotlib.packed_rotno_2_rotwell(packed_rotno, rotamer_vector);
		}
		rotamer_bin = rotlib.rotwell_2_rotno(rotamer_vector);

		PackedDunbrackRotamer< T, N, core::Real > interpolated_rotamer;
		rotlib.interpolate_rotamers(
			residue, scratch, packed_rotno, interpolated_rotamer);

		return true;

	}

};


class RotamerFeatures : public protocols::features::FeaturesReporter {
public:
	RotamerFeatures(){}

	RotamerFeatures(
		RotamerFeatures const & ) :
		FeaturesReporter()
	{}

	virtual ~RotamerFeatures(){}

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the residue_rotamers table schema
	void
	write_residue_rotamers_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);
};

} // namespace
} // namespace

#endif // include guard
