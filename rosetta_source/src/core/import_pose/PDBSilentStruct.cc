// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/PDBSilentStruct.cc
///
/// @brief Representation of PDB files in a silent-file format.
/// @author James Thompson

// C++ Headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <sstream>

// mini headers
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/import_pose/PDBSilentStruct.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/Pose.hh>

#include <core/io/pdb/pose_io.hh>
// Auto-header: duplicate removed #include <core/io/pdb/pdb_dynamic_reader.hh>

#include <core/conformation/Residue.fwd.hh>

// ObjexxFCL
#include <ObjexxFCL/FArray2D.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <execinfo.h>
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <ostream>
#include <set>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end

namespace core {
namespace import_pose {

using namespace ObjexxFCL;

static basic::Tracer tr("core.io.silent.PDBSilentStruct");

PDBSilentStruct::PDBSilentStruct(
	core::pose::Pose const & pose,
	std::string tag
) {
	fill_struct( pose, tag );
} // PDBSilentStruct

void PDBSilentStruct::print_header( std::ostream & out ) {
	print_score_header( out );
}

void PDBSilentStruct::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {
	decoy_tag( tag );
	if ( tag == "empty_tag" ) set_tag_from_pose( pose );

	energies_from_pose( pose );
	fd_.init_from_pose( pose );

	sequence( pose.sequence() );
}

bool PDBSilentStruct::init_from_lines(
	utility::vector1< std::string > const & lines,
	core::io::silent::SilentFileData & container
) {
	bool success( false );
	using std::string;
	using utility::vector1;
	using namespace core::io::silent;

	vector1< std::string > energy_names_;
	vector1< std::string >::const_iterator iter = lines.begin();
	if ( iter->substr(0,9) == "SEQUENCE:" ) iter++; // ignore sequence for now
	if ( iter->substr(0,6) != "SCORE:" ) {
		// get sequence and scorename data from the silent-file data object, because I don't have it!
		EnergyNamesOP enames = EnergyNamesOP(
			static_cast< EnergyNames * > ( container.get_shared_silent_data( energynames )() )
		);

		energy_names_ = enames->energy_names();
	} else {
		// get scorename data from the first two lines provided, put into container
		// for further use by other SilentStruct objects.

		EnergyNamesOP enames( new EnergyNames( *iter ) );
		container.set_shared_silent_data( energynames, enames  );
		energy_names_ = enames->energy_names();
	} // get header information

	std::string concatenated_pdb_info; // concatenated pdb information
	for ( vector1< string >::const_iterator end = lines.end(); iter != end;	++iter ) {
		string tag;
		std::istringstream line_stream( *iter );

		if ( iter->substr(0,7) == "SCORE: " ) { // SCORE: line with values from this structure.
			std::string tag;
			line_stream >> tag;
			if ( line_stream.fail() || tag != "SCORE:" ) {
				tr.Error << "bad format in first score line of silent file" << std::endl;
				tr.Error << "line = " << *iter << std::endl;
				tr.Error << "tag = " << tag << std::endl;
			}

			vector1< string >::const_iterator energy_iter;
			for ( energy_iter = energy_names_.begin();
						energy_iter != energy_names_.end(); ++energy_iter
			) {
				line_stream >> tag;
				if ( *energy_iter != "description" ) { // currently the only text-based field, might change in future.
					Real score_val = (Real) float_of( tag );
					add_energy( *energy_iter, score_val );
				} else {
					line_stream >> tag;
				}
			} // for ( energy_iter ... )
			decoy_tag( tag ); // decoy_tag should be last column of this line.
		} else if ( iter->substr(0,6) == "REMARK" ) {
			// do nothing!
		} else {
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			if ( option[ out::file::silent_preserve_H ]() || iter->substr(13,1) != "H" ) {
				// FileData needs \n's
				concatenated_pdb_info += iter->substr( 0, 79 ) + '\n';
				string temp_decoy_tag  = iter->substr( 80 );
			}
		} // else
	} // for ( iter ... )

	pdb_lines_ = concatenated_pdb_info;
	fd_ = core::io::pdb::PDB_DReader::createFileData( concatenated_pdb_info );

	success = true;
	return success;
} // init_from_lines

void PDBSilentStruct::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set
) const {
	core::import_pose::build_pose( fd_, pose, residue_set );
	core::import_pose::read_additional_pdb_data( pdb_lines_, pose, fd_ );
} // fill_pose

void PDBSilentStruct::fill_pose(
	core::pose::Pose & pose
) const {

	using namespace core::chemical;
	ResidueTypeSetCAP residue_set;
	residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	fill_pose( pose, *residue_set );
	finish_pose( pose );
} // fill_pose


void PDBSilentStruct::print_conformation( std::ostream & output ) const {
	using std::string;

	string data = core::io::pdb::PDB_DReader::createPDBData(fd_);
	output.write( data.c_str(), data.size() );
} // print_conformation

Real PDBSilentStruct::get_debug_rmsd() {
	tr.Error << "get_debug_rmsd stubbed out!" << std::endl;
	return 0.0;
}

ObjexxFCL::FArray2D< Real >
PDBSilentStruct::get_CA_xyz() const {
	tr.Error << "PDBSilentStruct::get_CA_xyz" << std::endl;
	return FArray2D< Real > ( 3, 1 );
}

PDBSilentStruct & PDBSilentStruct::operator= (
	PDBSilentStruct const &
)
{
	utility_exit_with_message( "called ProteinSilentStruct::operator=)" );
	exit(0); //  just to keep the compiler happy
}

} // namespace import_pose
} // namespace core
