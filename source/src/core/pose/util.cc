// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/util.cc
/// @brief  Pose class utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Rhiju Das, Steven Lewis, Vikram K. Mulligan


// Unit header
#include <core/pose/util.hh>

// Package headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/PositionConservedResiduesStore.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/rna/util.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/Exceptions.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.string.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
//#include <stdio.h>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <boost/functional/hash.hpp>
#include <boost/foreach.hpp>


namespace core {
namespace pose {

static basic::Tracer TR("core.pose.util");

void
append_pose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	bool new_chain
){
	append_subpose_to_pose(pose1, pose2, 1, pose2.total_residue(), new_chain);
}

void
append_subpose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::Size start_res,
	core::Size end_res,
	bool new_chain
){
	if(pose2.total_residue()>=start_res){
		pose1.append_residue_by_jump(pose2.residue(start_res), pose1.total_residue() , "", "", new_chain);
		for(core::Size i=start_res+1; i<=end_res; ++i){
			if(pose2.residue(i).is_lower_terminus()){
				if(i > 1 && pose2.chain(i) == pose2.chain(i-1)){
					pose1.append_residue_by_jump(pose2.residue(i), pose1.total_residue(), "","", false);
				}
				else{
					pose1.append_residue_by_jump(pose2.residue(i), pose1.total_residue(), "","", true);
				}
			}
			else{
				pose1.append_residue_by_bond(pose2.residue(i));
			}
		}
	}
}

void jumps_from_pose(const core::pose::Pose& pose, Jumps* jumps) {
	assert(jumps);
	for (Size i = 1; i <= pose.num_jump(); ++i) {
		jumps->insert(i);
	}
}

void remove_virtual_residues(core::pose::Pose* pose) {
  assert(pose);
  for (core::Size i = 1; i <= pose->total_residue(); ++i) {
    if (pose->residue_type(i).name() == "VRT")
      pose->conformation().delete_residue_slow(i);
  }
}

void swap_transform(Size jump_num, const kinematics::RT& xform, Pose* pose) {
  assert(pose);
  assert(jump_num <= pose->num_jump());

  const kinematics::FoldTree& tree = pose->fold_tree();

  int upstream = tree.upstream_jump_residue(jump_num);
  int downstream = tree.downstream_jump_residue(jump_num);
  const conformation::Residue& upstream_res = pose->residue(upstream);
  const conformation::Residue& downstream_res = pose->residue(downstream);

  if (upstream_res.natoms() < 3 || downstream_res.natoms() < 3) {
    TR.Warning << "Insufficient number of atoms for stub creation on one or more"
               << " jump residues-- " << upstream << ", " << downstream
               << std::endl;
    return;
  }

  id::StubID upstream_stub(
      id::AtomID(0, upstream_res.atom_index("N")),
      id::AtomID(1, upstream_res.atom_index("CA")),
      id::AtomID(2, upstream_res.atom_index("C")));

  id::StubID downstream_stub(
      id::AtomID(0, downstream_res.atom_index("N")),
      id::AtomID(1, downstream_res.atom_index("CA")),
      id::AtomID(2, downstream_res.atom_index("C")));

  assert(upstream_stub.valid());
  assert(downstream_stub.valid());

  pose->conformation().set_stub_transform(
      upstream_stub,
      downstream_stub,
      xform);
}

bool is_position_conserved_residue(const Pose& pose, core::Size residue) {
  using basic::datacache::BasicDataCache;
  using core::pose::datacache::PositionConservedResiduesStore;
  using core::pose::datacache::PositionConservedResiduesStoreCOP;

  assert(residue > 0);
  assert(residue <= pose.total_residue());

  const BasicDataCache& cache = pose.data();
  if (!cache.has(core::pose::datacache::CacheableDataType::POSITION_CONSERVED_RESIDUES))
    return false;

  PositionConservedResiduesStoreCOP store =
      static_cast<PositionConservedResiduesStore const *>(
          cache.get_const_ptr(core::pose::datacache::CacheableDataType::POSITION_CONSERVED_RESIDUES)());

  return store->is_conserved(residue);
}

void
create_subpose(
	Pose const & src,
	utility::vector1< Size > const & positions,
	kinematics::FoldTree const & f,
	Pose & pose
)
{
	Size const nres( f.nres() );
	assert( nres == positions.size() );

	pose.clear();

	for ( Size i=1; i<= nres; ++i ) {
		Size const seqpos( positions[i] );
		conformation::Residue const & rsd( src.residue( seqpos ) );
		// If the residue and the previous residue are bonded in the source pose, they should be bonded in the new pose
		if( i>1 && rsd.is_polymer_bonded( positions[ i-1 ] ) ) {
			pose.append_residue_by_bond( rsd );
		} else {
			pose.append_residue_by_jump( rsd, 1 );
		}
		if ( i>1 ) {
			// check if this residue should be in a new chain. not a perfect check...
			conformation::Residue const & prev_rsd( src.residue( positions[i-1] ) );
			if ( prev_rsd.is_upper_terminus() || rsd.is_lower_terminus() || prev_rsd.chain() != rsd.chain() ) {
				assert( pose.total_residue() == i );
				pose.conformation().insert_chain_ending( i-1 );
			}
		}
	}

	// now set the desired foldtree
	pose.fold_tree(f);

}

////////////////////////////////////////////////////////////////////////////
void
partition_pose_by_jump(
	pose::Pose const & src,
	int const jump_number,
	pose::Pose & partner1,
	pose::Pose & partner2
)
{
	Size const nres( src.total_residue() );

	// split src pose's foldtree
	kinematics::FoldTree f1,f2;
	src.fold_tree().partition_by_jump( jump_number, f1, f2 );

	TR << src.fold_tree() << '\n' << f1 << '\n' << f2 << '\n';

	// identify residues in the two partners
	ObjexxFCL::FArray1D_bool partner1_pos( nres, false ); // FARRAY! DOH!!
	src.fold_tree().partition_by_jump( jump_number, partner1_pos );

	utility::vector1< Size > partner1_pos_list, partner2_pos_list;
	for ( Size i=1; i<= nres; ++i ) {
		if ( partner1_pos(i) ) partner1_pos_list.push_back( i );
		else partner2_pos_list.push_back( i );
	}

	create_subpose( src, partner1_pos_list, f1, partner1 );

	create_subpose( src, partner2_pos_list, f2, partner2 );

/// SJF 11Sep13 why dump the poses?
//	src.dump_pdb( "complex.pdb" );
//	partner1.dump_pdb( "partner1.pdb" );
//	partner2.dump_pdb( "partner2.pdb" );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Crude way to guess secondary structure given a pose. This function
/// sets the sec_struct array of pose.conformation_ to the result of the
/// guesswork. This has been ported directly from rosetta++.
void
set_ss_from_phipsi(
	pose::Pose & pose
)
{
	// ss       :ss = 1 helix, ss = 2 sheet, ss = 3 other
	const int sstemp_offset=3;
	utility::vector1 < int > sstemp( sstemp_offset*2 + pose.total_residue() );
	utility::vector1 < int > ss( pose.total_residue() );

	sstemp[sstemp_offset-1] = 3; // assign loop to fictious residues at ends of chain
	sstemp[sstemp_offset+0] = 3;
	sstemp[sstemp_offset+pose.total_residue()+1] = 3;
	sstemp[sstemp_offset+pose.total_residue()+2] = 3;

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {

	// <murphp>
		if( !pose.residue_type(i).is_protein() ) { // make sure we don't inquire about the phi/psi of a non-protein residue
			sstemp[sstemp_offset+i] = 3;
		}
		else {
			// </murphp>

			if ( pose.phi(i) < -20.0 && pose.psi(i) > -90.0 && pose.psi(i) < -10.0 ) {
				sstemp[sstemp_offset+i] = 1;
			} else if ( pose.phi(i) < -20.0 && (pose.psi(i) > 20.0 || pose.psi(i) < -170.0) ) {
				sstemp[sstemp_offset+i] = 2;
			} else {
				sstemp[sstemp_offset+i] = 3;
			}

		}
	}

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( sstemp[sstemp_offset+i] == 2 ) {
			if ( sstemp[sstemp_offset+i-1] == 2 && sstemp[sstemp_offset+i+1] == 2 ) {
				ss[i] = 2;
			} else if ( sstemp[sstemp_offset+i-1] == 2 && sstemp[sstemp_offset+i-2] == 2 ) {
				ss[i] = 2;
			} else if ( sstemp[sstemp_offset+i+1] == 2 && sstemp[sstemp_offset+i+2] == 2 ) {
				ss[i] = 2;
			} else {
				ss[i] = 3;
			}
		} else if ( sstemp[sstemp_offset+i] == 1 ) {
			if ( sstemp[sstemp_offset+i-1] == 1 && sstemp[sstemp_offset+i+1] == 1 ) {
				ss[i] = 1;
			} else if ( sstemp[sstemp_offset+i-1] == 1 && sstemp[sstemp_offset+i-2] == 1 ) {
				ss[i] = 1;
			} else if ( sstemp[sstemp_offset+i+1] == 1 && sstemp[sstemp_offset+i+2] == 1 ) {
				ss[i] = 1;
			} else {
				ss[i] = 3;
			}
		} else {
			ss[i] = 3;
		}
	}

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( ss[i] == 1 ) {
			pose.set_secstruct(i,'H');
		} else if ( ss[i] == 2 ) {
			pose.set_secstruct(i,'E');
		} else {
			pose.set_secstruct(i,'L');
		}
	}
}

/// @details Adds a virtual residue to a pose as the root. Jump is to the
/// residue closest to <xyz>. If the pose is already rooted on a VRT res,
/// do nothing.
void addVirtualResAsRoot(const numeric::xyzVector<core::Real>& xyz, core::pose::Pose& pose) {
	int nres = pose.total_residue();

	// if already rooted on virtual residue, return
	if ( pose.residue( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		TR.Warning << "addVirtualResAsRoot() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return;
	}

	// return if the pose is empty (otherwise will segfault)
	if (nres == 0) {
		TR.Warning << "WARNING: addVirtualResAsRoot() called with empty pose!" << std::endl;
		return;
	}

	// check for terminal ligands
	int last_peptide_res = nres;
	while ( !pose.residue( last_peptide_res ).is_polymer() )
		last_peptide_res--;

	// try to avoid putting the vrt too close to termini
	int i_min = 1;

#ifdef WIN32
	int r_start = static_cast< int > ( std::floor(   static_cast< double > (last_peptide_res) /3. ) );
	int r_end   = static_cast< int > ( std::ceil ( 2.* static_cast< double > (last_peptide_res)/3. ) );
#else
	int r_start = static_cast< int > ( std::floor(   last_peptide_res/3 ) );
	int r_end   = static_cast< int > ( std::ceil ( 2*last_peptide_res/3 ) );
#endif

	// If less than three total residues, reset starting residue
	if (r_start == 0) {
		r_start = 1;
	}

	core::Real d_min = 99999, this_d;
	for ( int i=r_start; i<=r_end; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );

		if (rsd.aa() == core::chemical::aa_vrt) continue;
		if (!rsd.is_protein() ) continue;

		conformation::Atom const & atom( rsd.atom("CA") );
		this_d = (atom.xyz() - xyz).length();
		if (this_d < d_min) {
			d_min = this_d;
			i_min = i;
		}
	}
	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSetCAP const &residue_set(
								          core::chemical::ChemicalManager::get_instance()->residue_type_set
	                            ( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
	                                                    );
	core::chemical::ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("VRT") );
	if (rsd_type_list.size() == 0) {
		utility_exit_with_message("Cannot find residue type VRT" );
	}
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type_list[1] ) );

	// move to <xyz>
	for ( Size j=1; j<= new_res->natoms(); ++j ) {
		new_res->atom(j).xyz( new_res->atom(j).xyz()+xyz );
	}

	pose.append_residue_by_jump( *new_res , i_min );

	// update PDBinfo
	if (pose.pdb_info()) {
		pose.pdb_info()->chain( pose.total_residue(), 'z' );
		pose.pdb_info()->number( pose.total_residue(), 1 );
		pose.pdb_info()->obsolete( false );
	}

	// make the virt atom the root
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder( nres+1 );
	TR.Debug << "addVirtualResAsRoot() setting new fold tree to " << newF << std::endl;
	TR.Debug << "   i_min = " << i_min << "   d_min = " << d_min << std::endl;
	pose.fold_tree( newF );
}

/// @detail Get center of mass of a pose.
numeric::xyzVector< core::Real >
get_center_of_mass( core::pose::Pose const & pose ){
	numeric::xyzVector< core::Real > massSum( 0.0 );
	Size const & nres = pose.total_residue();
	Size nAtms=0;
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if (rsd.aa() == core::chemical::aa_vrt) continue;
		for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );
			massSum += atom.xyz();
			nAtms++;
		}
	}
	massSum /= nAtms;
	return massSum;
}

/// @detail Find residue closest to center-of-mass
void addVirtualResAsRoot( core::pose::Pose & pose ) {
	int nres = pose.total_residue();
	// return if the pose is empty (otherwise will segfault)
	if (nres == 0) {
		TR.Warning << "WARNING: addVirtualResAsRoot() called with empty pose!" << std::endl;
		return;
	}
	numeric::xyzVector< core::Real > massSum = get_center_of_mass( pose );
	return addVirtualResAsRoot(massSum, pose);
}

/// @brief  Reads the comments from the pdb file and adds it into comments
void
read_comment_pdb(
	std::string const &file_name,
	core::pose::Pose & pose
) {
	utility::io::izstream data(file_name);
  if ( !data ) {
		TR<<"ERROR! PDB SCAFOLD NAME NOT FOUND!!"<<std::endl;
		utility_exit();
  }
  std::string line;
  while( getline( data, line ) ) {
  	if( line != "##Begin comments##")
			continue;
			getline( data, line );
		while (line != "##End comments##"){
			//TR<<"Testing read comments! :"<<line<<std::endl;
			std::string const key;
			std::string const value;
			utility::vector1<std::string> comment_line(utility::string_split(line,' '));
			core::pose::add_comment(pose,comment_line[1],comment_line[2]);
			getline( data, line );
	}
}
}

void
dump_comment_pdb(
	std::string const &file_name,
	core::pose::Pose const& pose
) {
	std::ofstream out( file_name.c_str() );
	pose.dump_pdb(out);
	// verbose output
	out << "END\n";
	  out << "##Begin comments##" << std::endl;
	 	                        using namespace std;
	 	                        map< string, string > const comments = core::pose::get_all_comments(pose);
	 	                        for( std::map< string, string >::const_iterator i = comments.begin(); i != comments.end(); ++i ){
	 	                                out << i->first<<" "<<i->second << std::endl;
	 	                        }
	 	                        out << "##End comments##" << std::endl;
	out.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool getPoseExtraScores(
	core::pose::Pose const & pose,
	std::string const name,
	core::Real & value
)
{
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapCOP;

	// make sure that the pose has one of these.
	if( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ){
		return false;
	}

	CacheableStringFloatMapCOP data
		= dynamic_cast< CacheableStringFloatMap const * >
			( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
	assert( data.get() != NULL );

	std::map< std::string, float >::const_iterator it = data->map().find( name );
	if ( it == data->map().end() ) {
		return false;
	}
	value = it->second;
	return true;
}

bool getPoseExtraScores(
	core::pose::Pose const & pose,
	std::string const name,
	std::string & value
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	// make sure that the pose has one of these.
	if( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ){
		return false;
	}

	CacheableStringMapCOP data
		= dynamic_cast< CacheableStringMap const * >
			( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
	assert( data.get() != NULL );

	std::map< std::string, std::string >::const_iterator it = data->map().find( name );
	if ( it == data->map().end() ) {
		return false;
	}
	value = it->second;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void setPoseExtraScores(
	core::pose::Pose & pose,
	std::string name,
	core::Real value
)
{
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA,
			new CacheableStringFloatMap()
		);
	}

	CacheableStringFloatMapOP data
		=  dynamic_cast< CacheableStringFloatMap* >
			( pose.data().get_raw_ptr(core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA) );

	runtime_assert( data.get() != NULL );
	data->map()[name] = value;
}

void setPoseExtraScores(
	core::pose::Pose & pose,
	std::string const name,
	std::string value
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA,
			new CacheableStringMap()
		);
	}

	CacheableStringMapOP data
		=  dynamic_cast< CacheableStringMap* >
			( pose.data().get_raw_ptr(core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA) );

	runtime_assert( data.get() != NULL );
	data->map()[name] = value;
}

void add_comment(
	core::pose::Pose & pose,
	std::string const & key,
	std::string const & val
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;

	// make sure that the pose has a map of strings
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::STRING_MAP,
			new CacheableStringMap()
		);
	}

	CacheableStringMapOP data
		=  dynamic_cast< CacheableStringMap* >
			( pose.data().get_raw_ptr(core::pose::datacache::CacheableDataType::STRING_MAP) );

	runtime_assert( data.get() != NULL );
	data->map()[key] = val;
} // add_comment

void add_score_line_string(
	core::pose::Pose & pose,
	std::string const & key,
	std::string const & val
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;

	// make sure that the pose has a map of strings
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS,
			new CacheableStringMap()
		);
	}

	CacheableStringMapOP data
		=  dynamic_cast< CacheableStringMap* >
			( pose.data().get_raw_ptr(core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS) );

	runtime_assert( data.get() != NULL );
	data->map()[key] = val;
}

void clearPoseExtraScores(
	core::pose::Pose & pose
) {

	{
		using basic::datacache::CacheableStringFloatMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA,
			new CacheableStringFloatMap()
		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA,
			new CacheableStringMap()
		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::STRING_MAP,
			new CacheableStringMap()
		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS,
			new CacheableStringMap()
		);
	}
}

void clearPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name
)
{
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapOP;
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;

	if( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ){
		CacheableStringFloatMapOP data
			= dynamic_cast< CacheableStringFloatMap* >
				( pose.data().get_raw_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
		assert( data.get() != NULL );

		if ( data->map().find( name ) != data->map().end() ) data->map().erase( name );
	}

	if( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ){
		CacheableStringMapOP data
			= dynamic_cast< CacheableStringMap* >
				( pose.data().get_raw_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
		assert( data.get() != NULL );

		if ( data->map().find( name ) != data->map().end() ) data->map().erase( name );
	}
}

bool get_comment(
	core::pose::Pose const & pose,
	std::string const & key,
	std::string & val
) {
	using std::map;
	using std::string;
	map< string, string > comment_map = get_all_comments( pose );

	if ( comment_map.find( key ) == comment_map.end() ) {
		return false;
	}

	val = comment_map[ key ];
	return true;
}

bool get_score_line_string(
	core::pose::Pose const & pose,
	std::string const & key,
	std::string & val
) {
	using std::map;
	using std::string;
	map< string, string > score_line_strings_map = get_all_score_line_strings( pose );

	if ( score_line_strings_map.find( key ) == score_line_strings_map.end() ) {
		return false;
	}

	val = score_line_strings_map[ key ];
	return true;
}

void delete_comment(
	core::pose::Pose & pose,
	std::string const & key
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		CacheableStringMapOP data
			=  dynamic_cast< CacheableStringMap* >
				( pose.data().get_raw_ptr(core::pose::datacache::CacheableDataType::STRING_MAP) );
		std::map< std::string, std::string >::iterator it;
		it = data->map().find(key);
		if ( it != data->map().end() ) {
			data->map().erase(it);
		}
	}
}

std::map< std::string, std::string >
get_all_score_line_strings(
	core::pose::Pose const & pose
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	std::map< std::string, std::string > score_line_strings;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) ) {
		CacheableStringMapCOP data
			= dynamic_cast< CacheableStringMap const * >
				( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) );
		score_line_strings = data->map();
		runtime_assert( data.get() != NULL );
	}
	return score_line_strings;
}

std::map< std::string, std::string >
get_all_comments(
	core::pose::Pose const & pose
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	std::map< std::string, std::string > comments;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		CacheableStringMapCOP data
			= dynamic_cast< const CacheableStringMap* >
				( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::STRING_MAP ) );
		comments = data->map();
		runtime_assert( data.get() != NULL );
	}
	return comments;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< char > read_psipred_ss2_file( pose::Pose const & pose ) {
	using namespace basic::options;

	std::string filename( option[ OptionKeys::in::file::psipred_ss2 ]().name() );

	utility::vector1< char > secstructs;
	utility::io::izstream data( filename );

	if ( !data ) {
		TR.Warning << "Cannot open psipred_ss2 file " << filename << std::endl;
		return secstructs;
	}

	std::string line;
	Size count(0);
	while( getline( data, line ) ) {
		if( line[0] == '#' || line == "" )
			continue;
		std::istringstream line_stream( line );
		Size pos;
		char aa, sec;
		line_stream >> pos >> aa >> sec;
		count++;
		if ( sec != 'H' && sec != 'E' && sec != 'C' ) {
			TR.Warning << "unrecognized secstruct char : " << sec << " at seqpos " << count << std::endl;
		}
		if ( sec == 'C' ) {
			secstructs.push_back( 'L' );
		} else {
			secstructs.push_back( sec );
		}
	}

	// chu get number of protein residues
	Size nres=0;
	for ( Size i =1 ; i <= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() ) nres++;
	}

	assert( secstructs.size() == nres);
	if( secstructs.size() != nres )
		secstructs.clear();

	return secstructs;
}
////////////////////////////////////////////////////////////////////////////////////////////////////


void conf2pdb_chain_default_map( core::pose::Pose const & pose, std::map<int,char> & chainmap ) {
	chainmap.clear();
	char letter = 'A';
	for(core::Size i = 1; i <= pose.conformation().num_chains(); ++i){
		chainmap[i] = letter;
		if('Z'==letter) utility_exit_with_message("too many chains to map to letters!!!");
		letter = static_cast<char>(letter + 1);
	}
}

/// @brief get Conformation chain -> PDBInfo chain mapping
/// @remarks Any chains whose PDBInfo chain records are marked entirely as
///  PDBInfo::empty_record() will be mapped to that character.  Note that
///  Conformation -> PDBInfo is always unique, but the reverse may not be true.
/// @return the mapping if PDBInfo available and chains appear consistent,
///  otherwise returns an empty mapping
std::map< int, char > conf2pdb_chain( core::pose::Pose const & pose ) {
	using core::Size;
	using core::pose::PDBInfo;
	typedef std::map< int, char > Conf2PDB;

	Conf2PDB conf2pdb;

	if ( !pose.pdb_info().get() ) {
		TR.Warning << "WARNING: conf2pdb_chain(): PDBInfo does not exist, returning default map 1=A, 2=B, ..." << std::endl;
		conf2pdb_chain_default_map(pose,conf2pdb);
		return conf2pdb;
	}

	for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
		int const conf = pose.chain( i );
		char const pdb = pose.pdb_info()->chain( i );

		Conf2PDB::iterator c2p = conf2pdb.find( conf );
		if ( c2p != conf2pdb.end() ) { // must check if existing record inconsistent
			if ( c2p->second != pdb ) {

				// three cases:
				//  (1) replace an existing empty record
				//  (2) it's an unneeded empty record, so continue
				//  (3) there is an actual problem
				if ( pdb != PDBInfo::empty_record() && c2p->second == PDBInfo::empty_record() ) {
					// replace the record
					c2p->second = pdb;
				} else if ( pdb == PDBInfo::empty_record() ) {
					continue; // skip empty record
				} else {
					// something is inconsistent
					TR.Warning << "WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					TR.Warning << "WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HERE BE DRAGONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					TR.Warning << "WARNING: conf2pdb_chain(): chain mapping inconsistent, returning default map p 1=A, 2=B, ... ";
					TR.Warning << "WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					TR.Warning << "existing " << c2p->first << " -> " << c2p->second << "  |  ";
					TR.Warning << "new " << conf << " -> " << pdb << std::endl;
					conf2pdb_chain_default_map(pose,conf2pdb);
					return conf2pdb;
				}

			}
		} else { // record doesn't exist yet
			conf2pdb[ conf ] = pdb;
		}

	} // foreach residue

	assert( conf2pdb.size() == pose.conformation().num_chains() );
	return conf2pdb;
}


/// @brief renumber PDBInfo based on Conformation chains; each chain starts from 1
/// @param[in,out] pose The Pose to modify.
/// @param[in] fix_chains If true, the procedure will attempt to fix any empty record
///  characters it finds in the PDBInfo. (default true)
/// @param[in] start_from_existing_numbering If true, will attempt to start each
///  chain from the existing numbering in the PDBInfo.  E.g. if the first residue
///  of chain 2 in the Conformation is 27, then the renumbering of the chain in
///  PDBInfo will start from 27. (default true)
/// @param[in] keep_insertion_codes If true, will maintain insertion codes and
///  will not increment the pdb residue numbering for those residues.  This means
///  new numbering with insertion codes will only reflect properly if the
///  old numbering included the base numbering of the insertion code residues,
///  i.e. 100 100A 100B and not just 100A 100B (with 100 never appearing).
///  (default false)
/// @param[in] rotate_chain_ids If true, allows support for more than 26 pdb chains
///  by rotating [A,Z] continuously.  WARNING: This will break the assumption
///  made by the PDBPoseMap that each pdb chain id is unique, so make sure you
///  are not using the PDBPoseMap feature downstream in your code path without
///  corrections! (default false)
/// @remarks If fixing chains and there is only one chain and the PDBInfo exists
///  but all records are marked as empty, will renumber and set the PDBInfo chain
///  to 'A'.
/// @return true if renumbering successful, false otherwise
bool renumber_pdbinfo_based_on_conf_chains(
	core::pose::Pose & pose,
	bool fix_chains,
	bool const start_from_existing_numbering,
	bool const keep_insertion_codes,
	bool const rotate_chain_ids
)
{
	using core::Size;
	using core::pose::PDBInfo;
	typedef std::map< int, char > Conf2PDB;

	if ( !pose.pdb_info().get() ) {
		TR.Warning << "WARNING: renumber_pdbinfo_based_on_conf_chains(): no PDBInfo, returning" << std::endl;
		return false;
	}

	Conf2PDB conf2pdb = conf2pdb_chain( pose );

	if ( fix_chains ) {
		if ( conf2pdb.empty() ) { // something is wrong with chain consistency
			TR.Warning << "WARNING: renumber_pdbinfo_based_on_conf_chains(): Request to fix PDBInfo chains, but ";
			TR.Warning << "chain mapping is inconsistent, so that step will be skipped." << std::endl;
			fix_chains = false;
		} else { // Try to fill in any empty record characters.

			// two different schemes: rotating and fixed length
			// WARNING: Rotating will break assumption of unique chain ids
			// inside PDBPoseMap, so make sure you are not using the PDBPoseMap
			// feature after calling this function without correcting!
			// First either remove or rotate any existing chains to the end of
			// the list.
			std::string letters( "ABCDEFGHIJKLMNOPQRSTUVWXYZ" );
			for ( Conf2PDB::iterator i = conf2pdb.begin(), ie = conf2pdb.end(); i != ie; ++i ) {
				if ( i->second != PDBInfo::empty_record() ) {
					std::string::size_type const j = letters.find( i->second );
					if ( j != std::string::npos ) {
						if ( rotate_chain_ids ) { // rotating
							letters.push_back( letters.at( j ) );
						}
						letters.erase( j, 1 );
					}
				}
			}

			// Now fill in empty records.
			Size lidx = 0;
			for ( Conf2PDB::iterator i = conf2pdb.begin(), ie = conf2pdb.end(); i != ie; ++i ) {
				if ( i->second == PDBInfo::empty_record() ) {
					if ( rotate_chain_ids ) { // rotating
						i->second = letters.at( lidx % letters.size() );
					} else { // fixed length
						runtime_assert( lidx < letters.size() );
						i->second = letters.at( lidx );
					}
					++lidx;
				}
			}

		} // if conf2pdb.empty()
	} // if fix_chains

	PDBInfo & pdbinfo = *pose.pdb_info();

	// grab all the chain endings
	utility::vector1< Size > chain_endings = pose.conformation().chain_endings();
	chain_endings.push_back( pose.n_residue() ); // add the last res, which is not in the list

	Size res = 1;
	for ( utility::vector1< Size >::const_iterator i = chain_endings.begin(), ie = chain_endings.end(); i != ie; ++i ) {
		Size const chain_end = *i;
		int pdb_res = 0; // new chain, so reset pdb_res counter
		char chain;
		char icode;
		if ( start_from_existing_numbering && pdbinfo.chain( res ) != PDBInfo::empty_record() ) {
			pdb_res = pdbinfo.number( res ) - 1;
		}

		// find the corresponding pdb chain
		Conf2PDB::const_iterator c2p = conf2pdb.find( pose.chain( chain_end ) );
		assert( ( fix_chains && c2p != conf2pdb.end() ) || !fix_chains ); // otherwise something's very wrong

		for ( ; res <= chain_end; ++res ) {
			// handle the pdb chain only if necessary
			chain = pdbinfo.chain( res );
			if ( pdbinfo.chain( res ) == PDBInfo::empty_record() && fix_chains && c2p != conf2pdb.end() ) {
				chain = c2p->second;
			}

			// If keeping insertion codes, increment pdb_res counter only if
			// no insertion code or we're at position 1, in case there's an
			// insertion code at 1.
			icode = pdbinfo.icode( res );
			if ( keep_insertion_codes && ( pdbinfo.icode( res ) == ' ' || res == 1 ) ) {
				++pdb_res;
			} else if ( !keep_insertion_codes ) { // always increment and clear insertion code
				icode = ' ';
				++pdb_res;
			}

			// The new pdb info for this residue must be setup in one shot.
			// The way we're redoing the info in this function can cause the
			// pdb2pose map to become out-of-sync if we attempt to make the
			// changes by separately calling chain(), icode(), and number().
			pdbinfo.set_resinfo( res, chain, pdb_res, icode );
		}

	} // foreach chain ending

	// no point updating pdb_info if it's just thrown away
	pose.pdb_info()->obsolete( false );

	assert( res == pose.n_residue() + 1 );

	return true;
}

/// @brief checks if the pose geometry is ideal
/// @param[in] pose The Pose to check.
/// @return true if all pose positions have ideal bond lengths and angles
///  up to some very small epsilon
bool is_ideal_pose(
	core::pose::Pose const & pose
) {
	bool is_ideal=true;
	for ( core::Size i=1 ; i < pose.total_residue(); i++ ) {//leaving out last residue which always returns non-ideal for some reason
		if ( !is_ideal_position(i, pose) )	{
			is_ideal=false;
			break;
		}
	}
	return is_ideal;
}

/// @brief checks if the pose geometry is ideal in position seqpos
/// @param[in] pose The Pose to check.
/// @return true if position seqpos has ideal bond lengths and angles
///  up to some very small epsilon
bool is_ideal_position(
	core::Size seqpos,
	core::pose::Pose const & pose
)
{
	return conformation::is_ideal_position(
		seqpos,	pose.conformation()	);
}

/// @brief this function removes all residues from the pose which are not protein residues.
/// @details This removal includes, but is not limited to, metals, DNA, RNA, and ligands.
/// It will NOT remove ligands which are canonical residues (for example, if a protein binds an alanine monomer,
/// the monomer will be untouched).
void remove_nonprotein_residues( core::pose::Pose & pose )
{
	core::Size i(1);
	while(i <= pose.total_residue()) {
		if(!pose.residue_type(i).is_protein()) pose.conformation().delete_residue_slow(i);
		else ++i;
	}
}

/// @brief this function removes all residues with both UPPER and LOWER terminus types.
/// This is intended for removing ligands that are canonical residues.
void remove_ligand_canonical_residues( core::pose::Pose & pose )
{
	if(pose.total_residue() == 1) { //if we have only one residue, it cannot be removed, and this is going to crash
		throw utility::excn::EXCN_Msg_Exception("Pose utility remove_ligand_canonical_residues: I have received a pose with only one residue but cannot delete the last residue of the pose.");
	}

	core::Size i(1);
	while(i <= pose.total_residue()) {
		if(pose.residue_type(i).is_upper_terminus() && pose.residue_type(i).is_lower_terminus())
			pose.conformation().delete_residue_slow(i);
		else ++i;
	}
}


/// @details this function compares the 3-d coordinates of two poses.
/// Along the way it is forced to check for certain other (in)equalities to prevent vector overrruns, etc.
/// These include: pose length, ResidueType, and # atoms in residue.
/// Inequalities other than 3-d coordinates result in a warning message (you shouldn't have been comparing those poses!)
/// This is NOT a complete equality operator for a pose, but I think it does a good job with the coordinates.
/// Note that it performs floating-point coordinate comparisons by floor(X*10^n_dec_places) -
/// this may cause failures if your pose is a billion angstroms from 0,0,0.
/// This comparison is preferred to an epsilon comparison std::abs( a.x - b.x ) < epsilon because it can run into
/// situations where a == b and b == c, but a != c (thanks to APL for pointing this out).
/// The last argument, n_dec_places, is the number of decimal places of precision when comparing.
bool compare_atom_coordinates(core::pose::Pose const & lhs, core::pose::Pose const & rhs, core::Size const n_dec_places){

	//number of decimal places of precision - 3 (1000) is equivalent to PDB precision.
	core::Real const n_dec(pow(static_cast< Real > (10), static_cast< int > (n_dec_places))); //this is a Real to prevent premature rounding

	//first compare pose sizes; prerequisite to iterating through length
	core::Size const lhssize(lhs.total_residue()), rhssize(rhs.total_residue());

	if(lhssize != rhssize) {
		TR.Warning << "poses of different length in compare_atom_coordinates; doomed to fail!" << std::endl;
		return false;
	}

	//now iterate through residues and make comparisons
	for( core::Size i(1); i<=lhssize; ++i){
		//check equality of residue types
		core::chemical::ResidueType const & lhstype(lhs.residue_type(i)), & rhstype(rhs.residue_type(i));
		if(lhstype.name() != rhstype.name()) { //string matching is sufficient because ResidueType objects have unique names
			TR.Warning << "nonmatching ResidueTypes at " << i << " in compare_atom_coordinates" << std::endl;
			return false;
		}

		//get atoms vectors to compare
		core::conformation::Atoms const & lhsatoms(lhs.residue(i).atoms()), rhsatoms(rhs.residue(i).atoms());
		core::Size const lhsatmsize(lhsatoms.size()), rhsatmsize(rhsatoms.size());
		if(lhsatmsize != rhsatmsize) { //check vector length equality
			TR.Warning << "nonmatching numbers of atoms at residue " << i << " in compare_atom_coordinates" << std::endl;
			TR.Warning << "How did we even get here?  ResidueType comparison should have failed!" << std::endl;
			return false;
		}

		//iterate through atoms vector
		for( core::Size atm(1); atm <= lhsatmsize; ++atm){
			if(  (std::floor(lhsatoms[atm].xyz().x()*n_dec) != std::floor(rhsatoms[atm].xyz().x()*n_dec))
				|| (std::floor(lhsatoms[atm].xyz().y()*n_dec) != std::floor(rhsatoms[atm].xyz().y()*n_dec))
				|| (std::floor(lhsatoms[atm].xyz().z()*n_dec) != std::floor(rhsatoms[atm].xyz().z()*n_dec)) )
				return false; //no warning messages, this is the "expected" failure
		}//iterate through atoms vector
	}//over all residues

	return true; //whoo! we made it!
}//compare_atom_coordinates


bool
compare_binary_protein_silent_struct(
	Pose const & lhs,
	Pose const & rhs
) {
	core::io::silent::BinarySilentStruct lhs_silent_struct(lhs, "" );
	std::stringstream lhs_str;
	lhs_silent_struct.print_conformation(lhs_str);

	core::io::silent::BinarySilentStruct rhs_silent_struct(rhs, "" );
	std::stringstream rhs_str;
	rhs_silent_struct.print_conformation(rhs_str);

	return lhs_str.str() == rhs_str.str();
}

////////////////////////////////////////////////////////////////////////////////////////////////

id::NamedAtomID
atom_id_to_named_atom_id(
	core::id::AtomID const & atom_id,
	Pose const & pose
) {
	conformation::Residue const& rsd( pose.residue( atom_id.rsd() ) );
	return core::id::NamedAtomID( rsd.atom_name( atom_id.atomno() ), atom_id.rsd() );
}

///@details returns an AtomID corresponding to your NamedAtomID
/// check for a valid AtomID after this.
/// following conditions return invalid ID :
/// rsd > total_residue
/// atom not present in residue ( e.g., no CB in GLY )
id::AtomID
named_atom_id_to_atom_id(
	core::id::NamedAtomID const & named_atom_id,
	Pose const & pose,
	bool raise_exception /*default true*/
) {
	using namespace core::id;
	// work out the stubID
	if ( named_atom_id.valid() ) {
		if ( named_atom_id.rsd() <= pose.total_residue() ) { //if in range, ... otherwise leave atomno_ on 0 --> valid() == false
			chemical::ResidueType const& rt ( pose.residue_type ( named_atom_id.rsd() ) );
			if ( rt.has( named_atom_id.atom() ) ) {
				return AtomID( rt.atom_index( named_atom_id.atom() ), named_atom_id.rsd() );
			} else {
				// tr.Error << "Error: can't find atom " << named_atom_id.atom() << " in residue "
// 					<< rt.name() << ", residue has " << rt.nheavyatoms() << " heavy atoms." << std::endl;
// 				tr.Error << "atom names are: " << std::endl;
				//rt.show_all_atom_names( tr.Error );
				if ( raise_exception ) throw id::EXCN_AtomNotFound( named_atom_id );
				return id::BOGUS_ATOM_ID;
			}
		} else {
// 			tr.Error << "Error: can't find residue " << named_atom_id.rsd()
// 				<< " in pose (pose.total_residue() = ) "
// 				<< pose.total_residue() << std::endl;
			if ( raise_exception ) throw id::EXCN_AtomNotFound( named_atom_id );
 			return id::BOGUS_ATOM_ID;
		}
	} else {
		if ( raise_exception ) throw id::EXCN_AtomNotFound( named_atom_id );
		return id::BOGUS_ATOM_ID;
	}
}


id::NamedStubID
stub_id_to_named_stub_id( id::StubID const& stub_id, core::pose::Pose const& pose ) {
	using namespace core::id;
	if ( stub_id.center().valid() ) {
		return NamedStubID(
			NamedAtomID( atom_id_to_named_atom_id( stub_id.center() , pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 1 ), pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 2 ), pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 3 ), pose ) )
		);
	} else {
		return NamedStubID( // TODO does this make sense? if input stub is bad, what should be done?
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 1 ), pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 2 ), pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 3 ), pose ) )
		);
	}
}

id::StubID
named_stub_id_to_stub_id( id::NamedStubID const& named_stub_id, core::pose::Pose const& pose ) {
	using namespace core::id;
	if ( named_stub_id.center().valid() ) {
		return StubID(
			AtomID( named_atom_id_to_atom_id( named_stub_id.center() , pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 2 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 3 ), pose ) )
		);
	} else {
		return StubID( // TODO does this make sense? if input stub is bad, what should be done?
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 2 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 3 ), pose ) )
		);
	}
}


///////////////////////////////////////////////////////////////
std::string tag_from_pose( core::pose::Pose const & pose ) {
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;

	std::string tag( "empty_tag" );
	if ( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		tag =
			static_cast< basic::datacache::CacheableString const & >
			( pose.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	}

	return tag;
}

void tag_into_pose( core::pose::Pose & pose, std::string const & tag ) {
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;
	pose.data().set( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,  new CacheableString(tag) );
}

core::Real energy_from_pose(
	core::pose::Pose const & pose, core::scoring::ScoreType const & sc_type
) {
	return pose.energies().total_energies()[ sc_type ];
}

core::Real energy_from_pose(
	core::pose::Pose const & pose, std::string const & sc_type
) {
	return pose.energies().total_energies()[ core::scoring::score_type_from_name( sc_type ) ];
}

core::Real total_energy_from_pose( core::pose::Pose const & pose ) {
	Real total_energy = pose.energies().total_energy();
	// when initiated from a silent struct, total energy lives in "extra data".
	if ( total_energy == Real( 0.0 ) ) getPoseExtraScores( pose, "score", total_energy );
	return total_energy;
}

// criterion for sorting.
bool
sort_pose_by_score( core::pose::PoseOP const & pose1, core::pose::PoseOP const & pose2 ) {
	return ( total_energy_from_pose( *pose1 ) < total_energy_from_pose( *pose2 ) );
}

void
transfer_phi_psi( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose, core::Size ir, core::Size jr )
{
		core::Size tgtlength = tgtpose.total_residue();
		core::Size srclength = srcpose.total_residue();
		for( core::Size ires = ir; ires <= std::min( srclength, std::min( tgtlength, jr)) ; ires++){
			if (!srcpose.residue_type( ires ).is_protein() || !tgtpose.residue_type( ires ).is_protein()) continue;
			tgtpose.set_phi(   ires, srcpose.phi(ires) );
			tgtpose.set_psi(   ires, srcpose.psi(ires) );
			tgtpose.set_omega( ires, srcpose.omega(ires) );
		}
}

void
transfer_phi_psi( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose ){
    transfer_phi_psi( srcpose, tgtpose, 1, std::min( tgtpose.total_residue(), srcpose.total_residue() ) );
}

//fpd transfer the RB geometry from one pose to another
void
transfer_jumps( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose )
{
	core::kinematics::FoldTree f_tgt=tgtpose.fold_tree(), f_src=srcpose.fold_tree();

	// project fold tree
	core::Size srcjmps = srcpose.num_jump();
	for( core::Size ijmp = 1; ijmp <= srcjmps ; ijmp++){
		core::kinematics::Edge srcedge_i = srcpose.fold_tree().jump_edge(ijmp);
		core::Size jjmp = tgtpose.fold_tree().jump_nr( srcedge_i.start(), srcedge_i.stop() );
		if (jjmp != 0)
			tgtpose.set_jump( jjmp, srcpose.jump(ijmp) );
		else
			TR.Debug << "In transfer_jumps() unable to map jump: " << srcedge_i.start() << " ,  " << srcedge_i.stop() << std::endl;
	}
}


/// helper function for residue replacement/residuetype switching
/// these functions should probably move to pose/util.cc
/// @note  Will call new_rsd->fill_missing_atoms if the new residue has atoms
/// that the old one doesn't
void
replace_pose_residue_copying_existing_coordinates(
	pose::Pose & pose,
	Size const seqpos,
	core::chemical::ResidueType const & new_rsd_type
	)
{

	core::conformation::Residue const & old_rsd( pose.residue( seqpos ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( new_rsd_type ) );
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	pose.replace_residue( seqpos, *new_rsd, false );

}

///////////////////////////////////////////////////////////////////////////////
/// @brief removes variant from an existing residue
/// @details
/// @note
core::conformation::ResidueOP
remove_variant_type_from_residue(
	core::conformation::Residue const & old_rsd,
	core::chemical::VariantType const & variant_type,
	pose::Pose const & pose
	)
{
	if ( !old_rsd.has_variant_type( variant_type ) ) return old_rsd.clone();

	// the type of the desired variant residue
	core::chemical::ResidueTypeSet const & rsd_set( old_rsd.residue_type_set() );
	core::chemical::ResidueType const & new_rsd_type( rsd_set.get_residue_type_with_variant_removed( old_rsd.type(), variant_type ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( new_rsd_type, old_rsd, pose.conformation() ) );
	core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	if ( old_rsd.nchi() == new_rsd_type.nchi() ) {
		for ( Size chino=1; chino <= old_rsd.nchi(); ++chino ) {
			new_rsd->set_chi( chino, old_rsd.chi( chino ) );
		}
	} else {
		TR << "The chi angles will not be updated and your dunbrack score for this rotamer will be huge; this function is only meant to add a variant type to a residue of the same type" << std::endl;
	}

	return new_rsd;
}

///////////////////////////////////////////////////////////////////////////////
/// @brief construct a variant of an existing residue
/// @details
/// @note
conformation::ResidueOP
add_variant_type_to_residue(
	conformation::Residue const & old_rsd,
	chemical::VariantType const & variant_type,
	pose::Pose const & pose
	) {

	if ( old_rsd.has_variant_type( variant_type ) ) return old_rsd.clone();

	// the type of the desired variant residue
	chemical::ResidueTypeSet const & rsd_set( old_rsd.residue_type_set() );
	chemical::ResidueType const & new_rsd_type( rsd_set.get_residue_type_with_variant_added( old_rsd.type(), variant_type ) );
	conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( new_rsd_type, old_rsd, pose.conformation() ) );
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	if ( old_rsd.nchi() == new_rsd_type.nchi() ) {
		for ( Size chino=1; chino <= old_rsd.nchi(); ++chino ) {
			new_rsd->set_chi( chino, old_rsd.chi( chino ) );
		}
	} else {
		TR << "The chi angles will not be updated and your dunbrack score for this rotamer will be huge; this function is only meant to add a variant type to a residue of the same type" << std::endl;
	}

	return new_rsd;
}

///////////////////////////////////////////////////////////////////////////////
/// @brief construct a variant of an existing pose residue
/// @details eg make a terminus variant, and replace the orignal in pose.
/// @note this copies any atoms in common between old and new residues, rebuild the others
void
add_variant_type_to_pose_residue(
	pose::Pose & pose,
	chemical::VariantType const & variant_type,
	Size const seqpos
	)
{
	runtime_assert( seqpos != 0 );
	if ( pose.residue( seqpos ).has_variant_type( variant_type ) ) return;

	conformation::Residue const & old_rsd( pose.residue( seqpos ) );

	// the type of the desired variant residue
	chemical::ResidueTypeSet const & rsd_set( old_rsd.residue_type_set() );
	chemical::ResidueType const & new_rsd_type( rsd_set.get_residue_type_with_variant_added( old_rsd.type(), variant_type ) );

	core::pose::replace_pose_residue_copying_existing_coordinates( pose, seqpos, new_rsd_type );

    // update connections
    for (Size i_con=1; i_con<=pose.conformation().residue_type(seqpos).n_residue_connections(); ++i_con) {
        if (pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con) != 0) {
            Size connected_seqpos = pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con);
            Size connected_id = pose.residue(seqpos).connect_map(i_con).connid();
            pose.conformation().update_noncanonical_connection(seqpos, i_con, connected_seqpos, connected_id);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
/// @brief construct a non-variant of an existing pose residue
/// @details eg remove a terminus variant, and replace the orignal in pose.
/// @note this copies any atoms in common between old and new residues, rebuild the others
void
remove_variant_type_from_pose_residue(
	pose::Pose & pose,
	chemical::VariantType const & variant_type,
	Size const seqpos
	)
{
	if ( !pose.residue( seqpos ).has_variant_type( variant_type ) ) return;

	conformation::Residue const & old_rsd( pose.residue( seqpos ) );

	// the type of the desired variant residue
	chemical::ResidueTypeSet const & rsd_set( old_rsd.residue_type_set() );
	chemical::ResidueType const & new_rsd_type( rsd_set.get_residue_type_with_variant_removed( old_rsd.type(), variant_type ) );

	core::pose::replace_pose_residue_copying_existing_coordinates( pose, seqpos, new_rsd_type );

    // update connections
    for (Size i_con=1; i_con<=pose.conformation().residue_type(seqpos).n_residue_connections(); ++i_con) {
        if (pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con) != 0) {
            Size connected_seqpos = pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con);
            Size connected_id = pose.residue(seqpos).connect_map(i_con).connid();
            pose.conformation().update_noncanonical_connection(seqpos, i_con, connected_seqpos, connected_id);
        }
    }
}

void
add_lower_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
	)
{
	add_variant_type_to_pose_residue( pose, chemical::LOWER_TERMINUS, seqpos );
}

void
add_upper_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
	)
{
	add_variant_type_to_pose_residue( pose, chemical::UPPER_TERMINUS, seqpos );
}

void
remove_lower_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
	)
{
	core::pose::remove_variant_type_from_pose_residue( pose, chemical::LOWER_TERMINUS, seqpos );
}

void
remove_upper_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
	)
{
	core::pose::remove_variant_type_from_pose_residue( pose, chemical::UPPER_TERMINUS, seqpos );
}

///@brief returns a Distance
core::Real
pose_max_nbr_radius( Pose const & pose )
{
	core::Real maxrad( 0.0 );
	for ( 	core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( pose.residue( ii ).nbr_radius() > maxrad ) maxrad = pose.residue_type(ii).nbr_radius();
	}
	return maxrad;
}

///////////////////////////////////////////////////////////////////////////////
void
setup_dof_to_torsion_map(
	pose::Pose const & pose,
	id::DOF_ID_Map< id::TorsionID > & dof_map
	)
{
	Size const n_res( pose.n_residue() );

	// Set DOF mask size and initialize to an invalid TorsionID
	core::pose::initialize_dof_id_map( dof_map, pose, id::BOGUS_TORSION_ID );

	for ( Size i = 1; i <= n_res; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );

		// first the backbone torsion angles
		{
			// PHIL note the Residue-based helper functions you need for this
			// PHIL also note the pose.conformation() interface
			int const n_torsions( rsd.mainchain_atoms().size() );
			for ( int j=1; j<= n_torsions; ++j ) {
				id::TorsionID const tor_id( i, id::BB, j );
				id::DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( tor_id ) );
				if ( id.valid() ) {
					dof_map[ id ] = tor_id;
				}
			}
		}

		{
			// PHIL note the Residue-based helper functions you need for this
			// PHIL also note the pose.conformation() interface
			int const n_torsions( rsd.nchi() );
			for ( int j=1; j<= n_torsions; ++j ) {
				id::TorsionID const tor_id( i, id::CHI, j );
				id::DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( tor_id ) );
				if ( id.valid() ) {
					dof_map[ id ] = tor_id;
				}
			} // j=1,chi-torsions
		} // scope


	} // i=1,n_res

	for ( Size i = 1; i <= pose.num_jump(); ++i ) {
		for ( int j=1; j<= 6; ++j ) {
			id::TorsionID const tor_id(i,id::JUMP,j);
			id::DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( tor_id ) );
			dof_map[ id ] =tor_id;
		}
	} // i=1.num_jump
}


///////////////////////////////////////////////////////////////////////////////
// convert from allow-bb/allow-chi MoveMap to simple DOF_ID boolean mask needed by the minimizer
void
setup_dof_mask_from_move_map(
		kinematics::MoveMap const & mm,
		pose::Pose const & pose,
		id::DOF_ID_Mask & dof_mask
)
{
	using namespace id;

	bool const PHI_default( mm.get( PHI ) );

	// Set DOF mask size and initialize to false.
	initialize_dof_id_map( dof_mask, pose, false );

	// DOF_Type defaults
	// could/should do this with a loop over types?
	// currently ignoring rb types, set these individually by jump number
	dof_mask.set( PHI  , mm.get( PHI   ) );
	dof_mask.set( THETA, mm.get( THETA ) );
	dof_mask.set( D    , mm.get( D     ) );

	// Torsion angles
	Size const n_res( pose.n_residue() );
	for ( Size i = 1; i <= n_res; ++i ) {
		// PHIL note the Residue-based helper functions you need for this
		// PHIL also note the pose.conformation() interface

		conformation::Residue const & rsd( pose.residue(i));

		// first the backbone torsion angles
		Size const n_bb_torsions( rsd.mainchain_atoms().size() );

		// Note: In many (most?) cases, a ResidueType with a ring will have overlap between its internal ring torsions
		// and its main-chain torsions as defined by the AtomTree.  In such cases, the ring atoms may or may not be
		// considered part of the backbone.  For example, in the case of carbohydrates, phi, psi, and omega are to be
		// considered backbone torsions, but the other main-chain torsions should be ignored.  If one wants to modify
		// those other mainchain torsions, they should be treated as nu angles.  This block of code is specific to
		// carbohydrates, but if any other ResidueTypes have strict definitions of backbone vs ring angles, similar code
		// could be added here to "subtract out" main chain torsions already covered by nu torsions. ~Labonte
		Size n_cyclic_main_chain_torsions = 0;  // By default, a residue is acyclic.
		if (rsd.is_carbohydrate()) {
			if (rsd.carbohydrate_info()->is_acyclic()) {
				n_cyclic_main_chain_torsions = 0;
			} else if (rsd.carbohydrate_info()->has_exocyclic_linkage()) {
				n_cyclic_main_chain_torsions = n_bb_torsions - 3;  // minus PHI, PSI, & OMEGA, the actual BB torsions
			} else /* doesn't have an omega angle */ {
				n_cyclic_main_chain_torsions = n_bb_torsions - 2;  // minus PHI & PSI, the actual BB torsions
			}
		}

		for ( uint j=1; j<= n_bb_torsions; ++j ) {
			bool mm_setting;
			if (j <= n_cyclic_main_chain_torsions) {
				mm_setting = false;  // Do not move "backbone" torsions that are part of a ring.
			} else {
				mm_setting = mm.get(TorsionID(i, BB, j));
			}
			if ( mm_setting == PHI_default ) continue;
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id(TorsionID(i, BB, j)));
			if ( id.valid() ) {  // If not valid, it's probably just a terminal/chainbreak torsion.
				dof_mask[ id ] = mm_setting;
			}
		} // j=1, n_bb_torsions

		// then the side chain torsions
		Size const n_chi_torsions( rsd.nchi() );
		for ( uint j = 1; j <= n_chi_torsions; ++j ) {
			bool mm_setting = mm.get(TorsionID(i, CHI, j));
			if ( mm_setting == PHI_default ) continue;
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id(TorsionID(i, CHI, j)));
			if ( id.valid() ) {
				dof_mask[ id ] = mm_setting;
			} else {
				TR.Warning << "WARNING: Unable to find atom_tree atom for this " <<
						"Rosetta chi angle: residue " << i << " CHI " << j << std::endl;
			}
		} // j=1, n_chi_torsions

		// finally, the internal ring torsions
		Size const n_nu_torsions = rsd.n_nus();
		for ( uint j = 1; j <= n_nu_torsions; ++j ) {
			bool mm_setting = mm.get(TorsionID(i, NU, j));
			if ( mm_setting == PHI_default ) continue;
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id(TorsionID(i, NU, j)));
			if ( id.valid() ) {
				dof_mask[ id ] = mm_setting;
			} else {
				TR.Warning << "WARNING: Unable to find atom_tree atom for this " <<
						"Rosetta nu angle: residue " << i << " NU " << j << std::endl;
			}
		} // j=1, n_nu_torsions
	} // i=1, n_res

	// Jumps.
	for ( Size i=1; i<= pose.num_jump(); ++i ) {
		if ( mm.get_jump(i) ) {
			for ( int j=1; j<= 6; ++j ) {
				DOF_ID const & id( pose.conformation().dof_id_from_torsion_id(TorsionID(i, JUMP, j)));
				dof_mask[ id ] = true;
			}
		}
	} // i=1, num_jump

	/////////////////////////
	// DOFs set individually
	for ( kinematics::MoveMap::DOF_ID_Map::const_iterator it=mm.dof_id_begin(), it_end=mm.dof_id_end();
			it != it_end; ++it ) {
		dof_mask[ it->first ] = it->second;
	}
} // setup_dof_mask_from_move_map


bool
has_chain(std::string const & chain, core::pose::Pose const & pose){
	assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return has_chain(chain_char, pose);
}

bool
has_chain(char const & chain, core::pose::Pose const & pose){
	for(core::Size i=1; i <= pose.conformation().num_chains(); ++i){
		char this_char= get_chain_from_chain_id(i, pose);
		if(this_char == chain){
			return true;
		}
	}
	return false;
}

bool
has_chain(core::Size chain_id, core::pose::Pose const & pose){
	for(core::Size i=1; i <= pose.conformation().num_chains(); ++i){
		if (i == chain_id){
			return true;
		}
	}
	return false;
}

std::set<core::Size>
get_jump_ids_from_chain_ids(std::set<core::Size> const chain_ids, core::pose::Pose const & pose){
	std::set<core::Size> jump_ids;
	std::set<core::Size>::const_iterator chain_id= chain_ids.begin();
	for(; chain_id != chain_ids.end(); ++chain_id){
		core::Size jump_id= get_jump_id_from_chain_id(*chain_id, pose);
		jump_ids.insert(jump_id);
	}
	return jump_ids;
}

core::Size
get_jump_id_from_chain_id(core::Size const & chain_id,const core::pose::Pose & pose){
	for(core::Size jump_id=1; jump_id <= pose.num_jump(); jump_id++){
		core::Size ligand_residue_id= (core::Size) pose.fold_tree().downstream_jump_residue(jump_id);
		core::Size ligand_chain_id= pose.chain(ligand_residue_id);
		if(chain_id==ligand_chain_id){
			return jump_id;
		}
	}
	utility_exit();
	return 0;// this will never happen
}

core::Size
get_chain_id_from_chain(std::string const & chain, core::pose::Pose const & pose){
	assert(chain.size()==1);// chain is one char
	if( chain.size() > 1) utility_exit_with_message("Multiple chain_ids per chain! Are you using '-treat_residues_in_these_chains_as_separate_chemical_entities', and not using compatible movers?" );
	char chain_char= chain[0];
	return get_chain_id_from_chain(chain_char, pose);
}

core::Size
get_chain_id_from_chain(char const & chain, core::pose::Pose const & pose){
	utility::vector1<core::Size> chain_ids = get_chain_ids_from_chain(chain, pose);
	if(chain_ids.size() != 1)
	{
		throw utility::excn::EXCN_RangeError("chain_id "+utility::to_string(chain)+" does not exist");
	}
	return chain_ids[1];
}

utility::vector1<core::Size>
get_chain_ids_from_chain(std::string const & chain, core::pose::Pose const & pose){
	assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return get_chain_ids_from_chain(chain_char, pose);

}

utility::vector1<core::Size>
get_chain_ids_from_chain(char const & chain, core::pose::Pose const & pose){
	utility::vector1<core::Size> chain_ids;
	for(core::Size i=1; i <= pose.conformation().num_chains(); i++){
		char this_char= get_chain_from_chain_id(i, pose);
		if(this_char == chain){
			chain_ids.push_back(i);
		}
	}
	return chain_ids;
}

core::Size
get_jump_id_from_chain(std::string const & chain, core::pose::Pose const & pose){
	core::Size chain_id= get_chain_id_from_chain(chain, pose);
	return get_jump_id_from_chain_id(chain_id, pose);
}

core::Size
get_jump_id_from_chain(char const & chain, core::pose::Pose const & pose){
	core::Size chain_id= get_chain_id_from_chain(chain, pose);
	return get_jump_id_from_chain_id(chain_id, pose);
}

utility::vector1<core::Size>
get_jump_ids_from_chain(char const & chain, core::pose::Pose const & pose){
	utility::vector1<core::Size> jump_ids;
	utility::vector1<core::Size> chain_ids = get_chain_ids_from_chain(chain, pose);
	BOOST_FOREACH(core::Size chain_id, chain_ids){
		jump_ids.push_back( get_jump_id_from_chain_id(chain_id, pose));
	}
	return jump_ids;
}

utility::vector1<core::Size>
get_jump_ids_from_chain(std::string const & chain, core::pose::Pose const & pose){
	assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return get_jump_ids_from_chain(chain_char, pose);
}

core::Size get_chain_id_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose){
	core::Size ligand_residue_id= (core::Size) pose.fold_tree().downstream_jump_residue(jump_id);
	return pose.chain(ligand_residue_id);
}

char
get_chain_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose){
	core::Size chain_id= get_chain_id_from_jump_id(jump_id, pose);
	return get_chain_from_chain_id(chain_id, pose);
}

core::conformation::ResidueCOPs
get_chain_residues(core::pose::Pose const & pose, core::Size const chain_id){
	core::Size begin= pose.conformation().chain_begin(chain_id);
	core::Size const end= pose.conformation().chain_end(chain_id);
	core::conformation::ResidueCOPs residues;
	for(; begin <= end; ++begin){
		residues.push_back( new core::conformation::Residue(pose.residue(begin)) );
	}
	return residues;
}

char
get_chain_from_chain_id(core::Size const & chain_id, core::pose::Pose const & pose){
	core::Size first_chain_residue= pose.conformation().chain_begin( chain_id );
	return pose.pdb_info()->chain(first_chain_residue);
}

core::Size num_heavy_atoms(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
) {
	core::Size total_heavy_atoms = 0;
	for (; begin <= end; ++begin) {
		total_heavy_atoms += pose.residue(begin).nheavyatoms();
	}
	TR.Debug << "# of heavy atoms: "<< total_heavy_atoms << std::endl;
	return total_heavy_atoms;
}

core::Size num_atoms(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
) {
	core::Size total_atoms = 0;
	for (; begin <= end; ++begin) {
		total_atoms += pose.residue(begin).natoms();
	}
	TR.Debug << "# of heavy atoms: "<< total_atoms << std::endl;
	return total_atoms;
}

core::Size num_hbond_acceptors(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
) {
	core::Size total_hbond_acceptors = 0;
	for (; begin <= end; ++begin) {
		total_hbond_acceptors += pose.residue(begin).n_hbond_acceptors();
	}
	TR.Debug << "# of heavy atoms: "<< total_hbond_acceptors << std::endl;
	return total_hbond_acceptors;
}

core::Size num_hbond_donors(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
) {
	core::Size total_hbond_donors = 0;
	for (; begin <= end; ++begin) {
		total_hbond_donors += pose.residue(begin).n_hbond_donors();
	}
	TR.Debug << "# of heavy atoms: "<< total_hbond_donors << std::endl;
	return total_hbond_donors;
}

core::Size num_chi_angles(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
) {
	core::Size total_chi_angles = 0;
	for (; begin <= end; ++begin) {
		total_chi_angles += pose.residue(begin).nchi();
	}
	return total_chi_angles;
}

core::Real
mass(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
){
	core::Real mass = 0;
	for (; begin <= end; ++begin) {
		mass += pose.residue(begin).type().mass();
	}
	return mass;
}


core::Size get_hash_from_chain(char const & chain, core::pose::Pose const & pose)
{
	std::size_t hash = 0;

	core::Size chain_id = get_chain_id_from_chain(chain,pose);
	core::Size chain_begin = pose.conformation().chain_begin(chain_id);
	core::Size chain_end = pose.conformation().chain_end(chain_id);
	for(core::Size res_num = chain_begin; res_num <= chain_end; ++res_num)
	{
		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for(core::Size atom_num = 1; atom_num <= natoms; ++atom_num)
		{
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			boost::hash_combine(hash,current_xyz);
		}
	}

	return hash;
}

core::Size get_hash_excluding_chain(char const & chain, core::pose::Pose const & pose)
{
	std::size_t hash = 0;

	core::Size chain_id = get_chain_id_from_chain(chain,pose);

	for(core::Size res_num = 1; res_num <= pose.n_residue(); ++res_num)
	{
		if((int)chain_id == pose.chain(res_num))
		{
			continue;
		}
		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for(core::Size atom_num = 1; atom_num <= natoms; ++atom_num)
		{
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			boost::hash_combine(hash,current_xyz);
		}
	}

	return hash;
}

std::string get_sha1_hash_excluding_chain(char const & chain, core::pose::Pose const & pose)
{

	std::stringstream coord_stream;

	core::Size chain_id = get_chain_id_from_chain(chain,pose);

	for(core::Size res_num = 1; res_num <= pose.n_residue(); ++res_num)
	{
		if((int)chain_id == pose.chain(res_num))
		{
			continue;
		}

		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for(core::Size atom_num = 1; atom_num <= natoms; ++atom_num)
		{
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			coord_stream << numeric::truncate_and_serialize_xyz_vector(current_xyz,3);
		}
	}
	return utility::string_to_sha1(coord_stream.str());
}

void
initialize_disulfide_bonds(
	Pose & pose
) {

		// disulfides
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	// Fix disulfides if a file is given
	if ( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].user() ) {
		core::io::raw_data::DisulfideFile ds_file( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ]() );
		utility::vector1< std::pair<Size,Size> > disulfides;
		ds_file.disulfides(disulfides, pose);
		pose.conformation().fix_disulfides( disulfides );
	} else {
		if ( option[ in::detect_disulf ].user() ?
				option[ in::detect_disulf ]() : // detect_disulf true
				pose.is_fullatom() // detect_disulf default but fa pose
			) {
			pose.conformation().detect_disulfides();
		}
	}
}

std::string extract_tag_from_pose( core::pose::Pose &pose )
{
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;
	using basic::datacache::CacheableStringOP;

	if( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ){
			CacheableStringOP data =  dynamic_cast< CacheableString* > (  (pose.data().get_raw_ptr( ( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG  )) ));
			if( data.get() == NULL ) return std::string("UnknownTag");
			else               return data->str();
	}

	return std::string("UnknownTag");
}

core::id::SequenceMapping sequence_map_from_pdbinfo( Pose const & first, Pose const & second ) {
	core::id::SequenceMapping retval(first.total_residue(), second.total_residue());
	core::pose::PDBInfoCOP first_pdbinfo = first.pdb_info();
	core::pose::PDBInfoCOP second_pdbinfo = second.pdb_info();

	if ( first_pdbinfo && !first_pdbinfo->obsolete() && second_pdbinfo && !second_pdbinfo->obsolete() ) {
		for ( core::Size ii(1); ii<= first.total_residue(); ++ii ) {
			// pdb2pose returns 0 for "not found" - 0 is also used for "not found" for SequenceMapping.
			retval[ii] = second_pdbinfo->pdb2pose( first_pdbinfo->chain(ii), first_pdbinfo->number(ii), first_pdbinfo->icode(ii) );
		}
	} else {
		TR << "One or both poses do not have usable PDBInfo, using sequence alignment instead." << std::endl;
		retval = core::sequence::map_seq1_seq2( new core::sequence::Sequence(first), new core::sequence::Sequence(second) );
	}

	return retval;
}

core::Size canonical_residue_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for(core::Size i = 1; i <= pose.total_residue();++i)
	{
		core::conformation::Residue const & resi(pose.residue(i));
		if(resi.aa() <= core::chemical::num_canonical_aas)
		{
			++count;
		}
	}
	return count;
}

core::Size noncanonical_residue_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for(core::Size i = 1; i <= pose.total_residue();++i)
	{
		core::conformation::Residue const & resi(pose.residue(i));
		if(resi.aa() > core::chemical::num_canonical_aas)
		{
			++count;
		}
	}
	return count;
}

core::Size canonical_atom_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for(core::Size i = 1; i <= pose.total_residue();++i)
	{
		core::conformation::Residue const & resi(pose.residue(i));
		if(resi.aa() <= core::chemical::num_canonical_aas)
		{
			count += resi.natoms();
		}
	}
	return count;
}

core::Size noncanonical_atom_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for(core::Size i = 1; i <= pose.total_residue();++i)
	{
		core::conformation::Residue const & resi(pose.residue(i));
		if(resi.aa() > core::chemical::num_canonical_aas)
		{
			count += resi.natoms();
		}
	}
	return count;
}

core::Size noncanonical_chi_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for(core::Size i = 1; i <= pose.total_residue();++i)
	{
		core::conformation::Residue const & resi(pose.residue(i));
		if(resi.aa() > core::chemical::num_canonical_aas)
		{
			count += resi.nchi();
		}
	}
	return count;
}

////////////////////////////////////////////////////////////////////////////////////
/// @begin center_of_mass
///
/// @brief calculates the center of mass of a pose
/// @detailed
///				the start and stop positions (or residues) within the pose are used to
///				find the starting and finishing locations
///
/// @authors Monica Berrondo June 14 2007
///
/// @last_modified Javier Castellanos June 4 2012
/////////////////////////////////////////////////////////////////////////////////
numeric::xyzVector< core::Real>
center_of_mass(
	pose::Pose const & pose,
	int const start,
	int const stop
)
{
	Vector center( 0.0 );
	for ( int i=start; i<=stop; ++i ) {
		if( !pose.residue( i ).is_protein()) {
			Vector ca_pos( pose.residue( i ).nbr_atom_xyz() );
			center += ca_pos;
	 	} else {
			Vector ca_pos( pose.residue( i ).atom( "CA" ).xyz() );
			center += ca_pos;
			}
	}
	center /= (stop-start+1);

	return center;
}

////////////////////////////////////////////////////////////////////////////////////
/// @begin residue_center_of_mass
///
/// @brief calculates the center of mass of a pose
/// @detailed
///				the start and stop positions (or residues) within the pose are used to
///				find the starting and finishing locations
///
/// @authors Monica Berrondo June 14 2007
///
/// @last_modified Javier Castellanos June 4 2012
/////////////////////////////////////////////////////////////////////////////////
int
residue_center_of_mass(
	pose::Pose const & pose,
	int const start,
	int const stop
)
{
	Vector center = center_of_mass(pose, start, stop );
	return core::pose::return_nearest_residue( pose, start, stop, center );
}

////////////////////////////////////////////////////////////////////////////////////
/// @begin core::pose::return_nearest_residue
///
/// @brief finds the residue nearest some position passed in (normally a
///		center of mass)
/// @detailed
///				the start and stop positions (or residues) within the pose are used to
///				find the starting and finishing locations
///
/// @authors Monica Berrondo June 14 2007
///
/// @last_modified June 29 2007
/////////////////////////////////////////////////////////////////////////////////
int
return_nearest_residue(
	pose::Pose const & pose,
	int const begin,
	int const end,
	Vector center
)
{
	Real min_dist = 9999.9;
	int res = 0;
	for ( int i=begin; i<=end; ++i )
	{
		Vector ca_pos;
		if( !pose.residue( i ).is_protein() ){
			ca_pos = pose.residue( i ).nbr_atom_xyz();
			} else {
		//Vector ca_pos( pose.residue( i ).atom( "CA" ).xyz() );
			ca_pos = pose.residue( i ).atom( "CA" ).xyz() ;
			}

		ca_pos -= center;
		Real tmp_dist( ca_pos.length_squared() );
		if ( tmp_dist < min_dist ) {
			res = i;
			min_dist = tmp_dist;
		}
	}
	return res;
}

#ifdef USELUA
void lregister_util( lua_State * lstate ) {
	luabind::module(lstate, "core")
	[
		luabind::namespace_("pose")
		[
			luabind::def("getExtraScore", &getPoseExtraScores, luabind::pure_out_value(_3)),
			luabind::def("setExtraScore", &setPoseExtraScores),
			luabind::def("get_comment", &get_comment, luabind::pure_out_value(_3)),
			luabind::def("add_comment", &add_comment),
			luabind::def("getScore", (core::Real (*) (core::pose::Pose &, std::string const & ) ) &energy_from_pose),
			luabind::def("getTotalScore", (core::Real (*) (core::pose::Pose & ) ) &total_energy_from_pose)
		]
	];
}
#endif

// silly conversion from std::map< AtomID, AtomID> to rosetta's silly AtomID_Map class.
id::AtomID_Map< id::AtomID >
convert_from_std_map( std::map< id::AtomID, id::AtomID > const & atom_map,
											core::pose::Pose const & pose ){
	id::AtomID_Map< id::AtomID > atom_ID_map;
	initialize_atomid_map( atom_ID_map, pose, id::BOGUS_ATOM_ID );
	for ( std::map< id::AtomID, id::AtomID >::const_iterator it = atom_map.begin(); it != atom_map.end(); it++ ){
		atom_ID_map.set( it->first, it->second );
	}
	return atom_ID_map;
}

/// @brief Add cutpoint variants to all residues annotated as cutpoints in the pose.
void
correctly_add_cutpoint_variants( core::pose::Pose & pose ) {
	const core::kinematics::FoldTree& tree(pose.fold_tree());
	for (core::Size i = 1; i < pose.total_residue(); ++i) { // Less than because cutpoints are between i and i+1
		if ( tree.is_cutpoint(i) ) {
			correctly_add_cutpoint_variants( pose, i, false );
		}
	}
}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// try to unify all cutpoint addition into this function.
	void
	correctly_add_cutpoint_variants( core::pose::Pose & pose,
																	 Size const cutpoint_res,
																	 bool const check_fold_tree /* = true*/){

		using namespace core::chemical;

		runtime_assert( cutpoint_res < pose.total_residue() );
		if ( check_fold_tree ) runtime_assert( pose.fold_tree().is_cutpoint( cutpoint_res ) );

		remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS, cutpoint_res );
		remove_variant_type_from_pose_residue( pose, "THREE_PRIME_PHOSPHATE", cutpoint_res );
		remove_variant_type_from_pose_residue( pose, "C_METHYLAMIDATION", cutpoint_res );

		remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS, cutpoint_res + 1 );
		remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, cutpoint_res + 1 );
		remove_variant_type_from_pose_residue( pose, "FIVE_PRIME_PHOSPHATE", cutpoint_res + 1 );
		remove_variant_type_from_pose_residue( pose, "N_ACETYLATION", cutpoint_res + 1);

		if ( pose.residue_type( cutpoint_res ).is_RNA() )	 rna::correctly_position_cutpoint_phosphate_torsions( pose, cutpoint_res );

		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint_res   );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint_res + 1 );

		// important -- to prevent artificial penalty from steric clash.
		if ( pose.residue_type( cutpoint_res ).is_NA() ) {
			runtime_assert( pose.residue_type( cutpoint_res + 1 ).is_NA() );
			pose.conformation().declare_chemical_bond( cutpoint_res, " O3'", cutpoint_res+1, " P  " );
		}	else if ( pose.residue_type( cutpoint_res ).is_protein() ) {
			runtime_assert( pose.residue_type( cutpoint_res + 1 ).is_protein() );
			pose.conformation().declare_chemical_bond( cutpoint_res, " C  ", cutpoint_res+1, " N  " );
		}
	}


} // pose
} // core
