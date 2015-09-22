// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author ashworth
/// @brief

#include <protocols/dna/PDBOutput.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/dna/util.hh> // dna_full_name3
#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh> // make_hbBasetoAcc_unitvector
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

#include <utility/file/file_sys_util.hh> // file_exists, create_directory
#include <utility/file/gzip_util.hh>
#include <utility/string_util.hh>

#include <algorithm> // std::min
#include <iostream>
#include <fstream>

// option key includes

#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>
#include <utility/io/ozstream.hh>
#include <ObjexxFCL/format.hh>


namespace protocols {
namespace dna {

using utility::string_split;
using utility::vector1;
using namespace core;
using namespace core::pose;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::scoring;
using namespace core::pack;
using namespace basic::options;
using namespace ObjexxFCL::format;

using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static THREAD_LOCAL basic::Tracer TR( "protocols.dna.PDBOutput", t_info );
////////////////////////////////////////////////////////////////////////////////////////////////////

PDBOutput::PDBOutput()
: jd2::PDBJobOutputter(),
	pose_copy_(/* 0 */),
	reference_pose_(/* 0 */),
	chi_diff_threshold_(0.0001),
	mainchain_torsion_diff_threshold_(0.0001),
	enabled_(true)
{}

PDBOutput::~PDBOutput(){}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details pose is const here, so it must be scored already if score information is expected in output file
/// @author ashworth
void
PDBOutput::final_pose( JobOP job, Pose const & pose, std::string const & /*tag*/ )
{
	if ( !enabled_ ) return; // to allow easy overrides of excess pdb writing in higher-level code
	call_output_observers( pose, job );
	pose_copy_ = PoseOP( new Pose( pose ) );
	ozstream pdbout( extended_name(job) );
	if ( !pdbout.good() ) {
		utility_exit_with_message( "Unable to open file: " + extended_name(job) + "\n" );
	}
	pdbout << "REMARK BEGIN ROSETTA INFO\n";
	extract_data_from_Job( job, pdbout );
	output_pdb( pdbout );
	// append job info to pdb (parent virtual function call)
	pdbout.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details scores pdb
/// @author ashworth
void
PDBOutput::operator() (
	Pose const & pose,
	std::string const & name
)
{
	if ( !enabled_ ) return; // to allow easy overrides of excess pdb writing in higher-level code
	pose_copy_ = PoseOP( new Pose( pose ) );
	make_subdirs( name );
	ozstream pdbout( name );
	if ( !pdbout.good() ) {
		utility_exit_with_message( "Unable to open file: " + name + "\n" );
	}
	pdbout << "REMARK BEGIN ROSETTA INFO\n";
	output_info( pdbout );
	output_pdb( pdbout );
	pdbout.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details
/// @author ashworth
void
PDBOutput::output_pdb( ozstream & pdbout )
{
	if ( ! enabled_ ) return; // to allow easy overrides of excess pdb writing in higher-level code
	runtime_assert( pose_copy_ != 0 );
	if ( ! score_function_ ) score_function_ = get_score_function();
	( *score_function_ )( *pose_copy_ );
	// if the sparse_pdb_output option is used, this limits most output to residues which differ from the reference structure (if it exists)
	get_residue_indices_to_output();
	output_score_info( pdbout );
	output_design_tags( pdbout );
	output_hbond_info( pdbout );
	output_buried_unsatisfied_hbonds( pdbout );
	pose_copy_->dump_pdb( pdbout, res_indices_to_output_ );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void PDBOutput::starting_pose( Pose const & pose ) { reference_pose( pose ); }

void PDBOutput::reference_pose( Pose const & pose )
{
	// make hard copy to guarantee that the reference pose remains unchanged
	reference_pose_ = PoseCOP( PoseOP( new Pose( pose ) ) );
	designed_residues_.assign( reference_pose_->total_residue(), false );
}
pose::PoseCOP PDBOutput::reference_pose() const { return reference_pose_; }

void PDBOutput::score_function( ScoreFunction const & sf )
{
	score_function_ = sf.clone();
}
ScoreFunctionCOP PDBOutput::score_function() const { return score_function_; }

void PDBOutput::designed_residue(
	Size index,
	bool value /* = true */
)
{
	Size const s( designed_residues_.size() );
	int const ext( index - s );
	// not sure this auto-expansion is working...
	if ( ext > 0 ) designed_residues_.insert( designed_residues_.end(), ext, false );
	designed_residues_[ index ] = value;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// update designed_residues_ from a PackerTask
void
PDBOutput::note_designed_residues( PackerTaskCOP ptask )
{
	for ( Size index(1), end( ptask->total_residue() ); index < end; ++index ) {
		if ( ptask->design_residue(index) ||
				ptask->residue_task(index).has_behavior("TARGET") ) {
			designed_residue(index);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
PDBOutput::output_info( ozstream & pdbout )
{
	for ( StringsMap::const_iterator itr( info_map_.begin() ), end( info_map_.end() );
			itr != end; ++itr ) {
		// the 'key' or title for this particular set of information
		pdbout << itr->first << '\n';
		for ( Strings::const_iterator line( itr->second.begin() ), end2( itr->second.end() );
				line != end2; ++line ) {
			pdbout << *line << '\n';
		}
		pdbout << "REMARK\n";
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
PDBOutput::add_info(
	std::string const & key,
	Strings const & info,
	bool append /* = true */
)
{
	StringsMap::iterator finditer( info_map_.find( key ) );
	if ( finditer == info_map_.end() || !append ) {
		info_map_[ key ] = info;
	} else {
		Strings & existing_info( finditer->second );
		for ( Strings::const_iterator it( info.begin() ), end( info.end() ); it != end; ++it ) {
			existing_info.push_back( *it );
		}
	}
}

bool
PDBOutput::remove_info( std::string const & key )
{
	return info_map_.erase( key );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private methods

/// @brief compares identity, then internal degrees of freedom for between residues
bool
PDBOutput::residues_are_different(
	Residue const & res1,
	Residue const & res2
) const
{
	// check aa
	if ( res1.aa() != res2.aa() ) return true;
	// check number of heavy atoms
	if ( res1.nheavyatoms() != res2.nheavyatoms() ) return true;
	// check sidechain torsions
	vector1< Real > const & res1chi( res1.chi() ), res2chi( res2.chi() );
	for ( vector1< Real >::const_iterator it1( res1chi.begin() ), it2( res2chi.begin() ),
			end1( res1chi.end() ), end2( res2chi.end() ); it1 != end1 && it2 != end2; ++it1, ++it2 ) {
		if ( std::abs( *it1 - *it2 ) > chi_diff_threshold_ ) {
			//   TR(t_info) << "Chis differ: " << F(6,3,*it1) << " " << F(6,3,*it2) << std::endl;
			return true;
		}
	}
	// check mainchain torsions
	vector1< Real > const & res1mc( res1.mainchain_torsions() ), res2mc( res2.mainchain_torsions() );
	for ( vector1< Real >::const_iterator it1( res1mc.begin() ), it2( res2mc.begin() ),
			end1( res1mc.end() ), end2( res2mc.end() ); it1 != end1 && it2 != end2; ++it1, ++it2 ) {
		if ( std::abs( *it1 - *it2 ) > mainchain_torsion_diff_threshold_ ) {
			//   TR(t_info) << "mainchain dihedrals differ: " << F(6,3,*it1) << " " << F(6,3,*it2) << std::endl;
			return true;
		}
	}
	// 'not different'
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
PDBOutput::get_residue_indices_to_output()
{
	res_indices_to_output_.clear();
	if ( !reference_pose_ || !option[ OptionKeys::dna::design::sparse_pdb_output ]() ) {
		for ( Size i(1), end( pose_copy_->total_residue() ); i <= end; ++i ) {
			res_indices_to_output_.push_back(i);
		}
	} else {
		for ( Size i(1), endpose( pose_copy_->total_residue() ),
				endref( reference_pose_->total_residue() ); i <= endpose; ++i ) {
			if ( i <= endref &&
					!residues_are_different( pose_copy_->residue(i), reference_pose_->residue(i) ) ) continue;
			res_indices_to_output_.push_back(i);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
std::string
string_join(
	PDBOutput::Strings const & list,
	std::string sep = ", "
)
{
	std::string os;
	for ( std::list< std::string >::const_iterator it( list.begin() ), end( list.end() );
			it != end; ++it ) os += sep + *it;
	return os;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details outputs a set of tags describing how residues varied
/// @author ashworth
void
PDBOutput::output_design_tags( ozstream & pdbout ) const
{
	pdbout << "REMARK Residues varied in this design:\n";

	std::list< std::string > moved, moved_DNA, mutated, mutated_DNA;

	for ( Size index(1), nres( pose_copy_->total_residue() ); index <= nres; ++index ) {
		std::ostringstream os;
		if ( pose_copy_->pdb_info() ) {
			os << pose_copy_->pdb_info()->number( index ) << " "
				<< pose_copy_->pdb_info()->chain( index );
		} else os << index << " " << pose_copy_->chain( index );

		if ( reference_pose_ ) {
			// compare to reference pose if it exists
			// mutation
			if ( pose_copy_->residue_type(index).aa() != reference_pose_->residue_type(index).aa() ) {
				if ( pose_copy_->residue_type(index).is_DNA() ) mutated_DNA.push_back( os.str() );
				else mutated.push_back( os.str() );
				// no mutation, but structure changed
			} else if ( residues_are_different( pose_copy_->residue(index),
					reference_pose_->residue(index) ) ) {
				if ( pose_copy_->residue_type(index).is_DNA() ) moved_DNA.push_back( os.str() );
				else moved.push_back( os.str() );
			}
		}
	}
	pdbout << "REMARK MovedRes" << string_join( moved ) << '\n';
	pdbout << "REMARK MovedDNA" << string_join( moved_DNA ) << '\n';
	pdbout << "REMARK MutatedRes" << string_join( mutated ) << '\n';
	pdbout << "REMARK MutatedDNA" << string_join( mutated_DNA ) << '\n';
	pdbout << "REMARK\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
PDBOutput::output_score_info( ozstream & pdbout )
{
	// weights (non-zero)
	EnergyMap const & weights( score_function_->weights() );
	pdbout << "REMARK Non-zero ScoreFunction weights:\n";
	for ( int i(1); i <= n_score_types; ++i ) {
		if ( weights[ ScoreType(i) ] == 0 ) continue;
		std::string const name( ScoreTypeManager::name_from_score_type( ScoreType(i) ) );
		pdbout << "REMARK " << A(16,name) << " " << F( 6, 4, weights[ ScoreType(i) ] ) << '\n';
	}
	pdbout << "REMARK\n";
	// end weights

	// per-residue energies
	if ( option[ OptionKeys::score::output_residue_energies ]() ) {
		pdbout << "REMARK Weighted per-residue energies\n";
		// header
		int const colwidth(10);
		pdbout << "REMARK " << A( 4, "pdb" ) << A( 4, "type" ) << A( 2, "ch" );
		for ( int i(1); i <= n_score_types; ++i ) {
			if ( weights[ ScoreType(i) ] == 0 ) continue;
			pdbout << A( colwidth, ScoreTypeManager::name_from_score_type( ScoreType(i) ) );
		}
		pdbout << '\n';
		// end header
		// residues
		for ( vector1< Size >::const_iterator pos( res_indices_to_output_.begin() ),
				end( res_indices_to_output_.end() ); pos != end; ++pos ) {
			pdbout << "REMARK ";
			if ( pose_copy_->pdb_info() ) {
				pdbout << I( 4, pose_copy_->pdb_info()->number(*pos) )
					<< A( 4, pose_copy_->residue(*pos).type().name3() )
					<< A( 2, pose_copy_->pdb_info()->chain(*pos) );
			} else {
				pdbout << I( 4, *pos ) << A( 4, pose_copy_->residue(*pos).type().name3() )
					<< A( 2, pose_copy_->chain(*pos) );
			}
			EnergyMap const & res_energies( pose_copy_->energies().residue_total_energies(*pos) );
			for ( int st(1); st <= n_score_types; ++st ) {
				if ( weights[ ScoreType(st) ] == 0 ) continue;
				pdbout << F( colwidth, 3, res_energies[ ScoreType(st) ] * weights[ ScoreType(st) ] );
			}
			pdbout << '\n';
		}
		// end residues
		pdbout << "REMARK\n";
		// end per-residue energies
	}

	pdbout << "REMARK Weighted total energies:\n";
	EnergyMap const & total_energies( pose_copy_->energies().total_energies() );
	for ( int i(1); i <= n_score_types; ++i ) {
		Real const E( total_energies[ ScoreType(i) ] ),
			W( weights[ ScoreType(i) ] );
		//  if ( E == 0 || W == 0 ) continue;
		if ( W == 0 ) continue;
		Real const weighted_E( E * W );
		std::string const name( ScoreTypeManager::name_from_score_type( ScoreType(i) ) );
		pdbout << "REMARK " << A(16,name) << " " << F( 8, 2, weighted_E ) << '\n';
	}
	// "special" ('weight' value for total_score defaults to zero?)
	pdbout << "REMARK " << A(16,"total_score") << " " << F( 8,2,total_energies[total_score] ) << '\n';
	pdbout << "REMARK \n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Generates a table of hydrogen bond info and categorizes hbonds by type
/// @author ashworth
void
PDBOutput::output_hbond_info( ozstream & pdbout )
{
	if ( ! option[ OptionKeys::run::output_hbond_info ]() ) return;

	using namespace scoring::hbonds;
	HBondSet hbond_set;
	// omission can result in failed 'graph_state_ == GOOD' assertion, but pose is const...
	// pose.update_residue_neighbors();
	// this yields UNWEIGHTED energies
	fill_hbond_set( *pose_copy_, false, hbond_set );

	// first output column headers for file readability
	pdbout << "REMARK Loc, res, pos, pdb, chain, atom, x, y, z, res, pos, pdb, chain, atom, x, y, z, hbE, HAdis, xD, xH\n";

	// total protein-BASE hydrogen-bonding energy
	Real hb_spec( 0.0 );
	std::string category("NONE");

	DnaInterfaceFinder interface;
	interface.determine_protein_interface( *pose_copy_ );

	// expand list of residues to consider into nres bool vector for simple lookup
	vector1< bool > relevant_residue( pose_copy_->total_residue(), false );
	for ( vector1< Size >::const_iterator index( res_indices_to_output_.begin() ),
			end( res_indices_to_output_.end() ); index != end; ++index ) {
		relevant_residue[ *index ] = true;
	}

	for ( Size i(1); i <= hbond_set.nhbonds(); ++i ) {

		HBond const & hb( hbond_set.hbond(i) );
		Size const don_res_i( hb.don_res() ), acc_res_i( hb.acc_res() );
		Residue const & don_rsd( pose_copy_->residue( don_res_i ) ),
			acc_rsd( pose_copy_->residue( acc_res_i ) );
		// skip NA-NA
		if ( don_rsd.is_DNA() && acc_rsd.is_DNA() ) continue;
		// skip (protein) backbone-backbone hbonds
		if ( hb.don_hatm_is_protein_backbone() || hb.acc_atm_is_protein_backbone() ) continue;
		// skip hbonds that don't involve a relevant residue
		if ( !relevant_residue[ don_res_i ] && !relevant_residue[ acc_res_i ] ) continue;
		// skip protein residues that are not close to DNA
		if ( !don_rsd.is_DNA() && !interface.protein_neighbors().find( don_res_i )->second.close() ) {
			continue;
		}
		if ( !acc_rsd.is_DNA() && !interface.protein_neighbors().find( acc_res_i )->second.close() ) {
			continue;
		}

		Size const donH_i( hb.don_hatm() ), acc_i( hb.acc_atm() );
		std::string const don_hatm_name( don_rsd.atom_name( donH_i ) ),
			acc_atm_name( acc_rsd.atom_name( acc_i ) );

		Real weight(0.0);
		if ( acc_rsd.atom_is_backbone( acc_i ) || don_rsd.atom_is_backbone( donH_i ) ) {
			category = "sc_bb";
			weight = score_function_->weights()[ hbond_bb_sc ];
		} else {
			category = "sc_sc";
			weight = score_function_->weights()[ hbond_sc ];
		}
		if ( acc_rsd.is_DNA() || don_rsd.is_DNA() ) {
			category = ( category == "sc_bb" ? "dna_bb" : "dna_base" );
			// dna water adduct
			if ( acc_rsd.atom_type( acc_i ).name() == "HOH" ||
					don_rsd.atom_type( donH_i ).name() == "HOH"  ) {
				category = "dna_wat";
				weight = score_function_->weights()[ h2o_hbond ];
			}
		}
		Real const weighted_E( hb.energy() * weight );
		if ( category == "dna_base" || category == "dna_wat" ) hb_spec += weighted_E;

		// assume donH_i indexes a hydrogen whose only bonded neighbor is the donor heavy
		Size const don_i( don_rsd.bonded_neighbor( donH_i ).front() );

		// output information
		numeric::xyzVector< Real > const Hxyz( don_rsd.atom( donH_i ).xyz() ),
			Dxyz( don_rsd.atom( don_i ).xyz() ),
			Axyz( acc_rsd.atom( acc_i ).xyz() ),
			Bxyz( acc_rsd.atom( acc_rsd.atom_base( acc_i ) ).xyz() ),
			B2xyz( acc_rsd.atom( acc_rsd.abase2( acc_i ) ).xyz() );
		// AHdis, xD, xH : Acc-Hyd dist, angle A-H-D, angle B-A-H. used in scoring
		// copied from hbonds_geom::hb_energy_deriv
		Vector const HD( Dxyz - Hxyz ), AH( Hxyz - Axyz );
		Real const HDdis( HD.length() ), AHdis( AH.length() );
		Real const inv_HDdis( 1.0f / HDdis ), inv_AHdis( 1.0f / AHdis );
		Vector const HDunit( HD * inv_HDdis ), AHunit( AH * inv_AHdis );
		Vector BAunit, PBxyz;
		make_hbBasetoAcc_unitvector(
			hbond_set.hbond_options(),
			get_hbe_acc_hybrid(hb.eval_type()),
			Axyz, Bxyz, B2xyz, PBxyz, BAunit );
		Real const xD( dot( AHunit, HDunit ) ), xH( dot( BAunit, AHunit ) );

		int don_pdb(0), acc_pdb(0);
		if ( pose_copy_->pdb_info() ) {
			don_pdb = pose_copy_->pdb_info()->number( don_res_i );
			acc_pdb = pose_copy_->pdb_info()->number( acc_res_i );
		}
		pdbout << "REMARK " << A(8,category) << " ";
		pdbout << don_rsd.name3() << " " << I(4,don_res_i) << " " << I(4,don_pdb) << " ";
		if ( pose_copy_->pdb_info() ) pdbout << pose_copy_->pdb_info()->chain( don_res_i );
		else pdbout << pose_copy_->chain( don_res_i );
		pdbout << " " << don_hatm_name << " "
			// these coordinates are used by PyMOL plugins to create (and display) hydrogen bond cgo objects
			<< F(8,3,Hxyz[0]) << " " << F(8,3,Hxyz[1]) << " " << F(8,3,Hxyz[2]) << " ";
		pdbout << acc_rsd.name3() << " " << I(4,acc_res_i) << " " << I(4,acc_pdb) << " ";
		if ( pose_copy_->pdb_info() ) pdbout << pose_copy_->pdb_info()->chain( acc_res_i );
		else pdbout << pose_copy_->chain( acc_res_i );
		pdbout << " " << acc_atm_name << " "
			// these coordinates are used by PyMOL plugins to create (and display) hydrogen bond cgo objects
			<< F(8,3,Axyz[0]) << " " << F(8,3,Axyz[1]) << " " << F(8,3,Axyz[2]) << " ";
		pdbout << F(6,2, weighted_E ) << " " << AHdis << " " << xD << " " << xH << '\n';
	}
	// pdbout << "REMARK base-specific hbond energy: " << F(6,2,hb_spec) << '\n';
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details
/// @author ashworth
void
PDBOutput::output_buried_unsatisfied_hbonds( ozstream & pdbout )
{
	if ( !option[ OptionKeys::run::decoystats ]() ) return;

	using namespace protocols::toolbox::pose_metric_calculators;
	using id::AtomID_Map;
	BuriedUnsatisfiedPolarsCalculator bu_calc( "default", "default" );

	basic::MetricValue< AtomID_Map< bool > > map;
	bu_calc.get( "atom_bur_unsat", map, *pose_copy_ );
	for ( vector1< Size >::const_iterator pos( res_indices_to_output_.begin() ),
			end( res_indices_to_output_.end() ); pos != end; ++pos ) {
		for ( Size atomi(1); atomi <= map.value().n_atom(*pos); ++atomi ) {
			if ( map.value()(*pos,atomi) ) {
				ResidueType const & rt( pose_copy_->residue_type(*pos) );
				std::string name3( rt.name3() );
				if ( rt.is_DNA() ) name3 = dna_full_name3( name3 );

				pdbout << "REMARK Unsat ";
				if ( pose_copy_->pdb_info() ) {
					pdbout << pose_copy_->pdb_info()->chain(*pos) << '.'
						<< pose_copy_->pdb_info()->number(*pos);
				} else {
					pdbout << pose_copy_->chain(*pos) << '.' << *pos;
				}
				pdbout << '.' << name3 << " " << rt.atom_name(atomi) << std::endl;
			}
		}
	}

	// instantiation of BuriedUnsatisfiedPolarsCalculator looks up (or creates) instantiations of NumberHBondsCalculator and SasaCalculatorLegacy -- output of info from these too while we're at it
	// somehow pose has become aware of these calculators and will not recompute them(?)

	//  basic::MetricValue< id::AtomID_Map< Size > > atom_hbonds;
	//  pose_copy_->metric( bu_calc.name_of_hbond_calc(), "atom_Hbonds", atom_hbonds );
	basic::MetricValue< Size > nhbonds;
	pose_copy_->metric( bu_calc.name_of_hbond_calc(), "all_Hbonds", nhbonds );

	TR(t_info) << "pdb nhbonds " << nhbonds.value() << std::endl;
	pdbout << "REMARK nhbonds " << nhbonds.value() << '\n';

	//  basic::MetricValue< id::AtomID_Map< Real > > atom_sasa;
	//  pose_copy_->metric( bu_calc.name_of_sasa_calc(), "atom_sasa", atom_sasa );
	basic::MetricValue< Real > total_sasa;
	pose_copy_->metric( bu_calc.name_of_sasa_calc(), "total_sasa", total_sasa );

	TR(t_info) << "pdb total_sasa " << total_sasa.value() << std::endl;
	pdbout << "REMARK total_sasa " << total_sasa.value() << '\n';

	// packing statistical score
	PackstatCalculator calc_packstat;
	basic::MetricValue< Real > total_packstat;
	calc_packstat.get( "total_packstat", total_packstat, *pose_copy_ );
	TR(t_info) << "pdb total_packstat " << total_packstat.value() << std::endl;
	pdbout << "REMARK total_packstat " << total_packstat.value() << '\n';

}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
make_subdirs( std::string const & name )
{
	// for names that include directory paths: try to create any such directories that do not yet exist
	typedef utility::vector1< std::string > StringVec;
	StringVec const subdirs( string_split( name, '/' ) );
	for ( StringVec::const_iterator dir( subdirs.begin() ); dir != subdirs.end()-1; ++dir ) {
		if ( utility::file::file_exists( dir->c_str() ) ) continue;
		utility::file::create_directory( dir->c_str() );
	}
}

} // namespace dna
} // namespace protocols
