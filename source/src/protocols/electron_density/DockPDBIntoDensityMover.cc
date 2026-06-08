// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockPDBIntoDensityMover.cc
/// @brief protocols for folding into density
/// @details Re-tooled DockIntoDensityMover to support parallel execution, and more intermediate debugging steps
/// at the cost of only working with PDBs and not fragments.
/// @author Danny Farrell

#include <protocols/electron_density/DockPDBIntoDensityMover.hh>
#include <protocols/electron_density/DockIntoDensityUtils.hh>
#include <protocols/electron_density/util.hh>

#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>

#include <core/pose/init_id_map.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/AtomInformation.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/electron_density/xray_scattering.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/fourier/SHT.hh>

#include <ObjexxFCL/format.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>

#include <sys/stat.h>


using basic::Tracer;


namespace protocols {
namespace electron_density {

static basic::Tracer TR( "protocols.electron_density.DockPDBIntoDensityMover" );

//////////////////
//

MinimizePoseIntoDensityOptions
create_minimize_pose_into_density_options(core::pose::Pose const & pose, core::Real const dens_wt) {
	core::scoring::ScoreFunctionOP scorefxn_dens( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
	scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_allatom, 1.0);

	core::scoring::ScoreFunctionOP scorefxn_refine = core::scoring::get_score_function();
	scorefxn_refine->set_weight( core::scoring::elec_dens_fast, dens_wt);
	scorefxn_refine->set_weight( core::scoring::coordinate_constraint, 1.0);

	core::scoring::ScoreFunctionOP scorefxn_refine_rb( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
	scorefxn_refine_rb->set_weight( core::scoring::elec_dens_fast, dens_wt);
	scorefxn_refine_rb->set_weight( core::scoring::coordinate_constraint, 1.0);

	core::kinematics::MoveMapOP bbmm( utility::pointer::make_shared< core::kinematics::MoveMap >() );
	bbmm->set_bb( true ); bbmm->set_chi( true ); bbmm->set_jump( true );
	protocols::minimization_packing::MinMoverOP bbmin( utility::pointer::make_shared< protocols::minimization_packing::MinMover >(
		bbmm, scorefxn_refine, "lbfgs_armijo_nonmonotone", 1e-5, true ) );
	bbmin->max_iter(200); // make a parameter?

	// packer
	core::pack::task::TaskFactoryOP tf( utility::pointer::make_shared< core::pack::task::TaskFactory >() );
	tf->push_back( core::pack::task::operation::TaskOperationCOP( utility::pointer::make_shared< core::pack::task::operation::InitializeFromCommandline >() )); // get extra rotamer flags from command line
	tf->push_back( core::pack::task::operation::TaskOperationCOP( utility::pointer::make_shared< core::pack::task::operation::IncludeCurrent >() )); // include current rotamer by default
	tf->push_back( core::pack::task::operation::TaskOperationCOP( utility::pointer::make_shared< core::pack::task::operation::RestrictToRepacking >() )); // do not design
	protocols::minimization_packing::PackRotamersMoverOP packer( utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >() );
	packer->task_factory( tf );
	packer->score_function( scorefxn_refine );

	// Setup rigid-body movemap now! must allow jumps from virt(root) to move
	core::kinematics::MoveMapOP rbmm = core::kinematics::MoveMapOP( utility::pointer::make_shared< core::kinematics::MoveMap >() );
	rbmm->set_bb( false ); rbmm->set_chi( false ); rbmm->set_jump( false );
	int const root = pose.fold_tree().root();
	utility::vector1< core::kinematics::Edge > root_edges = pose.fold_tree().get_outgoing_edges(root);
	for ( core::Size i=1; i<=root_edges.size(); ++i ) rbmm->set_jump ( root_edges[i].label() , true );

	// Setup rigid-body min now!
	protocols::minimization_packing::MinMoverOP rbmin( utility::pointer::make_shared< protocols::minimization_packing::MinMover >( rbmm, scorefxn_refine_rb, "lbfgs_armijo_nonmonotone", 1e-5, true ) );
	rbmin->max_iter(200); // make a parameter?

	return MinimizePoseIntoDensityOptions{scorefxn_dens, scorefxn_refine, scorefxn_refine_rb, rbmin, bbmin, packer};
}
//

void
dump_to_silent( core::pose::Pose const & pose, std::string const & fname, std::string tag, bool overwrite ) {
	if ( overwrite ) remove( fname.c_str() );
	if ( tag.size() == 0 ) tag = "result_1";
	core::io::silent::SilentFileOptions sf_opt;
	core::io::silent::SilentFileData sfd( fname, false, false, "binary", sf_opt);
	core::io::silent::BinarySilentStruct silent_stream( sf_opt, pose, tag  );
	sfd.write_silent_struct( silent_stream, fname );
}


inline bool check_if_file_exists(const std::string& name) {
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}


void
DockPDBIntoDensityMover::set_nRsteps_from_pose( core::pose::Pose const & pose ) {
	numeric::xyzVector< core::Real > const com = core::pose::get_center_of_mass( pose );
	core::Real const extent = get_radius( pose);
	utility::vector1< core::Real > pose_1dspec;
	// grid spacing == delR
	core::Size ngrid = (core::Size) std::ceil( extent / delR_ + 2);
	pose_1dspec.resize(ngrid, 0.0);

	core::Real massSum=0.0;
	core::Size resstart = 1;
	core::Size resend = pose.size();
	if ( point_radius_ != 0 ) {
		core::Size midres = (pose.size()+1)/2;
		resstart = midres;
		resend = midres;
	}
	for ( core::Size i=resstart; i<=resend; ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			core::Real const binnum = ( atom.xyz()-com ).length() / delR_ + 1.0;
			core::Real const fpart = binnum - std::floor(binnum);
			core::Size const binint = (core::Size)std::floor(binnum);
			pose_1dspec[binint] += (1-fpart);
			pose_1dspec[binint+1] += (fpart);
			massSum += 1.0;
		}
	}

	// set nRsteps!
	core::Real const fracDens=fragDens_; // choose radius covering this fraction of density mass
	core::Real running_total=0.0;
	for ( int i=1; i<=(int)ngrid; ++i ) {
		running_total += pose_1dspec[i] / massSum;
		if ( running_total > fracDens ) {
			nRsteps_ = i;
			break;
		}
	}
}


void
DockPDBIntoDensityMover::setMultiNative( utility::vector1< core::pose::PoseOP > const & all_natives) {
	all_natives_ = all_natives;

	// first set COMs
	for ( core::Size i = 1; i <= all_natives_.size(); ++i ) {
		core::pose::addVirtualResAsRoot( (*(all_natives_[i])) );
		numeric::xyzVector< core::Real > current_native_com(numeric::xyzVector< core::Real >(0,0,0));
		core::Size N=0;
		for ( core::Size ii=1; ii <= all_natives_[i]->size(); ++ii ) {
			if ( all_natives_[i]->residue(ii).aa() == core::chemical::aa_vrt ) continue;
			if ( all_natives_[i]->residue(ii).is_protein() ) {
				current_native_com += all_natives_[i]->residue(ii).xyz(2);
				++N;
				if ( ii == ( all_natives_[i]->size() ) / 2 ) {
					all_native_seq_center_.push_back( all_natives_[i]->residue(ii).xyz(2) );
				}
			} else if ( all_natives_[i]->residue(ii).is_NA() ) {
				current_native_com += all_natives_[i]->residue(ii).xyz(1);
				++N;
			}
		}
		current_native_com /= N;
		all_native_com_.push_back( current_native_com );
	}

	// second set middle CAs. (middle ca is xyz of CA closest to COM)
	all_native_mca_.resize( all_native_com_.size(), numeric::xyzVector< core::Real >(0,0,0) );
	for ( core::Size i = 1; i <= all_natives_.size(); ++i ) {
		core::pose::addVirtualResAsRoot( (*(all_natives_[i])) );
		core::Real closest_distance_to_COM = 999;
		for ( core::Size ii=1; ii <= all_natives_[i]->size(); ++ii ) {
			if ( all_natives_[i]->residue(ii).aa() == core::chemical::aa_vrt ) continue;
			if ( all_natives_[i]->residue(ii).is_protein() ) {
				core::Real dist = all_natives_[i]->residue(ii).atom(2).xyz().distance( all_native_com_[i]);
				if (  dist < closest_distance_to_COM ) {
					all_native_mca_[i] = all_natives_[i]->residue(ii).atom(2).xyz(); // set mca
					closest_distance_to_COM = dist;
				}
			} else if ( all_natives_[i]->residue(ii).is_NA() ) {
				continue; // TODO in future, specify?
			}
		}
	}

	// All natives should be the same except orientation
	TR << "NOTE: distance from COM (in NATIVE) to nearest CA is " << all_native_mca_[1].distance( all_native_com_[1] ) << std::endl;
}


core::Real
DockPDBIntoDensityMover::get_radius( core::pose::Pose const & pose ) {
	// Not efficient because usually we want both radius & COM, but it's more useful/easier to read
	// if we just seperate the two.
	numeric::xyzVector< core::Real > pose_com = core::pose::get_center_of_mass( pose );
	core::Real maxSize = 0;
	for ( core::Size i=1; i<= pose.size(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			maxSize = std::max( maxSize, pose_com.distance( atom.xyz() ) );
		}
	}
	return maxSize; // this is the radius
}


void
DockPDBIntoDensityMover::predefine_search( utility::vector1< numeric::xyzVector<core::Real> > &pts_in ) {
	points_defined_ = true;

	points_to_search_.clear();

	for ( core::Size i=1; i<=pts_in.size(); ++i ) {
		numeric::xyzVector< core::Real > x_idx;
		core::scoring::electron_density::getDensityMap().cart2idx( pts_in[i] , x_idx );
		points_to_search_.push_back( x_idx );
	}
}


void
DockPDBIntoDensityMover::score_and_dump_natives() {
	using namespace core::io::silent;
	core::scoring::ScoreFunctionOP scorefxn_dens( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
	scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_allatom, 1.0);
	protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");
	std::string fname = silent_ +  "_native.silent";
	remove( fname.c_str() );
	SilentFileOptions sfopts;
	SilentFileData sfd( fname, false, false, "binary", sfopts );
	for ( core::Size i = 1; i <= all_natives_.size(); ++i ) {
		core::Real score = (*scorefxn_dens)(*all_natives_[i]);
		core::pose::Pose posecopy = (*all_natives_[i]);

		if ( centroid_silent_out_ ) {
			if ( !posecopy.is_centroid() ) {
				to_cen.apply( posecopy );
			}
		}

		std::string base_tag = utility::string_split( silent_, '/').back();
		std::string tag = base_tag + "_" + utility::to_string(1000 + i);

		BinarySilentStruct silent_stream( sfopts, posecopy, tag  );
		silent_stream.clear_energies();
		silent_stream.add_energy( "score", score ); // REQUIRED!! order is IMPORTANT
		silent_stream.add_energy( "dens_score", score );
		silent_stream.add_energy( "dens_rank", 1000 + i );
		silent_stream.add_energy( "rms", 0.0 );
		silent_stream.add_energy( "native_idx", i );
		sfd.write_silent_struct( silent_stream, fname );
	}
}


void
DockPDBIntoDensityMover::read_in_points_to_search() {
	points_to_search_.clear();
	std::string instring;
	// TODO: fix this to be an option!
	if ( points_to_search_fname_.size() != 0 ) { instring = points_to_search_fname_; }
	else if ( silent_.size() != 0 ) { instring = silent_ + ".pts"; }
	else { instring = "points_to_search.pts"; }

	if ( !check_if_file_exists( instring ) ) throw CREATE_EXCEPTION(utility::excn::BadInput, "Couldn't find the file " + instring);
	std::ifstream inpoints( instring.c_str() );
	numeric::xyzVector< core::Real > current_point_to_search(0.0, 0.0, 0.0 );
	core::Size current_point_idx;
	while ( inpoints >> current_point_idx >> current_point_to_search[0] >> current_point_to_search[1] >> current_point_to_search[2] ) {
		points_to_search_.push_back(current_point_to_search);
	}
	if ( points_to_search_.size() == 0 ) throw CREATE_EXCEPTION(utility::excn::BadInput, "Found 0 points to search. Please check the file " + instring );
}


core::Real
DockPDBIntoDensityMover::compare_and_align_poses( core::pose::Pose & query, core::pose::Pose const & native ) {
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, query, core::id::AtomID::BOGUS_ATOM_ID() );  // maps every atomid to bogus atom
	core::Size N = 0;
	for ( core::Size i = 1; i <= query.size(); ++i ) {
		int const query_resnum = query.pdb_info()->number(i);
		std::string const query_chain = query.pdb_info()->chain(i);
		core::Size const native_posenum = native.pdb_info()->pdb2pose( query_chain, query_resnum );
		//int r2_resnum = r2->pdb_info()->number(r2_posenum);
		if ( native_posenum == 0 ) continue;
		if ( query.residue(i).name1() != native.residue(native_posenum).name1() ) {
			TR.Fatal << "nat seq: " << native.sequence() << std::endl;
			TR.Fatal << "query seq: " << query.sequence() << std::endl;
			TR.Fatal << "nat chain: " << native.pdb_info()->chain(native_posenum) << " native pdbinfo res: " << native.pdb_info()->number(native_posenum) << std::endl;
			TR.Fatal << "query chain: " << query_chain << " query pdbinfo res: " << query_resnum << std::endl;
			TR.Fatal << "nat name: " << native.residue(native_posenum).name3() << " query name: " << query.residue(i).name3() << std::endl;
			TR.Fatal << "query pose num: " << i << " native pose num: " << native_posenum << std::endl;
			throw CREATE_EXCEPTION(utility::excn::BadInput, "Found that in compare_and_align_poses, poses did not align in residue numbering or sequence");
		}
		core::id::AtomID const id1( query.residue(i).atom_index("CA"), i );
		core::id::AtomID const id2( native.residue(native_posenum).atom_index("CA"), native_posenum );
		atom_map[ id1 ] = id2;
		++N;
	}
	if ( N == 0 ) {
		TR.Warning << "Was unable to align any residues between query pose and native. Maybe your native doesn't include this region?" << std::endl;
		return 1000;
	}
	core::Real const aligned_rms = core::scoring::superimpose_pose( query, native, atom_map );
	return aligned_rms;
}

RBfitResultDB
DockPDBIntoDensityMover::read_in_partial_search_results( utility::vector1< std::string > const & local_result_filenames ) const {
	if ( local_result_filenames.size() == 0 ) throw CREATE_EXCEPTION(utility::excn::BadInput, "Didn't find any local_result_filenames. Please check your -local_result_files flag");
	RBfitResultDB all_results( 10e8 );
	for ( auto &current_local_result_filename : local_result_filenames ) {
		std::ifstream inresults( current_local_result_filename.c_str() );
		// skip header
		std::string header;
		std::getline( inresults, header);
		// begin reading
		core::Real pose_idx, result_num, score, rot_xx, rot_xy, rot_xz, rot_yx, rot_yy;
		core::Real rot_yz, rot_zx, rot_zy, rot_zz, pre_x, pre_y, pre_z, post_x, post_y, post_z;
		while ( inresults >> pose_idx >>  result_num >>  score
				>>  rot_xx >> rot_xy >> rot_xz
				>>  rot_yx >> rot_yy >> rot_yz
				>>  rot_zx >> rot_zy >> rot_zz
				>>  pre_x >> pre_y >> pre_z
				>>  post_x >>  post_y >>  post_z  )
				{
			// generate variables
			numeric::xyzMatrix< core::Real > const rotation(
				numeric::xyzMatrix< core::Real >::cols(
				rot_xx, rot_yx, rot_zx,
				rot_xy, rot_yy, rot_zy,
				rot_xz, rot_yz, rot_zz) );
			numeric::xyzVector<core::Real> const pretrans( pre_x, pre_y, pre_z );
			numeric::xyzVector<core::Real> const posttrans( post_x, post_y, post_z );

			RBfitResult const single_result( pose_idx, score, rotation, pretrans, posttrans );
			all_results.add_element( single_result);
		}
	}
	if ( all_results.size() == 0 ) {
		std::stringstream ss;
		ss << "Found 0 results in the following files:" << std::endl;
		for ( auto & i : local_result_filenames ) {
			ss << i << " ";
		}
		ss << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput, ss.str() );
	}
	return all_results;
}


void
DockPDBIntoDensityMover::set_search_responsibilities() {
	// core_idx_[] 1 = current core, 2 is total cores
	if ( core_idx_.size() == 0 ) {
		point_search_start_ = 1;
		point_search_end_ = points_to_search_.size();
		TR << "Found no core idx settings, so running serially." << std::endl;
		return;
	}
	core::Size avg_workload = points_to_search_.size() / core_idx_[2]; // integers auto round down.
	point_search_start_ = ( core_idx_[1] - 1 ) * avg_workload + 1;
	point_search_end_ = ( core_idx_[1] ) * avg_workload;
	if ( core_idx_[1] == core_idx_[2] ) {
		point_search_end_ = points_to_search_.size();
	}
	TR << "This process will be responsible for searching from: " << point_search_start_ << " to " << point_search_end_ << std::endl;
}


void
DockPDBIntoDensityMover::dump_RBfitDB_to_pdbs( RBfitResultDB & resultDB, core::pose::PoseOP const & poseOP ) const {
	while ( resultDB.size() > 0 ) {
		core::Size const result_rank = resultDB.size();
		RBfitResult const sol_i = resultDB.pop();
		core::pose::PoseOP posecopy( utility::pointer::make_shared< core::pose::Pose >( *poseOP ) );
		apply_transform( *posecopy, sol_i );
		std::string const pdb_name = "result_" + std::to_string(result_rank) + ".pdb";
		posecopy->dump_pdb( pdb_name );
	}
}


void
DockPDBIntoDensityMover::import_refinement_silent_files( utility::vector1< std::string > const & refinement_result_filenames, RevRefinementResultDB & refinement_results ) {
	using namespace core::io::silent;
	for ( core::Size n=1; n<=refinement_result_filenames.size(); ++n ) {
		// load the domain placements into silent_file_data
		TR << "about to read in files from " << refinement_result_filenames[n] << std::endl;
		SilentFileOptions sfopts;
		SilentFileData silent_file_data( sfopts );
		silent_file_data.read_file( refinement_result_filenames[n] );

		// for every result in our silent file data
		for ( SilentFileData::iterator iter = silent_file_data.begin(), end = silent_file_data.end(); iter != end; ++iter ) {
			std::string tag = iter->decoy_tag();
			utility::vector1< std::string > fields = utility::string_split( tag, '_' );

			// if ( fields[1] == "empty" ) continue; // check?? TODO
			// if (fields[3] == "0001" ) continue; // this gets rid of jd2 output
			core::Real score = iter->get_energy( "dens_score" );
			core::pose::PoseOP imported_poseOP (utility::pointer::make_shared< core::pose::Pose >() );
			iter->fill_pose( *imported_poseOP );
			// core::Real rms = iter->get_energy( "rms" ); // TODO: maybe keep?
			// score, pose
			refinement_results.add_element( RefinementResult( score, imported_poseOP ) );
		}
	}
}


void
DockPDBIntoDensityMover::set_refinement_responsibilities( RBfitResultDB & resultDB ) {
	// core_idx_[] 1 = current core, 2 is total cores
	if ( core_idx_.size() == 0 ) {
		if ( refine_start_ == 0 || refine_end_ == 0 ) { // serial workload
			refine_start_ = 1;
			refine_end_ = resultDB.size();
		}
	} else { // auto generated workload
		core::Size avg_workload = resultDB.size() / core_idx_[2]; // integers auto round down.
		refine_start_ = ( core_idx_[1] - 1 ) * avg_workload + 1;
		refine_end_ = ( core_idx_[1] ) * avg_workload;
		if ( core_idx_[1] == core_idx_[2] ) {
			refine_end_ = resultDB.size();
		}
	} // if none of the above, we assume it was manually set
	TR << "This process will be responsible for refining from: " << refine_start_ << " to " << refine_end_ << std::endl;
}


void
DockPDBIntoDensityMover::check_for_existing_output_file( std::string const & mode ) const {
	// .user? ? maybe instead of overwrite true?
	if ( mode == "refine" || mode == "combine_refine" || mode == "cluster_silent" ) {
		if ( check_if_file_exists( silent_ ) ) {
			if ( !overwrite_ ) {
				std::string const error_text = [&]{
					std::stringstream tmp_error_text;
					tmp_error_text << "Found that file: "  << silent_ <<  " already exists!\nDid you forget to pass -overwrite true? ";
					return tmp_error_text.str();
				}();
				throw CREATE_EXCEPTION(utility::excn::IOError, error_text);
			} else {
				std::remove(silent_.c_str());
			}
		}
	} else {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Check for existing output file called a string that wasn't listed in the options" );
	}
}


void
DockPDBIntoDensityMover::minimize_poseOP_into_density(
	core::pose::PoseOP const posecopy,
	MinimizePoseIntoDensityOptions const & params,
	utility::vector1< core::Real > & scores_vector
) {

	core::pose::addVirtualResAsRoot( *posecopy );
	if ( constrain_refinement_ != 0 ) {
		core::scoring::constraints::ConstraintSetOP const cst_set( utility::pointer::make_shared< core::scoring::constraints::ConstraintSet >() );
		core::Size const root_i = posecopy->fold_tree().root();
		if ( posecopy->residue( root_i ).aa() == core::chemical::aa_vrt && root_i == posecopy->total_residue() ) {
			for ( core::Size idx = 1; idx <= posecopy->total_residue(); ++idx ) {
				core::conformation::Residue const & rsd( posecopy->residue( idx ) );
				for ( core::Size ii = 2; ii <= 2; ++ii ) {
					core::scoring::func::FuncOP fx( utility::pointer::make_shared< core::scoring::func::FlatHarmonicFunc >( 0.0, 0.05, constrain_refinement_));
					cst_set->add_constraint( core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP(
						utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >(
						core::id::AtomID( ii, idx ), core::id::AtomID( 1, posecopy->total_residue() ), rsd.xyz( ii ), fx ) ) ) );
				}
			}
		}
		posecopy->constraint_set( cst_set );
	}

	core::Real const scoreb = (*params.scorefxn_dens)(*posecopy);
	core::Real scorei=0;
	numeric::xyzVector< core::Real > start_com = core::pose::get_center_of_mass( *posecopy );
	scores_vector.push_back( scoreb );

	protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD );
	// rbmin/pack/fullmin
	try {
		if ( do_refine_ ) {
			if ( dump_inter_silent_ ) dump_to_silent( *posecopy, "prerb.silent", "prerb", true );
			if ( dump_inter_ ) posecopy->dump_pdb("prerb.pdb");
			params.rbmin->apply( *posecopy );
			if ( dump_inter_silent_ ) dump_to_silent( *posecopy, "postrb.silent", "postrb", true );
			if ( dump_inter_ ) posecopy->dump_pdb("postrb.pdb");

			if ( posecopy->is_centroid() ) {
				to_all_atom.apply( *posecopy );
			}
			if ( dump_inter_silent_ ) dump_to_silent( *posecopy, "_testmin.silent", "_testmin", true );

			for ( core::Size i=1; i<=ncyc_; ++i ) {
				params.packer->apply( *posecopy );
				scorei = (*params.scorefxn_dens)(*posecopy);

				if ( min_backbone_ ) {
					//scorefxn_refine->show( *posecopy );
					params.bbmin->apply( *posecopy );
					//scorefxn_refine->show( *posecopy );
				} else {
					params.rbmin->apply( *posecopy );
				}
			}
		}
	} catch (...) {
		TR << "ERROR FAILURE (2)! WARNING, LOSING A POSE TO MINIMIZATION ERROR" << std::endl;
		scores_vector.push_back( 10000.00 );
		scores_vector.push_back( 10000.00 );
		return;
	}

	numeric::xyzVector< core::Real > const end_com = core::pose::get_center_of_mass( *posecopy );
	// TODO THIS IS A HORRIBLE HACK
	// Sometimes poses are thrown to inifinity, this hopefully ignores those poses
	// https://github.com/RosettaCommons/main/pull/2728
	if ( start_com.distance(end_com) > 100 ) {
		TR << "ERROR FAILURE! WARNING, LOSING A POSE TO MINIMIZATION ERROR" << std::endl;
		scores_vector.push_back( 10000.00 );
		scores_vector.push_back( 10000.00 );
		return;
	}
	if ( false ) ( *params.scorefxn_refine)(*posecopy);
	scores_vector.push_back( scorei );
	// rescore (expensive sf)
	core::Real const scoref = (*params.scorefxn_dens)(*posecopy);
	scores_vector.push_back( scoref );
}


void
DockPDBIntoDensityMover::refine_RBfitResultDB(
	core::pose::PoseOP const &poseOP,
	RBfitResultDB & results_in,
	RevRefinementResultDB & results_out
) {
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	// do partial refinement
	// copy to vector and then act on only a range of the vector indicies
	// refine best results first :D
	core::Size const ntotal = results_in.size();
	utility::vector1< RBfitResult > results_in_vector(ntotal); // idk if this works

	while ( results_in.size() > 0 ) {
		RBfitResult sol_i = results_in.pop();
		results_in_vector[ ( results_in.size() + 1 ) ] = sol_i;
	}
	TR << "Refining " << refine_end_ - refine_start_ + 1 << " poses" << std::endl;
	if ( constrain_refinement_ != 0 ) {
		TR << "Note: We will be adding Coordinate Constraints" << std::endl;
	}

	for ( core::Size i = refine_start_; i <= refine_end_; ++i ) {
		RBfitResult const sol_i = results_in_vector[i];
		//TR << "[" << ntotal-results_in.size() << "/" << ntotal << "]" << "\r" << std::flush;

		core::pose::PoseOP posecopy ( utility::pointer::make_shared< core::pose::Pose >( *poseOP ) );

		if ( dump_inter_ ) posecopy->dump_pdb( utility::to_string(i) + "_pretrans.pdb");
		if ( dump_inter_silent_ ) dump_to_silent( *posecopy, utility::to_string(i) + "_pretrans.silent", utility::to_string(i) + "_pretrans", true );

		apply_transform( *posecopy, sol_i );
		core::pose::addVirtualResAsRoot( *posecopy );
		utility::vector1< core::Real > scores_vector;

		if ( dump_inter_ ) posecopy->dump_pdb( utility::to_string(i) + "_premin.pdb");
		if ( dump_inter_silent_ ) dump_to_silent( *posecopy, utility::to_string(i) + "_premin.silent", utility::to_string(i) + "_premin", true );

		// check RMS, maybe future COM check first before you get_rms, but usually don't have that many domains, so it's not really a big slowdown.
		core::Size best_native_index = 0;
		core::Real rms_b=9999.0, gdt_before=0;
		if ( all_natives_.size() != 0 ) {
			for ( core::Size native_i = 1; native_i <= all_natives_.size(); ++native_i ) {
				core::Real rms = get_rms( *posecopy, *all_natives_[ native_i ], symminfo_, true);
				if ( rms < rms_b ) {
					rms_b = rms;
					best_native_index = native_i;
				}
			}
			gdt_before = get_gdt( *posecopy, *all_natives_[best_native_index], symminfo_, true);
		}
		//////// MINIMIZE LOGIC IN HERE ///////////
		minimize_poseOP_into_density( posecopy, create_minimize_pose_into_density_options(*posecopy, dens_wt_), scores_vector );
		if ( dump_inter_ ) posecopy->dump_pdb( utility::to_string(i) + "_postmin.pdb");
		if ( dump_inter_silent_ ) dump_to_silent( *posecopy, utility::to_string(i) + "_postmin.silent", utility::to_string(i) + "_postmin", true );

		// only do rms check on our closest native from our last check
		if ( all_natives_.size() != 0 ) {
			if ( dump_inter_ ) posecopy->dump_pdb( utility::to_string(i) + "_xx1ostmin.pdb");
			if ( dump_inter_ ) all_natives_[best_native_index]->dump_pdb( utility::to_string(i) + "_xx2ostmin.pdb");
			core::Real rms = get_rms( *posecopy, *all_natives_[best_native_index], symminfo_, true );
			core::Real gdt_after = get_gdt( *posecopy, *all_natives_[best_native_index], symminfo_, true );
			TR << "[" << i << "/" << ntotal << "] " << scores_vector[1] << " : " << scores_vector[2] << " : " << scores_vector[3] << " rms before=" << rms_b << " rms after=" << rms <<  " gdt before: " << gdt_before << " gdt after: " << gdt_after << std::endl;
		} else {
			TR << "[" << i << "/" << ntotal << "] " << scores_vector[1] << " : " << scores_vector[2] << " : " << scores_vector[3] <<  std::endl;
		}
		// TODO CONTINUE HORRIBLE HACK
		if ( scores_vector[3] == 10000.00 ) {
			TR << "CONTINUE HORRIBLE HACK ON RESULT: " <<  i << std::endl;
			core::pose::PoseOP tmppose ( utility::pointer::make_shared< core::pose::Pose >( *poseOP ) );
			results_out.add_element( RefinementResult( scores_vector[3], tmppose ) );
			continue;
		}
		// store
		results_out.add_element( RefinementResult( scores_vector[3], posecopy ) );
	}
}


void
DockPDBIntoDensityMover::cluster_RevRefinementDB( RevRefinementResultDB & results, core::Size target_size = 0 ) {
	// dumb -- copy to vector
	utility::vector1< RefinementResult > results_sort;
	while ( results.size() > 0 )
			results_sort.push_back( results.pop() );
	utility::vector1< bool > selector(results_sort.size(), false);

	for ( int i=results_sort.size(); i>=1; --i ) {
		selector[i] = true;
		for ( int j=i+1; j<=(int)results_sort.size() && selector[i]; ++j ) {
			if ( selector[j] && get_rms(*results_sort[i].pose_, *results_sort[j].pose_, symminfo_, false) < cluster_radius_ ) {
				selector[i] = false;
			}
		}
		if ( selector[i] ) {
			results.add_element( results_sort[i] );
		}
		if ( results.size() == target_size && target_size != 0 ) return;
	}
}

void
DockPDBIntoDensityMover::analyze_RefinementDB(RevRefinementResultDB resultDB) {
	TR.Debug << "START analyze_refinementDB here" << std::endl;
	while ( resultDB.size() > 0 ) {
		core::Size rank = resultDB.size();
		RefinementResult sol_i = resultDB.pop();
		core::Real best_rms = 9999.0;
		core::Size best_rms_idx = 0;
		if ( all_natives_.size() != 0 ) {
			for ( core::Size i = 1; i <= all_natives_.size(); ++i ) {
				core::Real rms = get_rms( *sol_i.pose_, *all_natives_[i], symminfo_, true );
				if ( rms < best_rms ) {
					best_rms = rms;
					best_rms_idx = i;
				}
			}
		}
		TR.Debug << "rank: " << rank << " score: " << sol_i.score_ << " rms: " << best_rms << " nat idx: " << best_rms_idx << std::endl;
	}
	TR.Debug << "END analyze_refinmentDB" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Start main callable functions /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO: big- should make all innner functions private. or protected
// TODO: add a header check that will check to see if we have all points (just a warning, no error because we want to be able to do things partially.)


utility::vector1< numeric::xyzVector<core::Real> >
DockPDBIntoDensityMover::get_points_to_search( core::pose::PoseOP const poseOP ) {
	if ( score_natives_ ) {
		TR << "Scoring natives and dumping to silent..." << std::endl;
		score_and_dump_natives();
	}
	if ( cheat_native_com_ ) {
		TR << "Predefining search via native COMs..." << std::endl;
		predefine_search( all_native_com_ );
	}
	if ( cheat_native_mca_ ) {
		TR << "Predefining search via native MCAs..." << std::endl;
		predefine_search( all_native_mca_ );
	}
	if ( start_model_.size() != 0 ) {
		core::pose::PoseCOP poseOP(core::import_pose::pose_from_file(start_model_));
		protocols::electron_density::remove_occupied_density_from_density(*poseOP, core::scoring::electron_density::getDensityMap());
	}
	if ( !(cheat_native_com_ || cheat_native_mca_) ) {
		TR << "Selecting points to search..." << std::endl;
		points_to_search_ = select_density_points(
			*poseOP,
			SelectDensityPointsOptions{
			B_, nRsteps_, gridStep_, topNtrans_,
			delR_, fragDens_, point_radius_, laplacian_offset_,
			rot_middle_ca_, convolute_single_residue_,
			symminfo_, all_natives_, all_native_com_, all_native_mca_
			},
			nRsteps_
		);
	}
	TR << "Writing points to search..." << std::endl;
	dump_points_to_search_to_pdb_or_txt(points_to_search_, points_to_search_pdb_fname_, points_to_search_fname_);
	return points_to_search_;
}


RBfitResultDB
DockPDBIntoDensityMover::apply_search( core::pose::Pose & pose, core::Size const result_size ) {
	if ( points_to_search_.empty() ) read_in_points_to_search();
	set_search_responsibilities();
	RBfitResultDB results(result_size); // placeholder

	{
		utility::vector1< core::Real > pose_1dspec;
		nRsteps_ = get_spectrum( pose, pose_1dspec, delR_, fragDens_, convolute_single_residue_, rot_middle_ca_ );
	}

	density_grid_search(
		1, pose, results,
		points_to_search_,
		DensityGridSearchOptions{
		B_, nRsteps_, max_rot_per_trans_, point_search_start_, point_search_end_,
		delR_, cluster_radius_, laplacian_offset_, 5.0,
		rot_middle_ca_, convolute_single_residue_, true,
		point_search_results_fname_,
		symminfo_, all_natives_, all_native_com_, all_native_mca_
		}
	);

	return results;
}


void
DockPDBIntoDensityMover::combine_search( utility::vector1< std::string > const & local_result_filenames, core::pose::Pose const & pose, RBfitResultDB & all_results ) {
	if ( all_results.size() == 0 ) all_results = read_in_partial_search_results( local_result_filenames );
	TR << "Imported " << all_results.size() << " results from fast density search" << std::endl;
	cluster_RBfitResultDB_fast(
		all_results,
		delR_,
		nRsteps_,
		cluster_radius_,
		topNfilter_,
		true,
		core::scoring::electron_density::getDensityMap());
	TR << "After clustering we have " << all_results.size() << " results from fast density search" << std::endl;
	if ( all_natives_.size() != 0 ) {
		TR << "Comparing our DB results to the native... " << std::endl;
		compare_RBfitDB_to_native( all_results, pose, all_natives_, all_native_com_, all_native_mca_, symminfo_, rot_middle_ca_, 15.0);
	}

	std::string const myfile = [&]{
		if ( combined_search_results_fname_.empty() ) return std::string("results_for_refinement.srfr");
		else return combined_search_results_fname_;
	}();
	std::ofstream outresults( myfile.c_str() );
	write_RBfitResultDB( all_results, outresults );
}


void
DockPDBIntoDensityMover::search_results_to_pdb( utility::vector1< std::string > const & local_result_filenames, core::pose::PoseOP const & poseOP ) {
	RBfitResultDB all_results = read_in_partial_search_results( local_result_filenames );
	dump_RBfitDB_to_pdbs( all_results, poseOP );
}


RevRefinementResultDB
DockPDBIntoDensityMover::apply_refinement(
	utility::vector1< std::string > const & local_result_filenames,
	core::pose::PoseOP const poseOP,
	RBfitResultDB & results_to_refine,
	bool const write_to_file
) {
	if ( results_to_refine.size() == 0 ) results_to_refine = read_in_partial_search_results( local_result_filenames );
	while ( results_to_refine.size() > topNfilter_ ) results_to_refine.pop();
	RevRefinementResultDB refinement_results( topNfilter_ ); // keep same as results_to_refine size!
	set_refinement_responsibilities( results_to_refine );
	check_for_existing_output_file( "refine" );
	refine_RBfitResultDB( poseOP, results_to_refine, refinement_results );

	std::string const out_silent_file = [&]{
		if ( silent_.empty() ) {
			return  "final_result_inter_"
				+ ObjexxFCL::right_string_of(refine_start_, 6, '0') +
				"_" + ObjexxFCL::right_string_of(refine_end_, 6, '0') + ".ref_silent";
		} else { return silent_; }
	}();

	if ( write_to_file ) {
		dump_RefinementDB_to_silent(
			refinement_results,
			silent_,
			out_silent_file,
			final_chain_,
			/*centroid_output*/ true,
			/*appent_to_outfile*/ false,
			all_natives_,
			symminfo_,
			/* legacy_rms */ false);
	}
	return refinement_results;
}

void
DockPDBIntoDensityMover::manual_refine_pdb( core::pose::PoseOP const & poseOP ) {
	utility::vector1< core::Real > scores_vector;
	minimize_poseOP_into_density( poseOP, create_minimize_pose_into_density_options(*poseOP, dens_wt_), scores_vector );
	poseOP->dump_pdb("manual_refined.pdb");
}


void
DockPDBIntoDensityMover::combine_refinement( utility::vector1< std::string > const & refinement_result_filenames, RevRefinementResultDB & refinement_results ) {
	TR << "Running combine refinement, will output max " << topNfilter_ << " results!" << std::endl;
	check_for_existing_output_file( "combine_refine" );
	if ( refinement_results.size() == 0 ) {
		refinement_results = RevRefinementResultDB( topNfilter_ * 5 );
		import_refinement_silent_files( refinement_result_filenames, refinement_results);
	}

	cluster_RevRefinementDB( refinement_results, topNfilter_ );
	TR << "Just clustered and have " << refinement_results.size() << " results remaining!" << std::endl;
	// only dump top N
	while ( refinement_results.size() > topNfinal_ ) { refinement_results.pop(); }

	std::string const out_silent_file = silent_.empty() ? "final_result.silent" : silent_;

	dump_RefinementDB_to_silent(
		refinement_results,
		silent_,
		out_silent_file,
		final_chain_,
		/*centroid_output*/ true,
		/*appent_to_outfile*/ false,
		all_natives_,
		symminfo_,
		/* legacy_rms */ false);
}

void
DockPDBIntoDensityMover::cluster_silent( utility::vector1< std::string > const & final_result_filenames ) {
	TR << "Running cluster silent, will output max " << topNfilter_ << " results!" << std::endl;
	check_for_existing_output_file( "cluster_silent" );
	RevRefinementResultDB final_results( topNfilter_ * 5 );
	import_refinement_silent_files( final_result_filenames, final_results);
	analyze_RefinementDB(final_results);
	TR << "haven't clustered yet" << std::endl;
	cluster_RevRefinementDB( final_results, topNfilter_ );
	TR << "Just clustered and have " << final_results.size() << " results remaining!" << std::endl;
	// only dump top N
	while ( final_results.size() > topNfinal_ ) { final_results.pop(); }

	std::string const out_silent_file = [&]{
		if ( silent_.empty() ) return std::string("clustered.silent");
		else return silent_ + "_clustered.silent";
	}();

	dump_RefinementDB_to_silent(
		final_results,
		silent_,
		out_silent_file,
		final_chain_,
		/*centroid_output*/ true,
		/*appent_to_outfile*/ false,
		all_natives_,
		symminfo_,
		/* legacy_rms */ false);
}


void
DockPDBIntoDensityMover::run_aio( core::pose::PoseOP const poseOP ) {
	core_idx_ = utility::vector1< core::Size >({1, 1});
	get_points_to_search( poseOP );

	point_search_start_ = 1;
	point_search_end_ = points_to_search_.size();

	RBfitResultDB search_results = apply_search( *poseOP, 1000000 );
	combine_search(utility::vector1< std::string >(), *poseOP, search_results);
	RevRefinementResultDB refined_results = apply_refinement(utility::vector1< std::string >(), poseOP, search_results, false);
	combine_refinement(utility::vector1< std::string >(), refined_results);
	analyze_RefinementDB(refined_results);
}



void
DockPDBIntoDensityMover::apply( core::pose::Pose & pose ) {
	TR << "Todo: not implemented yet " << std::endl;
	TR << pose.sequence() << std::endl;
	// brief: run everything in serial
	// get filenames from previous functions, and
}



}
}

