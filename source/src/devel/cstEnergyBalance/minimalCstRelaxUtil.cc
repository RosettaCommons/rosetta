// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.AlignmentCluster

/// @file   apps/pilot/brunette/minimalCstRelaxUtil.cc
///
/// @brief  This is the utility for minimalCstRelax and minimalCstHomology
/// @author TJ Brunette (tjbrunette@gmail.com)

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <core/types.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
//#include <core/util/SwitchResidueTypeSet.hh>

#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/id/SequenceMapping.hh>
//#include <core/io/silent/SilentStruct.hh>
//#include <core/io/silent/SilentFileData.hh>
//#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//#include <core/import_pose/pose_stream/util.hh>
//#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/kinematics/FoldTree.hh>

// AUTO-REMOVED #include <protocols/comparative_modeling/util.hh>
// AUTO-REMOVED #include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/comparative_modeling/ThreadingMover.hh>
// AUTO-REMOVED #include <protocols/comparative_modeling/PartialThreadingMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MissingDensityToJumpMover.hh>
#include <protocols/simple_moves/MissingDensityToJumpMover.fwd.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/cst_util.hh>

//utilities
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>


#include <devel/cstEnergyBalance/minimalCstRelaxUtil.hh>
//#include <core/init/init.hh>
// option key includes
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/run.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>

#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <ostream>

#include <utility/vector1.hh>


namespace devel {
namespace cstEnergyBalance {

static thread_local basic::Tracer tr( "minimalCstRelaxUtil" );

using core::pose::Pose;
using core::Size;
using core::Real;
using std::string;
using std::map;
using utility::vector1;
using namespace core::sequence;

/// @brief Gets the x,y,z center of mass of a pose
numeric::xyzVector< Real > get_centerOfMass(const Pose& pose ){
	using namespace core::conformation;
	Size nAtms=0;
	numeric::xyzVector< Real > massSum(0.0,0.0,0.0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if (rsd.aa() == core::chemical::aa_vrt) continue;
		for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			massSum += atom.xyz();
			nAtms++;
		}
	}
	massSum /= nAtms;
	return massSum;
}

///@brief Gets the x,y,z center of mass from a partial thread.
numeric::xyzVector< Real > get_centerOfMass(const Pose& templatePose, const SequenceAlignment aln,const string query_sequence){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::sequence;
	using namespace protocols;
	numeric::xyzVector< core::Real > atomLocation(0.0,0.0,0.0);
	core::pose::Pose query_pose;
	core::pose::make_pose_from_sequence(
																			query_pose, query_sequence, *( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )) );
	comparative_modeling::ThreadingMover mover( aln, templatePose);
	mover.build_loops( false );
	mover.randomize_loop_coords( option[OptionKeys::nonlocal::randomize_missing]() );
	mover.repack_query( false );
	mover.apply( query_pose );

	addVirtualResAsRoot(query_pose);
	atomLocation = query_pose.residue(query_pose.fold_tree().root()).xyz(1);
	return(atomLocation);
}

/// @brief Generates a set of coordinate constraints that correspond to the
/// CA that are being constrained.
core::scoring::constraints::ConstraintSetOP convert_caAtomsToConstrain_to_coordCsts(std::set< Size > caAtomsToConstrain,const Pose& pose){
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::conformation;
	using core::id::SequenceMapping;
	using core::id::AtomID;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::set< Size >::iterator caAtomsToConstrain_start,caAtomsToConstrain_stop;
	Real const default_coord_sdev( option[ OptionKeys::relax::minirelax_sdev ]());
	static string atom_name( "CA" );
	ConstraintSetOP cst_set( new ConstraintSet );
	caAtomsToConstrain_start =caAtomsToConstrain.begin();
	caAtomsToConstrain_stop = caAtomsToConstrain.end();
	tr.Debug << "---------------------------------------------------" << std::endl;
	while(caAtomsToConstrain_start != caAtomsToConstrain_stop){
		Size residue_id = *caAtomsToConstrain_start;
		tr.Debug << "constraining " << residue_id << std::endl;
		Residue const & rsd( pose.residue(residue_id) );
		cst_set->add_constraint(new CoordinateConstraint(AtomID(rsd.atom_index("CA"),residue_id), AtomID(1,pose.total_residue()), rsd.xyz(atom_name),new core::scoring::func::HarmonicFunc(0.0,default_coord_sdev)));
		caAtomsToConstrain_start++;
	}
return(cst_set);
}

/// @brief Generates a list of core residues to constrain. The residues chosen are the number of residues to constrain closest to the center of mass
std::set<Size> get_coreResiduesToConstrain(const Size numbResiduesToConstrain,const Pose& pose){
	std::set< Size > residuesToConstrain;
	numeric::xyzVector< Real > massSum = get_centerOfMass(pose);
	std::multimap<Real,Size> distCenterMassCaPos;
	std::multimap<Real,Size>::iterator start_massCenterCaPos,stop_massCenterCaPos;
	static string atom_name( "CA" );
	for(Size ii = 1; ii <= pose.total_residue(); ++ii ){
		Real distance = pose.residue(ii).xyz(atom_name).distance(massSum);
		distCenterMassCaPos.insert(std::pair<Real,Size>(distance,ii));
	}
	Size added_numbCsts = 0;
	start_massCenterCaPos = distCenterMassCaPos.begin();
	stop_massCenterCaPos = distCenterMassCaPos.end();
	while((start_massCenterCaPos != stop_massCenterCaPos)&&(added_numbCsts<numbResiduesToConstrain)){
		residuesToConstrain.insert(start_massCenterCaPos->second);
		tr.Debug << "mass center csts being added" << std::endl;
		tr.Debug <<start_massCenterCaPos->first << "," << start_massCenterCaPos->second << std::endl;
		start_massCenterCaPos++;
		added_numbCsts++;
	}
	return residuesToConstrain;
}

/// @brief Generates a list of residues to constrain based on movement between relaxed pose and an unmodified pose.
std::set<Size> get_distDeviationResiduesToConstrain(const Real distDeviationThresh,Size gapBtwConstrainedResidues, Pose& relaxed_pose, const Pose& unmodified_pose){
	std::set< Size > residuesToConstrain;
	static string atom_name( "CA" );
	core::scoring::calpha_superimpose_pose(relaxed_pose, unmodified_pose);
	Size count = 0;
	tr.Debug << "---------------------------------------------------" << std::endl;
	for ( Size ii = 1; ii <= unmodified_pose.total_residue(); ++ii ) {
		Real distance = relaxed_pose.residue(ii).xyz(atom_name).distance( unmodified_pose.residue(ii).xyz(atom_name));
		tr.Debug << ii << "-" << distance << std::endl;
		if(distance >= distDeviationThresh){
			if(count == 0){
				tr.Debug << "adding cst " << ii << std::endl;
				residuesToConstrain.insert(ii);
				count = gapBtwConstrainedResidues; //counts down from gap size
			}
			else
				count--;
		}
	}
	return residuesToConstrain;
}

///@brief Generates a list of residues to constrain based on the center of the 3 residues clossest to the center of mass of contiguous regions
std::set<Size> get_coreDistDeviationResiduesToConstrain(const Real distDeviationThresh, Pose& relaxed_pose, const Pose& unmodified_pose){
	//get all residues to constrain
	std::set< Size >::iterator resSet_iter,tmp_resSet_iter;
	std::set< Size > final_residuesToConstrain;
	std::set< Size > all_residuesToConstrain = get_distDeviationResiduesToConstrain(distDeviationThresh,0,relaxed_pose,unmodified_pose);
	Size CORE_RES_TO_CONSTRAIN = 3;
	//get constrained regions longer than 3 residues
	if(all_residuesToConstrain.size() <= 3)
		return all_residuesToConstrain;
	else{
		Size firstContigRes,lastContigRes;
		std::set< Size > tmp_residuesToConstrain;
		resSet_iter = all_residuesToConstrain.begin();
		firstContigRes = *resSet_iter;
		lastContigRes = *resSet_iter;
		tmp_residuesToConstrain.insert(*resSet_iter);
		resSet_iter++;
		while(resSet_iter != all_residuesToConstrain.end()){
			if(*resSet_iter == lastContigRes+1){
				lastContigRes = *resSet_iter;
				tmp_residuesToConstrain.insert(*resSet_iter);
			}
			else{
				if(lastContigRes-firstContigRes >= 3){
					//The case where the residues were contiguous
					tmp_residuesToConstrain.clear();
					core::pose::PoseOP	subset_poseOP = new core::pose::Pose(unmodified_pose,firstContigRes,lastContigRes);
					tmp_residuesToConstrain = get_coreResiduesToConstrain(CORE_RES_TO_CONSTRAIN,*subset_poseOP);
					tmp_resSet_iter = tmp_residuesToConstrain.begin();
					while(tmp_resSet_iter != tmp_residuesToConstrain.end()){
						final_residuesToConstrain.insert(*tmp_resSet_iter+firstContigRes-1);
						tmp_resSet_iter++;
					}
				}
				else{
					//The case where the residues were not contiguous
					tmp_resSet_iter = tmp_residuesToConstrain.begin();
					while(tmp_resSet_iter != tmp_residuesToConstrain.end()){
						final_residuesToConstrain.insert(*tmp_resSet_iter);
						tmp_resSet_iter++;
					}
				}
				//after contiguous? region reset the position index
				firstContigRes = *resSet_iter;
				lastContigRes = *resSet_iter;
				tmp_residuesToConstrain.clear();
				tmp_residuesToConstrain.insert(*resSet_iter);
			}
			resSet_iter++;
		}
		//post loop wrap up the final two cases.
		if(lastContigRes-firstContigRes >= 3){
			//The case where the residues were contiguous
			tmp_residuesToConstrain.clear();
			core::pose::PoseOP	subset_poseOP = new core::pose::Pose(unmodified_pose,firstContigRes,lastContigRes);
			tmp_residuesToConstrain = get_coreResiduesToConstrain(CORE_RES_TO_CONSTRAIN,*subset_poseOP);
			tmp_resSet_iter = tmp_residuesToConstrain.begin();
			while(tmp_resSet_iter != tmp_residuesToConstrain.end()){
				final_residuesToConstrain.insert(*tmp_resSet_iter+firstContigRes-1);
				tmp_resSet_iter++;
			}
		}
		else{
			//The case where the residues were not contiguous
			tmp_resSet_iter = tmp_residuesToConstrain.begin();
			while(tmp_resSet_iter != tmp_residuesToConstrain.end()){
				final_residuesToConstrain.insert(*tmp_resSet_iter);
				tmp_resSet_iter++;
			}
		}
		return(final_residuesToConstrain);
	}
}

/// @brief This util uses the above utilities to create a list of which CA residues to constrain.*Pose is not const because gaps of more than 3.5 angstroms are changed to be jumps
/// First loop through it only constrains the core residues
std::set<Size> get_residuesToConstrain(const Size /*coordCstGapInitial*/, const Real gdtThresh, Pose& pose){
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace protocols;
	using namespace relax;
	const Real DIST_DEVIATION_THRESH_MIN = 1.0; //max dist deviation
	const Real DIST_DEVIATION_THRESH_START = 2.5;
	const Real DIST_DEVIATION_THRESH_STEP = 0.5;
	const Real COORDINATE_CST_WT = 1;
	const Size CORE_RESIDUES_TO_CONSTRAIN = 3;
	const Size RDS_OF_RELAX = 2; //default is 2
	// Size coordCstGap = coordCstGapInitial; // Unused variable causes warning.
	//Get input pdbs-------------------------------------------------------------
	std::set< Size > caAtomsToConstrain;
	//if missing density exists fix gaps.
	protocols::simple_moves::MissingDensityToJumpMoverOP fixMissingDensityMover (new protocols::simple_moves::MissingDensityToJumpMover());
	fixMissingDensityMover->apply(pose);
	core::scoring::ScoreFunctionOP scorefxn_ = get_score_function();
	core::scoring::ScoreFunctionOP scorefxn_w_csts_ = get_score_function();
	scorefxn_w_csts_->set_weight( coordinate_constraint, COORDINATE_CST_WT );
	core::pose::PoseOP relaxed_poseOP = new core::pose::Pose(pose);
	Real dist_deviation_thresh = DIST_DEVIATION_THRESH_START;
	//single round of relax--------------------
	FastRelax relaxer( scorefxn_,RDS_OF_RELAX,"NO CST RAMPING");
	relaxer.apply(*relaxed_poseOP);
	delete_virtual_residues(*relaxed_poseOP);
	//Determine CA deviation that occurs in single round------------------
	Real gdtmm(core::scoring::CA_gdtmm(*relaxed_poseOP,pose));
	if(gdtmm<gdtThresh){
		std::set< Size > coreResiduesToConstrain = get_coreResiduesToConstrain(CORE_RESIDUES_TO_CONSTRAIN,pose);
		std::set< Size > coreDistDeviationResiduesToConstrain = get_coreDistDeviationResiduesToConstrain(dist_deviation_thresh,*relaxed_poseOP,pose);
		std::set_union(coreResiduesToConstrain.begin(), coreResiduesToConstrain.end(), coreDistDeviationResiduesToConstrain.begin(),  coreDistDeviationResiduesToConstrain.end(), std::inserter(caAtomsToConstrain, caAtomsToConstrain.begin()));
		ConstraintSetOP cst_set = convert_caAtomsToConstrain_to_coordCsts(caAtomsToConstrain,pose);
		core::pose::PoseOP rd2relaxed_poseOP = new core::pose::Pose(pose);
		add_virtual_residue_to_cterm(*rd2relaxed_poseOP);
		rd2relaxed_poseOP->constraint_set(cst_set);
		FastRelax relaxer2( scorefxn_w_csts_,RDS_OF_RELAX,"NO CST RAMPING");
		relaxer2.apply(*rd2relaxed_poseOP);
		delete_virtual_residues(*rd2relaxed_poseOP);
		gdtmm = (core::scoring::CA_gdtmm(*rd2relaxed_poseOP,pose));
		//repeat relaxing while gdt is below .99 and the dist_deviation_thresh has not fallen below 1.0
		while((gdtmm < gdtThresh) && (dist_deviation_thresh >= DIST_DEVIATION_THRESH_MIN)){
			tr.Debug << "gdtmm" << gdtmm << std::endl;
			std::set< Size > tmp_caAtomsToConstrain;
			std::set< Size > distDeviationResiduesToConstrain = get_distDeviationResiduesToConstrain(dist_deviation_thresh,0,*rd2relaxed_poseOP, pose);//on the second rounds and later gap size is always set to 0. That way this doesn't limits the number of loops
			tr.Debug << "empty ?" << distDeviationResiduesToConstrain.empty() << std::endl;
			if(!distDeviationResiduesToConstrain.empty()){
				std::set_union(caAtomsToConstrain.begin(), caAtomsToConstrain.end(), distDeviationResiduesToConstrain.begin(), distDeviationResiduesToConstrain.end(), std::inserter(tmp_caAtomsToConstrain, tmp_caAtomsToConstrain.begin()));
				caAtomsToConstrain = tmp_caAtomsToConstrain;
				cst_set = convert_caAtomsToConstrain_to_coordCsts(caAtomsToConstrain,pose);
				rd2relaxed_poseOP = new core::pose::Pose(pose);
				add_virtual_residue_to_cterm(*rd2relaxed_poseOP);
				rd2relaxed_poseOP->constraint_set(cst_set);
				relaxer2.apply(*rd2relaxed_poseOP);
				delete_virtual_residues(*rd2relaxed_poseOP);
				gdtmm = core::scoring::CA_gdtmm(*rd2relaxed_poseOP,pose);
			}
			else{
				dist_deviation_thresh = dist_deviation_thresh - DIST_DEVIATION_THRESH_STEP;
			}
		}
	}
	return(caAtomsToConstrain);
}

/// @brief Same as above but with defaults for gap size and gdt_thresh
std::set<Size> get_residuesToConstrain(Pose& pose){
		return get_residuesToConstrain(0,.99,pose);
}
/// @brief Outputs coordinate contraints in standard coordinate constraint format
void output_coordCsts(const std::set< Size > caAtomsToConstrain,std::ostream & out, Pose& pose){
	using namespace protocols;
	using namespace relax;
	if(caAtomsToConstrain.size() != 0){//don't output coordinate constraints if no atoms to constrain
		std::set< Size >::const_iterator caAtomsToConstrain_start,caAtomsToConstrain_stop;
		int virtualRes_id = pose.total_residue();
		numeric::xyzVector< core::Real > atomLocation(0.0,0.0,0.0);
		caAtomsToConstrain_start = caAtomsToConstrain.begin();
		while(caAtomsToConstrain_start != caAtomsToConstrain.end()){
			Size residue_id = *caAtomsToConstrain_start;
			atomLocation = pose.residue(residue_id).xyz("CA");
			out << "CoordinateConstraint  CA " << residue_id <<" N " << virtualRes_id  << " " << atomLocation.x() << " " << atomLocation.y() << " " << atomLocation.z() << " HARMONIC 0 1" << std::endl;
			caAtomsToConstrain_start++;
		}
		numeric::xyzVector< Real > virtualResLocation =  pose.residue(virtualRes_id).xyz("N");
		out << "CoordinateConstraint N " << virtualRes_id << " N " << virtualRes_id << " " << virtualResLocation.x() << " " << virtualResLocation.y() << " " << virtualResLocation.z() << " HARMONIC 0 1" << std::endl;
	}
}

// @brief Outputs coordinate contraints in standard coordinate constraint
// uses the partial thread and places the constraint in the center of mass
// of mass of the partial thread.
void output_coordCsts(const std::set< Size > caAtomsToConstrain,std::ostream & out, Pose& pose, const SequenceAlignment aln,const string query_sequence, bool only_res_out){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	std::set< Size >::const_iterator caAtomsToConstrain_start,caAtomsToConstrain_stop;
	if(caAtomsToConstrain.size() > 0){
		if(only_res_out){
			caAtomsToConstrain_start = caAtomsToConstrain.begin();
			while(caAtomsToConstrain_start != caAtomsToConstrain.end()){
				Size templateResidue_id = aln.sequence_mapping(2,1)[*caAtomsToConstrain_start]; //The remaped CA position
				if(templateResidue_id !=0){
					out << templateResidue_id << std::endl;
				}
				caAtomsToConstrain_start++;
			}
		}
		else{
			numeric::xyzVector< core::Real > atomLocation(0.0,0.0,0.0);
			//getting center of mass
			Real virtualRes_id = query_sequence.size()+1;
			atomLocation = get_centerOfMass(pose,aln,query_sequence);
			out << "CoordinateConstraint ORIG " << virtualRes_id << " ORIG " << virtualRes_id << " " <<  atomLocation.x() << " " << atomLocation.y() << " " << atomLocation.z() << " HARMONIC 0 1" << std::endl;
			caAtomsToConstrain_start = caAtomsToConstrain.begin();
			while(caAtomsToConstrain_start != caAtomsToConstrain.end()){
				Size templateResidue_id = aln.sequence_mapping(2,1)[*caAtomsToConstrain_start]; //The remaped CA position
				if(templateResidue_id !=0){
					Size residue_id =*caAtomsToConstrain_start;
					atomLocation = pose.residue(residue_id).xyz("CA");
					out << "CoordinateConstraint  CA " << templateResidue_id <<" ORIG " << virtualRes_id  << " " << atomLocation.x() << " " << atomLocation.y() << " " << atomLocation.z() << " HARMONIC 0 1" << std::endl;
				}
				caAtomsToConstrain_start++;
			}
		}
	}
}



/// @brief input coordinate constraints from file in standard coordinate
/// constraint format saves the data as a list of CA to constrain.
///For this particular implementation of reading the constraints in I only
/// care about the residue id's
std::set< Size > input_coordCsts(const string inputFileName){
	std::set< Size > caAtomsToConstrain;
	utility::io::izstream data( inputFileName.c_str() );
	if ( !data ) {
		utility_exit_with_message( "[ERROR] Unable to open constraints file: "+ inputFileName );
	}
	std::string line;
	while( getline( data, line ) ) {
		Size res1, res2;
		std::string cstType, name1, name2;
		std::string func_type;
		std::string type;
		numeric::xyzVector< Real > xyz_target(0.0,0.0,0.0);
		std::istringstream line_stream( line );
		line_stream >> cstType
								>> name1 >> res1
								>> name2 >> res2
								>> xyz_target.x()
								>> xyz_target.y()
								>> xyz_target.z()
								>> func_type;
		if (cstType == "CoordinateConstraint"){
			if(name1 == "CA")
				caAtomsToConstrain.insert(res1);
			if(name2 == "CA")
				caAtomsToConstrain.insert(res2);
		}
	}
	return caAtomsToConstrain;
}

/// @brief outputs Fasta with virtual atoms attached.
void output_fastaWVirtual(const std::string fastaFileName, std::ostream & out){
	using namespace core::sequence;
	vector1< SequenceOP > fastaSequenceOP = read_fasta_file(fastaFileName);
	out <<">" << fastaSequenceOP[1]->id() << std::endl;
	out << fastaSequenceOP[1]->sequence() << "X[VRT]" << std::endl;
}

/// @brief inputs the sequence alignments and maps them to to either pdbid or input file name(alignment file name).
map<string,SequenceAlignment> input_alignmentsMapped(bool mapToPdbid){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  map<string,SequenceAlignment> alns;
  map<string,SequenceAlignment>::iterator location_alns;
	if(option[ in::file::alignment ].user()){
		vector1< std::string > align_fns = option[ in::file::alignment ]();
		for ( Size ii = 1; ii <= align_fns.size(); ++ii ) {
			vector1< SequenceAlignment > tmp_alns = core::sequence::read_aln(
				  option[ cm::aln_format ](), align_fns[ii]
				  );
			for ( Size jj = 1; jj <= tmp_alns.size(); ++jj ) {
				string mapToName;
				if(mapToPdbid == true){
					string aln_id = tmp_alns[jj].sequence(2)->id();
					string pdbid = aln_id.substr(0,5);
					mapToName = pdbid;
				}
				else{
					std::stringstream numbConvert;
					string alnBaseName = utility::file_basename(option[ in::file::alignment ]()[ii]);
					numbConvert << alnBaseName << "_" << jj;
					string aln_id = numbConvert.str();
					mapToName = aln_id;
				}
				alns.insert(std::pair<string,SequenceAlignment>(mapToName,tmp_alns[jj]));
			}
		}
	}
  return(alns);
}

///@brief inputs the sequence alignments and keeps them ordered the way they
/// were in the alingment.filt file.
vector1<SequenceAlignment> input_alignments(){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::sequence;
  vector1< std::string > align_fns = option[ in::file::alignment ]();
  vector1< SequenceAlignment > alns = core::sequence::read_aln(
				    option[ cm::aln_format ](), align_fns[1]
				  );
  return(alns);
}

///@brief gets the first and last residues from an alignment ignoring residues that are gapped alnIdx can be 1 or 2 depending if you want the first/last from the target or template.
void get_terminal_aln_res(const SequenceAlignment aln,const Size alnIdx, Size & firstRes, Size & lastRes){
	firstRes = 0;
	lastRes = 0;
	//gets the first and last residues
	for(Size ii=1; ii<=aln.length(); ii++){
		if(!aln.is_gapped(ii)){
			lastRes = aln.sequence_indices(ii)[alnIdx];
			if (firstRes ==0)
				firstRes = aln.sequence_indices(ii)[alnIdx];
		}
	}
}




} //namespace cstEnergyBalance
}//namespace devel
