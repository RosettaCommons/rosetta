// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/energy_based_clustering/EnergyBasedClusteringProtocol.cc
/// @brief Performs the work done by the energy_based_clustering app.  Uses an energy-biased cookie-cutter approach to
/// cluster a large number of structures without generating an all-by-all RMSD matrix.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Associated headers
#include <protocols/energy_based_clustering/EnergyBasedClusteringProtocol.hh>

// Core headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/variant_util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/TorsionID.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/rms_util.hh>

//Protocols headers
#include <protocols/relax/FastRelax.hh>
#include <protocols/cyclic_peptide/CycpepSymmetryFilter.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>

//Utility headers
#include <utility/vector1.hh>

//Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Numeric headers
#include <numeric/xyzVector.fwd.hh>
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/model_quality/rms.hh>

static basic::Tracer TR( "protocols.cluster.energy_based_clustering.EnergyBasedClusteringProtocol" );

#define CNCa_ANGLE 121.7
#define CNH_ANGLE 119.15
#define CaCN_ANGLE 116.2
#define OCN_ANGLE 123.01

namespace protocols {
namespace energy_based_clustering {

/// @brief Default constructor -- initializes from global options system.
EnergyBasedClusteringProtocol::EnergyBasedClusteringProtocol():
	utility::pointer::ReferenceCount(),
	options_(true), //Read options from global options system
	constraints_file_contents_(""),
	n_clusters_from_last_run_(0)
{}

/// @brief Options constructor -- avoids access to global options system.
EnergyBasedClusteringProtocol::EnergyBasedClusteringProtocol(
	EnergyBasedClusteringOptions const & options
) :
	utility::pointer::ReferenceCount(),
	options_( options ), //Don't read options from global options system; just copy them from input.
	constraints_file_contents_(""),
	n_clusters_from_last_run_(0)
{}

EnergyBasedClusteringProtocol::~EnergyBasedClusteringProtocol()= default;

EnergyBasedClusteringProtocolOP
EnergyBasedClusteringProtocol::clone() const {
	return EnergyBasedClusteringProtocolOP( new EnergyBasedClusteringProtocol( *this ) );
}

/// @brief Indicate which options-system options are relvant.
void
EnergyBasedClusteringProtocol::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys::cluster::energy_based_clustering;

	option.add_relevant( prerelax );
	option.add_relevant( relax_rounds );
	option.add_relevant( cluster_by );
	option.add_relevant( use_CB );
	option.add_relevant( cluster_radius );
	option.add_relevant( residues_to_ignore );
	option.add_relevant( chains_to_ignore );
	option.add_relevant( limit_structures_per_cluster );
	option.add_relevant( limit_clusters );
	option.add_relevant( cyclic );
	option.add_relevant( cyclic_symmetry );
	option.add_relevant( cyclic_symmetry_mirroring );
	option.add_relevant( cyclic_symmetry_threshold );
	option.add_relevant( cluster_cyclic_permutations );
	option.add_relevant( mutate_to_ala );
	option.add_relevant( disulfide_positions );
	option.add_relevant( homooligomer_swap );
	option.add_relevant( silent_output );
	option.add_relevant( cst_file );
	option.add_relevant( extra_rms_atoms );
	option.add_relevant( rebuild_all_in_dihedral_mode );
	option.add_relevant( basic::options::OptionKeys::out::prefix );
	option.add_relevant( basic::options::OptionKeys::in::file::s );
	option.add_relevant( basic::options::OptionKeys::in::file::l );
	option.add_relevant( basic::options::OptionKeys::in::file::silent );
}

/// @brief Perform the clustering, based on options set in the options object.
void
EnergyBasedClusteringProtocol::go() {

	//Set up the scorefunction to use:
	core::scoring::ScoreFunctionOP sfxn( set_up_scorefunction() );

	//Check that options were sensibly set:
	do_option_checks();

	//Parse the user-specified list of additional atoms to use in the RMSD calculation:
	utility::vector1<core::id::NamedAtomID> extra_atom_list;
	parse_extra_atom_list(extra_atom_list); //Does nothing if no list provided.

	core::Size count = 0;
	core::Real lowestE = 0; //Lowest energy encountered so far.
	core::Size lowestE_index = 0; //The number of the pose with the lowest energy encountered.
	utility::vector1 <core::Real> poseenergies; //Vector of energies of the poses.
	utility::vector1 < utility::vector1 <core::Real> > posedata; //Vector of vectors to store the data that will be used for clustering.
	utility::vector1 < utility::vector1 <core::Real> > dihedral_reconstruction_data; //Vector of vectors to store atom positions in dihedral mode.  Only populated if rebuild_all_in_dihedral_mode is true.
	utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > alignmentdata; //Vector of arrays of x,y,z coordinates of atoms to be used for alignment in Cartesian clustering.
	utility::vector1 < core::Size > cluster_assignments; //List of which cluster each structure is assigned to.
	utility::vector1 < core::Size > cluster_offsets; //List of offsets for cyclic permutations when clustering (measured in number of amino acid positions we've offset by).
	utility::vector1 < core::Size > cluster_oligomer_permutations; //List of permutations of oligomers if the user has specified the option to do that.

	TR << "Scoring energies of all input structures." << std::endl;

	core::pose::Pose firstpose;

	//Symmetry filter (used only if the cyclic_symmetry flag is used):
	protocols::cyclic_peptide::CycpepSymmetryFilter symmfilter;
	if ( options_.cyclic_symmetry_ > 1 ) {
		symmfilter.set_symm_repeats( options_.cyclic_symmetry_ );
		symmfilter.set_mirror_symm( options_.cyclic_symmetry_mirroring_ );
		symmfilter.set_angle_threshold( options_.cyclic_symmetry_threshold_ );
	}

	//Import and score all structures:
	do_initial_import_and_scoring( count, lowestE, lowestE_index, firstpose, symmfilter, poseenergies, posedata, alignmentdata, dihedral_reconstruction_data, cluster_assignments, cluster_offsets, cluster_oligomer_permutations, sfxn, extra_atom_list );

	TR << "Clustering, starting with lowest-energy structure as the center of the first cluster." << std::endl;

	core::Size unclustered_count( count ); //The number of structures that remain to be clustered
	core::Size cluster_count( 0 ); //The number of clusters created
	utility::vector1 < utility::vector1 <core::Real> > clustcenter_list; //Vector of vectors to store the data for the cluster centres
	utility::vector1 < utility::vector1 <core::Size> > clusterlist_sortedbyenergy; //A vector storing lists of the states assigned to each cluster, sorted by energy.

	while ( unclustered_count > 0 ) {
		//Make a new cluster
		cluster_count++;

		//The cumulative weighting of the points considered so far
		//core::Real weighting_accumulator = option[v_weightbyenergy]() ? exp(-lowestE/option[v_kbt]()) : 1.0;
		//core::Real currentweighting = 0.0;

		//Assign the lowest energy structure to the next cluster
		cluster_assignments[lowestE_index] = cluster_count;
		utility::vector1 <core::Size> cluster_sortedbyenergy;
		cluster_sortedbyenergy.push_back(lowestE_index);
		unclustered_count--;

		TR << "Started cluster " << cluster_count << " and added structure " << lowestE_index << " to it." << std::endl;
		utility::vector1 <core::Real> clustcenter = posedata[lowestE_index]; //Set the center of the current cluster

		//Make a list of unassigned candidate structures:
		utility::vector1 <core::Size> candidatelist;
		candidatelist.reserve( unclustered_count );
		for ( core::Size istruct=1; istruct<=count; istruct++ ) {
			if ( cluster_assignments[istruct]==0 ) candidatelist.push_back(istruct);
		}
		TR << "\tMade list of " << candidatelist.size() << " unassigned structures." << std::endl;

		//Assign members of the candidate list to the current cluster if they fall within R_cluster of the cluster center.
		core::Real currentdist(0.0);
		core::Size currentcyclicoffset(0); //Only used for calculating cyclic permutations
		core::Size current_oligomer_permutation(0); //Only used for calculating permutations when swapping around homodimers.
		for ( core::Size istruct(1), istructmax(candidatelist.size()); istruct<=istructmax; ++istruct ) {
			currentdist=calc_dist(clustcenter, posedata[candidatelist[istruct]], options_.cluster_by_, alignmentdata[lowestE_index], alignmentdata[candidatelist[istruct]], firstpose.size(), firstpose, currentcyclicoffset, current_oligomer_permutation);

			if ( currentdist < options_.cluster_radius_ ) {
				std::streamsize const old_precision(TR.precision());
				TR.precision(6);
				TR << "\tAdding structure " << candidatelist[istruct] << " (" << currentdist << " from the cluster centre)." << std::endl;
				TR.precision(old_precision);
				cluster_assignments[candidatelist[istruct]]=cluster_count;
				cluster_sortedbyenergy.push_back(candidatelist[istruct]);
				if ( options_.cluster_cyclic_permutations_ ) cluster_offsets[candidatelist[istruct]]=currentcyclicoffset;
				if ( options_.homooligomer_swap_ ) cluster_oligomer_permutations[candidatelist[istruct]]=current_oligomer_permutation;
			}
		} //for loop through all candidates

		sort_cluster_list(cluster_sortedbyenergy, poseenergies); //Sort the list of states in this cluster by energy
		clusterlist_sortedbyenergy.push_back(cluster_sortedbyenergy); //Add the list of states in this cluster to the list of states in each cluster

		//The following is commented-out code related to PCA analysis of clusters:
		/*utility::vector1 < utility::vector1 < core::Real > > pca_vector_list; //A place to store the PCA vectors
		utility::vector1 < core::Real > coeff_list; //A place to store coefficients (amplitudes) for each PCA vector
		if(!option[v_skip_PCA]()) {
		shift_center_and_PCA(clustcenter, pca_vector_list, coeff_list, lowestE_index, posedata, poseenergies,
		cluster_sortedbyenergy, cluster_count, options_.cluster_by_, firstpose, cluster_offsets,
		cluster_oligomer_permutations, extra_atom_list);
		output_PCA(pca_vector_list, coeff_list, cluster_count);
		clustcenter_list.push_back(clustcenter); //Store the current cluster center
		}*/

		bool no_unassigned(true);
		for ( core::Size i=1; i<=poseenergies.size(); i++ ) {
			if ( cluster_assignments[i]!=0 ) continue; //Continue if this structure has been assigned
			if ( no_unassigned || poseenergies[i]<lowestE ) { //If this is the first unassigned encountered OR the lowest energy unassigned encountered so far
				no_unassigned=false;
				lowestE=poseenergies[i];
				lowestE_index=i;
			}
		}

		if ( no_unassigned ) break; //If this is still true, there were no unassigned structures.
	}


	//Outputs:
	TR << "Cluster\tStructure\tFile_out" << std::endl;

	for ( core::Size i=1, imax(clusterlist_sortedbyenergy.size()); i<=imax; i++ ) {
		if ( options_.limit_clusters_ != 0 && i > static_cast<core::Size>( options_.limit_clusters_ ) ) {
			TR << "Maximum number of clusters for output reached." << std::endl;
			break;
		}
		core::pose::Pose pose1;
		for ( core::Size j=1, jmax(clusterlist_sortedbyenergy[i].size()); j<=jmax; j++ ) {
			core::pose::Pose temppose;
			if ( options_.limit_structures_per_cluster_==0 || j<=static_cast<core::Size>(options_.limit_structures_per_cluster_) ) {
				pose_from_posedata (firstpose, temppose, options_.cluster_by_, posedata[ clusterlist_sortedbyenergy[i][j] ], dihedral_reconstruction_data[clusterlist_sortedbyenergy[i][j]], options_.rebuild_all_in_dihedral_mode_);
				if ( j>1 && options_.homooligomer_swap_ ) swap_chains( temppose, cluster_oligomer_permutations[clusterlist_sortedbyenergy[i][j]] ); //Swap chains around.
				if ( j==1 ) pose1=temppose;
				else align_with_offset(temppose, pose1, (options_.cluster_cyclic_permutations_ ? cluster_offsets[ clusterlist_sortedbyenergy[i][j] ] : 0), extra_atom_list );
			}
			std::stringstream outfile1, outfile2;
			char curstructtag[128];
			std::string const prepended( options_.output_prefix_.empty() ? "" : options_.output_prefix_ + "_" );
			sprintf(curstructtag, "%sc.%lu.%lu", prepended.c_str(), i, j);
			if ( options_.silent_output_ ) {
				core::io::silent::SilentFileOptions opts;
				opts.read_from_global_options();
				core::io::silent::SilentFileData outsilentfiledata(opts);
				core::io::silent::SilentStructOP outsilentstruct ( core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts) );
				outsilentstruct->fill_struct(temppose, curstructtag);
				if ( j==1 ) {
					outfile1.clear();
					if ( !options_.output_prefix_.empty() ) {
						outfile1 << options_.output_prefix_ << "_";
					}
					outfile1 << "clusters_firstmember.out";;
					outsilentfiledata.write_silent_struct((*outsilentstruct), outfile1.str());
				}
				if ( options_.limit_structures_per_cluster_==0 || j<=options_.limit_structures_per_cluster_ ) {
					if ( options_.limit_structures_per_cluster_==0 ) outfile2 << "clusters_all_members.out";
					else {
						outfile2.clear();
						if ( !options_.output_prefix_.empty() ) {
							outfile2 << options_.output_prefix_ << "_";
						}
						outfile2 << "clusters_first_" << options_.limit_structures_per_cluster_ << "_members.out" << std::endl;
					}
					outsilentfiledata.write_silent_struct((*outsilentstruct), outfile2.str());
				}
			} else {
				outfile1.clear();
				outfile1 << curstructtag << ".pdb";
				if ( options_.limit_structures_per_cluster_==0 || j<=options_.limit_structures_per_cluster_ ) temppose.dump_pdb(outfile1.str());
			}
			TR << i << "\t" << clusterlist_sortedbyenergy[i][j] << "\t" << curstructtag << std::endl;
		}
	}

	/*printf("\nWriting cluster centers.\n"); fflush(stdout);
	for(core::Size i=1; i<=clustcenter_list.size(); i++) {
	core::pose::Pose cenpose;
	pose_from_posedata(firstpose, cenpose, clustermode, clustcenter_list[i]);
	char outfile[64];
	sprintf(outfile, "center_%lu.pdb", i);
	cenpose.dump_pdb(outfile);
	}*/

	n_clusters_from_last_run_ = clusterlist_sortedbyenergy.size();

	TR << "*****JOB COMPLETED*****" << std::endl;
	TR.flush();
}


/////////////////////// PRIVATE FUNCTIONS ///////////////////////

/// @brief Function to determine whether a value is in a list
bool
EnergyBasedClusteringProtocol::is_in_list (
	core::Size const val,
	utility::vector1 < core::Size > const &vallist
) const {
	if ( vallist.size()>0 ) {
		for ( core::Size i=1, listlength=vallist.size(); i<=listlength; ++i ) {
			if ( vallist[i]==val ) return true;
		}
	}
	return false;
}

/// @brief Align one pose to another with an offset in the residue count.
void
EnergyBasedClusteringProtocol::align_with_offset (
	core::pose::Pose &pose1,
	core::pose::Pose const &pose2, //the target pose -- doesn't change.
	core::Size const offset,
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) const {
	using namespace core::id;

	AtomID_Map< AtomID > amap;
	signed int ir;
	core::pose::initialize_atomid_map(amap,pose1, AtomID::BOGUS_ATOM_ID() );
	for ( core::Size ir2(1), nres(pose2.total_residue()); ir2 <= nres; ++ir2 ) {
		ir= ir2-offset;
		if ( ir<=0 ) ir+=pose1.total_residue();
		for ( core::Size ia(1), iamax(pose2.residue_type(ir2).nheavyatoms()); ia <= iamax; ++ia ) {
			if ( use_in_rmsd_offset(pose1,pose2,ir2,ia, offset, extra_atom_list) ) {
				amap[AtomID(ia,ir)] = AtomID(ia,ir2);
			}
		}
	}
	core::scoring::superimpose_pose( pose1, pose2, amap );
}

/// @brief Reconstruct a pose from posedata.
/// @details The reconstruction_data vector is unused for anything except dihedral-based building whe rebuild_all_in_dihedral_mode is true.  In that case,
/// this must be a full vector of coordinates of atoms.  Otherwise, it can be empty.
void
EnergyBasedClusteringProtocol::pose_from_posedata (
	core::pose::Pose const &inputpose,
	core::pose::Pose &outputpose,
	EBC_ClusterType const clustermode,
	utility::vector1<core::Real> const &posedata,
	utility::vector1<core::Real> const &reconstruction_data,
	bool const rebuild_all_in_dihedral_mode/*=false*/
) const {
	outputpose=inputpose; //Copy the input pose for initial pose setup.
	make_disulfides(outputpose);

	core::Size counter(0);
	numeric::xyzVector <core::Real> atomxyz;
	switch(clustermode) {
	case EBC_bb_cartesian :
		for ( core::Size ir(1), irmax(outputpose.size()); ir<=irmax; ++ir ) {
			for ( core::Size ia(1), iamax(outputpose.residue(ir).natoms()); ia<=iamax; ++ia ) {
				counter+=3;
				atomxyz.assign( posedata[counter-2], posedata[counter-1], posedata[counter] );
				outputpose.set_xyz(core::id::AtomID(ia, ir), atomxyz);
			}
		}
		break; //case
	case EBC_bb_dihedral : //Dihedral clustering case.
		if ( rebuild_all_in_dihedral_mode ) {
			for ( core::Size ir(1), irmax(outputpose.size()); ir<=irmax; ++ir ) {
				for ( core::Size ia(1), iamax(outputpose.residue(ir).natoms()); ia<=iamax; ++ia ) {
					counter+=3;
					atomxyz.assign( reconstruction_data[counter-2], reconstruction_data[counter-1], reconstruction_data[counter] );
					outputpose.set_xyz(core::id::AtomID(ia, ir), atomxyz);
				}
			}
		} else {
			for ( core::Size i(1), nres(outputpose.total_residue()); i<=nres; ++i ) { //Loop through all residues
				//Can't set DoFs of ligands, metals, or ignored residues or chains, unfortunately.
				if ( outputpose.residue_type(i).is_ligand() || outputpose.residue_type(i).is_metal() ) continue;
				if ( options_.chains_to_ignore_.size() != 0 && is_in_list(outputpose.chain(i), options_.chains_to_ignore_) ) continue;
				if ( options_.residues_to_ignore_.size() != 0 && is_in_list(i, options_.residues_to_ignore_) ) continue;

				for ( core::Size j(1), jmax(outputpose.residue(i).mainchain_torsions().size()); j<=jmax; ++j ) { //Loop through torsions of current residue
					if ( j == 1 && (outputpose.residue_type(i).is_lower_terminus() || outputpose.residue(i).connected_residue_at_lower() == 0) ) continue;
					if ( (j == jmax || j == jmax - 1) && ( outputpose.residue_type(i).is_upper_terminus() || outputpose.residue(i).connected_residue_at_upper() == 0 ) ) continue;
					++counter;
					outputpose.set_torsion( core::id::TorsionID( i, core::id::BB, j ), posedata[counter] );
				}
			}
		}
		break; //case
	}

	make_disulfides(outputpose);
	outputpose.update_residue_neighbors();
}

/// @brief Sort the list of states in a cluster by energies.
/// @details Inefficient selection sort used.  That's okay -- this sorts small lists.
void
EnergyBasedClusteringProtocol::sort_cluster_list(
	utility::vector1 <core::Size> &statelist,
	utility::vector1 <core::Real> &poseenergies
) const {
	//This won't be the world's most efficient sort.  I'll use a selection sort algorithm:
	core::Real lowestE = 0.0;
	core::Size lowestEentry = 0;
	core::Size buffer = 0;

	for ( core::Size i(1), imax(statelist.size()); i<imax; ++i ) {
		for ( core::Size j(i), jmax(statelist.size()); j<=jmax; ++j ) {
			if ( j==i || poseenergies[statelist[j]]<lowestE ) {
				lowestE = poseenergies[statelist[j]];
				lowestEentry=j;
			}
		}
		if ( lowestEentry!=i ) {
			buffer=statelist[i];
			statelist[i]=statelist[lowestEentry];
			statelist[lowestEentry]=buffer;
		}
	}
}

/// @brief Function to calculate the RMSD between two poses, based on whatever is in "posedata" for the second only.
/// This assumes that a pose already exists for the first, and that this is provided in refpose.
core::Real
EnergyBasedClusteringProtocol::calc_dist(
	utility::vector1 <core::Real> const &vect1, //Backbone dihedral vector 1
	utility::vector1 <core::Real> const &vect2, //Backbone dihedral vector 2
	EBC_ClusterType const clustmode,
	utility::vector1 < numeric::xyzVector < core::Real > > const &alignmentvect1, //Alignment vector 1 (for Cartesian clustering)
	utility::vector1 < numeric::xyzVector < core::Real > > const &alignmentvect2, //Alignment vector 2 (for Cartesian clustering)
	core::Size const nresidues,
	core::pose::Pose const &firstpose //Used for reference only!
) const {
	static const std::string errmsg( "Error in protocols::cluster::energy_based_clustering::EnergyBasedClusteringProtocol::calc_dist(): " );

	core::Real accumulator( 0.0 );

	runtime_assert_string_msg( alignmentvect1.size() == alignmentvect2.size(), errmsg + "Alignment vector size mismatch." );

	if ( clustmode==EBC_bb_cartesian ) {
		accumulator=numeric::model_quality::calc_rms(alignmentvect1, alignmentvect2); //Calculate the rms, if we're doing Cartesian clustering
	} else if ( clustmode==EBC_bb_dihedral ) { //Backbone dihedral clustering
		for ( core::Size ir(1),index(1); ir<=nresidues; ir++ ) {
			if ( firstpose.residue_type(ir).is_ligand() || firstpose.residue_type(ir).is_metal() ) continue; //Skip ligands and metals.

			//The number of torsion angles in the current residue:
			core::Size curres_torsioncount( firstpose.residue(ir).mainchain_torsions().size() );
			if ( firstpose.residue(ir).is_lower_terminus() || firstpose.residue(ir).connected_residue_at_lower() == 0 ) curres_torsioncount--;
			if ( firstpose.residue(ir).is_upper_terminus() || firstpose.residue(ir).connected_residue_at_upper() == 0 ) curres_torsioncount-=2;

			if ( options_.chains_to_ignore_.size() > 0 ) {
				bool continue_on(false);
				for ( core::Size j(1), jmax(options_.chains_to_ignore_.size()); j<=jmax; ++j ) {
					if ( firstpose.residue(ir).chain() == static_cast<core::Size>(options_.chains_to_ignore_[j]) ) {
						index+=curres_torsioncount; //Increment the dihedral angle index to skip this residue.
						continue_on=true;
						break; //Go on to the next residue.
					}
				}
				if ( continue_on ) continue; //Go on to the next residue.
			}

			if ( options_.residues_to_ignore_.size() > 0 ) {
				bool continue_on(false);
				for ( core::Size j(1), jmax(options_.residues_to_ignore_.size()); j<=jmax; ++j ) {
					if ( ir == static_cast<core::Size>(options_.residues_to_ignore_[j]) ) {
						index+=curres_torsioncount; //Increment the dihedral angle index to skip this residue.
						continue_on=true;
						break; //Go on to the next residue.
					}
				}
				if ( continue_on ) continue;
			}

			//If we get here, then the current residue is to be included in the calculation.
			for ( core::Size i=1; i<=curres_torsioncount; i++ ) {
				accumulator+=pow( numeric::principal_angle_degrees( vect1[index] - vect2[index] ), 2.0 );
				index++;
			}
		}
		accumulator=sqrt(accumulator);
	}

	return accumulator;
}

/// @brief Overloaded form of calc_dist for cyclic permutations and for homooligomer swapping.
core::Real
EnergyBasedClusteringProtocol::calc_dist(
	utility::vector1 <core::Real> const &vect1,
	utility::vector1 <core::Real> const &vect2,
	EBC_ClusterType const clustmode,
	utility::vector1< numeric::xyzVector < core::Real > > const &alignmentvect1,
	utility::vector1< numeric::xyzVector < core::Real > > const &alignmentvect2,
	core::Size const nresidues,
	core::pose::Pose const &firstpose, //Used for reference only
	core::Size &offset, //Used only for cyclic permutations (OUTPUT)
	core::Size &permutation //Used only for homooligomer permutations (OUTPUT)
) const {

	static std::string const errmsg( "Error in protocols::cluster::energy_based_clustering::EnergyBasedClusteringProtocol::calc_dist(): " );

	core::Real dist(0.0);

	if ( options_.cluster_cyclic_permutations_ ) {
		if ( clustmode==EBC_bb_cartesian ) { //Clustering by backbone atom Cartesian coordinates
			utility::vector1< numeric::xyzVector < core::Real > > avect2 (alignmentvect2.size());
			for ( core::Size ioffset(0); ioffset<nresidues; ioffset+=options_.cyclic_permutation_offset_ ) { //Loop through all possible offsets
				core::Size offset_index_1(1); //The index of the first atom (i.e. alignment atom 1 of the OLD first residue), when offset by ioffset residues.
				//Figure out starting value of offset_index_1.  This means that we're putting the atoms from the END of the sequence ahead of the first atom:
				if ( ioffset>0 ) {
					for ( core::Size ir=nresidues, irmin(nresidues-ioffset); ir>irmin; --ir ) {
						offset_index_1 += alignment_atoms_in_res(firstpose, ir);
					}
				}

				for ( core::Size ir(1), index(1); ir<=nresidues; ++ir ) { //For each offset, loop through all residues
					//The number of alignment atoms stored for this residue:
					core::Size atoms_in_this_res(alignment_atoms_in_res(firstpose,ir));

					//Offset the alignment vectors:
					for ( core::Size ia(1); ia<=atoms_in_this_res; ++ia ) {
						avect2[offset_index_1] = alignmentvect2[index];
						++index; //Move to next atom.
						++offset_index_1; //Move this with index.
					}
					if ( ir==(nresidues-ioffset) ) offset_index_1=1; //Reset this if we've reached the wrap-around point.
				} //Looping through all residues

				core::Real const curdist = calc_dist(vect1, vect2, clustmode, alignmentvect1, avect2, nresidues, firstpose);

				if ( ioffset==0 || curdist<dist ) { //If this offset has yielded the lowest RMS encountered so far.
					dist = curdist;
					offset = ioffset;
				}
			}
		} else if ( clustmode==EBC_bb_dihedral ) { //Clustering by backbone dihedrals
			//Make copies of vect1 and vect2
			utility::vector1 <core::Real> v2 = vect2;
			for ( core::Size ioffset(0); ioffset<nresidues; ioffset+=options_.cyclic_permutation_offset_ ) {
				core::Size offset_index_1(2);
				if ( ioffset>0 ) {
					for ( core::Size ir=nresidues, irmin(nresidues-ioffset); ir>irmin; --ir ) {
						offset_index_1 += alignment_torsions_in_res(firstpose, ir);
					}
				}

				for ( core::Size ir(1), index(1); ir<=nresidues; ++ir ) { //Loop through all residues
					core::Size const torsions_in_this_res( alignment_torsions_in_res(firstpose, ir) );
					for ( core::Size itors=1; itors<=torsions_in_this_res; ++itors ) {
						v2[offset_index_1] = vect2[index];
						++offset_index_1;
						++index;
						if ( offset_index_1 > v2.size() ) offset_index_1=1; //Reset this at the wrap-around point.
					}
				}

				core::Real const curdist( calc_dist(vect1, v2, clustmode, alignmentvect1, alignmentvect2, nresidues, firstpose) );

				if ( ioffset==0 || curdist<dist ) { //If this offset has yielded the lowest distance encountered so far.
					dist = curdist;
					offset = ioffset;
				}
			}
		}
	} else if ( options_.homooligomer_swap_ ) { //if(options_.cluster_cyclic_permutations_) //consider all permutations of chains
		core::Size numchains(firstpose.num_chains());
		for ( core::Size i(1), imax( numeric::factorial(numchains)); i<=imax; i++ ) {
			if ( clustmode==EBC_bb_cartesian ) {
				utility::vector1< numeric::xyzVector < core::Real > > alignmentvect2_swapped;
				swap_alignment_vector(alignmentvect2, alignmentvect2_swapped, firstpose, i);
				core::Real curdist ( calc_dist(vect1, vect2, clustmode, alignmentvect1, alignmentvect2_swapped, nresidues, firstpose) );
				if ( i==1 || curdist < dist ) {
					dist=curdist;
					permutation=i;
				}
			} else if ( clustmode==EBC_bb_dihedral ) {
				//Currently doesn't work
				utility_exit_with_message( errmsg + "Homooligomer swapping currently does not work with dihedral-based clustering." );
				/*
				utility::vector1 < core::Real > vect2_swapped;
				swapvector(vect2, vect2_swapped, firstpose, i); //TODO -- CURRENTLY DOESN'T WORK!
				core::Real curdist = calc_dist(vect1, vect2_swapped, clustmode, alignmentvect1, alignmentvect2, nresidues, firstpose);
				if(i==1 || curdist < dist) {
				dist=curdist;
				permutation=i;
				}
				*/
			}
		}

	} else { //else if ( options_.v_homooligomer_swap_ ) {
		dist = calc_dist(vect1, vect2, clustmode, alignmentvect1, alignmentvect2, nresidues, firstpose);
	}

	return dist;
}

/// @brief Swap chains in an alignment vector (vector of x,y,z atom coordinates).
/// @details Can only be called if options_.cluster_by_ == EBC_bb_cartesian.
void
EnergyBasedClusteringProtocol::swap_alignment_vector (
	utility::vector1 < numeric::xyzVector <core::Real > > const &parentvect,
	utility::vector1 < numeric::xyzVector <core::Real > >  &swappedvect,
	core::pose::Pose const &refpose,
	core::Size const permutation_number
) const {
	debug_assert( options_.cluster_by_ == EBC_bb_cartesian ); //Should be true.

	swappedvect = parentvect; //Copy the parent vector

	core::Size const nchain(refpose.num_chains());

	utility::vector1 < core::Size > chainlist(nchain); //List of chains, which will be rearranged in each permutation's order
	utility::vector1 < core::Size > startaa_list(nchain), endaa_list(nchain); //List of starting amino acids of chains
	utility::vector1 < core::Size > startindex_list(nchain), endindex_list(nchain); //List of the starting index in the parentvect for each chain

	//Populate and permute the chain list:
	for ( core::Size i(1); i<=nchain; i++ ) chainlist[i] = i;
	if ( permutation_number>1 ) {
		for ( core::Size i(2); i<=permutation_number; ++i ) {
			std::next_permutation(chainlist.begin(), chainlist.end());
		}
	}

	//Store the starting and ending amino acid numbers for each chain, in the perturbed chain order:
	for ( core::Size i=1; i<=nchain; i++ ) get_start_and_end_aa(i, refpose, startaa_list[i], endaa_list[i]);

	//Store the starting indices for each chain, in the perturbed chain order:
	for ( core::Size i=1; i<=nchain; i++ ) get_start_and_end_indices(i, refpose, startindex_list[i], endindex_list[i]);

	core::Size counter(1);
	for ( core::Size ichain(1); ichain<=nchain; ++ichain ) {
		core::Size counter2( startindex_list[ichain] ); //Starting index in the parent vector for this chain
		for ( core::Size ir(startaa_list[ichain]); ir<=endaa_list[ichain]; ++ir ) { //Loop through all residues in the current chain
			for ( core::Size iatom(1), natom(alignment_atoms_in_res(refpose, ir)); iatom<=natom; iatom++ ) { //Loop through all alignment atoms in the current residue
				swappedvect[counter] = parentvect[counter2];
				counter++; counter2++; //Increment both counters to point at the next atom.
			}
		}
	}
}

/// @brief Returns the number of alignment atoms in a given residue.
/// @details Some special-case logic here, since oxygens are counted in peptide-like building-blocks.
core::Size
EnergyBasedClusteringProtocol::alignment_atoms_in_res (
	core::pose::Pose const &pose,
	core::Size const position
) const {

	core::chemical::ResidueType const &type( pose.residue_type(position) );
	core::conformation::Residue const &res( pose.residue(position) );

	//Cases where there are no atoms in the residue in question:
	//Ligands and metals:
	if ( type.is_ligand() || type.is_metal() ) return 0;
	//Ignored chains:
	if ( options_.chains_to_ignore_.size() > 0 ) {
		if ( is_in_list( pose.chain(position), options_.chains_to_ignore_ ) ) return 0;
	}
	//Ignored residues:
	if ( options_.residues_to_ignore_.size() > 0 ) {
		if ( is_in_list(position, options_.residues_to_ignore_) ) return 0;
	}

	core::Size numatoms( type.mainchain_atoms().size() );
	if ( (type.is_alpha_aa() || type.is_beta_aa() || type.is_gamma_aa() || type.is_oligourea() || type.is_peptoid()) &&
			!(type.is_upper_terminus() || res.connected_residue_at_upper() == 0)
			) {
		++numatoms; //O used for clustering
	}
	if ( options_.use_CB_ && res.has( res.is_peptoid() ? "CA1" : "CB" ) ) ++numatoms; //Add one atom if "CB" is to be counted.
	return numatoms;
}

/// @brief Returns the number of alignment torsion angles stored for a given residue.
core::Size
EnergyBasedClusteringProtocol::alignment_torsions_in_res(
	core::pose::Pose const &pose,
	core::Size const position
) const {
	core::chemical::ResidueType const &type( pose.residue_type(position) );
	core::conformation::Residue const &res( pose.residue(position) );

	//Cases where there are no atoms in the residue in question:
	//Ligands and metals:
	if ( type.is_ligand() || type.is_metal() ) return 0;
	//Ignored chains:
	if ( options_.chains_to_ignore_.size() > 0 ) {
		if ( is_in_list( pose.chain(position), options_.chains_to_ignore_ ) ) return 0;
	}
	//Ignored residues:
	if ( options_.residues_to_ignore_.size() > 0 ) {
		if ( is_in_list(position, options_.residues_to_ignore_) ) return 0;
	}

	core::Size numtors( res.mainchain_torsions().size() );
	if ( type.is_lower_terminus() || res.connected_residue_at_lower() == 0 ) --numtors; //No phi for first position.
	if ( type.is_upper_terminus() || res.connected_residue_at_upper() == 0 ) numtors -= 2; //No psi or omega for last position.

	return numtors;
}




/// @brief Function containing the special-case logic for alpha, beta, and gamma-amino acids.
bool
EnergyBasedClusteringProtocol::use_this_atom(
	core::conformation::Residue const &res1,
	core::conformation::Residue const &res2,
	core::Size const atomno,
	std::string const &atname1,
	std::string const &atname2/*=""*/
) const {

	core::chemical::ResidueType const &type1( res1.type() );
	core::chemical::ResidueType const &type2( res2.type() );

	if ( atname2 != "" && atname1 != atname2 ) return false;
	if ( !type1.has(atname1) ) return false;

	//Special-case logic for peptoids and alpha-, beta-, and gamma-amino acids:
	if ( (type1.is_alpha_aa() || type1.is_peptoid()) && (type2.is_alpha_aa() || type2.is_peptoid()) ) {
		if ( atname1 == "N" ||
				atname1 == "CA" ||
				atname1 == "C" ||
				( atname1 == "O" && res1.has_upper_connect() && res1.connected_residue_at_upper()!=0 &&
				res2.has_upper_connect() && res2.connected_residue_at_upper()!=0
				) ||
				( options_.use_CB_ && atname1 == "CB" ) ||
				atname1 == "CA1"
				) {
			return true;
		}
	} else if ( type1.is_beta_aa() && type2.is_beta_aa() ) {
		if ( atname1 == "N" ||
				atname1 == "CA" ||
				atname1 == "CM" ||
				atname1 == "C" ||
				( atname1 == "O" && res1.has_upper_connect() && res1.connected_residue_at_upper()!=0 &&
				res2.has_upper_connect() && res2.connected_residue_at_upper()!=0
				) ||
				( options_.use_CB_ && atname1 == "CB" )
				) {
			return true;
		}
	} else if ( type1.is_oligourea() && type2.is_oligourea() ) {
		if ( atname1 == "N" ||
				atname1 == "CA" ||
				atname1 == "CM" ||
				atname1 == "NU" ||
				atname1 == "C" ||
				( atname1 == "O" && res1.has_upper_connect() && res1.connected_residue_at_upper()!=0 &&
				res2.has_upper_connect() && res2.connected_residue_at_upper()!=0
				) ||
				( options_.use_CB_ && atname1 == "CB" )
				) {
			return true;
		}
	} else if ( type1.is_gamma_aa() && type2.is_gamma_aa() ) {
		if ( atname1 == "N" ||
				atname1 == "C4" ||
				atname1 == "C3" ||
				atname1 == "C2" ||
				atname1 == "C" ||
				( atname1 == "O" && res1.has_upper_connect() && res1.connected_residue_at_upper()!=0 &&
				res2.has_upper_connect() && res2.connected_residue_at_upper()!=0
				) ||
				( options_.use_CB_ && atname1 == "CB" )
				) {
			return true;
		}
	} else {
		//General polymer case: use mainchain atoms.
		if ( is_in_list( atomno, type1.mainchain_atoms() ) && is_in_list( atomno, type2.mainchain_atoms() ) ) {
			return true;
		}
	}
	return false;
}

/// @brief Given a residue number, and atom number, and a pose, determine whether that atom gets used in the RMSD calculation.
bool
EnergyBasedClusteringProtocol::use_in_rmsd(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size resno,
	core::Size atomno,
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) const {
	//If there are chains to ignore:
	if ( options_.chains_to_ignore_.size()>0 ) {
		for ( core::Size i(1), imax(options_.chains_to_ignore_.size()); i<=imax; ++i ) {
			if ( pose1.chain(resno) == options_.chains_to_ignore_[i] || pose2.chain(resno) == options_.chains_to_ignore_[i] ) return false; //If this chain is to be ignored, return false.
		}
	}

	//If there are residues to ignore:
	if ( options_.residues_to_ignore_.size()>0 ) {
		for ( core::Size i(1), imax(options_.residues_to_ignore_.size()); i<=imax; ++i ) {
			if ( resno == options_.residues_to_ignore_[i] ) return false; //If this residue is to be ignored, return false.
		}
	}

	//The name of the atom:
	std::string atname1( utility::strip( pose1.residue(resno).atom_name(atomno), " " ) );
	if ( atname1 == "CN" && pose1.residue_type(resno).is_alpha_aa() && pose1.residue_type(resno).is_n_methylated() ) atname1 = "CA1"; //Allow N-methyl amino acids to be aligned to peptoids.
	if ( !pose2.residue_type(resno).has(atname1) ) {
		if ( atname1 == "CA1" && pose2.residue_type(resno).is_alpha_aa() && pose2.residue_type(resno).is_n_methylated() && pose2.residue_type(resno).has("CN") ) return true;
		return false;
	}

	if ( use_this_atom( pose1.residue(resno), pose2.residue(resno), atomno, atname1 ) ) return true;

	//Check extra atoms:
	for ( core::Size i(1), imax(extra_atom_list.size()); i<=imax; ++i ) {
		if ( resno!=extra_atom_list[i].rsd() ) continue;
		if ( pose1.residue(resno).has( extra_atom_list[i].atom() ) && pose2.residue(resno).has( extra_atom_list[i].atom() ) && pose1.residue(resno).atom_index( extra_atom_list[i].atom() )==atomno ) return true;
	}

	return false;
}

/// @brief Given two poses, a residue number in the second pose, and atom number, and and offset for circular permutation, determine whether the atom is to be used in RMSD calculation.
bool
EnergyBasedClusteringProtocol::use_in_rmsd_offset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size resno2,
	core::Size atomno,
	core::Size const pose1_offset,
	utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) const {
	signed int resno(resno2-pose1_offset);
	if ( resno<=0 ) resno+=pose1.total_residue();

	if ( options_.chains_to_ignore_.size()>0 ) {
		for ( core::Size i(1), imax(options_.chains_to_ignore_.size()); i<=imax; ++i ) {
			if ( pose2.chain(resno2)==options_.chains_to_ignore_[i] ) return false; //If this chain is to be ignored, return false.
		}
	}

	if ( options_.residues_to_ignore_.size()>0 ) {
		for ( core::Size i=1, imax(options_.residues_to_ignore_.size()); i<=imax; i++ ) {
			if ( resno2 == options_.residues_to_ignore_[i] ) return false; //If this residue is to be ignored, return false.
		}
	}

	std::string atname2( utility::strip( pose2.residue(resno2).atom_name(atomno), " " ) );

	if ( use_this_atom( pose1.residue(resno), pose2.residue(resno2), atomno, atname2, atname2 ) ) return true;

	//Check extra atoms:
	for ( core::Size i(1), imax(extra_atom_list.size()); i<=imax; ++i ) {
		if ( resno!=static_cast<signed int>(extra_atom_list[i].rsd()) ) continue;
		if (
				pose1.residue(resno).has( extra_atom_list[i].atom() ) &&
				pose2.residue(resno2).has( extra_atom_list[i].atom() )
				) {
			return true;
		}
	}

	return false;
}

/// @brief Check that symmetry options have been set sensibly.
void
EnergyBasedClusteringProtocol::do_option_checks() const {

	//Check whether the v_mutate_to_ala flag has been set:
	if ( options_.mutate_to_ala_ ) {
		TR << "The \"mutate_to_ala\" option was set.  Input structures will be mutated to a chain of (alpha-D-, alpha-L-, or beta-3-) alanines (with the exception of cysteine residues)." << std::endl;
	}

	//Check that the user hasn't specified the same chain multiple times in the chains_to_ignore option:
	if ( options_.chains_to_ignore_.size()>1 ) {
		for ( core::Size i(2), imax( options_.chains_to_ignore_.size() ); i<=imax; ++i ) {
			for ( core::Size j=1; j<=i; ++j ) {
				runtime_assert_string_msg( options_.chains_to_ignore_[i] != options_.chains_to_ignore_[j],
					"Error!  The same chain must not be specified multiple times with the \"chains_to_ignore\" option.  Crashing gracelessly." );
			}
		}
	}

	//Check that the user hasn't specified the same residue multiple times in -v_ignoreresidue:
	if ( options_.residues_to_ignore_.size()>1 ) {
		for ( core::Size i(2), imax(options_.residues_to_ignore_.size()); i<=imax; i++ ) {
			for ( core::Size j=1; j<i; j++ ) {
				runtime_assert_string_msg( options_.residues_to_ignore_[i] != options_.residues_to_ignore_[j],
					"Error!  The same residue must not be specified multiple times with the \"residues_to_ignore\" option.  Crashing.");
			}
		}
	}

	//Check whether the v_homooligomer_swap flag has been set:
	if ( options_.homooligomer_swap_ ) {
		runtime_assert_string_msg( options_.cluster_by_ == EBC_bb_cartesian, "Error!  The \"homooligomer_swap\" option can currently only be used with the \"bb_cartesian\" clustering mode!  Crashing.\n");
		runtime_assert_string_msg( !options_.cyclic_, "Error!  The \"homooligomer_swap\" option is not currently compatible with the \"cyclic\" flag!  Crashing.\n" );
		TR << "The \"homooligomer_swap\" option was specified.  When calculating RMSD values, all permutations of chains will be considered in the alignment." << std::endl;
	}

	//Check the clustering radius:
	runtime_assert_string_msg( options_.cluster_radius_ > 0.0,
		"Error!  The clustering radius must be greater than zero.  Crashing");
	{
		std::string unitsstring;
		if ( options_.cluster_by_ == EBC_bb_cartesian ) unitsstring = " Angstroms";
		else if ( options_.cluster_by_ == EBC_bb_dihedral ) unitsstring = " degrees";
		TR << "Using a cluster radius of " << options_.cluster_radius_ << unitsstring << "." << std::endl;
	}

	//Check kbt:
	/*runtime_assert_string_msg( options_.kbt_ > 0.0, "Error!  The value of k_B*T specified with the \"kbt\" option must be greater than zero.  Crashing.");
	if(options_.weight_by_energy_) {
	TR << "Weighting structures by energy when calculating cluster centers.  Setting k_B*T=" << options_.kbt_ << "." << std::endl;
	}*/

	if ( options_.cyclic_symmetry_ ) {
		runtime_assert_string_msg( options_.cyclic_, "Error!  The \"cyclic_symmetry\" option requires \"cyclic\" to be true." );
		runtime_assert_string_msg( options_.cyclic_symmetry_ > 1, "Error!  The \"cyclic_symmetry\" option requires a setting greater than 1." );
		if ( options_.cyclic_symmetry_ > 1 ) TR << "Setting cyclic symmetry to " << options_.cyclic_symmetry_ << "." << std::endl;
	}
	if ( options_.cyclic_symmetry_mirroring_ ) {
		runtime_assert_string_msg( options_.cyclic_, "Error!  The \"cyclic_symmetry_mirroring\" option requires \"cyclic\" to be true." );
		runtime_assert_string_msg( options_.cyclic_symmetry_ > 1 && options_.cyclic_symmetry_ % 2 == 0, "Error!  The \"cyclic_symmetry_mirroring\" option requires that the \"cyclic_symmetry\" option is also used, and set to a value divisible by 2." );
		TR << "Enabling mirror symmetry." << std::endl;
	}
	if ( options_.cyclic_symmetry_threshold_specified_ ) {
		runtime_assert_string_msg( options_.cyclic_, "Error!  The \"-cyclic_symmetry_threshold\" option requires \"cyclic\" to be true." );
		runtime_assert_string_msg( options_.cyclic_symmetry_ > 1, "Error!  The \"cyclic_symmetry_threshold\" option requires that the \"-cyclic_symmetry\" option is set to more than 1." );
		TR << "Setting cyclic symmetry threshold to " << options_.cyclic_symmetry_threshold_ << "." << std::endl;
	}

	//Checks releated to v_cluster_cyclic_permutations:
	if ( options_.cluster_cyclic_permutations_ ) {
		runtime_assert_string_msg( options_.cyclic_, "Error!  The \"cluster_cyclic_permutations\" option cannot be used without the \"cyclic\" option.  Crashing gracelessly.");
		TR << "Cyclic permutations will be tried when aligning structures during clustering (\"cluster_cyclic_permutations\" option)." << std::endl;
	}

}

/// @brief Get the global scorefunction, and set up constraints
/// energies appropriately.
core::scoring::ScoreFunctionOP
EnergyBasedClusteringProtocol::set_up_scorefunction() const {
	using namespace core::scoring;

	core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );

	if ( !options_.cst_files_.empty() || !constraints_file_contents_.empty() ) { //If a constraints file has been specified by the user, turn on the atom_pair, angle, and dihedral constraint weights unless otherwise on.
		TR << "Turning on constraint weights, since a constraint file was specified." << std::endl;
		if ( sfxn->get_weight(atom_pair_constraint) < 1.0e-6 ) sfxn->set_weight(atom_pair_constraint, 1.0);
		if ( sfxn->get_weight(angle_constraint) < 1.0e-6 ) sfxn->set_weight(angle_constraint, 1.0);
		if ( sfxn->get_weight(dihedral_constraint) < 1.0e-6 ) sfxn->set_weight(dihedral_constraint, 1.0);
	}

	//Alter the scoring function for v_cyclic:
	if ( options_.cyclic_ ) {
		TR << "Setting constraint weights for a peptide bond between the N- and C-termini (\"cyclic\" option)." << std::endl;
		if ( sfxn->get_weight(atom_pair_constraint) < 1.0e-6 ) sfxn->set_weight(atom_pair_constraint, 1.0);
		if ( sfxn->get_weight(angle_constraint) < 1.0e-6 ) sfxn->set_weight(dihedral_constraint, 1.0);
		if ( sfxn->get_weight(dihedral_constraint) < 1.0e-6 ) sfxn->set_weight(angle_constraint, 1.0);
	}

	return sfxn;
}

/// @brief Given a vector of strings of format resnum:atomname, parse this out into a vector of NamedAtomIDs.
void
EnergyBasedClusteringProtocol::parse_extra_atom_list (
	utility::vector1 < core::id::NamedAtomID > &extra_atom_list
) const {

	if ( options_.extra_rms_atoms_.empty() ) return; //Do nothing if there's no list of extra atoms.
	extra_atom_list.clear();

	utility::vector1 <std::string> const & tags = options_.extra_rms_atoms_;

	//Parse each tag:
	TR << "The following additional atoms will be used in the RMSD calculation:\n";
	for ( core::Size i(1), imax(tags.size()); i<=imax; ++i ) { //Loop through each tag
		core::Size const colonposition( tags[i].find( ':' ));
		std::string const resstring( tags[i].substr( 0, colonposition ) );
		auto const res( static_cast<core::Size>( atoi( resstring.c_str() ) ) ); //The residue number
		std::string const atomname( tags[i].substr( colonposition + 1) ); //The atom name
		TR << "\tResidue " << res << ", Atom " << atomname << "\n";
		extra_atom_list.push_back( core::id::NamedAtomID( atomname, res ) );
	}

	TR << std::endl;
}

/// @brief Function to add user-specified constraints (specified with a CST file) to a pose.
void
EnergyBasedClusteringProtocol::add_user_constraints (
	core::pose::Pose &mypose
) const {
	core::Size const n_files( options_.cst_files_.size() ); //Number of constraints files.

	if ( !n_files ) return; //Do nothing if no constraint file is specified.

	//Add each constraint file to the pose:
	for ( core::Size i(1); i<=n_files; ++i ) {

		protocols::constraint_movers::ConstraintSetMoverOP cst_maker ( new protocols::constraint_movers::ConstraintSetMover );
		std::string const & cstfile ( options_.cst_files_[i] );

		cst_maker->constraint_file(cstfile);
		cst_maker->add_constraints(true); //Add constraints to anything else already there.
		cst_maker->apply(mypose);
	}
}

/// @brief Swap chains around in a pose (circular permutation of chain indices).
void
EnergyBasedClusteringProtocol::swap_chains(
	core::pose::Pose &currentpose, //Input and output
	core::Size const permutation_number //Input -- the chain perturbation number
) const {
	if ( permutation_number == 1 ) return; //Do nothing if this is the first perturbation (no rearrangement of the list of chains).

	core::Size const chaincount( currentpose.num_chains() );
	if ( chaincount < 2 ) return; //Do nothing if there's only one chain.

	utility::vector1 < core::Size > chainlist( chaincount ); //List of chains -- order will be perturbed.
	utility::vector1 < core::Size > startaa_list( chaincount ); //List of starting amino acids of chains
	utility::vector1 < core::Size > endaa_list( chaincount ); //List of starting amino acids of chains
	utility::vector1 < core::Size > startindex_list; //List of the starting index in the parentvect for each chain

	for ( core::Size i=1; i<=chainlist.size(); i++ ) chainlist[i]=i; //Initialize the list

	//Reorder the list of chains:
	if ( permutation_number>1 ) {
		for ( core::Size i=2; i<=permutation_number; i++ ) {
			std::next_permutation(chainlist.begin(), chainlist.end());
		}
	}

	//Get the starting and ending aa index for each chain:
	for ( core::Size i=1; i<=chaincount; i++ ) {
		get_start_and_end_aa(chainlist[i], currentpose, startaa_list[i], endaa_list[i]);
	}

	//A pose into which residues will be copied:
	core::pose::Pose newpose;
	for ( core::Size i(1); i<=chaincount; ++i ) {
		core::pose::append_subpose_to_pose( newpose, currentpose, startaa_list[i], endaa_list[i], true );
	}

	currentpose=newpose;
	currentpose.update_residue_neighbors();
}

/// @brief Get the start and end indices of chain chain_index, and store these in startaa and endaa.
void
EnergyBasedClusteringProtocol::get_start_and_end_aa(
	core::Size const chain_index,
	core::pose::Pose const &pose,
	core::Size & startaa,
	core::Size & endaa
) const {
	bool in_chain(false);
	startaa = 0; endaa = 0;
	for ( core::Size i(1), nres(pose.total_residue()); i<=nres; ++i ) {
		if ( !in_chain && pose.chain(i) == chain_index ) {
			in_chain = true;
			startaa = i;
		} else if ( in_chain && (pose.chain(i) != chain_index || i == nres) ) {
			in_chain = false;
			if ( i == nres ) endaa = i;
			else endaa = i - 1;
		}
		if ( startaa!=0 && endaa!=0 ) return;
	}
	// We shouldn't reach here.
	utility_exit_with_message( "Error in protocols::cluster::energy_based_clustering::EnergyBasedClusteringProtocol::get_start_and_end_aa(): Failed to find start and end indices.  Is the requested chain in the pose?" );
}

/// @brief Get the start and end indices of chain chain_index, and store these in startaa and endaa.
void
EnergyBasedClusteringProtocol::get_start_and_end_indices(
	core::Size const chain_index,
	core::pose::Pose const &pose,
	core::Size &startindex,
	core::Size &endindex
) const {
	startindex = 0;
	endindex = 0;
	bool in_chain(false);
	core::Size counter(1);
	for ( core::Size ir(1), nres(pose.total_residue()); ir<=nres; ++ir ) {
		if ( !in_chain && pose.chain(ir) == chain_index ) {
			startindex = counter;
			in_chain = true;
			counter += alignment_atoms_in_res( pose, ir );
		} else if ( in_chain && ( pose.chain(ir) != chain_index || ir == nres ) ) {
			endindex = (ir == nres ? counter + alignment_atoms_in_res(pose, ir) - 1 : counter - 1);
			in_chain = false;
		}
		if ( startindex != 0 && endindex != 0 ) return;
	}
	// We shouldn't reach here.
	utility_exit_with_message( "Error in protocols::cluster::energy_based_clustering::EnergyBasedClusteringProtocol::get_start_and_end_indices(): Failed to find start and end indices.  Is the requested chain in the pose?" );
}

/// @brief Function to form the disulfides based on user-specified disulfide positions (if provided), or based
/// on automatic detection (if not).
void
EnergyBasedClusteringProtocol::make_disulfides (
	core::pose::Pose &mypose
) const {
	core::Size const n_disulf_positions( options_.disulfide_positions_.size() ); //Number of user-specified disulfide positions.
	if ( n_disulf_positions ) { //If the user has specified disulfides
		runtime_assert_string_msg( n_disulf_positions % 2 == 0, "Error in protocols::cluster::energy_based_clustering::EnergyBasedClusteringProtocol::make_disulfides(): The number of disulfide positions specified must be divisible by two." );
		utility::vector1 < std::pair < core::Size, core::Size > > disulfpairs;
		for ( core::Size i(1), imax( n_disulf_positions ); i<=imax; i+=2 ) {
			disulfpairs.push_back( std::pair< core::Size, core::Size > ( options_.disulfide_positions_[i], options_.disulfide_positions_[i+1] ) );
		}
		mypose.conformation().fix_disulfides(disulfpairs);
	} else { //If the user has not specified disulfides
		mypose.conformation().detect_disulfides();
	}
}

/// @brief Function to mutate a pose to a chain of alanines.
/// @details This does not mutate cysteine residues involved in disulfide bonds.  Note that this
/// assumes that disulfide bonds have already been built in the pose.
void
EnergyBasedClusteringProtocol::mutate_to_alanine(
	core::pose::Pose &mypose
) const {

	for ( core::Size ir(1), irmax(mypose.total_residue()); ir<=irmax; ir++ ) { //Loop through all residues
		if ( !mypose.residue_type(ir).is_beta_aa() && !mypose.residue_type(ir).is_alpha_aa() && !mypose.residue_type(ir).is_peptoid() ) continue; //Skip non-amino acid residues
		if ( mypose.residue_type(ir).is_disulfide_bonded() ) continue; //Skip cysteines involved in disulfides
		std::string aaname( "ALA" ); //By default, mutate to alanine.
		if ( mypose.residue_type(ir).is_alpha_aa() ) {
			if ( mypose.residue_type(ir).is_d_aa() ) {
				aaname = "DALA"; //If it's a D-amino acid, mutate to D-alanine
			}
		} else if ( mypose.residue_type(ir).is_beta_aa() ) { //If it's a beta-amino acid, mutate to beta-3-alanine.
			if ( mypose.residue_type(ir).is_d_aa() ) {
				aaname="DB3A";
			} else {
				aaname="B3A";
			}
		} else if ( mypose.residue_type(ir).is_peptoid() ) { //If it's a peptoid, mutate to sarcosine.
			aaname="GLY:N_Methylation"; //Sarcosine is "GLY:N_Methylation" in Rosetta.
		} else if ( mypose.residue_type(ir).is_oligourea() ) { //If it's an oligourea, mutate to OU3_ALA.
			aaname="OU3_ALA"; //An alanine-like oligourea.
		}
		protocols::simple_moves::MutateResidue mutres(ir, aaname); //Mutate residue mover
		mutres.apply(mypose); //Apply the mutation
	}

	mypose.update_residue_neighbors();
}

/// @brief Function to check whether two poses have matching classes of residues at matching positions:
void
EnergyBasedClusteringProtocol::check_backbones_match (
	core::pose::Pose const &pose1,
	core::pose::Pose const &pose2
) const {
	static const std::string errmsg("Error in protocols::cluster::energy_based_clustering::EnergyBasedClusteringProtocol::backbones_match(): ");

	for ( core::Size ir=1, nres=pose1.total_residue(); ir<=nres; ++ir ) {
		if ( is_in_list( pose1.chain(ir), options_.chains_to_ignore_ ) ) continue; //Skip this residue if it's in the chain ignore list.
		if ( is_in_list( ir, options_.residues_to_ignore_) ) continue; //Skip this residue if it's in the ignore list.

		core::chemical::ResidueType const &type1( pose1.residue_type(ir) );
		core::chemical::ResidueType const &type2( pose2.residue_type(ir) );
		if ( (type1.is_alpha_aa() || type1.is_peptoid()) && !(type2.is_alpha_aa() || type2.is_peptoid()) ) {
			utility_exit_with_message( errmsg + "Residue " + std::to_string(static_cast<unsigned long>(ir)) + " of first pose is an alpha-amino acid or peptoid, but the corresponding residue of a subsequent pose is not!" );
		} else if ( type1.is_beta_aa() && !type2.is_beta_aa() ) {
			utility_exit_with_message( errmsg + "Residue " + std::to_string(static_cast<unsigned long>(ir)) + " of first pose is a beta-amino acid, but the corresponding residue of a subsequent pose is not!" );
		} else if ( type1.is_gamma_aa() && !type2.is_gamma_aa() ) {
			utility_exit_with_message( errmsg + "Residue " + std::to_string(static_cast<unsigned long>(ir)) + " of first pose is a gamma-amino acid, but the corresponding residue of a subsequent pose is not!" );
		} else if ( type1.is_NA() && !type2.is_NA() ) {
			utility_exit_with_message( errmsg + "Residue " + std::to_string(static_cast<unsigned long>(ir)) + " of first pose is a nucleic acid residue, but the corresponding residue of a subsequent pose is not!" );
		} else if ( (type1.is_ligand() || type1.is_metal()) && !(type2.is_ligand() || type2.is_metal()) ) {
			utility_exit_with_message( errmsg + "Residue " + std::to_string(static_cast<unsigned long>(ir)) + " of first pose is a metal or ligand, but the corresponding residue of a subsequent pose is not!" );
		}
	}
}

/// @brief Given an input pose, store only the relevant data needed for clustering.
/// @details The "alignmentdata" array is only used for Cartesian clustering.
void
EnergyBasedClusteringProtocol::storeposedata(
	core::pose::Pose const &pose,
	utility::vector1 < core::Real > &posedata,
	utility::vector1< numeric::xyzVector< core::Real > > &alignmentdata, //Only for Cartesian clustering: x,y,z coordinates of atoms to be used for alignment.
	utility::vector1 < core::Real > &dihedral_mode_reconstruction_data, //Only for reconstructing pose in dihedral mode.
	EBC_ClusterType const clustermode,
	utility::vector1 < core::id::NamedAtomID > const &extra_atom_list
) const {
	posedata.clear();
	alignmentdata.clear();
	core::Size const nres( pose.total_residue() );

	//Count the number of atoms in the alignment data:
	if ( clustermode == EBC_bb_cartesian ) {
		core::Size alignmentdatasize( 0 );
		for ( core::Size ir=1; ir<=nres; ir++ ) {
			if ( is_in_list( pose.chain(ir), options_.chains_to_ignore_ ) ) continue; //Skip this residue if it's in the chain ignore list.
			if ( is_in_list( ir, options_.residues_to_ignore_) ) continue; //Skip this residue if it's in the ignore list.
			if ( pose.residue_type(ir).is_ligand() || pose.residue_type(ir).is_metal() ) continue; //Skip this if not an alpha or beta amino acid
			for ( core::Size ia=1,iamax=pose.residue(ir).natoms(); ia<=iamax; ia++ ) {
				if ( use_in_rmsd(pose,pose,ir,ia, extra_atom_list) ) alignmentdatasize++;
			}
		}
		alignmentdata.resize( alignmentdatasize );
		posedata.reserve( pose.total_atoms() );
	} else if ( clustermode == EBC_bb_dihedral ) {
		posedata.reserve( pose.total_residue() * 3 ); //Best initial guess at ultimate size of posedata vector: three torsions per residue.  Could be more or less, depending on details of residue types and connectivity.
		if ( options_.rebuild_all_in_dihedral_mode_ ) dihedral_mode_reconstruction_data.reserve(pose.total_atoms());
	}


	core::Size icount=0;
	for ( core::Size ir=1, nres=pose.size(); ir<=nres; ir++ ) { //Loop through all residues
		switch(clustermode) { //Depending on the clustering mode do different things
		case EBC_bb_cartesian : //backbone Cartesian coordinate clustering
			for ( core::Size ia=1, natom=pose.residue(ir).natoms(); ia<=natom; ia++ ) {
				posedata.push_back(pose.residue(ir).atom(ia).xyz().x());
				posedata.push_back(pose.residue(ir).atom(ia).xyz().y());
				posedata.push_back(pose.residue(ir).atom(ia).xyz().z());
				bool useresidue(true);
				if ( pose.residue_type(ir).is_ligand() || pose.residue_type(ir).is_metal() ) useresidue=false; //Skip if this isn't an alpha or beta amino acid
				else {
					if ( is_in_list(ir, options_.residues_to_ignore_) ) { useresidue=false; }
					if ( is_in_list(pose.chain(ir), options_.chains_to_ignore_ ) ) { useresidue=false; }
				}
				if ( useresidue && use_in_rmsd(pose, pose, ir, ia, extra_atom_list) ) { //Store data for Cartesian alignment:
					icount++;
					alignmentdata[icount].x( pose.residue(ir).atom(ia).xyz().x() );
					alignmentdata[icount].y( pose.residue(ir).atom(ia).xyz().y() );
					alignmentdata[icount].z( pose.residue(ir).atom(ia).xyz().z() );
				}
			}
			break;
		case EBC_bb_dihedral : //backbone dihedral clustering
			if ( options_.rebuild_all_in_dihedral_mode_ ) { //Store xyz coords for later rebuilding.
				for ( core::Size ia=1, natom=pose.residue(ir).natoms(); ia<=natom; ia++ ) {
					dihedral_mode_reconstruction_data.push_back( pose.residue(ir).atom(ia).xyz().x() );
					dihedral_mode_reconstruction_data.push_back( pose.residue(ir).atom(ia).xyz().y() );
					dihedral_mode_reconstruction_data.push_back( pose.residue(ir).atom(ia).xyz().z() );
				}
			}

			if ( pose.residue_type(ir).is_ligand() || pose.residue_type(ir).is_metal() ) continue;
			if ( is_in_list(ir, options_.residues_to_ignore_) ) continue;
			if ( is_in_list(pose.chain(ir), options_.chains_to_ignore_) ) continue;

			for ( core::Size itors(1), itors_max(pose.residue(ir).mainchain_torsions().size()); itors<=itors_max; ++itors ) {
				if ( !pose.residue_type(ir).is_polymer() ) continue; //Shouldn't be possible, but just to be sure.
				if ( itors == 1 ) {
					if ( !pose.residue(ir).has_lower_connect() || pose.residue(ir).connected_residue_at_lower() == 0 ) continue;
				}
				if ( itors == itors_max || itors == itors_max - 1 ) {
					if ( !pose.residue(ir).has_upper_connect() || pose.residue(ir).connected_residue_at_upper() == 0 ) continue;
				}
				posedata.push_back( numeric::principal_angle_degrees( pose.residue(ir).mainchain_torsion(itors) ) );
			}
			break;
		default :
			break;
		} //End switch
	} //End loop through all residues

	/*
	if (clustermode==EBC_bb_dihedral) { //Additional data to be stored in the backbone dihedral clustering case:
	for (core::Size ir=1; ir<=pose.size(); ir++) { //Store side chain dihedrals, too, for rebuilding structures later
	if(pose.residue(ir).nchi()==0) continue;
	for(core::Size ichi=1; ichi<=pose.residue(ir).nchi(); ichi++) {
	posedata.push_back( numeric::principal_angle_degrees( pose.residue(ir).chi(ichi) ) );
	}
	}
	}
	*/
}

/// @brief Import all structures and score them, setting up derived data.
/// @details Performs disk access.
void
EnergyBasedClusteringProtocol::do_initial_import_and_scoring(
	core::Size &count,
	core::Real &lowestE,
	core::Size &lowestE_index,
	core::pose::Pose &firstpose,
	protocols::cyclic_peptide::CycpepSymmetryFilter const &symmfilter,
	utility::vector1 <core::Real> &poseenergies,
	utility::vector1 < utility::vector1 <core::Real> > &posedata,
	utility::vector1 < utility::vector1< numeric::xyzVector< core::Real > > > &alignmentdata,
	utility::vector1 < utility::vector1 <core::Real> > &dihedral_reconstruction_data,
	utility::vector1 < core::Size > &cluster_assignments,
	utility::vector1 < core::Size > &cluster_offsets,
	utility::vector1 < core::Size > &cluster_oligomer_permutations,
	core::scoring::ScoreFunctionOP sfxn,
	utility::vector1 < core::id::NamedAtomID > const &extra_atom_list
) const {

	core::import_pose::pose_stream::MetaPoseInputStream input( core::import_pose::pose_stream::streams_from_cmd_line() );

	//Loop for initial import and scoring:
	while ( input.has_another_pose() ) { //For every input structure
		count++;

		core::pose::Pose pose; //Create the pose
		input.fill_pose( pose ); //Import it

		if ( options_.cyclic_ ) {
			add_cyclic_constraints(pose); //If this is a cyclic peptide, add constraints for the terminal peptide bond:
			if ( options_.cluster_cyclic_permutations_ && options_.cluster_by_ == EBC_bb_cartesian && options_.use_CB_ ) {
				for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ir++ ) {
					if ( !pose.residue(ir).has("CB") ) { //Error if we don't have the same set of atoms on which to cluster in each residue.
						utility_exit_with_message("Error!  When using Cartesian clustering and clustering cyclic permutations, all residues in the input structures must have the same number of atoms to be used in the alignment.  If beta carbons are included, the sequence cannot contain glycine.  Crashing.");
					}
				}
			}
		}

		//Filter by symmetry:
		if ( options_.cyclic_symmetry_ > 1 ) {
			if ( !symmfilter.apply( pose ) ) {
				if ( options_.cyclic_symmetry_mirroring_ ) {
					TR << "\nStructure is not S" << options_.cyclic_symmetry_ << " symmetric.  Skipping." << std::endl;
				} else {
					TR << "\nStructure is not C" << options_.cyclic_symmetry_ << " symmetric.  Skipping." << std::endl;
				}
				count--;
				continue;
			}
		}

		if ( count % 100 == 0 ) {
			TR << count << " structures loaded...";
		}

		cluster_assignments.push_back(0); //Initially, every structure is assigned to cluster 0 (unassigned).
		if ( options_.cluster_cyclic_permutations_ ) cluster_offsets.push_back(0); //Initially, we assume that each structure may be aligned without any cyclic permutation, so these are all zero.
		if ( options_.homooligomer_swap_ ) cluster_oligomer_permutations.push_back(1); //Initially, we assume the first permutation of homooligomer subunits for each structure.

		if ( pose.num_jump() > 0 && options_.cluster_by_ == EBC_bb_dihedral ) {
			utility_exit_with_message( "Error!  Backbone dihedral clustering is not currently compatible with multi-chain PDB files.  Crashing." ); //TODO -- make this compatible!
		}

		if ( pose.num_jump() > 0 && options_.homooligomer_swap_ ) {
			const std::string chain1seq = pose.chain_sequence(1);
			for ( core::Size i=2, imax=pose.num_jump()+1; i<=imax; i++ ) {
				if ( pose.chain_sequence(i) != chain1seq ) {
					printf("Error!  When using the v_homooligomer_swap option, the lengths and sequences of all chains in the input structures must be identical.  Crashing.\n");
					fflush(stdout); exit(1);
				}
			}
		}

		add_user_constraints(pose); //This checks for a user-specified CST file, and does nothing if there isn't one.  (Adds constraints if there is one.)

		make_disulfides(pose); //Add user-specified disulfide bonds.

		if ( options_.prerelax_ ) {
			protocols::relax::FastRelax frlx(sfxn, options_.relax_rounds_);
			frlx.apply(pose); //Relax the pose if the user has specified that it be relaxed.
		}

		if ( options_.mutate_to_ala_ ) mutate_to_alanine(pose); //Mutate the pose to a chain of alanines if necesssary.

		(*sfxn)(pose); //Score the input pose
		if ( count==1 ) {
			firstpose=pose; //Store the first pose.
		} else {
			check_backbones_match(firstpose, pose);
		}
		poseenergies.push_back(pose.energies().total_energy()); //Store the pose energy

		//Store the pose data that will be used for clustering:
		posedata.push_back( utility::vector1< core::Real >() );
		alignmentdata.push_back( utility::vector1< numeric::xyzVector< core::Real > >() );
		if ( options_.rebuild_all_in_dihedral_mode_ || dihedral_reconstruction_data.size() == 0 ) dihedral_reconstruction_data.push_back(utility::vector1<core::Real>());
		storeposedata(pose, posedata[posedata.size()], alignmentdata[alignmentdata.size()], dihedral_reconstruction_data[dihedral_reconstruction_data.size()], options_.cluster_by_, extra_atom_list);

		if ( count==1 || pose.energies().total_energy() < lowestE ) {
			lowestE=pose.energies().total_energy();
			lowestE_index = count;
		}
	} //Loop over all input structures.
} //do_initial_import_and_scoring()


/// @brief Function to add cyclic constraints to a pose.
void
EnergyBasedClusteringProtocol::add_cyclic_constraints (
	core::pose::Pose &mypose
) const {
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	//Strip off termini:
	if ( mypose.residue_type(1).is_lower_terminus() ) {
		core::pose::remove_lower_terminus_type_from_pose_residue(mypose, 1);
	}
	if ( mypose.residue_type(mypose.total_residue()).is_upper_terminus() ) {
		core::pose::remove_upper_terminus_type_from_pose_residue(mypose, mypose.total_residue());
	}

	//Declare a chemical bond connecting the termini:
	protocols::cyclic_peptide::DeclareBond decbond;
	decbond.set( 1, "N", mypose.size(), "C", false );
	decbond.apply(mypose);

	{//Peptide bond length constraint:
		FuncOP harmfunc1 ( new HarmonicFunc( 1.3288, 0.02) );
		ConstraintCOP distconst1 ( new AtomPairConstraint (
			AtomID( mypose.residue(mypose.size()).atom_index("C") , mypose.size() ) ,
			AtomID( mypose.residue(1).atom_index("N") , 1) ,
			harmfunc1
			) );
		mypose.add_constraint (distconst1);
	}

	{ //Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		FuncOP circharmfunc1 ( new CircularHarmonicFunc( numeric::constants::d::pi, 0.02) );
		ConstraintCOP dihedconst1 ( new DihedralConstraint (
			AtomID( mypose.residue(mypose.size()).atom_index( (mypose.residue_type( mypose.size()).is_beta_aa() ? "CM" : "CA" ) ) , mypose.size() ),
			AtomID( mypose.residue(mypose.size()).atom_index("C") , mypose.size() ),
			AtomID( mypose.residue(1).atom_index("N") , 1) ,
			AtomID( mypose.residue(1).atom_index("CA") , 1) ,
			circharmfunc1
			) );
		mypose.add_constraint (dihedconst1);
	}

	{ //Peptide bond angle constraints:
		FuncOP circharmfunc2a ( new CircularHarmonicFunc( CNCa_ANGLE/180.0*numeric::constants::d::pi, 0.02) );
		FuncOP circharmfunc2b ( new CircularHarmonicFunc( CaCN_ANGLE/180.0*numeric::constants::d::pi, 0.02) );

		ConstraintCOP angleconst1 ( new AngleConstraint (
			AtomID( mypose.residue(mypose.size()).atom_index("C") , mypose.size() ),
			AtomID( mypose.residue(1).atom_index("N") , 1) ,
			AtomID( mypose.residue(1).atom_index("CA") , 1) ,
			circharmfunc2a
			) );
		ConstraintCOP angleconst2 ( new AngleConstraint (
			AtomID( mypose.residue(mypose.size()).atom_index( (mypose.residue(mypose.size()).has("CM")?"CM":"CA") ) , mypose.size() ),
			AtomID( mypose.residue(mypose.size()).atom_index("C") , mypose.size() ),
			AtomID( mypose.residue(1).atom_index("N") , 1) ,
			circharmfunc2b
			) );

		mypose.add_constraint (angleconst1);
		mypose.add_constraint (angleconst2);
	}
}



/****************************************************************************************************

/////////////////// OLD FUNCTIONS RELATED TO PRINCIPAL COMPONENT ANALYSIS ///////////////////////////
/////////////////// (I'll reinstate these someday.)                       ///////////////////////////

//Function to trim a vector of dihedral angles, and to add tranlation and Euler angle vectors for all the jumps in a pose to it.
void trim_and_add_jump_data (
utility::vector1 < core::Real > &myvect,
const core::Size numdihedrals,
const core::pose::Pose &mypose
) {
myvect.resize(numdihedrals + 6*mypose.num_jump()); //Discard extra information in myvect
if(mypose.num_jump() > 0) { //If there are jumps, store translation and rotation vectors for each jump in mypose.
for(core::Size i=numdihedrals+1, counter=1; i<=myvect.size(); counter++) {
numeric::xyzVector < core::Real > transvect = mypose.jump(counter).get_translation();
numeric::EulerAngles < core::Real > euler_angles(mypose.jump(counter).get_rotation());
for(core::Size j=0; j<3; j++) myvect[i++] = transvect[j]; //Note that xyzVectors are zero-based (unlike everything else in Rosetta -- yay for consistency).
myvect[i++] = euler_angles.phi_degrees();
myvect[i++] = euler_angles.theta_degrees();
myvect[i++] = euler_angles.psi_degrees();
//The following is just for debugging:
//printf("T=(%.2f,%.2f,%.2f)\nR=(%.2f,%.2f,%.2f)\n", transvect[0], transvect[1], transvect[2], euler_angles.phi_degrees(), euler_angles.theta_degrees(), euler_angles.psi_degrees());
}
}

return;
}

//Fuction to shift the current center of the cluster, weighted by energies.
//This also performs PCA analysis, returning a matrix of PCA vectors.
//The PCA vectors correspond to backbone dihedral angles, with jumps (tx, ty, tz, rx, ry, rz) appended at the end.
//Angles are in degrees, transforms are in angstroms.
void shift_center_and_PCA(
utility::vector1 <core::Real> &clustcenter,
utility::vector1 < utility::vector1 < core::Real > > &pca_vector_list,
utility::vector1 < core::Real > &coeff_list,
const core::Size clustcenterindex,
const utility::vector1 < utility::vector1 <core::Real> > &posedata,
const utility::vector1 <core::Real> &poseenergies,
const utility::vector1 <core::Size> &statesincluster, //This is the list of states assigned to this cluster, sorted by energies!
const core::Size currentclusterindex,
const core::Size clustmode,
const core::pose::Pose &firstpose,
const utility::vector1 <core::Size> &cluster_offsets, //Only used if clustering cyclic permutations together
const utility::vector1 <core::Size> &cluster_oligomer_permutations, //Only used if clustering permutations of homooligomer chains
utility::vector1 <core::id::NamedAtomID> const &extra_atom_list
) {
using namespace basic::options;
using namespace basic::options::OptionKeys;

printf("\tShifting center of cluster %lu.\n", currentclusterindex);

core::Size ndihedrals = count_bb_dihedrals(firstpose); //The number of dihedral angles to use in PCA analysis

printf("\t\tStarting with structure %lu.\n", clustcenterindex); fflush(stdout);

//Clear the storage containers:
pca_vector_list.clear();
coeff_list.clear();

core::Real weighting_accumulator = option[v_weightbyenergy]() ? exp(-poseenergies[clustcenterindex]/option[v_kbt]()) : 1.0; //The weighting accumulator (the denominator) starts off with the partition weight of the cluster center.
utility::vector1 <core::Real> coeff;
coeff.push_back(weighting_accumulator); //coeff[1] is the coefficient for the first state

utility::vector1 <core::Real> newclustcenter; //Used for clustmode=1 clustering
utility::vector1 <core::Real> deltaclustcenter; //The shift in the cluster center

FArray2D <core::Real> dummyarray; //Need this to satisfy requirements of storeposedata function.

core::pose::Pose clustcenterpose=firstpose;

utility::vector1 < utility::vector1 < core::Real > > Dphipsiomega; //A matrix whose columns are the phi, psi, and omega values of each state in the cluster, with the cluster centre subtracted off.  Columns will then be multiplied by the weighting coefficients.
Dphipsiomega.resize(statesincluster.size()); //A column for each state.
for(core::Size i=1, imax=statesincluster.size(); i<=imax; ++i) {
Dphipsiomega[i].resize( ndihedrals+6*clustcenterpose.num_jump() ); //A row for each backbone dihedral angle and jump parameter.
}
//Dphipsiomega.setlength(statesincluster.size(), ndihedrals+6*clustcenterpose.num_jump()); //This should have a column for each state, and a row for each backbone dihedral angle and jump parameter.

if(clustmode==1) {
pose_from_posedata(firstpose, clustcenterpose, clustmode, clustcenter);
storeposedata(clustcenterpose, newclustcenter, dummyarray, 2, extra_atom_list); //If clustermode is 1, newclustcenter holds the vector of backbone (and side chain) dihedrals at this point.  NOTE THAT THERE'S EXTRA INFORMATION AT THE END OF NEWCLUSTCENTER.
} else if(clustmode==2) newclustcenter=clustcenter;

trim_and_add_jump_data (newclustcenter, ndihedrals, clustcenterpose);  //Discard extra data and add jump data.

utility::vector1 <core::Size> not_angle_list; //A list of entries in the newclustcenter vector that are NOT angles (translations instead of rotations).
if(clustcenterpose.num_jump()>0) {
for(core::Size i=ndihedrals+1, counter=1, njump=clustcenterpose.num_jump(); counter<=njump; counter++) {
not_angle_list.push_back(i);
not_angle_list.push_back(i+1);
not_angle_list.push_back(i+2);
i+=6;
}
}

deltaclustcenter.resize(ndihedrals + 6*clustcenterpose.num_jump());
for(core::Size i=1; i<=deltaclustcenter.size(); i++) deltaclustcenter[i]=0.0; //Initialize the array
for(core::Size i=1; i<=newclustcenter.size(); i++) Dphipsiomega[1][i] = newclustcenter[i]; //Initialize the first column of Dphipsiomega

if(statesincluster.size() > 1) { //If there is more than one state in this cluster.
for (core::Size i=2; i<=statesincluster.size(); i++) { //Loop through all other states in the cluster.

printf("\t\tCalculating influence of structure %lu.\n", statesincluster[i]); fflush(stdout);

coeff.push_back(option[v_weightbyenergy]() ? exp(-poseenergies[statesincluster[i]]/option[v_kbt]()) : 1.0);
weighting_accumulator+=coeff[i]; //The weighting accumulator will ultimately be the demoninator in the partition function.

if(clustmode==1) {
utility::vector1 <core::Real> currentdata;
core::pose::Pose currentpose;
pose_from_posedata(firstpose, currentpose, clustmode, posedata[statesincluster[i]]);
if(option[v_homooligomer_swap]()) swapchains( currentpose, cluster_oligomer_permutations[statesincluster[i]] ); //Swap chains around.
storeposedata(currentpose, currentdata, dummyarray, 2, extra_atom_list); //currentdata is a vector of backbone dihedrals at this point.  NOTE EXTRA DATA AT END -- SIDE-CHAIN DIHEDRALS.

//If we're clustering cyclic permutations, it's necessary to offset the backbone dihedral vector:
if(option[v_cluster_cyclic_permutations]()) slidearound(currentdata, cluster_offsets[statesincluster[i]], firstpose.size(), firstpose);
trim_and_add_jump_data (currentdata, ndihedrals, currentpose); //Discard extra data and add jump data.

//Add this structure to Dphipsiomega:
for(core::Size j=1; j<=currentdata.size(); j++) Dphipsiomega[i][j]=currentdata[j];

//Add this to deltaclustcenter:
for(core::Size j=1, diffmode=2; j<=deltaclustcenter.size(); j++) {
diffmode = 2;
if(j>ndihedrals && is_in_list(j, not_angle_list)) diffmode=1; //If we're past the dihedral list, check whether this entry is a translation or a rotation.  If it's a translation, diffmode=1.
deltaclustcenter[j] += (coeff[i]*elementdiff(currentdata[j], newclustcenter[j], diffmode));
}
} else if (clustmode==2) {
utility::vector1 <core::Real> currentdata = posedata[statesincluster[i]];
//If we're clustering cyclic permutations, it's necessary to offset the backbone dihedral vector:
if(option[v_cluster_cyclic_permutations]()) slidearound(currentdata, cluster_offsets[statesincluster[i]], firstpose.size(), firstpose);
for(core::Size j=1; j<=deltaclustcenter.size(); j++) deltaclustcenter[j] += (coeff[i]*elementdiff(currentdata[j], newclustcenter[j], clustmode));
//Add this structure to Dphipsiomega:
for(core::Size j=1; j<=ndihedrals; j++) Dphipsiomega[i][j]=currentdata[j];
}
} //Looping through states in the cluster

for(core::Size i=1; i<=deltaclustcenter.size(); i++) newclustcenter[i]+=(deltaclustcenter[i]/weighting_accumulator); //Shift the cluster center.

check_angle_bounds(newclustcenter, not_angle_list);

//Subtract clustcenter from Dphipsiomega, and multiply each resultant column by coeff:
for(core::Size i=1, imax=statesincluster.size(); i<=imax; i++) { //Loop through all states
for(core::Size j=1, jmax=newclustcenter.size(); j<=jmax; j++) { //Loop through all backbone dihedrals
Dphipsiomega[i][j] = coeff[i]/weighting_accumulator*elementdiff(Dphipsiomega[i][j],newclustcenter[j], (is_in_list(j, not_angle_list)?1:2));
//Dphipsiomega[i-1][j-1] = coeff[i]*(Dphipsiomega[i-1][j-1]-newclustcenter[j]);
//printf("%.4f\t", Dphipsiomega[i-1][j-1]); //DELETE ME
}
//printf("\n"); fflush(stdout); //DELETE ME
}

//PCA analysis:
printf("\tPerforming principal component analysis for cluster %lu.\n", currentclusterindex); fflush(stdout);
//alglib::ae_int_t PCAresult = 0;
//alglib::real_1d_array variances;
//variances.setlength(ndihedrals+6*clustcenterpose.num_jump());
//alglib::real_2d_array PCmatrix;
//PCmatrix.setlength(ndihedrals+6*clustcenterpose.num_jump(),ndihedrals+6*clustcenterpose.num_jump());
//alglib::pcabuildbasis(Dphipsiomega, statesincluster.size(), ndihedrals+6*clustcenterpose.num_jump(), PCAresult, variances, PCmatrix);
std::pair< utility::vector1< utility::vector1 <core::Real> >,  utility::vector1 <core::Real> > pcaresult = numeric::principal_components_and_eigenvalues_ndimensions( Dphipsiomega, false );

//Store the variances and principal component vectors for output from this function.  Only the first N-1 should be nonzero, where N is the number of states in the cluster.
for(core::Size i=1; i<std::min(statesincluster.size(), ndihedrals+6*clustcenterpose.num_jump()); i++) {
coeff_list.push_back( pcaresult.second[i] );
//printf("%.8f\n", variances[i-1]);fflush(stdout); //DELETE ME
utility::vector1 < core::Real > PCA_vector = pcaresult.first[i];
//for(core::Size j=1; j<=ndihedrals+6*clustcenterpose.num_jump(); j++) PCA_vector.push_back(PCmatrix[j-1][i-1]);
pca_vector_list.push_back(PCA_vector);
}

if(clustmode==1) {
pose_from_posedata(firstpose, clustcenterpose, 2, newclustcenter, false);
storeposedata(clustcenterpose, newclustcenter, dummyarray, clustmode, extra_atom_list);
}

clustcenter=newclustcenter;

}

return;
}

//Function to write out a PCA file:
//Note that a PCA file is a list corresponding to backbone dihedrals, followed by entries for each jump.  Each jump has six entries: translation(x,y,z) and rotation(phi, theta, psi).
void output_PCA (
const utility::vector1 < utility::vector1 < core::Real > > &pca_vector_list,
const utility::vector1 < core::Real > &coeff_list,
const core::Size &clustnumber
) {
FILE* outputfile;
char filename [128];
sprintf(filename, "PCA_%lu.txt", clustnumber);
outputfile = fopen(filename, "w");
//printf("\tOpened %s for write.\n", filename); fflush(stdout); //DELETE ME

if(coeff_list.size() > 0) {
fprintf (outputfile, "%lu\t%lu\n", pca_vector_list[1].size(), coeff_list.size());
for(core::Size i=1; i<=coeff_list.size(); i++) fprintf(outputfile, "%.14e\t", coeff_list[i]);
fprintf(outputfile, "\n");
for(core::Size i=1; i<=pca_vector_list.size(); i++) {
for(core::Size j=1; j<=pca_vector_list[i].size(); j++) fprintf(outputfile, "%.14e\t", pca_vector_list[i][j]);
fprintf(outputfile, "\n");
}
} else {
fprintf (outputfile, "0\t0\n");
}

fclose(outputfile);
printf("\tWrote %s.\n", filename); fflush(stdout);
return;
}

****************************************************************************************************/

} //energy_based_clustering
} //protocols
