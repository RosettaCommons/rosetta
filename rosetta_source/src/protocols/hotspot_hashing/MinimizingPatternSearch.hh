// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author 

#include <string>
#include <fstream>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <devel/init.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <basic/Tracer.hh>

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>

#include <protocols/cluster/APCluster.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>

#include <protocols/hotspot_hashing/SearchPattern.hh>

#ifndef INCLUDED_protocols_hotspot_hashing_MinimizingPatternSearch_hh
#define INCLUDED_protocols_hotspot_hashing_MinimizingPatternSearch_hh

namespace protocols {
namespace hotspot_hashing {

basic::Tracer tr( "protocols.hotspot_hashing.MinimizingPatternSearch.hh" );

template <class T>
class Cluster
{
	public:
    T exemplar;
		utility::vector1<T> members;
};

class ResidueClusteringAlgorithm
{
	bool fxnal_group_only;

	public:
		ResidueClusteringAlgorithm() :
			fxnal_group_only(true)
		{
		}
		
		void clusterResidues(utility::vector1<core::conformation::ResidueCOP> residues, utility::vector1< Cluster<core::Size> > & clusters)
		{
			protocols::cluster::APCluster apcluster( residues.size() );

			// Store residue-residue rmsd to determine minimum rmsd for APCluster self-similarity
			utility::vector1< core::Real > residue_rms;
			residue_rms.reserve( (residues.size() * (residues.size() - 1)) / 2 );

			// Iterate across all residue pairs setting cluster rmsd
			for(core::Size i=1; i < residues.size(); ++i)
			{
				for(core::Size j=i+1; j <= residues.size(); ++j)
				{
					core::Real ij_rmsd = core::scoring::residue_sc_rmsd_no_super( residues[i], residues[j], fxnal_group_only);

					residue_rms.push_back(ij_rmsd);

					// Using negative rmsd values for clustering
					// RMSD relation is symmetric, must be set explicitly
					apcluster.set_sim( i, j, 0-ij_rmsd ); 
					apcluster.set_sim( j, i, 0-ij_rmsd );
				}
			}

			// find median rmsd
			sort( residue_rms.begin(), residue_rms.end() );
			core::Real median_rmsd = *(residue_rms.begin() + residue_rms.size() / 2);

			for(core::Size i = 1; i <= residues.size(); ++i)
			{
				// Use median rmsd for self-similarity
				apcluster.set_sim( i, i, 0 - median_rmsd );
			}

			// these settings are relatively low. clustering may not converge completely, but should be OK
			core::Size maxits( 1000 );
			core::Size convits( 10 );
			core::Real lambda( 0.75 );

			apcluster.cluster( maxits, convits, lambda );

			// Extract all exemplars and cluster members
			utility::vector1< core::Size > exemplars;

			clusters.resize(0);
			clusters.reserve(exemplars.size());

			apcluster.get_all_exemplars( exemplars );

			foreach(core::Size i, exemplars)
			{
				Cluster<core::Size> cluster;

				cluster.exemplar = i;
				apcluster.get_cluster_for(i, cluster.members);

				clusters.push_back(cluster);
			}
		}
};

std::string pre_outfile_suffix = ".premin.pdb";
std::string post_outfile_suffix = ".postmin.pdb";
std::string scorefile_suffix = ".score.txt";
std::string clusterfile_suffix = ".clusters.txt";

class MinimizingPatternSearch : public utility::pointer::ReferenceCount
{
	public:
		std::string output_basename;

		core::pose::Pose targetPose;
		core::conformation::ResidueCOP targetResidue;

		SearchPattern& searchPattern;

		core::scoring::ScoreFunctionOP scorefxn;

		bool constrain_residue_;
		bool sc_only_;

		std::ofstream pre_outfile;
		std::ofstream scorefile;
		std::ofstream clusterfile;
		std::ofstream post_outfile;

		MinimizingPatternSearch(std::string output_basename, core::pose::Pose targetPose, core::conformation::ResidueCOP targetResidue, SearchPattern& searchPattern, core::scoring::ScoreFunctionOP scorefxn, bool constrain_residue) :
			output_basename(output_basename),
			targetPose(targetPose),
			targetResidue(targetResidue),
			searchPattern(searchPattern),
			scorefxn(scorefxn),
			constrain_residue_(constrain_residue),
			sc_only_(true)
		{
		};

		void execute()
		{
			pre_outfile.open((output_basename + pre_outfile_suffix).c_str());
			post_outfile.open((output_basename + post_outfile_suffix).c_str());
			scorefile.open((output_basename + scorefile_suffix).c_str());
			clusterfile.open((output_basename + clusterfile_suffix).c_str());

			scorefile << "model pre_score post_score" << std::endl;

			utility::vector1<TransformPair> searchpoints = searchPattern.Searchpoints();

			tr.Info << "Initialized search pattern. Points: " << searchpoints.size() << std::endl;

			utility::vector1<core::conformation::ResidueCOP> minimized_residues;
			minimized_residues.reserve(searchpoints.size());

			core::Size transformindex = 1;

			foreach(TransformPair transform, searchpoints)
			{
				// Clone target pose, residue is cloned internally pose.
				core::pose::Pose pose(targetPose);

				core::Size residuejumpindex;
				core::Size residueindex;

				placeResidueAtTransform(pose, targetResidue, transform, residuejumpindex, residueindex);

				(*scorefxn)(pose);

				// Add additional scoring parameters
				if (sc_only_ == true)
				{
					core::pose::add_variant_type_to_pose_residue( pose, "SHOVE_BB", residueindex );
				}

				core::Real prescore = scoreResidue(pose, residuejumpindex, residueindex);
				
				if (constrain_residue_ == true)
				{
					addResidueCoordinateConstraints(pose, residuejumpindex, residueindex);
				}

				logPreMinPose(transformindex, pose, residuejumpindex, residueindex);

				// Minimize 
				minimizePose(pose, residuejumpindex, residueindex);

				if (constrain_residue_ == true)
				{
					pose.remove_constraints();
				}

				(*scorefxn)(pose);

				core::Real postscore = scoreResidue(pose, residuejumpindex, residueindex);
				logPostMinPose(transformindex, pose, residuejumpindex, residueindex);

				scorefile << transformindex << " " << prescore << " " << postscore << std::endl;

				core::conformation::ResidueCOP minimized_residue( pose.residue( residueindex ));
				minimized_residues.push_back(minimized_residue);

				transformindex++;
			}

			utility::vector1< Cluster<core::Size> > clusters;
			ResidueClusteringAlgorithm clusterer;
			clusterer.clusterResidues(minimized_residues, clusters);

			foreach(Cluster<core::Size> cluster, clusters)
			{
				foreach(core::Size member, cluster.members)
				{
					clusterfile << cluster.exemplar << " " << member << std::endl;
				}
			}

			pre_outfile.close();
			post_outfile.close();
			scorefile.close();
			clusterfile.close();
		};

		static void placeResidueAtTransform( core::pose::Pose & pose, core::conformation::ResidueCOP sourceResidue, TransformPair & transform, core::Size & residuejumpindex, core::Size & residueindex)
		{
			tr.Debug << "Placing at transform: " << transform << std::endl;

			// Places residue at last jump & residue number
			placeResidueOnPose(pose, sourceResidue);
			residueindex = pose.total_residue();
			residuejumpindex = pose.num_jump();

			core::kinematics::Stub upstreamstub = pose.conformation().upstream_jump_stub(residuejumpindex);
			core::kinematics::Jump residuejump = pose.jump(residuejumpindex);
			core::kinematics::Jump newjump(residuejump);
			
			core::conformation::ResidueCOP const residue( pose.residue( residueindex ));

			// Place residue in "canonical" position
			// Calculate transforms to move target atom to 0,0,0 and rotate into canonical position
			TransformPair residuetransform = residueStubCentroidTransform(residue);

			if (residuetransform.translation.length() != 0)
			{
				newjump.translation_along_axis(upstreamstub, residuetransform.translation, residuetransform.translation.length());
			}

			newjump.rotation_by_matrix(upstreamstub, Vector(), residuetransform.rotation);

			// Apply target transformation
			newjump.rotation_by_matrix(upstreamstub, Vector(), transform.rotation);
			if (transform.translation.length() != 0)
			{
				newjump.translation_along_axis(upstreamstub, transform.translation, transform.translation.length());
			}

			pose.set_jump(residuejumpindex, newjump);

			core::conformation::ResidueCOP const finalresidue( pose.residue( residueindex ));
			tr.Debug << "Placed residue at anchor location: " << finalresidue->xyz(residue->atom_index("CA")) << std::endl;
		};

		static TransformPair residueCanonicalTransform(core::conformation::ResidueCOP const residue)
		{
			// Canonical transform aligns CA atom to <0, 0, 0>
			// CA->CB vector along <1,0,0>
			// CA->C vector on the XY plane (<CA->C> * <0,0,1> == 0)
			
			Vector position = -residue->xyz(residue->atom_index("CA"));
			
			Vector zunit = Vector(0, 0, 1);
			Vector xunit = Vector(1, 0, 0);

			Vector cacb_vector = residue->xyz(residue->atom_index("CB")) + position;
			Matrix cacb_rotation = rotation_matrix( cacb_vector.cross(xunit), angle_of(cacb_vector, xunit));

			Vector cac_vector_prime = cacb_rotation * (residue->xyz(residue->atom_index("C")) + position);
			Vector cac_vector_zyprojection = Vector(0, cac_vector_prime.y(), cac_vector_prime.z());
			Matrix cac_rotation = rotation_matrix( cac_vector_zyprojection.cross(zunit), angle_of(cac_vector_zyprojection, zunit));

			return TransformPair(position, cac_rotation * cacb_rotation);
		}

		static TransformPair residueStubCentroidTransform(core::conformation::ResidueCOP const residue)
		{
			// Canonical transform aligns CA atom to <0, 0, 0>
			// CA->SC heavyatom centroid vector along <1,0,0>
			// CA->C vector on the XY plane (<CA->C> * <0,0,1> == 0)
			
			Vector position = -residue->xyz(residue->atom_index("CA"));
			
			Vector xunit = Vector(1, 0, 0);
			Vector yunit = Vector(0, 1, 0);

			Vector cacentroid_vector = residueStubCentroid(residue) + position;
			Matrix cacentroid_rotation = rotation_matrix( cacentroid_vector.cross(xunit), angle_of(cacentroid_vector, xunit));

			Vector cac_vector_prime = cacentroid_rotation * (residue->xyz(residue->atom_index("C")) + position);
			Vector cac_vector_zyprojection = Vector(0, cac_vector_prime.y(), cac_vector_prime.z());
			Matrix cac_rotation = rotation_matrix( cac_vector_zyprojection.cross(yunit), angle_of(cac_vector_zyprojection, yunit));

			return TransformPair(position, cac_rotation * cacentroid_rotation);
		}

		static Vector residueStubCentroid(core::conformation::ResidueCOP const residue)
		{
			Vector centroid;
			centroid = 0;

			if (residue->first_sidechain_atom() > residue->nheavyatoms())
			{
				//TODO Generate pseudocentroid from mainchain atoms
				return centroid;
			}

			for (core::Size i = residue->first_sidechain_atom(); i <= residue->nheavyatoms(); ++i)
			{
				centroid += residue->xyz(i);
			}

			centroid /= (1 + residue->nheavyatoms() - residue->first_sidechain_atom());

			return centroid;
		}

		static void placeResidueOnPose(core::pose::Pose & pose, core::conformation::ResidueCOP residue)
		{
			pose.append_residue_by_jump( *residue, pose.total_residue(), "", residue->atom_name(residue->nbr_atom()), true );

			tr.Debug << "Appended residue on pose. Residue: " << pose.total_residue() << " Jump: " << pose.num_jump() << " Anchor atom: " << residue->atom_name(residue->nbr_atom()) << std::endl;
		};

		static void addResidueCoordinateConstraints(core::pose::Pose & pose, core::Size residuejumpindex, core::Size residueindex)
		{
			using core::scoring::constraints::CoordinateConstraint;
			using core::scoring::constraints::HarmonicFunc;
			using core::id::AtomID;

			core::conformation::ResidueCOP const residue( pose.residue( residueindex ));

			pose.add_constraint( new CoordinateConstraint( AtomID(residue->atom_index("CA"), residueindex), AtomID(2, 1), pose.xyz(AtomID(residue->atom_index("CA"), residueindex)), new HarmonicFunc(0, 1.0)));
			pose.add_constraint( new CoordinateConstraint( AtomID(residue->atom_index("CB"), residueindex), AtomID(2, 1), pose.xyz(AtomID(residue->atom_index("CB"), residueindex)), new HarmonicFunc(0, 1.0)));
			pose.add_constraint( new CoordinateConstraint( AtomID(residue->atom_index("C"), residueindex), AtomID(2, 1), pose.xyz(AtomID(residue->atom_index("C"), residueindex)), new HarmonicFunc(0, 1.0)));
		}

		void minimizePose(core::pose::Pose & pose, core::Size residuejumpindex, core::Size residueindex)
		{
			core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
			movemap->set_jump( true );
			movemap->set_chi( false );
			movemap->set_bb( false );

			for (core::Size j = 1; j <= pose.num_jump(); ++j)
			{
				if (j == residuejumpindex)
				{
					movemap->set_jump(j, true);
				}
				else
				{
					movemap->set_jump(j, false);
				}
			}

			protocols::simple_moves::MinMover dfpMinTightTol( movemap, scorefxn, "dfpmin_armijo_nonmonotone_atol", 0.02, false /*use_nblist*/ );
			
			dfpMinTightTol.apply(pose);
		};

		core::Real scoreResidue( core::pose::Pose const & pose, core::Size residuejumpindex, core::Size residueindex)
		{
			core::Real score = 0;
			core::scoring::EnergyMap const score_weights = pose.energies().weights();

			for ( core::graph::EdgeListConstIterator egraph_it = pose.energies().energy_graph().get_node( residueindex )->const_edge_list_begin(); egraph_it != pose.energies().energy_graph().get_node( residueindex )->const_edge_list_end(); ++egraph_it)
      {
				core::scoring::EnergyEdge const * edge = static_cast< core::scoring::EnergyEdge const * > (*egraph_it);
				score += edge->dot( score_weights );
			}

			return score;
		}

		virtual void logPreMinPose(core::Size transformindex, core::pose::Pose & pose, core::Size residuejumpindex, core::Size residueindex)
		{
			utility::vector1< core::Size > outindex;
			outindex.push_back(residueindex);

			pre_outfile << "MODEL     "+utility::to_string( transformindex )+"\n";
			pose.dump_pdb( pre_outfile, outindex );
			pre_outfile << "ENDMDL\n";
		};

		virtual void logPostMinPose(core::Size transformindex, core::pose::Pose & pose, core::Size residuejumpindex, core::Size residueindex)
		{
			utility::vector1< core::Size > outindex;
			outindex.push_back(residueindex);

			post_outfile << "MODEL    "+utility::to_string( transformindex )+"\n";
			pose.dump_pdb( post_outfile, outindex );
			post_outfile << "ENDMDL \n";
		};
};

}
}

#endif
