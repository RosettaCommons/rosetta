// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <string>
#include <fstream>
#include <vector>

#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>

#include <protocols/cluster/APCluster.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/xyz.functions.hh>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/lexical_cast.hpp>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

basic::Tracer tr( "apps.pilot.lsmgridprobe.cc" );

typedef numeric::xyzMatrix< core::Real > Matrix;
typedef numeric::xyzVector< core::Real > Vector;

class TransformPair {
public:
	TransformPair (Vector translation, Matrix rotation) :
		translation(translation),
		rotation(rotation)
	{
	};

	TransformPair() :
		translation(0, 0, 0),
		rotation(Matrix::diag(1, 1, 1))
	{
	};

	Vector translation;
	Matrix rotation;
};

class VectorPair
{
	public:
		VectorPair() :
			position(),
			direction()
		{
		};

		VectorPair(Vector position, Vector direction) :
			position(position),
			direction(direction)
		{
		};

		VectorPair(const VectorPair& src) :
			position(src.position),
			direction(src.direction)
		{
		};

		Vector position;
		Vector direction;
};

std::ostream& operator<<(std::ostream &strm, const Vector &vector) {
  return strm << "[" << vector.x() << "," << vector.y() << "," << vector.z() << "]";
}

std::ostream& operator<<(std::ostream &strm, const Matrix &matrix) {
  return strm << "[" << matrix.row_x() << "," << matrix.row_y() << "," << matrix.row_z() << "]";
}

std::ostream& operator<<(std::ostream &strm, const TransformPair &tp) {
  return strm << "(" << tp.translation << "," << tp.rotation << ")";
}

std::ostream& operator<<(std::ostream &strm, const VectorPair &vp) {
  return strm << "(" << vp.position << "," << vp.direction << ")";
}

class SearchPattern : public utility::pointer::ReferenceCount
{
	public:
		virtual utility::vector1<TransformPair> Searchpoints() = 0;
};

class ConstPattern : public SearchPattern
{
	public:
		virtual utility::vector1<TransformPair> Searchpoints()
		{
			utility::vector1<TransformPair> searchpoints;
			searchpoints.push_back(TransformPair(Vector(10, 0, 0), numeric::z_rotation_matrix_degrees((core::Real)90)));
			searchpoints.push_back(TransformPair(Vector(0, 0, 0), numeric::z_rotation_matrix_degrees((core::Real)0)));
			searchpoints.push_back(TransformPair(Vector(-10, 0, 0), numeric::z_rotation_matrix_degrees((core::Real)270)));

			return searchpoints;
		}
};

class TestPattern : public SearchPattern
{
	public:
		TestPattern()
		{
		}

		virtual utility::vector1<TransformPair> Searchpoints()
		{
			utility::vector1<TransformPair> searchpoints;

			core::Real x = 0;
			core::Real y = 0;
			core::Real z = 0;

			for (x = -1; x <= 1; x += 1)
			{
				for (y = -1; y <= 1; y += 1)
				{
					for (z = -1; z <= 1; z += 1 )
					{
						for (core::Real angle = 0; angle < 360; angle += 180)
						{
							Vector translation = Vector(x, y, z);
							Matrix rotation = numeric::z_rotation_matrix_degrees(angle);

							TransformPair tp(translation, rotation);
							searchpoints.push_back(tp);
						}
					}
				}
			}

			return searchpoints;
		};
};

class LSMSearchPattern : public SearchPattern
{
	public:
		LSMSearchPattern(
				VectorPair lsmspec,
				core::Real angle_sampling,
				core::Real translocation_sampling,
				core::Real max_radius,
				core::Real distance_sampling,
				core::Real max_distance) : 
			lsmspec_(lsmspec),
			angle_sampling_(angle_sampling),
			translocation_sampling_(translocation_sampling),
			max_radius_(max_radius),
			distance_sampling_(distance_sampling),
			max_distance_(max_distance)
		{

		}

		virtual utility::vector1<TransformPair> Searchpoints()
		{
			utility::vector1<TransformPair> searchpoints;

			Vector xunit = Vector(0, 0, 1);
			Matrix normal_rotation = rotation_matrix( lsmspec_.direction.cross(xunit), angle_of(lsmspec_.direction, xunit));
			
			for (core::Real x = -max_radius_; x <= max_radius_; x += translocation_sampling_)
			{
				for (core::Real y = -max_radius_; y <= max_radius_; y += translocation_sampling_)
				{
					if(sqrt(x*x + y*y) <= max_radius_)
					{
						for (core::Real z = 0; z <= max_distance_; z += distance_sampling_)
						{
							for (core::Real angle = 0; angle < 360; angle += angle_sampling_)
							{
								Vector translation = Vector(x, y, z);
								Matrix rotation = numeric::x_rotation_matrix_degrees(angle);

								TransformPair tp(translation + lsmspec_.position, rotation * normal_rotation);
								searchpoints.push_back(tp);
							}
						}
					}
				}
			}

			return searchpoints;
		}

	private:
		VectorPair lsmspec_;
		core::Real angle_sampling_;
		core::Real translocation_sampling_;
		core::Real max_radius_;
		core::Real distance_sampling_;
		core::Real max_distance_;
};

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

		void placeResidueAtTransform( core::pose::Pose & pose, core::conformation::ResidueCOP sourceResidue, TransformPair & transform, core::Size & residuejumpindex, core::Size & residueindex)
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
			TransformPair residuetransform = residueCanonicalTransform(residue);

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

		TransformPair residueCanonicalTransform(core::conformation::ResidueCOP const residue)
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

		void placeResidueOnPose(core::pose::Pose & pose, core::conformation::ResidueCOP residue)
		{
			pose.append_residue_by_jump( *residue, pose.total_residue(), "", residue->atom_name(residue->nbr_atom()), true );

			tr.Debug << "Appended residue on pose. Residue: " << pose.total_residue() << " Jump: " << pose.num_jump() << " Anchor atom: " << residue->atom_name(residue->nbr_atom()) << std::endl;
		};

		void addResidueCoordinateConstraints(core::pose::Pose & pose, core::Size residuejumpindex, core::Size residueindex)
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

core::conformation::ResidueOP initializeResidue(core::pose::Pose pose, std::string resname)
{
	core::chemical::ResidueTypeSet const & residue_set ( pose.residue(1).residue_type_set() );
	core::chemical::ResidueType const & restype( residue_set.name_map( resname ) );

	return core::conformation::ResidueFactory::create_residue( restype );
}

bool tryParseLSMSpec(std::string lsmstring, core::Size & id, VectorPair & lsmspec)
{
	std::vector< std::string > vectors;
	boost::split( vectors, lsmstring, boost::is_any_of(":") );

	if(vectors.size() != 3){
		return false;
	}

	std::vector< std::string > position;
	std::vector< std::string > direction;

	boost::split( position, vectors[1], boost::is_any_of(",") );
	boost::split( direction, vectors[2], boost::is_any_of(",") );

	if(position.size() != 3 || direction.size() != 3){
		return false;
	}

	using boost::lexical_cast;

	id = lexical_cast<core::Size>(vectors[0]);

	lsmspec = VectorPair(
			Vector(
				lexical_cast<core::Real>(position[0]), 
				lexical_cast<core::Real>(position[1]), 
				lexical_cast<core::Real>(position[2])),
			Vector(
				lexical_cast<core::Real>(direction[0]), 
				lexical_cast<core::Real>(direction[1]), 
				lexical_cast<core::Real>(direction[2])));
	
	return true;
}

namespace basic{
	namespace options{
		namespace OptionKeys{
			basic::options::StringOptionKey const lsmdef("lsmdef");

			basic::options::RealOptionKey const angle_sampling("angle_sampling");
			basic::options::RealOptionKey const translocation_sampling("translocation_sampling");
			basic::options::RealOptionKey const max_radius("max_radius");
			basic::options::RealOptionKey const distance_sampling("distance_sampling");
			basic::options::RealOptionKey const max_distance("max_distance");

			basic::options::BooleanOptionKey const constrain_grid("constrain_grid");
		}
	}
}

int main( int argc, char * argv [] )
{
	using basic::options::option;
  using namespace basic::options::OptionKeys;

  option.add( lsmdef, "lsmdef");
  option.add( angle_sampling, "angle_sampling");
  option.add( translocation_sampling, "translocation_sampling");
  option.add( max_radius, "max_radius");
  option.add( distance_sampling, "distance_sampling");
  option.add( max_distance, "max_distance");

  option.add( constrain_grid, "constrain_grid").def(false);

	devel::init( argc, argv );

	// Turn on extra rotamers
	option[ packing::ex1::ex1 ].value( true );
	option[ packing::ex2::ex2 ].value( true );
	option[ packing::extrachi_cutoff ].value( 0 );

	// Necessary to make sure NANs in hbonding don't cause an exit
	option[ in::file::fail_on_bad_hbond ].value( false );

	// Read target pose
	std::string targetFilename;
	core::pose::Pose targetPose;
	if ( option[hotspot::target].user() )
	{
		targetFilename = option[ hotspot::target ]();
		core::import_pose::pose_from_pdb( targetPose, targetFilename );
	}
	else
	{
		utility_exit_with_message("Must specify a target structure using -target <filename>");
	}

	// Read starting residue
	std::string resname;
	if (option[ hotspot::residue ].user() ) {
		resname = option[ hotspot::residue ]()[1];
	}
	else
	{
		utility_exit_with_message("You must specify probe residue identity using -residue <3Letter>");
	}

	// Initialize residue representation
	core::conformation::ResidueOP residue;
	core::chemical::ResidueTypeSet const & residue_set ( targetPose.residue(1).residue_type_set() );
	core::chemical::ResidueType const & restype( residue_set.name_map( resname ) );
	residue = core::conformation::ResidueFactory::create_residue( restype );

	// Read LSM definition
	std::string lsmstring;
	core::Size lsmid;
	VectorPair lsmspec;

	if (option[ lsmdef ].user() ) {
		lsmstring = option[ lsmdef ]();
	}
	else{
		utility_exit_with_message("You must specify lsm spec (id:location:normal) using -lsmdef <id>:<x>,<y>,<z>:<x>,<y>,<z>");
	}

	if (!tryParseLSMSpec(lsmstring, lsmid, lsmspec)){
		utility_exit_with_message("Invalid lsm spec, use format (id:location:normal)  <id>:<x>,<y>,<z>:<x>,<y>,<z>");
	}

	// Generate output filename
	utility::file::FileName output_file(targetPose.pdb_info()->name());
	output_file.ext(resname + "." + boost::lexical_cast<std::string>(lsmid));

	utility::file::PathName	output_directory( option[out::path::all]() );

	utility::file::FileName output_basename(output_file, output_directory);

	// Initialize score function
	core::scoring::ScoreFunctionOP scorefxn;
	if( option[ score::weights ].user() ) {
		scorefxn = core::scoring::getScoreFunction();
	}
	else {
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score13" );
		scorefxn->set_weight( core::scoring::envsmooth, 0 );
	}

	core::scoring::methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.hbond_options().use_hb_env_dep( option[ hotspot::envhb]() );
	scorefxn->set_energy_method_options( options );
	
	LSMSearchPattern lsm(
			lsmspec,
			option[ angle_sampling ],
			option[ translocation_sampling ],
			option[ max_radius ],
			option[ distance_sampling ],
			option[ max_distance]
			);
	
	SearchPattern& pattern = lsm;

	MinimizingPatternSearch search(output_basename.name(), targetPose, residue, pattern, scorefxn, option[ constrain_grid ] );
	search.execute();
}
