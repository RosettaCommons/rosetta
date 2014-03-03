// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 MembraneProteinFactory.cc
///
/// @brief 	 MembraneProteinFactory
/// @details The membrane protein factory creates a single pose from various membrane proteins
///			 loaded on the front end and initialized as membrane proteins. This single framework
///			 will then be passed off to the MembraneHub (which coordinates I/O) and sent back to the protocol it was
///			 called from (usually in pose loading)
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

/// @brief Class: Membrane Protein Factory
/// @details Construct a membrane protein in Rosetta
/// @author Rebecca Alford (rfalford12@gmail.com)
///
/// The membrane protein factory constructs a membrane protein in Rosetta through the following steps:
///     (1) Construct fold tree with pre-defined membrane fold tree topology
///     (2) Add a MembraneConformation which maintains info about the membrane protein
///     (3) Add a virtual residue to the root of the foldtree to define the membrane
///     (4) Add a new virtual residue representing the embedding of each protein chain in the membrane
///
/// The membrane protein factory provides three constructors - single chain construction, multi chain construction via the resource
/// manager and multi chain construction via the options system (existing system). The new membrane protein framework
/// is not backwards compatible with the previous RosettaMembrane code (Yarov-Yaravoy et al. 2006) and uses a different set of options
///
/// Maintains the following invariants:
///   (1) For an n chain pose, will return a pose with n+1 chains (n chains plus a membrane chain)
///   (2) Initializes a correct membrane foldtree topology (initial)
///   (3) Membrane chain is always the n+1 chain
///   (4) Membrane fold tree is a valid fold tree
///
/// Notes:
///     - This class is not a mover. It will not move the backbone whatsoever (see MembraneRigidInitialMover)
///

#ifndef INCLUDED_core_membrane_MembraneProteinFactory_cc
#define INCLUDED_core_membrane_MembraneProteinFactory_cc

// Unit Headers
#include <core/membrane/MembraneProteinFactory.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/properties/LipidAccInfo.hh>
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/membrane/geometry/MembraneResidueFactory.hh>
#include <core/membrane/geometry/EmbeddingFactory.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/conformation/Conformation.hh>
#include <core/membrane/MembraneConformation.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace core::membrane;
using namespace core::membrane::properties;
using namespace core::pose;

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.membrane.MembraneProteinFactory");

/// @brief      Membrane Protein Factory
/// @details    Initializes a pose as a membrane protein in Rosetta

namespace core {
namespace membrane {


    //// Constructors /////////////////////

    /// @brief  Default Constructor (Private
    /// @details Private Default constructor (DO NOT USE)
    MembraneProteinFactory::MembraneProteinFactory() :
        utility::pointer::ReferenceCount()
    {}

    /// @brief   Resource Manager Constructor (No Options Specified)
    /// @details Construct a Membrane Protein Factory loading all required resources
    ///          from the resource manager using the default option settings
    ///
    /// @throws EXCN_Illegal_Arguments
    ///             if file is "" - cannot construct a membrane protein without chain refs for the resource manager
    ///
    /// @return MembraneProteinFactory
    MembraneProteinFactory::MembraneProteinFactory( std::string membrane_chains ) :
        utility::pointer::ReferenceCount(),
        prefix_file_(membrane_chains),
        fullatom_(true),
        include_lips_(false)
    {
        initialize_chains();
        initialize_resources();
    }

    /// @brief   Resource Manager Constructor (Non-Default)
    /// @details Construct a Membrane Protein Factory loading all required resources
    ///          from the resource manager
    ///
    /// @param  include_lips
    ///             load and include lipid accessibility data in scoring
    /// @param  membrane_chains
    ///             text file storing references to membrane chains in the membrane protein
    /// @param  fullatom
    ///             specified fullatom residue typeset
    ///
    /// @return MembraneProteinFactory
    MembraneProteinFactory::MembraneProteinFactory(
                                                   bool include_lips,
                                                   std::string membrane_chains,
                                                   bool fullatom

    ) :
        utility::pointer::ReferenceCount(),
        prefix_file_(membrane_chains),
        fullatom_(fullatom),
        include_lips_(include_lips)
    {
        initialize_chains();
        initialize_resources();
    }

    /// @brief   Interactive Constructor (With Default MP Options)
    /// @details Construct a Membrane Protein Factory loading all required resources
    ///          from user specified inputs
    ///
    /// @param chians
    ///         list of pose chains to compose into a single pose
    /// @param topologies
    ///         list of spannign topology objects for each chain
    /// @param embeddings
    ///         list of embedding data objects for each chain
    /// @param lipid_acc
    ///         list of lipid accessibility data objects for each chain
    ///
    /// @return MembraneProteinFactory
    MembraneProteinFactory::MembraneProteinFactory(
                                                   utility::vector1< PoseOP > chains,
                                                   utility::vector1< SpanningTopologyOP > topologies,
                                                   utility::vector1< core::membrane::util::EmbedConfigInfoOP > embeddings,
                                                   utility::vector1< LipidAccInfoOP > lipid_acc
                                                   ) :
    utility::pointer::ReferenceCount(),
        prefix_file_(""),
        fullatom_(true),
        include_lips_(false)
    {
        chains_ = chains;
        topologies_ = topologies;
        embeddings_ = embeddings;
        lipid_acc_ = lipid_acc;
    }

    /// @brief   Interactive Constructor (with non default options)
    /// @details Construct a Membrane Protein Factory loading all required resources
    ///          from user specified inputs
    ///
    /// @param  include_lips
    ///             load and include lipid accessibility data in scoring
    /// @param  fullatom
    ///             specified fullatom residue typeset
    /// @param  chians
    ///             list of pose chains to compose into a single pose
    /// @param  topologies
    ///             list of spannign topology objects for each chain
    /// @param  embeddings
    ///             list of embedding data objects for each chain
    /// @param  lipid_acc
    ///             list of lipid accessibility data objects for each chain
    ///
    /// @return MembraneProteinFactory
    MembraneProteinFactory::MembraneProteinFactory(
                                                   bool include_lips,
                                                   bool fullatom,
                                                   utility::vector1< PoseOP > chains,
                                                   utility::vector1< SpanningTopologyOP > topologies,
                                                   utility::vector1< core::membrane::util::EmbedConfigInfoOP > embeddings,
                                                   utility::vector1< LipidAccInfoOP > lipid_acc
                                                   ) :
        utility::pointer::ReferenceCount(),
        prefix_file_(""),
        fullatom_(fullatom),
        include_lips_(include_lips)

    {
        chains_ = chains;
        topologies_ = topologies;
        embeddings_ = embeddings;
        lipid_acc_ = lipid_acc;
    }

    /// @brief    Default Destructor
    /// @details
    ///
    /// @note
    MembraneProteinFactory::~MembraneProteinFactory()
    {}

    //// Public Member Functions ///////////////////

    /// @brief 	 Create Membrane Protein
    /// @details Create a membrne proteins from a series of loaded membrane proteins
    ///
    /// @return  Pose (as starting structure)
    core::pose::PoseOP
    MembraneProteinFactory::create_membrane_pose() {

        // Create new pose, build and return
        core::pose::PoseOP pose = new core::pose::Pose();
        build_pose(pose);

        return pose;
    }

    //// Private Member Functions //////////////////

    /// @brief      Initialize Chains
    /// @details    Initialize membrane chains from initialized prefix file provided in the constructor
    ///
    /// @throws     Argument exception if chain list not specified (also this is well docuemnted, no excuse)
    void
    MembraneProteinFactory::initialize_chains() {

        using namespace basic::options;
        using namespace core::membrane::util;

        // Ensure Prefix file is specified
        if ( prefix_file_.compare("") == 0 ) {
                throw new EXCN_NonMembrane("Cannot construct a membrane pose if membrane chains file is not provided");
        }

        // Grab file and create stream
        std::string infile = prefix_file_;
        utility::io::izstream stream (infile);

        std::string line;
        std::string desc;
        stream.open(infile);

        if (stream) {

            // Grab the first line
            getline(stream, line);

            core::Size i = 1;
            while ( !stream.eof() ) {

                std::istringstream l(line);
                l >> desc;
                chains_map_.insert( std::pair< core::Size, std::string >( i, desc ) );
                getline(stream, line);
                i++;
            }

        } else {
            throw new EXCN_Illegal_Arguments("Cannot open file " + infile );
        }

    }

    /// @brief Build Pose
    /// @details Create pose containing membrane/embedding residues from multi-
    /// chain input.
    /// @throws EXCN_Resource_Manager, EXCN_Membrane_Bounds
    void
    MembraneProteinFactory::build_pose( core::pose::PoseOP pose ) {

        using namespace core::conformation;
        using namespace core::kinematics;
        using namespace core::membrane::geometry;
        using namespace core::membrane;
        using namespace core::pose;

        // Keep Track of Data for the foldtree
        utility::vector1< std::pair< int, int > > embres_map; // chain num, jump anchor, embedding residue
        int nchains = (int) chains_.size();
        embres_map.resize(nchains);

        TR << "Printing number of chains at the top of the build pose method" << nchains << std::endl;

        // Create multi chain pose and write teh chain map
        for ( int i = 1; i <= nchains; i++ ) {

            // Create a given pose from multiple chains
            core::pose::PoseOP temp_pose = chains_[i];
            append_pose_to_pose(*pose, *temp_pose, true);
            TR << "Printing total number of resiues during iteration " << i << "of the append: " << pose->total_residue() << std::endl;
        }

        // Create and add Membrane residue with defaults
        core::Vector center(0, 0, 0);
        core::Vector normal(0, 0, 1);
        core::Real depth = 30.0;

        mrf_.add_membrane_residue(center, normal, depth, *pose, fullatom_);
        int root = (int) pose->total_residue();

        // For each pose/prefix pair, create an embedding residue, add it to the existing pose
        // chain and then append that chain to the master pose
        for ( int i = 1; i <= nchains; i++ ) {

            // Calculate jump anchor residue from chain center of mass
            core::Size jump = residue_center_of_mass( *pose, pose->conformation().chain_begin(i), pose->conformation().chain_end(i));

            // Grab Resources for the given chain
            core::membrane::properties::SpanningTopologyOP topology = topologies_[i];
            core::membrane::util::EmbedConfigInfoOP def = embeddings_[i];

            // Create a factory and apply residue
            EmbeddingFactoryOP factory = new EmbeddingFactory(pose, def, topology);
            factory->create_and_add_embedding(fullatom_, jump);
            embres_map[ i ] = std::pair< int, int >( jump, pose->total_residue() );
        }

        // Construct and set a new conformation
        MembraneConformationOP mp_conf( new MembraneConformation( pose->conformation(), embres_map, root ) );

        // Initialize Spanning Topology
        initialize_topology( mp_conf );

        // Should initialize multi chain lipds here
        initialize_lips_exp( mp_conf );

        // Set final conf
        pose->set_new_conformation( mp_conf );

        // Done!
        return;
    }

    /// @brief Initialize Spanning Topology
    /// @details Initialize spanning topology in the final pose
    void
    MembraneProteinFactory::initialize_topology( MembraneConformationOP mp_conf  )
    {
        using namespace core::membrane::properties;

        for ( core::Size i = 1; i <= topologies_.size(); i++ )
        {
            // Make a new const ref and add to conf
            SpanningTopology const & topology( *topologies_[i] );
            mp_conf->add_topology_by_chain( topology, i );
        }
    }

    /// @brief Initialize Lipds Exposure Data
    /// @details Initialize lipid exposure data in the final pose
    void
    MembraneProteinFactory::initialize_lips_exp( MembraneConformationOP mp_conf ) {

        using namespace core::membrane::properties;

        if ( include_lips_ ) {

            for ( core::Size i = 1; i <= lipid_acc_.size(); i++ )
            {
                LipidAccInfo const & lipids( *lipid_acc_[i] );
                mp_conf->add_lips_by_chain(lipids, i );
            }
        }
    }

    /// @brief      Initialize Resources from the Resource Manager
    /// @details    Load required resources for initializing a membrane protein
    ///
    /// @note Precondition: Initialized prefix list
    /// @throws EXCN_Resource_Manager (missing reuqired resource)
    void
    MembraneProteinFactory::initialize_resources() {

        using namespace basic::resource_manager;
        using namespace core::membrane::util;

        // Resize maps based on chains_map size
        chains_.resize(chains_map_.size());
        topologies_.resize(chains_map_.size());
        embeddings_.resize(chains_map_.size());
        lipid_acc_.resize(chains_map_.size());

        for ( core::Size i = 1; i <= chains_map_.size(); i++ ) {

            // Prefix Based Resource Tags
            if( ! chains_map_.count(i) ) {
                utility_exit_with_message( "Initialization error: chain not found in MembraneProteinFactory." );
            }
            std::string base_desc = chains_map_[i]; // chains_map_ is non-const
            std::string topo_desc = base_desc + "_span";
            std::string embed_desc = base_desc + "_embed";
            std::string lipid_desc = base_desc + "_lips";

            // Get pose from resource manager
            if ( ! ResourceManager::get_instance()->has_resource_with_description( base_desc ) )
            {
                throw EXCN_Resource_Definition( "Cannot load chain with description " + base_desc );
                TR << "Cannot load any of my resource description things" << std::endl;
            }
            chains_[i] = basic::resource_manager::get_resource< core::pose::Pose >( base_desc );

            // Get topology from resource manager
            if ( ! ResourceManager::get_instance()->has_resource_with_description( topo_desc ) )
            {
                throw EXCN_Resource_Definition( "Cannot load topology with description " + topo_desc );
            }
            topologies_[i] = basic::resource_manager::get_resource< core::membrane::properties::SpanningTopology >( topo_desc );

            // Get embedding from resource manager
            if ( ! ResourceManager::get_instance()->has_resource_with_description( embed_desc) )
            {
                throw EXCN_Resource_Definition( "Cannot load membrane embedding with description " + embed_desc );
            }
            embeddings_[i] = basic::resource_manager::get_resource< core::membrane::util::EmbedConfigInfo >( embed_desc );

            // If user specified to include lips, load lipid acc data
            if ( include_lips_ ) {

                // Get lipids accessibility data from the resource manager
                if ( ! ResourceManager::get_instance()->has_resource_with_description( lipid_desc) )
                {
                    throw EXCN_Resource_Definition( "Cannot load lipid accessibility data with description " + lipid_desc );
                }
                lipid_acc_[i] = basic::resource_manager::get_resource< core::membrane::properties::LipidAccInfo >( lipid_desc );
            }
        }
        return;
    }

} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneProteinFactory_cc


