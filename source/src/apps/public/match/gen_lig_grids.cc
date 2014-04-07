// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file gen_lig_grids.cc
/// @brief generate matcher grids for a protein with a ligand
/// @author Yuan Liu
/// @detail flags example
/// -ignore_unrecognized_res
/// -grid_delta 0.5
/// -grid_lig_cutoff 4.0
/// -grid_bb_cutoff 2.25
/// -grid_active_res_cutoff 5.0

//std
#include <iostream>
#include <map>
#include <string>

#include <devel/init.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/io/pdb/file_data.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>

// AUTO-REMOVED #include <core/scoring/packstat/compute_sasa.hh>


// AUTO-REMOVED #include <utility/string_util.hh>
#include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>

//from comput_sasa.cc
// AUTO-REMOVED #include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <ObjexxFCL/format.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/packstat/types.hh>
#include <utility/vector1.hh>
#include <fstream>

#include <utility/excn/Exceptions.hh>

using namespace std;
using namespace core;
using namespace core::io::pdb;
using namespace core::conformation;
using namespace core::chemical;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring::packstat;

using namespace ObjexxFCL::format;
using namespace numeric;
using namespace utility;

OPT_KEY(Real, grid_delta)
OPT_KEY(Real, grid_lig_cutoff)
OPT_KEY(Real, grid_bb_cutoff)
OPT_KEY(Real, grid_active_res_cutoff)

basic::Tracer TR("apps.public.match.gen_lig_grids");

inline void assure(std::ifstream& in, const char* filename = "")
{
    using namespace std;
    if (!in)
    {
        TR << "Could not open file!" << endl;
        utility_exit_with_message(filename);
    }
}

inline void assure(std::ofstream& in, const char* filename = "")
{
    using namespace std;
    if (!in)
    {
        TR << "Could not open file!" << endl;
        utility_exit_with_message(filename);
    }
}

int main( int argc, char * argv [] )
{
    try {
	//normal init
    NEW_OPT(grid_delta, "Size of grids", 0.5);
    NEW_OPT(grid_lig_cutoff, "Grid range around ligand", 4.0);
    NEW_OPT(grid_bb_cutoff, "Backbone occupation size", 2.25);
    NEW_OPT(grid_active_res_cutoff, "Active residues around ligand", 5.0);

    devel::init( argc, argv );

    pose::Pose scaffold;

    if ( option[ in::file::s ].user() )
    {
        core::import_pose::pose_from_pdb( scaffold, option[ in::file::s ]()[1] );
    }
    else
    {
        TR.Error << "User did not specify the pdb file!" << endl;
        exit( EXIT_FAILURE );
    }

    //get the atom spheres and center of the residues
    //PosePackData pd_scaffold = pose_to_pack_data(scaffold);

    //get information from ligand
    string proname=option[ in::file::s ]()[1];
    string ligname=option[ in::file::s ]()[2];
    TR << "Ligand pdb file: " << ligname << endl;
    ifstream ligfp(ligname.c_str());
    assure(ligfp, ligname.c_str());

    string line;
    utility::vector1< xyzVector<core::Real> > ligxyz;
    core::Real maxligx=-9999, maxligy=-9999, maxligz=-9999;
    core::Real minligx=9999, minligy=9999, minligz=9999;

    //pick first residue
    int lastresndx=-9999;
    bool boxflag=true;

    while ( getline(ligfp, line) )
    {
		if (line.substr(0,6)!="HETATM") continue;
        core::Real x = atof(line.substr(30,8).c_str());
        core::Real y = atof(line.substr(38,8).c_str());
        core::Real z = atof(line.substr(46,8).c_str());
        int n = atoi(line.substr(22,4).c_str());
        //if (lastresndx>0 && n!=lastresndx) boxflag=false;
        //only take the first residue, in case there are many
        //what if it is a peptide? change the ligand pdb?
        if (lastresndx>0 && n!=lastresndx) break;
        lastresndx = n;
        //debug
        //TR << x << " " << y << " " << z << endl;
        if (boxflag)
        {
            minligx = minligx<x?minligx:x;
            minligy = minligy<y?minligy:y;
            minligz = minligz<z?minligz:z;
            maxligx = maxligx>x?maxligx:x;
            maxligy = maxligy>y?maxligy:y;
            maxligz = maxligz>z?maxligz:z;
        }
        xyzVector<core::Real> xyz(x,y,z);
        ligxyz.push_back(xyz);
    }
    TR << "Ligand atom number: " << ligxyz.size() << endl;
    if (ligxyz.size()<1)
    {
        utility_exit_with_message("No ligand atom found!");
    }

    //define the box
    const core::Real dx=option[grid_delta];
    const core::Real dy=option[grid_delta];
    const core::Real dz=option[grid_delta];             //width of the bin
    const core::Real cutoff = option[grid_lig_cutoff];  //cutoff: ligand atom
    const core::Real bbcutoff = option[grid_bb_cutoff]; //cutoff for bb atom

    const core::Real cenx = (maxligx+minligx)/2.0;
    const core::Real ceny = (maxligy+minligy)/2.0;
    const core::Real cenz = (maxligz+minligz)/2.0;
    const core::Real wdx = maxligx-minligx+cutoff*2.0;
    const core::Real wdy = maxligy-minligy+cutoff*2.0;
    const core::Real wdz = maxligz-minligz+cutoff*2.0;
    const core::Real basex=cenx-wdx/2.0;
    const core::Real basey=ceny-wdy/2.0;
    const core::Real basez=cenz-wdz/2.0;
    const core::Size nx = int(wdx/dx);
    const core::Size ny = int(wdy/dy);
    const core::Size nz = int(wdz/dz);

    //bbbox
    //const core::Real maxbbx = (wdx>15.0)?(wdx+5.0):20.0;
    //const core::Real maxbby = (wdy>15.0)?(wdy+5.0):20.0;
    //const core::Real maxbbz = (wdz>15.0)?(wdz+5.0):20.0;	//range of the box for bb grid

    //create the active res list
    utility::vector1< core::Size > active_res_ndx;
    core::Real threshold = option[grid_active_res_cutoff];
    //for each residue
    for ( core::Size j=1; j<=scaffold.n_residue(); j++)
    {
        core::conformation::ResidueCOP res( &scaffold.residue(j) );
        bool flag=true;

        //for each atom
        for ( core::Size k=1; k<=res->atoms().size() && flag; k++)
        {
            Vector const & xyzatom = res->atoms()[k].xyz();
            //for each cav ball
            for ( core::Size n=1; n<=ligxyz.size() && flag; n++ )
            {
                Vector const d = xyzatom - ligxyz[n];
                //cutoff 5.0A from will
                if ( d.length_squared() < threshold*threshold) flag=false;
            }
        }

        //save the neighbor number of this residue if it is the active res
        if (!flag)
        {
            active_res_ndx.push_back(j);
        }
    }

    // output the pos file
    //string posname(proname+".pos");
    string posname(proname+"_0.pos");
    ofstream out_pos(posname.c_str());
    assure(out_pos,posname.c_str());
    for ( core::Size i=1; i<=active_res_ndx.size(); i++ )
    {
        out_pos << active_res_ndx[i] << " ";
    }
    out_pos.close();

    // output the gridlig file
    //string gridname(proname+".gridlig");
    string gridname(proname+"_0.gridlig");
    ofstream out_grid(gridname.c_str());
    assure(out_grid,gridname.c_str());

    //Title
    out_grid << "NAME: gridlig" << endl;
    out_grid << "BASE: " << basex << " " << basey << " " << basez << endl;
    out_grid << "SIZE: " << nx << " "<< ny << " "<< nz << endl;
    out_grid << "LENGTH: " << dx << " " << dy << " " << dz << endl;

    for (core::Size i=1; i<=nx; i++)
    {
        for (core::Size j=1; j<=ny; j++)
        {
            for (core::Size k=1; k<=nz; k++)
            {
                //the grid's xyz
                core::Real xx = basex+dx*i;
                core::Real yy = basey+dy*j;
                core::Real zz = basez+dz*k;
                xyzVector<core::Real> xyz(xx,yy,zz);

                //test if this grid is empty
                bool flag=false;

                //lig occupy
                //check all ligxyz
                for ( core::Size m = 1; m <= ligxyz.size() && !flag; m++ )
                {

                    xyzVector<core::Real> d = xyz - ligxyz[m];
                    if (d.length_squared()<cutoff*cutoff) flag=true;

                }

                //remove res occupancy
                for ( core::Size m=1; m <= active_res_ndx.size() && flag; m++ )
                {
                    //for each active res
                    core::conformation::ResidueCOP res( &scaffold.residue(active_res_ndx[m]) );
                    //for each atom
                    for ( core::Size n=1; n <= res->atoms().size() && flag; n++ )
                    {
                        //backbone
                        if (res->atom_is_backbone(n))
                        {
                            xyzVector<core::Real> d = xyz - res->atoms()[n].xyz();
                            if ( d.length_squared() < (bbcutoff*bbcutoff) ) flag=false;
                        }
                    }
                }

                if (flag) out_grid<<"1";
                else out_grid<<"0";
                out_grid << " ";
            }
            out_grid << endl;
        }
        out_grid << endl;
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return EXIT_SUCCESS;
}

