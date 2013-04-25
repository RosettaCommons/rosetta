#! /bin/bash

# Regenerate each libraries' project file from its scons build definition
echo "Regenerating project files..."
for project in basic core.1 core.2 core.3 core.4 core.5 numeric ObjexxFCL utility protocols.1 protocols_a.2 protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4 protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4 protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5 protocols_f.5 protocols.6 protocols.7
do
	python makeFiles.py $project
done

# Remove all references to Visual Studio 2010+ project files from the solution.
# Unfortunately, this also removes project dependencies required to set build order.
echo "Removing references to vcxproj..."
sed -i 's/vcxproj/vcproj/g' boinc_mini.sln

# Rename executables and program databases
echo "Renaming executables..."
sed -i 's/minirosetta_beta_3.25_windows_intelx86.exe/minirosetta.exe/g' apps/minirosetta.vcproj
sed -i 's/minirosetta_graphics_3.25_windows_intelx86.exe/minirosetta_graphics.exe/g' apps/minirosetta_graphics.vcproj
sed -i 's/minirosetta_beta_3.25_windows_intelx86.pdb/minirosetta.pdb/g' apps/minirosetta.vcproj
sed -i 's/minirosetta_graphics_3.25_windows_intelx86.pdb/minirosetta_graphics.pdb/g' apps/minirosetta_graphics.vcproj

# Include additional libraries
echo "Including libraries..."
sed -i 's/basic.lib/basic.lib sqlite3.lib cppdb.lib/g' apps/minirosetta.vcproj
sed -i 's/basic.lib/basic.lib sqlite3.lib cppdb.lib/g' apps/minirosetta_graphics.vcproj

echo "Removing upgrade logs..."
rm -rf _UpgradeReport_Files UpgradeLog*.XML

echo "Done!"
