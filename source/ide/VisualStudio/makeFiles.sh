#! /bin/bash

# Regenerate each libraries' project file from its scons build definition
echo "Regenerating project files..."
for project in basic devel core.1 core.2 core.3 core.4 core.5 numeric ObjexxFCL utility protocols.1 protocols.3 protocols.4 protocols.7 protocols.8 protocols_a.2 protocols_a.5 protocols_a.6 protocols_b.2 protocols_b.5 protocols_b.6 protocols_c.5 protocols_c.6 protocols_d.5 protocols_d.6 protocols_e.5 protocols_e.6 protocols_f.5 protocols_g.5 protocols_h.5 
do
	python makeFiles.py $project
done

echo "Done!"
