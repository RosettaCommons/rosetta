#! /bin/bash

# Regenerate each libraries' project file from its scons build definition
echo "Regenerating project files..."
for project in basic core.1 core.2 core.3 core.4 core.5 numeric ObjexxFCL utility protocols.1 protocols_a.2 protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4 protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4 protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5 protocols_f.5 protocols.6 protocols.7
do
	python makeFiles.py $project
done

echo "Done!"
