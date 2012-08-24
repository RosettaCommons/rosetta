this directory is intended to house PyRosetta helper scripts etc.

Unfortunately, these seem to be quite diverse and potentially depend on
a great deal of outside software

Evan will be adding a bunch of method here soon but the plan with be
to have the following dubsirectories

download_tutorials
environments
methods
pymol

I will be adding methods here soon/move the existing tools here

more explicitely

download_tutorials
    contains plain text and .pdf files for downloading 3rd party software
    ex. download_pymol.pdf , download_psipred.pdf

environments
    simple python "scripts" which load all desired methods for a unique
    interpreter session (an "environment") or script environment
    I have found these extremely useful since it allows you to organize
    methods in the filesystem  however you like while not convoluting
    imports and preventing repetetive typing
    ex. start_pyrosetta.py will start pyrosetta and additionally load
        any extra toolbox methods you want
    start_bioinformatics_environment.py will load non-pyrosetta methods in
        toolbox which may be useful for fast investigations

methods
    the bulk tool package methods, these will rely on different applications
    but CAN all be combined into a single environment
    ...I would love this to be the raw toolbox but its a bit confusing with
    the other non-pythonic helper directories

pymol
    associated PyMOL scripts which (typically) cannot be loaded into a
    ipython session (since the GUI needs to be active)


