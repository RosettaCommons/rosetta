Several standards of formatting and language are employed within these sample
scripts to help you better understand them. Please consult the brief summary
of PyRosetta syntax and conventions found online if you have any questions.

Rosetta objects are typically referred with their proper object name (Mover,
FoldTree, MoveMap, etc.) except when used in abstract discussion or when
referencing a specific instance. This is usually reserved for Pose which
is usually an instance and thus "pose".

Comments, lists, or other instructions are indented (usually 4 spaces) if
    they continue onto successive lines.

Lists are typically setup with numbers, letters, dashes, and whatever else so:
    1. read this document
    2. investigate a sample script
        a. read its basic summary
        b. read the code in-depth
            -understand what each object does
        c. run the protocol
            -put in proper arguments
            -analyze the results
                >look at output files
                >view any output protein in PyMOL (or you favorite viewer)

"""
Comment blocks, sections of code surrounded by triple quotes (""") indicate
significant sections of text or a method definition. A comment block outside
of a method describes something general or addresses a slightly different topic.
"""

# individual comments should be brief
# sometimes, a comment gets too long and
#    continues onto the next line
# comments should match the step specified in the protocol outline

#### four comment symbols indicate a comment or region of the code that
####    is likely to change if you adjust the protocol

A necessary file will be referenced by its filename, a string. All filename
variables will end in "_filename" such as the common "pdb_filename". Nearly
all methods used which require information be read in from a file will
make use of a filename variable.

List comprehensions are usually written on multiple lines with a single line
explaining the output, what to loop over, or any conditional statements applied.

The word "create" implies a single object is created, constructed, initialized,
    or whatever you want to call it or the object is created and modified
The word "load" implies some additional file is required and some of its
    information is loaded into an object or into the script
The phrase "sets up" or "setup" implies a single process that may involve
    multiple objects, for example, setting up a PackRotamersMover requires
    a PackerTask, which requires a Pose, but instead of saying all this,
    just say "sets up PackRotamersMover"
The word "perform" implies a number of individual steps will occur
The word "obtain" implies that data is extracted from another object, often
    into a form which is easier to access or compatible with other objects

Each sample script here will (hopefully) contain
1.  A general explanation of the protocol performance and requirements
        a. summary
        b. instructions
        c. references
2.  The basic protocol overview outlining individual processes which
        are performed in the script
3.  The code itself written as a single method with several input arguments
4.  A guideline for interpreting results for the protocol
5.  Commandline compatibility code defining all appropriate options and then
        calling the single method (from 3.) with these values
6.  A section discussing alternate scenarios for similar code, this will
        hopefully include:
        a.  how to obtain any necessary "non-standard" files
        b.  a real example problem and appropriate parameter values
        c.  alternate sampling options
        d.  alternate scoring options
        e.  any other common changes or notable alternatives

The scripts:
    dna_interface.py
    docking.py
    folding.py
    ligand_interface.py
    refinement.py
all output a large number of structures to PyMOL (using the PyMOL_Mover) by
default. If these structures are not appearing, ensure that you are supplying
the proper IP address using the "--PyMOL_Mover_ip" option (or edit the sample
script to your liking and that your computer is not blocking the packets sent
by the PyMOL_Mover. to run these scripts without this output, specify
"PyMOL_Mover_ip off". To view an intuitive representation of the protocol,
I recommend the following commands in PyMOL:
    1. execute the command "remove hydro" to remove hydrogens
    2. execute the command "show cartoon"
    3. clock the "play" button on the bottom right of the viewer GUI

