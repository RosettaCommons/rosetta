#
# ui config file tempalte


#
# devel apps, use this to specifying list of pilot apps that is not yet ready for public release but should be nevertheless be always build.
# Each entry should be a tuple (path, project-name)
# where path is the relative path to project file from main/source/ui/apps/pilot
# and project-name is name of .pro file in this path without extension
# for example to add project main/source/ui/apps/pilot/sergey/a/b/c.pro specify ('sergey/a/b', 'c')
#
# Note: if you want to DISABLE building of particualr (or all) devel apps locally - please do this by specifying devel_apps in local `config.py` file
devel_apps = [
    ('awatkins/rna_denovo', 'rna_denovo'),
    ('vmullig/bundle_gui',  'parametric_design'),
    ('sergey/pose_viewer',  'pose_viewer'),
]
