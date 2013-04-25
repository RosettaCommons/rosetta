add_items(
    SOURCES
    vgrid.c
    vmgrid.c
    vopot.c
    vpmg.c
    vpmgp.c
)

add_items(
    EXTERNAL_HEADERS
    vgrid.h
    vmgrid.h
    vopot.h
    vpmg.h
    vpmgp.h
)

add_sublibrary(mg apbs_generic apbs_pmgc)
