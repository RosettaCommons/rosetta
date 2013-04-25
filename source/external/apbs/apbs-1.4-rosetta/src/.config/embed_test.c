// TODO: figure out what is going on with this

#define EMBED(rctag)

static const char* rctag;
static void* use_rcsid=(0 ? &use_rcsid : (void**)&rcsid);
EMBED(rcsid)
