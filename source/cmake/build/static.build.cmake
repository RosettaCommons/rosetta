SET( LINK_TYPE STATIC )
### Any external libraries to link statically
if( ${COMPILER} STREQUAL "cl" )
    SET( LINK_EXTERNAL_LIBS zlib ws2_32 )
else()
    SET( LINK_EXTERNAL_LIBS z.a libcppdb-static.a )
endif()
#SET(LINK_STATIC_LIBS boost_regex-gcc41-mt.a boost_thread-gcc41-mt.a GLEW.a ircclient.a openal.a vorbisfile.a vorbis.a ogg.a curl.a ssl.a crypto.a freetype.a png.a lua.a z.a)
