project( tinyxml2 )
set( TINYXML2_SOURCES tinyxml2.cpp ) # Setup the list of sources here.
add_library( tinyxml2 ${TINYXML2_SOURCES} )
# You can refer to other libraries which should be available from the root cmakelists.
target_link_libraries( tinyxml2 )
