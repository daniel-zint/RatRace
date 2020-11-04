set(BACKGROUNDGRID_INCLUDE ${SRC}/BackgroundGrid)
add_library(BackgroundGrid STATIC 
    ${BACKGROUNDGRID_INCLUDE}/BackgroundGrid.cpp 
    ${BACKGROUNDGRID_INCLUDE}/BackgroundGrid.h
    ${BACKGROUNDGRID_INCLUDE}/SizeGrid.cpp 
    ${BACKGROUNDGRID_INCLUDE}/SizeGrid.h
)
target_include_directories(BackgroundGrid PRIVATE 
    ${OPENMESH_INCLUDE_DIR}
)

source_group("BackgroundGrid" FILES 
    ${BACKGROUNDGRID_INCLUDE}/BackgroundGrid.cpp 
    ${BACKGROUNDGRID_INCLUDE}/BackgroundGrid.h
)
source_group("SizeGrid" FILES 
    ${BACKGROUNDGRID_INCLUDE}/SizeGrid.cpp 
    ${BACKGROUNDGRID_INCLUDE}/SizeGrid.h
)

set(BACKGROUNDGRID_LIBS glog::glog)
if(NOT MSVC)
	set(BACKGROUNDGRID_LIBS ${BACKGROUNDGRID_LIBS} X11)
endif()
target_link_libraries(BackgroundGrid ${BACKGROUNDGRID_LIBS})
target_compile_definitions(BackgroundGrid PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)
