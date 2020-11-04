set(PROJECT_SRC_DIR ${SRC}/${PROJECT_NAME})

#collect everything within the src folder
file(GLOB PROJECT_SRC_FILES 
    ${PROJECT_SRC_DIR}/*.cpp 
    ${PROJECT_SRC_DIR}/*.h
    ${PROJECT_SRC_DIR}/*.cu
    #${PROJECT_SRC_DIR}/*/*.cpp 
    #${PROJECT_SRC_DIR}/*/*.h
    #${PROJECT_SRC_DIR}/*/*.cu
)

set(WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
configure_file(workingDirectory.h.in workingDirectory.h)

# static libs
set(STOPWATCH_INCLUDE ${SRC}/Stopwatch)
add_library(Stopwatch STATIC 
    ${SRC}/Stopwatch/Stopwatch.cpp 
    ${STOPWATCH_INCLUDE}/Stopwatch.h
)
target_compile_definitions(Stopwatch PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

set(BLOSSOM_INCLUDE ${SRC}/Blossom5)
file(GLOB BLOSSOM_FILES
	${BLOSSOM_INCLUDE}/*.cpp
	${BLOSSOM_INCLUDE}/*.h
	${BLOSSOM_INCLUDE}/*/*.cpp
	${BLOSSOM_INCLUDE}/*/*.h
)
add_library(Blossom5 STATIC 
    ${BLOSSOM_FILES}
)
target_compile_definitions(Blossom5 PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

add_subdirectory(BackgroundGrid)

#add project executable
add_executable(${PROJECT_NAME} 
    main.cpp 
    ${PROJECT_SRC_FILES}
)

target_include_directories(${PROJECT_NAME} PRIVATE 
    #${SRC}
    ${PROJECT_SRC_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ### externals ###
    ${OPENMESH_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${LIBS}
    Stopwatch
	Blossom5
	BackgroundGrid
)

target_compile_definitions(${PROJECT_NAME} PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
endif()