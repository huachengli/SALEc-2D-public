cmake_minimum_required(VERSION 3.10)
project(sale_rebuild C)

set(CMAKE_C_STANDARD 99)

find_package(MPI REQUIRED)
if(MPI_FOUND)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
    message(ERROR "cannot find MPI")
endif()


set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -fopenmp)

set(PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
include_directories(${PROJECT_SOURCE_DIR}/mpi_vec)
include_directories(${PROJECT_SOURCE_DIR}/io)
include_directories(${PROJECT_SOURCE_DIR}/geometry)
include_directories(${PROJECT_SOURCE_DIR}/eos)
include_directories(${PROJECT_SOURCE_DIR})

link_libraries(m z)

set(SRCS src/mpi_vec/mpi_vec.c src/mpi_vec/mpi_vec_ght.c
        src/mpi_vec/linked_list.c
        src/mpi_vec/index_var.c
        src/io/InputParser.c src/io/output.c
        src/geometry/interpolate_fn.c
        src/eos/eos_state.c src/eos/aneos_sieden.c
        src/geometry/interpolate.c src/geometry/linear.c
        src/sale2d.c src/sale2d_back.c
        src/advect.c
        src/strength.c
        src/tracer.c
        src/ejecta.c
        src/write_vtk.c
        src/target_profile.c
        src/construct.c
        src/velocity.c
        src/geometry/interpolate_2d9.c)

add_executable(${PROJECT_NAME} src/salec2d.c ${SRCS} src/io/iSALEInputParser.c)
add_executable(test_till src/other/test_functions.c ${SRCS})
add_executable(merge_vts src/other/merge_vts.c ${SRCS})
add_executable(generate_input src/io/generate_input.c src/io/iSALEInputParser.c
        ${SRCS})
add_dependencies(${PROJECT_NAME} merge_vts test_till)
target_compile_options(${PROJECT_NAME} PRIVATE
        -O3
#        -DDEBUG_MODE=1
#        -DDEBUG_RANK=1
        -DOUTPUT_GHT
#        -DDEBUG_RANK_ONLY
#        -save-temps=obj
#        -DCOMM_OFF
#        -DCHECK_NAN
        )
target_compile_options(test_till PRIVATE
        -DZLIB_DEBUG_MODE=1
        )
target_compile_options(merge_vts PRIVATE
        -O3)
target_compile_options(generate_input PRIVATE
        -O0)

add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME} ${CMAKE_BINARY_DIR}/salec2d
)