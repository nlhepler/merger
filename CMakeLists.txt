cmake_minimum_required(VERSION 2.8)

# set(CMAKE_C_COMPILER "gcc-4.7")
# set(CMAKE_CXX_COMPILER "g++-4.7")
# set(CMAKE_EXE_LINKER_FLAGS "-g -pg -Wall")

project(merger)

set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

set(CMAKE_BUILD_TYPE DEBUG)
# set(CMAKE_BUILD_TYPE RELEASE)

set(
    COMMON_CXX_FLAGS
    "-Wall -Werror -Wno-unknown-pragmas -fstack-protector-all -fpie"
    )

set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g -pg ${COMMON_CXX_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "-O3 -g -pg -Wall")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${COMMON_CXX_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wall")

add_library(
    bam
    STATIC
	external/samtools/bgzf.c
    external/samtools/kstring.c
    external/samtools/bam_aux.c
    external/samtools/bam.c
    external/samtools/bam_import.c
    external/samtools/sam.c
    external/samtools/bam_index.c
    external/samtools/bam_pileup.c
    external/samtools/bam_lpileup.c
    external/samtools/bam_md.c
    external/samtools/razf.c
    external/samtools/faidx.c
    external/samtools/bedidx.c
    external/samtools/knetfile.c
    external/samtools/bam_sort.c
    external/samtools/sam_header.c
    external/samtools/bam_reheader.c
    external/samtools/kprobaln.c
    external/samtools/bam_cat.c
    )

set_property(
    TARGET bam
    PROPERTY COMPILE_DEFINITIONS _FILE_OFFSET_BITS=64 _LARGEFILE64_SOURCE _USE_KNETFILE _CURSES_LIB=1
    )

add_executable(
    merger
    src/merger/main.cpp
    src/merger/args.cpp
    src/aligned.cpp
    src/bamfile.cpp
    src/merge.cpp
    src/util.cpp
    )

target_link_libraries(merger bam pthread m z)

# do not remove -Wall and -Werror: please fix the errors instead of being lazy
set_property(
    TARGET merger
    APPEND PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src/merger
    )

add_executable(
    puncher
    src/puncher/main.cpp
    src/puncher/args.cpp
    src/aligned.cpp
    src/bamfile.cpp
    src/coverage.cpp
    src/rateclass.cpp
    src/util.cpp
    )

target_link_libraries(puncher bam pthread m z)

set_property(
    TARGET puncher
    APPEND PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src/puncher
    )

# add_executable(
#     sampler
#     src/sampler/main.cpp
#     src/sampler/args.cpp
#     src/aligned.cpp
#     src/bamfile.cpp
#     src/merge.cpp
#     src/util.cpp
#     )
#
# target_link_libraries(sampler bam pthread m z)
#
# set_property(
#     TARGET sampler
#     APPEND PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src/sampler
#     )

add_executable(
    variants
    src/variants/main.cpp
    src/variants/args.cpp
    src/aligned.cpp
    src/bamfile.cpp
    src/coverage.cpp
    src/rateclass.cpp
    src/util.cpp
    )

target_link_libraries(variants bam pthread m z)

set_property(
    TARGET variants
    APPEND PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src/variants
    )

set_property(
    TARGET merger puncher variants # sampler variants
    APPEND PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

add_executable(
    binmix
    src/binmix/main.cpp
    src/rateclass.cpp
    )

target_link_libraries(binmix m)

set_property(
    TARGET binmix
    APPEND PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

install(
    TARGETS merger puncher variants # sampler variants
	RUNTIME DESTINATION bin
	OPTIONAL
    )

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUXX)
    find_package(OpenMP)

    if (${OPENMP_FOUND})
        set_property(
            TARGET binmix merger puncher variants
            APPEND PROPERTY COMPILE_FLAGS "${OpenMP_CXX_FLAGS}"
        )
        set_property(
            TARGET binmix merger puncher variants
            APPEND PROPERTY LINK_FLAGS "${OpenMP_CXX_FLAGS}"
            )
    endif (${OPENMP_FOUND})
endif (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUXX)
