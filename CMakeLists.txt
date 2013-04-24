cmake_minimum_required(VERSION 2.8)
project(merger)

# set(CMAKE_C_COMPILER "gcc-4.7")
# set(CMAKE_CXX_COMPILER "g++-4.7")
# set(CMAKE_EXE_LINKER_FLAGS "-g -pg -Wall")
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

set(CMAKE_CONFIGURATION_TYPES Release)

include_directories(
    src/
    samtools/
    )

add_library(
    bam
    STATIC
	samtools/bgzf.c
    samtools/kstring.c
    samtools/bam_aux.c
    samtools/bam.c
    samtools/bam_import.c
    samtools/sam.c
    samtools/bam_index.c
    samtools/bam_pileup.c
    samtools/bam_lpileup.c
    samtools/bam_md.c
    samtools/razf.c
    samtools/faidx.c
    samtools/bedidx.c
    samtools/knetfile.c
    samtools/bam_sort.c
    samtools/sam_header.c
    samtools/bam_reheader.c
    samtools/kprobaln.c
    samtools/bam_cat.c
    )

set_target_properties(
    bam
    PROPERTIES
    COMPILE_FLAGS "-O3 -Wall"
    )

set_property(
    TARGET
    bam
    PROPERTY 
    COMPILE_DEFINITIONS _FILE_OFFSET_BITS=64 _LARGEFILE64_SOURCE _USE_KNETFILE _CURSES_LIB=1 
    )

add_executable(
    merger
    src/main.cpp
    src/argparse.cpp
    src/merge.cpp
    src/bamfile.cpp
    )

target_link_libraries(merger bam pthread m z)

# do not remove -Wall and -Werror: please fix the errors instead of being lazy
set_target_properties(
    merger
    PROPERTIES
    COMPILE_FLAGS "-O3 -Wall -Werror"
    LINK_FLAGS "-O3 -Wall -Werror"
    # COMPILE_FLAGS "-O3 -g -pg -Wall -Werror"
    )


install(
    TARGETS merger
	RUNTIME DESTINATION bin
	OPTIONAL
    )

find_program(
    ASTYLE
    astyle
    PATHS ENV PATH
    )

if(ASTYLE)
    add_custom_target(
        style
        COMMAND ${ASTYLE} --options=astyle.txt --suffix=none --recursive 'src/*.cpp' 'src/*.hpp' 'src/*.tpp'
        DEPENDS astyle.txt
        )
endif(ASTYLE)
