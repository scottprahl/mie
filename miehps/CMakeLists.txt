cmake_minimum_required (VERSION 2.8)
project (MIHPS)

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.

# Create a library called "Hello" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
# add_library (neib neib.cpp)
# add_library (chargeconfig chargeconfig.cpp

include_directories (${MIHPS_SOURCE_DIR})

# Make sure the linker can find the Hello library once it is built.
link_directories (${MIHPS_BINARY_DIR})



set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_Cc_COMPILER cpp)

# set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.

add_executable (timingtestclass miehpstimingtestclass.cpp miehps.cpp)
add_executable (oonewcallhpsmie oonewcallhpsmie.cpp miehps.cpp)
add_executable (oowaterdroplet oowaterdroplet.cpp miehps.cpp)
add_executable (oogold oogold.cpp miehps.cpp)

# Link the executable to the Hello library.
target_link_libraries (timingtestclass gmp gmpxx)
target_link_libraries (oonewcallhpsmie gmp gmpxx)
target_link_libraries (oowaterdroplet gmp gmpxx)
target_link_libraries (oogold gmp gmpxx)
