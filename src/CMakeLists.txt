# Locate ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Generate library from source files
add_library(Spectrum)
target_sources(Spectrum PRIVATE Spectrum.cpp)


# Build executables
add_executable(ToyEndpoint ToyEndpoint.cpp)
target_link_libraries(ToyEndpoint PRIVATE Spectrum ${ROOT_LIBRARIES})
target_include_directories(ToyEndpoint PRIVATE ${PROJECT_SOURCE_DIR}/TBetaGenerator)