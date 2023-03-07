# Locate ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Build executables
add_executable(ToyEndpoint ToyEndpoint.cpp)
target_link_libraries(ToyEndpoint PRIVATE ${ROOT_LIBRARIES})
target_include_directories(ToyEndpoint PRIVATE ${PROJECT_SOURCE_DIR}/TBetaGenerator)