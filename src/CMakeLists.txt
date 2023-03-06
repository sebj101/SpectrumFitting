# Build executables
add_executable(ToyEndpoint ToyEndpoint.cpp)
target_include_directories(ToyEndpoint PRIVATE ${PROJECT_SOURCE_DIR}/TBetaGenerator)