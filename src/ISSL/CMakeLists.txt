cmake_minimum_required(VERSION 3.8)

project(ISSL_offtarget_scoring VERSION 1.0)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams)

add_library(tutorial_compiler_flags INTERFACE)
target_compile_features(tutorial_compiler_flags INTERFACE cxx_std_17)

add_subdirectory(penalties)

add_executable(ISSLCreateIndex ISSLCreateIndex.cpp)
target_link_libraries(ISSLCreateIndex PUBLIC 
                      tutorial_compiler_flags 
                      Boost::boost 
                      Boost::iostreams
                      )

target_include_directories(ISSLCreateIndex PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/include"
                          )

add_executable(ISSLScoreOfftargets ISSLScoreOfftargets.cpp)
target_link_libraries(ISSLScoreOfftargets PUBLIC 
                      tutorial_compiler_flags 
                      OpenMP::OpenMP_CXX
                      Boost::boost 
                      Boost::iostreams
                      penalties
                      )

target_include_directories(ISSLScoreOfftargets PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/include"
                          )

# NOTE: This version was created to simply report off-targets for results section
#       REMOVE later                          
add_executable(ISSLGetOfftargets ISSLGetOfftargets.cpp)
target_link_libraries(ISSLGetOfftargets PUBLIC 
                     tutorial_compiler_flags 
                     OpenMP::OpenMP_CXX
                     Boost::boost 
                     Boost::iostreams
                     penalties
                     )

target_include_directories(ISSLGetOfftargets PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           "${PROJECT_SOURCE_DIR}/include"
                           )

# install(TARGETS tutorial_compiler_flags DESTINATION lib)

set(include_dirs 
   "${PROJECT_SOURCE_DIR}/include/fmt"
   "${PROJECT_SOURCE_DIR}/include/libpopcnt"
   "${PROJECT_SOURCE_DIR}/include/phmap"
   )
# install(DIRECTORY ${include_dirs} DESTINATION include)

set(installable_bins ISSLCreateIndex ISSLScoreOfftargets ISSLGetOfftargets)
install(TARGETS ${installable_bins} DESTINATION .)
