add_library(penalties penalties.cpp)

target_include_directories(penalties INTERFACE 
                          "${CMAKE_CURRENT_SOURCE_DIR}" 
                          )

target_include_directories(penalties PRIVATE 
                          "${PROJECT_SOURCE_DIR}/include"
                          )

target_link_libraries(penalties PUBLIC tutorial_compiler_flags)

# install(TARGETS penalties DESTINATION lib)
# install(FILES penalties.hpp DESTINATION include)