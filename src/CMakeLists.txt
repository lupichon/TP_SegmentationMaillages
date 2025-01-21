add_subdirectory(geomAlgoLib)

add_executable(mainApp "main.cpp")

target_link_libraries(mainApp PRIVATE geomAlgoLib)
