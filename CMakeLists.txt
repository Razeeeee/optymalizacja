cmake_minimum_required(VERSION 3.16)
project(Optymalizacja)

# Ustaw standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dodaj flagi kompilacji
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Jeśli nie określono typu budowania, ustaw na Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Znajdź wszystkie pliki źródłowe
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "src/*.h")

# Stwórz wykonywany plik
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Dodaj katalog z nagłówkami
target_include_directories(${PROJECT_NAME} PRIVATE src)

# Ustaw ścieżkę wyjściową dla pliku wykonywalnego
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Wyświetl informacje o konfiguracji
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Source files: ${SOURCES}")
