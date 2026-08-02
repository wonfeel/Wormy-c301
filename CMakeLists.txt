cmake_minimum_required(VERSION 3.18)
project(TesseraWorm VERSION 1.0 LANGUAGES CXX C)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ===== Движок - подтягивается по URL/тегу, не как ручной клон =====
# GIT_TAG указывает на ветку по умолчанию; закрепите на конкретный тег/коммит,
# как только в TesseraEngine появится первый релиз - плавающая ветка удобна на
# старте, но не для воспроизводимой сборки.
include(FetchContent)
FetchContent_Declare(
    TesseraEngine
    GIT_REPOSITORY https://github.com/wonfeel/TesseraEngine.git
    GIT_TAG main
)
FetchContent_MakeAvailable(TesseraEngine)

# ===== Копирование ассетов =====
add_custom_target(copy_assets ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${tesseraengine_SOURCE_DIR}/Shaders" "${CMAKE_BINARY_DIR}/Shaders"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/Shaders" "${CMAKE_BINARY_DIR}/Shaders"
    COMMENT "Copying engine + worm shaders"
)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/demo/worm/data")
    add_custom_command(TARGET copy_assets POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/demo/worm/data" "${CMAKE_BINARY_DIR}/worm_data"
        COMMENT "Copying worm connectome data"
    )
endif()

# Собрать одно приложение из main.cpp (+ соседние .cpp/.c в той же папке).
function(fe_add_app MAIN_CPP TARGET_NAME)
    get_filename_component(APP_DIR "${MAIN_CPP}" DIRECTORY)
    file(GLOB_RECURSE APP_EXTRA "${APP_DIR}/*.cpp" "${APP_DIR}/*.c")
    list(REMOVE_ITEM APP_EXTRA "${MAIN_CPP}")
    add_executable(${TARGET_NAME} "${MAIN_CPP}" ${APP_EXTRA})
    add_dependencies(${TARGET_NAME} copy_assets)
    target_link_libraries(${TARGET_NAME} PRIVATE EngineLib)
    target_include_directories(${TARGET_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${tesseraengine_SOURCE_DIR}
        ${TESSERA_ENGINE_GLFW_INCLUDE_DIR}
        ${TESSERA_ENGINE_GLAD_INCLUDE_DIR}
        ${TESSERA_ENGINE_GLM_INCLUDE_DIR}
    )
    if(TESSERA_ENGINE_IMGUI_AVAILABLE)
        target_include_directories(${TARGET_NAME} PRIVATE ${TESSERA_ENGINE_IMGUI_DIR})
    endif()
    message(STATUS "Added ${TARGET_NAME} (${MAIN_CPP})")
endfunction()

# Обойти дерево, собрать приложения из каждого main.cpp. Имя цели — из пути
# относительно базовой папки: demo/worm -> Demo_worm, tests/worm_regression -> Test_worm_regression.
function(fe_discover_apps BASE_DIR PREFIX)
    file(GLOB_RECURSE MAINS "${BASE_DIR}/*.cpp")
    foreach(MAIN ${MAINS})
        get_filename_component(FN "${MAIN}" NAME)
        if(FN STREQUAL "main.cpp")
            get_filename_component(MDIR "${MAIN}" DIRECTORY)
            file(RELATIVE_PATH REL "${BASE_DIR}" "${MDIR}")
            string(REPLACE "/" "_" NAME "${REL}")
            fe_add_app("${MAIN}" "${PREFIX}_${NAME}")
        endif()
    endforeach()
endfunction()

fe_discover_apps("${CMAKE_CURRENT_SOURCE_DIR}/demo"  "Demo")
fe_discover_apps("${CMAKE_CURRENT_SOURCE_DIR}/tests" "Test")

# ===== CTest-регистрация (только настоящие, argument-free, PASS/FAIL тесты) =====
# fe_discover_apps выше собирает ВСЕ tests/*/main.cpp как исполняемые цели,
# включая десятки одноразовых калибровочных харнессов (свои CLI-аргументы,
# свой формат вывода, не задуманы как pass/fail regression) - регистрировать
# их ВСЕ в ctest бессмысленно и вводит в заблуждение. Регистрируем явно только
# то, что реально устроено как regression-тест: без аргументов, чёткий
# PASS/FAIL, ненулевой exit code при провале.
if(TARGET Test_worm_locomotion)
    add_test(NAME Test_worm_locomotion COMMAND Test_worm_locomotion)
endif()
if(TARGET Test_worm_regression)
    add_test(NAME Test_worm_regression COMMAND Test_worm_regression)
endif()
