find_program(DOXYGEN_PATH doxygen REQUIRED)

if(DOXYGEN_PATH_NOTFOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation. Please install it on your system")
else()
    message(STATUS "Doxygen found.")

    # Download the Doxygen Awesome CSS files
    file(DOWNLOAD https://raw.githubusercontent.com/jothepro/doxygen-awesome-css/v2.3.1/doxygen-awesome.css
        ${CMAKE_CURRENT_LIST_DIR}/doxygen-awesome.css)
    file(DOWNLOAD https://raw.githubusercontent.com/jothepro/doxygen-awesome-css/v2.3.1/doxygen-awesome-sidebar-only.css
        ${CMAKE_CURRENT_LIST_DIR}/doxygen-awesome-sidebar-only.css)

    # Define the output file generated by Doxygen
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_BINARY_DIR}/html)
    set(DOXYGEN_OUTPUT_FILE ${DOXYGEN_OUTPUT_DIR}/index.html)

    # Add a custom command to build the documentation
    add_custom_command(
        OUTPUT ${DOXYGEN_OUTPUT_FILE}
        COMMAND ${DOXYGEN_PATH} ${CMAKE_CURRENT_LIST_DIR}/Doxyfile
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating Doxygen documentation..."
        VERBATIM
    )

    # Ensure the documentation is built as part of the default build process
    add_custom_target(generate_docs ALL
        DEPENDS ${DOXYGEN_OUTPUT_FILE}
    )

    # Optional message after build
    add_custom_command(
        TARGET generate_docs POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E echo "Documentation successfully generated. You can preview it at: ${DOXYGEN_OUTPUT_FILE}"
    )
endif()
