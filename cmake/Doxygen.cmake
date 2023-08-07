find_package(Doxygen
        REQUIRED dot)

set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_HTML_OUTPUT
        ${PROJECT_SOURCE_DIR}/RFL_docs)
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE
        "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

# Settings for Inheritance Diagram etc.
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_CLASS_DIAGRAMS YES)
set(DOXYGEN_HIDE_UNDOC_RELATIONS NO)
set(DOXYGEN_CALL_GRAPH YES)
set(DOXYGEN_CALLER_GRAPH YES)
set(DOXYGEN_HAVE_DOT             YES)
set(DOXYGEN_CLASS_GRAPH          YES)
set(DOXYGEN_COLLABORATION_GRAPH  NO)
set(DOXYGEN_UML_LOOK             YES)
set(DOXYGEN_UML_LIMIT_NUM_FIELDS 50)
set(DOXYGEN_TEMPLATE_RELATIONS   YES)
set(DOXYGEN_DOT_GRAPH_MAX_NODES  100)
set(DOXYGEN_MAX_DOT_GRAPH_DEPTH  0)
set(DOXYGEN_DOT_TRANSPARENT      YES)
set(DOXYGEN_INTERACTIVE_SVG  YES)
set(DOXYGEN_USE_MATHJAX YES)
set(DOXYGEN_SOURCE_BROWSER YES)

# Exclude tests from documentation
set(DOXYGEN_EXCLUDE_PATTERNS */tests/*)


function(Doxygen input output)
    if (NOT DOXYGEN_FOUND)
        add_custom_target(doxygen COMMAND false
                COMMENT "Doxygen not found")
        return()
    endif()
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_HTML_OUTPUT
            ${PROJECT_BINARY_DIR}/${output})


    doxygen_add_docs(${output}
            ${input}
            COMMENT "Generate HTML documentation"
            )
endfunction()
