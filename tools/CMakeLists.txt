SET(TOOLS_SRC
   generateTestDomain
   mergeObj
)

SET(TOOLS_RENDER_HEADERS
    TreeImageRenderer.h
)

SET(TOOLS_RENDER_SRC
    TreeImageRenderer.cpp
)

FOREACH(FILE ${TOOLS_SRC})
  add_executable(${FILE} ${FILE}.cpp ${SOURCE_BASE_FILES} ${SOURCE_BASE_HEADER_FILES} )
  target_link_libraries (${FILE}  ${DGTAL_LIBRARIES} )
ENDFOREACH(FILE)

add_executable(graph2render graph2render.cpp ${SOURCE_BASE_FILES} ${TOOLS_RENDER_SRC} ${SOURCE_BASE_HEADER_FILES} ${TOOLS_RENDER_SRC})
target_link_libraries(graph2render ${DGTAL_LIBRARIES})
