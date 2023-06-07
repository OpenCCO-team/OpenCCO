SET(TOOLS_SRC
   generateTestDomain
   mergeObj
)

FOREACH(FILE ${TOOLS_SRC})
  add_executable(${FILE} ${FILE}.cpp ${SOURCE_BASE_FILES} ${SOURCE_BASE_HEADER_FILES} )
  target_link_libraries (${FILE}  ${DGTAL_LIBRARIES} )
ENDFOREACH(FILE)

ADD_SUBDIRECTORY(renderer)
