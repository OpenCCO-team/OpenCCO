SET(TESTS_SRC
  testCompCCO
  testGenerateTree
  testGeom
  testKamiya
  testSchreiner
)


FOREACH(FILE ${TESTS_SRC})
  add_executable(${FILE} ${FILE}.cpp ${SOURCE_BASE_FILES} ${SOURCE_BASE_HEADER_FILES} )
  target_link_libraries (${FILE} ${CERES_LIBRARIES} ${DGTAL_LIBRARIES} )
ENDFOREACH(FILE)
target_compile_definitions(testGenerateTree PUBLIC SAMPLE_DIR="${PROJECT_SOURCE_DIR}/samples/")
target_compile_definitions(testGeom PUBLIC SAMPLE_DIR="${PROJECT_SOURCE_DIR}/samples/")
