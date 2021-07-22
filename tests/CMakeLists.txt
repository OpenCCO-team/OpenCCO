SET(TESTS_SRC
  testCCO
  testGeom
  testKamiya
  testSchreiner
)

FOREACH(FILE ${TESTS_SRC})
  add_executable(${FILE} ${FILE}.cpp  ${SOURCE_BASE_HEADER_FILES} ${SOURCE_BASE_FILES})
  target_link_libraries (${FILE} Ceres::ceres ${DGTAL_LIBRARIES})
ENDFOREACH(FILE)

