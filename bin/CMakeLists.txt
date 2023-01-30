SET(BIN_SRC
  generateTree2D
  generateTree3D
)

FOREACH(FILE ${BIN_SRC})
  add_executable(${FILE} ${FILE}.cpp ${SOURCE_BASE_FILES} ${SOURCE_BASE_HEADER_FILES} )
  target_link_libraries (${FILE} Ceres::ceres ${DGTAL_LIBRARIES} )
ENDFOREACH(FILE)
target_compile_definitions(testGeom PUBLIC SAMPLE_DIR="${PROJECT_SOURCE_DIR}/samples/")
