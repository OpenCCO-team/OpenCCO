
SET(CCO_TESTS
   testGeomTools
   testCoronaryArteryTree
   testRescale
   TestGeneration
   TestGeneration3D
   TEST
   TestGenerationTreshold
   TestImage
)




FOREACH(FILE ${CCO_TESTS})
  add_executable(${FILE} ${FILE} ${CCO_SRC})
   target_link_libraries (${FILE} ${DGTAL_LIBRARIES}  )
ENDFOREACH(FILE)


