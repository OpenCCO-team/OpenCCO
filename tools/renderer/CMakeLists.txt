SET(TOOLS_RENDERER_SRC
	graph2render
	animationRenderer
)

SET(TOOLS_RENDERER_HEADERS
	TreeImageRenderer.h
)

SET(TOOLS_RENDER_SRC
	TreeImageRenderer.cpp
)

SET(SVG_EXPORT_HEADERS
	src_svg_export/streamable.h
	src_svg_export/svg_attributes.h
	src_svg_export/svg_elements.h
)

SET(SVG_EXPORT_SRC
    src_svg_export/streamable.cpp
	src_svg_export/svg_attributes.cpp
	src_svg_export/svg_elements.cpp
)

# build a lib so that the svg_export sources are not compiled once per executable
ADD_LIBRARY(svg_export_lib ${SVG_EXPORT_SRC} ${SVG_EXPORT_HEADERS})

FOREACH(FILE ${TOOLS_RENDERER_SRC})
	add_executable(${FILE} ${FILE}.cpp ${TOOLS_RENDER_SRC} ${SOURCE_BASE_HEADER_FILES} ${TOOLS_RENDERER_HEADERS})
	target_link_libraries (${FILE} PRIVATE ${CERES_LIBRARIES} ${DGTAL_LIBRARIES} svg_export_lib)
ENDFOREACH(FILE)
