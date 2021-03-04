# Loop over all available xml files
FILE(GLOB xmls [a-zA-Z\-]*.xml)

FOREACH(xml ${xmls})
  GET_FILENAME_COMPONENT(xml ${xml} NAME)
  CONFIGURE_FILE(${xml} ${xml})
ENDFOREACH()
