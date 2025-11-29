#-----------------------------------------------------------------------------
set(MODULE_NAME AliMuwahedModule)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${AliMuwahedModule}.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${AliMuwahedModule}.png
  Resources/UI/${AliMuwahedModule}.ui
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${AliMuwahedModule}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${AliMuwahedModule}.py)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
