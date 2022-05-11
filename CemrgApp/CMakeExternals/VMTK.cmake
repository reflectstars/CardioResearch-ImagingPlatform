# ----------
# VMTK
# ----------
set(proj VMTK)
set(proj_DEPENDENCIES ITK VTK)
set(${proj}_DEPENDS ${proj})

if(MITK_USE_VMTK)
  # 'Sanity' checks
  if(DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR})
    message(FATAL_ERROR "${proj}_DIR variable is defined but corresponds to non-existing directory")
  endif()

  if(NOT DEFINED ${proj}_DIR)
    set(additional_cmake_args)

    if(CTEST_USE_LAUNCHERS)
      list(APPEND additional_cmake_args
      -DCMAKE_PROJECT_VTK_VMTK_INCLUDE:FILEPATH=${CMAKE_ROOT}/Modules/CTestUseLaunchers.cmake
    )
    endif()

    if(MITK_USE_Python)
      list(APPEND additional_cmake_args
        -DVT