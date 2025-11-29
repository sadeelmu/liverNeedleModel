cmake_minimum_required(VERSION 3.13.4)

project(AliMuwahedExtension)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/AliMuwahedExtension")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Caroline Essert (Univ. Strasbourg), Kinan Ali (HealthTech), and Sadel Muwahed(HealthTech)")
set(EXTENSION_DESCRIPTION "This extension was created by Kinan Ali and Sadel Muwahed")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/AliMuwahedExtension.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/AliMuwahedExtension/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(AnotherModule)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
