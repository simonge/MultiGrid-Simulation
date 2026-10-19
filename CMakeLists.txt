################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

cmake_minimum_required(VERSION 3.30..)

project(MG VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Request the Geant4 features used by this application: GDML geometry,
# interactive UI, and visualization. ui_all/vis_all enable the drivers
# available in the selected Geant4 installation.
find_package(Geant4 REQUIRED COMPONENTS gdml ui_all vis_all)

# Set NCrystalGeant4_DIR when configuring if the package is not on CMake's
# default search path (setup.sh does this for the repository-local venv).
find_package(NCrystalGeant4 0.0.1 CONFIG REQUIRED)

# Link only the ROOT components used by the application. ROOT's imported
# targets carry their include paths and transitive link dependencies.
find_package(ROOT CONFIG REQUIRED COMPONENTS Core RIO Tree Hist)

file(GLOB MG_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cc")
file(GLOB MG_HEADERS CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/*.hh")

add_executable(MG MG.cc ${MG_SOURCES} ${MG_HEADERS})
target_include_directories(MG PRIVATE
  "${PROJECT_SOURCE_DIR}/include"
  ${Geant4_INCLUDE_DIRS}
)
target_link_libraries(MG PRIVATE
  NCrystalGeant4::NCrystalGeant4
  ${Geant4_LIBRARIES}
  ROOT::Core
  ROOT::RIO
  ROOT::Tree
  ROOT::Hist
)

include(GNUInstallDirs)
install(TARGETS MG RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
