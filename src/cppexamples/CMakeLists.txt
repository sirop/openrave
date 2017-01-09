cmake_minimum_required (VERSION 2.6.0)
project (cppexamples)

set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../..) # only used for cppexamples!
find_package(OpenRAVE REQUIRED)

set(OpenBLAS_DIR /opt/OpenBLAS/lib/cmake/openblas)
find_package(OpenBLAS REQUIRED)
set(OpenBLAS_LIBRARY_DIRS /opt/OpenBLAS/lib)
set(OpenBLAS_LIBRARIES /opt/OpenBLAS/lib/libopenblas.so)

#find_library(GFORTRAN_LIBRARY NAMES gfortran)

if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -Wall")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

find_package(Boost ${OpenRAVE_Boost_VERSION} EXACT COMPONENTS iostreams python thread system)

include_directories(${OpenRAVE_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${OpenBLAS_INCLUDE_DIRS})
if( Boost_INCLUDE_DIRS )
  include_directories(${Boost_INCLUDE_DIRS})
endif()

link_directories(${OpenRAVE_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS} ${OpenBLAS_LIBRARY_DIR})

macro(build_openrave_executable name)
  add_executable(${name} ${name}.cpp)
  set_target_properties(${name} PROPERTIES COMPILE_FLAGS "${OpenRAVE_CXX_FLAGS}")
  set_target_properties(${name} PROPERTIES LINK_FLAGS "${OpenRAVE_LINK_FLAGS}")
  target_link_libraries(${name} ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${OpenBLAS_LIBRARIES})
  install(TARGETS ${name} DESTINATION . )
endmacro(build_openrave_executable)

macro(build_openrave_plugin name)
add_library(${name} SHARED ${name}.cpp)
set_target_properties(${name} PROPERTIES COMPILE_FLAGS "${OpenRAVE_CXX_FLAGS}")
set_target_properties(${name} PROPERTIES LINK_FLAGS "${OpenRAVE_LINK_FLAGS}")
target_link_libraries(${name} ${OpenRAVE_LIBRARIES} ${OpenBLAS_LIBRARIES})
install(TARGETS ${name} DESTINATION . )
endmacro(build_openrave_plugin)

build_openrave_plugin(plugincpp)
#build_openrave_plugin(customreader)

#build_openrave_executable(orcollision)
#build_openrave_executable(orconveyormovement)
#build_openrave_executable(orloadviewer)
#build_openrave_executable(ikfastloader)
#build_openrave_executable(orikfilter)
build_openrave_executable(jacobian)
build_openrave_executable(jacobianTime)
#build_openrave_executable(ormulticontrol)
#build_openrave_executable(ormultithreadedplanning)
#build_openrave_executable(orpr2turnlever)
#build_openrave_executable(orplanning_module)
#build_openrave_executable(orplanning_multirobot)
#build_openrave_executable(orplanning_planner)
#build_openrave_executable(orplanning_door)
#build_openrave_executable(orplanning_ik)
#build_openrave_executable(orshowsensors)
#build_openrave_executable(ortrajectory)



