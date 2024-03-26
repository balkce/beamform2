cmake_minimum_required(VERSION 3.8)
project(beamform2)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter)
endif()

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(jack_msgs REQUIRED)
find_package(rosjack REQUIRED)

find_package(rosidl_default_generators REQUIRED)

set(msg_files
  "msg/InterfTheta.msg"
)

rosidl_generate_interfaces(${PROJECT_NAME}
  ${msg_files}
  DEPENDENCIES std_msgs
)

INCLUDE(FindPkgConfig)
pkg_check_modules (GLIB2 glib-2.0)
#message( "Found glib-2.0: " ${GLIB2_FOUND} )
pkg_check_modules (FFTW fftw3)
#message( "Found fftw3: " ${FFTW_FOUND} )
pkg_check_modules (EIGEN eigen3)
#message( "Found eigen: " ${EIGEN_FOUND} )
pkg_check_modules (JSON jsoncpp)
#message( "Found jsoncpp: " ${JSON_FOUND} )

include_directories(
  ${GLIB2_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIRS}
  ${FFTW_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
)

rosidl_get_typesupport_target(cpp_typesupport_target
  ${PROJECT_NAME} rosidl_typesupport_cpp)

set(dependencies
  rclcpp std_msgs jack_msgs rosjack 
)

add_executable(das src/das.cpp)
ament_target_dependencies(das ${dependencies})
target_link_libraries(das
  ${FFTW_LIBRARIES}
)

add_executable(mvdr src/mvdr.cpp)
ament_target_dependencies(mvdr ${dependencies})
target_link_libraries(mvdr
  ${FFTW_LIBRARIES}
)

add_executable(gsc src/gsc.cpp)
ament_target_dependencies(gsc ${dependencies})
target_link_libraries(gsc
  ${FFTW_LIBRARIES}
)

add_executable(lcmv src/lcmv.cpp)
ament_target_dependencies(lcmv ${dependencies})
target_link_libraries(lcmv
  "${cpp_typesupport_target}"
  ${FFTW_LIBRARIES}
)

add_executable(gss src/gss.cpp)
ament_target_dependencies(gss ${dependencies})
target_link_libraries(gss
  "${cpp_typesupport_target}"
  ${FFTW_LIBRARIES}
)

add_executable(phase src/phase.cpp)
ament_target_dependencies(phase ${dependencies})
target_link_libraries(phase
  ${FFTW_LIBRARIES}
)

add_executable(phasempf src/phasempf.cpp)
ament_target_dependencies(phasempf ${dependencies})
target_link_libraries(phasempf
  ${FFTW_LIBRARIES}
)

add_executable(rosjack_read src/jack_read.cpp)
ament_target_dependencies(rosjack_read ${dependencies})

add_executable(rosjack_write src/jack_write.cpp)
ament_target_dependencies(rosjack_write ${dependencies})

add_executable(rosjack_ref src/jack_ref.cpp)
ament_target_dependencies(rosjack_ref ${dependencies})
target_link_libraries(rosjack_ref
  ${FFTW_LIBRARIES}
)

add_executable(mcra src/mcra.cpp)
ament_target_dependencies(mcra ${dependencies})
target_link_libraries(mcra
  ${FFTW_LIBRARIES}
)

install(TARGETS das
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS mvdr
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS gsc
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS lcmv
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS gss
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS phase
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS phasempf
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS rosjack_read
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS rosjack_write
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS rosjack_ref
  DESTINATION lib/${PROJECT_NAME}
)
install(TARGETS mcra
  DESTINATION lib/${PROJECT_NAME}
)

install(DIRECTORY
  launch
  DESTINATION share/${PROJECT_NAME}
)
install(DIRECTORY
  config
  DESTINATION share/${PROJECT_NAME}
)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_export_dependencies(rosidl_default_runtime)
ament_package()
