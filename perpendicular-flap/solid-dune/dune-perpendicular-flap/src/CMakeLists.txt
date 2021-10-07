find_package(precice REQUIRED CONFIG)

set(programs
		solid-dune)

#include precice
find_package(precice)

foreach(_program ${programs})
	add_executable(${_program} ${_program}.cc)
	target_link_libraries(${_program} PRIVATE precice::precice)
endforeach()
