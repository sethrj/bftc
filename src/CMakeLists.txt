add_custom_target(exec)

#include blitz
find_package(Trilinos COMPONENTS epetra teuchos triutils REQUIRED)
#find_package(Trilinos COMPONENTS epetra aztecoo)

include_directories(${Trilinos_INCLUDE_DIR})

set(EXEC_NAMES
  ex15
  )

# link each target to appropriate libs
foreach(TARGET_NAME ${EXEC_NAMES})
  add_executable(${TARGET_NAME}
    "${TARGET_NAME}.cpp"
    )
  target_link_libraries(${TARGET_NAME}
    ${Trilinos_LIBRARIES}
    )
endforeach(TARGET_NAME)

# new target: build examples builds these
add_dependencies(exec ${EXEC_NAMES})
