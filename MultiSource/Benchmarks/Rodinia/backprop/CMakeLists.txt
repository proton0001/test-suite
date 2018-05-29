set(PROG backprop)
list(APPEND LDFLAGS -lm)

if(SMALL_PROBLEM_SIZE)
  list(APPEND CFLAGS -DSMALL_DATASET)
endif()

set(WORKDIR ${CMAKE_CURRENT_BINARY_DIR})
llvm_multisource()
