set(PROG hotspot)

if(SMALL_PROBLEM_SIZE)
  list(APPEND CFLAGS -DSMALL_DATASET)
endif()
set(RUN_OPTIONS ${CMAKE_CURRENT_SOURCE_DIR}/temp_file.txt ${CMAKE_CURRENT_SOURCE_DIR}/power_file.txt)

set(WORKDIR ${CMAKE_CURRENT_BINARY_DIR})
llvm_multisource()
