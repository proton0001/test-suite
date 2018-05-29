set(PROG bfs)
set(RUN_OPTIONS ${CMAKE_CURRENT_SOURCE_DIR}/graph.txt)
set(WORKDIR ${CMAKE_CURRENT_BINARY_DIR})
llvm_multisource()
