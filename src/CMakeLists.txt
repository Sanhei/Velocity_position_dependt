

# problem 5.1
add_executable(D_trans
  main.cpp
)

target_link_libraries(D_trans
  ${Boost_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

install(TARGETS D_trans
  DESTINATION bin
)

