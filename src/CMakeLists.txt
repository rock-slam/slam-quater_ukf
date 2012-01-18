rock_library(quater_ukf
    SOURCES Dummy.cpp
    HEADERS Dummy.hpp)

rock_executable(quater_ukf_bin Main.cpp
    DEPS quater_ukf)

