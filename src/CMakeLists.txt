rock_find_pkgconfig(Eigen3 eigen3 REQUIRED)

rock_library(quater_ukf
    SOURCES ukf.cpp
    HEADERS ukf.hpp)

rock_executable(quater_ukf_bin Main.cpp
    DEPS quater_ukf)

