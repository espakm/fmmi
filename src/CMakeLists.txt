add_library(fmmi INTERFACE)
target_sources(fmmi INTERFACE fmmi/smatrix.hpp fmmi/dmatrix.hpp fmmi/fmmi.hpp)
target_include_directories(fmmi INTERFACE .)
