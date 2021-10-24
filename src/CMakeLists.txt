set(FMMI_SOURCES
    fmmi/dmatrix.hpp
    fmmi/dmatrix_rec.hpp
    fmmi/smatrix.hpp
    fmmi/smatrix_rec.hpp
)

add_library(fmmi INTERFACE)
target_sources(fmmi INTERFACE ${FMMI_SOURCES})
target_include_directories(fmmi INTERFACE .)
