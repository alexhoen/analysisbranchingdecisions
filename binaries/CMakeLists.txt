include(GenerateExportHeader)

set(CMAKE_CXX_STANDARD 14)

include_directories(${SCIP_INCLUDE_DIRS} ${SOPLEX_INCLUDE_DIRS} ${PAPILO_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})
find_package(Boost REQUIRED iostreams serialization program_options)

set(CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${Boost_LIBRARIES})
check_cxx_source_compiles(
   "#include <boost/iostreams/filter/gzip.hpp>
    int main() { auto decomp = boost::iostreams::gzip_decompressor(); (void)decomp; return 0; }"
   EXACT_USE_BOOST_IOSTREAMS_WITH_ZLIB )
check_cxx_source_compiles(
   "#include <boost/iostreams/filter/bzip2.hpp>
    int main() { auto decomp = boost::iostreams::bzip2_decompressor(); (void)decomp; return 0; }"
   EXACT_USE_BOOST_IOSTREAMS_WITH_BZIP2 )


add_executable(exact-executable ${CMAKE_CURRENT_LIST_DIR}/../src/exact.cpp)

set_target_properties(exact-executable PROPERTIES OUTPUT_NAME exact RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
target_link_libraries(exact-executable exact-core  ${SCIP_LIBRARIES} ${SOPLEX_LIBRARIES} ${PAPILO_IMPORTED_TARGETS} ${Boost_LIBRARIES})
target_compile_definitions(exact-executable PRIVATE EXACT_USE_EXTERN_TEMPLATES EXACT_MPS_WRITER)

