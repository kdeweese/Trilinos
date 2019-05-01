# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    ImageBased
    SOURCES Driver_ImageBased.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    ImageBased
    NAME "ImageBased_Laplace2D"
    ARGS "--linAlgebra=Tpetra --xml=image_based.xml --matrixType=Laplace2D --nx=10 --ny=10 --greyscale=auto"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    ImageBased
    NAME "ImageBased_smoothed_Laplace2D"
    ARGS "--linAlgebra=Tpetra --xml=image_based_smoothed.xml --matrixType=Laplace2D --nx=10 --ny=10 --greyscale=auto"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )


  TRIBITS_COPY_FILES_TO_BINARY_DIR(ImageBased_cp
    SOURCE_FILES image_based.xml image_based_smoothed.xml
    )


ENDIF()
