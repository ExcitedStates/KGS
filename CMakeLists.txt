cmake_minimum_required (VERSION 2.8.8)

project(KGS)

set(KGS_MAJOR_VERSION 0)
set(KGS_MINOR_VERSION 4)
set(KGS_PATCH_VERSION 2)
set(KGS_VERSION ${KGS_MAJOR_VERSION}.${KGS_MINOR_VERSION}.${KGS_PATCH_VERSION})




# Enable C++11 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

add_library( libKGS

		core/Atom.h
		core/Chain.h
		Color.h
		core/Configuration.h
		core/Coordinate.h
		CTKTimer.h
		DisjointSets.h
		core/Grid.h
		HbondIdentifier.h
		IO.h
		JacobianRelated.h
		Logger.h
		math/MathUtility.h
		core/Molecule.h
		core/Bond.h
		core/ProteinHBond.h
		core/Residue.h
		ResidueProfiles.h
		core/Rigidbody.h
		core/graph/KinGraph.h
		RunFirst.h
		math/SVD.h
		SamplingOptions.h
		planners/SamplingPlanner.h
		planners/DihedralRRT.h
		planners/RRTPlanner.h
		planners/PoissonPlanner.h
		planners/BidirectionalMovingFront.h
		Selection.h
		SiteID.h
		core/Transformation.h
		Util.h
		KGSExceptions.h
		freeMKL.h
		moves/Move.h
		moves/RawMove.h
		moves/CompositeMove.h
		moves/NullspaceMove.h
		moves/RebuildMove.h
		moves/ClashAvoidingMove.h
		directions/Direction.h
		directions/MixedDirection.h
		directions/RandomDirection.h
		directions/DihedralDirection.h
		directions/MSDDirection.h
		directions/LSNullspaceDirection.h
		directions/VDWDirection.h
		math/infnan.h
		math/math.h
		math3d/primitives.h


		core/Atom.cpp
		core/Chain.cpp
		Color.cpp
		core/Configuration.cpp
		core/Coordinate.cpp
		CTKTimer.cpp
		DisjointSets.cpp
		core/Grid.cpp
		HbondIdentifier.cpp
		IO.cpp
		JacobianRelated.cpp
		Logger.cpp
		math/MathUtility.cpp
		core/Molecule.cpp
		core/Bond.cpp
		core/ProteinHBond.cpp
		core/Residue.cpp
		core/Rigidbody.cpp
		core/graph/KinGraph.cpp
		RunFirst.cpp
		math/SVD.cpp
		SamplingOptions.cpp
		planners/SamplingPlanner.cpp
		planners/DihedralRRT.cpp
		planners/RRTPlanner.cpp
		planners/PoissonPlanner.cpp
		planners/BidirectionalMovingFront.cpp
		Selection.cpp
		core/Transformation.cpp
		Util.cpp
		KGSExceptions.cpp
		freeMKL.cpp
		moves/Move.cpp
		moves/RawMove.cpp
		moves/CompositeMove.cpp
		moves/NullspaceMove.cpp
		moves/RebuildMove.cpp
		moves/ClashAvoidingMove.cpp
		directions/Direction.cpp
		directions/MixedDirection.cpp
		directions/RandomDirection.cpp
		directions/DihedralDirection.cpp
		directions/MSDDirection.cpp
		directions/LSNullspaceDirection.cpp
		directions/VDWDirection.cpp
		metrics/Dihedral.cpp
		metrics/RMSD.cpp
		math3d/primitives.cpp

		math/GSLSVD.cpp
		math/GSLSVD.h
		math/MKLSVD.cpp
		math/MKLSVD.h
		math/gsl_helpers.cpp
		math/gsl_helpers.h
		math/Nullspace.cpp
		math/Nullspace.h
		loopclosure/ExactIK.cpp
		loopclosure/ExactIK.h
		core/dofs/DOF.cpp
		core/dofs/DOF.h
        core/dofs/TorsionDOF.cpp
        core/dofs/TorsionDOF.h
        core/dofs/GlobalTranslateDOF.cpp
        core/dofs/GlobalTranslateDOF.h
		core/dofs/GlobalRotateDOF.cpp
		core/dofs/GlobalRotateDOF.h
		core/graph/KinEdge.cpp
		core/graph/KinEdge.h
		core/graph/KinVertex.cpp
		core/graph/KinVertex.h
		core/graph/KinTree.cpp
		core/graph/KinTree.h
		)

set_target_properties( libKGS PROPERTIES PREFIX "" )


find_package(GSL)
if( NOT ${GSL_FOUND})
	message(STATUS "GSL not found, downloading and compiling from source")

	include(ExternalProject)
	set(GSLVersion 2.1)
	ExternalProject_Add(GSL
			PREFIX extern/GSL
			URL ftp://ftp.gnu.org/gnu/gsl/gsl-${GSLVersion}.tar.gz
			CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/extern --enable-shared=no --with-pic=yes
			BUILD_COMMAND           make -j8
			INSTALL_COMMAND         make install
			BUILD_IN_SOURCE 1)
	set(GSL_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib/ )
	set(GSL_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include/)
	set(GSL_LIBRARIES -lgsl -lgslcblas -lm)
	add_dependencies(libKGS GSL)
endif()
link_directories(   ${GSL_LIBRARY_DIR} )
include_directories(${GSL_INCLUDE_DIR} )
#link_libraries( gsl gslcblas m )


# Attempt to find MKL
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
find_package(MKL REQUIRED)
if(${MKL_FOUND})
	include_directories(${MKL_INCLUDE_DIR})
	#list(APPEND CMAKE_EXE_LINKER_FLAGS -lmkl_intel_lp64)
	link_libraries( mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread )
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__INTEL_COMPILER -D__INTEL_MKL" )
endif()

include_directories(${CMAKE_SOURCE_DIR})		#Absolute header paths
#include_directories(${CMAKE_SOURCE_DIR}/utils/) #LoopTK headers

add_executable( kgs             applications/MainKGS.cpp)
add_executable( rigidity        applications/MainRigidity.cpp )
add_executable( print_torsions  applications/MainTorsionPrinter.cpp )
add_executable( ensemble_pca    applications/MainEnsemblePCA.cpp )
add_executable( rmsd		    applications/MainRmsd.cpp )
add_executable( testik		    applications/MainTestExactIK.cpp )
add_executable( mytest  	    applications/MainTest.cpp )

target_link_libraries( kgs              libKGS ${GSL_LIBRARIES} )
target_link_libraries( rigidity         libKGS ${GSL_LIBRARIES} )
target_link_libraries( print_torsions   libKGS ${GSL_LIBRARIES} )
target_link_libraries( ensemble_pca     libKGS ${GSL_LIBRARIES} )
target_link_libraries( rmsd			    libKGS ${GSL_LIBRARIES} )
target_link_libraries( testik		    libKGS ${GSL_LIBRARIES} )
target_link_libraries( mytest  		    libKGS ${GSL_LIBRARIES} )

set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR} )



