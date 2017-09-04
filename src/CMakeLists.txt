cmake_minimum_required (VERSION 2.8.8)

project(KGS)

set(KGS_MAJOR_VERSION 0)
set(KGS_MINOR_VERSION 5)
set(KGS_PATCH_VERSION 1)
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
	message (STATUS "${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS}")
add_library( libKGS

		core/Atom.h
		core/Chain.h
		core/Configuration.h
		core/Coordinate.h
		Color.h
		CTKTimer.h
		DisjointSets.h
		core/Grid.h
		HbondIdentifier.h
		IO.h
		Logger.h
		core/Molecule.h
		core/Bond.h
		core/HBond.h
		core/DBond.h
		core/Residue.h
		ResidueProfiles.h
		core/Rigidbody.h
		core/graph/KinGraph.h
		applications/options/ApplicationOptions.h
		applications/options/ExploreOptions.h
		applications/options/TransitionOptions.h
		applications/options/RelativeTransitionOptions.h
		applications/options/RigidityOptions.h
		applications/options/AlignOptions.h
		applications/options/HierarchyOptions.h
		planners/SamplingPlanner.h
		planners/DihedralRRT.h
		planners/RRTPlanner.h
		planners/MCMCPlanner.h
		planners/PoissonPlanner.h
		planners/PoissonPlanner2.h
		planners/BidirectionalMovingFront.h
		Selection.h
		SiteID.h
		Util.h
		KGSExceptions.h
		freeMKL.h
		moves/Move.h
		moves/RawMove.h
		#        moves/CompositeMove.h
		moves/NullspaceMove.h
		moves/RebuildMove.h
		moves/FastClashAvoidingMove.h
		moves/ClashAvoidingMove.h
		moves/DecreaseStepMove.h
		directions/Direction.h
		directions/BlendedDirection.h
		directions/RandomDirection.h
		directions/DihedralDirection.h
		directions/MSDDirection.h
		directions/LSNullspaceDirection.h
		directions/VDWDirection.h
		math/MathUtility.h
		math/SVD.h
		math/infnan.h
		math/math.h
		math3d/primitives.h
		metrics/Metric.h
		metrics/RMSD.h
		metrics/RMSDnosuper.h

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
        #JacobianRelated.cpp
        Logger.cpp
        math/MathUtility.cpp
        core/Molecule.cpp
        core/Bond.cpp
        core/HBond.cpp
        core/DBond.cpp
        core/Residue.cpp
        core/Rigidbody.cpp
        core/graph/KinGraph.cpp
        math/SVD.cpp
        applications/options/ApplicationOptions.cpp
		applications/options/ExploreOptions.cpp
		applications/options/TransitionOptions.cpp
		applications/options/RelativeTransitionOptions.cpp
		applications/options/RigidityOptions.cpp
		applications/options/AlignOptions.cpp
		applications/options/HierarchyOptions.cpp
        planners/SamplingPlanner.cpp
        planners/DihedralRRT.cpp
        planners/RRTPlanner.cpp
		planners/MCMCPlanner.cpp
        planners/PoissonPlanner.cpp
        planners/PoissonPlanner2.cpp
        planners/BidirectionalMovingFront.cpp
        Selection.cpp
        #core/Transformation.cpp
        Util.cpp
        KGSExceptions.cpp
        freeMKL.cpp
        moves/Move.cpp
        moves/RawMove.cpp
#        moves/CompositeMove.cpp
        moves/NullspaceMove.cpp
        moves/RebuildMove.cpp
        moves/FastClashAvoidingMove.cpp
        moves/ClashAvoidingMove.cpp
        moves/DecreaseStepMove.cpp
        directions/Direction.cpp
        directions/BlendedDirection.cpp
        directions/RandomDirection.cpp
        directions/DihedralDirection.cpp
        directions/MSDDirection.cpp
        directions/LSNullspaceDirection.cpp
        directions/VDWDirection.cpp
        metrics/Dihedral.cpp
        metrics/Metric.cpp
        metrics/RMSD.cpp
        metrics/RMSDnosuper.cpp
        math3d/primitives.cpp

        math/SVDGSL.cpp
        math/SVDGSL.h
        math/SVDMKL.cpp
        math/SVDMKL.h
        math/CudaSVD.cpp
        math/CudaSVD.h
        math/gsl_helpers.cpp
        math/gsl_helpers.h
        math/Nullspace.cpp
        math/Nullspace.h
        math/NullspaceSVD.cpp
        math/NullspaceSVD.h
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
        core/dofs/FixedLink.cpp
        core/dofs/FixedLink.h
        core/graph/KinEdge.cpp
        core/graph/KinEdge.h
        core/graph/KinVertex.cpp
        core/graph/KinVertex.h
        core/graph/KinTree.cpp
        core/graph/KinTree.h
		directions/LSNrelativeDirection.cpp
        directions/LSNrelativeDirection.h
        core/graph/Constraint.cpp
        core/graph/Constraint.h
        core/ConstraintSet.cpp
        core/ConstraintSet.h )

set_target_properties( libKGS PROPERTIES PREFIX "" )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(GSL REQUIRED)
#find_package(GSL)
#if( NOT ${GSL_FOUND})
#	message(STATUS "GSL not found, downloading and compiling from source")
#
#    include(ExternalProject)
#	set(GSLVersion 2.1)
#	ExternalProject_Add(GSL
#			PREFIX extern/GSL
#			URL ftp://ftp.gnu.org/gnu/gsl/gsl-${GSLVersion}.tar.gz
#			CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/extern --enable-shared=no --with-pic=yes
#			BUILD_COMMAND           make -j8
#			INSTALL_COMMAND         make install
#			BUILD_IN_SOURCE 1)
#	set(GSL_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib/ )
#	set(GSL_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include/)
#       set(GSL_LIBRARIES -lgsl -lgslcblas -lm)
#    add_dependencies(libKGS GSL)
#endif()
link_directories(   ${GSL_LIBRARY_DIR} )
include_directories(${GSL_INCLUDE_DIR} )

option(ForceGSL "ForceGSL" OFF)

if( NOT ${ForceGSL} )
  # Attempt to find MKL
  find_package(MKL REQUIRED)
  if(${MKL_FOUND})
    include_directories(${MKL_INCLUDE_DIR})
    #list(APPEND CMAKE_EXE_LINKER_FLAGS -lmkl_intel_lp64)
    link_libraries( mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__INTEL_MKL" )
  endif()
endif()

#find_package(CUDA)
#if(${CUDA_FOUND})
#    link_libraries(${CUDA_LIBRARIES} ${CUDA_cusolver_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
#    include_directories(${CUDA_INCLUDE_DIRS})
#    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__CUDA_SVD" )
#endif()


include_directories(${CMAKE_SOURCE_DIR})		#Absolute header paths
#include_directories(${CMAKE_SOURCE_DIR}/utils/) #LoopTK headers

add_executable( kgs_explore         applications/MainExplore.cpp)
add_executable( kgs_rigidity        applications/MainRigidity.cpp )
#add_executable( kgs_print_torsions  applications/MainTorsionPrinter.cpp )
add_executable( kgs_ensemble_pca    applications/MainEnsemblePCA.cpp )
add_executable( kgs_rmsd            applications/MainRmsd.cpp )
add_executable( kgs_align            applications/MainAlign.cpp )
add_executable( kgs_transition      applications/MainTransition.cpp )
add_executable( kgs_loopexplore     applications/MainLoopExplore.cpp )
add_executable( kgs_hierarchy       applications/MainHierarchy.cpp )
#add_executable( kgs_mcmc            applications/MainMCMC.cpp )
add_executable( kgs_test            applications/MainTest.cpp )
add_executable( kgs_relative_transition            applications/MainRelativeTransition.cpp )

target_link_libraries( kgs_explore          libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_rigidity         libKGS ${GSL_LIBRARIES} )
#target_link_libraries( kgs_print_torsions   libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_ensemble_pca     libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_rmsd             libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_align             libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_transition       libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_loopexplore      libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_hierarchy        libKGS ${GSL_LIBRARIES} )
#target_link_libraries( kgs_mcmc             libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_test        libKGS ${GSL_LIBRARIES} )
target_link_libraries( kgs_relative_transition             libKGS ${GSL_LIBRARIES} )

set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR} )



