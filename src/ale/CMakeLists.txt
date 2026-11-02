
set(kalerax_SOURCES ale.cpp
  AleArguments.cpp
  AleEvaluator.cpp
  AleOptimizer.cpp
  AleState.cpp
  Highways.cpp
  OptimizationClasses.cpp
  TrimFamilies.cpp
  )

# kalerax: the WGD/LORe fork of AleRax. The binary is named distinctly from the
# upstream AleRax executable.
add_program(kalerax "${kalerax_SOURCES}")
install(TARGETS kalerax DESTINATION bin)

# WGD regression test (not installed): _hasWGD all-false must reproduce the
# original DL likelihood. See tests/wgd_regression.cpp.
add_program(wgd_regression "${CMAKE_SOURCE_DIR}/tests/wgd_regression.cpp")

# LORe resolution-branch marginal readout test (not installed): the U-aware
# backtrace posterior under the fitted global r. See tests/lore_marginal.cpp.
add_program(lore_marginal "${CMAKE_SOURCE_DIR}/tests/lore_marginal.cpp")
# WGD/LORe-aware scenario sampler vs likelihood identity (synteny probe).
# See tests/export_consistency.cpp.
add_program(export_consistency "${CMAKE_SOURCE_DIR}/tests/export_consistency.cpp")

# LORe AORe-vs-LORe model-selection sanity check (the old r in {0,1} profile
# proxy, kept as a coarse sanity check). See tests/lore_modelselection_sanity.cpp.
add_program(lore_modelselection_sanity
            "${CMAKE_SOURCE_DIR}/tests/lore_modelselection_sanity.cpp")

