<a id="developers-run-clang-tidy-locally"></a>

# How to run the clang-tidy linter

WarpX’s CI tests include several checks performed with the [clang-tidy](https://clang.llvm.org/extra/clang-tidy/) linter.
The complete list of checks performed is defined in the `.clang-tidy` configuration file.

### clang-tidy configuration file

```yaml
Checks: '
    -*,
    bugprone-*,
        -bugprone-easily-swappable-parameters,
        -bugprone-implicit-widening-of-multiplication-result,
        -bugprone-misplaced-widening-cast,
        -bugprone-unchecked-optional-access,
        -bugprone-unused-local-non-trivial-variable,
    cert-*,
        -cert-err58-cpp,
        -cert-int09-c,
    clang-analyzer-*,
        -clang-analyzer-optin.performance.Padding,
        -clang-analyzer-optin.mpi.MPI-Checker,
        -clang-analyzer-osx.*,
        -clang-analyzer-optin.osx.*,
    clang-diagnostic-*,
    cppcoreguidelines-*,
        -cppcoreguidelines-avoid-c-arrays,
        -cppcoreguidelines-avoid-magic-numbers,
        -cppcoreguidelines-avoid-non-const-global-variables,
        -cppcoreguidelines-init-variables,
        -cppcoreguidelines-macro-usage,
        -cppcoreguidelines-missing-std-forward,
        -cppcoreguidelines-narrowing-conversions,
        -cppcoreguidelines-non-private-member-variables-in-classes,
        -cppcoreguidelines-owning-memory,
        -cppcoreguidelines-pro-*,
    google-build-explicit-make-pair,
    google-build-namespaces,
    google-explicit-constructor,
    google-global-names-in-headers,
    misc-*,
        -misc-no-recursion,
        -misc-non-private-member-variables-in-classes,
        -misc-include-cleaner,
        -misc-use-internal-linkage,
    modernize-*,
        -modernize-avoid-c-arrays,
        -modernize-return-braced-init-list,
        -modernize-use-trailing-return-type,
        -modernize-use-constraints,
        -modernize-use-designated-initializers,
        -modernize-use-std-numbers,
        -modernize-use-ranges,
        -modernize-use-integer-sign-comparison,
        -modernize-use-starts-ends-with,
        -modernize-use-std-format,
    mpi-*,
    performance-*,
        -performance-enum-size,
    portability-*,
    readability-*,
        -readability-avoid-nested-conditional-operator,
        -readability-convert-member-functions-to-static,
        -readability-else-after-return,
        -readability-enum-initial-value,
        -readability-function-cognitive-complexity,
        -readability-identifier-length,
        -readability-implicit-bool-conversion,
        -readability-isolate-declaration,
        -readability-magic-numbers,
        -readability-math-missing-parentheses,
        -readability-named-parameter,
        -readability-uppercase-literal-suffix,
        -readability-redundant-casting,
        -readability-container-contains,
        -readability-redundant-inline-specifier,
        -readability-redundant-member-init,
        -readability-reference-to-constructed-temporary
    '

CheckOptions:
- key:          bugprone-narrowing-conversions.WarnOnIntegerToFloatingPointNarrowingConversion
  value:        "false"
- key:          misc-definitions-in-headers.HeaderFileExtensions
  value:        "H,"
- key:          modernize-pass-by-value.ValuesOnly
  value:        "true"
- key:          misc-use-anonymous-namespace.HeaderFileExtensions
  value:        "H,"
- key:          performance-move-const-arg.CheckTriviallyCopyableMove
  value:        "false"

HeaderFilterRegex: 'Source[a-z_A-Z0-9\/]+\.H$'

# TODO Modernize this after switching to C++20:
# -modernize-use-constraints
# -modernize-use-designated-initializers
# -modernize-use-std-numbers
# -modernize-use-ranges
# -modernize-use-integer-sign-comparison
# -modernize-use-starts-ends-with
# -modernize-use-std-format
# -readability-container-contains

# TODO Consider enabling the following checks:
# -bugprone-unused-local-non-trivial-variable
# -misc-use-internal-linkage
# -performance-enum-size
# -readability-container-contains
# -readability-enum-initial-value
# -readability-math-missing-parentheses
# -readability-redundant-inline-specifier
# -readability-redundant-member-init
# -readability-reference-to-constructed-temporary
```

Under [Tools/Linter](https://github.com/BLAST-WarpX/warpx/blob/development/Tools/Linter), the script `runClangTidy.sh` can be used to run the clang-tidy linter locally.

### clang-tidy local run script

```bash
#!/usr/bin/env bash
#
# Copyright 2024 Luca Fedeli
#
# This file is part of WarpX.
#

# This script is a developer's tool to perform the
# checks done by the clang-tidy CI test locally.
#
# Note: this script is only tested on Linux

echo "============================================="
echo
echo "This script is a developer's tool to perform the"
echo "checks done by the clang-tidy CI test locally"
echo "_____________________________________________"

# Check source dir
REPO_DIR=$(cd $(dirname ${BASH_SOURCE})/../../ && pwd)
echo
echo "Your current source directory is: ${REPO_DIR}"
echo "_____________________________________________"

# Set number of jobs to use for compilation
PARALLEL="${WARPX_TOOLS_LINTER_PARALLEL:-4}"
echo
echo "${PARALLEL} jobs will be used for compilation."
echo "This can be overridden by setting the environment"
echo "variable WARPX_TOOLS_LINTER_PARALLEL, e.g.: "
echo
echo "$ export WARPX_TOOLS_LINTER_PARALLEL=8"
echo "$ ./Tools/Linter/runClangTidy.sh"
echo "_____________________________________________"

# Check clang version
export CC="${CLANG:-"clang"}"
export CXX="${CLANGXX:-"clang++"}"
export CTIDY="${CLANGTIDY:-"clang-tidy"}"
echo
echo "The following versions of the clang compiler and"
echo "of the clang-tidy linter will be used:"
echo
echo "clang version:"
which ${CC}
${CC} --version
echo
echo "clang++ version:"
which ${CXX}
${CXX} --version
echo
echo "clang-tidy version:"
which ${CTIDY}
${CTIDY} --version
echo
echo "This can be overridden by setting the environment"
echo "variables CLANG, CLANGXX, and CLANGTIDY e.g.: "
echo "$ export CLANG=clang-19"
echo "$ export CLANGXX=clang++-19"
echo "$ export CLANGTIDY=clang-tidy-19"
echo "$ ./Tools/Linter/runClangTidy.sh"
echo
echo "******************************************************"
echo "* Warning: clang v19 is currently used in CI tests.  *"
echo "* It is therefore recommended to use this version.   *"
echo "* Otherwise, a newer version may find issues not     *"
echo "* currently covered by CI tests while older versions *"
echo "* may not find all the issues.                       *"
echo "******************************************************"
echo "_____________________________________________"

# Prepare clang-tidy wrapper
echo
echo "Prepare clang-tidy wrapper"
echo "The following wrapper ensures that only source files"
echo "in WarpX/Source/* are actually processed by clang-tidy"
echo
cat > ${REPO_DIR}/clang_tidy_wrapper << EOF
#!/bin/bash
REGEX="[a-z_A-Z0-9\/]*WarpX\/Source[a-z_A-Z0-9\/]+.cpp"
if [[ \$4 =~ \$REGEX ]];then
  ${CTIDY} \$@
fi
EOF
chmod +x ${REPO_DIR}/clang_tidy_wrapper
echo "clang_tidy_wrapper: "
cat ${REPO_DIR}/clang_tidy_wrapper
echo "_____________________________________________"

# Compile Warpx using clang-tidy
echo
echo "*******************************************"
echo "* Compile Warpx using clang-tidy          *"
echo "* Please ensure that all the dependencies *"
echo "* required to compile WarpX are met       *"
echo "*******************************************"
echo

rm -rf ${REPO_DIR}/build_clang_tidy

cmake -S ${REPO_DIR} -B ${REPO_DIR}/build_clang_tidy \
  -DCMAKE_CXX_CLANG_TIDY="${REPO_DIR}/clang_tidy_wrapper;--system-headers=0;--config-file=${REPO_DIR}/.clang-tidy" \
  -DCMAKE_VERBOSE_MAKEFILE=ON  \
  -DWarpX_DIMS="1;2;3;RZ"      \
  -DWarpX_MPI=ON               \
  -DWarpX_COMPUTE=OMP          \
  -DWarpX_FFT=ON               \
  -DWarpX_QED=ON               \
  -DWarpX_QED_TABLE_GEN=ON     \
  -DWarpX_OPENPMD=ON           \
  -DWarpX_PRECISION=SINGLE

cmake --build ${REPO_DIR}/build_clang_tidy -j ${PARALLEL} 2> ${REPO_DIR}/build_clang_tidy/clang-tidy.log

cat ${REPO_DIR}/build_clang_tidy/clang-tidy.log
echo
echo "============================================="
```

It is a prerequisite that WarpX is compiled following the instructions that you find in our [Users](../install/users.md#install-methods-cmake) or [Developers](../install/cmake.md#install-build-cmake) sections.

The script generates a wrapper to ensure that clang-tidy is only applied to WarpX source files and compiles WarpX in 1D, 2D, 3D, and RZ geometry, using such wrapper.

By default WarpX is compiled in single precision with PSATD solver, QED module, QED table generator and embedded boundary in order to ensure broader coverage with the clang-tidy tool.

Few optional environment variables can be set to tune the behavior of the script:

* `WARPX_TOOLS_LINTER_PARALLEL`: set the number of cores used for compilation;
* `CLANG`, `CLANGXX`, and `CLANGTIDY`: set the version of the compiler and the linter.

For continuous integration we currently use clang version 19 and it is recommended to use this version locally as well.
A newer version may find issues not currently covered by CI tests (checks are opt-in), while older versions may not find all the issues.

Here’s an example of how to run the script after setting the appropriate environment variables:

```bash
export WARPX_TOOLS_LINTER_PARALLEL=12
export CLANG=clang-19
export CLANGXX=clang++-19
export CLANGTIDY=clang-tidy-19

./Tools/Linter/runClangTidy.sh
```
