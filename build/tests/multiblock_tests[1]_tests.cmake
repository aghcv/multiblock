add_test([=[IOUtils.CanReadSimpleSTL]=]  /Users/aghorban/code/multiblock/build/tests/multiblock_tests [==[--gtest_filter=IOUtils.CanReadSimpleSTL]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[IOUtils.CanReadSimpleSTL]=]  PROPERTIES WORKING_DIRECTORY /Users/aghorban/code/multiblock/build/tests SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  multiblock_tests_TESTS IOUtils.CanReadSimpleSTL)
