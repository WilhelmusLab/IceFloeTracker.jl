using TestItemRunner
@run_package_tests filter = ti -> (:e2e in ti.tags)