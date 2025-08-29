using TestItemRunner
@run_package_tests filter = ti -> (:smoke in ti.tags)