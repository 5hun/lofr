library(RUnit)
library(lofr)
test.suite <- defineTestSuite("all", 
    dirs=file.path(getwd(), "runit"))
runTestSuite(test.suite)
