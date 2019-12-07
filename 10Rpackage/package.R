setwd("F:/Rkeyan")
library(devtools)
create_package("balala")
#ctrl+shift+alt+'R':write annotation for functions
setwd("F:/Rkeyan/balala")
document()
setwd("F:/Rkeyan")
install("balala")