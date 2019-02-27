# Define a new class - finding
finding <- setClass("finding", slots=c(rationale="character",
                               libs="character",
                               data="list",
                               vis = "list",
                               conclusion ="character",
                               discussion = "character"))

# Example
# F2<- finding(rationale = "Testing the new class",
#              libs = "purrr",
#              data = list(test="a"),
#              vis = list(plot="1"),
#              conclusion = "This is dumb test",
#              discussion = "test could be ok?")
#
# F2@rationale <- "This finding is obwiously going to be based on pure hypothetical shooting"
