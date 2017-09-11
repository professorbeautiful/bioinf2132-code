tempList = list()
tempChangeList = function() {
  tempList$a = 1
  cat("I changed to ")
  print(tempList)
}
tempChangeList()
tempList   ### tempList did not change. Local variable

tempList = list()
tempChangeList = function() tempList$a = 1
tempChangeList()
tempList   
### tempList did not change. Local variable
tempChangeList = function(myList=tempList) myList$a = 1
tempChangeList()
tempList   
### tempList did not change. Arg is call-by-value, not call-by-reference.


### Now use the global assignment operator
tempChangeList = function() tempList$a <<- 1 
tempChangeList()
tempList  # With the "<<-" operator, it does change.
# You can also you the function assign().

#### Now, the same for Environments.
tempEnv = new.env()
tempChangeEnv = function() tempEnv$a = 1
tempChangeEnv()
tempEnv$a     
### Aha, it DID change.  You don't need the global assignment.
tempChangeEnv = function(myEnv=tempEnv) myEnv$a = 2
tempChangeEnv()
tempEnv$a
### Aha, it changed again. Environments are ONLY references.

tempList2 = tempList
tempList2
tempList2$a = 2
tempList2
tempList


tempEnv2 = tempEnv
tempEnv2
get("a", envir=tempEnv2)
tempEnv2$a
tempEnv2$a = 4
tempEnv2$a
tempEnv$a

temp = function() {
     if(is.null(environment(temp)$x)) x = 0; 
     x<<-x+1; x
   }
temp()
temp()
temp()

