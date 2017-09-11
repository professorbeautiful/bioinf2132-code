bios2063names = scan(what=character(0), sep="=")
Jon
Amie
Raf
Caleb
Joyeeta
Ying
Kevin


length(bios2063names)   ###  

pickAStudent = function() {
	if( identical(find("bios2063chosen"), character(0))) bios2063chosen = character(0)
	if( length(bios2063names) == length(bios2063chosen)) bios2063chosen=character(0)
	name = sample(setdiff(y=bios2063chosen, bios2063names), 1)
	bios2063chosen <<- c(bios2063chosen, name)
	name
}

pickAStudent()