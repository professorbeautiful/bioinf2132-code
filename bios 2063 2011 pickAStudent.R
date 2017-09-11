bios2063names = scan(what=character(0), sep="=")
Thomas Blaze
Olive Buhule
Ian Cook
Jiang Yingda
Modhurima Moitra
Megan Olson Hunt
John Pleis 
WenjingQi
Emmanuel Sampene
Zhang Yang
Xia Feihong
Ye Ye
Ling Yun
Amy Kao
Liu Qing
Tanya Bekhuis
Andrea Aldrich
Chimwemwe
Ching-Wen Lee
Drew Michanowicz
Chen Yi-Fan
Kevin McDade

length(bios2063names)   ###  22
##  Also Kevin McDade, 


#Lopa,Samia H.
#Singhabahu,Dilrukshika Manori

pickAStudent = function() {
	if( identical(find("bios2063chosen"), character(0))) bios2063chosen = character(0)
	if( length(bios2063names) == length(bios2063chosen)) bios2063chosen=character(0)
	name = sample(setdiff(y=bios2063chosen, bios2063names), 1)
	bios2063chosen <<- c(bios2063chosen, name)
	name
}

pickAStudent()