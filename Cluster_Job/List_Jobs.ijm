#@String (label="Username", value="username", description="Username on the host") username
#@String (label="Host", value="hex", description="hostname to which command are send") hostname
#@Boolean (label="All users", value=false) allusers
/*
 * List cluster slurm jobs 
 */

if (isOpen("Results")) {selectWindow("Results");run("Close");}
if (allusers) {
	ret = exec("ssh", username+"@"+hostname, "squeue");	
} else {
	ret = exec("ssh", username+"@"+hostname, "squeue -u "+username);
}

lines = split(ret,"\n");
if (lines.length == 1) {
	print("No jobs are running");
}
for (i = 1; i < lines.length; i++) {
	elem = split(lines[i]," ");	
	setResult("JOBID",i-1,elem[0]);
	setResult("PARTITION",i-1,elem[1]);
	setResult("NAME",i-1,elem[2]);
	setResult("USER",i-1,elem[3]);
	setResult("STATUS",i-1,elem[4]);
	setResult("TIME",i-1,elem[5]);
}


