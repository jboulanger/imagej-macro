#@String (label="Username", value="username", description="Username on the host") username
#@String (label="Host", value="hex", description="hostname to which command are send") hostname
#@String (label="Job ID", value="hex", description="job id") jobid
/*
 * Cancel a slurm job
 */
 
ret = exec("ssh", username+"@"+hostname, "scancel "+jobid);

