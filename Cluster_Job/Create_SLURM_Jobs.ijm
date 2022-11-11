#@File   (label="Input file", description="Select the file to process") input_file
#@String (label="Username", value="username", description="Username on the host") username
#@String (label="Host", value="hex", description="hostname to which command are send") hostname
#@File   (label="Template file", description="Script to send to the host") local_template_path
#@File   (label="Local share",description="Local mounting point of the network share", style="directory") local_share
#@String (label="Remote share",description="Remote mounting point of the network share", value="/cephfs/") remote_share
#@String (label="Action", choices={"Submit job","List","Cancel"}) action

/*
 * Single script to interact with a SLURM cluster
 * 
 */

delay = 500; // artificial delays 
remote_path = replace(convert_slash(input_file), convert_slash(local_share), remote_share);
template_name = File.getName(local_template_path);
remote_jobs_dir = remote_share + "/jobs";
local_jobs_dir = local_share + File.separator + "jobs";
remote_template_path = remote_jobs_dir + '/' + template_name;
local_template_cpy_path = local_jobs_dir + File.separator + template_name;


if (matches(action, "Submit job")) {
	processSingle();
} else if (matches(action, "Multiple jobs")) {
	processMulti();
} else if (matches(action, "List")) {
	listJobs();
} else if (matches(action, "Cancel")) {
	cancelJobs();
} 

function processSingle() {
	/*
	 * Process a single file
	 */
	 
	initJobsFolder();
		
	// copy the template file
	print("Copy template " + template_name + " file to " + local_template_cpy_path + ".");
	File.copy(local_template_path, local_template_cpy_path);
	
	// convert encoding of the file, otherwise sbatch will complain
	ret = exec("ssh", username+"@"+hostname, "dos2unix", remote_template_path);
	wait(delay);
	
	// start the job with a ssh command
	print("Start job " + template_name + " on " + hostname);
	ret = exec("ssh", username+"@"+hostname, "sbatch", "--chdir", remote_jobs_dir, template_name, "\""+remote_path+"\"");
	print(" >>" + ret);
	wait(delay);
}

function processMulti() {
	/*
	 * Process the file listed in a table
	 */
	 
	initJobsFolder();
	
	// rename the table
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	table_name = File.getNameWithoutExtension(Table.title);
	Table.rename(Table.title, table_name);
	table_name += "-"+year+""+month+""+dayOfMonth+""+hour+""+minute+""+second+".csv";
		
	// save the table
	local_path = local_jobs_dir + File.separator + table_name;
	remote_path = remote_jobs_dir + "/" + table_name;
	print("Save the table "+Table.title+" to file " + local_path);
	Table.save(local_path);
	n = Table.size;
	
	// copy the template file
	print("Copy template " + template_name + " file to " + local_template_cpy_path + ".");
	File.copy(local_template_path, local_template_cpy_path);
	
	// convert encoding of the file, otherwise sbatch will complain
	ret = exec("ssh", username+"@"+hostname, "dos2unix", remote_template_path);
	wait(delay);
	
	// start the job with a ssh command
	print("Start job " + template_name + " on " + hostname);
	ret = exec("ssh", username+"@"+hostname, "sbatch", "--chdir", remote_jobs_dir, template_name, "\""+remote_path+"\"");
	print(" >>" + ret);
	wait(delay);
}

function listJobs() {
	/*
	 * List the running jobs on the cluster calling squeue
	 */
	print("[ List jobs ]");	
	if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
	ret = exec("ssh", username+"@"+hostname, "squeue -u "+username);
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
}

function initJobsFolder() {	
	// create a job folder if needed
	if (File.exists(local_jobs_dir) != 1) {
		print("Create a jobs folder in " + local_share);
		File.makeDirectory(local_jobs_dir);
	} else {
		print("Jobs folder already present in " + local_share);
	}
}

function cancelJob() {
	/*
	 * Cancel job by runing scancel
	 */
	ret = exec("ssh", username+"@"+hostname, "squeue -u "+username);
	lines = split(ret,"\n");
	if (lines.length == 1) {
		print("No jobs are running");
	}
	ids = newArray(lines.length-1);
	for (i = 1; i < lines.length; i++) {
		elem = split(lines[i]," ");
		ids[i-1] = elem[0]
	}
	Dialog.create("Select jobs to cancel");
	for (i = 0; i < ids.length; i++) {
		Dialog.addCheckbox(ids[i], default);
	}
	Dialog.show();	
	str = "";
	for (i = 0; i < ids.length; i++) {
		if (Dialog.getCheckbox()) {
			str += ids[i] + " "
		}
	}	
	ret = exec("ssh", username+"@"+hostname, "scancel "+jobid);
	print(ret);
}

function convert_slash(src) {
	// convert windows file separator to unix file separator if needed
	if (File.separator != "/") {
		a = split(src,File.separator);
		dst = "";
		for (i=0;i<a.length-1;i++) {
			dst += a[i] + "/";
		}
		dst += a[a.length-1];
		return dst;
	} else {
		return src;
	}
}