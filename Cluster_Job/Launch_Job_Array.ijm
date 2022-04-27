#@String (label="Username", value="username", description="Username on the host") username
#@String (label="Host", value="hex", description="hostname to which command are send") hostname
#@File (label="Job file", description="Script to send to the host") local_template_path
#@File (label="Local share",description="Local mounting point of the network share", style="directory") local_share
#@String (label="Remote share",description="Remote mounting point of the network share", value="/cephfs/") remote_share
#@String (label="Data root location",description="Local path to the dataset", value="/cephfs/") local_data_path


/*
 * Launch job array on a slurm managed cluster
 * 
 * The list of file is defined by the opened table.
 * To create a list of file use the macro File_Conversion/Parse_folder.ijm
 * 
 * Prepare before running the macro:
 * - map remote drives 
 * - copy data to remote host
 * - list jobs to run in a csv table and open it in imagej
 * - create a job script unless you want to reuse an existing one 
 * 
 * The macro will
 * - Creates a job folder in "Remote share" and copy the script "template file" and the "opened table"
 * - Convert file to unix file ending
 * - send a sbatch job.sh via ssh
 * 
 * Jerome Boulanger 2021
 */

delay = 500; // artificial delays that prevent denial of service

getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
table_name = File.getNameWithoutExtension(Table.title);
Table.rename(Table.title, table_name);
table_name += "-"+year+""+month+""+dayOfMonth+""+hour+""+minute+""+second+".csv";

template_name = File.getName(local_template_path);
remote_jobs_dir = remote_share + "/jobs";
local_jobs_dir = local_share + File.separator + "jobs";
remote_template_path = remote_jobs_dir + '/' + template_name;
local_template_cpy_path = local_jobs_dir + File.separator + template_name;
local_path = local_jobs_dir + File.separator + table_name;
remote_path = remote_jobs_dir + "/" + table_name;
remote_data_path = replace(convert_slash(local_data_path), convert_slash(local_share), remote_share);

// create a job folder if needed
if (File.exists(local_jobs_dir) != 1) {
	print("Create a jobs folder in " + local_share);
	File.makeDirectory(local_jobs_dir);
} else {
	print("Jobs folder already present in " + local_share);
}

// save the table
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
print("Remote data path: "+ remote_data_path);
print("Remote path: "+ remote_path);
print("Remote job dir: "+ remote_jobs_dir);
ret = exec("ssh", username+"@"+hostname, "sbatch", "--chdir", remote_jobs_dir, "--array=1-"+n, template_name, remote_path, remote_data_path);
wait(delay);
print(ret);
print("Ready\n");




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