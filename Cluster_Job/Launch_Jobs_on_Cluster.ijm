#@File (label="Input file",description="Select the file to process") local_path
#@String (label="Username", value="username", description="Username on the host") username
#@String (label="Host", value="hex", description="hostname to which command are send") hostname
#@File (label="Template file", description="Script to send to the host") local_template_path
#@File (label="Local path",description="Local mounting point of the network share", style="directory") local_share
#@String (label="Remote path",description="Remote mounting point of the network share", value="/cephfs/") remote_share

/*
 * Launch a slurm job with the file filename as input using a template script
 * 
 * Jerome Boulanger 2021 
 */

remote_path = replace(convert_slash(local_path), convert_slash(local_share), remote_share);
template_name = File.getName(local_template_path);
remote_jobs_dir = remote_share + "/jobs";
remote_template_path = remote_jobs_dir + '/' + template_name;
print("Create a job folder on " + hostname);
ret = exec("ssh", username+"@"+hostname, "mkdir", "-p", remote_jobs_dir);
print(ret);
wait(200);

print("Send template " + template_name + " file to " + hostname + ".");
ret = exec("scp", local_template_path, username+"@"+hostname+":", remote_template_path);
print(ret);
wait(200);

print("Start job " + File.getName(local_template_path) + " on " + hostname);
ret = exec("ssh", username+"@"+hostname, "sbatch", "--chdir", remote_jobs_dir, template_name, "\""+remote_path+"\"");
print(ret);
print("");


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