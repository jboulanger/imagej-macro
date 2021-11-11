#@File (label="Input file") local_path
#@String (label="Username", value="username") username
#@String (label="Host", value="hex") hostname
#@File (label="Template file") local_template_path
#@File (label="Local path",description="Local mounting point of the network share", style="directory") local_share
#@String (label="Remote path",description="Remote mounting point of the network share", value="/cephfs/") remote_share

/*
 * Launch a slurm job with the file filename as input using a template script
 * 
 * Jerome Boulanger 2021 
 */

template_name = File.getName(local_template_path);
remote_path = replace(local_path, local_share, remote_share);
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

