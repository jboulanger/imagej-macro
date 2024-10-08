#@File (label="Input file",description="Select the file to process") input_file
#@String(label="patch",choices={"5x5","7x7","9x9","11x11","13x13","3x3x3","5x5x3","5x5x5","5x5x1x3","7x7x1x3","9x9x1x3"}) patch_str
#@String(label="mode",choices={"2D","3D","3D+time","2D+time"}) mode_str
#@Float(label="smoothing",value=2,description="converted to a pvalue =10^-smoothing") smoothing
#@Integer(label="iterations",value=4) niter
#@Float(label="sensor gain",value=-1, description="if negative, noise parameter are estimated") gain
#@Float(label="sensor offset",value=0, description="offset") offset
#@Float(label="readout noise",value=0, description="readout in gray levels") readout
#@String (label="Username", value="username", description="Username on the host") username
#@String (label="Host", value="hex", description="hostname to which command are send") hostname
#@File (label="Local share",description="Local mounting point of the network share", style="directory") local_share
#@String (label="Remote share",description="Remote mounting point of the network share", value="/cephfs/") remote_share
#@String(label="Action",choices={"Process","List Jobs","Cancel Jobs"}) action


/*
 * Launch a slurm job with the file filename as input using a template script
 * 
 * Jerome Boulanger 2021 
 */

delay = 500; // artificial delays that prevent denial of service 
remote_path = replace(convert_slash(input_file), convert_slash(local_share), remote_share);
template_name = "denoise.sh";
remote_jobs_dir = remote_share + "/jobs";
local_jobs_dir = local_share + File.separator + "jobs";
remote_template_path = remote_jobs_dir + '/' + template_name;
local_template_path = local_jobs_dir + File.separator + template_name;
patch = parsePatchStr(patch_str);
mode = parseModeStr(mode_str);
noise = parseNoise(gain,offset,readout);
pval = pow(10,-smoothing);
var jobid = 0;	
// create a job folder if needed
if (File.exists(local_jobs_dir) != 1) {
	print("Create a jobs folder in " + local_share);
	File.makeDirectory(local_jobs_dir);
} else {
	print("Jobs folder already present in " + local_share);
}

if (matches(action, "Process")) {
	process();
} else if (matches(action, "List Jobs")) {
	listjobs();
} else if (matches(action, "Cancel Jobs")) {
	canceljob();
}

function process() {
	print("Processing");
	str = "#!/bin/tcsh\n";
	str += "#SBATCH --mail-type=END,FAIL\n";
	str += "#SBATCH --nodes=1\n";
	str += "#SBATCH --ntasks-per-node=1\n";
	str += "#SBATCH --cpus-per-task=112\n";
	str += "#SBATCH --time=05:00:00\n";
	str += "set ifile = \"$1\"\n";
	str += "set ofile = `echo $ifile | sed s/\\.tif/_denoised\\.tif/g` \n";
	str += "/lmb/home/jeromeb/bin/ndsafir -i \"$ifile\" -o \"$ofile\" "+patch+mode+" -iter "+niter+" -pval "+pval + " " + noise + "\n";
	File.saveString(str, local_template_path);
	ret = exec("ssh", username+"@"+hostname, "dos2unix", remote_template_path);		
	print(ret);
	ret = exec("ssh", username+"@"+hostname, "sbatch", "--chdir", remote_jobs_dir, template_name, "\""+remote_path+"\"");
	print(ret);
	parts  = split(ret," ");
	jobid = parseInt(parts[parts.length - 1]);
	print("Job ID:");
	print(jobid);	
	print("Job log file:");
	print(local_jobs_dir + File.separator + "slurm-"+jobid+".out");
	
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


function listjobs() {
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
		setResult("LOG",i-1,local_jobs_dir + File.separator + "slurm-"+elem[0]+".out");
	}
}

function canceljob() {
	jobid = getNumber("Job ID", jobid);
	ok = getBoolean("Do you want to cancel job " + jobid, "Yes", "No");
	if (ok) {
		ret = exec("ssh", username+"@"+hostname, "scancel "+jobid);
	}
}

//choices={"5x5","7x7","9x9","3x3x3","5x5x3","5x5x5"})
function parsePatchStr(str) {
	if (str=="5x5") {
		patch = newArray(2,2,0,0,0);
	} else if (str=="7x7") {
		patch = newArray(3,3,0,0,0);
	} else if (str=="9x9") {
		patch = newArray(4,4,0,0,0);
	} else if (str=="11x11") {
		patch = newArray(5,5,0,0,0);
	} else if (str=="13x13") {
		patch = newArray(5,5,0,0,0);
	} else if (str=="7x7") {
		patch = newArray(1,1,1,0,0);
	} else if (str=="5x5x3") {
		patch = newArray(2,2,1,0,0);
	} else if (str=="5x5x5") {
		patch = newArray(2,2,2,0,0);
	} else if (str=="5x5x1x3") {
		patch = newArray(2,2,0,0,1);
	} else if (str=="7x7x1x3") {
		patch = newArray(3,3,0,0,1);
	} else if (str=="9x9x1x3") {
		patch = newArray(5,5,0,0,1);
	} else {
		print("Undefined patch shape" + str);
		patch = newArray(4,4,0,0,0);
	}	
	return " -patch "+patch[0]+","+patch[1]+","+patch[2]+","+patch[3]+","+patch[4];
}

// {"2D","3D","3D+time","2D+time"}
function parseModeStr(str) {
	if (str=="2D") {
		mode = " -3D 0 -time 0";
	} else if (str=="3D") {
		mode = " -3D 1 -time 0";
	} else if (str=="3D+time") {
		mode = " -3D 1 -time 1";
	} else if (str=="2D+time") {
		mode = " -3D 0 -time 1";
	}
	return mode;
}

function parseNoise(_gain,_offset,_readout) {
	if (_gain > 0) {
		return " -noise " + _gain + "," + _offset + "," + _readout;
	} else {
		return "";
	}
}