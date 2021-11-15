# SLURM clusters jobs

This macro helps running SLURM jobs https://slurm.schedmd.com on a cluster with out installing additional software on a pc other than ImageJ. The macro uses predefined template scripts taking as input the filename. See template0.sh as an example. The script is responsible for defining the task to run and will be fed to the command ``sbatch`` on the host to initialize jobs. The template file will be copied to a 'job' folder in the path. 

Typically, the following command is issued to create a job:
```
ssh username@hostname sbatch template.sh filename
```


## Configuration
### Password less configuration
To enable the connection via ssh without password, ssh keys need to be generated and added to the file .ssh/authorized_keys on the host. In the following replace username and hostname by your username on the remote server and hostname by the name of the remote server.
1. Generate the keys using the following command:
```bash
ssh-keygen
```
2. Add the key to the list of recognized client on the host:
  - On unix systems (Mac/Linux) using ssh-copy-id:
    ```bash
    ssh-copy-id -i ~/.ssh/id_rsa.pub username@hostname
    ```
   - On windows system using the following command in the terminal (powershell)
     ```bash
     type %USERPROFILE%\.ssh\id_rsa.pub | ssh username@hostname "cat >> .ssh/authorized_keys"
     ```
3. Finally, test the access using
   ```
   ssh username@hostname
   ```
   no password should be required for the connection. You can also check the added line in the file .ssh/authorized_keys

### Creating template jobs
In order to define what processing to apply to the input file, one need to write a template script which can be reused for several files. The template script is written in the shell (command line interpreter, bash or tcsh for example) running used on the cluster. It can include directives for the slurm manager that are equivalent to the command line arguments for the command sbatch https://slurm.schedmd.com/sbatch.html. For example:
```
#SBATCH --partition=cpu
#SBATCH --mail-type=END
```

## SLURM job
You can run the macro in normal mode or in batch mode.
1. Map or connect the network drive where data are stored and accessible by the cluster.
2. Open the macro and press "run" or "batch"
3. Select the file you want to process that are already on the network drive
4. Set the your username and the name of the server you want to connect to.
5. Select the path to the network drive (on windows this will be the letter of the drive, for example Z:\)
7. Write the path to the drive as accessible from the cluster. This should be starting by / (absolute path)
8. Press ok

## JOB array
It might be more convenient to use a job array instead of a job. To define the array, a table is used with each line defining a job. This enable to define filenames but also various parameters that can be used by the macro. Use [Parse_Folder.ijm](https://raw.githubusercontent.com/jboulanger/imagej-macro/main/File_Conversion/Parse_Folders.ijm) to create a list of files you want to process, optionally edit this file in excel for example and save it as csv. Open the table and run the macro. We assume again that the network share is accessible both on the local computer and the remote computer.
The parameters are:
- Username: your username on the cluster
- Host : the name of the server to connect to.
- Template File: the script that will run on the host. See templates-array0.sh as an example of script
- Local share: local mounting point of the network share (eg X:\ on windows)
- Remote share: Remote mounting point of the network share (eg /bignetworkdrive/$USER/)
- Data root location: the root (absolute path) to the data on the locally mapped drive. (X:\data)




