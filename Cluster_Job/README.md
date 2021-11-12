# Clusters jobs

This macro helps running jobs on a cluster with out installing additional software on a pc other than
ImageJ.

## Configuration
## Password less configuration
To enable the connection via ssh without password, ssh keys need to be generated and copied to the host.
- first generate the keys using the following command:
``ssh-keygen``
- then copy the key to the host:
  - On unix systems (Mac/Linux) using ssh-copy-id:
    ``ssh-copy-id -i ~/.ssh/id_rsa.pub username@hostname``
   - On windows system using the following command in the terminal (powershell)
     ``type %USERPROFILE%\.ssh\id_rsa.pub | ssh username@hostname "cat >> .ssh/authorized_keys"
- now test the access using
-   ``ssh username@hostname``

## Creating template jobs
