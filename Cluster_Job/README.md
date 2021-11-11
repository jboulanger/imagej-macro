# Clusters jobs

This macro helps running jobs on a cluster with out installing additional software on a pc other than
ImageJ.

## Configuration
## Password less configuration
To enable the connection via ssh without password, ssh keys need to be generated and copied to the host.
- first generate the keys using the following command:
``sh ssh-keygen``
- then copy the key to the host:
``sh ssh-copy-id -i ~/.ssh/id_rsa.pub username@hostname``

## Creating template jobs
