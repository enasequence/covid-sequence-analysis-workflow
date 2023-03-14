sudo yum update -y
sudo amazon-linux-extras install docker
sudo yum install docker
sudo service docker start

sudo groupadd docker
sudo usermod -aG docker $(whoami)
#reboot instance after add user & rerun `sudo service docker start`
sudo service docker start 
