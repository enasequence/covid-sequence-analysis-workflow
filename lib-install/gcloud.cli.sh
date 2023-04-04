echo "Installing gcloud cli"
curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-421.0.0-linux-x86_64.tar.gz
tar -xf ./google-cloud-cli-421.0.0-linux-x86_64.tar.gz
mv google-cloud-sdk /opt/google-cloud-sdk
ln -s /opt/google-cloud-sdk/bin/* /bin/
rm google-cloud-cli-422.0.0-linux-x86_64.tar.gz
# mv ./google-cloud-sdk ~/
# rm ./google-cloud-cli-421.0.0-linux-x86_64.tar.gz
# $HOME/google-cloud-sdk/install.sh
# $HOME/google-cloud-sdk/bin/gcloud init
gcloud auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
gcloud config set project "prj-int-dev-covid19-nf-gls"
echo "Installed gcloud cli"