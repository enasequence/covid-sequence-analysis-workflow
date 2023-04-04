echo "Installing gcloud cli"
export $(grep -v '^#' .env | xargs)
curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-421.0.0-linux-x86_64.tar.gz
tar -xf ./google-cloud-cli-421.0.0-linux-x86_64.tar.gz
mv ./google-cloud-sdk ~/
rm ./google-cloud-cli-421.0.0-linux-x86_64.tar.gz
$HOME/google-cloud-sdk/install.sh
$HOME/google-cloud-sdk/bin/gcloud init
echo "Installed gcloud cli"