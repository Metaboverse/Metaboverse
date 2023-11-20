VERSION=0.11.1

# Install ubuntu programs
apt-get install wget jq

# Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -O /tmp/miniconda3.sh
bash /tmp/miniconda3.sh -b -p ~/miniconda3
rm /tmp/miniconda3.sh
conda init bash
source ~/.bashrc

# Install node 
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.34.0/install.sh | bash
source ~/.bashrc
nvm install node

# Run build
git clone https://github.com/Metaboverse/Metaboverse.git
cd Metaboverse/resources
bash ./__build-nix__.sh $VERSION --build

