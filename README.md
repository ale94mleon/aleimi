Here I need to provied a correct documentation.

Coded with  MOPAC2016 (Version: 21.329L)

# Setting the environment
 * Install MOPAC2016 (Version: 21.329L). Go to the web and look for the installation
 * Now you will need a lot of different libraries to work with. Must of them are reached through conda and others through pip. Therefore you have to install conda in order to retrieve the corresponded modules.

```
conda create -n aleimi
conda activate aleimi
conda config --add channels conda-forge
conda install rdkit cython
conda install openbabel
pip install rmsd numpy pandas 
```

I think that the rest of the library are standard python libraries. If you get some missing library, install it.
# Package
To clone the repository you must have permission.
```
mkdir ~/GITLAB
cd GITLAB
git clone git@gitlab.com:md_tools/aleimi.git
```
In the `~/.bashrc` add the following lines
```
alias aleimi='conda activate aleimi'
export PYTHONPATH=$PYTHONPATH:/home/$USER/GITLAB/
export PATH=$PATH:/home/$USER/GITLAB/aleimi
```
Now you should be able to use the command line without problems:

`aleimi -h`

And get:
```
positional arguments:
  suppl                 The path were the molecule(s) is(are)

optional arguments:
  -h, --help            show this help message and exit
  -p PARAMS, --params PARAMS
                        Parameters to run ALEIMI
```
# How to use it
It should be simple. Let's assume that you have in `~/my_project/mols` a batch of molecules (and nothing more); it doesn't matter the format (.pdb, .mol, .pdbqt, .smi, etc.), `aleimi` will try to handled. And also you would like to pass some parameters because you don't like the default ones. Then you have a text file like `~/my_project/my_awesome_params.txt`.

Then it should as easy as call
```
aleimi ~/my_project/mols -p ~/my_project/my_awesome_params.txt
```
And all the magic will be done!

## -p, --params option
This option will control the parameters of the functions:
```
aleimi.tools.mopac
aleimi.boltzmann.main
aleimi.extractor.main
```
Here I had to add more description about it.
