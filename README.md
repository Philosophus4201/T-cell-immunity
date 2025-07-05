# T cell cross-recactivity predictor pipeline


## 1) Installaton tutorial


To work, we recommend creating two conda environments named "IEPAPI" and "BATMAN", but if you wish, you can specify the path to different Python interpreters directly in the script. Python environments must use python==3.7.9 and python==3.11.11, as well as the dependencies listed in "IEPAPI_dependencies.txt" and "BATMAN_dependencies.txt". To create all the necessary environments, you can simply run this code:
```console
cd ./T-cell-immunity
conda env create -f BATMAN_env.yml
conda env create -f IEPAPI_env.yml
```

The script must be run from T-cell-immunity/ directory and IEPAPI need to be installed. You can clone IEPAPI repository using this command:
```console
git clone https://github.com/ddd9898/IEPAPI.git
```
## 1) Run code


Now you can use Cross_reaact.py as module in your scrips. For convenience, we have provided a notebook (Example_notebook.ipynb) that shows how you can run the code on your data. You must run the notebook in the BATMAN environment. 

