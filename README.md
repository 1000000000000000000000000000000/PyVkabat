# PyVkabat
## About PyVkabat
PyVkabat is a bioinformatics/computational chemistry tool used to calculated secondary structure variability at each residue position in a protein of interest. To reduce run time, each of the requests are made in parallel.

### Vkabat equation
Vkabat = k * N / n1

where:
<ul><li>k is the number of secondary structure classes predicted for a residue</li>
    <li>N is the total number of predictions made for a given residue</li>
    <li>n1 is how often the most predicted secondary structure is observed for a residue</li></ul>
    
### PyVkabat web servers used
PyVkabat submits the user's protein sequence of interest to 15 different algorithms:
<ul>
  <li>gor1</li>
  <li>gor3</li>
  <li>dpm</li>
  <li>sopm</li>
  <li>predator</li>
  <li>hnn</li>
  <li>dsc</li>
  <li>mlrc</li>
  <li>yaspin</li>
  <li>jpred</li>
  <li>phd</li>
  <li>prof</li>
  <li>sspro</li>
  <li>jnet</li>
  <li>psipred</li>
</ul>

## Installation Guide
### Prerequisite Software
<ul>
<li>Miniconda or Annaconda</li>
<li>Windows or Linux or Mac</li>
</ul>

### Instructions
After downloading the repository and successfully installing the prerequisite software on your machine, navigate to the project directory. In the command prompt, enter the following:
```
conda create -n PyVkabat python=3.10 -y
```

Next, enter the following:
```
conda activate PyVkabat
```

Now, enter the following (make sure you are in the directory where pyvkabat_requirements.txt is located):
```
pip install -r pyvkabat_requirements.txt
```

If the previous command did not return any errors, you have successfully installed python and the required python packages to run the program. Now just clean up by entering the following command:
```
conda deactivate
```

That's it, you're done!

## Running the program
To run the program, open the command prompt. You can navigate to the project folder by entering the following:
```
cd <PATH TO PROJECT FOLDER>
```

Now that you are in the project folder, enter the following command:
```
conda activate PyVkabat && python ./PyVkabat.py
```
