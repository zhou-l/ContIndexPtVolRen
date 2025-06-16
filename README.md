# Continuous Indexed Points for Multivariate Volume Visualization
The project includes the visualization tool for exporing processed multivariate volumes, and the Matlab code that preprocesses the multivariate volumes and indexed points.

More information can be found in the paper "Continuous Indexed Points for Multivariate Volume Visualization" by Liang Zhou, Xinyi Gou, and Daniel Weiskopf. (The CVM journal).
If you use our approach, please cite the following paper.
```
@ARTICLE{Zhou2025cvm,
  author={Zhou, Liang and Gou, Xinyi and Weiskopf, Daniel},
  title={Continuous Indexed Points for Multivariate Volume Visualization}, 
  journal={Computational Visual Media}, 
  year={},
  volume={},
  number={},
  pages={}
  }
````
![Continuous indexed points for volume visualization](/images/cvm25.png)

# Installations
To install the visualization tool, unzip "bin.zip" and "opencv_world450.zip" in the vistool directory. 

Place the unzipped "opencv_world450.dll" to the unzipped bin directory. 

[Binary files are placed separately due to the file size limit of github.]

# How to Run?

Run the visualization tool using "QImprovedContPlots.exe" in the unzipped bin directory.

# Example Datasets
A number of synthetic datasets are included in the "data" directory. 
The user can open a ".ipcproj" file through "File"->"Open Project File" to load one dataset at a time. 
You will have to close the application and restart it to load another dataset.
![The synthetic data of 'cvmappendix_synthetic3D.ipcproj' (the dataset of Fig. 14 of the paper ) ](/images/cvm_synth.png)


# Documentation
Find the documentation "Instruction.pdf" in the vistool directory.

# The Matlab Code
The Matlab code concern the preprocessing of multivariate volumes and converting discrete indexed points to continuous.
