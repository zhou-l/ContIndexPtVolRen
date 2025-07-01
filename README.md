# Continuous Indexed Points for Multivariate Volume Visualization
The project includes the visualization tool for exporing processed multivariate volumes, and the Matlab code that preprocesses the multivariate volumes and indexed points.

More information can be found in the paper "Continuous Indexed Points for Multivariate Volume Visualization" by Liang Zhou, Xinyi Gou, and Daniel Weiskopf. 
If you use our approach, please cite as follows.
```
@ARTICLE{Zhou2025cvm,
  author={Zhou, Liang and Gou, Xinyi and Weiskopf, Daniel},
  title={Continuous Indexed Points for Multivariate Volume Visualization}, 
  journal={Computational Visual Media}, 
  year={in press},
  volume={},
  number={},
  pages={},
  doi={10.26599/CVM.2025.9450496 }
  }
````
The arXiv version of the paper can be found at: [arXiv:2506.19400](https://arxiv.org/abs/2506.19400).

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


# Documentation and notes
Find the documentation "Instruction.pdf" in the vistool directory.

In addition, key buttons are used to control the viewers.

T-Key: press "T" with the volume renderer being the current view to recompile shaders.

E-Key: press "E" with the volume renderer being the current view to cycle through compare modes for spatial and data domain indexed points computation.

F-Key: press "F" with the parallel coordinates being the current view to cycle through compare modes for spatial and data domain indexed points computation. 

The associations of paper figures and datasets:
Fig1-DTI (1-flat); Fig3-fig3_synth3D, fig3_planeSynth; Fig4-Isabel; Fig5-CT tooth; Fig6,7-BraTs MRI; Fig 8-fig3_synth3D, fig3_planeSynth; Fig 9-Isabel; Fig 10-DTI(2-flat); Fig12,13,14 - synthVol3D; Fig15-CT tooth (val vs. gradient mag); Fig16 - BraTs (T1 vs. gradient mag). 

# Preprocessing and the Matlab Code
The indexed volumes can be computed in the visualization software when loading the original data volumes and enable the correlation computation in the spatial neighborhood. For correct results, please set the sampling rate for computation to 1 (100%). This can also be done with the Matlab code which is slower.
While the computation gives the indexed point volumes, it also yield a large .txt file that contains the 2D positions of indexed points on the image plane of the parallel coordinates that can be converted into a continuous density representation.

The Matlab code concern the preprocessing of multivariate volumes and converting discrete indexed points to a continuous density.
