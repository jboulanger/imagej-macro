# Cellpose in Fiji

This macro runs the [cellpose](https://www.nature.com/articles/s41592-020-01018-x) segmentation algorithm  on the current image and load the regions as regions in the ROI manager.

It requires to install cellpose first, then to point to the installation path.

## Installation

On windows, install pytorch and cellpose in an environment. Check the version of your cuda drivers using `nvidia-smi`. Installing torchvision and torchaudio at the same time than pytorch made conda unable to resolve dependencies. It might be needed to install qt and pyqt separately to have access to the graphical interface.

```
conda create -n cellpose python=3.9
conda activate cellpose
conda install pytorch pytorch-cuda=11.7 -c pytorch -c nvidia
conda install -c conda-forge qt pyqt
pip install cellpose
```

