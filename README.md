# Feature Graph Learning for 3D Point Cloud Denoising

This repository is the official MATLAB implementation of the following paper:

Wei Hu, Xiang Gao, Gene Cheung, Zongming Guo, "Feature Graph Learning for 3D Point Cloud Denoising," accepted to IEEE Transactions on Signal Processing (TSP), March, 2020.

# Usage

This code is tested on `MATLAB R2018a` under Windows 10 x64 platform. You can also run the following command in the MATLAB command line on Linux/macOS:

```
MATLAB> main
```

**Note:** you can directly run the script from the MATLAB command window.

To support more datasets and different noise levels, you can run the `add_gaussian_noise` script to generate your own data and place them in the `models` folder.

# Reference

Please cite our paper if you use any part of the code from this repository:

```
@article{hu2019feature,
  title={Feature Graph Learning for 3D Point Cloud Denoising},
  author={Hu, Wei and Gao, Xiang and Cheung, Gene and Guo, Zongming},
  journal={IEEE Transactions on Signal Processing (TSP)},
  month={March},
  year={2020}
}
```
