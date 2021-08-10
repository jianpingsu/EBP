# EBP
This is an implementation of the paper "Efficient Bijective Parameterizations".

### Dependencies
* [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/download/)
* [Eigen](http://eigen.tuxfamily.org/)
* [Triangle](http://www.cs.cmu.edu/~quake/triangle.html)
* [MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html#gs.8lvw52)
* [Qt](http://download.qt.io/archive/qt/) optional for GUI

### Usage

```
shell.exe input_mesh.obj
```
Note that: the input mesh should be a disk-topology tri-mesh. The parameter `num_procs` should be set according to your CPU cores number.
