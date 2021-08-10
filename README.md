# EBP
This is an implementation of the paper "Efficient Bijective Parameterizations".

### Dependencies
* [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/download/)
* [Eigen](http://eigen.tuxfamily.org/)
* [Triangle](http://www.cs.cmu.edu/~quake/triangle.html)
* [CGAL](https://www.cgal.org/download.html)
* [PARDISO](https://pardiso-project.org/)
* [Qt](http://download.qt.io/archive/qt/) optional for GUI

### Usage

```
shell.exe input_mesh.obj
```
Note that: the input mesh should be a disk-topology tri-mesh. The parameter `num_procs` should be set according to your CPU cores number.
