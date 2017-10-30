# About
This project was written during my graduation in Mechatronics Engineering.
It has code to solve Poisson equation in a 1x1 square, as well as code to solve
Navier Stokes equations in the same square.
Initially it was written in C++, but after I learned functional programming (OCaml), 
I thought that maybe the right tool for this task would be a programming language that
was more expressive. At the same time, I didn't want to use Python cause I was afraid of
having performance issues due to the several loops that I had to implement. This way I found
Julia and was pretty happy with it. In 2 or 3 days I converted one year of work in C++ to Julia
after having spent little time learning it.
There is one bug in the project that I wasn't able to find out. When I compared values in the mesh
with references from the literature, the values didn't match exactly. There was an error or around 0.5% that I still don't know the cause.

<a name="Vortex"/>
## Vortex
This project is part of the vortex group
- **Old website:** <http://www.vortex.unb.br>
- **New Website:** <http://www.vortexresearchgroup.com/>

<a name="Ferrofluidos"/>
## Ferrofluidos

This is a Scientific Junior Research Project in the area of magnetic fluids currently in development during my graduation (my = Ataias) at University of Brasília. The supervisor is professor [Yuri Dumaresq](http://yuri.mat.unb.br/). The goal is to simulate magnetic fluids in a couple different scenarios. For that, mathematical equations are discretized using finite differences then executed on a computer. The main equations are: Laplace, Poisson and Navier Stokes. They are working well on a cavity and the current work is on the magnetization. The domain of the simulation is a cavity, here this is a square mesh of nondimensional size equals to one that has n points.

<a name="Dependencies"/>
## Dependencies

- **Julia:** <http://julialang.org>
	- [ArgParse](https://github.com/carlobaldassi/ArgParse.jl)
	- [HDF5](https://github.com/JuliaLang/HDF5.jl)
- **Python 3:** <https://www.python.org>
	- [Matploblib](http://matplotlib.org)
	- [NumPy](http://www.numpy.org)
	- [h5py](http://www.h5py.org/)
- **LaTeX:** <http://latex-project.org>

If you wish to run the project in a virtual machine, it can be easily set up checking the vagrant repository available at <https://github.com/ataias/ubuntu-vagrant>. If you are using ubuntu, you can check what programs you actually need to install verifying the shell scripts there.

<a name="How to compile and run"/>
## Compile and run

The simulation code is written in Julia and the graphics are made using [matplotlib](http://matplotlib.org/). Julia compiles it the first time that you use the code, using compilation Just-In-Time. As for Python, code is interpreted.

First off, clone the repository:

	git clone https://github.com/ataias/ferrofluidos.git

Before running, the line with the repository path should be added to your `~/.juliarc.jl`. If the repository was cloned in your home folder, you should do:

	push!(LOAD_PATH,ENV["HOME"] * "/ferrofluidos/src")

After that, cd `ferrofluidos/src` to be in the folder with the source code. A simple simulation can be done running `frames.jl`. For that, go to a terminal and run in the folder:

	julia frames.jl 50 10.0 2.5 0.5 0.8 0.8 0.0 -0.05

Details are available in the code and you can also type `julia frames.jl --help` to learn more. The output file is HDF5. Once the program finishes, you can process it calling `read_and_plot.py` in the same folder where the simulation files, .h5~ and .txt~ files are. If they are in the `src` folder, you can run

	chmod +x read_and_plot.py
	./read_and_plot.py

If you have any questions, don't hesitate to contact me.
