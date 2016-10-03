This repository simulate how real-time electronic states evolve after unsetting
and initial bias voltage between left and right side of a nano-capacity / 
junction like structure.

The following heatmap shows this discharging current of a bilayer graphene 
nano-capacitor. Graphene layer A is located at z = 10 A while graphene layer
B is located at z = 13.5 A. At time = 0 atto, an equilibrium ground state is
achieved with the left layer gated at -Vbias/2 = -0.5 V and the right layer
gated at Vbias/2 = 0.5 V. After unsetting this bias voltage, a discharging 
current flows back to the right electrode. Please see the RT-tddft README for
complete info.

![Alt text](rt-tddft/doc/mapCur.png?raw=true)

# RT-tddft README

This repository contains the Real Time-Time Dependent Density Functional
Theory (rt-tddft) extension and the Quantum Espresso code on top of it.

- **Installation**

	To install rt-tddft, please compile the quantum espresso code (QE-5.2.0)
	first.

		unzip espresso.zip
		cd espresso
		./configure --prefix=[YOUR QE PREFIX]
		make pw

	Parallelization is automatically enabled in QE, however you can choose
	to disable it. Please see ./configure --help for more info. The RT-tddft
	code inherits QE's parallelization option by including make.sys.

		cd rt-tddft
		autoconf
		./configure --with-qepath=[YOUR QE PREFIX]
		make build

- **Ground State Calculation**

	Biased ground state are computed by the pwscf routine, fictious dipole
	moment is cancelled with a mirror image of the computational cell. please
	see the graphene example GrBias.pw-in for how to setup the bias voltage
	and the mirror image.

	![Alt text](rt-tddft/doc/gsRho.png?raw=true)

- **Time Propagation**

	The initial bias voltage can be unset with a linear decay formula. Time 
	integration is carried out with the following schemes which all strictly
	conserve the charge. For this time-dependent Hamiltonian:

		'CN'    : Crank-Nicolson, O(dt^1) local error, unconditionally stable
		'CN2'   : second order Crank-Nicolson, O(dt^2) local error
		'CN-mid': mid point Crank-Nicolson, O(dt^2) local error, unconditionally stable

	The integration schemes are driven by two different linear solvers:

		'itsolver': Iterative solver
		'cgsolver': Conjugate Gradient square solver, implemented by Xiaofeng Qian at MIT

- **Restart**

	To restart a terminated calculation, please specify *init_step* as the latest
	step number output from last calculation. Please also make sure 
	that other parameters including the bias potential, solver, method are all
	the same except *num_step*.

- **Dump**

	The input parameter *dump* and *dump_dir* can be used to dump real-time charge
	density and Kohn-Sham potential along z direction to a specified directory. This
	is enabled by default and can be unset with *dump*=.false. The dump files are loaded
	back into a python library for post-processing.

- **Post-Processing**

	Post-processing are handled with C++ and then wrapped into a python module --- tdpost.
	To install the post-processing tool please have Cython and numpy installed first:
	
		cd tools
		./install.sh

# Quantum Espresso README

This is the distribution of the Quantum ESPRESSO suite of codes (ESPRESSO: 
opEn-Source Package for Research in Electronic Structure, Simulation, 
and Optimization), promoted by the IOM-DEMOCRITOS National Simulation Center 
of the Italian CNR (http://www.democritos.it). 

("make" alone prints a list of acceptable targets). Binaries go in bin/.
For more information, see the general documentation in directory Doc/, 
package-specific documentation in Doc/, and the web site
http://www.quantum-espresso.org/

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
